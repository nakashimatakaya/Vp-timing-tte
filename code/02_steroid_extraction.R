################################################################################
# 02_steroid_extraction.R
#
# Extract systemic corticosteroid exposure from MIMIC-IV (eMAR, inputevents)
# and eICU-CRD (medication, infusiondrug).
# Creates per-patient steroid flags relative to time zero.
#
# Requires: pt_vp6, miiv, eicu (ricu sources attached)
# Outputs:  pt_vp6 with steroid_pre and hourly steroid flags appended
################################################################################

library(dplyr)
library(stringr)
library(tidyr)
library(readr)

# ==============================================================================
# 1. Prepare pt_vp6 with row identifiers
# ==============================================================================

pt_base <- pt_vp6 %>%
  ungroup() %>%
  mutate(
    row_id     = row_number(),
    src_std    = case_when(
      str_detect(str_to_lower(as.character(src)), "miiv|mimic") ~ "miiv",
      str_detect(str_to_lower(as.character(src)), "eicu")       ~ "eicu",
      TRUE ~ str_to_lower(as.character(src))
    ),
    all_id_chr = as.character(all_id),
    t0_hr      = as.numeric(t0_hr),
    can_eval   = !is.na(all_id_chr) & !is.na(t0_hr) & src_std %in% c("miiv", "eicu")
  )

pt_ref <- pt_base %>%
  filter(can_eval) %>%
  transmute(row_id, src_std, all_id = all_id_chr, t0_hr)

# ==============================================================================
# 2. Target IDs for each database
# ==============================================================================

target_ids_mimic <- miiv$icustays %>%
  as_tibble() %>%
  semi_join(
    pt_ref %>% filter(src_std == "miiv") %>%
      transmute(stay_id = as.integer(all_id)) %>% distinct(),
    by = "stay_id"
  ) %>%
  transmute(
    all_id = as.character(subject_id),
    pt_all_id = as.character(stay_id),
    subject_id, hadm_id, stay_id, intime
  ) %>%
  distinct()

target_ids_eicu <- pt_ref %>%
  filter(src_std == "eicu") %>%
  distinct(all_id) %>%
  filter(!is.na(all_id)) %>%
  mutate(pt_all_id = all_id)

# ==============================================================================
# 3. Drug pattern helpers
# ==============================================================================

steroid_pat <- regex(
  paste(c("hydrocortisone", "solu[- ]?cortef",
          "methylpred", "methylprednisolone", "solu[- ]?medrol",
          "prednisone", "prednisolone",
          "dexameth", "dexamethasone", "decadron",
          "betameth", "betamethasone",
          "triamcinolone", "fludrocortisone", "cortisone"),
        collapse = "|"),
  ignore_case = TRUE
)

nonsystemic_pat <- regex(
  paste(c("TOPIC", "TOPICAL", "OPH", "OPHTH", "EYE",
          "OTIC", "EAR", "NASAL", "INHAL", "NEB",
          "SPRAY", "DROP", "CREAM", "OINT", "LOTION",
          "GEL", "FOAM", "VAGIN", "RECTAL", "ENEMA"),
        collapse = "|"),
  ignore_case = TRUE
)

keep_systemic <- function(drug_txt = NA_character_, route_txt = NA_character_) {
  txt <- str_squish(str_to_upper(
    paste(coalesce(as.character(drug_txt), ""),
          coalesce(as.character(route_txt), ""))
  ))
  txt == "" | !str_detect(txt, nonsystemic_pat)
}

is_admin_status <- function(x) {
  z <- str_squish(str_to_upper(coalesce(as.character(x), "")))
  z != "" &
    !str_detect(z, "NOT GIVEN|NOT ADMIN|REFUS|HOLD|MISS|DELAY|CANCEL|REMOVE|STOP|RETURN") &
    str_detect(z, "ADMIN|APPL|CONFIRM|GIVEN|BOLUS")
}

# ==============================================================================
# 4. MIMIC-IV steroid events
# ==============================================================================

# 4a. eMAR
mimic_steroid_emar <- miiv$emar %>%
  as_tibble() %>%
  semi_join(target_ids_mimic %>% select(subject_id, hadm_id) %>% distinct(),
            by = c("subject_id", "hadm_id"), copy = TRUE) %>%
  select(subject_id, hadm_id, emar_id, charttime, medication, event_txt) %>%
  inner_join(miiv$emar_detail %>% as_tibble() %>% select(emar_id, product_description),
             by = "emar_id") %>%
  collect() %>%
  inner_join(target_ids_mimic %>% select(subject_id, hadm_id, pt_all_id, stay_id, intime),
             by = c("subject_id", "hadm_id")) %>%
  mutate(
    drug_txt = str_squish(paste(coalesce(product_description, ""), coalesce(medication, ""))),
    event_offset_hr = as.numeric(difftime(charttime, intime, units = "hours"))
  ) %>%
  filter(str_detect(drug_txt, steroid_pat),
         keep_systemic(drug_txt),
         is_admin_status(event_txt),
         !is.na(event_offset_hr)) %>%
  transmute(src_std = "miiv", all_id = pt_all_id, source_table = "emar",
            drug_name = coalesce(product_description, medication), event_offset_hr)

# 4b. inputevents
mimic_steroid_input <- miiv$inputevents %>%
  as_tibble() %>%
  semi_join(target_ids_mimic %>% select(stay_id) %>% distinct(),
            by = "stay_id", copy = TRUE) %>%
  select(stay_id, starttime, itemid, amount, rate) %>%
  inner_join(miiv$d_items %>% as_tibble() %>% select(itemid, label), by = "itemid") %>%
  collect() %>%
  inner_join(target_ids_mimic %>% select(stay_id, pt_all_id, intime), by = "stay_id") %>%
  mutate(
    drug_txt = str_squish(coalesce(label, "")),
    event_offset_hr = as.numeric(difftime(starttime, intime, units = "hours"))
  ) %>%
  filter(str_detect(drug_txt, steroid_pat),
         keep_systemic(drug_txt),
         !is.na(event_offset_hr),
         coalesce(as.numeric(amount), 0) > 0 | coalesce(as.numeric(rate), 0) > 0) %>%
  transmute(src_std = "miiv", all_id = pt_all_id, source_table = "inputevents",
            drug_name = label, event_offset_hr)

# ==============================================================================
# 5. eICU-CRD steroid events
# ==============================================================================

# 5a. medication table
eicu_steroid_med <- eicu$medication %>%
  as_tibble() %>%
  semi_join(target_ids_eicu %>% transmute(patientunitstayid = as.integer(all_id)) %>% distinct(),
            by = "patientunitstayid", copy = TRUE) %>%
  select(patientunitstayid, drugorderoffset, drugstartoffset,
         drugordercancelled, drugname, routeadmin) %>%
  collect() %>%
  mutate(
    drug_txt = str_squish(coalesce(drugname, "")),
    route_txt = coalesce(routeadmin, ""),
    event_offset_hr = as.numeric(coalesce(drugstartoffset, drugorderoffset)) / 60
  ) %>%
  filter(str_detect(drug_txt, steroid_pat),
         keep_systemic(drug_txt, route_txt),
         str_to_upper(coalesce(drugordercancelled, "NO")) != "YES",
         !is.na(event_offset_hr)) %>%
  transmute(src_std = "eicu", all_id = as.character(patientunitstayid),
            source_table = "medication", drug_name = drugname, event_offset_hr)

# 5b. infusiondrug table
eicu_steroid_inf <- eicu$infusiondrug %>%
  as_tibble() %>%
  semi_join(target_ids_eicu %>% transmute(patientunitstayid = as.integer(all_id)) %>% distinct(),
            by = "patientunitstayid", copy = TRUE) %>%
  select(patientunitstayid, infusionoffset, drugname) %>%
  collect() %>%
  mutate(
    drug_txt = str_squish(coalesce(drugname, "")),
    event_offset_hr = as.numeric(infusionoffset) / 60
  ) %>%
  filter(str_detect(drug_txt, steroid_pat),
         keep_systemic(drug_txt),
         !is.na(event_offset_hr)) %>%
  transmute(src_std = "eicu", all_id = as.character(patientunitstayid),
            source_table = "infusiondrug", drug_name = drugname, event_offset_hr)

# ==============================================================================
# 6. Combine and restrict to relevant time window
# ==============================================================================

need_window <- pt_ref %>%
  group_by(src_std, all_id) %>%
  summarise(need_lo = min(t0_hr, na.rm = TRUE) - 24,
            need_hi = max(t0_hr, na.rm = TRUE) + 6,
            .groups = "drop")

steroid_events <- bind_rows(
  mimic_steroid_emar, mimic_steroid_input,
  eicu_steroid_med, eicu_steroid_inf
) %>%
  distinct(src_std, all_id, source_table, drug_name, event_offset_hr) %>%
  inner_join(need_window, by = c("src_std", "all_id")) %>%
  filter(event_offset_hr >= need_lo, event_offset_hr <= need_hi) %>%
  select(-need_lo, -need_hi)

# ==============================================================================
# 7. Compute per-row steroid features
# ==============================================================================

steroid_features <- pt_ref %>%
  select(row_id, src_std, all_id, t0_hr) %>%
  left_join(steroid_events %>% select(src_std, all_id, event_offset_hr),
            by = c("src_std", "all_id")) %>%
  mutate(t_rel = event_offset_hr - t0_hr) %>%
  group_by(row_id) %>%
  summarise(
    steroid_pre24h_any = as.integer(any(!is.na(t_rel) & t_rel >= -24 & t_rel < 0)),
    steroid_pre24h_closest_hr = {
      x <- -t_rel[!is.na(t_rel) & t_rel >= -24 & t_rel < 0]
      if (length(x) == 0L) NA_real_ else min(x, na.rm = TRUE)
    },
    steroid_post_0_1h = as.integer(any(!is.na(t_rel) & t_rel >= 0 & t_rel < 1)),
    steroid_post_1_2h = as.integer(any(!is.na(t_rel) & t_rel >= 1 & t_rel < 2)),
    steroid_post_2_3h = as.integer(any(!is.na(t_rel) & t_rel >= 2 & t_rel < 3)),
    steroid_post_3_4h = as.integer(any(!is.na(t_rel) & t_rel >= 3 & t_rel < 4)),
    steroid_post_4_5h = as.integer(any(!is.na(t_rel) & t_rel >= 4 & t_rel < 5)),
    steroid_post_5_6h = as.integer(any(!is.na(t_rel) & t_rel >= 5 & t_rel <= 6)),
    steroid_post_0_6h_any = as.integer(any(!is.na(t_rel) & t_rel >= 0 & t_rel <= 6)),
    .groups = "drop"
  )

# ==============================================================================
# 8. Attach to pt_vp6
# ==============================================================================

flag_cols <- c("steroid_pre24h_any", "steroid_post_0_1h", "steroid_post_1_2h",
               "steroid_post_2_3h", "steroid_post_3_4h", "steroid_post_4_5h",
               "steroid_post_5_6h", "steroid_post_0_6h_any")

pt_vp6_steroid <- pt_base %>%
  left_join(steroid_features, by = "row_id")

# Ensure flag columns exist even if steroid_features was empty
for (col in flag_cols) {
  if (!col %in% names(pt_vp6_steroid)) {
    pt_vp6_steroid[[col]] <- NA_integer_
  }
}

pt_vp6_steroid <- pt_vp6_steroid %>%
  mutate(
    across(all_of(flag_cols),
           ~ if_else(can_eval, coalesce(.x, 0L), NA_integer_))
  ) %>%
  select(-row_id, -src_std, -all_id_chr, -can_eval)

pt_vp6 <- pt_vp6_steroid %>%
  mutate(
    steroid_at_t0 = as.integer(coalesce(steroid_pre24h_any, 0L)),
    steroid_0_6h  = as.integer(coalesce(steroid_post_0_6h_any, 0L)),
    steroid_pre   = as.integer(coalesce(as.logical(steroid_at_t0), FALSE))
  )

cat("steroid_pre distribution:\n")
print(table(pt_vp6$steroid_pre, useNA = "always"))
