source(paste0(dirname(dirname(getwd())),'/map.r'))

source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "clinical_help.r"))

library(gtsummary)
library(gt)

clin <- fread(paste0(READY_DIR, "clinical_ready.csv"))

cohorts <- 
fread("/mnt/bioinfnas2/immunocomp/shared_reference_data/cohorts/cohorts_ready.csv") %>% 
 se(sampleId, cohort) %>% 
 mu(cohort = ifelse( cohort %in% c("Colon", "Rectum"), "Colorectum", cohort))

locations <- 
fread(paste0(SHARE_DIR, "treatments.csv")) %>% 
 fi(cohort == "Pan-Cancer") %>%
 tm(sampleId, 
    location = primaryTumorLocation)#,
    #treatment = derived_treatmentName, 
    #mechanism = derived_treatmentMechanism, 
    #type = derived_treatmentType)

meta <- fread( paste0( META_DIR, "metadata_update_feb21_2025.csv"))
pretreatments <- meta %>% tm(sampleId, pretreated = tolower(hasSystemicPreTreatment))

rna <- 
fread(paste0(SHARE_DIR, "biomarkers_base.csv")) %>%
 tm(sampleId, rna = as.numeric(!is.na(rna_geneset_KEGG_N_GLYCAN_BIOSYNTHESIS)))

go <- fread(paste0(SHARE_DIR, "go.csv")) %>% se(sampleId, group, cohortGo, non_response)

treatments <- 
go %>%
 fi(group == "treatment", grepl("Pan-Cancer", cohortGo)) %>%
 mu(treatment = gsub("Pan-Cancer / ", "", cohortGo))

mechanisms <- 
go %>%
 fi(group == "mechanism", grepl("Pan-Cancer", cohortGo)) %>%
 mu(mechanism = gsub("Pan-Cancer / ", "", cohortGo)) %>%
 se(-cohortGo, -non_response, -group)

types <- 
go %>%
 fi(group == "type", grepl("Pan-Cancer", cohortGo)) %>%
 mu(type = gsub("Pan-Cancer / ", "", cohortGo)) %>%
 se(-cohortGo, -non_response, -group)

#biopsies_go <- 
#c("Liver (Metastases)", "Lymph node (Metastases)", "Primary Tumor Location", "Bone (Metastases)",
#  "Peritoneum/omentum (Metastases)", "Skin (subcutaneous) (Metastases)", "Lung (Metastases)")
biopsies_go <- 
c("Liver (Metastases)", "Lymph node (Metastases)", "Primary Tumor Location")

locations_go <- 
c("Breast", "Lung", "Skin", "Prostate", 
  "Soft tissue", "Pancreas", "Ovary", "Bladder", "Kidney", 
  "Colorectum", "Gastroensophageal")

treatments_go <- 
go %>% 
 gb(cohortGo, group) %>% 
 su(ct = n(), no_dcb = sum(non_response), dcb = ct - no_dcb) %>% 
 fi(ct > 30, no_dcb >=  15, dcb >=  15) %>% 
 ug() %>%
 fi(group %in% c("treatment", "mechanism", "type"), grepl("Pan-Cancer", cohortGo)) %>%
 mu(treatment = gsub("Pan-Cancer / ", "", cohortGo)) %>%
 pu(treatment)

step_a <- 
clin %>%
 tm(sampleId, 
    study = clin_study, 
    age = clin_age,
    sex = clin_sex, 
    pretreatments = pretreatment_lines_group, 
    biopsy = biopsyStructure, 
    biopsy_type = clin_biopsy_type, 
    durableClinicalBenefit) %>% 
 ij(locations, by = "sampleId") %>%
 lj(treatments, by = "sampleId") %>%
 lj(mechanisms, by = "sampleId") %>% 
 lj(types, by = "sampleId") %>% 
 lj(cohorts, by = "sampleId") %>%
 lj(pretreatments, by = "sampleId") %>%
 lj(rna, by = "sampleId")

step_b <- 
step_a %>%
 mu(
 biopsy_go = 
 case_when( 
  biopsy_type == 0 ~ "Primary Tumor Location", 
  location == "Skin" & grepl("Skin", biopsy) ~ "Primary Tumor Location",
  TRUE ~ paste0(biopsy, " (Metastases)")),
 biopsy_metastases =   
 case_when( 
  biopsy_type == 0 ~ "Primary Tumor Location", 
  location == "Skin" & grepl("Skin", biopsy) ~ "Primary Tumor Location",
  TRUE ~ "Metastases")) %>% 
 mu(biopsy_go = ifelse(biopsy_go %in% biopsies_go, biopsy_go, "Other location (Metastases)")) %>%
 mu(
  location_go = 
  case_when(
    location %in% c("Colon", "Rectum") ~ "Colorectum", 
    location == "Esophagus/gastroesophageal junction" ~ "Gastroesophageal",   
    location %in% locations_go ~ location,
    TRUE ~ "Other location")) %>%
 mu(treatment_go = ifelse(treatment %in% treatments_go, treatment, "Other treatment")) %>%
 mu(mechanism_go = ifelse(mechanism %in% treatments_go, mechanism, "Other mechanism")) %>%
 mu(type_go = ifelse(type %in% treatments_go, type, "Other type")) %>%
 se(-cohortGo, -treatment, -biopsy, -location, -non_response, -cohort, -group, -biopsy_type)

table_base <- 
rbind(step_b, step_b %>% mu(study = "Pooled")) %>%
 gb(study) %>% 
 mu(study_nice = paste0(study, " (n = ", n(), ")"), total = n()) %>% 
 ug() %>% 
 se(-study)

age_report <- 
table_base %>%
 gb(study_nice) %>%
 su( q1_age = round(quantile(age, c(.25), na.rm = TRUE)),
     q2_age = round(quantile(age, c(.5), na.rm = TRUE)),
     q3_age = round(quantile(age, c(.75), na.rm = TRUE))) %>%
 mu(age_report = paste0( q2_age , " (",q1_age, "-", q3_age,")")) %>%
 se(study_nice, age_report) %>%
 sp(study_nice, age_report) %>%
 mu(gp = "age", val = "Age, median (IQR)") %>%
 relocate(gp)

other_report <- 
table_base %>%
 se(-age) %>% 
 ga(gp, val, -sampleId, -study_nice, -total) %>%
 gb(study_nice, total, gp, val) %>%
 su(ct = n(), pct = ct/total, show = paste0(ct, " (", round(100*pct), ")")) %>%
 ug() %>% 
 unique() %>% 
 se(-total, -pct,-ct) %>%
 sp(study_nice, show)

order_gps <- 
c("durableClinicalBenefit", "age", "sex", "pretreated", "rna", 
  "biopsy_go", "location_go", "treatment_go", "mechanism_go", "type_go")

counter <- function(x) as.numeric(strsplit(x, " \\(")[[1]][1])

together_report <- 
rbind(age_report, other_report) %>% 
 relocate(val) %>% 
 relocate(gp) %>%
 fi(gp != "biopsy_metastases", 
    !(gp == "pretreated" & val == "no"), 
    !(gp == "durableClinicalBenefit" & val == 0),
    !((gp == "sex" & val == 0)), 
    gp != "pretreatments", 
    !(gp == "rna" & val == 0)) %>%
 mu(gp = factor(gp, levels = order_gps)) %>%
 rw() %>% mu(count = counter(`Pooled (n = 2594)`)) %>% ug() %>%
 ar(gp, desc(count)) 

together_report %>%
 mu( val = case_when(
    gp == "durableClinicalBenefit" ~ "Durable Clinical Benefit, n (%)",
    gp == "age" ~ "Age, median (IQR)",
    gp == "sex" ~ "Female, n (%)", 
    gp == "pretreated" ~ "Systemic Pretreatment, n (%)", 
    gp == "rna" ~ "RNA Available, n (%)", 
    TRUE ~ val)) %>%
 drop_na(gp) %>%
 se(-gp, -count) %>%
 mutate(across(everything(), ~ replace_na(.x, "0 (0)")))
