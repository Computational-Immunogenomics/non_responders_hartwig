source(paste0(dirname(dirname(getwd())),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))

library(patchwork)

go <- readRDS(paste0(SHARE_DIR, "3_ready.rds")) 

share_output <- 
go %>% 
 tm(cohort = cohortGo,
    feature = clean_name,
    non_responders_given_event = e_nr,
    responders_given_event = e_r,
    non_responders_given_no_event = ne_nr,
    reponders_given_no_event = ne_r,
    total_patients,
    odds_ratio = or, 
    odds_ratio_ci95_high = ci_high,
    odds_ratio_ci95_low = ci_low,
    fisher_exact_pvalue = fisher_pval,
    fisher_exact_pvalue_fdr_adjusted = p_fdr_fisher,
    pfs_hazard_est = surv_est,
    pfs_hazard_ci95_high = surv_high,
    pfs_hazard_ci95_low = surv_low,
    cox_ph_pvalue = surv_pval, 
    cox_ph_pvalue_fdr_adjusted = p_fdr_surv,
    treatment, 
    mechanism) %>%
 ar(responders_given_event, desc(non_responders_given_event))

fwrite(share_output, paste0(FIG_DIR,  "Supplemental_Table_3.csv"))
