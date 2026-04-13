source(paste0(dirname(dirname(getwd())),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))
source(paste0(HELP_DIR, "fisher.r"))

library(survminer)
library(rlang)
library(cowplot)
library(patchwork)
library(ggpubr)

pathways <- fread("/home/josephusset@vhio.org/projects/non_responders_share/2_biomarkers/pathways.csv")

#pathways

ddr_drivers <- 
pathways %>% 
 fi(pathway == "DDR") %>%
 mu(tmp = paste0("driver_", gene)) %>%
 pu(tmp)

#ddr_drivers

go <- fread(paste0(SHARE_DIR, "fisher_base.csv"))

examples <- 
fread(paste0(TMP_DIR, "pfs_highlights.csv")) %>% 
 se(cohortGo, feature, example_base) %>% 
 fi(!grepl("Arachidonic", example_base))

#head(go)

ready <- 
go %>%
 rename(trt = derived_treatmentName, val = drivers_pathway_DDR) %>% 
 fi(cohortGo == "Lung NSCLC / Chemotherapy") %>%
 se(sampleId, cohortGo, contains("Pfs"), contains("best"), contains("durable"), val, trt, any_of(ddr_drivers)) %>%
 mu(gp = case_when(
    trt %in% c("Pemetrexed", "Gemcitabine", "Docetaxel") ~ "Non-Platinum Single Agents (Pemetrexed/Gemcitabine/Docetaxel)",
    grepl("Carboplatin ##", trt) ~ "Platinum Doublets (Carboplatin-based)",
    grepl("Cisplatin ##", trt) ~ "Platinum Doublets (Cisplatin-based)"),
    gp = ifelse( (grepl("Cisplatin", trt) | grepl("Carboplatin", trt)), "Platinum Chemotherapies", "Non-Platinum Chemotherapies"))

drivers <- 
ready %>%
 se(gp, contains("driver")) %>% 
 ga(driver, val, -gp) %>%
 fi(val == 1) %>% 
 gb(gp) %>% 
 su(drivers = paste0(sort(gsub("driver_", "", driver)), collapse = ","))

go <- ready %>% lj(drivers, by = "gp") %>% mu(gp2 = paste0(gp, "\nDDR Drivers (", drivers, ")"))

go %>% 
 se(sampleId, contains("driver"), val) %>%
 ar(desc(driver_RAD51B)) %>%
 relocate(val)

data_prep <- function(gp = "Platinum Doublets (Carboplatin-based)") {
    data <- go %>% fi(gp2 == !!gp ) %>% ug()
    surv_formula <- expr(Surv(daysToPfsEvent, pfsEvent) ~ val)
    fits <- survfit(eval(surv_formula, envir = environment()), data = data)
    pval <- signif(survdiff(eval(surv_formula, envir = environment()), data = data)$pvalue,2)
    list("data" = data, "fits" = fits, "surv_formula" = surv_formula, "pval" = pval)
}

figurer <- function(example, data, fits, pval, pad = 30 ) {
        
    max_time <- max(data$daysToPfsEvent, na.rm = TRUE)
    xmax <- min(max_time, 1100)
    
    oo <- 
    ggsurvplot(
     fits,
     data = data,
     palette = c("#7AABD3", "#e52f28"),
     conf.int = TRUE, 
     risk.table = TRUE, 
     pval.coord = c(700, .95), 
     xlim = c(0, xmax), 
     break.time.by = 300, 
     ggtheme = theme_minimal(),
     xlab = "Days", 
     ylab = "Progression Free Survival Probability", 
     title = example) 

    oo$plot <- 
    oo$plot + 
     annotate("text", x = 500, y = 0.9, label = paste0("Log-rank        \np-value = ", pval), size = 5) + 
     theme(plot.title = element_text(hjust = 0.5)) + 
     guides(color = guide_legend(nrow = 1, byrow = TRUE)) 

    as_ggplot(gridExtra::arrangeGrob(oo$plot, oo$table, layout_matrix = matrix(c(1,1,1,1,1,1,1,2,2)))) + 
    theme(plot.margin = margin(t = pad, r = pad, b = pad, l = pad))
}

cohorts <- unique(go %>% pu(gp2))

survival_figures <- list()
for( i in cohorts){
    data_ready <- data_prep(i)
    data <- data_ready$data 
    fits <- data_ready$fits
    pval <- data_ready$pval
    surv_formula <- data_ready$surv_formula
    survival_figures[[i]] <- figurer( i, data = data, fits = fits, pval = pval)
}

options(repr.plot.height = 8, repr.plot.width =16) 
plt_share <- 
wrap_plots(survival_figures, ncol = 2) + 
plot_annotation(
    title = "Lung NSCLC: Progression free survival plots by DDR driver status",
    theme = theme(
      plot.title = element_text(size = 22, face = "bold", hjust = .5)
    )
  )

plt_share
ggsave(paste0(FIG_DIR, "lung_nsclc_ddr_pfs.png"), plt_share, width = 16, height = 8)

base <- fread(paste0(SHARE_DIR, "fisher_base.csv"))

data <- 
base %>%
 fi(grepl("NSCLC", primaryTumorType), grepl("Immunotherapy", derived_treatmentType)) %>% 
 #fi(#grepl("Skin", primaryTumorLocation), grepl("Immunotherapy", derived_treatmentType)) %>% 
 gb(sampleId) %>% mu(rk = row_number()) %>% fi(rk == 1) %>% 
 se(sampleId, cohortGo, contains("Pfs"), derived_treatmentType, derived_treatmentName, primaryTumorLocation) %>% ug() %>%
 mu(cohortGo = "Lung NSCLC / Anti-PD(L)-1")

surv_object <- Surv(time = data$daysToPfsEvent, event = data$pfsEvent)
fit <- survfit(surv_object ~ 1, data = data)

options(repr.plot.width = 6, repr.plot.height = 6)
oo <- ggsurvplot(
  fit,
  data = data,
  conf.int = TRUE,          # optional: show confidence interval
  risk.table = TRUE,        # optional: show risk table
  ggtheme = theme_minimal(), # optional: a clean theme
   break.time.by = 300,    
  xlab = "Days", 
  ylab = "Progression Free Survival Probability", 
  title = "Lung NSCLC / Anti-PD(L)-1 Progression Free Survival"
)

pad = 30 

oo$plot <- 
oo$plot + 
 #annotate("text", x = 500, y = 0.9, label = paste0("Log-rank        \np-value = ", pval), size = 5) + 
 theme(plot.title = element_text(hjust = 0.5)) + 
 guides(color = guide_legend(nrow = 1, byrow = TRUE)) 

oo_go <- as_ggplot(gridExtra::arrangeGrob(oo$plot, oo$table, layout_matrix = matrix(c(1,1,1,1,1,1,1,2,2)))) + 
theme(plot.margin = margin(t = pad, r = pad, b = pad, l = pad))

paste0(FIG_DIR, "lung_nsclc_pd1_pfs.png")

oo_go
ggsave(paste0(FIG_DIR, "lung_nsclc_pd1_pfs.png"), oo_go, width = 6, height = 6)
