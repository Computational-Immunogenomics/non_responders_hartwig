source(paste0(dirname(dirname(getwd())),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))

library(patchwork)
library(scales)
library(ggrepel)
library(ggcorrplot)
library(corrplot)
library(survminer)
library(cowplot)

extra_theme <- 
theme(axis.text.x = element_text(angle = 0, size = 12), 
      axis.text.y = element_text(size = 12), 
      plot.title = element_text(size = 16),
      plot.margin = unit(c(1, 1, 1, 0), "cm")) 

type_mapping <- c(
  "Anti-PD-1" = "Immunotherapy",
  "Pyrimidine (ant)agonist" = "Chemotherapy",
  "Hormonal therapy" = "Hormonal therapy",
  "Chemotherapy" = "Chemotherapy",
  "Multikinase inhibitor" = "Targeted therapy",
  "Alkaloid / Platinum" = "Chemotherapy",
  "Aromatase inhibitor" = "Hormonal therapy",
  "Platinum / Pyrimidine (ant)agonist" = "Chemotherapy",
  "Targeted therapy" = "Targeted therapy",
  "Anti-AR" = "Hormonal therapy",
  "Immunotherapy" = "Immunotherapy")

base <- 
readRDS(paste0(SHARE_DIR, "3_ready.rds")) %>% 
 rw() %>%
 mu( treatment_type = factor(
     ifelse(low_responder & fisher_pval < .01, type_mapping[mechanism], "Not highlighted"),
     levels = c("Hormonal therapy", "Immunotherapy", "Targeted therapy", "Chemotherapy", "Not highlighted"))) %>%
 ug() %>%
 mu(`FDR adjusted significance` = gp2) %>% 
 mu(`Response Rate < 5%` = low_responder) 

my_colors <- c("#F04437", "#E81F64", "#903E97", "#65499E", "#4356A5", "#478FCC", "#34A4DD", "#00BCD4", "#009889", "#4BB04F", "#8BC34C", "#CCDA3A", "#FCED3A", "#FFC10E", "#F8991D", "#F1592C", "#7A5649", "#9F9E9E", "#607F8C")

fill_map <- 
list(
'Not highlighted' = my_colors[18],
'Immunotherapy' = my_colors[1], 
'Hormonal therapy' = my_colors[15],    
'Targeted therapy' = my_colors[5],        
'Chemotherapy' = my_colors[10]
) 

shapes_map <- c("Fisher Significant" = 24, "PFS Significant" = 21, "Rest" = 22)

highlights_main <- 
base %>%
 fi(cohortGo %in% c("Pan-Cancer / Anti-PD-1", "Pan-Cancer / Targeted therapy")) %>%
 fi(feature %in% c("purity_tmbStatus_low", "purity_tmbStatus_high", "hotspot_BRAF_V600E", "purity_msStatus_high")) %>%
 ar(prob_response) %>%
 se(cohortGo, feature, p_fdr_fisher, prob_response) %>%
 gb(feature) %>% mu(rk = row_number(p_fdr_fisher)) %>% fi(rk == 1) %>% se(-rk) %>%
 mu(example = 
     case_when(feature == "hotspot_BRAF_V600E" ~ "BRAF V600E\nPan-Cancer / Targeted Therapy",
               feature == "purity_tmbStatus_high" ~ "TMB High\nPan-Cancer / Anti-PD-1", 
               feature == "purity_tmbStatus_low" ~ "TMB Low\nPan-Cancer / Anti-PD-1",
               feature == "purity_msStatus_high" ~ "MSI High\nPan-Cancer / Anti-PD-1"))

main_figure <- 
 base %>% 
   ggplot( 
    aes(x = prob_response, 
        y = -log10(p_fdr_fisher),
        alpha = `Response Rate < 5%`, 
        color = `Response Rate < 5%`,
        size = `Response Rate < 5%`,
        shape = `FDR adjusted significance`, 
        fill = treatment_type, 
        label = example)) + 
   geom_point() + 
   scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = .2)) + 
   scale_size_manual(values = c("TRUE" = 4, "FALSE" = 2)) +
   scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey")) +
   scale_fill_manual(values = fill_map) +
   scale_shape_manual(values = shapes_map) +
   go_theme + 
   extra_theme + 
   geom_hline(yintercept = 1, alpha = .7, color = "red") + 
   geom_vline(xintercept = .05, alpha = .7, size = .1) + 
   labs(x = "Estimated Probability of Response", 
        y = "Fisher's Test\n-Log10 (FDR Adjusted p-value)", 
        title = "Systematic Analysis Results") + 
   guides(alpha = "none", 
          shape = guide_legend(override.aes = list(size = 4)),
          fill = guide_legend(override.aes = list(shape = 21, size = 4))
         ) + 
   scale_x_continuous( labels = percent_format(accuracy = 1), limits = c(-.01,1)) +
   annotate("text", x = 0, y = 4.5, label = "Non-\nresponse\nfeatures", size = 2.7, fontface = "bold")

main_figure_with_highlights <- 
main_figure +
 geom_text(data = highlights_main, 
            aes(x = prob_response, y =  -log10(p_fdr_fisher) + .15, label = example),inherit.aes = FALSE, color = "black", size = 3)

options(repr.plot.width = 10, repr.plot.height = 7)
main_figure_with_highlights
ggsave( paste0(FIG_DIR, "volcano_highlight_positive.png"), plot = main_figure_with_highlights, width = 10, height = 7)
ggsave( paste0(FIG_DIR, "volcano_highlight_positive.pdf"), plot = main_figure_with_highlights, width = 10, height = 7)

remove <- c("Pan-Cancer / Pazopanib\nRNA Arachidonic Acid Metabolism Very Low",
            "Lung NSCLC / Chemotherapy\nDrivers Pathway DDR", 
            "Driver KRAS Pan-Cancer / Targeted therapy")

zoom_df <- base %>% fi(selected_example) %>% ar(gp) %>% fi(!example %in% remove, feature != "driver_KRAS")

highlight <- 
main_figure + 
 scale_x_continuous( labels = percent_format(accuracy = 1), breaks = c(0,.05), limits = c(-.03,.051))  + 
 ylim(.47,1.41) + 
 theme(axis.title.y = element_blank(), legend.key.size = unit(1.5, "lines") ) + 
 labs(x = "Estimated Probability of Response", 
      y = "Fisher's Test of Odds Ratio\n-Log10 (p-value)", 
      title = "Non-response Features") +
 geom_text_repel(data = zoom_df, 
                 aes(label = example), size = 3,  nudge_y = .07,
                    force = 1,               # increase repulsion force (default = 1)
                    force_pull = 0.01,        # reduce attraction toward anchor point (default = 0.1)
                    box.padding = .7,
                    point.padding = 1,
                 max.overlaps = Inf) 

go <- plot_grid(main_figure + 
                theme(legend.position = "none"), 
                highlight, ncol = 2, 
                rel_widths = c(1, 1.25), 
                labels = c("a", "b"), 
                label_size = 18)

options(repr.plot.width = 15, repr.plot.height = 7)
share <- 
go + 
plot_annotation(
    title = "Non-Response - Systematic Testing Results",
    subtitle = "2,663 biomarkers tested for Non-Response association (Fisher's Exact, Cox-PH) across across 56 cohorts",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      plot.caption = element_text(size = 10, face = "italic"),
      plot.tag = element_text(face = "bold", size = 14)  
    )
  )
share

ggsave( paste0(FIG_DIR, "volcano_main.png"), plot = share, width = 15, height = 7)

saveRDS(share, paste0(FIG_DIR, "volcano_main.rds"))
