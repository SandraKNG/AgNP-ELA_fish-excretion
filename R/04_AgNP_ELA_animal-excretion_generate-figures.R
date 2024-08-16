  # *** Whole-lake silver nanoparticles addition promotes phosphorus and 
  # silver excretion by yellow perch (Perca flavescens)  ***
  # this code was developped by S. Klemet-N'Guessan in 2020-2023
  
  # load libraries ----
  library(gt) # to make tables
  library(ggpubr) # to add multiple graphs on the same page + stats
  library(fishflux) # to do fish stoich model
  library(parallel) # for mcmapply function for sensitivity analysis
  library(fishualize) # for colours based on fish colours
  library(tidybayes)
  library(patchwork) # to align multiple plots
  library(cowplot) # to align multiple plots
  library(scales) # for log scale formatting
  
  # make tables 
  # Table S1 ----
  combined_anova %>%  
    gt(groupname_col = "groupname") %>% 
    cols_label(
      Predictors = "",
      Df = md("**df**"),
      'Sum Sq' = md("**SS**"),
      'Mean Sq' = md("**MS**"),
      'F value' = md("**F**"),
      'Pr(>F)' = md("***p***")
    ) %>% 
    cols_align(
      align = "center",
      columns = c(Df, 'Sum Sq', 'Mean Sq', 'F value', 'Pr(>F)')
    ) %>% 
    fmt_number(
      columns = c('Sum Sq', 'Mean Sq', 'F value'),
      decimals = 3
    ) %>% 
    fmt_number(
      columns = c('Pr(>F)'),
      decimals = 2
     ) %>% 
  tab_style(
    style = cell_borders(
      sides = c("top", "bottom"),
      style = 'hidden'
    ),
    locations = cells_body(
      columns = everything(),
      rows = everything()
    )
  ) %>% 
  gtsave("tables_figures/final-tables_figures/tableS1.rtf")
  
  # Table S2 ----
  contrasts %>% 
    gt(groupname_col = "test") %>% 
    fmt_number(
      columns = everything(),
      decimals = 3
    ) %>% 
    cols_align(
      align = "center"
    ) %>% 
    cols_align(
      align = "left",
      columns = contrast
    ) %>% 
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    ) %>% 
    cols_label(
      lower.CL = 'lower CI',
      upper.CL = "upper CI",
      t.ratio = "t ratio",
      p.value = md("*p*")
    ) %>% 
    gtsave("tables_figures/final-tables_figures/tableS2.rtf")
  
  # make tables 
  # Table S1 ----
  combined_anova %>%  
    gt(groupname_col = "groupname") %>% 
    cols_label(
      Predictors = "",
      Df = md("**df**"),
      'Sum Sq' = md("**SS**"),
      'Mean Sq' = md("**MS**"),
      'F value' = md("**F**"),
      'Pr(>F)' = md("***p***")
    ) %>% 
    cols_align(
      align = "center",
      columns = c(Df, 'Sum Sq', 'Mean Sq', 'F value', 'Pr(>F)')
    ) %>% 
    fmt_number(
      columns = c('Sum Sq', 'Mean Sq', 'F value'),
      decimals = 3
    ) %>% 
    fmt_number(
      columns = c('Pr(>F)'),
      decimals = 2
     ) %>% 
  tab_style(
    style = cell_borders(
      sides = c("top", "bottom"),
      style = 'hidden'
    ),
    locations = cells_body(
      columns = everything(),
      rows = everything()
    )
  ) %>% 
  gtsave("tables_figures/final-tables_figures/tableS1.rtf")
  
  # Table S2 ----
  contrasts %>% 
    gt(groupname_col = "test") %>% 
    fmt_number(
      columns = everything(),
      decimals = 3
    ) %>% 
    cols_align(
      align = "center"
    ) %>% 
    cols_align(
      align = "left",
      columns = contrast
    ) %>% 
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    ) %>% 
    cols_label(
      lower.CL = 'lower CI',
      upper.CL = "upper CI",
      t.ratio = "t ratio",
      p.value = md("*p*")
    ) %>% 
    gtsave("tables_figures/final-tables_figures/tableS2.rtf")
  
  # Table S3 ----
  MDD_results %>% 
    gt() %>% 
    gtsave("tables_figures/final-tables_figures/tableS3.rtf")
  
  # set up plot parameters ----
  lake.labels <- c('Experimental', 'Reference')
  exp.labels <- c('Pre-addition', 'Year 1', 'Year 2', 'Post-addition')
  exp_red.labels <- c('Pre-addition', 'Year 2', 'Post-addition')
  lake.colors <- c("gray60", "black")
  point.size = 1.5
  point.size2 = .5
  point.alpha = .5
  CI.alpha = .25
  line.size = .75
  jitter.position = position_jitterdodge(jitter.width = 0.2)
  
  get_legend <- function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
    }
  
  plot_emm <- function(df, xinter) {
    ggplot(df, aes(x = Year, y = response, color = Lake, group = Lake)) +
      geom_point(size = point.size, position = position_dodge(0.5)) +
      scale_x_discrete(labels = exp_red.labels) +
      geom_errorbar(aes(ymax = lower.CL, ymin = upper.CL), 
                    width = 0.2, lwd = 0.5, position = position_dodge(0.5)) +
      geom_vline(xintercept = xinter, linetype = 'dashed', linewidth = .2) +
      theme_classic(base_size = 10) +
      scale_colour_manual(name = 'Lake',
                          labels = lake.labels,
                          values = lake.colors) 
  }
  
  plot_effsize <- function(df) {
    ggplot(df, aes(x = Year, y = exp(estimate))) +
    geom_pointrange(aes(ymax = exp(lower.CL), ymin = exp(upper.CL)), lwd = 0.4) +
      labs(x = '',
           y = 'Effect size') +
      scale_x_discrete(labels = exp_red.labels) +
      theme_classic(base_size = 10) 
  }
  
  plot_wilcox <- function(y){
    ggplot(excr.Tag, 
           aes(x = Year, y = y)) +
      geom_boxplot(outlier.shape = NA, size = 0.5, colour = "gray60") +
      geom_jitter(size = point.size, position = position_jitter(0.2),
                  alpha = 0.5, colour = "gray60") +
      stat_compare_means(label.x = 0.9, label = 'p.format', 
                         size = 3) +
      theme_classic(base_size = 10) +
      scale_x_discrete(labels = c('Year 1', 'Year 2')) 
  }
  
  plot_mS <- function(y, lim1, lim2){
    ggplot(group_by(iter_yp, iter), aes(x = dry.mass, y = y)) +
      stat_lineribbon(alpha = 0.8,
                      linewidth = line.size,
                      show.legend = F) +
      scale_fill_brewer() +
      theme_classic(base_size = 10) +
      scale_colour_manual(name = 'Lake',
                          labels = lake.labels,
                          values = lake.colors) +
      scale_shape_manual(
        labels = exp.labels,
        values = c(16, 8, 13, 17),
        na.translate = F
      ) +
      scale_y_log10(
        breaks = trans_breaks("log10", function(x)
          10 ^ x),
        labels = trans_format("log10", math_format(10 ^ .x)),
        limits = c(lim1, lim2)
      )
      
      # coord_trans(y = 'log10', ylim = c(lim1, lim2))
  }
  
  plot_mVM <- function(y, y_se, lim1, lim2){
    ggplot(group_by(iter_yp_VM, iter), 
           aes(x = Mass, y = y)) +
      geom_line(size = line.size) +
      geom_ribbon(aes(ymin = y - y_se,
                      ymax = y + y_se), alpha = .1) +
      theme_classic(base_size = 10) +
      scale_colour_manual(name = 'Lake',
                          labels = lake.labels,
                          values = lake.colors) +
      scale_shape_manual(labels = exp_red.labels,
                         values = c(16, 8, 13, 17), na.translate = F) +
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x)), 
                    limits = c(lim1, lim2))
      # coord_trans(y = 'log10', ylim = c(lim1, lim2))
  }
  
  # Figure 1 ----
  # N excretion
  Nexcr.p <- plot_emm(pwcN[["emmeans"]], 2.5) +
    geom_point(data = NPexcr %>% filter(Year != '2014'),
               aes(x = Year, y = massnorm.N.excr),
               size = point.size2, alpha = point.alpha,
               position = jitter.position) +
    labs(x = '',
         y = 'Mass-specific \n N excretion (μg N/g/h)') 
  Nexcr.p
  
  eff_sizeN.p <- plot_effsize(pwcN[["effect_sizes"]]) 
  eff_sizeN.p
  
  # P excretion
  Pexcr.p <- plot_emm(pwcP[["emmeans"]], 3.5) +
    geom_point(data = NPexcr,
               aes(x = Year, y = massnorm.P.excr),
               size = point.size2, alpha = point.alpha,
               position = position_jitterdodge(jitter.width = .2)) +
    scale_x_discrete(labels = exp.labels) +
    labs(x = '',
         y = 'Mass-specific \n P excretion (μg P/g/h)') 
  Pexcr.p
  
  eff_sizeP.p <- ggplot(pwcP[["effect_sizes"]], aes(x = Year, y = exp(-estimate))) +
    geom_pointrange(aes(ymax = exp(-lower.CL), ymin = exp(-upper.CL)), lwd = 0.4) +
    scale_x_discrete(labels = exp.labels) +
    labs(x = '',
         y = 'Effect size') +
    theme_classic(base_size = 10) 
  eff_sizeP.p
  
  # N:P excretion
  # N excretion
  NPexcr.p <- plot_emm(pwcNP[["emmeans"]], 2.5) +
    geom_point(data = NPexcr %>% filter(Year != '2014'),
               aes(x = Year, y = massnorm.NP.excr),
               size = point.size2, alpha = point.alpha,
               position = position_jitterdodge(jitter.width = .2)) +
    labs(x = '',
         y = 'Mass-specific \n N:P excretion (molar)') 
  NPexcr.p
  
  eff_sizeNP.p <- plot_effsize(pwcNP[["effect_sizes"]]) +
    labs(x = '',
         y = 'Effect size') 
  eff_sizeNP.p
  
  # combine plots ----
  ggarrange(Nexcr.p, eff_sizeN.p, Pexcr.p, eff_sizeP.p, NPexcr.p, 
            eff_sizeNP.p, ncol = 2, nrow = 3, 
            labels = c("(a)", "(b)", "(c)", "(d)","(e)", "(f)"),
            font.label = list(size = 10), label.x = 0.2, label.y = 1,
            legend = 'top', common.legend = T, align = 'v')
  ggsave('tables_figures/final-tables_figures/Fig1.tiff', 
         width = 7, height = 7, 
         units = 'in', dpi = 600, compression = 'lzw')
  
  # Figure 2 ----
  # Tag
  TAgexcr.p <- plot_wilcox(excr.Tag$massnorm.Tag.excr) +
  labs(x = '', 
       y = 'Mass-specific \n Ag excretion (μg Ag/g/h)') +
    theme(axis.text.x = element_blank())
  TAgexcr.p
  # N:Ag
  NAgexcr.p <- plot_wilcox(excr.Tag$massnorm.NAg.excr)+
    labs(x = '', 
         y = 'Mass-specific \n N:Ag excretion (molar)')+
    theme(axis.text.x = element_blank())
  NAgexcr.p
  # P:Ag
  PAgexcr.p <- plot_wilcox(excr.Tag$massnorm.PAg.excr) +
    labs(x = '', 
         y = 'Mass-specific \n P:Ag excretion (molar)')
  PAgexcr.p
  
  # combine plots ----
  ggarrange(TAgexcr.p, NAgexcr.p, PAgexcr.p, ncol = 1, nrow = 3, 
            labels = c("(a)", "(b)", "(c)"), 
            font.label = list(size = 10), label.x = 0.2, label.y = 1,
            legend = 'none', common.legend = T, align = 'v')
  ggsave('tables_figures/final-tables_figures/Fig2.tiff', 
         width = 3.33, height = 7, 
         units = 'in', dpi = 300)
  
  # ..Figure 3 ----
  # Schiettekate et al model ----
  # N excretion 
  Nexcr_mS.p <- plot_mS(iter_yp$N.excretion.mod, 10^0.5, 10^3.3) +
    geom_point(aes(x = Mass, y = N.excretion.rate, color = Lake, shape = Year), 
               data = NPexcr %>%  filter(Year != '2014'),
               size = point.size, alpha = point.alpha + .2) +
    labs(x = "", y = "N excretion (μg N/ind/h)") +
    scale_shape_manual(labels = exp_red.labels,
                       values = c(16, 8, 17), na.translate = F) +
    theme(axis.text.x = element_blank())
  Nexcr_mS.p
  
  # P excretion
  Pexcr_mS.p <- plot_mS(iter_yp$P.excretion.mod, 10^-2.5, 10^3) +
    geom_point(aes(x = Mass, y = P.excretion.rate, color = Lake, shape = Year), 
               data = NPexcr,
               size = point.size, alpha = point.alpha + .2) +
     labs(x = "", y = "P excretion (μg P/ind/h)") 
  Pexcr_mS.p
  
  # V&M model ----
  # N excretion
  Nexcr_mVM.p <- plot_mVM(iter_yp_VM$N.excretion.m, 
                          iter_yp_VM$N.excretion.m.se, 10^0.5, 10^3.3) +
    geom_point(aes(x = Mass, y = N.excretion.rate, color = Lake, shape = Year), 
               data = NPexcr %>%  filter(Year != '2014'),
               size = point.size, alpha = point.alpha + .2) +
    scale_shape_manual(labels = exp_red.labels,
                       values = c(16, 8, 17), na.translate = F) +
    labs(x = "", 
         y = "") +
    theme(axis.text.x = element_blank())
  Nexcr_mVM.p
  
  # P excretion
  Pexcr_mVM.p <- plot_mVM(iter_yp_VM$P.excretion.m, 
                          iter_yp_VM$P.excretion.m.se, 10^-2.5, 10^3) +
    geom_point(aes(x = Mass, y = P.excretion.rate, color = Lake, shape = Year), 
               data = NPexcr,
               size = point.size, alpha = point.alpha + .2) +
    labs(x = "", 
         y = "") 
  Pexcr_mVM.p
  
  # combine plots ----
  fig3.legend <- get_legend(Pexcr_mS.p)
  fig3 <- ggarrange(Nexcr_mS.p, Nexcr_mVM.p, Pexcr_mS.p,  Pexcr_mVM.p, 
            nrow = 2, ncol = 2,
            legend = 'right', common.legend = F, align = 'hv', 
            labels = c('(a)', '(b)', '(c)', '(d)'), 
            label.x = 0.17, label.y = 1, font.label = list(size = 8), 
            legend.grob = fig3.legend)
<<<<<<< HEAD
  annotate_figure(fig3, 
                  bottom = text_grob('Dry mass (g)', size = 10, y = 1))
  
  ggsave('tables_figures/final-tables_figures/Fig3.tiff', 
         width = 7, height =  5, 
         units = 'in', dpi = 600, compression = 'lzw', bg = 'white')
=======
  ggsave('tables_figures/final-tables_figures/Fig3.tiff', 
         width = 7, height =  4, 
         units = 'in', dpi = 600, compression = 'lzw')
>>>>>>> fc56b5e280398a99aa1b7c8f1a834384b58e39e3
  
  # Figure 4 ----
  # plots ----
  lim_ref.p <- 
    ggplot() +
    geom_segment(aes(x = 4.7, y = 4.7/Snp, xend = Dnm, yend = Dpm), size = 0.5) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = Dnm, yend = 6.7/Snp), size = 0.5) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = 6.7, yend = Dpm), size = 0.5) +
    geom_text(aes(x = 5.25, y = 1.25, label = "N"), size = 13) +
    geom_text(aes(x = 6.25, y = 1.3, label = "C"), size = 13) +
    geom_text(aes(x = 6, y = 1.1, label = "P"), size = 13) +
    labs(x = "Diet N (%)", y = "Diet P (%)") +
    theme_bw(base_size = 8) +
    theme(aspect.ratio = 1) 
  lim_ref.p
  
  Cing_ref.p <-
    ggplot(simd) +
    geom_raster(aes(x = Dn, y = Dp, fill = C.ingestion.mod)) +
    scale_fill_fish(option = "Trimma_lantana", trans = "sqrt") +
    geom_segment(aes(x = 4.7, y = 4.7/Snp, xend = Dnm, yend = Dpm), size = 0.5) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = Dnm, yend = 6.7/Snp), size = 0.5) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = 6.7, yend = Dpm), size = 0.5) +
    theme_bw(base_size = 8) +
    theme(aspect.ratio = 1) +
    labs(x = "Diet N (%)", y = "Diet P (%)", 
         fill = "C ingestion \n (μg C/ind/h)")
  Cing_ref.p
  
  Nexcr_ref.p <-
    ggplot(simd) +
    geom_raster(aes(x = Dn, y = Dp, fill = N.excretion.mod)) +
    geom_segment(aes(x = 4.7, y = 4.7/Snp, xend = Dnm, yend = Dpm), size = 0.5) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = Dnm, yend = 6.7/Snp), size = 0.5) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = 6.7, yend = Dpm), size = 0.5) +
    scale_fill_fish(option = "Trimma_lantana", trans = "sqrt") +
    theme_bw(base_size = 8) +
    theme(aspect.ratio = 1)  +
    labs(x = "Diet N (%)", y = "Diet P (%)", 
         fill = 'N excretion \n (μg N/ind/h)')
  Nexcr_ref.p
  
  Pexcr_ref.p <-
    ggplot(simd) +
    geom_raster(aes(x = Dn, y = Dp, fill = P.excretion.mod)) +
    scale_fill_fish(option = "Trimma_lantana", trans = "sqrt") +
    geom_segment(aes(x = 4.7, y = 4.7/Snp, xend = Dnm, yend = Dpm), size = 0.5) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = Dnm, yend = 6.7/Snp), size = 0.5) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = 6.7, yend = Dpm), size = 0.5) +
    theme_bw(base_size = 8) +
    theme(aspect.ratio = 1) +
    labs(x = "Diet N (%)", y = "Diet P (%)", 
         fill = 'P excretion \n (μg P/ind/h)') 
  Pexcr_ref.p
  
  # combine plots ----
  plot_grid(lim_ref.p, Cing_ref.p, Nexcr_ref.p, Pexcr_ref.p, 
            ncol = 2, nrow = 2, align = 'hv', axis = "r", 
            labels = c("(a)", "(b)", "(c)", "(d)"),
            label_size = 8)
  ggsave('tables_figures/final-tables_figures/Fig4.tiff', 
         width = 7, height = 4, 
         units = 'in', dpi = 600)
  
  # Figure S1 ----
  Mexcr.p <- ggplot(NPexcr, 
                   aes(x = Year, y = Mass, color = Lake)) +
    geom_jitter(size = 1.5, aes(x = Year, y = Mass), 
                position = position_jitterdodge(jitter.width = 0.2),
                show.legend = FALSE, alpha = 0.5) +
    theme_classic(base_size = 10) +
    geom_boxplot(outlier.shape = NA, size = 0.5, fill = NA) +
    scale_x_discrete(labels = exp.labels) +
    labs(x = '',
         y = 'Dry mass (g)') +
    scale_colour_manual(name = 'Lake',
                        labels = lake.labels,
                        values = lake.colors)
  Mexcr.p
  
  ggsave('tables_figures/final-tables_figures/FigS1.tiff', 
         width = 7, height = 3.5, 
         units = 'in', dpi = 600, compression = 'lzw')
  
  # Figure S2 ----
  lim_mod.p <- 
    ggplot(limitation, aes(x = dry.mass, y = prop_lim)) +
    geom_point(aes(x = dry.mass, y = prop_lim, color = nutrient), size = 1.5) +
    geom_line(aes(x = dry.mass, y = prop_lim, color = nutrient), size = 0.5) +
    labs(x = "Dry mass (g)", y = "Proportion of iterations", 
         color = "Limiting element") +
    theme_classic(base_size = 10) +
    scale_color_grey(start = 0.8, end = 0.1, labels = c("C", "N", "P")) 
  lim_mod.p
  
  ggsave('tables_figures/final-tables_figures/FigS2.tiff', 
         width = 7, height = 4, 
         units = 'in', dpi = 300, compression = 'lzw')
  
  
  # Figure S3 ----
  Cing_mod.p <- 
    ggplot(group_by(iter_yp, iter), aes(x = dry.mass, y = C.ingestion.mod)) +
    stat_lineribbon(alpha = 0.8, show.legend = T) +
    scale_fill_brewer(name = 'Confidence Interval',
                      labels = c('95%', '80%', '50%')) +
    labs(x = "Dry mass (g)", y = "C ingestion rate (μg C/ind/h)") +
    theme_classic(base_size = 10) 
  Cing_mod.p
  
  ggsave('tables_figures/final-tables_figures/FigS3.tiff', 
         width = 7, height = 3.5, 
         units = 'in', dpi = 600)
  
  # export final tables excel files ----
  write_csv(NPexcr, 'output/excr_final.csv')
  write_csv(NPexcr.ss1, 'output/excr_smry_lake_year.csv')
  