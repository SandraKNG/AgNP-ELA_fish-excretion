  # *** Silver nanoparticles addition promotes phosphorus and silver excretion by 
  # yellow perch (Perca flavescens) in a boreal lake ***
  # this code was developped by S. Klemet-N'Guessan in 2020 and 2021
  
  # load libraries ----
  library(ggpubr) # to add multiple graphs on the same page + stats
  library(fishflux) # to do fish stoich model
  library(parallel) # for mcmapply function for sensitivity analysis
  library(fishualize) # for colours based on fish colours
  library(tidybayes)
  library(patchwork) # to align multiple plots
  library(cowplot) # to align multiple plots
  library(scales) # for log scale formatting
  
  # set up plot parameters ----
  lake.labels <- c('AgNPs L222', 'Reference L239')
  exp.labels <- c('Pre-addition', 'Year 1', 'Year 2', 'Post-addition')
  exp_red.labels <- c('Pre-addition', 'Year 2', 'Post-addition')
  lake.colors <- c("gray60", "black")
  point.size = 1.5
  point.size2 = .5
  point.alpha = .5
  
  get_legend <- function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
    }
  
  plot_emm <- function(df) {
    ggplot(df, 
                      aes(x = Year, y = response, color = Lake, group = Lake)) +
      geom_point(size = point.size, position = position_dodge(0.5)) +
      scale_x_discrete(labels = exp_red.labels) +
      geom_errorbar(aes(ymax = lower.CL, ymin = upper.CL), 
                    width = 0.2, lwd = 0.5, position = position_dodge(0.5)) +
      theme_classic(base_size = 10) +
      scale_colour_manual(name = 'Lake',
                          labels = lake.labels,
                          values = lake.colors) 
  }
  
  plot_effsize <- function(df) {
    ggplot(df, aes(x = Year, y = exp(estimate))) +
    geom_pointrange(aes(ymax = exp(lower.CL), ymin = exp(upper.CL)), lwd = 0.4) +
      scale_x_discrete(labels = exp_red.labels) +
      theme_classic(base_size = 10) 
  }
  
  plot_boxplot <- function(y){
    ggplot(excr.Tag, 
           aes(x = Year, y = y)) +
      geom_boxplot(outlier.shape = NA, size = 0.5) +
      geom_jitter(size = point.size, position = position_jitter(0.2),alpha = 0.5) +
      stat_compare_means(label.x = 0.9, label = 'p.format', size = 3) +
      theme_classic(base_size = 10) +
      labs(x = '', 
           y = 'Mass-specific \n Ag excretion (μg Ag/g/h)') +
      scale_x_discrete(labels = c('Year 1', 'Year 2')) 
  }
  
  plot_mod <- function(y){
    ggplot(group_by(iter_yp, iter), aes(x = dry.mass, y = y)) +
      stat_lineribbon(alpha = 0.8, linewidth = 0.6, show.legend = F) +
      scale_fill_brewer() +
      theme_classic(base_size = 8) +
      scale_colour_manual(name = 'Lake',
                          labels = lake.labels,
                          values = lake.colors) +
      scale_shape_manual(labels = exp.labels,
                         values = c(16, 8, 13, 17), na.translate = F) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))
  }
  
  # Figure 1 ----
  # N excretion
  position = position_jitterdodge(jitter.width = 0.2)
  Nexcr.p <- plot_emm(pwcN[["emmeans"]]) +
    geom_point(data = NPexcr %>% filter(Year != '2014'),
               aes(x = Year, y = massnorm.N.excr),
               size = point.size2, alpha = point.alpha,
               position = position_jitterdodge(jitter.width = .2)) +
    labs(x = '',
         y = 'Mass-specific \n N excretion (μg N/g/h)') 
  Nexcr.p
  
  eff_sizeN.p <- plot_effsize(pwcN[["effect_sizes"]]) +
    labs(x = '',
         y = 'Effect size') 
  effsizeN.p
  
  # P excretion
  Pexcr.p <- plot_emm(pwcP[["emmeans"]]) +
    geom_point(data = NPexcr,
               aes(x = Year, y = massnorm.P.excr),
               size = point.size2, alpha = point.alpha,
               position = position_jitterdodge(jitter.width = .2)) +
    scale_x_discrete(labels = exp.labels) +
    labs(x = '',
         y = 'Mass-specific \n P excretion (μg P/g/h)') 
  Pexcr.p
  
  eff_sizeP.p <-  ggplot(pwcP[["effect_sizes"]], aes(x = Year, y = exp(-estimate))) +
    geom_pointrange(aes(ymax = exp(-lower.CL), ymin = exp(-upper.CL)), lwd = 0.4) +
    theme_classic(base_size = 10) +
    labs(x = '',
         y = 'Effect size') +
    scale_x_discrete(labels = exp.labels) 
  effsizeP.p
  
  # N:P excretion
  # N excretion
  NPexcr.p <- plot_emm(pwcNP[["emmeans"]]) +
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
  ggsave('figures/final-figures/Fig1.tiff', 
         width = 7, height = 7, 
         units = 'in', dpi = 600, compression = 'lzw')
  
  # Figure 2 ----
  # Tag
  TAgexcr.p <- plot_boxplot(excr.Tag$massnorm.Tag.excr)
  TAgexcr.p
  # N:Ag
  NAgexcr.p <- plot_boxplot(excr.Tag$massnorm.NAg.excr)
  NAgexcr.p
  # P:Ag
  PAgexcr.p <- plot_boxplot(excr.Tag$massnorm.PAg.excr)
  PAgexcr.p
  
  # combine plots ----
  ggarrange(TAgexcr.p, NAgexcr.p, PAgexcr.p, ncol = 1, nrow = 3, 
            labels = c("(a)", "(b)", "(c)"), 
            font.label = list(size = 10), label.x = 0.2, label.y = 1,
            legend = 'none', common.legend = T, align = 'v')
  ggsave('figures/final-figures/Fig2.tiff', 
         width = 3.33, height = 7, 
         units = 'in', dpi = 300)
  
  # ..Figure 3 ----
  # Schiettekate et al model ----
  # N excretion 
  Nexcr_mS.p <- plot_mod(iter_yp$N.excretion.mod) +
    geom_point(aes(x = Mass, y = N.excretion.rate, color = Lake, shape = Year), 
               data = NPexcr %>%  filter(Year != '2014'),
               size = point.size, alpha = point.alpha + .2) +
    labs(x = "", y = "N excretion (μg N/ind/h)") +
    scale_shape_manual(labels = exp_red.labels,
                       values = c(16, 8, 17), na.translate = F) 
  Nexcr_mS.p
  
  # P excretion
  Pexcr_mS.p <- plot_mod(iter_yp$P.excretion.mod) +
    geom_point(aes(x = Mass, y = P.excretion.rate, color = Lake, shape = Year), 
               data = NPexcr,
               size = point.size, alpha = point.alpha + .2) +
    labs(x = "", y = "P excretion (μg P/ind/h)") 
  Pexcr_mS.p
  
  # V&M model ----
  # N excretion
  Nexcr_mVM.p <- ggplot(group_by(iter_yp_VM, iter), 
                        aes(x = Mass, y = N.excretion.m)) +
    geom_line() +
    geom_ribbon(aes(ymin = N.excretion.m - N.excretion.m.se,
                    ymax = N.excretion.m + N.excretion.m.se), alpha = .1) +
    theme_classic(base_size = 8) +
    geom_point(aes(x = Mass, y = N.excretion.rate, color = Lake, shape = Year), 
               data = NPexcr %>%  filter(Year != '2014'),
               size = point.size, alpha = point.alpha + .2) +
    labs(x = "Dry mass (g)", 
         y = "N excretion (μg N/ind/h)") +
      scale_colour_manual(name = 'Lake',
                          labels = lake.labels,
                          values = lake.colors) +
      scale_shape_manual(labels = exp_red.labels,
                         values = c(16, 8, 17), na.translate = F) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))
  Nexcr_mVM.p
  
  # P excretion
  Pexcr_mVM.p <- ggplot(group_by(iter_yp_VM, iter), 
                        aes(x = Mass, y = P.excretion.m)) +
    geom_line() +
    geom_ribbon(aes(ymin = P.excretion.m - P.excretion.m.se,
                    ymax = P.excretion.m + P.excretion.m.se), alpha = .1) +
    theme_classic(base_size = 8) +
    geom_point(aes(x = Mass, y = P.excretion.rate, color = Lake, shape = Year), 
               data = NPexcr %>%  filter(Year != '2014'),
               size = point.size, alpha = point.alpha + .2) +
    labs(x = "Dry mass (g)", 
         y = "P excretion (μg N/ind/h)") +
    scale_colour_manual(name = 'Lake',
                        labels = lake.labels,
                        values = lake.colors) +
    scale_shape_manual(labels = exp_red.labels,
                       values = c(16, 8, 17), na.translate = F) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))
  Pexcr_mVM.p
  
  # combine plots ----
  fig3.legend <- get_legend(Pexcr_mS.p)
  ggarrange(Nexcr_mS.p, Pexcr_mS.p, Nexcr_mVM.p, Pexcr_mVM.p, 
            nrow = 2, ncol = 2,
            legend = 'right', common.legend = F, align = 'hv', 
            labels = c('(a)', '(b)', '(c)', '(d)'), 
            label.x = 0.15, label.y = 1, font.label = list(size = 8), 
            legend.grob = fig3.legend)
  ggsave('figures/final-figures/Fig3.tiff', 
         width = 7, height =  4, 
         units = 'in', dpi = 600, compression = 'lzw')
  
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
  ggsave('figures/final-figures/Fig4.tiff', 
         width = 7, height = 4, 
         units = 'in', dpi = 600)
  
  # Figure S1
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
  
  ggsave('figures/final-figures/FigS1.tiff', 
         width = 7, height = 3.5, 
         units = 'in', dpi = 600, compression = 'lzw')
  
  # Figure S2
  Cing_mod.p <- 
    ggplot(group_by(iter_yp, iter), aes(x = dry.mass, y = C.ingestion.mod)) +
    stat_lineribbon(alpha = 0.8, show.legend = T) +
    scale_fill_brewer(name = 'Confidence Interval',
                      labels = c('95%', '80%', '50%')) +
    labs(x = "Dry mass (g)", y = "C ingestion rate (μg C/ind/h)") +
    theme_classic(base_size = 10) 
  Cing_mod.p
  
  ggsave('figures/final-figures/FigS3.tiff', 
         width = 7, height = 3.5, 
         units = 'in', dpi = 600)
  
  # Figure S3 ----
  lim_mod.p <- 
    ggplot(limitation, aes(x = dry.mass, y = prop_lim)) +
    geom_point(aes(x = dry.mass, y = prop_lim, color = nutrient), size = 1.5) +
    geom_line(aes(x = dry.mass, y = prop_lim, color = nutrient), size = 0.5) +
    labs(x = "Dry mass (g)", y = "Proportion of iterations", 
         color = "Limiting element") +
    theme_classic(base_size = 10) +
    scale_color_grey(start = 0.8, end = 0.1, labels = c("C", "N", "P")) 
  lim_mod.p
  
  ggsave('figures/final-figures/FigS2.tiff', 
         width = 7, height = 4, 
         units = 'in', dpi = 300, compression = 'lzw')
  