  # NAME of paper
  # this code was developped by S. Klemet-N'Guessan in 2020 and 2021

  # load libraries and read silver nanoparticles (NP) dataset ----
  
  library(tidyverse)
  library(car)
  library(RColorBrewer)
  library(emmeans) # post-hoc pairwise comparisons using lm, glm, levDev models
  library(effectsize) # other way to calculate effect size
  #library(effects) # effects from a fitted model
  library(ggpubr) # to add multiple graphs on the same page + stats
  library(fishflux) # to do fish stoich model
  library(parallel) # for mcmapply function for sensitivity analysis
  library(fishualize) # for colours based on fish colours
  library(tidybayes)
  library(patchwork) # to align multiple plots
  library(cowplot) # to align multiple plots
  library(datawizard) # to do summary statistics

  NPer <- read.csv('20 04 21 AgNP ELA lakes fish excretion.csv',
                   stringsAsFactors = F, na.strings = c("", "NA", "."), 
                   strip.white = TRUE, sep = ",")
  Pmd <- read.csv('21 05 17 AgNP fish excretion model_R.csv',
                    stringsAsFactors = F, na.strings = c("", "NA", "."), 
                    strip.white = TRUE, sep = ",")
  paramm <- read.csv('21 09 20 param_modelsp_FishStoich.csv',
                     stringsAsFactors = F, na.strings = c("", "NA", "."), 
                     strip.white = TRUE, sep = ",")
  
  str(paramm)
  
  str(NPer)
  head(NPer)
  
  # clean, rename, and add variables to the dataset ----
  NPexcr <- NPer %>% 
    rename(Year = Sampling.year,
           N.excretion = N.excretion.rate..Âµg.ind.h.,
           P.excretion = P.excretion.rate..Âµg.ind.h.,
           C.excretion = C.excretion..mg.C.ind.h.,
           Tag.excretion = Tag.excretion.rate..Âµg.ind.h.,
           Excreted.CN = Excreted.C.N..molar.,
           Excreted.NP = Excreted.N.P..molar.,
           Excreted.CP = Excreted.C.P..molar.,
           Mass = Dry.mass..g.,
           Temperature = Temperature.Â.,
           FishID = ï..Fish.ID) %>% 
    mutate(Log.mass = log(Mass),
           Log10.mass = log10(Mass),
           Log.N.excretion = log(N.excretion),
           Log.P.excretion = log(P.excretion),
           Log.Tag.excretion = log(Tag.excretion),
           Log10.N.excretion = log10(N.excretion),
           Log10.P.excretion = log10(P.excretion),
           Log10.C.excretion = log10(C.excretion),
           P.excretion.d = (P.excretion*24)/10^6, # P excretion rate from ug/h to g/d
           N.excretion.d = (N.excretion*24)/10^6,
           Total.length.cm = Total.length..mm./10,
           Lake = as.factor(Lake),
           Year = as.factor(Year)) %>% 
    filter(Mass!= c(3.75, 2.25))
  
  str(NPexcr)
  
  
  str(Pmd)
  Pmod <- Pmd %>% rename(Ingestion.rate = ï..Ingestion.rate,
                         ModP.excretion = Modelled.P.excretion.rate,
                         MeasP.excretion = Measured.P.excretion) %>% 
    mutate(FoodCP = as.factor(Food.C.P),
           Year = as.factor(Year),
           Lake = as.factor(Lake))
  str(Pmod)
  
  parammod <- paramm %>% rename(ac_m = ï..ac_m)
  str(parammod)
  

  ############################### TEST #######################################
  
  # ..log10 mass vs log10 N or P excretion for all the data by lake ----
  # visualize Log10 N excretion vs. Log10 mass by Lake
  windows(width = 14, height = 7)
  ggplot(NPexcr, aes(x = Log10.mass, y = Log10.N.excretion)) +
    geom_point(size = 5) +
    theme_classic(base_size = 26) +
    facet_wrap(~Lake) +
    theme(axis.text = element_text(face = 'bold'),
          axis.line = element_line(size = 1),
          panel.background = element_rect(colour = 'black'),
          panel.grid = element_blank(),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'),
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') 
  # very bad relationship, almost flat for P excretion +
  # not convinced it is actually linear for N excretion
  # probably due to the low variance in mass among individuals (from 0.4 to 8g) 
  # for the entire dataset
  
  # get coeff of variation log10 N excretion vs. log10 mass
  b.N222 <- lm(Log10.N.excretion ~ Log10.mass, 
     data = NPexcr %>% filter(Lake == '222'))$coefficients["Log10.mass"]
  # b coeff = -0.06
  b.N239 <- lm(Log10.N.excretion ~ Log10.mass, 
     data = NPexcr %>% filter(Lake == '239'))$coefficients["Log10.mass"]
  # b coeff = 0.58
  lm(Log10.N.excretion ~ Log10.mass, data = NPexcr) # b coeff = 0.40
  
  # visualize Log10 P excretion vs. Log10 mass by Lake
  ggplot(NPexcr, aes(x = Log10.mass, y = Log10.P.excretion)) +
    geom_point(size = 5) +
    theme_classic(base_size = 26) +
    facet_wrap(~Lake) +
    theme(axis.text = element_text(face = 'bold'),
          axis.line = element_line(size = 1),
          panel.background = element_rect(colour = 'black'),
          panel.grid = element_blank(),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'),
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') 
  
  # get coeff of variation log10 P excretion vs. log10 mass
  b.P222 <- lm(Log10.P.excretion ~ Log10.mass,
     data = NPexcr %>% filter(Lake == '222'))$coefficients["Log10.mass"]
  # b coeff = -0.06
  b.P239 <- lm(Log10.P.excretion ~ Log10.mass,
    data = NPexcr %>% filter(Lake == '239'))$coefficients["Log10.mass"] 
  # b coeff = 0.75
  lm(Log10.P.excretion ~ Log10.mass, data = NPexcr) # b coeff = 0.06
  
  # get coeff of variation log10 P excretion vs. log10 mass
  windows(width = 14, height = 7)
  ggplot(NPexcr, aes(x = Log10.mass, y = Log10.C.excretion)) +
    geom_point(size = 5) +
    theme_classic(base_size = 26) +
    facet_wrap(~Lake) +
    theme(axis.text = element_text(face = 'bold'),
          axis.line = element_line(size = 1),
          panel.background = element_rect(colour = 'black'),
          panel.grid = element_blank(),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'),
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right')
  
  lm(Log10.C.excretion ~ Log10.mass,
     data = NPexcr %>% filter(Lake == '222')) # b coeff = -3.54
  lm(Log10.C.excretion ~ Log10.mass,
     data = NPexcr %>% filter(Lake == '239')) # b coeff = -1.34
  
  # ..do mass-normalized N and P excretion ----
  # get coeffs of variation for allometric  relationship
  NPmasscorr.combined <- data.frame(Lake = c('222', '239', '222', '239'),
                                    b.coeff.N.excr = c(b.N222, b.N239),
                                    b.coeff.P.excr = c(b.P222, b.P239))
  
  # mass-correct the excretion WITHOUT coeff of variation 
  NPexcr <- NPexcr %>% 
    mutate(massnorm.N.excr = N.excretion/Mass,
           massnorm.P.excr = P.excretion/Mass,
           massnorm.C.excr = C.excretion/Mass,
           massnorm.Tag.excr = Tag.excretion/Mass,
           Log10.massnorm.Tag.excr = log10(massnorm.Tag.excr),
           Log10.massnorm.N.excr = log10(massnorm.N.excr),
           Log10.massnorm.P.excr = log10(massnorm.P.excr),
           massnorm.NP.excr = (massnorm.N.excr/massnorm.P.excr)/(14/31),
           Log10.massnorm.NP.excr = log10(massnorm.NP.excr),
           massnorm.CN.excr = (massnorm.C.excr/massnorm.N.excr)/(12/14),
           massnorm.CP.excr = (massnorm.C.excr/massnorm.P.excr)/(12/31),
           massnorm.NP.excr2 = Excreted.NP/Mass,
           massnorm.NAg.excr = massnorm.N.excr/massnorm.Tag.excr/(14/107),
           massnorm.PAg.excr = massnorm.P.excr/massnorm.Tag.excr/(12/107))
  
  
  # take the average of mass-corrected excretion
  NPexcr.ts <- NPexcr %>% group_by(Lake, Year) %>% 
    summarize(massnorm.N.excr.av = mean(massnorm.N.excr, na.rm = TRUE),
              massnorm.N.excr.sd = sd(massnorm.N.excr, na.rm = TRUE),
              massnorm.P.excr.av = mean(massnorm.P.excr, na.rm = TRUE),
              massnorm.P.excr.sd = sd(massnorm.P.excr, na.rm = TRUE),
              massnorm.C.excr.av = mean(massnorm.C.excr, na.rm = TRUE),
              massnorm.C.excr.sd = sd(massnorm.C.excr, na.rm = TRUE),
              massnorm.NP.excr.av = mean(massnorm.NP.excr, na.rm = TRUE),
              massnorm.NP.excr.sd = sd(massnorm.NP.excr, na.rm = TRUE),
              massnorm.CN.excr.av = mean(massnorm.CN.excr, na.rm = TRUE),
              massnorm.CN.excr.sd = sd(massnorm.CN.excr, na.rm = TRUE),
              massnorm.CP.excr.av = mean(massnorm.CP.excr, na.rm = TRUE),
              massnorm.CP.excr.sd = sd(massnorm.CP.excr, na.rm = TRUE),
              massnorm.Tag.excr.av = mean(massnorm.Tag.excr, na.rm = TRUE),
              massnorm.Tag.excr.sd = sd(massnorm.Tag.excr, na.rm = TRUE),
              massnorm.NAg.excr.av = mean(massnorm.NAg.excr, na.rm = TRUE),
              massnorm.NAg.excr.sd = sd(massnorm.NAg.excr, na.rm = TRUE),
              massnorm.PAg.excr.av = mean(massnorm.PAg.excr, na.rm = TRUE),
              massnorm.PAg.excr.sd = sd(massnorm.PAg.excr, na.rm = TRUE))
  
  # ..summary statistics info for paper methods & discussion ----
  
  NPexcr.ss1 <- NPexcr %>% group_by(Lake) %>%
    summarize(Mass.min = min(Mass, na.rm = TRUE),
              Mass.max = max(Mass, na.rm = TRUE),
              Mass.av = mean(Mass, na.rm = TRUE),
              Mass.sd = sd(Mass, na.rm = TRUE),
              massnorm.N.excr.min = min(massnorm.N.excr, na.rm = TRUE),
              massnorm.N.excr.max = max(massnorm.N.excr, na.rm = TRUE),
              massnorm.N.excr.av = mean(massnorm.N.excr, na.rm = TRUE),
              massnorm.N.excr.sd = sd(massnorm.N.excr, na.rm = TRUE),
              massnorm.P.excr.min = min(massnorm.P.excr, na.rm = TRUE),
              massnorm.P.excr.max = max(massnorm.P.excr, na.rm = TRUE),
              massnorm.P.excr.sd = sd(massnorm.P.excr, na.rm = TRUE))
  
  NPexcr.ss1 <- NPexcr %>% group_by(Lake, Year) %>%
    select(c('Mass', 'massnorm.N.excr', 'massnorm.P.excr')) %>% 
    describe_distribution() 
  
  NPexcr.ss2 <- NPexcr %>% group_by(Lake, Year) %>%
    summarize(Tag.excr.min = min(massnorm.Tag.excr, na.rm = T),
              Tag.excr.max = max(massnorm.Tag.excr, na.rm = T),
              Tag.excr.av = mean(massnorm.Tag.excr, na.rm = T),
              Tag.excr.sd = sd(massnorm.Tag.excr, na.rm = T),
              P.excretion.av = mean(P.excretion, na.rm = T),
              WMass.av = mean(Mass..g., na.rm = T),
              DMass.av = mean(Mass, na.rm = T),
              Temp.av = mean(Temperature, na.rm = T))
           
  ################################ ASSUMPTIONS CHECK ########################
  
  # ..normality and variance of the data in 2012, 2014, 2015 ----
  histkdnc(NPexcr$massnorm.N.excr) # normal distribution
  quantile(NPexcr$massnorm.N.excr, c(0,0.1,0.5), na.rm = T)
  plot((density(log(NPexcr$massnorm.N.excr),na.rm= T)))
  histkdnc(NPexcr$P.excretion) # normal distribution
  shapiro.test(NPexcr$Log.massnorm.N.excr)
  boxplot(N.excretion ~ Lake, 
          data = NPexcr)
  leveneTest(N.excretion ~ Year, 
             data = NPexcr %>% filter(Year != '2014'))
  # plot mean P/N excretion vs. SD to see if there is a positive relationship between the two. If yes, a log-transformation is possible. Here both mean P and N excretions have a good positive relationship with their respective SD
  ggplot(NPexcr.ts, aes(massnorm.P.excr.sd, massnorm.P.excr.av)) +
    geom_point()
  ggplot(NPexcr.ts, aes(massnorm.N.excr.sd, massnorm.N.excr.av)) +
    geom_point()
  # test for heterogeneity of variances using Levene's test BUT old-school, not really telling where the heterogeneity is
  leveneTest(Log.massnorm.N.excr ~ Lake*Year, data = NPexcr)
  # variances significantly diff
  leveneTest(Log.massnorm.P.excr ~ Lake*Year, data = NPexcr) 
  # variances significantly diff
  
  # ..normality of the residuals 
  windows(width = 14, height = 7)
  qqPlot(lm(Log10.massnorm.NP.excr ~ Lake*Year, 
             data = NPexcr %>% filter(Year != '2014'))) 
  skewness(NPexcr$N.excretion, na.rm = T) 
   
  # skewness = 2.9 > 0.5 threshold (Webster and Oliver, 2007)
  qqPlot(lm(massnorm.P.excr ~ Lake*Year, data = NPexcr)) 
  skewness(NPexcr$P.excretion, na.rm = T)
  # skewness = 1.5 > 0.5 threshold (Webster and Oliver, 2007)
  qqPlot(lm(massnorm.C.excr ~ Lake, 
            data = NPexcr))
  qqPlot(lm(Tag.excretion ~ Year, 
            data = NPexcr))
  
  # ..plot standardized residuals of the lm against predicted values ----
  # and make informed decision by comparing it to plotted residuals 
  # of randomly generated datasets
  # Tag or N excretion
  z = lm(massnorm.Tag.excr ~ Year, 
          data = NPexcr$residual)
  n = length(z)
  xz = cbind(matrix(rnorm(12*n), nrow = n),z,matrix(rnorm(12*n), nrow = n))
  colnames(xz) = c(letters[1:12],"Z",letters[13:24])
  
  opar = par()
  par(mfrow=c(5,5));
  par(mar=c(0.5,0.5,0.5,0.5))
  par(oma=c(1,1,1,1));
  
  ytpos = (apply(xz,2,min)+3*apply(xz,2,max))/4
  cn = colnames(xz)
  
  for(i in 1:25) {
    qqnorm(xz[,i],axes=FALSE,ylab= colnames(xz)[i],xlab="",main="")
    qqline(xz[,i],col=2,lty=2)
    box("figure", col="darkgreen")
    text(-1.5,ytpos[i],cn[i])
  }
  par(opar)
  
  # P excretion
  z = glm(massnorm.P.excr ~ Lake*Year, 
          data = NPexcr, 
          family = Gamma(link = "log"))$residual
  n = length(z)
  xz = cbind(matrix(rnorm(12*n), nrow = n),z,matrix(rnorm(12*n), nrow = n))
  colnames(xz) = c(letters[1:12],"Z",letters[13:24])
  
  opar = par()
  par(mfrow=c(5,5));
  par(mar=c(0.5,0.5,0.5,0.5))
  par(oma=c(1,1,1,1));
  
  ytpos = (apply(xz,2,min)+3*apply(xz,2,max))/4
  cn = colnames(xz)
  
  for(i in 1:25) {
    qqnorm(xz[,i],axes=FALSE,ylab= colnames(xz)[i],xlab="",main="")
    qqline(xz[,i],col=2,lty=2)
    box("figure", col="darkgreen")
    text(-1.5,ytpos[i],cn[i])
  }
  par(opar)
  
  # whether using N.excretion, log or log10 transformed N excretion, 
  # residuals vs. predicted values still look quite different from 
  # those generated randomly: light/heavy-tailed on both ends (N excretion)
  # or on the right only (log and log10 transformation)
  
  
  ############################ DATA ANALYSIS & FIGURES ############################
  
  # ..fish mass ----
  
  lm.mass <- lm(Mass ~ Lake*Year,
                  data = NPexcr)
  summary(lm.mass)
  Anova(lm.mass)
  
  lm.mass_emmeans <- emmeans(lm.mass, ~ Lake:Year)
  lm.mass_emmeans
  pairs(lm.mass_emmeans)
  pairwise.emmc(lm.mass_emmeans)
  
  # Figure S1 ----
  # .....fish mass across years and lakes ----
  windows(width = 14, height = 7)
  mexcrp <- ggplot(NPexcr, 
         aes(x = Year, y = Mass, color = Lake)) +
    geom_jitter(size = 2, aes(x = Year, y = Mass), 
                position = position_jitterdodge(jitter.width = 0.2),
                show.legend = FALSE, alpha = 0.5) +
    theme_classic(base_size = 22) +
    # stat_compare_means(method = 'anova') +
    geom_boxplot(outlier.shape = NA, size = 1, fill = NA) +
    scale_x_discrete(labels = c('Pre-addition',  'Year 1', "Year 2")) +
    labs(x = '',
         y = 'Dry mass (g)') +
    theme_classic(base_size = 22) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.title.x = element_blank(),
          # axis.text.x = element_blank(),
          # axis.title.y = element_text(vjust = +1),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') +
    scale_colour_manual(name = 'Lake',
                        labels = c('AgNPs L222', 'Reference L239'),
                        values = c("black","gray60")) +
    annotate("text", x = 0.81, y = 1.6, label = 'a', 
             size = 7, fontface = 'bold') +
    annotate("text", x = 1.19, y = 0.36, label = 'b', 
             size = 7, fontface = 'bold', color = 'gray60') +
    annotate("text", x = 1.81, y = 1.40, label = 'a', 
             size = 7, fontface = 'bold') +
    annotate("text", x = 2.2, y = 0.5, label = 'b', 
             size = 7, fontface = 'bold', color = 'gray60') +
    annotate("text", x = 2.81, y = 0.45, label = 'c', 
             size = 7, fontface = 'bold') +
    annotate("text", x = 3.18, y = 0.30, label = 'bc', 
             size = 7, fontface = 'bold', color = 'gray60')
  mexcrp
  
  ggsave('final figures/Fig S1.png', 
         width = 14, height = 7, 
         units = 'in', dpi = 300)
  
  # Figure 1 ----
  # ..N excretion anova + figure ----
  
  # testing linear model with log-transformed data on raw column 
  # so that can revert back to response when doing emmeans
  lm.Nxnorm <- lm(log(massnorm.N.excr) ~ Lake*Year, 
     data = NPexcr %>% filter(Year != '2014'))
  anova(lm.Nxnorm)
  # which levels are different? Lakes from 2012 vs. 2015
  TukeyHSD(aov(lm.Nxnorm), ordered = F)
  
  # pairwise comparisons using emmeans
  lm.Nxnorm_emmeans <- emmeans(lm.Nxnorm, ~ Lake|Year, type = 'response')
  pairs(lm.Nxnorm_emmeans)
  emmip(lm.Nxnorm, Lake ~ Year)
  lm.Nxnorm_emmeans
  plot(lm.Nxnorm_emmeans, comparison = T)
  
  # extracting effects from emmeans
  lm.Nxnorm_emmeansdf <- as.data.frame(summary(lm.Nxnorm_emmeans))
  
  # calculating effect sizes + extracting effects
  eff_sizeNx <- eff_size(lm.Nxnorm_emmeans, sigma = sigma(lm.Nxnorm), 
                         edf = df.residual((lm.Nxnorm)))
  eff_sizeNx
  eff_sizeNxdf <- as.data.frame(summary(eff_sizeNx))
  eff_sizeNxdf2 <- eff_sizeNxdf %>% mutate(effect.size =
                                             ifelse(effect.size < 0,
                                                    effect.size*(-1), 0.8784323),
                                          lower.CL =
                                            ifelse(lower.CL < -0.8,
                                                   -0.4318393, -0.03054899),
                                          upper.CL =
                                            ifelse(upper.CL < 1,
                                                   0.97387891, 1.7874136))

  # ...data visualization 
  # .....time series 2012 to 2015 using model ----
  windows(width = 14, height = 7)
  Nexcr.p <- ggplot(lm.Nxnorm_emmeansdf, 
         aes(x = Year, y = response, color = Lake,
             group = Lake)) +
    geom_point(size = 4, position = position_dodge(0.5)) +
    geom_errorbar(aes(ymax = lower.CL, ymin = upper.CL), 
                  width = 0.2, lwd = 1.5, position = position_dodge(0.5)) +
    geom_jitter(data = NPexcr %>% filter(Year != '2014'),
                aes(x = Year, y = massnorm.N.excr),
                size = 1.5, position = position_jitterdodge(jitter.width = 0.2),
                show.legend = FALSE, alpha = 0.5) +
    scale_x_discrete(labels = c('Pre-addition',  'Year 2')) +
    theme_classic(base_size = 26) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.title.x = element_blank(),
          # axis.text.x = element_blank(),
          # axis.title.y = element_text(vjust = +1),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'top') +
    scale_colour_manual(name = 'Lake',
                      labels = c('AgNPs L222', 'Reference L239'),
                      values = c("black","gray60")) +
    scale_y_continuous(name = 'Mass-corrected \n N excretion (μg N/g/h)') +
    # annotate("text", x = 0.5, y = 9500, label = 'A)', 
    #          size = 9, fontface = 'bold') +
    annotate("text", x = 0.87, y = 995, label = 'a', 
             size = 7, fontface = 'bold') +
    annotate("text", x = 1.12, y = 1095, label = 'a', 
             size = 7, fontface = 'bold', color = 'gray60') +
    annotate("text", x = 1.87, y = 4195, label = 'b', 
             size = 7, fontface = 'bold') +
    annotate("text", x = 2.12, y = 2695, label = 'b', 
             size = 7, fontface = 'bold', color = 'gray60')
  Nexcr.p
  
  # .....effect sizes from 2012 2015 using model ----
  windows(width = 14, height = 7)
  eff_sizeN.p <- ggplot(eff_sizeNxdf2,
                        aes(x = Year, y = exp(effect.size))) +
    geom_point(size = 4) +
    geom_pointrange(aes(ymax = exp(lower.CL), ymin = exp(upper.CL)),lwd = 1.2) +
    #coord_flip() +
    theme_classic(base_size = 26) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face = "bold"),
          axis.title.x = element_blank(),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') +
    scale_y_continuous(name = 'Effect size') +
    scale_x_discrete(labels = c('Pre-addition', 'Year 1', 'Year 2')) 
  eff_sizeN.p
  
  # ..P excretion anova + figure ----
  # testing linear model with log-transformed data on raw column 
  # so that can revert back to response when doing emmeans
  lm.Pxnorm <- lm(log(massnorm.P.excr) ~ Lake*Year, 
                  data = NPexcr)
  Anova(lm.Pxnorm)
  # which levels are different? Lakes from 2014 vs. 2015
  TukeyHSD(aov(lm.Pxnorm), ordered = F)
  
  # emmeans
  lm.Pxnorm_emmeans <- emmeans(lm.Pxnorm, ~ Lake|Year, type = 'response')
  pairs(lm.Pxnorm_emmeans)
  emmip(lm.Pxnorm, Lake ~ Year)
  lm.Pxnorm_emmeans
   plot(lm.Pxnorm_emmeans, comparison = T)
  
  # extracting effects
  lm.Pxnorm_emmeansdf <- as.data.frame(summary(lm.Pxnorm_emmeans))
  
  # calculating effect sizes
  eff_sizePx <- eff_size(lm.Pxnorm_emmeans, sigma = sigma(lm.Pxnorm), 
           edf = df.residual((lm.Pxnorm)))
  eff_sizePx
  eff_sizePxdf <- as.data.frame(summary(eff_sizePx))
  
  # .....time series 2012 to 2015 using model ----
  windows(width = 14, height = 7)
  Pexcr.p <- ggplot(lm.Pxnorm_emmeansdf, 
         aes(x = Year, y = response, color = Lake,
             group = Lake)) +
    geom_point(size = 4, position = position_dodge(0.5)) +
    geom_errorbar(aes(ymax = lower.CL, ymin = upper.CL),
                  width = 0.2, lwd = 1.5, position = position_dodge(0.5)) +
    geom_jitter(data = NPexcr,
                aes(x = Year, y = massnorm.P.excr),
                size = 1.5, position = position_jitterdodge(jitter.width = 0.2),
                show.legend = FALSE, alpha = 0.5) +
    scale_x_discrete(labels = c('Pre-addition', 'Year 1', 'Year 2')) +
    theme_classic(base_size = 26) +
    theme(text = element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.title.x = element_blank(),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'top') +
    scale_y_continuous(name = 'Mass-corrected \n P excretion (μg P/g/h)') +
    scale_colour_manual(name = 'Lake',
                        labels = c('AgNPs L222', 'Reference L239'),
                        values = c("black","gray60")) +
    annotate("text", x = 0.87, y = 39, label = 'a', 
             size = 7, fontface = 'bold') +
    annotate("text", x = 1.12, y = 147, label = 'b', 
             size = 7, fontface = 'bold', color = 'gray60') +
    annotate("text", x = 1.87, y = 36, label = 'a', 
             size = 7, fontface = 'bold') +
    annotate("text", x = 2.12, y = 100, label = 'b', 
             size = 7, fontface = 'bold', color = 'gray60') +
    annotate("text", x = 2.87, y = 99, label = 'c', 
             size = 7, fontface = 'bold') +
    annotate("text", x = 3.12, y = 123, label = 'bc', 
             size = 7, fontface = 'bold', color = 'gray60')
  Pexcr.p
  
  # .....effect sizes from 2012-2015 using model ----
  windows(width = 14, height = 7)
  eff_sizeP.p <- ggplot(eff_sizePxdf,
         aes(x = Year, y = exp(-effect.size))) +
    geom_point(size = 4) +
    geom_pointrange(aes(ymax = exp(-lower.CL), ymin = exp(-upper.CL)),lwd = 1.2) +
    #coord_flip() +
    theme_classic(base_size = 26) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face = "bold"),
          axis.title.x = element_blank(),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') +
    scale_y_continuous(name = 'Effect size') +
    scale_x_discrete(labels = c('Pre-addition', 'Year 1', 'Year 2')) 
  eff_sizeP.p
  
  #..N:P excretion (molar) ----
  lm.NPxnorm <- lm(log(massnorm.NP.excr) ~ Lake*Year, 
                   data = NPexcr %>% filter(Year != '2014'))
  summary(lm.NPxnorm)
  Anova(lm.NPxnorm)
  # which levels are different? Lakes from 2014 vs. 2015
  TukeyHSD(aov(lm.NPxnorm), ordered = F)
  
  # emmeans
  lm.NPxnorm_emmeans <- emmeans(lm.NPxnorm, ~ Lake|Year, type = 'response')
  lm.NPxnorm_emmeans
  pairs(lm.NPxnorm_emmeans)
  
  # extracting effects
  lm.NPxnorm_marginal_means <- as.data.frame(summary(lm.NPxnorm_emmeans))
  
  # effect size
  eff_sizeNPx <- eff_size(lm.NPxnorm_emmeans, sigma = sigma(lm.NPxnorm), 
                          edf = df.residual((lm.NPxnorm)))
  eff_sizeNPx
  eff_sizeNPxdf <- as.data.frame(summary(eff_sizeNPx))
  
  # ...data visualization ----
  # .....time series 2012 to 2015 using model ----
  windows(width = 14, height = 7)
  NPexcr.p <- ggplot(lm.NPxnorm_marginal_means, 
                    aes(x = Year, y = response, color = Lake,
                        group = Lake)) +
    geom_point(size = 4, position = position_dodge(0.5)) +
    geom_errorbar(aes(ymax = lower.CL, ymin = upper.CL),
                  width = 0.2, lwd = 1.5, position = position_dodge(0.5)) +
    geom_jitter(data = NPexcr %>% filter(Year != '2014'),
                aes(x = Year, y = massnorm.NP.excr),
                size = 1.5, position = position_jitterdodge(jitter.width = 0.2),
                show.legend = FALSE, alpha = 0.5) +
    theme_classic(base_size = 26) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.title.x = element_blank(),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') +
    scale_fill_brewer(palette = "Accent") +
    scale_y_continuous(name = 'Mass-corrected \n N:P excretion (molar)') +
    scale_x_discrete(labels = c('Pre-addition', 'Year 2')) +
    scale_colour_manual(name = 'Lake',
                        labels = c('AgNPs L222', 'Reference L239'),
                        values = c("black","gray60")) +
    # annotate("text", x = 0.5, y = 450, label = 'E)', 
    #          size = 9, fontface = 'bold') +
    annotate("text", x = 0.87, y = 80, label = 'ac', 
             size = 7, fontface = 'bold') +
    annotate("text", x = 1.12, y = 40, label = 'b', 
             size = 7, fontface = 'bold', color = 'gray60') +
    annotate("text", x = 1.87, y = 190, label = 'c', 
             size = 7, fontface = 'bold') +
    annotate("text", x = 2.12, y = 115, label = 'c', 
             size = 7, fontface = 'bold', color = 'gray60') 
  NPexcr.p
  
  # .....effect sizes from 2012 2015 using model ----
  windows(width = 14, height = 7)
  eff_sizeNP.p <- ggplot(eff_sizeNPxdf,
                         aes(x = Year, y = exp(effect.size))) +
    geom_point(size = 4) +
    geom_pointrange(aes(ymax = exp(lower.CL), ymin = exp(upper.CL)),lwd = 1.2) +
    #coord_flip() +
    theme_classic(base_size = 26) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face = "bold"),
          axis.title.x = element_blank(),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') +
    scale_y_continuous(name = 'Effect size') +
    scale_x_discrete(labels = c('Pre-addition', 'Year 2')) 
  eff_sizeNP.p
  
  # combine all graphs into figure 1 ----
  windows(width = 14, height = 7)
  ggarrange(Nexcr.p, eff_sizeN.p, Pexcr.p, eff_sizeP.p, NPexcr.p, 
            eff_sizeNP.p, ncol = 2, nrow = 3, 
            labels = c("(a)", "(b)", "(c)", "(d)","(e)", "(f)"),
            font.label = list(size = 26), label.x = 0.2, label.y = 1,
            legend = 'top', common.legend = T, align = 'v')
  ggsave('final figures/Fig 1.png', 
         width = 18, height = 18, 
         units = 'in', dpi = 300)

  # Figure 2 ----
  # ..Tag excretion t-test + figure ----
  # Ag
  hist(NPexcr$massnorm.Tag.excr)
  gghistogram(NPexcr, x = 'massnorm.Tag.excr', rug = T, 
              fill = 'Year', add = 'mean')
  var.test(massnorm.Tag.excr ~ Year, NPexcr, alternative = 'two.sided')
  t.test(Log10.massnorm.Tag.excr ~ Year, NPexcr, var.equal = F)
  w.Tag <- wilcox.test(massnorm.Tag.excr ~ Year, NPexcr)
  w.Tag
  
  # N:Ag
  hist(NPexcr$massnorm.NAg.excr)
  var.test(massnorm.NAg.excr ~ Year, NPexcr, alternative = 'two.sided')
  w.NAg <- wilcox.test(massnorm.NAg.excr ~ Year, NPexcr)
  w.NAg
  
  #P:Ag
  hist(NPexcr$massnorm.PAg.excr)
  var.test(massnorm.PAg.excr ~ Year, NPexcr, alternative = 'two.sided')
  w.PAg <- wilcox.test(massnorm.PAg.excr ~ Year, NPexcr)
  w.PAg
  
  # ...TAg 
  windows(width = 14, height = 7)
  TAgexcrp <- ggplot(NPexcr %>% filter(Year != 2012), 
         aes(x = Year, y = massnorm.Tag.excr)) +
    geom_boxplot(outlier.shape = NA, size = 1) +
    geom_jitter(size = 4, position = position_jitter(0.2),alpha = 0.5) +
    stat_compare_means(label.x = 0.7, label = 'p.format', size = 5) +
    theme_classic(base_size = 18) +
    theme(text = element_text(family = "Tahoma"),
          axis.title = element_text(face = "bold"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') +
    scale_y_continuous(name = 'Mass-corrected \n Ag excretion (μg Ag/g/h)') +
    scale_x_discrete(labels = c('Year 1', 'Year 2')) 
  TAgexcrp
  
  # N:Ag
  windows(width = 14, height = 7)
  NAgexcrp <- ggplot(NPexcr %>% filter(Year != 2012), 
                     aes(x = Year, y = massnorm.NAg.excr)) +
    geom_boxplot(outlier.shape = NA, size = 1) +
    geom_jitter(size = 4, position = position_jitter(0.2),alpha = 0.5) +
    stat_compare_means(label.x = 0.7, label = 'p.format', size = 5) +
    theme_classic(base_size = 18) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face = "bold"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') +
    scale_y_continuous(name = 'Mass-corrected \n N:Ag excretion (molar)') +
    scale_x_discrete(labels = c('Year 1', 'Year 2'))
  NAgexcrp
  
  # P:Ag
  PAgexcrp <- ggplot(NPexcr %>% filter(Year != 2012), 
                     aes(x = Year, y = massnorm.PAg.excr)) +
    geom_boxplot(outlier.shape = NA, size = 1) +
    geom_jitter(size = 4, position = position_jitter(0.2),alpha = 0.5) +
    stat_compare_means(label.x = 0.7, size = 5, label = 'p.format') +
    theme_classic(base_size = 18) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face = "bold"),
          axis.title.x = element_blank(),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'),
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') +
    scale_y_continuous(name = 'Mass-corrected \n P:Ag excretion (molar)') +
    scale_x_discrete(labels = c('Year 1', 'Year 2'))
  PAgexcrp
  
  # combine all graphs into figure 2 ----
  ggarrange(TAgexcrp, NAgexcrp, PAgexcrp, ncol = 1, nrow = 3, 
            labels = c("(a)", "(b)", "(c)"), 
            font.label = list(size = 18), label.x = 0.11, label.y = 1,
            legend = 'none', common.legend = T, align = 'v')
  ggsave('final figures/Fig 2.png', 
         width = 12, height = 12, 
         units = 'in', dpi = 300)
  
  
  # FishStoichModel ----
  parameters <- model_parameters("Perca flavescens", otolith = F,
                                           family = "Percidae", temp = 18, mirror = "se")
  # find length-weight relationship and caudal fin aspect ratio
  find_lw("Perca flavescens", mirror = "se")
  aspect_ratio("Perca flavescens")
  # find growth parameters
  growth_params <- growth_params("Perca flavescens", otolith = FALSE)
  # find trophic level
  trophic_level("Perca flavescens")
  
  # ..Run model ----
  # Yellow perch
  total.length <- exp((log(NPexcr$Mass..g./parammod$lwa_m)/parammod$lwb_m))
  mass <- unique(NPexcr$Mass..g.)
  mass <- NPexcr$Mass..g.
  FStoichm <- cnp_model_mcmc(TL = total.length, param = parammod)
  output <- fishflux::extract(FStoichm, c("Fn","Fp", "Ic", "Gp", "lim", 
                                          "Sc", "Sn", "Sp"))
  
  # get mass instead of total length
  get_iter <- function(x){
    get <- t(plyr::ldply(x))
    colnames(get) <- get[1,]
    get <- data.frame(apply(get[-1,],2,as.numeric))
    get$iter <- 1:nrow(get)
    return(get)
  }
  iter_yp <- (lapply(FStoichm$stanfit, FUN = function(x){
    rstan::extract(x, c("Fn", "Fp", "w1", "Ic", "IN"))})) %>%
    lapply(FUN = get_iter) %>%
    dplyr::bind_rows()
  iter_yp <- iter_yp %>% mutate(tl = rep(total.length, each = 500),
                                mass = rep(NPexcr$Mass..g., each = 500),
                                dry.mass = mass*0.25,
                                N.excretion.mod = Fn*10^6/24,
                                P.excretion.mod = Fp*10^6/24,
                                C.ingestion.mod = Ic*10^6/24)
  
  # ..Figure 3 ----
  # plot N excretion
  Nexcr_mod.p <- 
    ggplot(group_by(iter_yp, iter), aes(x = dry.mass, y = N.excretion.mod)) +
    stat_lineribbon(alpha = 0.8, show.legend = F) +
    scale_fill_brewer() +
    theme_classic(base_size = 22) +
    geom_point(aes(x = Mass, y = N.excretion, color = Lake, shape = Year), 
               data = NPexcr %>%  filter(Year != '2014'),
               size = 4) +
    labs(x = "Dry mass (g)", y = "Mass-uncorrected \n N excretion (μg N/ind/h)") +
    scale_colour_manual(name = 'Lake',
                        labels = c('AgNPs 222', 'Reference 239'),
                        values = c("black","gray60")) +
    scale_shape_manual(labels = c('Pre-addition', 'Year 2'),
                       values = c(16, 15, 17), na.translate = F) +
    theme(text = element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'none') +
    annotate("text", x = 0.15, y = 2300, label = '(a)',
             size = 8, fontface = 'bold')
  Nexcr_mod.p
  
  # plot P excretion
  Pexcr_mod.p <- 
    ggplot(group_by(iter_yp, iter), aes(x = dry.mass, y = P.excretion.mod)) +
    stat_lineribbon(alpha = 0.8, show.legend = F) +
    scale_fill_brewer() +
    theme_classic(base_size = 22) +
    geom_point(aes(x = Mass, y = P.excretion, color = Lake, shape = Year), 
               data = NPexcr, size = 4) +
    labs(x = "Dry mass (g)", 
         y = "Mass-uncorrected \n P excretion (μg P/ind/h)") +
    scale_colour_manual(name = 'Lake',
                        labels = c('AgNPs 222', 'Reference 239'),
                        values = c("black","gray60")) +
    scale_shape_manual(labels = c('Pre-addition', 'Year 1', 'Year 2'),
                       values = c(16, 17, 15), na.translate = F) +
    theme(text = element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') +
    annotate("text", x = 0.15, y = 56, label = '(b)',
             size = 8, fontface = 'bold')
  Pexcr_mod.p
  
  # combine all graphs into figure 3 ----
  # get legend function
  get_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  
  fig3.legend <- get_legend(Pexcr_mod.p)
  
  windows(width = 14, height = 7)
  ggarrange(Nexcr_mod.p, Pexcr_mod.p, nrow = 2,
            legend = 'right', common.legend = F, align = 'v', 
            legend.grob = fig3.legend)
  ggsave('final figures/Fig 3.png', 
         width = 14, height =  12, 
         units = 'in', dpi = 300)
  
  # ..Figure 4 ----
  # Predict P excretion increase with diet increase
  # only averages for this simulation
  params <- parammod[grep("_m" ,names(parammod))]  
  
  mod <- cnp_model_mcmc(4:5, params)
  
  ex <- fishflux::extract(mod, c("Ic", "In", "Ip", "Fn", "Fp", 
                                 "Sc", "Sn", "Sp",
                                 "st_cn", "st_np")) %>% dplyr::filter(TL < 4.5)
  
  # TERs
  Scn <- ex$st_cn_median
  Snp <- ex$st_np_median
  
  # fix Dc at 45
  params$Dc_m <- 45
  Dcm <- params$Dc_m
  Dnm <- Dcm/Scn
  Dpm <- Dnm/Snp
  
  
  DN <- seq(from = 4.7, to = 6.7, by = 0.1)
  DP <- DN/Snp
  DNP <- expand.grid(DN, DP)
  colnames(DNP) <- c("Dn", "Dp")
  
  # Run model for multiple Dn and Dp values
  sim <- 
    lapply(1:nrow(DNP), function(i){
      params$Dn_m <- DNP$Dn[i]
      params$Dp_m <- DNP$Dp[i]
      mod <- fishflux::cnp_model_mcmc(TL = 4:5, param = params)
      ex <- fishflux::extract(mod, c("Ic", "In", "Ip", "Fn", "Fp", "Sc", "Sn", "Sp", 
                                     "lim", "st_cn", "st_np")) %>% 
        dplyr::filter(TL < 4.5)
      result <- cbind(DNP[i,], ex)
      print(i/nrow(DNP))
      return(result)
    })
  
  simd <- plyr::ldply(sim)
  simd <- simd %>% mutate(N.excretion.mod = Fn_median*10^6/24,
                          P.excretion.mod = Fp_median*10^6/24,
                          C.ingestion.mod = Ic_median*10^6/24)
  
  lim_ref.p <- 
    ggplot()+
    geom_segment(aes(x = 4.7, y = 4.7/Snp, xend = Dnm, yend = Dpm), size = 2) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = Dnm, yend = 6.7/Snp), size = 2) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = 6.7, yend = Dpm), size = 2) +
    geom_text(aes(x = 5.25, y = 1.25, label = "N"), size = 25) +
    geom_text(aes(x = 6.25, y = 1.3, label = "C"), size = 25) +
    geom_text(aes(x = 6, y = 1.1, label = "P"), size = 25) +
    labs(x = "Diet N (%)", y = "Diet P (%)") +
    theme_bw(base_size = 26) +
    theme(aspect.ratio = 1, 
          title = element_text(size = 18), 
          text = element_text(family = "Tahoma", size = 16),
          axis.title = element_text(face="bold"),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines')) 
  lim_ref.p
  
  Cing_ref.p <-
    ggplot(simd) +
    geom_raster(aes(x = Dn, y = Dp, fill = C.ingestion.mod)) +
    scale_fill_fish(option = "Trimma_lantana", trans = "sqrt") +
    geom_segment(aes(x = 4.7, y = 4.7/Snp, xend = Dnm, yend = Dpm), size = 2) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = Dnm, yend = 6.7/Snp), size = 2) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = 6.7, yend = Dpm), size = 2) +
    theme_bw(base_size = 26) +
    theme(aspect.ratio = 1, 
          title = element_text(size = 18), 
          text = element_text(family = "Tahoma", size = 16),
          axis.title = element_text(face="bold"),
          legend.title = element_text(face = 'bold', hjust = 0),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines')) +
    labs(x = "Diet N (%)", y = "Diet P (%)", 
         fill = "C ingestion \n (μg C/ind/h)")
  Cing_ref.p
  
  Nexcr_ref.p <-
    ggplot(simd) +
    geom_raster(aes(x = Dn, y = Dp, fill = N.excretion.mod)) +
    geom_segment(aes(x = 4.7, y = 4.7/Snp, xend = Dnm, yend = Dpm), size = 2) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = Dnm, yend = 6.7/Snp), size = 2) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = 6.7, yend = Dpm), size = 2) +
    scale_fill_fish(option = "Trimma_lantana", trans = "sqrt") +
    theme_bw(base_size = 26) +
    theme(aspect.ratio = 1, 
          title = element_text(size = 18), 
          text = element_text(family = "Tahoma", size = 16),
          axis.title = element_text(face="bold"),
          legend.title = element_text(face = 'bold', hjust = 0),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'))  +
    labs(x = "Diet N (%)", y = "Diet P (%)", 
         fill = 'N excretion \n (μg N/ind/h)')
  Nexcr_ref.p
  
  Pexcr_ref.p <-
    ggplot(simd) +
    geom_raster(aes(x = Dn, y = Dp, fill = P.excretion.mod)) +
    scale_fill_fish(option = "Trimma_lantana", trans = "sqrt") +
    geom_segment(aes(x = 4.7, y = 4.7/Snp, xend = Dnm, yend = Dpm), size = 2) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = Dnm, yend = 6.7/Snp), size = 2) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = 6.7, yend = Dpm), size = 2) +
    theme_bw(base_size = 26) +
    theme(aspect.ratio = 1, 
          title = element_text(size = 18), 
          text = element_text(family = "Tahoma", size = 16),
          axis.title = element_text(face="bold"),
          legend.title = element_text(face = 'bold', hjust = 0),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines')) +
    labs(x = "Diet N (%)", y = "Diet P (%)", 
         fill = 'P excretion \n (μg P/ind/h)') 
  Pexcr_ref.p
  
  # combine graphs
  plot_grid(lim_ref.p, Cing_ref.p, Nexcr_ref.p, Pexcr_ref.p, 
            ncol = 2, nrow = 2, align = 'hv', axis = "r", labels = c("(a)", "(b)", "(c)", "(d)"),
            label_size = 18)
  ggsave('final figures/Fig 4.png', 
         width = 12, height = 7, 
         units = 'in', dpi = 300)
  
  # ..Figure S2 & S3: Ingestion and N/P/C limitation ----
  # plot C ingestion
  Cing_mod.p <- 
    ggplot(group_by(iter_yp, iter), aes(x = dry.mass, y = C.ingestion.mod)) +
    stat_lineribbon(alpha = 0.8, show.legend = T) +
    scale_fill_brewer(name = 'Confidence Interval',
                      labels = c('95%', '80%', '50%')) +
    labs(x = "Dry mass (g)", y = "C ingestion rate (μg C/ind/h)") +
    theme_classic(base_size = 22) + 
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face = "bold"),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right')  
  Cing_mod.p
  
  ggsave('final figures/Fig S3.png', 
         width = 14, height = 7, 
         units = 'in', dpi = 300)
  
  # plot N/P/C limitation
  windows(width = 14, height = 7)
  limitation <- limitation(FStoichm, plot = F)
  a <- parammod$lwa_m
  b <- parammod$lwb_m
  limitation <- limitation %>% mutate(mass = a*tl^b,
                                      dry.mass = mass*0.25)
  
  lim_mod.p <- 
    ggplot(limitation, aes(x = dry.mass, y = prop_lim)) +
    geom_point(aes(x = dry.mass, y = prop_lim, color = nutrient), size = 2) +
    geom_line(aes(x = dry.mass, y = prop_lim, color = nutrient), size = 1) +
    labs(x = "Dry mass (g)", y = "Proportion of iterations", color = "Limiting element") +
    theme_classic(base_size = 22) +
    scale_color_grey(start = 0.8, end = 0.1, labels = c("C", "N", "P")) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face="bold"),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right')  
  lim_mod.p
  
  ggsave('final figures/Fig S2.png', 
         width = 14, height = 7, 
         units = 'in', dpi = 300)
  
  
  
  sessionInfo() 
  #################################### END OF CODE ##################################
  
             