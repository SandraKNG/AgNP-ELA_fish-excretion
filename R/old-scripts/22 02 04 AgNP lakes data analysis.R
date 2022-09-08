  # NAME of paper
  # this code was developped by S. Klemet-N'Guessan in 2020 and 2021

  # load libraries and read silver nanoparticles (NP) dataset ----
  
  library(tidyverse)
  library(car)
  library(RColorBrewer)
  library(emmeans) # post-hoc pairwise comparisons using lm, glm, levDev models
  #library(effects) # effects from a fitted model
  library(ggpubr) #to add multiple graphs on the same page
  library(fishflux) # to do fish stoich model
  library(parallel) # for mcmapply function for sensitivity analysis
  library(fishualize)
  library(tidybayes)

  NPer <- read.csv('20 07 24 AgNP ELA lakes fish excretion.csv',
                   stringsAsFactors = F, na.strings = c("", "NA", "."), 
                   strip.white = TRUE, sep = ",")
  Pmd <- read.csv('21 05 17 AgNP fish excretion model_R.csv',
                    stringsAsFactors = F, na.strings = c("", "NA", "."), 
                    strip.white = TRUE, sep = ",")
  paramm <- read.csv('22 02 11 param_modelsp_FishStoich.csv',
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
  
  # summary statistics info for paper methods & discussion ----
  
  NPexcr.ss1 <- NPexcr %>% group_by(Lake) %>%
    summarize(Mass.min = min(Mass, na.rm = TRUE),
              Mass.max = max(Mass, na.rm = TRUE),
              Mass.av = mean(Mass, na.rm = TRUE),
              Mass.sd = sd(Mass, na.rm = TRUE),
              massnorm.N.excr.min = min(massnorm.N.excr, na.rm = TRUE),
              massnorm.N.excr.max = max(massnorm.N.excr, na.rm = TRUE),
              massnorm.N.excr.sd = sd(massnorm.N.excr, na.rm = TRUE),
              massnorm.P.excr.min = min(massnorm.P.excr, na.rm = TRUE),
              massnorm.P.excr.max = max(massnorm.P.excr, na.rm = TRUE),
              massnorm.P.excr.sd = sd(massnorm.P.excr, na.rm = TRUE))
  
  NPexcr.ss2 <- NPexcr %>% group_by(Lake, Year) %>%
    summarize(Tag.excr.min = min(massnorm.Tag.excr, na.rm = T),
              Tag.excr.max = max(massnorm.Tag.excr, na.rm = T),
              Tag.excr.sd = sd(massnorm.Tag.excr, na.rm = T),
              P.excretion.av = mean(P.excretion, na.rm = T),
              WMass.av = mean(Mass..g., na.rm = T),
              DMass.av = mean(Mass, na.rm = T),
              Temp.av = mean(Temperature, na.rm = T))
  
  # NPexcr %>% summarize(min(Fork.length..mm., na.rm = TRUE))
  # NPexcr %>% summarize(sd(Fork.length..mm., na.rm = TRUE))
  # fork.2012 <- NPexcr %>% filter(Year == '2012',
  #                                Lake == '239') 
  # fork.2014 <- NPexcr %>% filter(Year == '2014') 
  # fork.2015 <- NPexcr %>% filter(Year == '2015') 
  # quantile(fork.2012$Fork.length..mm., c(0,0.1,0.5,0.75), na.rm = T)
  # NPexcr %>% filter (Year == '2015') %>% count (massnorm.N.excr)

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
           Log.massnorm.N.excr = log(massnorm.N.excr),
           Log.massnorm.P.excr = log(massnorm.P.excr),
           massnorm.NP.excr = (massnorm.N.excr/massnorm.P.excr)/(14/31),
           massnorm.CN.excr = (massnorm.C.excr/massnorm.N.excr)/(12/14),
           massnorm.CP.excr = (massnorm.C.excr/massnorm.P.excr)/(12/31),
           massnorm.NP.excr2 = Excreted.NP/Mass)
  
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
              massnorm.CP.excr.sd = sd(massnorm.CP.excr, na.rm = TRUE))

  
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
  qqPlot(lm(massnorm.N.excr ~ Lake*Year, 
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
  
  # ..N excretion anova + figure ----
  # testing linear model with log-transformed data
  lm.Nxnorm <- lm(log(massnorm.N.excr) ~ Lake*Year, 
     data = NPexcr %>% filter(Year != '2014'))
  Anova(lm.Nxnorm)
  
  # pairwise comparisons using emmeans
  lm.Nxnorm_emmeans <- emmeans(lm.Nxnorm, ~ Lake|Year, type = 'response',
                               bias.adj = T)
  pairs(lm.Nxnorm_emmeans)
  emmip(lm.Nxnorm, Lake ~ Year)
  lm.Nxnorm_emmeans
  plot(lm.Nxnorm_emmeans, comparison = T)
  
  # extracting effects from emmeans
  lm.Nxnorm_emmeansdf <- as.data.frame(summary(lm.Nxnorm_emmeans))
  
  # calculating effect sizes + extracting effects
  eff_sizeNx <- eff_size(lm.Nxnorm_emmeans, sigma = sigma(lm.Nxnorm), 
                         edf = df.residual((lm.Nxnorm)))
  eff_sizeNxdf <- as.data.frame(summary(eff_sizeNx))
  eff_sizeNxdf <- eff_sizeNxdf %>% mutate(effect.size = 
                                             ifelse(effect.size < 0, 
                                                    effect.size*(-1), 0.8784323),
                                          lower.CL = 
                                            ifelse(lower.CL < -0.8,
                                                   lower.CL*(-1), -0.03054899),
                                          upper.CL = 
                                            ifelse(upper.CL < 0.5,
                                                   upper.CL*(-1), 1.7874136
))
  
  # ...data visualization ----
  # .....N excretion vs. fish mass from 2012 to 2015 ----
  # Lake 222
  windows(width = 16, height = 20)
  ggplot(NPexcr %>%  filter(Year != '2014'), 
         aes(x = Mass, y = N.excretion)) +
    geom_point(size = 2) +
    theme_classic(base_size = 20) +
    facet_grid(rows = vars(Lake), cols = vars(Year)) +
    theme(axis.text = element_text(face = 'bold'),
          axis.line = element_line(size = 1),
          panel.background = element_rect(colour = 'black'),
          panel.grid = element_blank(),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'),
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') +
    #scale_fill_manual(name = 'Lake',
                      #labels = c('AgNPs 222', 'Reference 239'),
                      #values = c("black","gray60")) +
    scale_y_continuous(name = 'TDN excretion (μg N/g/h)') +
    scale_x_continuous(name = 'Mass (g)') 
  
  
  # .....time series 2012 to 2015 using model ----
  windows(width = 14, height = 7)
  Nexcrp <- ggplot(lm.Nxnorm_emmeansdf, 
         aes(x = Year, y = response, color = Lake,
             group = Lake)) +
    geom_point(size = 4, position = position_dodge(0.5)) +
    geom_errorbar(aes(ymax = lower.CL, ymin = upper.CL), 
                  width = 0.2, lwd = 1.5, position = position_dodge(0.5)) +
    geom_jitter(data = NPexcr %>% filter(Year != '2014'),
                aes(x = Year, y = massnorm.N.excr),
                size = 1.5, position = position_jitterdodge(jitter.width = 0.2),
                show.legend = FALSE, alpha = 0.5) +
    # scale_x_discrete('', labels = c('', '')) +
    theme_classic(base_size = 26) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face="bold"),
          # axis.title.x = element_blank(),
          # axis.text.x = element_blank(),
          # axis.title.y = element_text(vjust = +1),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'top') +
    scale_colour_manual(name = 'Lake',
                      labels = c('AgNPs 222', 'Reference 239'),
                      values = c("black","gray60")) +
    scale_y_continuous(name = 'TDN excretion rate (μg N/g/h)') +
    annotate("text", x = 0.5, y = 9500, label = 'A', 
             size = 7, fontface = 'bold')
  Nexcrp
  
  # .....effect sizes from 2012 2015 using model ----
  windows(width = 14, height = 7)
  eff_sizeNxp <- ggplot(eff_sizeNxdf,
                        aes(x = Year, y = exp(effect.size))) +
    geom_point(size = 4) +
    geom_pointrange(aes(ymax = exp(lower.CL), ymin = exp(upper.CL)),lwd = 1.2) +
    #coord_flip() +
    theme_classic(base_size = 26) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face = "bold"),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') +
    scale_y_continuous(name = 'Effect size') +
    annotate("text", x = 0.5, y = 6, label = 'B', 
             size = 7, fontface = 'bold')
  eff_sizeNxp
  
  windows(width = 14, height = 7)
  ggarrange(Nexcrp, eff_sizeNxp, ncol = 2, widths = c(1.2, 1),
            legend = 'top', common.legend = T, align = 'h')
  ggsave('final figures/Fig 1.png', 
         width = 14, height = 7, 
         units = 'in', dpi = 300)
  
  # ..P excretion anova + figure ----
  # testing linear model with log-transformed data
  lm.Pxnorm <- lm(log(massnorm.P.excr) ~ Lake*Year, 
                  data = NPexcr)
  Anova(lm.Pxnorm)
  
  lm.Pxnorm_emmeans <- emmeans(lm.Pxnorm, ~ Lake|Year, type = 'response', 
                               bias.adj = TRUE)
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
  
  # ...data visualization ----
  # .....P excretion vs. fish mass from 2012 to 2015 ----
  # Lake 222
  windows(width = 16, height = 20)
  ggplot(NPexcr, 
         aes(x = Mass, y = P.excretion)) +
    geom_point(size = 2) +
    theme_classic(base_size = 20) +
    facet_grid(rows = vars(Lake), cols = vars(Year)) +
    theme(axis.text = element_text(face = 'bold'),
          axis.line = element_line(size = 1),
          panel.background = element_rect(colour = 'black'),
          panel.grid = element_blank(),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'),
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') +
    scale_y_continuous(name = 'TDP excretion (μg P/g/h)') +
    scale_x_continuous(name = 'Mass (g)') 
  
  # boxplot of P excretion in 2012, 2014, 2015
  windows(width = 14, height = 7)
  ggplot(NPexcr, aes(x = Lake, y = massnorm.P.excr, fill = Lake)) +
    theme_classic(base_size = 26) +
    geom_boxplot(alpha = 0.7, outlier.size = 3, outlier.colour = NULL) +
    #geom_jitter() +
    theme_bw(base_size = 20) +
    theme(text = element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          legend.position = 'none') +
    scale_fill_brewer(palette = "Accent") +
    scale_y_continuous(name = 'TDP excretion (μg P/g/h)') +
    scale_x_discrete(labels = c('AgNPs 222', 'Reference 239'))
  
  # boxplot P excretion vs. Lake and Year coded (for 2014/2015)
  windows(width = 14, height = 7)
  ggplot(NPexcr, aes(x = Lake, y = massnorm.P.excr, fill = Year)) +
    theme_classic(base_size = 26) +
    geom_boxplot(alpha = 0.7, outlier.size = 3, outlier.colour = NULL) +
    #geom_jitter() +
    theme_bw(base_size = 20) +
    theme(text = element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(name = 'TDP excretion (μg P/ind/h)') +
    scale_x_discrete(labels = c('AgNPs 222', 'Reference 239'))
  
  
  # .....time series 2012 to 2015 using model ----
  windows(width = 14, height = 7)
  Pexcrp <- ggplot(lm.Pxnorm_emmeansdf, 
         aes(x = Year, y = response, color = Lake,
             group = Lake)) +
    geom_point(size = 4, position = position_dodge(0.5)) +
    geom_errorbar(aes(ymax = lower.CL, ymin = upper.CL),
                  width = 0.2, lwd = 1.5, position = position_dodge(0.5)) +
    geom_jitter(data = NPexcr,
                aes(x = Year, y = massnorm.P.excr),
                size = 1.5, position = position_jitterdodge(jitter.width = 0.2),
                show.legend = FALSE, alpha = 0.5) +
    # scale_x_discrete('', labels = c('', '', '')) +
    theme_classic(base_size = 26) +
    theme(text = element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'top') +
    scale_y_continuous(name = 'TDP excretion rate (μg P/g/h)') +
    scale_colour_manual(name = 'Lake',
                        labels = c('AgNPs 222', 'Reference 239'),
                        values = c("black","gray60")) +
    annotate("text", x = 0.55, y = 230, label = 'A', 
             size = 7, fontface = 'bold')
  Pexcrp
  
  # .....effect sizes from 2012-2015 using model ----
  windows(width = 14, height = 7)
  eff_sizePxp <- ggplot(eff_sizePxdf,
         aes(x = Year, y = exp(-effect.size))) +
    geom_point(size = 4) +
    geom_pointrange(aes(ymax = exp(-lower.CL), ymin = exp(-upper.CL)),lwd = 1.2) +
    #coord_flip() +
    theme_classic(base_size = 26) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face = "bold"),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') +
    scale_y_continuous(name = 'Effect size') +
    annotate("text", x = 0.55, y = 39, label = 'B', 
             size = 7, fontface = 'bold')
  eff_sizePxp
    
  windows(width = 14, height = 7)
  ggarrange(Pexcrp, eff_sizePxp, ncol = 2, heights = c(1.2, 1),
            legend = 'top', common.legend = T, align = 'h')
  ggsave('final figures/Fig 2.png', 
         width = 14, height = 7, 
         units = 'in', dpi = 300)
  
  # .....modelled yellow perch excretion with P release data ----
  windows(width = 14, height = 7)
  MeasP.excretion <- Pmod$MeasP.excretion
  Pexcrmod.p <- ggplot(Pmod,aes(x = Ingestion.rate, y = ModP.excretion,
                                color = Lake, group = FoodCP)) +
    geom_line(aes(linetype = FoodCP), colour = 'black', 
              lwd = 1.2, show.legend = F) +
    labs(x = 'Ingestion rate (mg C/mg C/d)',
         y = 'P excretion rate (mg P/mg C/d)') +
    xlim(0, 0.15) +
    ylim(0, 0.005) +
    theme_classic(base_size = 26) +
    theme(text = element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.title.y = element_text(vjust = +1),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') +
    scale_linetype_manual(values = c('solid', 'solid')) +
    annotate("text", x = 0.13, y = 0.0011, label = 'Food P:C min',
             size = 6, fontface = 'bold') +
    annotate("text", x = 0.13, y = 0.0006, label = 'Food P:C max',
             size = 6, fontface = 'bold') +
    annotate("text", x = 0, y = 0.005, label = 'A', 
             size = 7, fontface = 'bold') +
    annotate("text", x = 0.15, y = 0.005, label = 'B', 
             size = 7, fontface = 'bold')
  Pexcrmod.p
  
  Pexcrmeas.p <- ggplot(Pmod, aes(x = Ingestion.rate, y = MeasP.excretion,
                                 color = Lake)) +
    geom_point(data = Pmod %>%  filter(Lake == 'L222'),
               aes(x = 0, shape = Year, color = Lake),
                               size = 7.5) +
    geom_point(data = Pmod %>%  filter(Lake == 'L239'),
               aes(x = 0, shape = Year, color = Lake),
               size = 5) +
    scale_color_manual(labels = c('AgNPs 222', 'Reference 239'),
                       values = c("black","gray60"), na.translate = F) +
    scale_shape_manual(labels = c('2012', '2014', '2015'),
                       values = c(15, 16, 17), na.translate = F) +
    theme_classic(base_size = 22) +
    theme(text = element_text(family = "Tahoma"),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.25, .25, .25, .25, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right',
          legend.spacing.y = unit(0, 'cm'),
          legend.background = element_blank(),
          legend.box.background = element_rect(linetype = 2, colour = "black"),
          legend.box.spacing = unit(2, 'lines')) 
  Pexcrmeas.p
  
  windows(width = 20, height = 9)
  ggarrange(Pexcrmod.p, Pexcrmeas.p, ncol = 2, widths = c(0.55, 0.45),
            legend = 'right', common.legend = F, align = 'h')
  ggsave('final figures/Fig 3.png', 
         width = 14, height = 7, 
         units = 'in', dpi = 300)
  
  # ..Tag excretion anova + figure ----
  aov.Tag <- aov(massnorm.Tag.excr ~ Year, 
               data = NPexcr)
  summary(aov.Tag)
  
  glm.Tagx <- glm(Tag.excretion ~ Year*Mass, 
                 data = NPexcr, family = Gamma(link = "log"))
  
  summary(glm.Tagx)
  
  # ...data visualization ----
  windows(width = 14, height = 7)
  ggplot(NPexcr.ts %>% filter(Year != 2012), 
         aes(x = Year, y = Tag.excr.av)) +
    geom_errorbar(aes(ymin = Tag.excr.av - Tag.excr.sd,
                      ymax = Tag.excr.av + Tag.excr.sd), 
                  width = 0.2, lwd = 1.5, position = position_dodge(0.5)) +
    geom_bar(stat = "identity", width = 0.5, color="black") +
    # geom_point(size = 5, position = position_dodge(0.5)) +
    # geom_boxplot(outlier.shape = NA) +
    geom_jitter(data = NPexcr %>%  filter(Year != 2012),
                aes(x = Year, y = massnorm.Tag.excr),
                size = 4, position = position_jitter(0.2),
                show.legend = FALSE, alpha = 0.5) +
    theme_classic(base_size = 22) +
    theme(text = element_text(family = "Tahoma"),
          axis.title = element_text(face = "bold"),
          legend.position = 'none') +
    #scale_fill_manual(values = c("black","gray60")) +
    scale_y_continuous(name = 'TAg excretion (μg TAg/g/h)') +
    scale_x_discrete(labels = c('2014', '2015'))
  
  # ..excreted N:P (molar) ----
  aov.NPx <- aov(massnorm.NP.excr ~ Lake*Year, 
                data = NPexcr %>% filter(Year!= '2014'))
  summary(aov.NPx)
  
  #LM
  lm.NPxnorm <- lm(log(massnorm.NP.excr) ~ Lake*Year, 
                     data = NPexcr %>% filter(Year != '2014'))
  summary(lm.NPxnorm)
  Anova(lm.NPxnorm)
  
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
  windows(width = 14, height = 7)
  ggplot(NPexcr %>%  filter(Year =='2015'), 
         aes(x = Lake, y = massnorm.NP.excr, fill = Lake)) +
    theme_classic(base_size = 26) +
    geom_boxplot(alpha = 0.7, outlier.size = 3, outlier.colour = NULL) +
    #geom_jitter() +
    theme_bw(base_size = 20) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face = "bold"),
          legend.position = 'none') +
    scale_fill_brewer(palette = "Accent") +
    scale_y_continuous(name = 'Excreted N:P (molar)') +
    scale_x_discrete(labels = c('AgNPs 222', 'Reference 239'))
  
  # .....time series 2012 to 2015 using model ----
  windows(width = 14, height = 7)
  NPexcrp <- ggplot(lm.NPxnorm_marginal_means, 
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
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') +
    scale_fill_brewer(palette = "Accent") +
    scale_y_continuous(name = 'Excreted N:P (molar)') +
    scale_colour_manual(name = 'Lake',
                          labels = c('AgNPs 222', 'Reference 239'),
                          values = c("black","gray60"))+
    annotate("text", x = 0.5, y = 440, label = 'A', 
             size = 7, fontface = 'bold')
  NPexcrp
  
  # .....effect sizes from 2012 2015 using model ----
  windows(width = 14, height = 7)
  eff_sizeNPxp <- ggplot(eff_sizeNPxdf,
                        aes(x = Year, y = exp(effect.size))) +
    geom_point(size = 4) +
    geom_pointrange(aes(ymax = exp(lower.CL), ymin = exp(upper.CL)),lwd = 1.2) +
    #coord_flip() +
    theme_classic(base_size = 26) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face = "bold"),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') +
    scale_y_continuous(name = 'Effect size') +
    annotate("text", x = 0.5, y = 6.7, label = 'B', 
             size = 7, fontface = 'bold') 
  eff_sizeNPxp
  
  windows(width = 14, height = 7)
  ggarrange(NPexcrp, eff_sizeNPxp, ncol = 2, heights = c(1.2, 1),
            legend = 'top', common.legend = T, align = 'h')
  ggsave('final figures/Fig 4.png', 
         width = 14, height = 7, 
         units = 'in', dpi = 300)
  
  # ..excreted C:N (molar) in 2015 ----
  aov.CN <- aov(massnorm.CN.excr ~ Lake, 
               data = NPexcr)
  summary(aov.CN)
  
  windows(width = 14, height = 7)
  ggplot(NPexcr %>%  filter(Year =='2015'), 
         aes(x = Lake, y = Excreted.CN, fill = Lake)) +
    theme_classic(base_size = 26) +
    geom_boxplot(alpha = 0.7, outlier.size = 3, outlier.colour = NULL) +
    #geom_jitter() +
    theme_bw(base_size = 20) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face = "bold"),
          legend.position = 'none') +
    scale_fill_brewer(palette = "Accent") +
    scale_y_continuous(name = 'Excreted C:N (molar)') +
    scale_x_discrete(labels = c('AgNPs 222', 'Reference 239'))
  
  # ..excreted C:P (molar) in 2015 ----
  aov.CP <- aov(massnorm.CP.excr ~ Lake, 
                data = NPexcr)
  summary(aov.CP)
  
  windows(width = 14, height = 7)
  ggplot(NPexcr %>%  filter(Year =='2015'), 
         aes(x = Lake, y = Excreted.CP, fill = Lake)) +
    theme_classic(base_size = 26) +
    geom_boxplot(alpha = 0.7, outlier.size = 3, outlier.colour = NULL) +
    #geom_jitter() +
    theme_bw(base_size = 20) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face = "bold"),
          legend.position = 'none') +
    scale_fill_brewer(palette = "Accent") +
    scale_y_continuous(name = 'Excreted C:P (molar)') +
    scale_x_discrete(labels = c('AgNPs 222', 'Reference 239'))
  
  # FishStoichModel ----
  parameters <- model_parameters("Perca flavescens", otolith = F,
                                           family = "Percidae", temp = 18, mirror = "se")
  # find length-weight relationship and caudal fin aspect ratio
  find_lw("Perca flavescens", mirror = "se")
  aspect_ratio("Perca flavescens")
  # find growth parameters
  growth_params("Perca flavescens", otolith = FALSE)
  
  # Run model
  # Yellow perch
  total.length <- exp((log(NPexcr$Mass..g./parammod$lwa_m)/parammod$lwb_m))
  # total.length <- unique(NPexcr$Total.length.cm)
  # tl <- total.length[!is.na(total.length)]
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
  iter_yp$tl <- rep(total.length, each = 500)
  iter_yp$mass <- rep(NPexcr$Mass..g., each = 500)
  iter_yp <- iter_yp %>% mutate(N.excretion.mod = Fn*10^6/24,
                                P.excretion.mod = Fp*10^6/24,
                                C.ingestion.mod = Ic*10^6/24)
  
  # ..Figure 3: plot excretion ----
  # plot N excretion
  plot_N <- 
    ggplot(group_by(iter_yp, iter), aes(x = mass, y = N.excretion.mod)) +
    stat_lineribbon(alpha = 0.8, show.legend = F) +
    scale_fill_brewer() +
    theme_classic(base_size = 26) +
    geom_point(aes(x = Mass..g., y = N.excretion, color = Lake, shape = Year), 
               data = NPexcr %>%  filter(Year != '2014'),
               size = 4) +
    labs(x = "Mass (g)", y = "N excretion rate (μg N/ind/h)") +
    scale_colour_manual(name = 'Lake',
                        labels = c('AgNPs 222', 'Reference 239'),
                        values = c("black","gray60")) +
    scale_shape_manual(labels = c('2012', '2015'),
                       values = c(16, 15, 17), na.translate = F) +
    theme(text = element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'none') 
  plot_N
  
  # plot P excretion
  plot_P <- 
    ggplot(group_by(iter_yp, iter), aes(x = mass, y = P.excretion.mod)) +
    stat_lineribbon(alpha = 0.8, show.legend = F) +
    scale_fill_brewer() +
    theme_classic(base_size = 26) +
    geom_point(aes(x = Mass..g., y = P.excretion, color = Lake, shape = Year), 
               data = NPexcr, size = 4) +
    labs(x = "Mass (g)", y = "P excretion rate (μg P/ind/h)") +
    scale_colour_manual(name = 'Lake',
                        labels = c('AgNPs 222', 'Reference 239'),
                        values = c("black","gray60")) +
    scale_shape_manual(labels = c('2012', '2014', '2015'),
                       values = c(16, 17, 15), na.translate = F) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face="bold"),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') 
  plot_P
  
  # get legend function
  get_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  
  fig4.legend <- get_legend(plot_P)
  
  windows(width = 14, height = 7)
  ggarrange(plot_N, plot_P, ncol = 2,
            legend = 'right', common.legend = F, align = 'h', 
            legend.grob = fig4.legend)
  ggsave('final figures/Fig 4.1.png', 
         width = 16, height = 7, 
         units = 'in', dpi = 300)
  
  # ..Figure S3 ----
  # plot sensitivity to changes in diet on CNP limitation
  sensitivity <- function(TL, param, iter = 1000, par,
                          out = c("Ic", "In", "Ip", "Gc",
                                  "Gn", "Gp", "Fc", "Fn",
                                  "Fp", "Wc", "Wn", "Wp"), ...) {
    
    #parameter SD's and means
    pm <- c("lt_m", "ac_m", "an_m", "ap_m", "Dc_m",
            "Dn_m", "Dp_m", "linf_m", "k_m", "t0_m",
            "theta_m", "r_m", "h_m", "lwa_m", "lwb_m",
            "mdw_m", "v_m", "F0nz_m", "F0pz_m", "Qc_m",
            "Qn_m", "Qp_m", "a_m", "f0_m")
    psd <- c("lt_sd", "ac_sd", "an_sd", "ap_sd", "Dc_sd",
             "Dn_sd", "Dp_sd", "linf_sd", "k_sd", "t0_sd",
             "theta_sd", "r_sd", "h_sd", "lwa_sd", "lwb_sd",
             "mdw_sd", "v_sd", "F0nz_sd", "F0pz_sd", "Qc_sd",
             "Qn_sd", "Qp_sd", "a_sd", "f0_sd")
    
    parm <- par[par %in% pm]
    parsd <- par[par %in% psd]
    
    param_m <- param[parm]
    param_sd <- param[parsd]
    
    sd_low <- 0.000000001
    
    param_sdl <- param_sd
    param_sdl[seq_len(length(parsd))] <- sd_low
    param_msdl <- append(param_m, param_sdl)
    
    #run cnp_model for all sd's with rest very low
    sd <- parsd
    res_sd <- as.data.frame(
      mcmapply(sd, FUN = function(x) {
        param_msdl[x] <- param_sd[x]
        mod <- cnp_model_mcmc(TL, param_msdl, iter, ...)$summary
        ext <- mod[match(out, mod$variable), "Q_97.5"] - mod[match(out, mod$variable), "Q_2.5"]
        ext
      }))
    
    row.names(res_sd) <- sapply(out, function(x) {
      rn <- paste(x, "_CI", sep = "")
      rn
    })
    
    res_sd <- as.data.frame(t(res_sd))
    
    #plot
    res <- res_sd
    res$input_sd <- row.names(res)
    res <- gather(res, "key", "value", - input_sd)
    sum <- summarise(group_by(res, key), sum = sum(value))
    res <- res %>% left_join(sum) %>%
      mutate(scale = value / sum)
    plot <- ggplot(res) +
      geom_tile(aes(x = key, y = input_sd, fill = scale)) +
      scale_fill_fish(option = "Trimma_lantana", end = 0.9) +
      geom_text(aes(x = key, y = input_sd,
                    label = formatC(value, format = "e", digits = 1))) +
      labs(x = "", y = "", fill = "Relative width 95% CI") +
      theme_bw(base_size = 18) +
      scale_color_fish_d(option = "Trimma_lantana", end = 0.8) +
      theme(text = element_text( family = "Tahoma"),
            axis.title = element_text(face="bold"),
            legend.title = element_text(face = 'bold'),
            legend.margin = margin(.15, .15, .15, .15, 'cm'),
            legend.key.height = unit(2, 'lines'), 
            legend.key.width = unit(3, 'lines'),
            legend.position = 'bottom') 
    print(plot)
    res_sd
  }
    
  windows(width = 14, height = 7)
  plot_sens <- sensitivity(TL = 6, param = parammod,
              par = c("Dn_sd","Dc_sd","Dp_sd"), out = c("Fn", "Fp", "Ic"))
  ggsave('final figures/Fig S4.png', 
         width = 7, height = 7, 
         units = 'in', dpi = 300)
  
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
  Dpm <- Dnm / Snp
  
  
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
  simd <- simd %>% mutate(P.excretion.mod = Fp_median*10^6/24)
  
  p1 <- 
    ggplot()+
    geom_segment(aes(x = 1.2, y = 1.2/Snp, xend = Dnm, yend = Dpm), size = 2) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = Dnm, yend = 3.2/Snp), size = 2) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = 3.2, yend = Dpm), size = 2) +
    geom_text(aes(x = 1.7, y = 1, label = "N"), size = 15) +
    geom_text(aes(x = 2.7, y = 1.1, label = "C"), size = 15) +
    geom_text(aes(x = 2.4, y = 0.7, label = "P"), size = 15) +
    labs( x = "Dn (%)", y = "Dp (%)") +
    theme_bw() +
    theme(aspect.ratio=1, title = element_text(size = 16), 
          text = element_text(size = 14)) 
  
  p2 <-
    ggplot(simd) +
    geom_raster(aes(x = Dn, y = Dp, fill = Ic_median)) +
    scale_fill_fish(option = "Trimma_lantana", trans = "sqrt") +
    geom_segment(aes(x = 1.2, y = 1.2/Snp, xend = Dnm, yend = Dpm), size = 2) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = Dnm, yend = 3.2/Snp), size = 2) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = 3.2, yend = Dpm), size = 2) +
    theme_bw() +
    theme(aspect.ratio=1, title = element_text(size = 16), 
          text = element_text(size = 14)) +
    labs(x = "Dn (%)", y = "Dp (%)", fill = "Ic (g/day)") 
  
  p3 <-
    ggplot(simd) +
    geom_raster(aes(x = Dn, y = Dp, fill = Fn_median)) +
    # geom_segment(aes(x = 1.2, y = 1.2/Snp, xend = Dnm, yend = Dpm), size = 2) +
    # geom_segment(aes(x = Dnm, y = Dpm, xend = Dnm, yend = 3.2/Snp), size = 2) +
    # geom_segment(aes(x = Dnm, y = Dpm, xend = 3.2, yend = Dpm), size = 2) +
    scale_fill_fish(option = "Trimma_lantana", trans = "sqrt") +
    theme_bw() +
    theme(aspect.ratio=1, title = element_text(size = 16), 
          text = element_text(size = 14))  +
    labs(x = "Dn (%)", y = "Dp (%)", fill = "Fn (g/day)") 
  p3
  
  p4 <-
    ggplot(simd) +
    geom_raster(aes(x = Dn, y = Dp, fill = P.excretion.mod)) +
    scale_fill_fish(option = "Trimma_lantana", trans = "sqrt") +
    # geom_segment(aes(x = 1.2, y = 1.2/Snp, xend = Dnm, yend = Dpm), size = 2) +
    # geom_segment(aes(x = Dnm, y = Dpm, xend = Dnm, yend = 3.2/Snp), size = 2) +
    # geom_segment(aes(x = Dnm, y = Dpm, xend = 3.2, yend = Dpm), size = 2) +
    theme_bw() +
    theme(aspect.ratio=1, title = element_text(size = 16), 
          text = element_text(size = 14)) +
    labs(x = "Dn (%)", y = "Dp (%)", fill = "P excretion (ug/ind/h)") 
  p4
  
  # plot N/P/C limitation, N/P excretion, Ingestion rate, Growth rate
  windows(width = 14, height = 7)
  
  plot_lim <- function(mod){
    
    lim <- lapply(mod$stanfit, function(x){
      ee <- rstan::extract(x,"lim")[[1]]
      c <- length(which(ee==1))/length(ee)
      n <- length(which(ee==2))/length(ee)
      p <- length(which(ee==3))/length(ee)
      return(data.frame(c = c,
                        n = n,
                        p = p))
    }) %>%
      
      dplyr::bind_rows() %>%
      dplyr::mutate(tl = unique(mod$summary$TL)) %>%
      
      tidyr::gather("nutrient", "prop_lim", - tl)
    
    plot <- ggplot(lim) +
      geom_point(aes(x = tl, y = prop_lim, color = nutrient), size = 2) +
      geom_line(aes(x = tl, y = prop_lim, color = nutrient), size = 1) +
      labs(x = "Total length (cm)", y = "Proportion of iterations", color = "Limiting element") +
      theme_classic(base_size = 18) +
      scale_color_fish_d(option = "Centropyge_loricula", end = 0.8) +
      theme(text = element_text( family = "Tahoma"),
            axis.title = element_text(face="bold"),
            legend.title = element_text(face = 'bold'),
            legend.margin = margin(.15, .15, .15, .15, 'cm'),
            legend.key.height = unit(2, 'lines'), 
            legend.key.width = unit(3, 'lines'),
            legend.position = 'bottom') 
    print(plot)
    
    return(plot)
    
  }
  
  plot_lim <- plot_lim(FStoichm)
  ggarrange(plot_lim, plot_sens, ncol = 2,
            legend = 'right', common.legend = F, align = 'h')
  ggsave('final figures/Fig S3.png', 
         width = 7, height = 5, 
         units = 'in', dpi = 300)
  
  # ..Figure S2: Ingestion and N/P/C limitation ----
  # plot C ingestion
  plot_Ic <- 
    ggplot(group_by(iter_yp, iter), aes(x = mass, y = C.ingestion.mod)) +
    stat_lineribbon(alpha = 0.8, show.legend = F) +
    scale_fill_brewer() +
    theme_classic(base_size = 18) + 
    labs(x = "Mass (g)", y = "C ingestion rate (μg C/ind/h)") 
  plot_Ic
  
  ggsave('final figures/Fig S2.png', 
         width = 7, height = 5, 
         units = 'in', dpi = 300)
  
  
  plot_cnp(FStoichm,  y = "Fp", x = "tl", probs = c(0.5, 0.8, 0.95)) +
    geom_point(data = NPexcr, aes(x = Total.length.cm, 
                y = P.excretion.d, color = Year, 
                shape = Lake), size = 3)
  
  
  fishflux::plot_cnp(FStoichm,  y = "Gp", x = "tl", probs = c(0.5, 0.8, 0.95))
  fishflux::plot_cnp(FStoichm,  y = "Ic", x = "tl", probs = c(0.5, 0.8, 0.95))
  fishflux::plot_cnp(FStoichm,  y = c("Fp", "Gp"), x = "tl", probs = c(0.5))
  
  
  
  sessionInfo() 
  #################################### END OF CODE ##################################
  
             