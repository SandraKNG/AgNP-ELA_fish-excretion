  # NAME of paper
  # this code was developped by S. Klemet-N'Guessan in 2020 and 2021

  # load libraries and read silver nanoparticles (NP) dataset ----
  
  library(tidyverse)
  library(car)
  library(RColorBrewer)
  library(emmeans) # post-hoc pairwise comparisons using lm, glm, levDev models
  #library(effects) # effects from a fitted model
  library(ggpubr) #to add multiple graphs on the same page
  library(fishflux)

  NPer <- read.csv('20 07 24 AgNP ELA lakes fish excretion.csv',
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
           P.excretion.d = (P.excretion*24)/10^6,
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
  
  NPexcr %>% summarize(min(Fork.length..mm., na.rm = TRUE))
  NPexcr %>% summarize(sd(Fork.length..mm., na.rm = TRUE))
  fork.2012 <- NPexcr %>% filter(Year == '2012',
                                 Lake == '239') 
  fork.2014 <- NPexcr %>% filter(Year == '2014') 
  fork.2015 <- NPexcr %>% filter(Year == '2015') 
  quantile(fork.2012$Fork.length..mm., c(0,0.1,0.5,0.75), na.rm = T)
  NPexcr %>% filter (Year == '2015') %>% count (massnorm.N.excr)

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
  # example
  param_zebsco <- param_zebsco
  model <- cnp_model_mcmc(TL = 5:20, param = param_zebsco)
  fishflux::extract(model, c("Fn","Fp"))
  # Yellow perch
  FStoichm <- cnp_model_mcmc(TL = 4:11, param = parammod)
  output <- fishflux::extract(FStoichm, c("Fn","Fp", "Ic", "Gp"))
  
  # plot N/P/C limitation, N/P excretion, Ingestion rate, Growth rate
  windows(width = 12, height = 7)
  limitation(FStoichm)
  plot_cnp(FStoichm,  y = "Fp", x = "tl", probs = c(0.5, 0.8, 0.95)) +
    #ylim(0,0.00125) +
    geom_point(data = NPexcr, aes(x = Total.length.cm, 
                y = P.excretion.d, color = Year, 
                shape = Lake), size = 2)
  fishflux::plot_cnp(FStoichm,  y = "Fn", x = "tl", probs = c(0.5, 0.8, 0.95)) +
    geom_point(data = NPexcr, aes(x = Total.length.cm, 
                                  y = N.excretion.d, 
                                  color = Year, shape = Lake), size = 2)
  fishflux::plot_cnp(FStoichm,  y = "Gp", x = "tl", probs = c(0.5, 0.8, 0.95))
  fishflux::plot_cnp(FStoichm,  y = "Ic", x = "tl", probs = c(0.5, 0.8, 0.95))
  fishflux::plot_cnp(FStoichm,  y = c("Fp", "Gp"), x = "tl", probs = c(0.5))
  
  # P excretion vs total length
  ggplot(output, aes(x = Total.length.cm, 
                     y = P.excretion.d)) +
    geom_point(NPexcr, aes(x = Total.length.cm, 
                           y = P.excretion.d), size = 2) +
    labs(x = "Total length (cm)",
         y = "P excretion (g/day)") +
    xlim(3.8,11) +
    theme_bw()
  
  # N excretion vs total length
  ggplot(NPexcr, aes(x = Total.length.cm, 
                     y = N.excretion.d)) +
    geom_point(size = 2) +
    labs(x = "Total length (cm)",
         y = "N excretion (g/day)") +
    xlim(3.8,11) +
    theme_bw()
  
  
  sessionInfo() 
  #################################### END OF CODE ##################################
  
             