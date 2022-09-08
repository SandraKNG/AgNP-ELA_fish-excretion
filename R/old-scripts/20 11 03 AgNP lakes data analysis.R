  
  # load libraries and read silver nanoparticles (NP) dataset ----
  
  library(tidyverse)
  #library(nortest)
  #library(descr)
  library(car)
  library(RColorBrewer)
  #library(moments)
  #library(MASS)
  #library(AEDForecasting)
  library(emmeans) # post-hoc pairwise comparisons using lm, glm, levDev models
  library(effects)
  
  ### Below I have written a general purpose function to generate Levene's deviates
  LeveneDeviates <- function(y, group, med=TRUE, log_trans=TRUE) {
    
    #log transform data?
    if (log_trans==TRUE)
      y = log(y)
    
    # allows for use of mean or median as measure of central tendency
    if (med==TRUE)
      meds <- tapply(y, group, median, na.rm=TRUE)
    else 
      meds <- tapply(y, group, mean, na.rm=TRUE) 
    
    # calculates deviates for each observation from a measure of central tendency for a "group"
    abs(y - meds[group])}

  NPer <- read.csv('20 07 24 AgNP ELA lakes fish excretion.csv',
                 stringsAsFactors = F, na.strings = c("", "NA", "."), 
                 strip.white = TRUE, sep = ",")  

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
           Lake = as.factor(Lake),
           Year = as.factor(Year))
  
  str(NPexcr)
  
  NPexcr.ts <- NPexcr %>% group_by(Lake, Year) %>%
    summarize(Tag.excr.av = mean(massnorm.Tag.excr, na.rm = T),
              Tag.excr.sd = sd(massnorm.Tag.excr, na.rm = T))
  
  NPexcr %>% summarize(max(Total.length..mm., na.rm = TRUE))
  NPexcr %>% summarize(sd(Mass, na.rm = TRUE))
  fork.2012 <- NPexcr %>% filter(Year == '2012') 
  fork.2014 <- NPexcr %>% filter(Year == '2014') 
  fork.2015 <- NPexcr %>% filter(Year == '2015') 
  quantile(fork.2015$Fork.length..mm., c(0,0.1,0.5,0.75), na.rm = T)

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
           Log.massnorm.N.excr = massnorm.N.excr,
           Log.massnorm.P.excr = massnorm.P.excr,
           massnorm.NP.excr = (massnorm.N.excr/massnorm.P.excr)/(14/31),
           massnorm.CN.excr = (massnorm.C.excr/massnorm.N.excr)/(12/14),
           massnorm.CP.excr = (massnorm.C.excr/massnorm.P.excr)/(12/31))
  
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
  
  
  ################################### DATA ANALYSIS #################################
  
  # ..fish size ----
  
  lm.size <- lm(Total.length..mm. ~ Lake*Year,
                  data = NPexcr)
  summary(lm.size)
  
  lm.size_emmeans <- emmeans(lm.size, ~ Lake:Year)
  pairs(lm.size_emmeans)
  
  # ..N excretion anova + figure ----
  # post-addition N excretion in 2015 vs. Lake
  aov.Nnx <- aov(massnorm.N.excr ~ Lake, 
               data = NPexcr %>% filter(Year == '2015'))
  # not pooling 2014 + 2015 together may be a difference between both
  # not looking at year effect (2014 vs. 2015) 
  # because in 2014, 1 and 6 data points for L239 and L222 respectively
  summary(aov.Nnx)
  aov.Nx <- aov(massnorm.N.excr ~ Lake*Year, 
               data = NPexcr %>% filter(Year != '2014')) 
  summary(aov.Nx)
  # (pre in 2012 vs. post-addition in 2015 N excretion)*/+(Lake)
  # excluding 2014 because of low sample sizes
  
  # testing linear model with log-transformed data
  lm.Nxnorm <- lm(log(massnorm.N.excr) ~ Lake*Year, 
     data = NPexcr %>% filter(Year != '2014'))
  
  Anova(lm.Nxnorm)
  
  lm.Nxnorm_emmeans <- emmeans(lm.Nxnorm, ~ Lake:Year, type = 'response')
  pairs(lm.Nxnorm_emmeans)
  emmip(lm.Nxnorm, Lake ~ Year)
  lm.Nxnorm_emmeans
  plot(lm.Nxnorm_emmeans, comparison = T)
  
  # extracting effects
  #lm.Nxnorm_marginal_means <- as.data.frame(Effect(c("Lake", "Year"),
                                                    #lm.Nxnorm, type = 'response'))
  lm.Nxnorm_marginal_means <- as.data.frame(summary(lm.Nxnorm_emmeans))
  
  # testing Generalized Linear Models with Gamma distribution
  glm.Nxnorm <- glm(massnorm.N.excr ~ Lake*Year, 
                data = NPexcr %>% filter(Year != '2014'), 
                family = Gamma(link = "log"))
  summary(glm.Nxnorm)
  Anova(glm.Nxnorm)
  
  # extracting effects
  glm.Nxnorm_marginal_means <- as.data.frame(Effect(c("Lake", "Year"),
                                                    glm.Nxnorm))
  
  glm.Nxnorm_emmeans <- emmeans(glm.Nxnorm, ~ Year|Lake, type = 'response')
  glm.Nxnorm_emmeans_pairwise <- emmeans(glm.Nxnorm, pairwise ~ Lake:Year, 
                                         type = 'response')
  glm.Nxnorm_emmeans
  glm.Nxnorm_emmeans_pairwise$contrasts
  pairs(glm.Nxnorm_emmeans)
  summary(glm.Nxnorm_emmeans)
  windows(width = 14, height = 7)
  plot(glm.Nxnorm_emmeans, comparisons = T)
  emmip(glm.Nxnorm, Lake ~ Year)
  # Confidence intervals for EMMs are based on the statistical properties of the individual EMMs, whereas comparison arrows are based on the statistical properties of differences of EMMs.

  # use Gamma distribution because that model is a good start for positive continuous response variables.
  # log link is important here as the default inverse link doesn't make any sense
  # You need to think about range of possible values when choosing a distribution. Gaussian to t distributions allow for negative N excretion which makes no physical sense. You want a distribution for positive reals if we can assume none of the observations is censored... 
  # A Gamma distribution or inverse Gaussian are choices in the GLM world for positive real values responses. If 0 is allowed (not here) then Tweedie might work. If you have censored data, y_i < Level of Detection, then you need to fit a censored model, say censored Gamma 
  # The key thing to understand is that the natural log link function is a function of the mean of y, not the y values themselves.
  # The log link exponentiates the linear predictors. It does not log transform the outcome variable.
  # When we transform the data in a linear model, we are no longer claiming that y is normally distributed around a mean, given the x values — we are claiming that our new outcome variable, ln(yi), is normally distributed.
  # In fact, we often make this transformation specifically because the values of y do not appear to be normally distributed around their average.
  
  # Generate Levene's deviate 
  # for each group, calculates median, 
  # then look at deviation of each observation from median
  NPexcr$LevDevN <- with(NPexcr, 
                         LeveneDeviates(massnorm.N.excr,
                                        group = interaction(Lake, Year, drop = TRUE),
                                        med  = TRUE, log_trans = TRUE))
  # add 10% of the smallest value
  quantile(NPexcr$LevDevN, probs = c(0, .01, .05, 0.025, .5), na.rm = TRUE)
  # smallest value = 0.002616777
  NPexcr$LevDevN2 <- NPexcr$LevDevN + 0.0003
  glm.LevDevN_gamma <- glm(LevDevN2 ~ Lake*Year,
                          data = NPexcr %>% filter(Year != '2014'), 
                          family = "Gamma")
  summary(glm.LevDevN_gamma)
  Anova(glm.LevDevN_gamma)
  
  glm.LevDevN_emmeans <- emmeans(glm.LevDevN_gamma, ~ Lake|Year, type = 'response')
  plot(estimate_modelN)
  
  # testing with linear model
  lm.LevDevN <- lm(LevDevN ~ Lake*Year,
                           data = NPexcr %>% filter(Year != '2014'))
  summary(lm.LevDevN)
  Anova(lm.LevDevN)
  
  # after looking at both GLM and LM for Levene's deviates, not really concerned about heterogeneity of group variances..why??
  
  # ...data visualization ----
  # boxplot of N excretion in 2015 
  windows(width = 14, height = 7)
  ggplot(NPexcr %>% filter(Year != '2012'), 
         aes(x = Lake, y = massnorm.N.excr, fill = Lake)) +
    theme_classic(base_size = 26) +
    geom_boxplot(alpha = 0.7, outlier.size = 3, outlier.colour = NULL) +
    theme_bw(base_size = 20) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face="bold"),
          legend.position = 'none') +
    scale_fill_brewer(palette = "Accent") +
    scale_y_continuous(name = 'TDN excretion (μg N/g/h)') +
    scale_x_discrete(labels = c('AgNPs 222', 'Reference 239'))
  
  # .....time series 2012 to 2015 ----
  windows(width = 14, height = 7)
  ggplot(NPexcr.ts %>% filter(Year != '2014'), 
         aes(x = Year, y = massnorm.N.excr.av, color = Lake,
                        group = Lake)) +
    geom_point(size = 5) +
    geom_line(size = 2) +
    theme_classic(base_size = 26) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face="bold"),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') +
    scale_fill_brewer(palette = "Accent") +
    scale_y_continuous(name = 'TDN excretion (μg N/g/h)') +
    scale_colour_discrete(name = 'Lake',
                          labels = c('AgNPs 222', 'Reference 239'))
  
  # .....time series 2012 to 2015 using model ----
  windows(width = 14, height = 7)
  ggplot(lm.Nxnorm_marginal_means, 
         aes(x = Year, y = emmean, color = Lake,
             group = Lake)) +
    geom_point(size = 5, position = position_dodge(0.5)) +
    geom_errorbar(aes(ymax = lower.CL, ymin = upper.CL),
                  width = 0.2, lwd = 1.5, position = position_dodge(0.5)) +
    geom_jitter(data = NPexcr %>% filter(Year != '2014'),
                aes(x = Year, y = log(massnorm.N.excr)),
                size = 3, position = position_jitterdodge(jitter.width = 0.2),
                show.legend = FALSE, alpha = 0.5) +
    theme_classic(base_size = 26) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face="bold"),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') +
    scale_colour_manual(name = 'Lake',
                      labels = c('AgNPs 222', 'Reference 239'),
                      values = c("black","gray60")) +
    scale_y_continuous(name = 'TDN excretion (μg N/g/h)') 
  
  # ..P excretion anova + figure ----
  # post-addition N excretion in 2015 vs. Lake
  aov.Pnx <- aov(massnorm.P.excr ~ Lake, 
               data = NPexcr  %>% filter(Year != '2012'))
  summary(aov.Pnx)
  aov.Px <- aov(massnorm.P.excr ~ Lake*Year, 
               data = NPexcr)
  summary(aov.Px)
  # where is the difference coming from?
  TukeyHSD(aov.Px, conf.level=0.95, which = 'Lake:Year') 
  # (pre in 2012 vs. post-addition in 2015 N excretion)*/+(Lake)
  # log-transformed because residuals not normally distributed
  
  # testing linear model with log-transformed data
  lm.Pxnorm <- lm(log(massnorm.P.excr) ~ Lake*Year, 
                  data = NPexcr)
  Anova(lm.Pxnorm)
  
  lm.Pxnorm_emmeans <- emmeans(lm.Pxnorm, ~ Lake:Year, type = 'response')
  pairs(lm.Pxnorm_emmeans)
  emmip(lm.Pxnorm, Lake ~ Year)
  lm.Pxnorm_emmeans
  plot(lm.Pxnorm_emmeans, comparison = T)
  
  # extracting effects
  lm.Pxnorm_marginal_means <- as.data.frame(summary(lm.Pxnorm_emmeans))
  
  #glm
  glm.Px <- glm(massnorm.P.excr ~ Lake*Year, 
                data = NPexcr)
  summary(glm.Px)
  glm.Pxnorm <- glm(massnorm.P.excr ~ Lake*Year, 
                    data = NPexcr, 
                    family = Gamma(link = "log")) 
  summary(glm.Pxnorm)
  Anova(glm.Pxnorm)
  
  #extract effects
  glm.Pxnorm_marginal_means <- as.data.frame(Effect(c("Lake", "Year"),
                                                    glm.Pxnorm))
  
  glm.Pxnorm_emmeans <- emmeans(glm.Pxnorm, ~ Lake|Year, type = 'response')
  glm.Pxnorm_emmeans
  windows(width = 14, height = 7)
  plot(glm.Pxnorm_emmeans)
  pairs(glm.Pxnorm_emmeans)
  
  # Generate Levene's deviate 
  NPexcr$LevDevP <- with(NPexcr, 
                         LeveneDeviates(massnorm.P.excr,
                                        group = interaction(Lake, Year, drop = TRUE),
                                        med  = TRUE, log_trans = TRUE))
  # add 10% of the smallest value
  NPexcr %>% summarize(min(NPexcr$LevDevP, na.rm = TRUE)) 
  # smallest value = 0.0008525152
  NPexcr$LevDevP2 <- NPexcr$LevDevP + 0.00009
  glm.LevDevP_gamma <- glm(LevDevP2 ~ Lake*Year,
                          data = NPexcr, 
                          family = "Gamma")
  
  summary(glm.LevDevP_gamma)
  Anova(glm.LevDevP_gamma)
  glm.LevDevP_emmeans <- emmeans(glm.LevDevP_gamma, ~ Lake|Year, type = 'response')
  plot(glm.LevDevP_emmeans)
  
  # testing with linear model
  lm.LevDevP <- lm(LevDevP ~ Lake*Year,
                   data = NPexcr)
  summary(lm.LevDevP)
  Anova(lm.LevDevP)
  
  # ...data visualization ----
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
  
  # .....time series 2012 to 2015 ----
  windows(width = 14, height = 7)
  ggplot(NPexcr.ts, aes(x = Year, y = massnorm.P.excr.av, color = Lake,
                        group = Lake)) +
    geom_point(size = 5) +
    geom_line(size = 2) +
    theme_classic(base_size = 26) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face="bold"),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') +
    scale_colour_brewer(palette = "Accent") +
    scale_y_continuous(name = 'TDP excretion (μg P/g/h)') +
    scale_colour_discrete(name = 'Lake',
                          labels = c('AgNPs 222', 'Reference 239')) +
    guides(colour = guide_legend(override.aes = list(size = 4))) 
  
  # .....time series 2012 to 2015 using model ----
  windows(width = 14, height = 7)
  ggplot(lm.Pxnorm_marginal_means, 
         aes(x = Year, y = emmean, color = Lake,
             group = Lake)) +
    geom_point(size = 5, position = position_dodge(0.5)) +
    geom_errorbar(aes(ymax = lower.CL, ymin = upper.CL),
                  width = 0.2, lwd = 1.5, position = position_dodge(0.5)) +
    geom_jitter(data = NPexcr,
                aes(x = Year, y = log(massnorm.P.excr)),
                size = 2, position = position_jitterdodge(jitter.width = 0.2),
                show.legend = FALSE, alpha = 0.5) +
    theme_classic(base_size = 26) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face="bold"),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') +
    scale_y_continuous(name = 'TDP excretion (μg P/g/h)') +
    scale_colour_manual(name = 'Lake',
                        labels = c('AgNPs 222', 'Reference 239'),
                        values = c("black","gray60"))
  
  # ..C excretion anova + figure ----
  aov.C <- aov(massnorm.C.excr ~ Lake, 
               data = NPexcr)
  summary(aov.C)
  
  glm.Cxnorm <- glm(massnorm.C.excr ~ Lake, 
               data = NPexcr)
  summary(glm.Cxnorm)
  glm.Cxnorm_marginal_means <- as.data.frame(Effect(c("Lake"),
                                                    glm.Cxnorm))
  

  # ...data visualization ----
  windows(width = 14, height = 7)
  ggplot(NPexcr, aes(x = Lake, y = massnorm.C.excr, fill = Lake)) +
    theme_classic(base_size = 26) +
    geom_boxplot(alpha = 0.7, outlier.size = 3, outlier.colour = NULL) +
    #geom_jitter() +
    theme_bw(base_size = 20) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face = "bold"),
          legend.position = 'none') +
    scale_fill_brewer(palette = "Accent") +
    scale_y_continuous(name = 'DOC excretion (mg C/g/h)') +
    scale_x_discrete(labels = c('AgNPs 222', 'Reference 239'))
  
  # .....time series using model ----
  windows(width = 14, height = 7)
  ggplot(glm.Cxnorm_marginal_means, 
         aes(x = Lake, y = fit, color = Lake)) +
    geom_point(size = 5, position = position_dodge(0.5)) +
    geom_errorbar(aes(ymax = lower, ymin = upper),
                  width = 0.2, lwd = 1.5, position = position_dodge(0.5)) +
    geom_line(size = 2, position = position_dodge((0.5))) +
    geom_jitter(data = NPexcr %>% filter(Year != '2014'),
                aes(x = Lake, y = massnorm.C.excr),
                size = 2, position = position_jitterdodge(jitter.width = 0.2),
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
    scale_y_continuous(name = 'DOC excretion (mg C/g/h)') +
    scale_colour_discrete(name = 'Lake',
                          labels = c('AgNPs 222', 'Reference 239'))
  
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
    theme_classic(base_size = 26) +
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
  lm.NPxnorm_marginal_means <- as.data.frame(Effect(c("Lake", "Year"),
                                                     lm.NPxnorm))
  lm.NPxnorm_marginal_means <- as.data.frame(summary(lm.NPxnorm_emmeans))

  #GLM 
  glm.NPxnorm <- glm(massnorm.NP.excr ~ Lake*Year, 
                    data = NPexcr %>% filter(Year != '2014'), 
                    family = Gamma(link = "log"))
  summary(glm.NPxnorm)
  Anova(glm.Nxnorm)
  
  # extracting effects
  glm.NPxnorm_marginal_means <- as.data.frame(Effect(c("Lake", "Year"),
                                                    glm.NPxnorm))
  
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
  
  # .....time series 2012 to 2015 ----
  windows(width = 14, height = 7)
  ggplot(NPexcr.ts %>%  filter(Year !='2014'), 
         aes(x = Year, y = massnorm.NP.excr.av, color = Lake,
                        group = Lake)) +
    geom_point(size = 5) +
    geom_line(size = 2) +
    theme_classic(base_size = 26) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face="bold"),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'), 
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') +
    scale_y_continuous(name = 'Excreted N:P (molar)') +
    scale_colour_discrete(name = 'Lake',
                          labels = c('AgNPs 222', 'Reference 239'),
                          )
  
  # .....time series 2012 to 2015 using model ----
  windows(width = 14, height = 7)
  ggplot(lm.NPxnorm_marginal_means, 
         aes(x = Year, y = emmean, color = Lake,
             group = Lake)) +
    geom_point(size = 5, position = position_dodge(0.5)) +
    geom_errorbar(aes(ymax = lower.CL, ymin = upper.CL),
                  width = 0.2, lwd = 1.5, position = position_dodge(0.5)) +
    geom_jitter(data = NPexcr %>% filter(Year != '2014'),
                aes(x = Year, y = log(massnorm.NP.excr)),
                size = 2, position = position_jitterdodge(jitter.width = 0.2),
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
                          values = c("black","gray60"))
  
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

  
             