  
  # load libraries and read silver nanoparticles (NP) dataset ----
  
  library(tidyverse)
  library(nortest)
  library(descr)
  library(car)
  library(RColorBrewer)
  library(moments)
  library(MASS)

  NPer <- read.csv('20 06 18 AgNP ELA lakes fish excretion.csv',
                 stringsAsFactors = F, na.strings = c("", "NA", "."), 
                 strip.white = TRUE, sep = ",")  

  str(NPer)
  head(NPer)
  
  # clean, rename, and add variables to the dataset
  NPexcr <- NPer %>% 
    rename(Year = Sampling.year,
           N.excretion = N.excretion.rate..Âµg.ind.h.,
           P.excretion = P.excretion.rate..Âµg.ind.h.,
           C.excretion = C.excretion..mg.C.ind.h.,
           Tag.excretion = Tag.excretion.rate..Âµg.ind.h.,
           Excreted.CN = Excreted.C.N..molar.,
           Excreted.NP = Excreted.N.P..molar.,
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
  # unique(NPexcr$Year)
  # tapply(NPexcr$N.excretion, INDEX=NPexcr$Year, FUN=length, na.rm = T) 
  
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
  NPmasscorr.combined <- data.frame(Lake = c('222', '239', '222', '239'),
                                    b.coeff.N.excr = c(b.N222, b.N239),
                                    b.coeff.P.excr = c(b.P222, b.P239))
  
  # NPmasscorr.combined <- bind_rows(NPmasscorr.222, NPmasscorr.239)
  NPexcr <- left_join(NPexcr, NPmasscorr.combined, by = "Lake") %>% 
    distinct() %>% 
    mutate(massnorm.N.excr = N.excretion/(Mass^(b.coeff.N.excr)),
           massnorm.P.excr = P.excretion/(Mass^(b.coeff.P.excr)),
           Log.massnorm.N.excr = massnorm.N.excr,
           Log.massnorm.P.excr = massnorm.P.excr)

  
  ################################ ASSUMPTIONS CHECK ########################
  
  # ..normality and variance of the data in both 2014 and 2015 ----
  histkdnc(NPexcr$massnorm.N.excr) # normal distribution
  histkdnc(NPexcr$P.excretion) # normal distribution
  shapiro.test(NPexcr$Log.massnorm.N.excr)
  boxplot(N.excretion ~ Lake, 
          data = NPexcr)
  leveneTest(Log.N.excretion ~ Lake, data = NPexcr) # variances not significantly diff
  leveneTest(P.excretion ~ Lake, data = NPexcr) # variances not significantly diff
  leveneTest(C.excretion ~ Lake, data = NPexcr) # variances not significantly diff
  
  # ..normality of the residuals 
  windows(width = 14, height = 7)
  qqPlot () # Log N excretion turned out to be better than N excretion
  skewness(NPexcr$N.excretion, na.rm = T) 
   
  # skewness = 2.9 > 0.5 threshold (Webster and Oliver, 2007)
  qqPlot(lm(P.excretion ~ Lake*Year, 
  skewness(NPexcr$P.excretion, na.rm = T)))
  # skewness = 1.5 > 0.5 threshold (Webster and Oliver, 2007)
  qqPlot(lm(C.excretion ~ Lake, 
            data = NPexcr))
  qqPlot(lm(Tag.excretion ~ Year, 
            data = NPexcr))
  # for massnormalized excretion
  plot(glm(massnorm.N.excr ~ Lake*Year, 
            data = NPexcr %>% filter(Year != '2014'), 
            family = Gamma(link = "log")))
  
  # ..plot standardized residuals of the lm against predicted values ----
  # and make informed decision by comparing it to plotted residuals 
  # of randomly generated datasets
  # Tag or N excretion
  z = glm(massnorm.N.excr ~ Lake*Year, 
          data = NPexcr %>% filter(Year != '2014'), 
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
  
  # P excretion
  z = lm(P.excretion ~ Lake*Year, data = NPexcr)$residual
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
  
  
  ################################### DATA ANALYSIS ##################################
  
  # ..dataset prep for time series
  NPexcr.ts <- NPexcr %>% group_by(Lake, Year) %>% 
    summarize(P.excretion.av = mean(P.excretion, na.rm = TRUE),
              P.excretion.sd = sd(P.excretion, na.rm = TRUE),
              N.excretion.av = mean(N.excretion, na.rm = TRUE),
              N.excretion.sd = sd(N.excretion, na.rm = TRUE),
              Excreted.NP.av = mean(Excreted.NP, na.rm = TRUE),
              Excreted.NP.sd = sd(Excreted.NP, na.rm = TRUE),
              massnorm.P.excr.av = mean(massnorm.P.excr, na.rm = TRUE),
              massnorm.P.excr.sd = sd(massnorm.P.excr, na.rm = TRUE))
  
  # ..N excretion anova + figure ----
  # post-addition N excretion in 2015 vs. Lake
  aov.Nnx <- aov(Log.N.excretion ~ Lake, 
               data = NPexcr %>% filter(Year == '2015'))
  # not pooling 2014 + 2015 together may be a difference between both
  # not looking at year effect (2014 vs. 2015) 
  # because in 2014, 1 and 6 data points for L239 and L222 respectively 
  aov.Nx <- aov(Log.N.excretion ~ Lake*Year, 
               data = NPexcr %>% filter(Year != '2014')) 
  # (pre in 2012 vs. post-addition in 2015 N excretion)*/+(Lake)
  # excluding 2014 because of low sample sizes
  glm.Nx <- glm(N.excretion ~ Lake*Year*Mass, 
                data = NPexcr %>% filter(Year != '2014'),
                family = Gamma(link = "log"))
  glm.Nxnorm <- glm(massnorm.N.excr ~ Lake*Year, 
                data = NPexcr %>% filter(Year != '2014'), 
                family = Gamma(link = "log")) 
  # use Gamma distribution because that model is a good start for positive continuous response variables.
  # log link is important here as the default inverse link doesn't make any sense
  # You need to think about range of possible values when choosing a distribution. Gaussian to t distributions allow for negative N excretion which makes no physical sense. You want a distribution for positive reals if we can assume none of the observations is censored... 
  # A Gamma distribution or inverse Gaussian are choices in the GLM world for positive real values responses. If 0 is allowed (not here) then Tweedie might work. If you have censored data, y_i < Level of Detection, then you need to fit a censored model, say censored Gamma 
  
  summary(aov.Nnx)
  summary(aov.Nx)
  summary(glm.Nx)
  summary(glm.Nxnorm)
    
  # data visualization   
  windows(width = 14, height = 7)
  ggplot(NPexcr %>% filter(Year != '2012'), 
         aes(x = Lake, y = N.excretion, fill = Lake)) +
    theme_classic(base_size = 26) +
    geom_boxplot(alpha = 0.7, outlier.size = 3, outlier.colour = NULL) +
    #geom_jitter() +
    theme_bw(base_size = 20) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face="bold"),
          legend.position = 'none') +
    scale_fill_brewer(palette = "Accent") +
    scale_y_continuous(name = 'TDN excretion (μg N/ind/h)') +
    scale_x_discrete(labels = c('AgNPs 222', 'Reference 239'))
  
  # ...time series 2012 to 2015 ----
  windows(width = 14, height = 7)
  ggplot(NPexcr.ts %>% filter(Year != '2014'), 
         aes(x = Year, y = N.excretion.av, color = Lake,
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
    scale_y_continuous(name = 'TDN excretion (μg N/ind/h)') +
    scale_colour_discrete(name = 'Lake',
                          labels = c('AgNPs 222', 'Reference 239'))
  

  
  # ..P excretion anova + figure ----
  # post-addition N excretion in 2015 vs. Lake
  aov.Pnx <- aov(P.excretion ~ Lake, 
               data = NPexcr  %>% filter(Year != '2012'))
  aov.Px <- aov(Log.P.excretion ~ Lake*Year, 
               data = NPexcr)
  # (pre in 2012 vs. post-addition in 2015 N excretion)*/+(Lake)
  # log-transformed because residuals not normally distributed
  
  glm.Px <- glm(P.excretion ~ Lake*Year*Mass, 
                data = NPexcr,
                family = Gamma(link = "log"))
  glm.Pxnorm <- glm(massnorm.P.excr ~ Lake*Year*Mass, 
                    data = NPexcr, 
                    family = Gamma(link = "log")) 
  
  TukeyHSD(aov.Px, conf.level=0.95, which = 'Lake:Year') 
  # where is the difference coming from?
  
  summary(aov.Pnx)
  summary(aov.Px)
  summary(glm.Px)
  summary(glm.Pxnorm)
  
  
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
    scale_y_continuous(name = 'TDP excretion (μg P/ind/h)') +
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
  
  # ...time series 2012 to 2015 ----
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
    scale_y_continuous(name = 'TDP excretion (μg P/ind/h)') +
    scale_colour_discrete(name = 'Lake',
                          labels = c('AgNPs 222', 'Reference 239')) +
    guides(colour = guide_legend(override.aes = list(size = 4))) 
  
  # ..C excretion anova + figure ----
  aov.C <- aov(C.excretion ~ Lake, 
               data = NPexcr)
  summary(aov.C)


  windows(width = 14, height = 7)
  ggplot(NPexcr, aes(x = Lake, y = C.excretion, fill = Lake)) +
    theme_classic(base_size = 26) +
    geom_boxplot(alpha = 0.7, outlier.size = 3, outlier.colour = NULL) +
    #geom_jitter() +
    theme_bw(base_size = 20) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face = "bold"),
          legend.position = 'none') +
    scale_fill_brewer(palette = "Accent") +
    scale_y_continuous(name = 'DOC excretion (mg C/ind/h)') +
    scale_x_discrete(labels = c('AgNPs 222', 'Reference 239'))
  
  
  # ..Tag excretion anova + figure ----
  aov.Tag <- aov(Log.Tag.excretion ~ Year, 
               data = NPexcr)
  summary(aov.Tag)
  
  
  windows(width = 14, height = 7)
  ggplot(NPexcr, aes(x = Year, y = Log.Tag.excretion, fill = Year)) +
    theme_classic(base_size = 26) +
    geom_boxplot(alpha = 0.7, outlier.size = 3, outlier.colour = NULL) +
    #geom_jitter() +
    theme_bw(base_size = 20) +
    theme(text = element_text( family = "Tahoma"),
          axis.title = element_text(face = "bold"),
          legend.position = 'none') +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(name = 'Log Tag excretion (μg Tag/ind/h)') +
    scale_x_discrete(labels = c('2014', '2015'))
  
  # ..excreted N:P (molar) ----
  windows(width = 14, height = 7)
  ggplot(NPexcr %>%  filter(Year =='2015'), 
         aes(x = Lake, y = Excreted.NP, fill = Lake)) +
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
  
  # ...time series 2012 to 2015 ----
  windows(width = 14, height = 7)
  ggplot(NPexcr.ts %>%  filter(Year !='2014'), 
         aes(x = Year, y = Excreted.NP.av, color = Lake,
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
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(name = 'Excreted N:P (molar)') +
    scale_colour_discrete(name = 'Lake',
                          labels = c('AgNPs 222', 'Reference 239'))
  
  # ..excreted C:N (molar) ----
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
  
  # ..

  
             