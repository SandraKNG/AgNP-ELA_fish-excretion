  # *** Whole-lake silver nanoparticles addition promotes phosphorus and 
  # silver excretion by yellow perch (Perca flavescens)  ***
  # this code was developped by S. Klemet-N'Guessan in 2020-2023
  
  # load libraries ----
  library(car)
  library(lmerTest) # for lmer
  
  # ..log10 mass vs log10 N or P excretion for all the data by lake ----
  # visualize Log10 N excretion vs. Log10 mass by Lake
  windows(width = 14, height = 7)
  ggplot(NPexcr, aes(x = Log10.mass, y = Log10.N.excretion, col = Year)) +
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
  # probably due to the high variance in mass among individuals (from 0.4 to 8g) 
  # and small fish size 
  
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
  # plot mean P/N excretion vs. SD to see if there is a positive relationship between the two. 
  # If yes, a log-transformation is possible. Here both mean P and N excretions have a good positive relationship with their respective SD
  ggplot(NPexcr.ts, aes(massnorm.P.excr.sd, massnorm.P.excr.av)) +
    geom_point()
  ggplot(NPexcr.ts, aes(massnorm.N.excr.sd, massnorm.N.excr.av)) +
    geom_point()
  # test for heterogeneity of variances using Levene's test BUT old-school, not really telling where the heterogeneity is
  leveneTest(log(massnorm.N.excr) ~ Lake*Year, data = NPexcr)
  leveneTest(log(massnorm.N.excr) ~ Lake, data = NPexcr)
  # variances significantly diff
  leveneTest(log(massnorm.P.excr) ~ Lake*Year, data = NPexcr) 
  leveneTest(log(massnorm.NP.excr) ~ Lake* Year, data = NPexcr)
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
  
  # Ag excretion
  # Ag
  gghistogram(excr.Tag, x = 'massnorm.Tag.excr', rug = T, 
              fill = 'Year', add = 'mean')
  ggqqplot(excr.Tag, x = "massnorm.Tag.excr", facet.by = "Year")
  excr.Tag %>%  group_by(Year) %>% shapiro_test(massnorm.Tag.excr)
  
  # N:Ag
  hist(NPexcr$massnorm.NAg.excr)
  
  # P:Ag
  hist(NPexcr$massnorm.PAg.excr)
  
  # Alternative models tests: ANCOVA and repeated measured ANOVA ----
  # testing whether using N.excretion, log or log10 transformed N excretion, 
  # residuals vs. predicted values still look quite different from 
  # those generated randomly: light/heavy-tailed on both ends (N excretion)
  # or on the right only (log and log10 transformation)
  
  # N excretion
  # ANCOVA
  lm.Nx <- lm(log10(N.excretion.rate) ~ log10(Mass)*Year, 
              data = NPexcr %>% filter(Year != '2014'))
  lm.Nx <- lm(log10(N.excretion.rate) ~ log10(Mass)*Lake, 
              data = NPexcr %>% filter(Year != '2014'))
  Anova(lm.Nx, type = 'III')
  # which levels are different? Lakes from 2012 vs. 2015
  TukeyHSD(aov(lm.Nxnorm), ordered = F)
  
  # repeated measures ANOVA
  # using site class (i.e., treatment and control) and 
  # period (i.e., before, during, and after) as fixed factors with an interaction effect, 
  # and lake and year as random factors 
  lmN <- lmer(log(massnorm.N.excr) ~ Period * SiteClass + (1 |Lake:Year), data = NPexcr %>% filter(Year != '2014'))
  anova(lmN)
  anova(lmN, ddf = "Kenward-Roger")
  summary(lmN)
  
  # P excretion
  # ANCOVA
  lm.Px <- lm(log10(P.excretion.rate) ~ log10(Mass)*Year, 
              data = NPexcr)
  lm.Px <- lm(log10(P.excretion.rate) ~ log10(Mass)*Lake, 
              data = NPexcr)
  Anova(lm.Px, type = 'III')
  summary(lm.Px)
  # which levels are different? Lakes from 2014 vs. 2015
  postHocs <- glht(lm.Px, linfct = mcp(Year = "Tukey"))
  
  # repeated measures ANOVA
  lmP <- lmer(log(massnorm.P.excr) ~ Period * SiteClass + (1 |Lake:Year), data = NPexcr)
  anova(lmP)
  anova(lmP, ddf = "Kenward-Roger")
  summary(lmP)
  
  # P excretion in 2022 only
  NPexcr_22 <- NPexcr_22 %>% mutate(massnorm.P.excr = P.excretion.rate/Mass)
  ggplot(NPexcr_22, aes(x = Lake, y = massnorm.P.excr,
                        group = Lake)) +
    geom_boxplot() +
    geom_point(aes(colour = Lake))
  
  # correlation between N, P, Ag excretion rates
  ggplot(NPexcr %>% filter(Year == 2015), 
         aes(x = massnorm.N.excr, y = massnorm.P.excr)) +
    geom_point()
  
  ggplot(NPexcr %>% filter(Year == 2015), 
         aes(x = massnorm.N.excr, y = massnorm.Tag.excr)) +
    geom_point()
  
  ggplot(NPexcr %>% filter(Year == 2015), 
         aes(x = massnorm.P.excr, y = massnorm.Tag.excr)) +
    geom_point()
  
  
  # Fish Stoich model ----
  # set up parameters
  parameters <- model_parameters("Perca flavescens", otolith = F,
                                 family = "Percidae", temp = 18, mirror = "se")
  # find length-weight relationship and caudal fin aspect ratio
  find_lw("Perca flavescens", mirror = "se")
  aspect_ratio("Perca flavescens")
  # find growth parameters
  growth_params <- growth_params("Perca flavescens", otolith = FALSE)
  # find trophic level
  # trophic_level("Perca flavescens")