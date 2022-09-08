  
  # load libraries and read silver nanoparticles (NP) dataset ----
  
  library(tidyverse)
  library(nortest)
  library(descr)
  library(car)

  NPer <- read.csv('20 06 02 AgNP ELA lakes fish excretion 2014_2015.csv',
                 stringsAsFactors = F, na.strings = c("", "NA", "."), 
                 strip.white = TRUE, sep = ",")  

  str(NPer)
  head(NPer)
  
  # clean, rename, and add variables to the dataset
  NPexcr <- NPer %>% 
    rename(Year = Sampling.year,
           N.excretion = N.excretion.rate..Âµg.ind.h.,
           P.excretion = P.excretion.rate..Âµg.ind.h.,
           Mass = Mass..g.,
           Temperature = Temperature.Â.) %>% 
    mutate(Log.mass = log(Mass),
           Log.N.excretion = log(N.excretion),
           Log.P.excretion = log(P.excretion),
           Log10.N.excretion = log10(N.excretion),
           Lake = as.factor(Lake),
           Year = as.factor(Year))
  
  str(NPexcr)
  
  
  ############################### TEST #######################################
  
  # data visualization 
  # mass vs N or P excretion for all the data
  windows(width = 14, height = 7)
  ggplot(NPexcr, aes(x = Log.mass, y = N.excretion)) +
           geom_point() 
  #very bad relationship, almost flat for P excretion +
  #not convinced it is actually linear for N excretion
  #probably due to the low variance in mass among individuals (from 0.4 to 8g)
  
  ggplot(NPexcr, aes(x = Log.mass, y = Log.N.excretion)) +
           geom_point() +
    facet_wrap(~Year)
  
  lm(Log.P.excretion ~ Log.mass, data = NPexcr)
  
  ################################ ASSUMPTIONS CHECK ########################
  
  # ..normality and variance of the data in both 2014 and 2015
  histkdnc(NPexcr$Log.N.excretion) # Poisson distribution
  hist(NPexcr$Log.N.excretion)
  boxplot(N.excretion ~ Lake, 
          data = NPexcr)
  
  # ..normality and variance of the data in 2014 ----
  NPexcr.2014 <- NPexcr %>% filter(Year == '2014')
  histkdnc(NPexcr.2014$N.excretion) # looks like Poisson distribution
  boxplot(N.excretion ~ Lake, 
          data = NPexcr %>% filter(Year == '2014'))
  leveneTest(N.excretion ~ Lake, data = NPexcr %>% filter(Year == '2014')) 
  
  histkdnc(NPexcr.2014$P.excretion) # normal distribution
  boxplot(P.excretion ~ Lake, 
          data = NPexcr %>% filter(Year == '2014'))
  leveneTest(P.excretion ~ Lake, data = NPexcr %>% filter(Year == '2014'))
  # check for equality of variances --> they are equal for both N and P excretions!
  
  # ..normality and variance of the data in 2015 ----
  NPexcr.2015 <- NPexcr %>% filter(Year == '2015')
  histkdnc(NPexcr.2015$N.excretion) # Poisson distribution
  boxplot(N.excretion ~ Lake, 
          data = NPexcr %>% filter(Year == '2015'))
  leveneTest(N.excretion ~ Lake, data = NPexcr %>% filter(Year == '2015')) 
  
  histkdnc(NPexcr.2015$P.excretion)
  boxplot(P.excretion ~ Lake, 
          data = NPexcr %>% filter(Year == '2015'))
  leveneTest(P.excretion ~ Lake, data = NPexcr %>% filter(Year == '2015'))
  # check for equality of variances --> they are equal for both N and P excretions!
  
  # ..normality of the residuals ----
  qqPlot(lm(N.excretion ~ Lake+Temperature, 
            data = NPexcr))
  
  # plot standardized residuals of the lm against predicted values 
  # and make informed decision by comparing it to plotted residuals 
  # of randomly generated datasets
  z = lm(Log10.N.excretion ~ Lake+Temperature, data = NPexcr)$residual
  n = length(z)
  xz = cbind(matrix(rnorm(12*n), nr = n),z,matrix(rnorm(12*n), nr = n))
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
  # those generated randomly ()
  
  
  ################################### DATA ANALYSIS ##################################
  
  # ..glm for N excretion ----
  GLM.nx <- glm(N.excretion ~ Lake, 
                data = NPexcr, family = lognormal) # GLM without interaction
  GLM.x <- glm(N.excretion ~ Lake+Temperature, 
      data = NPexcr, family = quasipoisson)
  aov.nx <- aov(N.excretion ~ Lake, 
      data = NPexcr)
  aov.nx <- aov(N.excretion ~ Lake+Temperature, 
                data = NPexcr) #not lookign at year effect because very low sample sizes in each (5-7)
  summary(GLM.nx)
  summary(GLM.x)
  summary(aov.nx)
  Anova(GLM.nx)
  anova()
  # GLM with interaction
  
  interaction.plot(x.factor = Year, 
                   trace.factor = Lake, 
                   response = N.excretion, data = NPexcr)

    
             