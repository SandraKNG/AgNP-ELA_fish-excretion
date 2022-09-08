  
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
           Lake = as.factor(Lake),
           Year = as.factor(Year))
  
  str(NPexcr)
  
  
  ############################### TEST #######################################
  
  # data visualization 
  # mass vs N or P excretion for all the data
  windows(width = 14, height = 7)
  ggplot(NPexcr, aes(x = Log.mass, y = Log.N.excretion)) +
           geom_point() 
  #very bad relationship, almost flat for P excretion +
  #not convinced it is actually linear for N excretion
  #probably due to the low variance in mass among individuals (from 0.4 to 8g)
  
  lm(Log.P.excretion ~ Log.mass, data = NPexcr)
  
  ################################ ASSUMPTIONS CHECK ########################
  
  # ..normality and variance in 2014 ----
  histkdnc(NPexcr$N.excretion) # looks like Poisson distribution
  boxplot(N.excretion ~ Lake, 
          data = NPexcr %>% filter(Year == '2014'))
  leveneTest(N.excretion ~ Lake, data = NPexcr %>% filter(Year == '2014')) 
  
  histkdnc(NPexcr$P.excretion) # normal distribution
  boxplot(P.excretion ~ Lake, 
          data = NPexcr %>% filter(Year == '2014'))
  leveneTest(P.excretion ~ Lake, data = NPexcr %>% filter(Year == '2014'))
  # check for equality of variances --> they are equal for both N and P excretions!
  
  # ..normality and variance in 2015 ----
  histkdnc(NPexcr$N.excretion) # Poisson distribution
  boxplot(N.excretion ~ Lake, 
          data = NPexcr %>% filter(Year == '2015'))
  leveneTest(N.excretion ~ Lake, data = NPexcr %>% filter(Year == '2015')) 
  
  histkdnc(NPexcr$P.excretion %>% filter(Year == '2014'))
  boxplot(P.excretion ~ Lake, 
          data = NPexcr %>% filter(Year == '2015'))
  leveneTest(P.excretion ~ Lake, data = NPexcr %>% filter(Year == '2015'))
  # check for equality of variances --> they are equal for both N and P excretions!
  
  ################################### DATA ANALYSIS ##################################
  
  # ..glm for N excretion ----
  GLM.nx <- glm(N.excretion ~ Lake, data = NPexcr) # GLM without interaction
  GLM.x <- glm(N.excretion ~ Lake+Year+Temperature, 
      data = NPexcr, family = Poisson)
  # GLM with interaction
  
  interaction.plot(x.factor = Year, 
                   trace.factor = Lake, 
                   response = N.excretion, data = NPexcr)

    
             