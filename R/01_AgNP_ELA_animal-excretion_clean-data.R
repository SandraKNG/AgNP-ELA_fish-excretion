  # *** Whole-lake silver nanoparticles addition promotes phosphorus and 
  # silver excretion by yellow perch (Perca flavescens)  ***
  # this code was developped by S. Klemet-N'Guessan in 2020-2023
  
  # load libraries ----
  library(tidyverse)
  library(datawizard) # to do summary statistics
  
  # read silver nanoparticles (NP) dataset ----
  NPer <- read_csv('data/2020-04-21_AgNP-ELA-lakes_fish-excretion.csv')
  param <- read_csv('data/2021-09-20_param_modelsp_FishStoich.csv')
  NPer_22 <- read_csv('data/2022-07-11_11-DOC-lakes_Mastersheet_YP.csv')
  
  str(NPer)
  head(NPer)
  str(param)
  parammod <- as.list(param)
  
  # clean, rename, and add variables to the dataset ----
  NPexcr <- NPer %>% 
    rename(Year = `Sampling year`,
           N.excretion.rate = `N excretion rate (µg/ind/h)`,
           P.excretion.rate = `P excretion rate (µg/ind/h)`,
           C.excretion.rate = `C excretion (µg C/ind/h)`,
           Tag.excretion.rate = `Tag excretion rate (µg/ind/h)`,
           Excreted.CN = `Excreted C:N (molar)`,
           Excreted.NP = `Excreted N:P (molar)`,
           Excreted.CP = `Excreted C:P (molar)`,
           Mass = `Dry mass (g)`,
           Wet.mass = `Mass (g)`,
           Temperature = `Temperature ©`,
           FishID = 'Fish ID') %>% 
    mutate(
      Total.length.cm = `Total length (mm)`/10,
      Fork.length = `Fork length (mm)`,
      Lake = as.factor(Lake),
      Year = as.factor(Year),
      Period = factor(ifelse(Year == '2012', 'Before', 'After')),
      SiteClass = factor(ifelse(Lake == '239','Control', 'Impact')),
      Period = fct_relevel(Period, 'After', after = Inf)) 
  
  str(NPexcr)
  
  NPexcr_22 <- NPer_22 %>% 
    select(ID, Site.name, Species.code, N.excretion.rate, P.excretion.rate, 
           C.excretion.rate, Dry.mass, Wet.mass) %>% 
    rename(
      Lake = Site.name,
      Mass = Dry.mass
    ) %>% 
    filter(
      !(Species.code %in% c('CTL1', 'CTL2'))
    ) %>% 
    dplyr::mutate(
      Year = factor(2022),
      Lake = if_else(Lake == 'L239', '239', '222')
    ) %>% 
    dplyr::filter(Mass <= 2.25) %>% 
    group_by(Lake) %>% 
    filter(case_when(
      Lake == '239' ~ Mass < .5,
      TRUE ~ TRUE
    ))
  
  # join two datasets
  NPexcr <- bind_rows(NPexcr, NPexcr_22)
  
  # mass-correct the excretion ----
  NPexcr <- NPexcr %>% 
    mutate(massnorm.N.excr = N.excretion.rate/Mass,
           massnorm.P.excr = P.excretion.rate/Mass,
           massnorm.C.excr = C.excretion.rate/Mass,
           massnorm.Tag.excr = Tag.excretion.rate/Mass,
           massnorm.NP.excr = (massnorm.N.excr/massnorm.P.excr)/(14/31),
           massnorm.CN.excr = (massnorm.C.excr/massnorm.N.excr)/(12/14),
           massnorm.CP.excr = (massnorm.C.excr/massnorm.P.excr)/(12/31),
           massnorm.NP.excr2 = Excreted.NP/Mass,
           massnorm.NAg.excr = massnorm.N.excr/massnorm.Tag.excr/(14/107),
           massnorm.PAg.excr = massnorm.P.excr/massnorm.Tag.excr/(12/107))
  
  # ..summary statistics ----
  
  NPexcr.ss1 <- NPexcr %>% group_by(Lake, Year) %>%
    select(c('Mass', 'Fork.length', 'massnorm.N.excr', 'massnorm.P.excr',
             'massnorm.Tag.excr')) %>% 
    describe_distribution() 
  NPexcr.ss1
  
  NPexcr.ss2 <- NPexcr %>% 
    select(c('Mass', 'massnorm.N.excr', 'massnorm.P.excr', 'massnorm.Tag.excr')) %>% 
    describe_distribution() 
  NPexcr.ss2
  