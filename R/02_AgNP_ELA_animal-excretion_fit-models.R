  # *** Silver nanoparticles addition promotes phosphorus and silver excretion by 
  # yellow perch (Perca flavescens) in a boreal lake ***
  # this code was developped by S. Klemet-N'Guessan in 2020 and 2021
  
  # load libraries ----
  library(emmeans) # post-hoc pairwise comparisons using lm, glm, levDev models
  library(car)
  
  # LM & emmeans ----
  emm <- function(lm){
    # pairwise comparisons using emmeans
    m1 <- emmeans(lm, pairwise ~ Lake:Year, type = 'response')
    m2 <- emmeans(lm, pairwise ~ Lake|Year, type = 'response')
    contrasts <- m1$contrasts %>% summary(infer = T) %>% 
      as_tibble() %>%  arrange(p.value)
    print(contrasts, n = 40)
    
    # Extracting effects from emmeans
    emm_df <- as_tibble(m1$emmeans)
    
    # Calculating effect sizes + extracting effects
    eff_size_df <- summary(eff_size(m2, sigma = sigma(lm), 
                                    edf = df.residual(lm))) %>% 
      as_tibble()
    print(eff_size_df)
    
    # Return a list containing the results
    result_list <- list(contrasts = contrasts, 
                        emmeans = emm_df, effect_sizes = eff_size_df)
    return(result_list)
  }
  
  # N excretion
  lmN <- lm(log10(massnorm.N.excr) ~ Lake*Year, 
            data = NPexcr %>% filter(Year != '2014'))
  anova(lmN)
  pwcN <- emm(lmN)
  
  # P excretion
  lmP <- lm(log10(massnorm.P.excr) ~ Lake*Year, 
            data = NPexcr)
  anova(lmP)
  pwcP <- emm(lmP)
  
  # P excretion
  lmNP <- lm(log10(massnorm.NP.excr) ~ Lake*Year, 
            data = NPexcr %>% filter(Year != '2014'))
  anova(lmNP)
  pwcNP <- emm(lmNP)
  
  # make contrast table based on all models
  contrasts <- bind_rows(pwcN[['contrasts']], pwcP[['contrasts']], 
                         pwcNP[['contrasts']], .id = 'column_label')
  contrasts <- contrasts %>% 
    mutate(test = if_else(column_label == 1, 'massnorm.N.excr',
                          if_else(column_label == 2, 'massnorm.P.excr', 
                                  'massnorm.NP.excr')))
  
  write.csv(contrasts, 'output/emmeans_contrasts.csv')
  
  # wilcox test ----
  excr.Tag <- NPexcr %>% 
    dplyr::filter(!(Year %in% c('2012', '2022')))
  # Ag
  hist(NPexcr$massnorm.Tag.excr)
  gghistogram(NPexcr, x = 'massnorm.Tag.excr', rug = T, 
              fill = 'Year', add = 'mean')
  var.test(massnorm.Tag.excr ~ Year,  excr.Tag, alternative = 'two.sided')
  w.Tag <- wilcox.test(massnorm.Tag.excr ~ Year,  excr.Tag)
  w.Tag
  
  # N:Ag
  hist(NPexcr$massnorm.NAg.excr)
  var.test(massnorm.NAg.excr ~ Year, NPexcr, alternative = 'two.sided')
  w.NAg <- wilcox.test(massnorm.NAg.excr ~ Year,  excr.Tag)
  w.NAg
  
  # P:Ag
  hist(NPexcr$massnorm.PAg.excr)
  var.test(massnorm.PAg.excr ~ Year, NPexcr, alternative = 'two.sided')
  w.PAg <- wilcox.test(massnorm.PAg.excr ~ Year,  excr.Tag)
  w.PAg
  
  # simulation models ----
  # ..V&M universal model ----
  set.seed(1)
  VM.m <- NPexcr %>%  select(Lake, Year, Temperature, Mass)
  VM.m <- VM.m %>% group_by(Lake, Year) %>% 
    mutate(Temperature = replace_na(Temperature, sample(16:20, 1)),
           N.excretion.m = 10^(1.461 + 0.684 * log10(Mass) + 
                                 0.0246 * Temperature - 0.2013 + 0.7804),
           N.excretion.m.se = 10^(0.0897 + 0.0177 * log10(Mass) + 
                                    0.0014 * Temperature - 0.0771 + 0.0655),
           P.excretion.m = 10^(0.6757 + 0.5656 * log10(Mass) + 
                                 0.0194 * Temperature - 0.248 + 0.7504),
           P.excretion.m.se = 10^(0.0992 + 0.0205 * log10(Mass) +
                                    0.002 * Temperature - 0.0922 + 0.0768),
    ) %>% 
    ungroup() %>% 
    dplyr::filter(!is.na(Mass))
  
  # Select columns and repeat values
  VM.m_sub <- VM.m %>% 
    select(starts_with(c('N', 'P', 'M'))) 
  
  iter_yp_VM <- data.frame(
    N.excretion.m = rep(VM.m_sub$N.excretion.m, each = 2500),
    N.excretion.m.se = rep(VM.m_sub$N.excretion.m.se, each = 2500),
    P.excretion.m = rep(VM.m_sub$P.excretion.m, each = 2500),
    P.excretion.m.se = rep(VM.m_sub$P.excretion.m.se, each = 2500),
    Mass = rep(VM.m_sub$Mass, each = 2500),
    iter = rep(1:2500, each = 131)
  ) 
  
  # ..Schiettekatte FishStoich model ----
  # set up parameters
  # parameters <- model_parameters("Perca flavescens", otolith = F,
  #                                family = "Percidae", temp = 18, mirror = "se")
  # # find length-weight relationship and caudal fin aspect ratio
  # find_lw("Perca flavescens", mirror = "se")
  # aspect_ratio("Perca flavescens")
  # # find growth parameters
  # growth_params <- growth_params("Perca flavescens", otolith = FALSE)
  # # find trophic level
  # trophic_level("Perca flavescens")
  
  # ..Run model ----
  # Yellow perch
  total.length <- exp((log(NPexcr$Wet.mass[!is.na(NPexcr$Wet.mass)]/
                             parammod$lwa_m)/parammod$lwb_m))
  wet.mass <- NPexcr$Wet.mass[!is.na(NPexcr$Wet.mass)]
  FStoichm <- fishflux::cnp_model_mcmc(TL = total.length, param = parammod, 
                                       iter = 5000)
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
    lapply( FUN = get_iter) %>%
    dplyr::bind_rows()
  
  # get dataset ready for excretion simulation models
  iter_yp <- iter_yp %>% 
    mutate(
      tl = rep(total.length, each = 2500),
      wet.mass = rep(NPexcr$Wet.mass[!is.na(NPexcr$Wet.mass)], each = 2500),
      dry.mass = rep(NPexcr$Mass[!is.na(NPexcr$Mass)], each = 2500),
      N.excretion.mod = Fn * 10 ^ 6 / 24,
      P.excretion.mod = Fp * 10 ^ 6 / 24,
      C.ingestion.mod = Ic * 10 ^ 6 / 24
    )
  
  # Predict P excretion increase with diet increase ----
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
  
  
  # NPC limitation
  limitation <- limitation(FStoichm, plot = F)
  a <- parammod$lwa_m
  b <- parammod$lwb_m
  limitation <- limitation %>% mutate(mass = a*tl^b,
                                      dry.mass = mass*0.25)
  
  
  