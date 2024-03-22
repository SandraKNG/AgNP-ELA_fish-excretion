  # *** Whole-lake silver nanoparticles addition promotes phosphorus and 
  # silver excretion by yellow perch (Perca flavescens)  ***
  # this code was developped by S. Klemet-N'Guessan in 2020-2023
  
  # load libraries ----
  library(lmerTest) # for lmer
  library(fishflux)
  library(emmeans) # post-hoc pairwise comparisons using lm, glm, levDev models
  library(car)
  library(performance)
  # library(pwr) # for power analysis for MDD
  # library(arm)
  
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
    print(eff_size(m2, sigma = sigma(lm), 
                   edf = df.residual(lm)))
    print(eff_size_df)
    
    # Return a list containing the results
    result_list <- list(contrasts = contrasts, 
                        emmeans = emm_df, effect_sizes = eff_size_df)
    return(result_list)
  }
  
  
  # Mass
  hist(NPexcr$Mass)
  lmMass <- lm(log(Mass) ~ Lake*Year, data = NPexcr)
  check_model(lmMass)
  anova(lmMass)
  summary(lmMass)
  pwcMass <- emm(lmMass)
  
  # N excretion
  # Treatment is the fixed effect of interest.
  # (1 | Lake) specifies random intercepts for each lake, allowing for individual variation in intercepts.
  # (1 | Lake:Year) specifies random intercepts for the interaction between Lake and Year, allowing for variation in intercepts across different combinations of Lake and Year.
  # This model accounts for repeated measures within each lake and allows for individual variation between lakes as well as variation in the effect of treatment across different lakes and years.
  NPexcr.N <- NPexcr %>% filter(Year != '2014')
  log.massnorm.N.excr <- log10(NPexcr$massnorm.N.excr)
  lmN <- lm(log(massnorm.N.excr) ~ Lake * Year, data = NPexcr.N)
  lmN <- lmer(log(massnorm.N.excr) ~ SiteClass * Period + (1 | Lake),
            data = NPexcr.N)
  # lmN <- lmer(log10(massnorm.N.excr) ~ Treatment * Year + (1 + Year | Lake), 
  #             data = NPexcr.N)
  anova_N <- anova(lmN)
  summary(lmN)
  # m2 <- emmeans(lmN, specs = "Treatment", type = 'response')
  # eff_size_df <- summary(eff_size(m2, sigma = sigma(lmN), 
  #                                 edf = df.residual(lmN))) %>% 
  #   as_tibble()
  
  NPexcr_c <- NPexcr %>%  
    filter(Lake == '239', Year == 2015) 
  transc <- NPexcr_c$massnorm.N.excr
  NPexcr_t <- NPexcr %>%  
    filter(Lake == '222', Year == 2015) 
  transt <- NPexcr_t$massnorm.N.excr
  
  # attempt 1 ----
  # Function to calculate MDD for ANOVA
  calculate_MDD_ANOVA <- function(anova_result, alpha = 0.05) {
    # Extract degrees of freedom for error and total
    df_error <- anova_result$"Df"[2]  # Degrees of freedom for error
    df_total <- sum(anova_result$"Df")  # Total degrees of freedom
    
    # Calculate the critical F-value
    critical_f_value <- qf(1 - alpha, df_error, df_total - df_error)
    
    # Calculate the MDD effect size
    num_groups <- 2 * 3  # Assuming 2 Lake groups and 3 Year groups
    mdd_effect_size <- sqrt(critical_f_value / num_groups)
    
    return(mdd_effect_size)
  }
  
  # Calculate MDD for the provided ANOVA result
  log_mdd_effect_size <- calculate_MDD_ANOVA(anova_N)
  # Back-transform MDD effect size
  mdd_effect_size <- 10^log_mdd_effect_size
  
  # Calculate the upper CI for ANOVA
  upper_CI_ANOVA <- mean(transc) - mean(transt) +
    qf(0.95, df_error, df_total - df_error) *
    sqrt(var(c(transc - mean(transc), transt - mean(transt)))) *
    sqrt(2 / N)
  
  
  # attempt 2 ----
  # Define parameters
  alpha <- 0.05  # Significance level
  power <- 0.80  # Desired power
  effect_size <- 2  # Desired effect size (Cohen's d)
  
  # Extract MSE
  mse <- tail(anova_N$"Mean Sq", 1)
  
  # Degrees of freedom
  # Find the row corresponding to the interaction term(s)
  interaction_row <- grep("\\*", rownames(anova_N))
  # Degrees of freedom associated with the interaction term(s)
  df_interaction <- anova_N[interaction_row, "Df"] 
  df_residual <- lmN$df.residual
  
  # Calculate critical F-value
  critical_f_value <- qf(1 - alpha, df_interaction, df_residual)
  
  # Sample sizes (if different by group)
  sample_sizes <- c(20, 15, 18, 7, 8, 4)  # Sample sizes per group or combination of levels
  
  # Calculate MDD
  # The MDD represents the smallest difference between group means that you can 
  # detect with the specified sample size, power, and significance level.
  mdd_log <- sqrt((2 * (1 / mean(sample_sizes)) * df_residual) / critical_f_value)
  mdd <- 10^mdd_log
  mdd

  
  pwcN <- emm(lmN)
  
  # P excretion
  lmP <- lmer(log10(massnorm.P.excr) ~ Treatment + (1 | Lake) + (1 | Lake:Year), 
            data = NPexcr)
  lmPalt <- lmer(log(massnorm.P.excr) ~ Period * SiteClass + (1 |Lake), data = NPexcr)
  lmP <- lm(log(massnorm.P.excr) ~ Lake*Year, 
             data = NPexcr)
  
  anova(lmP, ddf="Kenward-Roger")
  anova(lmPalt, ddf="Kenward-Roger")
  summary(lmPalt)
  pwcP <- emm(lmP)
  
  # N:P excretion
  lmNP <- lm(log(massnorm.NP.excr) ~ Lake*Year, 
            data = NPexcr %>% filter(Year != '2014'))
  nova(lmNP)
  summary(lmNP)
  pwcNP <- emm(lmNP)
  
  # make lm table based on all model ----
  lm_models <- list(
    "Mass" = lmMass,
    "Mass-normalized N excretion" = lmN,
    "Mass-normalized P excretion" = lmP,
    "Mass-normalized N:P excretion" = lmNP
  )
  combined_anova <- bind_rows(lapply(names(lm_models), function(model_name) {
    model <- lm_models[[model_name]]
    anova_result <- rownames_to_column(anova(model), 
                                       var = "Predictors") %>% as_tibble
    anova_result$groupname <- model_name
    return(anova_result)
  }))
  
  # make contrast table based on all models
  contrasts <- bind_rows(pwcMass[['contrasts']],
                         pwcN[['contrasts']], pwcP[['contrasts']], 
                         pwcNP[['contrasts']], .id = 'column_label')
  
  contrasts <- contrasts %>% 
    mutate(test = if_else(column_label == 1,'Mass', 
                          if_else(column_label == 2,  "Mass-normalized N excretion",
                                  if_else(column_label == 3,  "Mass-normalized P excretion", 
                                          "Mass-normalized N:P excretion")))) %>% 
    select(-column_label)
  
  # wilcox test ----
  excr.Tag <- NPexcr %>% 
    dplyr::filter(!(Year %in% c('2012', '2022'))) 
  # perform bootstrapping with 2000 replications
  set.seed(1)
  reps <- boot(data = excr.Tag, statistic = coef_function, 
               R = 2000, formula = massnorm.Tag.excr ~ Year)
  plot(reps)
  
  # Ag
  gghistogram(excr.Tag, x = 'massnorm.Tag.excr', rug = T, 
              fill = 'Year', add = 'mean')
  ggqqplot(excr.Tag, x = "massnorm.Tag.excr", facet.by = "Year")
  excr.Tag %>%  group_by(Year) %>% shapiro_test(massnorm.Tag.excr)
  var.test(massnorm.Tag.excr ~ Year,  excr.Tag, alternative = 'two.sided')
  w.Tag <- wilcox.test(massnorm.Tag.excr ~ Year, excr.Tag, var.equal = F)
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
    iter = rep(1:2500, each = 132)
  ) 
  
  # ..Schiettekatte FishStoich model ----
  # Yellow perch
  total.length <- exp((log(NPexcr$Wet.mass[!is.na(NPexcr$Wet.mass)]/
                             parammod$lwa_m)/parammod$lwb_m))
  wet.mass <- NPexcr$Wet.mass[!is.na(NPexcr$Wet.mass)]
  # ..Run model 
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
  
  # ..get dataset ready for excretion simulation models ----
  iter_yp <- iter_yp %>% 
    mutate(
      tl = rep(total.length, each = 2500),
      wet.mass = rep(NPexcr$Wet.mass[!is.na(NPexcr$Wet.mass)], each = 2500),
      dry.mass = rep(NPexcr$Mass[!is.na(NPexcr$Mass)], each = 2500),
      N.excretion.mod = Fn * 10 ^ 6 / 24,
      P.excretion.mod = Fp * 10 ^ 6 / 24,
      C.ingestion.mod = Ic * 10 ^ 6 / 24
    )
  
  # ..predict P excretion increase with diet increase ----
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
  
  
  # ..NPC limitation ----
  limitation <- limitation(FStoichm, plot = F)
  a <- parammod$lwa_m
  b <- parammod$lwb_m
  limitation <- limitation %>% mutate(mass = a*tl^b,
                                      dry.mass = mass*0.25)
  
  
  