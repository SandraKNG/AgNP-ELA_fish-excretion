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
  lmN <- lm(log(massnorm.N.excr) ~ Lake * Year, data = NPexcr %>% filter(Year != '2014'))
  anova_N <- anova(lmN)
  summary(lmN)
  
  
  # lm models
  
  lmN <- lm(log(massnorm.N.excr) ~ Lake*Year, 
            data = NPexcr %>% filter(Year != '2014'))
  anova(lmN)
  summary(lmN)
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
  anova(lmP)
  summary(lmP)
  pwcP <- emm(lmP)
  
  # N:P excretion
  lmNP <- lm(log(massnorm.NP.excr) ~ Lake*Year, 
             data = NPexcr %>% filter(Year != '2014'))
  anova(lmNP)
  P <- anova(lmNP)
  summary(lmNP)
  pwcNP <- emm(lmNP)
  
  # MDD Function ----
  MDD <- function(N1, N2, variance1, variance2 = NULL, alpha = 0.05, 
                  two.sided = TRUE, var.equal = TRUE) {
    m <- list()
    m$n1 <- N1
    m$n2 <- N2
    m$variance1 <- variance1
    m$variance2 <- ifelse(is.null(variance2), variance1, variance2)
    m$alpha <- alpha
    m$two.sided <- two.sided
    
    if (var.equal) {
      df <- N1 + N2 - 2
      nu <- sqrt(variance1) * sqrt(1 / N1 + 1 / N2)
    } else {
      df <- (variance1 / N1 + m$variance2 / N2)^2 / 
        ((variance1 / N1)^2 / (N1 - 1) + (m$variance2 / N2)^2 / (N2 - 1))
      nu <- sqrt(variance1 / N1 + m$variance2 / N2)
    }
    
    m$df <- df
    m$nu <- nu
    qt_val <- qt(alpha / ifelse(two.sided, 2, 1), df = df, lower.tail = FALSE)
    m$mdd <- qt_val * nu
    
    class(m) <- "MDD"
    return(m)
  }
  
  # Initialize Output DataFrame ----
  columns_to_analyze <- c("massnorm.N.excr", "massnorm.P.excr", "massnorm.NP.excr")
  output_size <- length(unique(NPexcr$Year)) * length(columns_to_analyze)
  out <- data.frame(
    Year = rep(NA, output_size),
    Variable = rep(NA, output_size),
    Nc = rep(NA, output_size),
    Nt = rep(NA, output_size),
    control_mean = rep(NA, output_size),
    treatment_mean = rep(NA, output_size),
    p = rep(NA, output_size), 
    pMDD = rep(NA, output_size), 
    pMDE = rep(NA, output_size), 
    pCI = rep(NA, output_size),
    mdd_exc = rep(NA, output_size),
    mde_exc = rep(NA, output_size),
    upperCI_exc = rep(NA, output_size)
  )
  
  # Data Processing and Analysis ----
  i <- 1
  for (y in unique(NPexcr$Year)) {
    for (col in columns_to_analyze) {

      
      # Filter data for control and treatment groups
      control <- NPexcr %>% 
        filter(Year == y, Lake == 239, !is.na(.data[[col]])) %>% 
        pull(.data[[col]])
      
      treatment <- NPexcr %>% 
        filter(Year == y, Lake == 222, !is.na(.data[[col]])) %>% 
        pull(.data[[col]])
      
      # cat("Control: ", y, col, paste(unlist(control), collapse = ", "), "\n")
      # cat("Treatment ", y, col, paste(unlist(treatment), collapse = ", "), "\n")
      
      # Skip iteration if data is insufficient
      if (length(control) < 2 | length(treatment) < 2) next
      
      Nc <- length(unique(control))
      Nt <- length(unique(treatment))
      N_av <- (Nc + Nt) / 2
      
      transc <- log(control + 1)
      transt <- log(treatment + 1)
      
      # cat("Control: ", y, col, paste(unlist(transc), collapse = ", "), "\n")
      # cat("Treatment ", y, col, paste(unlist(transc), collapse = ", "), "\n")
      
      # Combined Variance Calculation
      combined_variance <- var(c(transc - mean(transc), transt - mean(transt)))
      
      # MDD Calculation
      if (col == 'massnorm.P.excr') {
        mdd_ln <- MDD(N1 = Nc, N2 = Nt, variance1 = (var(transc) + var(transt)) / 2, 
                    alpha = 0.05, two.sided = T, var.equal = F)
      } else {
        mdd_ln <- MDD(N1 = Nc, N2 = Nt, variance1 = (var(transc) + var(transt)) / 2, 
                      alpha = 0.05, two.sided = T, var.equal = T)
      }
      
      
      # Post Hoc Power Calculation
      postpower <- power.t.test(n = N_av, delta = NULL, 
                                sd = sqrt(combined_variance), 
                                sig.level = 0.05, alternative = "two.sided", 
                                type = "two.sample", power = 0.8)
      
      # Upper Confidence Interval Calculation
      if (mean(transc) > mean(transt)) {
        upperCI <- mean(transc) - mean(transt) +
          qt(0.05, df = Nc + Nt - 2, lower.tail = FALSE) *
          sqrt(combined_variance) * sqrt(2 / (Nc + Nt))
      } else {
        upperCI <- mean(transt) - mean(transc) +
          qt(0.05, df = Nc + Nt - 2, lower.tail = FALSE) *
          sqrt(combined_variance) * sqrt(2 / (Nc + Nt))
      }
      
      # Backtransformation MDD, MDE and upper bound of the CI:
      mdd_exc <- (exp(mean(transc)) - exp(mean(transc) - mdd_ln$mdd))
      mde_exc <- (exp(mean(transc)) - exp(mean(transc) - postpower$delta))
      upperCI_exc <- (exp(mean(transc)) - exp(mean(transc) - upperCI))
      
      # Relation to Control Mean
      control_mean <- round(exp(mean(transc) - 1), 2)
      treatment_mean <- round(exp(mean(transt) - 1), 2)
      pMDD <- round(100 * mdd_exc / control_mean - 1, 2)
      pMDE <- round(100 * mde_exc / control_mean - 1, 2)
      pCI <- round(100 * upperCI_exc / control_mean - 1, 2)
      
      # Store Results in Output DataFrame
      if (col == 'massnorm.P.excr') {
      test <- t.test(transc, transt, alternative = "two.sided", var.equal = F)
      } else {
        test <- t.test(transc, transt, alternative = "two.sided", var.equal = T)
      }
      out[i,] <- c(y, col, Nc, Nt, control_mean, treatment_mean, round(test$p.value, 2), pMDD, pMDE, pCI, 
                   round(mdd_exc, 2), round(mde_exc, 2), round(upperCI_exc, 2))
      
      i <- i + 1
    }
  }

  MDD_results <- out %>% 
    filter(!is.na(Year)) %>% 
    select(Year, Variable, Nc, Nt, control_mean, treatment_mean, p, pMDD, mdd_exc, pCI, upperCI_exc) %>% 
    rename(N_reference = Nc,
           N_experimental = Nt,
           reference_mean = control_mean,
           experimental_mean = treatment_mean,
           MDD = mdd_exc,
           upperCI = upperCI_exc) %>% 
    mutate(Year = if_else(Year == 2012, "Pre-addition",
                          if_else(Year == 2014, "Year 1",
                                  if_else(Year == 2015, "Year 2", "Post-addition"))),
           Variable = case_when(
             Variable == "massnorm.N.excr" ~ "mass-specific N excretion",
             Variable == "massnorm.P.excr" ~ "mass-specific P excretion",
             Variable == "massnorm.NP.excr" ~ "mass-specific N:P excretion",
             TRUE ~ Variable  # Keep any other variables unchanged
             )) %>% 
    mutate(Year = fct_relevel(Year, "Pre-addition", "Year 1", "Year 2", "Post-addition")) %>%
    arrange(Variable, Year)
  

  lm_models <- list(
    "Mass" = lmMass,
    "Mass-specific N excretion" = lmN,
    "Mass-specific P excretion" = lmP,
    "Mass-specific N:P excretion" = lmNP
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
  
  
  