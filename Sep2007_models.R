# fit model with data from Sepulveda et al 2007

install.packages("lme4")
library(lme4)
library(here)
  library(tidyverse)

  Sep2007 <- readxl::read_xlsx("C:/Users/Claudia/Desktop/Doctorado/MetabolicRate_WhiteShark/Sepulveda_data_2007.xlsx")
  
  # clean up col names and log transform MR, mass, temp
  Sep2007 <- Sep2007 %>%
    rename(
      rel_MR = MR_mgO2_kg_h,
      U = U_TLs,
    ) %>%
    mutate(
      mass_g = Mass_kg * 1000,
      logMassKG = log10(Mass_kg),
      logMassG = log10(mass_g),
      MR_units = "mgO2_kg_h",
      U_units = "TL_s",
      abs_MR = rel_MR * Mass_kg,
      log_rel_MR = log10(rel_MR),
      log_abs_MR = log10(abs_MR)
    )
  
  # convert temp to arrhenius
  
  #Boltzmann constant
  Bol <- 8.617343e-05

  #convert MR temp into kelvin
  Sep2007 <- Sep2007 %>% 
    mutate(Temp_K = Temp_C + 273.15,
           inv_temp = (1/(Bol*Temp_K)))
  
  #standardize temperature associated with MR
  Sep2007$inv_temp_s <- as.vector(scale(Sep2007$inv_temp))
  

  # plot
  ggplot(Sep2007) +
    geom_point(aes(x = U, y = log_rel_MR, color = Temp_C))
  
  # plot
  ggplot(Sep2007) +
    geom_point(aes(x = logMassKG, y = log_abs_MR, color = Temp_C))
  
   ggplot(Sep2007) +
    geom_point(aes(x = Mass_kg, y = abs_MR, color = Temp_C))
 
  # model 1: MR ~ U
  mod1 <- lme4::lmer(log_abs_MR ~ U + (1|Individual), data = Sep2007)
  summary(mod1)  
  coef(mod1)
  AIC(mod1)
  
  mod_no_re <- lm(log_abs_MR ~ U, data = Sep2007)
  summary(mod_no_re)  
  AIC(mod_no_re)
  
  # model 2: MR ~ U + Mass
  mod2 <- lme4::lmer(log_abs_MR ~ U + logMassKG + (1|Individual), data = Sep2007)
  summary(mod2)
  AIC(mod2)
  
  # model 2 with no random effect
  mod2.1 <- lm(log_abs_MR ~ U + logMassKG, data = Sep2007)
  summary(mod2.1)
  AIC(mod2.1)
  
  # model 3: MR ~ U + Mass + Temp
  mod3 <- lme4::lmer(log_abs_MR ~ U + logMassKG + inv_temp_s + 
                       (1|Individual), data = Sep2007)
  summary(mod3)
  AIC(mod3)
  
  # model 3 with no random effect
  mod3.1 <- lm(log_abs_MR ~ U + logMassKG + inv_temp_s, data = Sep2007)
  summary(mod3.1)
  AIC(mod3.1)
  
  # model with interaction 
  mod4 <- lm(log_abs_MR ~ U + logMassKG + U:logMassKG, data = Sep2007)
  summary(mod4)
  AIC(mod4)
  
  # model with only mass 
  mod5 <- lm(log_abs_MR ~ logMassKG, data = Sep2007)
  summary(mod5)
  AIC(mod5)
  
  # model with only mass 
  mod6 <- lm(log_abs_MR ~ logMassKG + inv_temp_s, data = Sep2007)
  summary(mod6)
  AIC(mod6)
 
  
  
  # lowest AIC for most parsimonius model is model 2.1 - MR ~ speed + mass (with no random effect of individual)
  
  save(mod2.1, file = here("./output/refit_Sep2007_mod.rds"))
  
  ## predict with Claudia's data
  #dat_speed <- readRDS(here("data", "datostransformados.rds")) #%>%
  # # dplyr::select(
  # #   -Date, -Time, -Depth, 
  # #   -DT, -Sex, -Moon, -Tide, 
  # #   -ShoreDistance, -Thermocline, 
  # #   -ShipsNumber, -Hora, -Bathy,
  # #   -contains("NES")
  # # )
  #
  #full_dat <- data.table::fread(file = here("./data/full_dat.csv")) %>%
  # dplyr::select(
  #    -Date, -Time, -Depth, 
  #    -DT, -Sex, -Moon, -Tide, 
  #    -ShoreDistance, -Thermocline, 
  #    -ShipsNumber, -Hora, -Bathy,
  #    -contains("NES")
  #  )
 #
  #full_dat_sum <- full_dat %>%
  #  distinct_at(vars(Length, Weight)) 
  #
  #dat <- left_join(dat_speed, full_dat_sum, by = "Length")
  #
  ## predict MR from mod2 - NEED TO WORK ON THIS
  #new_dat <- data.frame(
  #  U = seq(from = min(dat$speed), to = max(dat$speed), length.out = 
  #),
  #fits <- predict(mod2, newdata = new_dat, re.form = NA))
  #
  ####
  #
  
 # t <- Sep2007 %>% map( ~ tidy(lm(logMR ~ U + logMass)))
  
  
  
  
  
  
  # 1. Obtener predicciones de ambos modelos
  Sep2007$pred_no_mass <- predict(mod_no_re)       # modelo sin masa
  Sep2007$pred_with_mass <- predict(mod2.1)        # modelo con masa
  
  # 2. Calcular sobreestimación en escala log10
  Sep2007$overestimation_log <- Sep2007$pred_no_mass - Sep2007$pred_with_mass
  
  # 3. Convertir a escala original (mgO2/h)
  Sep2007$overestimation_pct <- (10^Sep2007$pred_no_mass - 10^Sep2007$pred_with_mass) / 
    (10^Sep2007$pred_with_mass) * 100
  
  # 4. Calcular promedio de sobreestimación
  mean_overestimation <- mean(Sep2007$overestimation_pct, na.rm = TRUE)
  print(paste("Sobreestimación promedio:", round(mean_overestimation, 2), "%"))
  
  