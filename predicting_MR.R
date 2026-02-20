# calculating body speed, converting length to mass, and calculating metabolic rate

# load necessary packages



### IDEAS

# figure 1 - fit lines for all equations (Ezcurra, Sep, Graham, refit Sep)
# figure 2 - predicting MR with these different equations

  library(tidyverse)
  library(lubridate)
  library(fuzzyjoin)
  library(here)
  library(nlme)
  library(broom)
  library(purrr)
  library(tinter)

  # read in data
  datos <- readRDS(file = here("./data/datostransformados.rds"))

  # convert time 
  datos <- datos %>%
    mutate(Date = force_tz(with_tz(Date, tzone = "America/Los_Angeles"), tzone ="America/Los_Angeles"),
           Time = force_tz(with_tz(Time, tzone = "America/Los_Angeles"), tzone = "America/Los_Angeles"),
           DT = force_tz(with_tz(DT, tzone = "America/Los_Angeles"), tzone = "America/Los_Angeles"))

  # calculate speed 
  datos.filtro <- datos %>% drop_na(speed) %>% filter(!is.infinite(speed))

  df <- datos.filtro %>% group_by(SharkID) %>%
    summarise(diffDT.min = quantile(diffDT, probs = 0.025),
              diffDT.max = quantile(diffDT, probs = 0.975),
              Distance.min = quantile(Distance, probs = 0.025),
              Distance.max = quantile(Distance, probs = 0.975),
              Speed.min = quantile(speed, probs = 0.025),
              Speed.max = quantile(speed, probs = 0.975))
  
  # calculate corrected speed
  datos.filtro <- datos.filtro %>% 
    filter(between(speed, min(df$Speed.min), max(df$Speed.max))) %>%
    filter(between(diffDT, min(df$diffDT.min), max(df$diffDT.max))) %>%
    filter(between(Distance, min(df$Distance.min), max(df$Distance.max)))
  
  # convert body length to mass using length-weight regressions
  # regression equation for Mass to TL in Mollet & Cailliet 2006 is: M = 7.914 * TL^3.0958, where mass is in kg and length in m
  # in log space, this equation is: log M = log(7.914) + 3.0958 * log(TL)
  
  lwr_func <- function(Length, a = 7.914, b = 3.0958){
    
    M = a * (Length ^ b)
  }
  
  # way 1
  all_lengths <- datos.filtro$Length
  all_weights <- sapply(all_lengths, lwr_func)
  
  datos.filtro$bodyMass_kg <- all_weights
  
  # function to predict the mean speed by hour for each individual based on the loess model
  mean_speed_fun <- function(SharkID){

    df <- datos.filtro[datos.filtro$SharkID == SharkID ,]

    df <- df %>% mutate(horario = round(((second(Time) + (minute(Time)*60) + (hour(Time)*3600))/3600), digits = 2))

    df <- df %>% filter(speed > 0)

    p1 <- ggplot(data = df, aes(x = horario, y = speed)) + 
      geom_smooth(method = "loess", span = 0.1)

    mean_speed <- ggplot_build(p1)$data[[1]]

    mean_speed <- mean_speed %>% dplyr::select(x, y, se)
    colnames(mean_speed) <- c("horario", "mean_speed", "se_speed")

    df <- df %>% difference_inner_join(mean_speed, by = "horario", max_dist = 0.1)
  
    
    }

    IDs <- unique(datos.filtro$SharkID)

    datos_speed <- map(IDs, mean_speed_fun) %>% bind_rows()
    
    datos_speed$horario <- datos_speed$horario.y


  ############################################################################
  # Calculating metabolic rate from published equations in the literature ####
  ############################################################################

  # We will use the following equations to calculate metabolic rate from speed and mass


  # Sepulveda et al. 2007: LogMO 2= 2.0937(+-0.06) + 0.97(+-0.13)U
  # 2nd equation in Sepulveda: 2.3716(+-0.03) + 0.58(+-0.12)U
  # Graham et al. 1990 (and Semmens et al. 2013 and Anderson et al. 2022): log * MO 2 = (0.58*U) + (log*246)
  # "New" equation from Sepulveda et al. 2007 that was refit with a random effect and included body mass

  # mean_speed is in m/s and we need to convert to body length/s to use speed to predict MR
    
  datos_speed <- datos_speed %>%
    mutate(mean_speed_cm_s = 100 * mean_speed,
           length_cm = Length * 100,
           bl_s = mean_speed_cm_s/length_cm,
           log_mass_kg = log10(bodyMass_kg))


  # Ezcurra et al. 2012 - absolute metabolic rate is response variable
  datos_speed$AMR_E_abs <- 458.5*(datos_speed$bodyMass_kg^0.79)
  datos_speed$AMR_E_abs_log <- log10(datos_speed$AMR_E_abs) # log this to be comparable with other estimates of MR below
  datos_speed$AMR_E_rel <- datos_speed$AMR_E_abs/datos_speed$bodyMass_kg
  datos_speed$AMR_E_rel_log <- log10(datos_speed$AMR_E_rel)


  ## Ezcurra et al. 2012 - absolute metabolic rate is response variable
  #datos_speed <- data.frame(datos_speed, AMR_E_raw =0)
  #datos_speed$AMR_E_raw <- 458.5*(datos_speed$bodyMass_kg^0.79)
  #datos_speed$AMR_E <- log10(datos_speed$AMR_E_raw) # log this to be comparable with other estimates of MR below
  #datos_speed$AMR_E_rel <- datos_speed$AMR_E/datos_speed$log_mass_kg

  # no speed in equation
  
  # Sepulveda et al. 2007 - equation Sepulveda fit for paper, log relative MR
  datos_speed$AMR_S1_rel_log <- 2.0937 + (0.97 * datos_speed$bl_s)

  # first, backtransform from log space and then calculate absolute MR
  datos_speed$AMR_S1_abs <- (10^(datos_speed$AMR_S1_rel_log)) * datos_speed$bodyMass_kg
  
  # log absolute MR
   datos_speed$AMR_S1_abs_log <- log10(datos_speed$AMR_S1_abs)

  
  # Sepulveda et al. 2007 - equation from Graham (1990), relative MR
  datos_speed$AMR_S2_rel_log <- 2.3716 + (0.58 * datos_speed$bl_s)
  # speed is in TL/s where TL is in cm
  
  # backtransform and convert to absolute MR
  datos_speed$AMR_S2_abs <- (10^(datos_speed$AMR_S2_rel_log)) * datos_speed$bodyMass_kg
  
  # log absolute MR
  datos_speed$AMR_S2_abs_log <- log10(datos_speed$AMR_S2_abs)

  
  # Graham et al. 1990 directly - very similar to above, relative MR
  datos_speed$AMR_G_rel_log <- (0.58 * datos_speed$bl_s) + log10(246)
  # speed is in TL/s where TL is in cm
  
  # backtransform and convert to absolute MR
  datos_speed$AMR_G_abs <- (10^(datos_speed$AMR_G_rel_log)) * datos_speed$bodyMass_kg
 
 # log absolute MR
  datos_speed$AMR_G_abs_log <- log10(datos_speed$AMR_G_abs)

  # refit data from Sepulveda paper 
  
  #here, we fit a model with body mass as a covariate, instead of using mass-specific or relative
  #MR as in the Sepulveda equations
  # output is absolute MR in log10 space
  load(file = here("./output/refit_Sep2007_mod.rds"))

  # mod2 is lowest model 
    #mod2.1 <- lm(logMR ~ U + logMass, data = Sep2007)
    #summary(mod2.1)
    #AIC(mod2.1)

# predict absolute MR for individuals in Claudia's df
  
  nd <- tibble(
    U = datos_speed$bl_s,
    logMassKG = datos_speed$log_mass_kg
  
  )
  
  # use predict() function in R
  p <- predict(mod2.1, newdata = nd)
  
  datos_speed$p <- p
  # this is absolute metabolic rate in log10 space!
  
  # back-transform 
  datos_speed$raw_p <- 10^(datos_speed$p)
  # this is absolute MR 
  
  datos_speed <- datos_speed |>
    mutate(rel_mr = raw_p/bodyMass_kg,
           log_rel_mr = log10(rel_mr),
           log_abs_mr = p)
  
  
  # plots #####


  # select only columns needed at the moment and remove raw predictions from Ezcurra et al. 2001
  datos_speed_trim <- datos_speed %>%
    select(SharkID, bodyMass_kg, bl_s, contains("log")) 
 

  # convert df into long format for plotting
  long_df <- datos_speed_trim %>%
    pivot_longer(
      cols = contains("log") & !contains("mass"),
      names_to = "metabolic_types",
      values_to = "metabolic_rate"
    )
  
  # plot
  long_df |>
  filter(str_detect(metabolic_types, "abs")) |>
  ggplot() +
    geom_line(aes(x = bl_s, y = metabolic_rate, group = metabolic_types, color = metabolic_types)) +
    facet_wrap(~ SharkID, scales = "free") 
  
  long_df |>
  filter(str_detect(metabolic_types, "rel")) |>
  ggplot() +
    geom_line(aes(x = bl_s, y = metabolic_rate, group = metabolic_types, color = metabolic_types)) +
    facet_wrap(~ SharkID, scales = "free") 
  
  long_df |>
  filter(str_detect(metabolic_types, "abs")) |>
   ggplot() +
    geom_line(aes(x = bl_s, y = metabolic_rate, color = SharkID)) +
    facet_wrap(~ metabolic_types, scales = "free") 

 
# add sepulveda data points on this figure


  Sep2007 <- readxl::read_xlsx(here("./data/Sepulveda_data_2007.xlsx"))

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
  
 
  long_df |>
  filter(str_detect(metabolic_types, "rel")) |>
  ggplot() +
  geom_line(aes(x = bl_s, y = metabolic_rate, group = metabolic_types, color = metabolic_types)) +
  facet_wrap(~ SharkID, scales = "free") +
  geom_point(data = Sep2007, aes(x = U, y = log_rel_MR), color = "black")
  
  # does geom line make sense - not prediction line/model, individual fits
 
    
  long_df <- long_df |>
    mutate(metabolic_type_name = case_when(
      metabolic_types == "AMR_E_abs_log" ~ "Ezcurra et al. 2012,\nlog absolute MR",
      metabolic_types == "AMR_E_rel_log" ~ "Ezcurra et al. 2012,\nlog relative MR",
      metabolic_types == "AMR_G_rel_log" ~ "Graham 2001,\nlog relative MR",
      metabolic_types == "AMR_G_abs_log" ~ "Graham 2001,\nlog absolute MR",
      metabolic_types == "AMR_S1_abs_log" ~ "Sepulveda et al. 2007 Graham equation,\nmlog absolute MR",
      metabolic_types == "AMR_S1_rel_log" ~ "Sepulveda et al. 2007 Graham equation,\nmlog relative MR",
      metabolic_types == "AMR_S2_rel_log" ~ "Sepulveda et al. 2007,\nmlog relative MR",
      metabolic_types == "AMR_S2_abs_log" ~ "Sepulveda et al. 2007,\nmlog absolute MR",
      metabolic_types == "log_abs_mr" ~ "refit model from Sepulveda et al. 2007,\nmlog absolute MR",
      metabolic_types == "log_rel_mr" ~ "refit model from Sepulveda et al. 2007,\nmlog relative MR"
    ))
  
  pred_MR_ind <- 
    long_df |>
    filter(str_detect(metabolic_types, "abs")) |>
    ggplot() +
    geom_line(aes(x = bl_s, y = metabolic_rate, group = metabolic_type_name, color = metabolic_type_name)) +
    facet_wrap(~ SharkID) +
      ylab("Predicted (log10) metabolic rate") +
      xlab("Swimming speed, body lengths/s") +
    ggsidekick::theme_sleek()
  
  ggsave(pred_MR_ind, file = here("./output/MR_prediction_figure_ind.png"),
         height = 7, width = 10, units = "in")

  
 
 pred_MR_mrt <- 
    long_df |>
    filter(str_detect(metabolic_types, "abs")) |>
    ggplot() +
    geom_line(aes(x = bl_s, y = metabolic_rate, group = SharkID, color = SharkID)) +
    facet_wrap(~ metabolic_type_name) +
      ylab("Predicted (log10) metabolic rate") +
      xlab("Swimming speed, body lengths/s") +
    ggsidekick::theme_sleek()
  
  ggsave(pred_MR_mrt, file = here("./output/MR_prediction_figure_mr_types.png"),
         height = 7, width = 10, units = "in")


 
 
    # assign colors to each species
    
    MR_type <- unique(long_df$metabolic_types)
    IDs <- unique(long_df$SharkID)

    plot_cols <- crossing(
      SharkID = IDs,
      metabolic_types = MR_type
    )    

    meta_d <- long_df |>
      select(metabolic_types, metabolic_type_name) |>
      distinct_all()

    plot_cols <- left_join(plot_cols, meta_d)    

    plot_cols_abs <- plot_cols |>
      filter(str_detect(metabolic_types, "abs")) 
    
  # colors
  colfunc <- colorRampPalette(c("#B3E5FC", "#01579B"))
  colfunc(3)

  # pick colors
  color_palette <- colorRampPalette(c("red", "yellow", "blue", "green", "purple"))(11)


  # Create a function to generate palettes for each group
  generate_palette <- function(color_palette) {
    # Create an empty list to store the color palettes
    color_palettes <- list()
  
   # Loop over each base color and generate a palette of 3 colors
   for (i in 1:length(color_palette)) {

      color_palettes[[i]] <- tinter(color_palette[i], steps = 3)   
     
       }
   
   return(color_palettes)
  }

# Generate the palettes based on the base colors
palettes <- generate_palette(color_palette)

#palettes <-  lapply(palettes, function(x) head(x, -1))

# Print the generated color palettes
palettes <- unlist(palettes, recursive = TRUE)


plot_cols_abs$color <- palettes

plot_cols_abs <- plot_cols_abs |>
  mutate(color_set_up = paste0(SharkID, "_", metabolic_types))

t <- long_df |>
  filter(str_detect(metabolic_types, "abs")) |>
  mutate(color_set_up = paste0(SharkID, "_", metabolic_types))

t <-  left_join(t, plot_cols_abs) |>
rename(hex_code = color)


   pred_MR3 <-  
    t |>
    ggplot() +
    geom_point(aes(x = bl_s, y = metabolic_rate, color = color_set_up), alpha = 0.7) +
   # facet_wrap(~ SharkID, scales = "free") +
      ylab("Predicted (log10) metabolic rate") +
      xlab("Swimming speed, body lengths/s") +
      scale_color_manual(values = palettes) +
    ggsidekick::theme_sleek() +
     theme(legend.title = element_blank())
   
    ggsave(pred_MR3, file = here("./output/MR_prediction_figure_all_together.png"),
         height = 7, width = 10, units = "in")


  #### plot separate panels for each MR prediction ####
  
        
    
  long_df |>
        filter(str_detect(metabolic_types, "abs")) |>
    ggplot() +
    geom_point(aes(x = bl_s, y = metabolic_rate, color = SharkID)) +
    facet_wrap(~ metabolic_type_name, scales = "free") +
      ylab("Predicted (log10) metabolic rate") +
      xlab("Swimming speed, body lengths/s") +
    ggsidekick::theme_sleek()
  
  
    long_df |>
  filter(str_detect(metabolic_types, "abs")) |>
      ggplot() +
      geom_line(aes(x = bl_s, y = metabolic_rate, color = SharkID), alpha = 0.5) +
      facet_wrap(~ metabolic_type_name, scales = "free") +
        ylab("Predicted (log10) metabolic rate") +
        xlab("Swimming speed, body lengths/s") +
      ggsidekick::theme_sleek()
  
    
    