rm(list=ls(all=TRUE))

setwd("C:/Users/bened/GIT/Foodweb_Ecology_Mesograzer_Herbivory")
library(tibble)
library(readxl)
library(dplyr)
library(lme4)
library(performance)
library(car)
library(DHARMa)
library(report)
library(readxl)
library(ggplot2)
library(tidyverse)
library(lubridate)
library(dplyr)
library(tidyr)
library(ggbeeswarm)
library(car)
library(performance)
library(multcomp)
library(multcompView)
library(broom)
library(rstatix)
library(lmerTest)
library(here)
library(lme4)
library(broom.mixed)
library(emmeans)
library(lattice)
library(nlme)

T_A1 <- data.frame(read_xlsx("Matrices_HW.xlsx", sheet = "A1"))
T_A2 <- data.frame(read_xlsx("Matrices_HW.xlsx", sheet = "A2"))
T_B2 <- data.frame(read_xlsx("Matrices_HW.xlsx", sheet = "B2"))
T_C1 <- data.frame(read_xlsx("Matrices_HW.xlsx", sheet = "C1"))
T_C2 <- data.frame(read_xlsx("Matrices_HW.xlsx", sheet = "C2"))
T_D1 <- data.frame(read_xlsx("Matrices_HW.xlsx", sheet = "D1"))
T_D2 <- data.frame(read_xlsx("Matrices_HW.xlsx", sheet = "D2"))
T_E1 <- data.frame(read_xlsx("Matrices_HW.xlsx", sheet = "E1"))
T_E2 <- data.frame(read_xlsx("Matrices_HW.xlsx", sheet = "E2"))
T_F1 <- data.frame(read_xlsx("Matrices_HW.xlsx", sheet = "F1"))
T_F2 <- data.frame(read_xlsx("Matrices_HW.xlsx", sheet = "F2"))

biomass_phy <- read_xlsx("macrophytes-biomass.xlsx")
biomass_phy <- biomass_phy[-3, ]
biomass_phy$Fv_WW_g <- biomass_phy$Fv_WW_g*1000
biomass_phy$Fa_DW_g <- biomass_phy$Fa_DW_g*1000
biomass_phy$tot <- rowSums(biomass_phy[, c("Zm_C_mg", "Fv_WW_g", "Fa_DW_g")], na.rm = TRUE)
biomass_det <- read_csv("POC_mean_Aug_Sept_2015.csv")
biomass_det <- biomass_det[-3, ]
biomass_zoo <- read_xlsx("Biomass_HW.xlsx")

T_A1 <- as.data.frame(lapply(T_A1, as.numeric))
T_A2 <- as.data.frame(lapply(T_A2, as.numeric))
T_B2 <- as.data.frame(lapply(T_B2, as.numeric))
T_C1 <- as.data.frame(lapply(T_C1, as.numeric))
T_C2 <- as.data.frame(lapply(T_C2, as.numeric))
T_D1 <- as.data.frame(lapply(T_D1, as.numeric))
T_D2 <- as.data.frame(lapply(T_D2, as.numeric))
T_E1 <- as.data.frame(lapply(T_E1, as.numeric))
T_E2 <- as.data.frame(lapply(T_E2, as.numeric))
T_F1 <- as.data.frame(lapply(T_F1, as.numeric))
T_F2 <- as.data.frame(lapply(T_F2, as.numeric))

summarize_tank <- function(tank_matrix, tank_id, treatment, side) {
  Ga_herb <- sum(tank_matrix[1:3, 4])
  Ib_herb <- sum(tank_matrix[1:3, 5])
  Ll_herb <- sum(tank_matrix[1:3, 6])
  tot_herb <- Ga_herb + Ib_herb + Ll_herb
  
  Ga_det <- tank_matrix[7, 4]
  Ib_det <- tank_matrix[7, 5]
  Ll_det <- tank_matrix[7, 6]
  tot_det <- Ga_det + Ib_det + Ll_det
  
  tibble(
    ID = tank_id,
    treatment = treatment,
    side = side,
    Ga_herb, Ib_herb, Ll_herb, tot_herb,
    Ga_det, Ib_det, Ll_det, tot_det
  )
}

tank_list <- list(
  summarize_tank(T_A1, "A1", "0HW", 1),
  summarize_tank(T_A2, "A2", "0HW", 2),
  summarize_tank(T_B2, "B2", "3HW", 2),
  summarize_tank(T_C1, "C1", "1HW", 1),
  summarize_tank(T_C2, "C2", "1HW", 2),
  summarize_tank(T_D1, "D1", "0HW", 1),
  summarize_tank(T_D2, "D2", "0HW", 2),
  summarize_tank(T_E1, "E1", "3HW", 1),
  summarize_tank(T_E2, "E2", "3HW", 2),
  summarize_tank(T_F1, "F1", "1HW", 1),
  summarize_tank(T_F2, "F2", "1HW", 2)
)


tanks <- bind_rows(tank_list)

tanks

tanks$ratio_herb_det <- tanks$tot_herb/tanks$tot_det

tanks %>%
  group_by(treatment) %>%
  summarise(mean_herbivory = mean(tot_herb, na.rm = TRUE))

tanks %>%
  group_by(treatment) %>%
  summarise(mean_detrivory = mean(tot_det, na.rm = TRUE))

biomass_zoo$Biomass_mgC_m2 <- as.numeric(biomass_zoo$Biomass_mgC_m2)
tanks$Ga_zoo <- biomass_zoo$Biomass_mgC_m2[biomass_zoo$ID == "Ga"]
tanks$Ib_zoo <- biomass_zoo$Biomass_mgC_m2[biomass_zoo$ID == "Ib"]
tanks$Ll_zoo <- biomass_zoo$Biomass_mgC_m2[biomass_zoo$ID == "Ll"]
tanks$zoo_tot <- rowSums(tanks[, c("Ga_zoo", "Ib_zoo", "Ll_zoo")], na.rm = TRUE)

###############################################################################
###############################################################################
######## Anova für ratio phytoplakton #########################################
###############################################################################
###############################################################################

#  Herbivory 
tanks$herb <- tanks$tot_herb /  biomass_phy$tot 


tanks_clean <- na.omit(tanks[, c("herb", "treatment")])

shapiro_test(tanks$herb)

# ANOVA 
model_anova <- aov(herb ~ treatment, data = tanks_clean)
summary(model_anova)

# Tukey HSD-Posthoc-Test
tukey_result <- TukeyHSD(model_anova)
print(tukey_result)


tukey_letters <- multcompLetters4(model_anova, tukey_result)
letters_df <- as.data.frame.list(tukey_letters$treatment)
letters_df$treatment <- rownames(letters_df)


max_vals <- tanks_clean %>%
  group_by(treatment) %>%
  summarise(max_herb = max(herb, na.rm = TRUE))


max_vals$label <- letters_df[match(max_vals$treatment, letters_df$treatment), "Letters"]

treatment_colors <- c(
  "0HW" = "#add8e6", 
  "1HW" = "#f08080",  
  "3HW" = "#8b0000" 
)

# Plot
ggplot(tanks_clean, aes(x = treatment, y = herb, fill = treatment)) +
  geom_boxplot() +
  geom_text(data = max_vals, aes(x = treatment, y = max_herb + 0.00008, label = label), 
            color = "black", size = 5) +
  scale_fill_manual(values = treatment_colors) +
  labs(
    x = "Treatment",
    y = expression('Herbivory / PP biomass'~"["*day^{-1}*"]"),
    fill = "Treatment"
  ) +
  #coord_cartesian(ylim = c(0, 0.0005)) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.0001)) +
  theme_minimal()

###############################################################################
###############################################################################
##################### Anova herbivory #########################################
###############################################################################
###############################################################################

library(ggplot2)
library(multcompView)
library(dplyr)


tanks_clean <- na.omit(tanks[, c("tot_herb", "treatment")])

shapiro_test(tanks$tot_herb)

# ANOVA
model_anova <- aov(tot_herb ~ treatment, data = tanks_clean)
summary(model_anova)

# Tukey HSD-Posthoc-Test
tukey_result <- TukeyHSD(model_anova)
print(tukey_result)


tukey_letters <- multcompLetters4(model_anova, tukey_result)
letters_df <- as.data.frame.list(tukey_letters$treatment)
letters_df$treatment <- rownames(letters_df)


max_vals <- tanks_clean %>%
  group_by(treatment) %>%
  summarise(max_herb = max(tot_herb, na.rm = TRUE))


max_vals$label <- letters_df[match(max_vals$treatment, letters_df$treatment), "Letters"]

treatment_colors <- c(
  "0HW" = "#add8e6",  
  "1HW" = "#f08080",   
  "3HW" = "#8b0000"   
)

# Plot
ggplot(tanks_clean, aes(x = treatment, y = tot_herb, fill = treatment)) +
  geom_boxplot() +
  geom_text(data = max_vals, aes(x = treatment, y = max_herb + 0.0008, label = label), 
            color = "black", size = 5) +
  scale_fill_manual(values = treatment_colors) +
  labs(
    x = "Treatment",
    y = expression(Herbivory~"["*mg~C~m^{-2}~day^{-1}*"]"),
    fill = "Treatment"
  )+
  coord_cartesian(ylim = c(0, 200)) +
  theme_minimal()



################################################################################
######## ratio detritus ########################################################
################################################################################

tanks$ratio_tot_det <- tanks$tot_det / biomass_det$`mean(biomass_mgC_m2)`

shapiro_test(tanks$ratio_tot_det)

# ANOVA-Modell
anova_model_det <- aov(ratio_tot_det ~ treatment, data = tanks)

# Post-hoc Test (Tukey)
tukey_det <- TukeyHSD(anova_model_det)
print(tukey_det)

tukey_letters_det <- multcompLetters4(anova_model_det, tukey_det)
letters_df <- as.data.frame.list(tukey_letters_det$treatment)
letters_df$treatment <- rownames(letters_df)
colnames(letters_df)[1] <- "label"

mean_max_det <- tanks %>%
  group_by(treatment) %>%
  summarise(mean_ratio = mean(ratio_tot_det, na.rm = TRUE),
            max_ratio = max(ratio_tot_det, na.rm = TRUE)) %>%
  left_join(letters_df, by = "treatment")

# Plot 
ggplot(tanks, aes(x = treatment, y = ratio_tot_det, fill = treatment)) +
  geom_boxplot() +
  geom_text(data = mean_max_det, aes(x = treatment, y = max_ratio * 1.05, label = label),
            color = "black", size = 5) +
  scale_fill_manual(values = treatment_colors) +
  scale_y_continuous(limits = c(0, 0.3)) +
  labs(
    x = "Treatment",
    y = expression(Detritivory~"["*day^{-1}*"]"),
    fill = "Treatment"
  ) +
  theme_minimal()

################################################################################
################################################################################
############################ detritivory ANOVA #################################
################################################################################
################################################################################

shapiro_test(tanks$tot_det)

# ANOVA-Modell
anova_model_det <- aov(tot_det ~ treatment, data = tanks)

# Post-hoc Test (Tukey)
tukey_det <- TukeyHSD(anova_model_det)
print(tukey_det)

tukey_letters_det <- multcompLetters4(anova_model_det, tukey_det)
letters_df <- as.data.frame.list(tukey_letters_det$treatment)
letters_df$treatment <- rownames(letters_df)
colnames(letters_df)[1] <- "label"

mean_max_det <- tanks %>%
  group_by(treatment) %>%
  summarise(mean_ratio = mean(tot_det, na.rm = TRUE),
            max_ratio = max(tot_det, na.rm = TRUE)) %>%
  left_join(letters_df, by = "treatment")

# Plot
ggplot(tanks, aes(x = treatment, y = tot_det, fill = treatment)) +
  geom_boxplot() +
  geom_text(data = mean_max_det, aes(x = treatment, y = max_ratio * 1.05, label = label),
            color = "black", size = 5) +
  scale_fill_manual(values = treatment_colors) +
  scale_y_continuous(limits = c(0, 80)) +
  labs(
    x = "Treatment",
    y = expression(Detritivory~"["*mg~C~m^{-2}~day^{-1}*"]"),
    fill = "Treatment"
  ) +
  theme_minimal()

################################################################################
################################################################################
################### ratio Herbivory / Detritivory ANOVA ########################
################################################################################
################################################################################

shapiro.test(tanks$ratio_herb_det) # no normal distribution

tanks$treatment <- factor(tanks$treatment,
                          levels = c("0HW", "1HW", "3HW"))

# GLM 
glm_model <- glm(ratio_herb_det ~ treatment, 
                 data = tanks, 
                 family = Gamma(link = "log"))

summary(glm_model)

# Anova-Test (Type II) 
Anova(glm_model, type = "II")

# Post-hoc Test 
posthoc <- glht(glm_model, linfct = mcp(treatment = "Tukey"))
summary(posthoc)
cld_results <- cld(posthoc)

means <- tanks %>%
  group_by(treatment) %>%
  summarise(mean_ratio = mean(ratio_herb_det, na.rm = TRUE),
            max_ratio = max(ratio_herb_det, na.rm = TRUE)) %>%
  left_join(tibble(treatment = names(cld_results$mcletters$Letters),
                   label = cld_results$mcletters$Letters),
            by = "treatment")

# Plot
ggplot(tanks, aes(x = treatment, y = ratio_herb_det, fill = treatment)) +
  geom_boxplot() +
  geom_text(data = means, aes(x = treatment, y = max_ratio + 0.2, label = label),
            color = "black", size = 6) +
  scale_fill_manual(values = treatment_colors) +
  scale_y_continuous(limits = c(1, 4.5)) +
  labs(
    x = "Treatment",
    y = "Herbivory / Detritivory",
    fill = "Treatment"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

################################################################################
################################################################################
########################## ratio PP/Detritus ANOVA #############################
################################################################################
################################################################################

shapiro_test(tanks$ratio_herb_det)

tanks$ratio_herb_det <- biomass_phy$tot/biomass_det$`mean(biomass_mgC_m2)`

tanks$treatment <- factor(tanks$treatment,
                          levels = c("0HW", "1HW", "3HW"))

# ANONVA
anova_model <- aov(ratio_herb_det ~ treatment, data = tanks)
summary(anova_model)

# Tukey HSD Test
tukey <- glht(anova_model, linfct = mcp(treatment = "Tukey"))
tukey_cld <- cld(tukey)

means <- tanks %>%
  group_by(treatment) %>%
  summarise(mean_herb_det = mean(ratio_herb_det, na.rm = TRUE),
            max_herb_det = max(ratio_herb_det, na.rm = TRUE)) %>%
  left_join(tibble(treatment = names(tukey_cld$mcletters$Letters),
                   label = tukey_cld$mcletters$Letters),
            by = "treatment")

# 3. Plot
ggplot(tanks, aes(x = treatment, y = ratio_herb_det, fill = treatment)) +
  geom_boxplot() +
  geom_text(data = means, aes(x = treatment, y = max_herb_det + 0.2, label = label), 
            color = "black", size = 6) +
  scale_fill_manual(values = treatment_colors_alt) +
  scale_y_continuous(limits = c(0, 2500)) +
  labs(
    x = "Treatment",
    y = "Herbivory / Detritivory",
    fill = "Treatment"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


################################################################################
################################################################################
################################################################################

#### all plots together:

tanks_long <- tanks_clean %>%
  select(treatment, tot_herb) %>%
  mutate(
    metric = "Herbivory~'['*mg~C~m^{-2}~day^{-1}*']'",
    value = tot_herb
  ) %>%
  select(treatment, metric, value) %>%
  bind_rows(
    tanks %>%
      select(treatment, tot_det) %>%
      mutate(
        metric = "Detritivory~'['*mg~C~m^{-2}~day^{-1}*']'",
        value = tot_det
      ) %>%
      select(treatment, metric, value),
    
    tanks %>%
      select(treatment, ratio_herb_det) %>%
      mutate(
        metric = "'Herbivory / Detritivory'",
        value = ratio_herb_det
      ) %>%
      select(treatment, metric, value),
    
    tanks %>%
      select(treatment, herb) %>%
      mutate(
        metric = "'Herbivory / PP biomass'~'['*day^{-1}*']'",
        value = herb
      ) %>%
      select(treatment, metric, value)
  )



tanks_long <- tanks_clean %>%
  select(treatment, tot_herb) %>%
  mutate(metric = "Herbivory~'['*mg~C~m^{-2}~day^{-1}*']'", value = tot_herb) %>%
  select(treatment, metric, value) %>%
  bind_rows(
    tanks %>%
      select(treatment, tot_det) %>%
      mutate(metric = "Detritivory~'['*mg~C~m^{-2}~day^{-1}*']'", value = tot_det) %>%
      select(treatment, metric, value),
    
    tanks %>%
      select(treatment, ratio_herb_det) %>%
      mutate(metric = "'Herbivory / Detritivory'", value = ratio_herb_det) %>%
      select(treatment, metric, value),
    tanks %>%
      select(treatment, herb) %>%
      mutate(metric = "'Herbivory / PP biomass'~'['*day^{-1}*']'", value = herb) %>%
      select(treatment, metric, value)
  )%>%
  mutate(metric = factor(metric, levels = c(
    "Herbivory~'['*mg~C~m^{-2}~day^{-1}*']'",
    "Detritivory~'['*mg~C~m^{-2}~day^{-1}*']'",
    "'Herbivory / Detritivory'",
    "'Herbivory / PP biomass'~'['*day^{-1}*']'"
  )))

ggplot(tanks_long, aes(x = treatment, y = value, fill = treatment)) +
  geom_boxplot() +
  facet_wrap(
    ~ metric,
    scales = "free_y",
    ncol = 2,  # <- zwei Plots nebeneinander
    strip.position = "left",
    labeller = label_parsed
  ) +
  scale_fill_manual(values = treatment_colors) +
  labs(
    x = "Treatment",
    y = NULL,
    fill = "Treatment"
  ) +
  theme_minimal() +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    legend.position = "right",
    panel.border = element_rect(color = "grey70", fill = NA, size = 0.5)
  )
