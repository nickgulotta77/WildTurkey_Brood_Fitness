require(parallel)
require(remotes)
require(ggplot2)
require(knitr)
require(readr)
require(ggeffects)
require(dplyr)
require(tidybayes)
require(coda)
require(lubridate)
require(posterior)
require(hms)
require(MCMCglmm)
require(brms)
require(bayestestR)

## theme for plots
theme_turkey<- function(){
  theme_bw() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size=12),
          axis.line.x = element_line(color = "black"),
          axis.line.y = element_line(color = "black"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(size = 15),
          legend.text = element_text(size = 10),
          legend.title = element_blank(),
          #legend.position = c(0.9, 0.9),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.background = element_rect(color = "black",
                                           fill = "transparent",
                                           size = 4, linetype = "blank"))
}


# load data frame
data1<-read.csv("Final_data/daily_landcover_brood_2014-2025_FINAL.csv")

# filter so we keep first nesting attempts only 
data1<-data1%>% filter(attempt == 1) 

# filter for only adults since juveniles rarely have a successful brood
data1<-data1%>% filter(age == 1) #adults = 1; juveniles = 0

## summary data
library(dplyr)

## Summarize total days per individual
summary_total_days <- data1%>%
  group_by(ID) %>%
  summarise(
    Site        = first(Site),
    age         = first(age),
    Total_days  = max(Total_days, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Total_days))  # optional: sort by total days
# 
# # View the summary
print(summary_total_days, n = 20)


### Transform and standardized response/covariates ###
###                                                ###

# center and scale day so it's easier to fit into model
data1$day_c<-scale(data1$Brood_day, scale = T, center = T)
hist(data1$day_c)

# center and scale number of gps points so it's easier to fit into model
data1$points_c<-scale(data1$n_gps_points, scale = T, center = T)
hist(data1$points_c)

# make numeric
require(psych)
data1$Shrub<-as.numeric(data1$Shrub)
data1$Developed<-as.numeric(data1$Developed)
data1$Open<-as.numeric(data1$Open)
data1$Hardwoods<-as.numeric(data1$Hardwoods)

# check histograms for each response
hist(data1$Pine)
describe(data1$Pine)


hist(data1$Open)
describe(data1$Open)

hist(data1$Hardwoods)
describe(data1$Hardwoods)


# Transform response variable
data1$Open_cu<-(data1$Open)^(1/3)
data1$Hard_cu<-sqrt(data1$Hardwoods)
data1$Pine_cu<-(data1$Pine)^(1/3)

# Check transformations
hist(data1$Open_cu)
describe(data1$Open_cu)
qqnorm(data1$Open_cu)
qqline(data1$Open_cu)


hist(data1$Hard_cu)
describe(data1$Hard_cu)
qqnorm(data1$Hard_cu)
qqline(data1$Hard_cu)


hist(data1$Pine_cu)
describe(data1$Pine_cu)
qqnorm(data1$Pine_cu)
qqline(data1$Pine_cu)


# Center/scale transformations
data1$Open_c<-scale(data1$Open_cu, scale = T, center = T)#risk-taking
hist(data1$Open_c)

data1$Hard_c<-scale(data1$Hard_cu, scale = T, center = T)#risk-taking
hist(data1$Hard_c)

data1$Pine_c<-scale(data1$Pine_cu, scale = T, center = T)#risk-taking
hist(data1$Pine_c)



# how many individuals in data frame
unique(data1$ID) # n = 207


# Bayes stuff
# Use the cmdstanr backend for Stan because it's faster and more modern than
# the default rstan. You need to install the cmdstanr package first
# (https://mc-stan.org/cmdstanr/) and then run cmdstanr::install_cmdstan() to
# install cmdstan on your computer.
options(mc.cores = 12,  # Use 4 cores
        brms.backend = "cmdstanr")
bayes_seed <- 1234


#1. Pine ------
require(brms)
require(cmdstanr)
require(bayestestR)
pine_brm<-brm(Pine_c ~  day_c + I(day_c^2) + Site +(1|ID)+
                ar(time = Brood_day, gr = ID, p = 1),
              data= data1, family="gaussian",
              warmup=2000, iter=4000,
              chains=4, init="random",thin = 2,
              cores = 4,
              seed=bayes_seed,
              backend = "cmdstanr",
              control = list(adapt_delta = 0.99, max_treedepth=15)) 
saveRDS(pine_brm, file = "pine_daily_ARFINAL.rds")
#pine_brm<-readRDS("pine_daily_ARFINAL.rds")
summary(pine_brm)
describe_posterior(pine_brm)
plot(conditional_effects(pine_brm))

# check model fit
plot(pine_brm)
pp_check(pine_brm, type = "error_hist")
pp_check(pine_brm, ndraws = 100)



# load packages
library(bayesplot)
require(bayestestR)
require(posterior)

### Calculate Repeatability ###
###                         ###

# variance for ID
var.ID<- posterior_samples(pine_brm)$"sd_ID__Intercept"^2
mean(var.ID);HPDinterval(as.mcmc(var.ID),0.95)#this computes CrI

# residual variance
var.res<-posterior_samples(pine_brm)$"sigma"^2
mean(var.res);HPDinterval(as.mcmc(var.res),0.95)#this computes CrI
# repeatability estimate
repeatability<-var.ID/(var.ID + var.res )
describe_posterior(repeatability)



### Fitness analysis ###
###                  ###

#extract posterior samples
require(dplyr)
variables(pine_brm)
draws_df <- brms::as_draws_df(pine_brm)

# select by BLUP values
blup_cols <- grep("^r_ID\\[", names(draws_df), value = TRUE)
meta_cols <- intersect(c(".chain", ".draw", ".iteration"), names(draws_df))

# pull draws from fitted model
o1 <- draws_df %>%
  dplyr::select(dplyr::all_of(c(blup_cols, meta_cols))) %>%
  #dplyr::slice_sample(n = 2000) %>%
  dplyr::select(-dplyr::any_of(c(".chain", ".draw", ".iteration"))) 


# 1a. Hen mortality during brooding----
s_ga <- read.csv("brood_survival_ref_2014_2025.csv")

# Pull the IDs that are actually in the model (from the BLUP column names)
model_IDs <- unique(sub("^r_ID\\[([^,]+),.*$", "\\1", blup_cols))

# Filter s_ga to keep only those IDs
n_before <- nrow(s_ga)
s_ga <- s_ga %>% filter(ID %in% model_IDs)
n_after <- nrow(s_ga)

message("s_ga rows before filtering: ", n_before)
message("s_ga rows after filtering:  ", n_after)
message("Individuals dropped:        ", n_before - n_after)

# Optional sanity check: any model IDs missing from s_ga?
missing_IDs <- setdiff(model_IDs, s_ga$ID)
if (length(missing_IDs) > 0) {
  message("Warning: these model IDs are not in s_ga: ",
          paste(missing_IDs, collapse = ", "))
}


require(bayestestR)

#for loop that iterates through entire distribution of blup values
#results in 2000 unique glms that correspond to each BLUP value extracted from univariate models
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(
    as_tibble(o1),
    contains("r_ID")
  )[iter, ]
  blups_FR <- tibble(
    ID = sub("^r_ID\\[([^,]+),.*$", "\\1", colnames(blups_FR2)),
    blups_FR = as.numeric(blups_FR2))
  
  data_FR2_lm <- merge(blups_FR, s_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(hen_mort~ blups_FR, data = data_FR2_lm, family = binomial)
}

#extract coefficients from 2000 models
harv_ga_da <- as.mcmc(
  sapply(
    FR2_lm_list,
    function(x) {
      summary(x)$coefficients[2, 1]
    }
  )
)



#### Log-odd ratio for the effect of risk-taking on adult mortality----
plot(harv_ga_da) #plot posterior distribution of log-odd ratios
harv_ga_da_sur<-data.frame(harv_ga_da) #make coefficients data frame so they can viewed using function below
describe_posterior(harv_ga_da_sur) # results presented here are sign reversed to produce plots.


# Mean BLUP per individual for observed points
FR <- as_tibble(o1) %>%
  pivot_longer(
    everything(),
    names_to = "ID",
    values_to = "value"
  ) %>%
  mutate(
    ID = sub("^r_ID\\[([^,]+),.*$", "\\1", ID)
  ) %>%
  group_by(ID) %>%
  summarise(BLUP = mean(value), .groups = "drop")

# Fixed prediction grid on the BLUP scale
x_grid <- tibble(
  blups_FR = seq(
    min(FR$BLUP, na.rm = TRUE),
    max(FR$BLUP, na.rm = TRUE),
    length.out = 207
  )
)

# Predict from every glm using the same x-grid
predictions_df <- data.frame()

for (i in seq_along(FR2_lm_list)) {
  
  fit_i <- FR2_lm_list[[i]]  # this is the glm from posterior draw i
  
  pred_i <- x_grid %>%
    mutate(
      draw = i,
      prediction = predict(fit_i, newdata = x_grid, type = "response")
    )
  
  predictions_df <- rbind(predictions_df, pred_i)
}

df_summary <- predictions_df %>%
  group_by(blups_FR) %>%
  summarise(
    mean = mean(prediction, na.rm = TRUE),
    ci_lower = quantile(prediction, 0.025, na.rm = TRUE),
    ci_upper = quantile(prediction, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# observed points for each individual
points_df <- FR %>%
  left_join(s_ga, by = "ID") %>%
  mutate(
    BLUP_raw1 = BLUP *2.87+3.74,
    BLUP_raw2 = BLUP_raw1^3
  )

# make sure prediction summary is ordered
df_summary <- df_summary %>%
  mutate(
    BLUP_raw1 = blups_FR *2.87+3.74,
    BLUP_raw2 = BLUP_raw1^3
  ) %>%
  arrange(BLUP_raw2)

p33 <- ggplot(df_summary, aes(x = BLUP_raw2, y = mean)) +
  geom_ribbon(
    aes(ymin = ci_lower, ymax = ci_upper),
    fill = "#A98F64",
    linetype = 0,
    alpha = 0.25,
    linewidth = 0
  ) +
  geom_line(
    color = "#A98F64",
    linewidth = 5
  ) +
  geom_point(
    data = points_df,
    aes(x = BLUP_raw2, y =hen_mort),
    color = "#A98F64",
    alpha = 0.6,
    size = 8,
    inherit.aes = FALSE
  ) +
  theme_turkey() +
  xlab("Distance to pine landcover (m)") +
  ylab("Female Survival") +
  scale_y_continuous(n.breaks = 5) +
  theme(
    legend.title = element_text(size = 74, face = "bold"),
    legend.text = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 70, face = "bold"),
    plot.title = element_text(size = 72),
    axis.ticks.length = unit(.5, "cm"),
    axis.ticks = element_line(linewidth = 3, color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2.05),
    axis.text.y = element_text(size = 68, colour = "black"),
    axis.text.x = element_text(size = 68, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(r = 55, l = 20, b = 35, t = 35)
  )

p33




# 1a. Hen mortality during brooding----
#extact posterior samples
require(dplyr)
variables(pine_brm)
draws_df <- brms::as_draws_df(pine_brm)


# select by blup blup values
blup_cols <- grep("^r_ID\\[", names(draws_df), value = TRUE)
meta_cols <- intersect(c(".chain", ".draw", ".iteration"), names(draws_df))

# pull draws from fitted model
o1 <- draws_df %>%
  dplyr::select(dplyr::all_of(c(blup_cols, meta_cols))) %>%
  #dplyr::slice_sample(n = 2000) %>%
  dplyr::select(-dplyr::any_of(c(".chain", ".draw", ".iteration"))) 


# 1a. Hen mortality during brooding----
s_ga <- read.csv("brood_survival_ref_2014_2025.csv")

# Pull the IDs that are actually in the model (from the BLUP column names)
model_IDs <- unique(sub("^r_ID\\[([^,]+),.*$", "\\1", blup_cols))

# Filter s_ga to keep only those IDs
n_before <- nrow(s_ga)
s_ga <- s_ga %>% filter(ID %in% model_IDs)
n_after <- nrow(s_ga)

message("s_ga rows before filtering: ", n_before)
message("s_ga rows after filtering:  ", n_after)
message("Individuals dropped:        ", n_before - n_after)

# Optional sanity check: any model IDs missing from s_ga?
missing_IDs <- setdiff(model_IDs, s_ga$ID)
if (length(missing_IDs) > 0) {
  message("Warning: these model IDs are not in s_ga: ",
          paste(missing_IDs, collapse = ", "))
}


#for loop that iterates through entire distribution of blup values
#results in 2000 unique glms that correspond to each BLUP value extracted from univariate models
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(
    as_tibble(o1),
    contains("r_ID")
  )[iter, ]
  blups_FR <- tibble(
    ID = sub("^r_ID\\[([^,]+),.*$", "\\1", colnames(blups_FR2)),
    blups_FR = as.numeric(blups_FR2))
  
  data_FR2_lm <- merge(blups_FR, s_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm( brood_fate~ blups_FR, data = data_FR2_lm, family = binomial)
}

#extract coefficients from 2000 models
harv_ga_da <- as.mcmc(
  sapply(
    FR2_lm_list,
    function(x) {
      summary(x)$coefficients[2, 1]
    }
  )
)



#### Log-odd ratio for the effect of risk-taking on brood fate----
plot(harv_ga_da) #plot posterior distribution of log-odd ratios
harv_ga_da_sur<-data.frame(harv_ga_da) #make coefficients data frame so they can viewed using function below
describe_posterior(harv_ga_da_sur) # results presented here are sign reversed to produce plots.


# Mean BLUP per individual for observed points
FR <- as_tibble(o1) %>%
  pivot_longer(
    everything(),
    names_to = "ID",
    values_to = "value"
  ) %>%
  mutate(
    ID = sub("^r_ID\\[([^,]+),.*$", "\\1", ID)
  ) %>%
  group_by(ID) %>%
  summarise(BLUP = mean(value), .groups = "drop")

# Fixed prediction grid on the BLUP scale
x_grid <- tibble(
  blups_FR = seq(
    min(FR$BLUP, na.rm = TRUE),
    max(FR$BLUP, na.rm = TRUE),
    length.out = 207
  )
)

# Predict from every glm using the same x-grid
predictions_df <- data.frame()

for (i in seq_along(FR2_lm_list)) {
  
  fit_i <- FR2_lm_list[[i]]  # this is the glm from posterior draw i
  
  pred_i <- x_grid %>%
    mutate(
      draw = i,
      prediction = predict(fit_i, newdata = x_grid, type = "response")
    )
  
  predictions_df <- rbind(predictions_df, pred_i)
}

df_summary <- predictions_df %>%
  group_by(blups_FR) %>%
  summarise(
    mean = mean(prediction, na.rm = TRUE),
    ci_lower = quantile(prediction, 0.025, na.rm = TRUE),
    ci_upper = quantile(prediction, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# observed points for each individual
points_df <- FR %>%
  left_join(s_ga, by = "ID") %>%
  mutate(
    BLUP_raw1 = BLUP *2.87+3.74,
    BLUP_raw2 = BLUP_raw1^3
  )

# make sure prediction summary is ordered
df_summary <- df_summary %>%
  mutate(
    BLUP_raw1 = blups_FR *2.87+3.74,
    BLUP_raw2 = BLUP_raw1^3
  ) %>%
  arrange(BLUP_raw2)

p33 <- ggplot(df_summary, aes(x = BLUP_raw2, y = mean)) +
  geom_ribbon(
    aes(ymin = ci_lower, ymax = ci_upper),
    fill = "#647EA9",
    linetype = 0,
    alpha = 0.25,
    linewidth = 0
  ) +
  geom_line(
    color = "#647EA9",
    linewidth = 5
  ) +
  geom_point(
    data = points_df,
    aes(x = BLUP_raw2, y =hen_mort),
    color = "#647EA9",
    alpha = 0.6,
    size = 8,
    inherit.aes = FALSE
  ) +
  theme_turkey() +
  xlab("Distance to pine landcover (m)") +
  ylab("Brood Survival") +
  scale_y_continuous(n.breaks = 5) +
  theme(
    legend.title = element_text(size = 74, face = "bold"),
    legend.text = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 70, face = "bold"),
    plot.title = element_text(size = 72),
    axis.ticks.length = unit(.5, "cm"),
    axis.ticks = element_line(linewidth = 3, color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2.05),
    axis.text.y = element_text(size = 68, colour = "black"),
    axis.text.x = element_text(size = 68, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(r = 55, l = 20, b = 35, t = 35)
  )

p33




#2. Open ------
require(brms)
require(cmdstanr)
require(bayestestR)
open_brm<-brm(Open_c ~  day_c + I(day_c^2) + Site + (1|ID)+
                ar(time = Brood_day, gr = ID, p = 1),
              data= data1, family="gaussian",
              warmup=2000, iter=4000,
              chains=4, init="random",thin = 2,
              cores = 4,
              seed=bayes_seed,
              backend = "cmdstanr",
              control = list(adapt_delta = 0.99, max_treedepth=15)) 
#saveRDS(open_brm, file = "Open_daily_ARFINAL.rds")
#open_brm<-readRDS("Open_daily_ARFINAL.rds")
summary(open_brm)
describe_posterior(open_brm)
plot(conditional_effects(open_brm))

plot(open_brm)
pp_check(open_brm, type = "error_hist")
pp_check(open_brm, ndraws = 100)



# load packages
library(bayesplot)
require(bayestestR)
require(posterior)

### Calculate Repeatability ###
###                         ###

# variance for ID
var.ID<- posterior_samples(open_brm)$"sd_ID__Intercept"^2
mean(var.ID);HPDinterval(as.mcmc(var.ID),0.95)#this computes CrI
# residual variance
var.res<-posterior_samples(open_brm)$"sigma"^2
mean(var.res);HPDinterval(as.mcmc(var.res),0.95)#this computes CrI
# repeatability estimate
repeatability<-var.ID/(var.ID + var.res)
describe_posterior(repeatability)



### Fitness analysis ###
###                  ###

#extract posterior samples
require(dplyr)
variables(open_brm)
draws_df <- brms::as_draws_df(open_brm)


# select by blup values
blup_cols <- grep("^r_ID\\[", names(draws_df), value = TRUE)
meta_cols <- intersect(c(".chain", ".draw", ".iteration"), names(draws_df))

# pull draws from fitted model
o1 <- draws_df %>%
  dplyr::select(dplyr::all_of(c(blup_cols, meta_cols))) %>%
  #dplyr::slice_sample(n = 2000) %>%
  dplyr::select(-dplyr::any_of(c(".chain", ".draw", ".iteration"))) 


# 1a. Hen mortality during brooding----
s_ga <- read.csv("brood_survival_ref_2014_2025.csv")

# Pull the IDs that are actually in the model (from the BLUP column names)
model_IDs <- unique(sub("^r_ID\\[([^,]+),.*$", "\\1", blup_cols))

# Filter s_ga to keep only those IDs
n_before <- nrow(s_ga)
s_ga <- s_ga %>% filter(ID %in% model_IDs)
n_after <- nrow(s_ga)

message("s_ga rows before filtering: ", n_before)
message("s_ga rows after filtering:  ", n_after)
message("Individuals dropped:        ", n_before - n_after)

# Optional sanity check: any model IDs missing from s_ga?
missing_IDs <- setdiff(model_IDs, s_ga$ID)
if (length(missing_IDs) > 0) {
  message("Warning: these model IDs are not in s_ga: ",
          paste(missing_IDs, collapse = ", "))
}


require(bayestestR)

#for loop that iterates through entire distribution of blup values
#results in 2000 unique glms that correspond to each BLUP value extracted from univariate models
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(
    as_tibble(o1),
    contains("r_ID")
  )[iter, ]
  blups_FR <- tibble(
    ID = sub("^r_ID\\[([^,]+),.*$", "\\1", colnames(blups_FR2)),
    blups_FR = as.numeric(blups_FR2))
  
  data_FR2_lm <- merge(blups_FR, s_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(hen_mort~ blups_FR, data = data_FR2_lm, family = binomial)
}

#extract coefficients from 2000 models
harv_ga_da <- as.mcmc(
  sapply(
    FR2_lm_list,
    function(x) {
      summary(x)$coefficients[2, 1]
    }
  )
)



#### Log-odd ratio for the effect of risk-taking on adult mortality ----
plot(harv_ga_da) #plot posterior distribution of log-odd ratios
harv_ga_da_sur<-data.frame(harv_ga_da) #make coefficients data frame so they can viewed using function below
describe_posterior(harv_ga_da_sur) # results presented here are sign reversed to produce plots.


# Mean BLUP per individual for observed points
FR <- as_tibble(o1) %>%
  pivot_longer(
    everything(),
    names_to = "ID",
    values_to = "value"
  ) %>%
  mutate(
    ID = sub("^r_ID\\[([^,]+),.*$", "\\1", ID)
  ) %>%
  group_by(ID) %>%
  summarise(BLUP = mean(value), .groups = "drop")

# Fixed prediction grid on the BLUP scale
x_grid <- tibble(
  blups_FR = seq(
    min(FR$BLUP, na.rm = TRUE),
    max(FR$BLUP, na.rm = TRUE),
    length.out = 207
  )
)

# Predict from every glm using the same x-grid
predictions_df <- data.frame()

for (i in seq_along(FR2_lm_list)) {
  
  fit_i <- FR2_lm_list[[i]]  # this is the glm from posterior draw i
  
  pred_i <- x_grid %>%
    mutate(
      draw = i,
      prediction = predict(fit_i, newdata = x_grid, type = "response")
    )
  
  predictions_df <- rbind(predictions_df, pred_i)
}

df_summary <- predictions_df %>%
  group_by(blups_FR) %>%
  summarise(
    mean = mean(prediction, na.rm = TRUE),
    ci_lower = quantile(prediction, 0.025, na.rm = TRUE),
    ci_upper = quantile(prediction, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# observed points for each individual
points_df <- FR %>%
 left_join(s_ga, by = "ID") %>%
  mutate(
    BLUP_raw1 = BLUP * 2.05 + 5.5,
    BLUP_raw2 = BLUP_raw1^3
  )

# make sure prediction summary is ordered
df_summary <- df_summary %>%
  mutate(
    BLUP_raw1 = blups_FR * 2.05 +5.5,
    BLUP_raw2 = BLUP_raw1^3
  ) %>%
  arrange(BLUP_raw2)

p33 <- ggplot(df_summary, aes(x = BLUP_raw2, y = mean)) +
  geom_ribbon(
    aes(ymin = ci_lower, ymax = ci_upper),
    fill = "#A98F64",
    linetype = 0,
    alpha = 0.25,
    linewidth = 0
  ) +
  geom_line(
    color = "#A98F64",
    linewidth = 5
  ) +
  geom_point(
    data = points_df,
    aes(x = BLUP_raw2, y = hen_mort),
    color = "#A98F64",
    alpha = 0.6,
    size = 8,
    inherit.aes = FALSE
  ) +
  theme_turkey() +
  xlab("Distance to open landcover (m)") +
  ylab("Female Survival") +
  scale_y_continuous(n.breaks = 5) +
  theme(
    legend.title = element_text(size = 74, face = "bold"),
    legend.text = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 70, face = "bold"),
    plot.title = element_text(size = 72),
    axis.ticks.length = unit(.5, "cm"),
    axis.ticks = element_line(linewidth = 3, color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2.05),
    axis.text.y = element_text(size = 68, colour = "black"),
    axis.text.x = element_text(size = 68, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(r = 55, l = 20, b = 35, t = 35)
  )

p33




# 1a. Hen mortality during brooding----
#extact posterior samples
require(dplyr)
variables(open_brm)
draws_df <- brms::as_draws_df(open_brm)


# select by blup values
blup_cols <- grep("^r_ID\\[", names(draws_df), value = TRUE)
meta_cols <- intersect(c(".chain", ".draw", ".iteration"), names(draws_df))

# pull draws from fitted model
o1 <- draws_df %>%
  dplyr::select(dplyr::all_of(c(blup_cols, meta_cols))) %>%
  #dplyr::slice_sample(n = 2000) %>%
  dplyr::select(-dplyr::any_of(c(".chain", ".draw", ".iteration"))) 


# 1a. Hen mortality during brooding----
s_ga <- read.csv("brood_survival_ref_2014_2025.csv")

# Pull the IDs that are actually in the model (from the BLUP column names)
model_IDs <- unique(sub("^r_ID\\[([^,]+),.*$", "\\1", blup_cols))

# Filter s_ga to keep only those IDs
n_before <- nrow(s_ga)
s_ga <- s_ga %>% filter(ID %in% model_IDs)
n_after <- nrow(s_ga)

message("s_ga rows before filtering: ", n_before)
message("s_ga rows after filtering:  ", n_after)
message("Individuals dropped:        ", n_before - n_after)

# Optional sanity check: any model IDs missing from s_ga?
missing_IDs <- setdiff(model_IDs, s_ga$ID)
if (length(missing_IDs) > 0) {
  message("Warning: these model IDs are not in s_ga: ",
          paste(missing_IDs, collapse = ", "))
}


#for loop that iterates through entire distribution of blup values
#results in 2000 unique glms that correspond to each BLUP value extracted from univariate models
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(
    as_tibble(o1),
    contains("r_ID")
  )[iter, ]
  blups_FR <- tibble(
    ID = sub("^r_ID\\[([^,]+),.*$", "\\1", colnames(blups_FR2)),
    blups_FR = as.numeric(blups_FR2))
  
  data_FR2_lm <- merge(blups_FR, s_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm( brood_fate~ blups_FR, data = data_FR2_lm, family = binomial)
}

#extract coefficients from 2000 models
harv_ga_da <- as.mcmc(
  sapply(
    FR2_lm_list,
    function(x) {
      summary(x)$coefficients[2, 1]
    }
  )
)



#### Log-odd ratio for the effect risk-taking on brood fate----
plot(harv_ga_da) #plot posterior distribution of log-odd ratios
harv_ga_da_sur<-data.frame(harv_ga_da) #make coefficients data frame so they can viewed using function below
describe_posterior(harv_ga_da_sur) # results presented here are sign reversed to produce plots.



# Mean BLUP per individual for observed points
FR <- as_tibble(o1) %>%
  pivot_longer(
    everything(),
    names_to = "ID",
    values_to = "value"
  ) %>%
  mutate(
    ID = sub("^r_ID\\[([^,]+),.*$", "\\1", ID)
  ) %>%
  group_by(ID) %>%
  summarise(BLUP = mean(value), .groups = "drop")

# Fixed prediction grid on the BLUP scale
x_grid <- tibble(
  blups_FR = seq(
    min(FR$BLUP, na.rm = TRUE),
    max(FR$BLUP, na.rm = TRUE),
    length.out = 207
  )
)

# Predict from every glm using the same x-grid
predictions_df <- data.frame()

for (i in seq_along(FR2_lm_list)) {
  
  fit_i <- FR2_lm_list[[i]]  # this is the glm from posterior draw i
  
  pred_i <- x_grid %>%
    mutate(
      draw = i,
      prediction = predict(fit_i, newdata = x_grid, type = "response")
    )
  
  predictions_df <- rbind(predictions_df, pred_i)
}

df_summary <- predictions_df %>%
  group_by(blups_FR) %>%
  summarise(
    mean = mean(prediction, na.rm = TRUE),
    ci_lower = quantile(prediction, 0.025, na.rm = TRUE),
    ci_upper = quantile(prediction, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# observed points for each individual
points_df <- FR %>%
  left_join(s_ga, by = "ID") %>%
  mutate(
    BLUP_raw1 = BLUP * 2.05 + 5.5,
    BLUP_raw2 = BLUP_raw1^3
  )

# make sure prediction summary is ordered
df_summary <- df_summary %>%
  mutate(
    BLUP_raw1 = blups_FR * 2.05 + 5.5,
    BLUP_raw2 = BLUP_raw1^3
  ) %>%
  arrange(BLUP_raw2)

p33 <- ggplot(df_summary, aes(x = BLUP_raw2, y = mean)) +
  geom_ribbon(
    aes(ymin = ci_lower, ymax = ci_upper),
    fill = "#647EA9",
    linetype = 0,
    alpha = 0.25,
    linewidth = 0
  ) +
  geom_line(
    color = "#647EA9",
    linewidth = 5
  ) +
  geom_point(
    data = points_df,
    aes(x = BLUP_raw2, y = brood_fate),
    color = "#647EA9",
    alpha = 0.6,
    size = 8,
    inherit.aes = FALSE
  ) +
  theme_turkey() +
  xlab("Distance to open landcover (m)") +
  ylab("Brood Survival") +
  scale_y_continuous(n.breaks = 5) +
  theme(
    legend.title = element_text(size = 74, face = "bold"),
    legend.text = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 70, face = "bold"),
    plot.title = element_text(size = 72),
    axis.ticks.length = unit(.5, "cm"),
    axis.ticks = element_line(linewidth = 3, color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2.05),
    axis.text.y = element_text(size = 68, colour = "black"),
    axis.text.x = element_text(size = 68, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(r = 55, l = 20, b = 35, t = 35)
  )

p33


#3. Hardwood ------
require(brms)
require(cmdstanr)
require(bayestestR)
hardwood_brm<-brm(Hard_c ~  day_c + I(day_c^2) + Site + (1|ID)+
                ar(time = Brood_day, gr = ID, p = 1),
              data= data1, family="gaussian",
              warmup=2000, iter=4000,
              chains=4, init="random",thin = 2,
              cores = 4,
              seed=bayes_seed,
              backend = "cmdstanr",
              control = list(adapt_delta = 0.99, max_treedepth=15)) 
#saveRDS(hardwood_brm, file = "hardwood_daily_ARFINAL.rds")
#hardwood_brm<-readRDS('hardwood_daily_ARFINAL.rds')
summary(hardwood_brm)
describe_posterior(hardwood_brm)
plot(conditional_effects(hardwood_brm))

plot(hardwood_brm)
pp_check(hardwood_brm, type = "error_hist")
pp_check(hardwood_brm, ndraws = 100)



# load packages
library(bayesplot)
require(bayestestR)
require(posterior)

### Calculate Repeatability ###
###                         ###
variables(hardwood_brm)
# variance for ID
var.ID<- posterior_samples(hardwood_brm)$"sd_ID__Intercept"^2
mean(var.ID);HPDinterval(as.mcmc(var.ID),0.95)#this computes CrI
# residual variance
var.res<-posterior_samples(hardwood_brm)$"sigma"^2
mean(var.res);HPDinterval(as.mcmc(var.res),0.95)#this computes CrI
# repeatability estimate
repeatability<-var.ID/(var.ID + var.res)
describe_posterior(repeatability)



### Fitness analysis ###
###                  ###

#extract posterior samples
require(dplyr)
variables(hardwood_brm)
draws_df <- brms::as_draws_df(hardwood_brm)


# select by blup values
blup_cols <- grep("^r_ID\\[", names(draws_df), value = TRUE)
meta_cols <- intersect(c(".chain", ".draw", ".iteration"), names(draws_df))

# pull draws from fitted model
o1 <- draws_df %>%
  dplyr::select(dplyr::all_of(c(blup_cols, meta_cols))) %>%
  #dplyr::slice_sample(n = 2000) %>%
  dplyr::select(-dplyr::any_of(c(".chain", ".draw", ".iteration"))) 


# 1a. Hen mortality during brooding----
s_ga <- read.csv("brood_survival_ref_2014_2025.csv")

# Pull the IDs that are actually in the model (from the BLUP column names)
model_IDs <- unique(sub("^r_ID\\[([^,]+),.*$", "\\1", blup_cols))

# Filter s_ga to keep only those IDs
n_before <- nrow(s_ga)
s_ga <- s_ga %>% filter(ID %in% model_IDs)
n_after <- nrow(s_ga)

message("s_ga rows before filtering: ", n_before)
message("s_ga rows after filtering:  ", n_after)
message("Individuals dropped:        ", n_before - n_after)

# Optional sanity check: any model IDs missing from s_ga?
missing_IDs <- setdiff(model_IDs, s_ga$ID)
if (length(missing_IDs) > 0) {
  message("Warning: these model IDs are not in s_ga: ",
          paste(missing_IDs, collapse = ", "))
}


require(bayestestR)

#for loop that iterates through entire distribution of blup values
#results in 2000 unique glms that correspond to each BLUP value extracted from univariate models
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(
    as_tibble(o1),
    contains("r_ID")
  )[iter, ]
  blups_FR <- tibble(
    ID = sub("^r_ID\\[([^,]+),.*$", "\\1", colnames(blups_FR2)),
    blups_FR = as.numeric(blups_FR2))
  
  data_FR2_lm <- merge(blups_FR, s_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(hen_mort~ blups_FR, data = data_FR2_lm, family = binomial)
}

#extract coefficients from 2000 models
harv_ga_da <- as.mcmc(
  sapply(
    FR2_lm_list,
    function(x) {
      summary(x)$coefficients[2, 1]
    }
  )
)



#### Log-odd ratio for the effect of risk-taking on adult mortality ----
plot(harv_ga_da) #plot posterior distribution of log-odd ratios
harv_ga_da_sur<-data.frame(harv_ga_da) #make coefficients data frame so they can viewed using function below
describe_posterior(harv_ga_da_sur) # results presented here are sign reversed to produce plots.


# Mean BLUP per individual for observed points
FR <- as_tibble(o1) %>%
  pivot_longer(
    everything(),
    names_to = "ID",
    values_to = "value"
  ) %>%
  mutate(
    ID = sub("^r_ID\\[([^,]+),.*$", "\\1", ID)
  ) %>%
  group_by(ID) %>%
  summarise(BLUP = mean(value), .groups = "drop")

# Fixed prediction grid on the BLUP scale
x_grid <- tibble(
  blups_FR = seq(
    min(FR$BLUP, na.rm = TRUE),
    max(FR$BLUP, na.rm = TRUE),
    length.out = 207
  )
)

# Predict from every glm using the same x-grid
predictions_df <- data.frame()

for (i in seq_along(FR2_lm_list)) {
  
  fit_i <- FR2_lm_list[[i]]  # this is the glm from posterior draw i
  
  pred_i <- x_grid %>%
    mutate(
      draw = i,
      prediction = predict(fit_i, newdata = x_grid, type = "response")
    )
  
  predictions_df <- rbind(predictions_df, pred_i)
}

df_summary <- predictions_df %>%
  group_by(blups_FR) %>%
  summarise(
    mean = mean(prediction, na.rm = TRUE),
    ci_lower = quantile(prediction, 0.025, na.rm = TRUE),
    ci_upper = quantile(prediction, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# observed points for each individual
points_df <- FR %>%
  left_join(s_ga, by = "ID") %>%
  mutate(
    BLUP_raw1 = BLUP * 5.34 + 7.88,
    BLUP_raw2 = BLUP_raw1^2
  )

# make sure prediction summary is ordered
df_summary <- df_summary %>%
  mutate(
    BLUP_raw1 = blups_FR * 5.34 +7.88,
    BLUP_raw2 = BLUP_raw1^2
  ) %>%
  arrange(BLUP_raw2)

p33 <- ggplot(df_summary, aes(x = BLUP_raw2, y = mean)) +
  geom_ribbon(
    aes(ymin = ci_lower, ymax = ci_upper),
    fill = "#A98F64",
    linetype = 0,
    alpha = 0.25,
    linewidth = 0
  ) +
  geom_line(
    color = "#A98F64",
    linewidth = 5
  ) +
  geom_point(
    data = points_df,
    aes(x = BLUP_raw2, y = hen_mort),
    color = "#A98F64",
    alpha = 0.6,
    size = 8,
    inherit.aes = FALSE
  ) +
  theme_turkey() +
  xlab("Distance to hardwood landcover (m)") +
  ylab("Female Survival") +
  scale_y_continuous(n.breaks = 5) +
  theme(
    legend.title = element_text(size = 74, face = "bold"),
    legend.text = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 70, face = "bold"),
    plot.title = element_text(size = 72),
    axis.ticks.length = unit(.5, "cm"),
    axis.ticks = element_line(linewidth = 3, color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2.05),
    axis.text.y = element_text(size = 68, colour = "black"),
    axis.text.x = element_text(size = 68, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(r = 55, l = 20, b = 35, t = 35)
  )

p33




# 1a. Hen mortality during brooding----
#extact posterior samples
require(dplyr)
variables(hardwood_brm)
draws_df <- brms::as_draws_df(hardwood_brm)


# select by blup valeus
blup_cols <- grep("^r_ID\\[", names(draws_df), value = TRUE)
meta_cols <- intersect(c(".chain", ".draw", ".iteration"), names(draws_df))

# pull draws from fitted model
o1 <- draws_df %>%
  dplyr::select(dplyr::all_of(c(blup_cols, meta_cols))) %>%
  #dplyr::slice_sample(n = 2000) %>%
  dplyr::select(-dplyr::any_of(c(".chain", ".draw", ".iteration"))) 


# 1a. Hen mortality during brooding----
s_ga <- read.csv("brood_survival_ref_2014_2025.csv")

# Pull the IDs that are actually in the model (from the BLUP column names)
model_IDs <- unique(sub("^r_ID\\[([^,]+),.*$", "\\1", blup_cols))

# Filter s_ga to keep only those IDs
n_before <- nrow(s_ga)
s_ga <- s_ga %>% filter(ID %in% model_IDs)
n_after <- nrow(s_ga)

message("s_ga rows before filtering: ", n_before)
message("s_ga rows after filtering:  ", n_after)
message("Individuals dropped:        ", n_before - n_after)

# Optional sanity check: any model IDs missing from s_ga?
missing_IDs <- setdiff(model_IDs, s_ga$ID)
if (length(missing_IDs) > 0) {
  message("Warning: these model IDs are not in s_ga: ",
          paste(missing_IDs, collapse = ", "))
}


#for loop that iterates through entire distribution of blup values
#results in 2000 unique glms that correspond to each BLUP value extracted from univariate models
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(
    as_tibble(o1),
    contains("r_ID")
  )[iter, ]
  blups_FR <- tibble(
    ID = sub("^r_ID\\[([^,]+),.*$", "\\1", colnames(blups_FR2)),
    blups_FR = as.numeric(blups_FR2))
  
  data_FR2_lm <- merge(blups_FR, s_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm( brood_fate~ blups_FR, data = data_FR2_lm, family = binomial)
}

#extract coefficients from 2000 models
harv_ga_da <- as.mcmc(
  sapply(
    FR2_lm_list,
    function(x) {
      summary(x)$coefficients[2, 1]
    }
  )
)



#### Log-odd ratio for the effect of risk-taking on brood fate----
plot(harv_ga_da) #plot posterior distribution of log-odd ratios
harv_ga_da_sur<-data.frame(harv_ga_da) #make coefficients data frame so they can viewed using function below
describe_posterior(harv_ga_da_sur) # results presented here are sign reversed to produce plots.



# Mean BLUP per individual for observed points
FR <- as_tibble(o1) %>%
  pivot_longer(
    everything(),
    names_to = "ID",
    values_to = "value"
  ) %>%
  mutate(
    ID = sub("^r_ID\\[([^,]+),.*$", "\\1", ID)
  ) %>%
  group_by(ID) %>%
  summarise(BLUP = mean(value), .groups = "drop")

# Fixed prediction grid on the BLUP scale
x_grid <- tibble(
  blups_FR = seq(
    min(FR$BLUP, na.rm = TRUE),
    max(FR$BLUP, na.rm = TRUE),
    length.out = 207
  )
)

# Predict from every glm using the same x-grid
predictions_df <- data.frame()

for (i in seq_along(FR2_lm_list)) {
  
  fit_i <- FR2_lm_list[[i]]  # this is the glm from posterior draw i
  
  pred_i <- x_grid %>%
    mutate(
      draw = i,
      prediction = predict(fit_i, newdata = x_grid, type = "response")
    )
  
  predictions_df <- rbind(predictions_df, pred_i)
}

df_summary <- predictions_df %>%
  group_by(blups_FR) %>%
  summarise(
    mean = mean(prediction, na.rm = TRUE),
    ci_lower = quantile(prediction, 0.025, na.rm = TRUE),
    ci_upper = quantile(prediction, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# observed points for each individual
points_df <- FR %>%
  left_join(s_ga, by = "ID") %>%
  mutate(
    BLUP_raw1 = BLUP * 5.34 + 7.88,
    BLUP_raw2 = BLUP_raw1^2
  )

# make sure prediction summary is ordered
df_summary <- df_summary %>%
  mutate(
    BLUP_raw1 = blups_FR * 5.34 + 7.88,
    BLUP_raw2 = BLUP_raw1^2
  ) %>%
  arrange(BLUP_raw2)

p33 <- ggplot(df_summary, aes(x = BLUP_raw2, y = mean)) +
  geom_ribbon(
    aes(ymin = ci_lower, ymax = ci_upper),
    fill = "#647EA9",
    linetype = 0,
    alpha = 0.25,
    linewidth = 0
  ) +
  geom_line(
    color = "#647EA9",
    linewidth = 5
  ) +
  geom_point(
    data = points_df,
    aes(x = BLUP_raw2, y = brood_fate),
    color = "#647EA9",
    alpha = 0.6,
    size = 8,
    inherit.aes = FALSE
  ) +
  theme_turkey() +
  xlab("Distance to open landcover (m)") +
  ylab("Brood Survival") +
  scale_y_continuous(n.breaks = 5) +
  theme(
    legend.title = element_text(size = 74, face = "bold"),
    legend.text = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 70, face = "bold"),
    plot.title = element_text(size = 72),
    axis.ticks.length = unit(.5, "cm"),
    axis.ticks = element_line(linewidth = 3, color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2.05),
    axis.text.y = element_text(size = 68, colour = "black"),
    axis.text.x = element_text(size = 68, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(r = 55, l = 20, b = 35, t = 35)
  )

p33