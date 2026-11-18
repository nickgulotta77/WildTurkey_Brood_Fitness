
# summarize AR parameter
posterior_summary(
  open_brm,
  variable = variables(open_brm)[grepl("^ar", variables(open_brm))]
)

library(dplyr)

# Make sure data are ungrouped
data1_acf <- data1 %>%
  ungroup()

# Extract residuals from the AR model
res_ar <- as.data.frame(residuals(open_brm, summary = TRUE))

# Add residuals back to the data
data1_acf <- data1_acf %>%
  mutate(resid_open_ar = res_ar$Estimate) %>%
  arrange(ID, Brood_day)

# Standard pooled ACF plot
acf(
  data1_acf$resid_open_ar,
  lag.max = 10,
  main = "Hardwood landcover-ACF of residuals"
)
