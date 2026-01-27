# ---- packages ----
library(lme4)
library(lmerTest)
library(broom.mixed)
library(performance)
library(ggplot2)
library(dplyr)

# ---- read + prep ----
df <- read.csv("recom_ccl3_tnfa.csv", stringsAsFactors=FALSE) %>%
  mutate(
    donor = factor(donor),
    dose  = as.numeric(dose)
  )

# ---- fit mixed model ----
mod <- lmer(tnfa ~ dose + (1 | donor), data = df)

# ---- extract stats ----
tidy_mod  <- tidy(mod, effects="fixed")
slope_row <- tidy_mod %>% filter(term=="dose")
slope_est <- round(slope_row$estimate, 2)
slope_se  <- round(slope_row$std.error, 2)
pval      <- signif(slope_row$p.value, 3)

# conditional R2
r2c <- round(r2(mod)$R2_conditional, 2)

# ---- compose just p-value and R² text ----


# ---- build the plot ----
p <- ggplot(df, aes(x = dose, y = tnfa, color = donor)) +
  geom_point(size = 3) +
  stat_smooth(
    aes(group = 1), method = "lm", se = TRUE,
    color = "black", linetype = "dashed"
  ) +
  scale_x_continuous(name = "CCL3 dose (ng/mL)") +
  scale_y_continuous(name = expression(TNF~alpha~"(pg/mL)")) +
  theme_minimal(base_size = 14) +
  theme(
    # restore legend
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) 

print(p)

ggsave(
  filename = "CCL3_TNFalpha_dose_response.png",
  plot     = p,               # your ggplot object
  width    = 6,               # in inches
  height   = 4,               # in inches
  dpi      = 600,             # high resolution
  bg       = "transparent"    # transparent background
)

