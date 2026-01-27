#part 1- impact of treatment on cytokine
library(performance)
library(lme4)      
library(emmeans)   
library(dplyr)     
library(tidyr)  
library(ggplot2)
library(multcompView)

# read in
modz <- read.csv("modz_with_ctrls.csv", stringsAsFactors = FALSE)

#pivot
long_modz <- modz %>%
  pivot_longer(cols= starts_with("donor_"),
    names_to= "donor_raw",
    values_to= "mod_z") %>%
  mutate(donor= factor(sub("donor_(\\d+)_modz", "\\1", donor_raw)),
    cell_line= factor(treatment),
    analyte= factor(cytokine)) %>%
  select(donor, cell_line, analyte, mod_z)

# model
analytes <- unique(long_modz$analyte)
omnibus_tab <- data.frame(analyte = analytes, p = NA_real_, adj_p = NA_real_)
pairwise <- vector("list", length(analytes))
names(pairwise) <- analytes

for(a in analytes) {
  d <- filter(long_modz, analyte == a)
  m0 <- lmer(mod_z ~ (1 | donor),            data = d, REML = FALSE)
  m1 <- lmer(mod_z ~ cell_line + (1 | donor), data = d, REML = FALSE)
  
  omnibus_tab$p[omnibus_tab$analyte == a] <- anova(m0, m1)$`Pr(>Chisq)`[2]
  
  pw <- pairs(emmeans(m1, "cell_line"), adjust = "tukey")
  pairwise[[a]] <- summary(pw)}

# BH correction
omnibus_tab$adj_p <- p.adjust(omnibus_tab$p, method = "BH")
print(omnibus_tab)

#diagnostics (change name for each analyte)
m_vegf <- lmer(mod_z ~ cell_line + (1 | donor),
  data = filter(long_modz, analyte == "vegf"),
  REML = FALSE)

performance::check_model(m_vegf)

pw_vegf <- pairwise[["vegf"]]

#pull out p <0.05
sig_others <- pw_vegf %>% filter(grepl(" - Jones$", contrast), 
                                 p.value < 0.05) %>%
  mutate(other = sub(" - Jones", "", contrast), 
         other = gsub("^\\(|\\)$", "", other)) %>%
  pull(other)

#get mean +/- SE
summary_df <- long_modz %>%
  filter(analyte == "vegf") %>%
  group_by(cell_line) %>%
  summarise(mean_modz = mean(mod_z), se_modz = sd(mod_z) / sqrt(n()),
    .groups = "drop") %>%
  mutate(status = case_when(cell_line == "Jones" ~ "Jones",
      cell_line %in% sig_others ~ "Sig vs Jones",
      TRUE ~ "Not sig"))

# plot it!
ggplot(summary_df, aes(x = reorder(cell_line, mean_modz), y = mean_modz,
  color = status)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_modz - se_modz, ymax = mean_modz + se_modz), 
                width = 0.2) +
  scale_color_manual(values = c(
    "Jones"        = "red",
    "Sig vs Jones" = "blue",
    "Not sig"      = "grey")) +
  coord_flip() +
  labs(title= "VEGF: Mean Modified Z-Score by Cell Line",
    subtitle = "Red = Jones; Blue = significantly lower vs Jones (p<0.05)",
    x = "Cell Line",
    y = "Mean mod-z") +
  theme_minimal()
