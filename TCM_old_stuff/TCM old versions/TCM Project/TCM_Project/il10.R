library(ggplot2)
library(ggbreak)
library(dplyr)
library(tidyr)

# Create your sample data
df_raw <- read.table(text = "D17	IL.10	93.55	12.4	122.59	Osteosarcoma
1771	IL.10	132.39	16.86	50.26	Leukemia/Lymphoma
17CM98	IL.10	191.62	223.35	142.79	Melanoma
Abrams	IL.10	203.46	293.56	231.83	Osteosarcoma
Bliley	IL.10	171.31	10.82	129.65	TCC
CLL1390	IL.10	69.04	14.94	87.9	Leukemia/Lymphoma
CML-10C2	IL.10	83.14	421.87	88.75	Melanoma
CML-6M	IL.10	148.94	1271.97	163.33	Melanoma
CMT12	IL.10	319.63	396	297.42	Mammary
CMT27	IL.10	269.66	156.72	263.84	Mammary
CTAC	IL.10	63.92	45.88	75.41	Thyroid carcinoma
DEN	IL.10	98.3	94.38	97.9	Hemangiosarcoma
DH82	IL.10	107.32	355.97	200.14	Histiocytic sarcoma
Gracie	IL.10	71.46	227.66	186.79	Osteosarcoma
HMPOS	IL.10	156.91	19.65	67.59	Osteosarcoma
Jones	IL.10	42.85	44.8	75.77	Melanoma
McKinley	IL.10	100	24.14	63.35	Osteosarcoma
Nike	IL.10	186.54	4057.11	224	Histiocytic sarcoma
Moresco	IL.10	103.58	6.56	95.77	Osteosarcoma
OS2.4	IL.10	64.43	37.43	58.57	Osteosarcoma
OSA8	IL.10	116.61	4057.11	140.54	Osteosarcoma
Parks	IL.10	84.7	89.96	100	Melanoma
STSA-1	IL.10	51.93	157.2	161.73	Soft tissue sarcoma
Vogel	IL.10	27.42	100	94.95	Osteosarcoma
Yamane	IL.10	42.04	145.16	84.97	Osteosarcoma", 
                     header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Name the columns appropriately
colnames(df_raw) <- c("Cell_line", "Cytokine", "Rep1", "Rep2", "Rep3", "TumorType")

# Pivot the replicate columns into a single "Value" column
df_tidy <- df_raw %>%
  pivot_longer(cols = starts_with("Rep"), names_to = "Replicate", values_to = "Value")

# Filter for IL-10 (all rows are IL.10 in this case)
df_il10 <- df_tidy %>% filter(Cytokine == "IL.10")

# Summarize by cell line and tumor type: compute mean and SEM
summary_df <- df_il10 %>%
  group_by(Cell_line, TumorType) %>%
  summarise(
    mean_value = mean(Value, na.rm = TRUE),
    se = sd(Value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_value)) %>%
  mutate(Cell_line = factor(Cell_line, levels = unique(Cell_line)))

# Create the plot with a log-transformed y-axis
p <- ggplot(summary_df, aes(x = Cell_line, y = mean_value, fill = TumorType)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_value - se, ymax = mean_value + se), width = 0.2) +
  geom_jitter(data = df_il10, 
              aes(x = Cell_line, y = Value),
              width = 0.2, height = 0, shape = 21, size = 3,
              fill = "white", color = "black") +
  labs(title = "Mean IL-10 Levels by Cell Line", fill = "Tumor Type") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_y_log10()  # This applies the log transformation

# Display the plot
print(p)
