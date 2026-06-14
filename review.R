library(tidyverse)
library(readxl)
library(patchwork)
scopus = read_excel("scopus.xlsx", sheet = "studies")
str(scopus)

fig4_A <- scopus |>
  arrange(Year) |>
  mutate(cumulative = row_number()) |>
  ggplot(aes(x = Year, y = cumulative)) +
  geom_line(color = "#2e7d4f", linewidth = 1.5) +
  scale_x_continuous(breaks = seq(1995, 2025, by = 2)) +
  scale_y_continuous(breaks = seq(0, 60, by = 10)) +
  labs(
    x = NULL,
    y = "Cumulative number of studies"
  ) +
  theme_classic(15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig4_A

fig4_B = scopus |>
  filter(was_seasonality_studied %in% c("No", "Yes")) |>
  separate_longer_delim(environment, delim = " & ") |>
  ggplot(aes(x = was_seasonality_studied, 
             fill = environment)) +
  geom_bar(alpha = 0.6, width = 0.6, color = "black") +
  scale_fill_brewer(palette = "Spectral") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "Seasonality",
    y = "Number of studies",
    fill = "Environment"
  ) +
  theme_classic(base_size = 15)

fig4_B

fig4_C = scopus |>
  separate_longer_delim(fish_family, delim = " & ") |>
  count(fish_family, sort = TRUE) |>
  mutate(fish_family = fct_reorder(fish_family, n)) |>
  ggplot(aes(x = n, y = fish_family)) +
  geom_bar(stat = "identity", fill = "#2e7d4f",
           color = "black", alpha = 0.7) +
  labs(
    x = "Number of studies",
    y = NULL
  ) +
  theme_classic(base_size = 15)+
  scale_x_continuous(expand = c(0,0), limits = c(0,15),
                     breaks = seq(0,15, by = 5))

fig4_C

fig4 <- ((fig4_A / fig4_B) | fig4_C) +
  plot_annotation(tag_levels = "A")

ggsave("fig4.png", plot = fig4, dpi = 300, width = 14, height = 8)



