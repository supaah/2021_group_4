# Clear workspace ---------------------------------------------------------
rm(list = ls())


# Load libraries ----------------------------------------------------------
library("tidyverse")
library("broom")
library("purrr")


# Define functions --------------------------------------------------------
#source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
gravier_clean_aug <- read_tsv(file = "data/03_my_data_clean_aug.tsv")


# Wrangle data ------------------------------------------------------------
gravier_clean_aug_long <- gravier_clean_aug %>%
  select(-value) %>%
  pivot_longer(cols = contains("g"),
               names_to = "gene",
               values_to = "log2_expr_level")

gravier_clean_aug_long <- gravier_clean_aug_long %>%
  group_by(gene) %>%
  nest() %>%
  ungroup(gene)

set.seed(7)
gravier_clean_aug_long_nested <- gravier_clean_aug_long %>%
  sample_n(size = 100)


# Model data
gravier_clean_aug_long_nested <- gravier_clean_aug_long_nested %>% 
  mutate(mdl = map(data, ~glm(response ~ log2_expr_level, 
                              data = .x,
                              family = binomial(link = "logit"))))

gravier_clean_aug_long_nested <- gravier_clean_aug_long_nested %>%
  mutate(mdl_tidy = map(mdl, tidy, conf.int = TRUE)) %>%
  unnest(mdl_tidy)

gravier_clean_aug_long_nested <- gravier_clean_aug_long_nested %>%
  filter(term == "log2_expr_level")

gravier_clean_aug_long_nested <- gravier_clean_aug_long_nested %>%
  mutate(identified_as = case_when(p.value < 0.05 ~ "significant",
                                   p.value >= 0.05 ~ "non_significant"))

gravier_data_wide = gravier_clean_aug %>%
  select(response, pull(gravier_clean_aug_long_nested, gene))


pca_fit <- gravier_data_wide %>% 
  select(where(is.numeric)) %>% 
  prcomp(scale = TRUE)



# Visualise data ----------------------------------------------------------

pca_fit %>%
  augment(gravier_data_wide) %>% # add original dataset back in
  mutate(response = factor(response)) %>%
  ggplot(aes(.fittedPC1, .fittedPC2, color = response)) + 
  geom_point(size = 1.5) +
  theme_classic()+
  theme(legend.position = "bottom")+
  labs(title = "PCA plot")

# Write data --------------------------------------------------------------
write_tsv(x = gravier_data_wide,
          path = "data/04_gravier_data_wide.tsv.gz")
ggsave(filename = "results/04_PCA_plot.png", width = 16, height = 9, dpi = 72)
