################################################################################

# READ THE COMMENTS
# WORK THROUGH THE CODE SEQUENTIALLY, STEP-BY-STEP

################################################################################

###################################################
# Crop growth model for the APSIM-wheat data      #
#                                                 #
# Log-logistic growth model                       #
#                                                 #
# Shota Gugushvili                                #
#                                                 #
# 3 April 2025                                    #
###################################################

################################################################################

## Packages

library(tidyverse) # Data manipulation
library(forcats) # Factor manipulation
library(readr) # Import data
library(tidytext) # reorder_within() https://juliasilge.com/blog/reorder-within/

library(hrbrthemes) # Plotting themes
library(latex2exp) # LaTeX equations with ggplot2
library(ggplot2) # Plotting
library(bayesplot) # ggplot2-based package to plot MCMC output
library(viridis) # Colour schemes
library(ggrepel) # Repel labels
library(GGally) # For pairs plot
library(plotly) # 3-d plotting
library(ggbump) # Bump chart

library(brms) # Bayesian fitting
library(shinystan) # Shiny app for Stan
library(tidybayes) # Access and manipulate Stan output in tidy format
library(parallel) # Parallelisation
library(tictoc) # Timing

library(GDAtools) # medoids
library(factoextra) # dist

## Plotting options

# hrbrthemes requires these fonts:
# https://github.com/hrbrmstr/hrbrthemes/tree/master/inst/fonts/roboto-condensed

theme_set(hrbrthemes::theme_ipsum_rc(grid = "XY",
                                     base_size = 14,
                                     plot_title_size = 24,
                                     plot_title_margin = 12,
                                     subtitle_size = 16,
                                     subtitle_margin = 18,
                                     strip_text_size = 14,
                                     caption_size = 12,
                                     caption_margin = 12,
                                     axis_title_size = 12
                                     )
          )

bayesplot_theme_set(hrbrthemes::theme_ipsum_rc(grid = "XY"))

color_scheme_set("viridisC") # viridis schemes are colourblind-friendly

## brms option

options(mc.cores = 4,  # Use 4 cores
        brms.backend = "cmdstanr" # Backend: assumes cmdstanr is installed
        )

## Working directory

# Change accordingly on your computer

setwd("~/Documents/GitHub/conceptual_framework_ii")

wd <- getwd()
wd

## Data

# Download data from: https://doi.org/10.6084/m9.figshare.24843705
# Put data in the working directory
# Load data

data <- read.delim("SimulatedData124Envs_Biomass.txt")

# Quick check and conversion of some variables to factors; biomass scaled with 1000

glimpse(data)

data <- data %>%
  mutate(year = as_factor(year),
         loc = as_factor(loc),
         Env = as_factor(Env),
         geno = as_factor(geno)) %>%
  mutate(biomass = biomass / 1000)

################################################################################

## Subset of data

# Environments are indicated in env_classes_15may.csv

set.seed(12345678) # Seed for reproducibility

geno_indices <- sample(1:199, size = 48, replace = FALSE) # Original subset

geno_names <- levels(data$geno)[geno_indices] # Genotypes: random sample

# PCA on genotypic parameters used in APSIM
# PCA components are in geno_pc_param.csv (provided by Daniela Bustos-Korts)

pca_comps <- read_delim(file = "geno_pc_param.csv")

glimpse(pca_comps)
str(pca_comps)

pca_comps <- pca_comps %>%
  mutate(geno_copy = geno)

pca_comps$geno_copy <- ""
ix_label <- geno_indices
pca_comps$geno_copy[ix_label] <- pca_comps$geno[ix_label]

ggplot(data = pca_comps, aes(x = PC1, y = PC2, label = geno_copy)) +
  geom_label_repel(max.overlaps = 25) +
  geom_point(color = ifelse(pca_comps$geno_copy == "", "grey50", "red")) +
  xlab("PC1 (23.39%)") +
  ylab("PC2 (14.35%)")

ggsave(filename = "genotypes_pca.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 18,
       units = "cm",
       device = "png",
       bg = "white")

# Extract data for genotypes in geno_names

data_small <- data %>%
  select(Env, loc, year, geno, das, biomass) %>%
  filter(geno %in% geno_names)

# Plot data

data_small %>%
  filter(geno == geno_names[1]) %>%
  ggplot() +
  geom_line(aes(x = das, y = biomass, colour = year)) +
  facet_wrap(~loc, nrow = 2) +
  scale_colour_viridis_d(option = "plasma")

## Data subsampling

# Original APSIM data are daily
# Subsample data to bi-weekly intervals, starting from das = 15
# Scale time: dscaled = das / 10 - 1.4

subsampled_ts <- seq(from = 15, to = 240, by = 14) # 240 is more than the growing season length

data_small_subsampled <- data_small %>%
  filter(das %in% subsampled_ts) %>%
  mutate(dscaled = das / 10 - 1.4)

################################################################################

## brms fit

# Non-linear hierarchical models in brms: https://paul-buerkner.github.io/brms/articles/brms_nonlinear.html

# Set up a model formula

# Four-parameter log-logistic function
# y = asym + (Asym - asym) / (1 + (x/xmid)^(-slope) )

# Asym: right asymptote
# asym: left asymptote (taken 0)
# xmid: x-value for which the response is (Asym+asym)/2
# slope: slope parameter on x-axis

# Work with log-parameters

# Separate parameters a priori uncorrelated

formula_brms <- bf(biomass ~ mubiomass,
                   nlf(mubiomass ~ Asym / (1 + (dscaled / tmid)^(-slope))),
                   nlf(Asym ~ exp(logAsym)),
                   nlf(tmid ~ exp(logtmid)),
                   nlf(slope ~ exp(logslope)),
                   logAsym ~ 1  + (1 | Env) + (1 | geno) + (1 | geno:Env),
                   logtmid ~ 1  + (1 | Env) + (1 | geno) + (1 | geno:Env),
                   logslope ~ 1  + (1 | Env) + (1 | geno) + (1 | geno:Env),
                   nl = TRUE
                   )

# Check priors

get_prior(formula = formula_brms,
          family = brmsfamily("gaussian", link_sigma = "identity"),
          data = data_small_subsampled
          )

# Priors: log parameters

prior_brms <- prior_string("normal(0, 2)", class = "sigma") +
  prior_string("normal(2, 2)", class = "b", nlpar = "logAsym", coef = "Intercept") +
  prior_string("normal(2, 2)", class = "b", nlpar = "logtmid", coef = "Intercept") +
  prior_string("normal(2, 2)", class = "b", nlpar = "logslope", coef = "Intercept") +
  prior_string("exponential(1)", class = "sd", group = "geno", nlpar = "logAsym") +
  prior_string("exponential(1)", class = "sd", group = "geno", nlpar = "logtmid") +
  prior_string("exponential(1)", class = "sd", group = "geno", nlpar = "logslope") +
  prior_string("exponential(1)", class = "sd", group = "Env", nlpar = "logAsym") +
  prior_string("exponential(1)", class = "sd", group = "Env", nlpar = "logtmid") +
  prior_string("exponential(1)", class = "sd", group = "Env", nlpar = "logslope") +
  prior_string("exponential(1)", class = "sd", group = "geno:Env", nlpar = "logAsym") +
  prior_string("exponential(1)", class = "sd", group = "geno:Env", nlpar = "logtmid") +
  prior_string("exponential(1)", class = "sd", group = "geno:Env", nlpar = "logslope")

# Fit

# Threading in bmrs: https://cran.r-project.org/web/packages/brms/vignettes/brms_threading.html
# With 8 cores: 4 chains used and thus cores = 4 and threads = threading(2)

# Computations take days on a laptop
# tic() toc() for timing
# Instead of fitting, you can also load the brmsfit object that I saved after fitting the model
# In that case skip the tic() ... toc() part below

tic()
fit_brms <- brm(formula_brms, # Model formula
                data = data_small_subsampled, # Data
                prior = prior_brms, # Prior
                control = list(adapt_delta = 0.85, max_treedepth = 15), # NUTS parameters
                warmup = 500, # Warmup samples
                thin = 1, # No thinning
                chains = 4, # 4 chains are enough for diagnostics
                cores = 4, # 4 cores
                threads = threading(2), # Threading
                backend = "cmdstanr", # Backend
                iter = 2500 # Samples per chain
                )
toc()

## Save fit to RDS

# saveRDS(fit_brms, file = "fit_brms.RDS")

## Import the saved brmsfit object

# fit_brms <- readRDS("fit_brms.RDS")

## Summary

summary(fit_brms)

# Export summary to txt file

sink("fit_brms_summary.txt")
print(fit_brms)
sink()

## Parameter estimates

fixef_brms <- as.data.frame(fixef(fit_brms))
fixef_brms

# Creata the output directory in your working directory beforehand
# Then write csv files

write.csv(x = fixef_brms,
         file = paste0(wd, "/output/", "fixef.csv")
         )

ranef_brms_geno_logAsym <- as.data.frame(ranef(fit_brms, group = "geno")$geno[ , , "logAsym_Intercept"])
ranef_brms_geno_logtmid <- as.data.frame(ranef(fit_brms, group = "geno")$geno[ , , "logtmid_Intercept"])
ranef_brms_geno_logslope <- as.data.frame(ranef(fit_brms, group = "geno")$geno[ , , "logslope_Intercept"])

write.csv(x = ranef_brms_geno_logAsym,
          file = paste0(wd, "/output/", "ranef_geno_logAsym.csv")
          )

write.csv(x = ranef_brms_geno_logtmid,
          file = paste0(wd, "/output/", "ranef_geno_logtmid.csv")
          )

write.csv(x = ranef_brms_geno_logslope,
          file = paste0(wd, "/output/", "ranef_geno_logslope.csv")
          )

## Visualise interactions via tidybayes

# http://mjskay.github.io/tidybayes/articles/tidy-brms.html

# Get parameter names

# get_variables(fit_brms)

# Compute credible intervals for interaction parameters and pick significant interactions

logAsym_geno_Env_intervals <- fit_brms %>%
  gather_draws(`r_geno:Env__logAsym`[condition,term]) %>%
  median_qi() %>%
  filter(!((.lower < 0) & (.upper > 0))) %>%
  select(condition, .value, .lower, .upper)

write.csv(x = logAsym_geno_Env_intervals,
          file = paste0(wd, "/output/", "logAsym_geno_Env_intervals.csv")
          )

logtmid_geno_Env_intervals <- fit_brms %>%
  gather_draws(`r_geno:Env__logtmid`[condition,term]) %>%
  median_qi() %>%
  filter(!((.lower < 0) & (.upper > 0))) %>%
  select(condition, .value, .lower, .upper)

write.csv(x = logtmid_geno_Env_intervals,
          file = paste0(wd, "/output/", "logtmid_geno_Env_intervals.csv")
          )

logslope_geno_Env_intervals <- fit_brms %>%
  gather_draws(`r_geno:Env__logslope`[condition,term]) %>%
  median_qi() %>%
  filter(!((.lower < 0) & (.upper > 0))) %>%
  select(condition, .value, .lower, .upper)

write.csv(x = logslope_geno_Env_intervals,
          file = paste0(wd, "/output/", "logslope_geno_Env_intervals.csv")
          )

## Plot intervals

fit_brms %>%
  gather_draws(`sd_Env__logAsym_Intercept`,
               `sd_geno__logAsym_Intercept`,
               `sd_geno:Env__logAsym_Intercept`,
               `sd_Env__logtmid_Intercept`,
               `sd_geno__logtmid_Intercept`,
               `sd_geno:Env__logtmid_Intercept`,
               `sd_Env__logslope_Intercept`,
               `sd_geno__logslope_Intercept`,
               `sd_geno:Env__logslope_Intercept`
  ) %>%
  median_qi(.width = c(.95, .67)) %>%
  separate(.variable, into = c("group", "parameter"), sep = "__", remove = FALSE) %>%
  separate(parameter, into = c("parameter", NA), sep = "_", remove = TRUE) %>%
  mutate(parameter = as_factor(parameter)) %>%
  mutate(parameter = fct_recode(parameter,
                                "Asymptote" = "logAsym",
                                "Midpoint" = "logtmid",
                                "Slope" = "logslope")
  ) %>%
  mutate(group = as_factor(group)) %>%
  mutate(group = fct_recode(group,
                            "sd_geno" = "sd_geno",
                            "sd_trial" = "sd_Env",
                            "sd_geno:trial" = "sd_geno:Env")
  ) %>%
  mutate(group = fct_relevel(group, rev)) %>%
  ggplot(aes(y = group, x = .value, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  facet_wrap(~ parameter, ncol = 2) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

# Create folder called figures_brms in your working directory beforehand
# Then save figures

ggsave(filename = "interaction_intervals.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 18,
       units = "cm",
       device = "png",
       bg = "white")

fit_brms %>%
  gather_draws(`sd_Env__logAsym_Intercept`,
               `sd_geno__logAsym_Intercept`,
               `sd_geno:Env__logAsym_Intercept`,
               `sd_Env__logtmid_Intercept`,
               `sd_geno__logtmid_Intercept`,
               `sd_geno:Env__logtmid_Intercept`,
               `sd_Env__logslope_Intercept`,
               `sd_geno__logslope_Intercept`,
               `sd_geno:Env__logslope_Intercept`
  ) %>%
  separate(.variable, into = c("group", "parameter"), sep = "__", remove = FALSE) %>%
  separate(parameter, into = c("parameter", NA), sep = "_", remove = TRUE) %>%
  ggplot(aes(y = group, x = .value)) +
  stat_halfeye() +
  facet_wrap(~ parameter, nrow = 2) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

ggsave(filename = "interaction_densities.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 18,
       units = "cm",
       device = "png",
       bg = "white")

## Standard deviation ratios geno:Env / geno

fit_brms %>%
  spread_draws(`sd_geno__logAsym_Intercept`,
               `sd_geno:Env__logAsym_Intercept`,
               `sd_geno__logtmid_Intercept`,
               `sd_geno:Env__logtmid_Intercept`,
               `sd_geno__logslope_Intercept`,
               `sd_geno:Env__logslope_Intercept`
  ) %>%
  mutate(logAsym = `sd_geno:Env__logAsym_Intercept` / `sd_geno__logAsym_Intercept`,
         logtmid = `sd_geno:Env__logtmid_Intercept` / `sd_geno__logtmid_Intercept`,
         logslope = `sd_geno:Env__logslope_Intercept` / `sd_geno__logslope_Intercept`) %>%
  gather_draws(logAsym, logtmid, logslope) %>%
  median_qi(.width = c(.95, .67)) %>%
  ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  xlab("ratio") +
  ylab("parameter")

ggsave(filename = "sd_ratios_intervals.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 18,
       units = "cm",
       device = "png",
       bg = "white")

fit_brms %>%
  spread_draws(`sd_geno__logAsym_Intercept`,
               `sd_geno:Env__logAsym_Intercept`,
               `sd_geno__logtmid_Intercept`,
               `sd_geno:Env__logtmid_Intercept`,
               `sd_geno__logslope_Intercept`,
               `sd_geno:Env__logslope_Intercept`
  ) %>%
  mutate(Asymptote = `sd_geno:Env__logAsym_Intercept` / `sd_geno__logAsym_Intercept`,
         Midpoint = `sd_geno:Env__logtmid_Intercept` / `sd_geno__logtmid_Intercept`,
         Slope = `sd_geno:Env__logslope_Intercept` / `sd_geno__logslope_Intercept`) %>%
  gather_draws(Slope, Midpoint, Asymptote) %>%
  mutate(.variable = as_factor(.variable)) %>%
  mutate(.variable = fct_relevel(.variable, rev)) %>%
  ggplot(aes(y = .variable, x = .value)) +
  stat_halfeye() +
  xlab("ratio") +
  ylab("parameter")

ggsave(filename = "sd_ratios_densities.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 18,
       units = "cm",
       device = "png",
       bg = "white")

## Variance ratios geno:Env / geno

fit_brms %>%
  spread_draws(`sd_geno__logAsym_Intercept`,
               `sd_geno:Env__logAsym_Intercept`,
               `sd_geno__logtmid_Intercept`,
               `sd_geno:Env__logtmid_Intercept`,
               `sd_geno__logslope_Intercept`,
               `sd_geno:Env__logslope_Intercept`
  ) %>%
  mutate(logAsym = (`sd_geno:Env__logAsym_Intercept` / `sd_geno__logAsym_Intercept`)^2,
         logtmid = (`sd_geno:Env__logtmid_Intercept` / `sd_geno__logtmid_Intercept`)^2,
         logslope = (`sd_geno:Env__logslope_Intercept` / `sd_geno__logslope_Intercept`)^2) %>%
  gather_draws(logAsym, logtmid, logslope) %>%
  median_qi(.width = c(.95, .67)) %>%
  ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  xlab("ratio") +
  ylab("parameter")

ggsave(filename = "var_ratios_intervals.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 18,
       units = "cm",
       device = "png",
       bg = "white")

fit_brms %>%
  spread_draws(`sd_geno__logAsym_Intercept`,
               `sd_geno:Env__logAsym_Intercept`,
               `sd_geno__logtmid_Intercept`,
               `sd_geno:Env__logtmid_Intercept`,
               `sd_geno__logslope_Intercept`,
               `sd_geno:Env__logslope_Intercept`
  ) %>%
  mutate(logAsym = (`sd_geno:Env__logAsym_Intercept` / `sd_geno__logAsym_Intercept`)^2,
         logtmid = (`sd_geno:Env__logtmid_Intercept` / `sd_geno__logtmid_Intercept`)^2,
         logslope = (`sd_geno:Env__logslope_Intercept` / `sd_geno__logslope_Intercept`)^2) %>%
  gather_draws(logAsym, logtmid, logslope) %>%
  ggplot(aes(y = .variable, x = .value)) +
  stat_halfeye() +
  xlab("ratio") +
  ylab("parameter")

ggsave(filename = "var_ratios_densities.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 18,
       units = "cm",
       device = "png",
       bg = "white")

################################################################################

## Genotypic parameter estimates

# Compute summaries

summaries_geno_logAsym <- fit_brms %>%
  spread_draws(`r_geno__logAsym`[condition,term]) %>%
  median_qi(.width = c(.95))

summaries_geno_logtmid <- fit_brms %>%
  spread_draws(`r_geno__logtmid`[condition,term]) %>%
  median_qi(.width = c(.95))

summaries_geno_logslope <- fit_brms %>%
  spread_draws(`r_geno__logslope`[condition,term]) %>%
  median_qi(.width = c(.95))

summaries_geno_logAsym <- rename(summaries_geno_logAsym, genotype = condition)
summaries_geno_logtmid <- rename(summaries_geno_logtmid, genotype = condition)
summaries_geno_logslope <- rename(summaries_geno_logslope, genotype = condition)

## Import physiological parameters

# Import APSIM parameters (provided by Daniela Bustos-Korts)

Parameters_D1 <- read_csv("Parameters_D1.csv", col_names = TRUE, quote = "\'")

glimpse(Parameters_D1)

Parameters_D2_pars <- c("Genotype",
                        "fr_lf_sen_rate",
                        "grains_per_gram_stem",
                        "ll_modifier",
                        "max_grain_size",
                        "photop_sens",
                        "potential_grain_filling_rate",
                        "tt_floral_initiation",
                        "vern_sens")

Parameters_D2 <- Parameters_D1[, Parameters_D2_pars]

Parameters_D2 <- rename(Parameters_D2, genotype = Genotype)

Parameters_D2 <- Parameters_D2 %>%
  filter(genotype %in% geno_names)

## Compute correlations

# Placeholders

Parameters_D2_pars <- Parameters_D2_pars[-1]

logAsym_corr <- data.frame(phys_par = Parameters_D2_pars,
                           mod_par = rep("Asymptote", length(Parameters_D2_pars)),
                           corr = NA,
                           lower = NA,
                           upper = NA)

logtmid_corr <- data.frame(phys_par = Parameters_D2_pars,
                           mod_par = rep("Midpoint", length(Parameters_D2_pars)),
                           corr = NA,
                           lower = NA,
                           upper = NA)

logslope_corr <- data.frame(phys_par = Parameters_D2_pars,
                            mod_par = rep("Slope", length(Parameters_D2_pars)),
                            corr = NA,
                            lower = NA,
                            upper = NA)

# Loops

for (l in 1:length(Parameters_D2_pars)){
  cor_test_res <- cor.test(pull(Parameters_D2, Parameters_D2_pars[l]),
                           pull(summaries_geno_logAsym, "r_geno__logAsym")
                           )
  
  logAsym_corr[l, 3] <- cor_test_res$estimate
  logAsym_corr[l, 4] <- cor_test_res$conf.int[1]
  logAsym_corr[l, 5] <- cor_test_res$conf.int[2]
}

for (l in 1:length(Parameters_D2_pars)){
  cor_test_res <- cor.test(pull(Parameters_D2, Parameters_D2_pars[l]),
                           pull(summaries_geno_logtmid, "r_geno__logtmid")
  )
  
  logtmid_corr[l, 3] <- cor_test_res$estimate
  logtmid_corr[l, 4] <- cor_test_res$conf.int[1]
  logtmid_corr[l, 5] <- cor_test_res$conf.int[2]
}

for (l in 1:length(Parameters_D2_pars)){
  cor_test_res <- cor.test(pull(Parameters_D2, Parameters_D2_pars[l]),
                           pull(summaries_geno_logslope, "r_geno__logslope")
  )
  
  logslope_corr[l, 3] <- cor_test_res$estimate
  logslope_corr[l, 4] <- cor_test_res$conf.int[1]
  logslope_corr[l, 5] <- cor_test_res$conf.int[2]
}

# Results

par_corr <- rbind(logAsym_corr, logtmid_corr, logslope_corr) 

par_corr <- par_corr %>%
  mutate(phys_par = as_factor(phys_par), mod_par = as_factor(mod_par)) %>%
  mutate(mod_par = fct_relevel(mod_par, "Midpoint", after = 2))

par_corr %>%
  ggplot(aes(x = phys_par, y = corr)) + 
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  coord_flip() +
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
  facet_wrap( ~ mod_par, ncol = 2) +
  xlab("parameter") +
  ylab("correlation") 
  
ggsave(filename = "par_corr.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 18,
       units = "cm",
       device = "png",
       bg = "white")

################################################################################

## ET classification

# Compute posterior summaries

summaries_Env_logAsym <- fit_brms %>%
  spread_draws(`r_Env__logAsym`[condition,term]) %>%
  median_qi(.width = c(.95))

summaries_Env_logtmid <- fit_brms %>%
  spread_draws(`r_Env__logtmid`[condition,term]) %>%
  median_qi(.width = c(.95))

summaries_Env_logslope <- fit_brms %>%
  spread_draws(`r_Env__logslope`[condition,term]) %>%
  median_qi(.width = c(.95))

summaries_Env_logAsym <- rename(summaries_Env_logAsym, trial = condition)
summaries_Env_logtmid <- rename(summaries_Env_logtmid, trial = condition)
summaries_Env_logslope <- rename(summaries_Env_logslope, trial = condition)

# Combine parameters in a single data frame

Env_pars0 <- data.frame(logAsym = summaries_Env_logAsym$r_Env__logAsym,
                        logtmid = summaries_Env_logtmid$r_Env__logtmid,
                        logslope = summaries_Env_logslope$r_Env__logslope)

row.names(Env_pars0) <- summaries_Env_logAsym$trial

# Clustering via k-means

set.seed(5678123) # Seed for reproducibility

kmeans3 <- kmeans(scale(Env_pars0), centers = 3, iter.max = 1000, nstart = 100)

# Relabel clusters
# WARNING: this relabeling is tied to the above run of kmeans

cluster_labels <- kmeans3$cluster # Placeholder

cluster_labels[kmeans3$cluster == 1] <- 3 # Replacements
cluster_labels[kmeans3$cluster == 2] <- 1
cluster_labels[kmeans3$cluster == 3] <- 2

## Visualise clusters

# Plot pairs

loc_year <- str_split(summaries_Env_logAsym$trial, "_", simplify = TRUE)
head(loc_year)

Env_pars3 <- Env_pars0 %>%
  dplyr::rename(all_of(c(Asymptote = "logAsym", Midpoint = "logtmid", Slope = "logslope"))) %>%
  mutate(cluster = as_factor(cluster_labels), # as_factor(kmeans3$cluster),
         trial = summaries_Env_logAsym$trial,
         loc = loc_year[, 1],
         year = loc_year[, 2])

cluster_legend <- ggplot(data = Env_pars3) +
  geom_point(aes(x = Asymptote, y = Slope, colour = cluster)) +
  scale_colour_viridis_d(option = "turbo") +
  theme(legend.position = "bottom")

cluster_pairs <- ggpairs(Env_pars3,
                         mapping = aes(color = cluster),
                         columns = c("Asymptote", "Midpoint", "Slope"),
                         upper = list(continuous = "points"),
                         diag = "blank",
                         labeller = "label_parsed",
                         legend = grab_legend(cluster_legend)
                         ) +
  scale_colour_viridis_d(option = "turbo")

cluster_pairs

ggsave(filename = "Env_pars_pairs_shota.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 18,
       units = "cm",
       device = "png",
       bg = "white")

## 3-d scatterplot of clusters

Env_pars_scatter <- Env_pars0 %>%
  dplyr::rename(all_of(c(log_Mmax = "logAsym", log_tmid = "logtmid", log_r = "logslope"))) %>%
  mutate(cluster = as_factor(cluster_labels), # as_factor(kmeans3$cluster),
         trial = summaries_Env_logAsym$trial,
         loc = loc_year[, 1],
         year = loc_year[, 2]
         )

env_scatter3d <- plot_ly(Env_pars_scatter, x = ~log_Mmax, y = ~log_tmid, z = ~log_r,
                         type = "scatter3d",
                         mode = "markers+text",
                         text = ~trial,
                         color = ~cluster,
                         colors = viridis_pal(option = "H")(3))

env_scatter3d

htmlwidgets::saveWidget(as_widget(env_scatter3d),
                        "figures_brms/env_scatter3d.html",
                        title = "Clustering results")

## Boxplots of parameters

Env_pars_scatter_long <- Env_pars_scatter %>%
  pivot_longer(cols = c(log_Mmax, log_tmid, log_r), names_to = "parameter", values_to = "value")

jitter <- position_jitter(width = 0.01, height = 0, seed = 42)

Env_pars_scatter_long %>%
  ggplot(mapping = aes(loc, value, colour = loc)) +
  stat_pointinterval(colour = "black", position = position_nudge(x = 0.15)) +
  geom_point(position = jitter) +
  facet_wrap(~parameter, nrow = 3) +
  scale_colour_viridis_d(option = "plasma", name = "location") +
  xlab("location") +
  ylab("value")

Env_pars_scatter_long %>%
  ggplot(mapping = aes(cluster, value, colour = cluster)) +
  stat_pointinterval(colour = "black", position = position_nudge(x = 0.15)) +
  geom_point(position = jitter) +
  facet_wrap(~parameter, nrow = 3) +
  scale_colour_viridis_d(option = "plasma", name = "cluster") +
  xlab("cluster") +
  ylab("value")

## Bump chart

# Add location and year

Env_pars3_bump <- Env_pars3 %>%
  mutate(loc = loc_year[ , 1], year = loc_year[ , 2]) %>%
  mutate(year = as.numeric(year), cluster = as.numeric(cluster))

# Plot

Env_pars3_bump %>% 
  ggplot(aes(year, cluster)) +
  geom_point() +
  geom_bump() +
  facet_wrap(~ loc) +
  scale_y_continuous(breaks = c(1, 2, 3))

ggsave(filename = "Env_bump.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 18,
       units = "cm",
       device = "png",
       bg = "white")

################################################################################

## Shinystan for some diagnostics

# launch_shinystan(fit_brms)

## Plot some posterior intervals

mcmc_areas(fit_brms,
           pars = c("b_logAsym_Intercept", "b_logtmid_Intercept", "b_logslope_Intercept"),
           prob = 0.67, # 67% intervals
           prob_outer = 0.99, # 99%
           point_est = "median"
           ) 

ggsave(filename = "pop_effects.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 18,
       units = "cm",
       device = "png",
       bg = "white")

mcmc_areas(fit_brms,
           pars = c("sd_geno__logAsym_Intercept",
                    "sd_geno__logtmid_Intercept",
                    "sd_geno__logslope_Intercept",
                    "sd_geno:Env__logAsym_Intercept",
                    "sd_geno:Env__logtmid_Intercept",
                    "sd_geno:Env__logslope_Intercept"
           ),
           prob = 0.67, # 67% intervals
           prob_outer = 0.99, # 99%
           point_est = "median"
           ) 

ggsave(filename = "group_effects.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 18,
       units = "cm",
       device = "png",
       bg = "white")

################################################################################

## Plot some fits

## Some cleanup to save memory

# rm(data, data_small)

## Plot fit

# Make combinations of geno and Env

conditions_brms <- make_conditions(fit_brms, vars = c("geno", "Env"))

## Summarise interactions by location and genotype

# logAsym

logAsym_geno_Env_intervals2 <- logAsym_geno_Env_intervals %>%
  separate(condition, into = c("genotype", "location", "year"), sep = "_", remove = FALSE)

logAsym_geno_Env_intervals2 %>%
  group_by(location) %>%
  summarize(n = n())

logAsym_geno_Env_intervals2 %>%
  group_by(genotype) %>%
  summarize(n = n()) %>%
  ungroup()  %>%
  mutate(genotype = as.factor(genotype)) %>%
  mutate(genotype = fct_reorder(genotype, n)) %>%
  ggplot(aes(genotype, n)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "frequency",
       x = "genotype")

ggsave(filename = "logAsym_geno_interaction_frequencies.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 18,
       units = "cm",
       device = "png",
       bg = "white")

logAsym_geno_Env_intervals2 %>%
  group_by(genotype) %>%
  summarize(n = n()) %>%
  ungroup()  %>%
  mutate(genotype = as.factor(genotype)) %>%
  mutate(genotype = fct_reorder(genotype, n)) %>%
  ggplot(aes(genotype, 100 * n / 124)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "percentage",
       x = "genotype")

ggsave(filename = "logAsym_geno_interaction_frequencies2.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 18,
       units = "cm",
       device = "png",
       bg = "white")

logAsym_geno_Env_intervals2 %>%
  group_by(location, year) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(location = as.factor(location)) %>%
  ggplot(aes(year, n, fill = location)) +
  geom_col() +
  scale_fill_viridis_d(option = "plasma") +
  coord_flip() +
  labs(y = "frequency",
       x = "year") +
  theme(legend.position = "right") +
  facet_wrap(~ location)

ggsave(filename = "logAsym_year_interaction_frequencies.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 24,
       units = "cm",
       device = "png",
       bg = "white")

logAsym_geno_Env_intervals2 %>%
  group_by(location, year) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(location = as.factor(location)) %>%
  ggplot(aes(year, 100 * n / 48)) +
  # geom_col(position = "dodge") +
  geom_col() +
  # scale_fill_viridis_d(option = "plasma") +
  coord_flip() +
  labs(y = "percentage",
       x = "year") +
  # theme(legend.position = "right") +
  facet_wrap(~ location, nrow = 1)

ggsave(filename = "logAsym_year_interaction_frequencies2.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 24,
       units = "cm",
       device = "png",
       bg = "white")

# logtmid

logtmid_geno_Env_intervals2 <- logtmid_geno_Env_intervals %>%
  separate(condition, into = c("genotype", "location", "year"), sep = "_", remove = FALSE)

summary(logtmid_geno_Env_intervals2)

logtmid_geno_Env_intervals2 %>%
  group_by(location) %>%
  summarize(n = n())

logtmid_geno_Env_intervals2 %>%
  group_by(genotype) %>%
  summarize(n = n()) %>%
  ungroup  %>%
  mutate(genotype = as.factor(genotype)) %>%
  mutate(genotype = fct_reorder(genotype, n)) %>%
  ggplot(aes(genotype, n)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "frequency",
       x = "genotype")

## Combinations with interaction: g012 and g199 

# rm(data, data_small) # Some cleanup

conditions_brms2 <- conditions_brms %>%
  filter(geno %in% c("g012", "g199"))

conditions_brms2

cond_eff_brms2 <- conditional_effects(fit_brms,
                                      conditions = conditions_brms2,
                                      re_formula = NULL,
                                      method = "posterior_epred",
                                      resolution = 25
                                      )

# Plot fitted curves via ggplot: extract fitted values

cond_eff_brms_df2 <- data.frame(fitted = cond_eff_brms2$dscaled$estimate__,
                               dscaled = cond_eff_brms2$dscaled$dscaled,
                               trial = as.character(cond_eff_brms2$dscaled$Env),
                               geno = as_factor(cond_eff_brms2$dscaled$geno)
                               )

data_small_subsampled2 <- data_small_subsampled %>%
  rename(trial = Env) %>%
  filter(geno %in% c("g012", "g199"), trial %in% c("Emerald_1983", "Yanco_1983"))

cbbPalette <- c("#000000", "#D55E00")
  
cond_eff_brms_df2 %>%
  filter(trial %in% c("Emerald_1983", "Yanco_1983")) %>%
  ggplot() +
  geom_line(aes(x = dscaled, y = fitted, colour = geno)) +
  geom_point(data = data_small_subsampled2,
             aes(x = dscaled, y = biomass, colour = geno),
             # alpha = 0.5,
             show.legend = FALSE) +
  facet_wrap(~ trial, nrow = 2) +
  scale_colour_manual(values = cbbPalette) +
  theme(legend.position = "bottom") +
  xlab("time") +
  ylab("biomass")

ggsave(filename = "fits_interaction.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 18,
       units = "cm",
       device = "png",
       bg = "white")

# Check credible intervals

logAsym_geno_Env_intervals2 %>%
  filter(genotype %in% c("g012", "g199"),
         location %in% c("Emerald", "Yanco"),
         year %in% c("1983"))

# Plot full time series

data_small %>%
  filter(geno %in% c("g012", "g199"),
         loc %in% c("Emerald", "Yanco"),
         year %in% c("1983")) %>%
  ggplot() +
  geom_line(aes(x = das, y = biomass, colour = geno)) +
  facet_wrap(~ loc, nrow = 2) +
  scale_colour_viridis_d(option = "plasma", name = "genotype") +
  theme(legend.position = "bottom") +
  xlab("time") +
  ylab("biomass")

## PP checks

pp_check(fit_brms, ndraws = 50)

ggsave(filename = "ppc.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 18,
       units = "cm",
       device = "png",
       bg = "white")

################################################################################

## Some fitted curves

conditions_brms3 <- conditions_brms %>%
  filter(Env %in% c("Emerald_1985", "Merredin_1992", "Narrabri_1999", "Yanco_2006"))

conditions_brms3

cond_eff_brms3 <- conditional_effects(fit_brms,
                                      conditions = conditions_brms3,
                                      re_formula = NULL,
                                      method = "posterior_epred",
                                      resolution = 25
                                      )

# Plot fitted curves via ggplot: extract fitted values

cond_eff_brms_df3 <- data.frame(fitted = cond_eff_brms3$dscaled$estimate__,
                                dscaled = cond_eff_brms3$dscaled$dscaled,
                                trial = as.character(cond_eff_brms3$dscaled$Env),
                                geno = as_factor(cond_eff_brms3$dscaled$geno)
                                )

data_small_subsampled3 <- data_small_subsampled %>%
  rename(trial = Env) %>%
  filter(trial %in% c("Emerald_1985", "Merredin_1992", "Narrabri_1999", "Yanco_2006"))

cond_eff_brms_df3 %>%
  ggplot() +
  geom_line(aes(x = dscaled, y = fitted, colour = geno)) +
  geom_point(data = data_small_subsampled3,
             aes(x = dscaled, y = biomass, colour = geno),
             alpha = 0.5,
             show.legend = FALSE) +
  facet_wrap(~ trial, nrow = 2) +
  scale_colour_viridis_d(option = "plasma", name = "genotype") +
  theme(legend.position = "none") +
  xlab("time") +
  ylab("biomass")

ggsave(filename = "fits_4loc.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 18,
       units = "cm",
       device = "png",
       bg = "white")

################################################################################

## Visualisation of ETs

# Medoids of clusters

et_medoid_indices <- medoids(D = factoextra::get_dist(Env_pars0), cl = as.vector(cluster_labels))

et_medoids <- rownames(Env_pars0)[et_medoid_indices]
et_medoids

# Extract conditions

conditions_brms4 <- conditions_brms %>%
  filter(Env %in% et_medoids)

conditions_brms4

# Compute fits

cond_eff_brms4 <- conditional_effects(fit_brms,
                                      conditions = conditions_brms4,
                                      re_formula = NULL,
                                      method = "posterior_epred",
                                      resolution = 25
                                      )

# Plot fitted curves via ggplot: extract fitted values

cond_eff_brms_df4 <- data.frame(fitted = cond_eff_brms4$dscaled$estimate__,
                                dscaled = cond_eff_brms4$dscaled$dscaled,
                                trial = as.character(cond_eff_brms4$dscaled$Env),
                                geno = as_factor(cond_eff_brms4$dscaled$geno)
                                )

data_small_subsampled4 <- data_small_subsampled %>%
  rename(trial = Env) %>%
  filter(trial %in% et_medoids)

cond_eff_brms_df4 %>%
  ggplot() +
  geom_line(aes(x = dscaled, y = fitted, colour = geno)) +
  geom_point(data = data_small_subsampled4,
             aes(x = dscaled, y = biomass, colour = geno),
             alpha = 0.5,
             show.legend = FALSE) +
  facet_wrap(~ trial, nrow = 2) +
  scale_colour_viridis_d(option = "plasma", name = "genotype") +
  theme(legend.position = "none") +
  xlab("time") +
  ylab("biomass")

ggsave(filename = "fits_3et_medoids.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 18,
       units = "cm",
       device = "png",
       bg = "white")

data_small %>%
  filter(Env %in% et_medoids) %>%
  ggplot() +
  geom_line(aes(x = das, y = biomass, colour = geno), show.legend = FALSE) +
  facet_wrap(~ Env, nrow = 2) +
  scale_colour_viridis_d(option = "turbo", name = "genotype") +
  theme(legend.position = "none") +
  xlab("time") +
  ylab("biomass")

## ETs for medoid genotype

# Genotypic parameter estimates

geno_pars0 <- data.frame(logAsym = summaries_geno_logAsym$r_geno__logAsym,
                         logtmid = summaries_geno_logtmid$r_geno__logtmid,
                         logslope = summaries_geno_logslope$r_geno__logslope
                         )

geno_pars0

# Medoid genotype

geno_medoid_index <- medoids(D = factoextra::get_dist(geno_pars0), cl = rep(1, 48))
geno_medoid <- summaries_geno_logAsym$genotype[geno_medoid_index]
geno_medoid

# Extract conditions

conditions_brms5 <- conditions_brms %>%
  filter(Env %in% et_medoids, geno == geno_medoid)

conditions_brms5

# cond_eff_brms5 <- conditional_effects(fit_brms,
#                                       conditions = conditions_brms5,
#                                       re_formula = NULL,
#                                       method = "posterior_epred",
#                                       resolution = 25
#                                       )

cond_eff_brms5 <- conditional_effects(fit_brms,
                                      conditions = conditions_brms5,
                                      re_formula = NULL,
                                      method = "posterior_predict",
                                      resolution = 25
                                      )

# Plot fitted curves via ggplot: extract fitted values

str(cond_eff_brms5)

cond_eff_brms_df5 <- data.frame(fitted = cond_eff_brms5$dscaled$estimate__,
                                lower = cond_eff_brms5$dscaled$lower__,
                                upper = cond_eff_brms5$dscaled$upper__,
                                dscaled = cond_eff_brms5$dscaled$dscaled,
                                trial = as.character(cond_eff_brms5$dscaled$Env),
                                geno = as_factor(cond_eff_brms5$dscaled$geno)
                                )

data_small_subsampled5 <- data_small_subsampled %>%
  rename(trial = Env) %>%
  filter(trial %in% et_medoids, geno == geno_medoid)

cond_eff_brms_df5 %>%
  ggplot() +
  geom_ribbon(aes(x = dscaled, ymin = lower, ymax = upper, fill = trial),
              alpha = 0.55,
              show.legend = FALSE) +
  geom_line(aes(x = dscaled, y = fitted, linetype = trial)) +
  geom_point(data = data_small_subsampled5,
             aes(x = dscaled, y = biomass, shape = trial),
             show.legend = FALSE) +
  scale_fill_viridis_d(option = "viridis", name = "trial") +
  theme(legend.position = "bottom") +
  xlab("time") +
  ylab("biomass")

ggsave(filename = "fits_3et_medoids.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 18,
       units = "cm",
       device = "png",
       bg = "white")

data_small %>%
  filter(Env %in% et_medoids) %>%
  ggplot() +
  geom_line(aes(x = das, y = biomass, colour = geno), show.legend = FALSE) +
  facet_wrap(~ Env, nrow = 2) +
  scale_colour_viridis_d(option = "plasma", name = "genotype") +
  theme(legend.position = "none") +
  xlab("time") +
  ylab("biomass")

## Another version with genotype effects integrated out

# 2.5 and 97.5% quantiles

q_low <- function(x){
  quantile(x,
           probs = 0.025,
           na.rm = FALSE,
           names = TRUE,
           type = 7,
           digits = 2
  )
}

q_hi <- function(x){
  quantile(x,
           probs = 0.975,
           na.rm = FALSE,
           names = TRUE,
           type = 7,
           digits = 2
  )
}

# Compute predictions

cond_eff_brms_df6 <- fit_brms %>%
  add_epred_draws(newdata = expand.grid(Env = et_medoids, geno = NA, dscaled = seq(from = 0.1, to = 16.9, by = 0.7)),
                  re_formula = ~ (1 | Env),
                  allow_new_levels = TRUE,
                  sample_new_levels = "uncertainty"
  ) %>%
  summarize(across(.epred, lst(mean, q_low, q_hi), .names = "{.fn}"))

cond_eff_brms_df6 %>%
  ggplot() +
  # geom_point(data = data_small_subsampled5,
  #            aes(x = dscaled, y = biomass, colour = trial),
  #            alpha = 0.5,
  #            show.legend = FALSE) +
  geom_ribbon(aes(x = dscaled, ymin = q_low, ymax = q_hi, colour = Env), fill = "grey90", show.legend = FALSE) +
  geom_line(aes(x = dscaled, y = mean, colour = Env)) +
  scale_colour_viridis_d(option = "plasma", name = "trial") +
  theme(legend.position = "bottom") +
  xlab("time") +
  ylab("biomass")

## Location medoids with medoid genotype

# Medoids of locations

loc_medoid_indices <- medoids(D = factoextra::get_dist(Env_pars0), cl = rep(1:4, each = 31))

loc_medoids <- rownames(Env_pars0)[loc_medoid_indices]

loc_medoids

# Extract conditions

conditions_brms7 <- conditions_brms %>%
  filter(Env %in% loc_medoids, geno == geno_medoid)

conditions_brms7

# cond_eff_brms7 <- conditional_effects(fit_brms,
#                                       conditions = conditions_brms7,
#                                       re_formula = NULL,
#                                       method = "posterior_epred",
#                                       resolution = 25
#                                       )

cond_eff_brms7 <- conditional_effects(fit_brms,
                                      conditions = conditions_brms7,
                                      re_formula = NULL,
                                      method = "posterior_predict",
                                      resolution = 25
                                      )

# Plot fitted curves via ggplot: extract fitted values

cond_eff_brms_df7 <- data.frame(fitted = cond_eff_brms7$dscaled$estimate__,
                                lower = cond_eff_brms7$dscaled$lower__,
                                upper = cond_eff_brms7$dscaled$upper__,
                                dscaled = cond_eff_brms7$dscaled$dscaled,
                                trial = as.character(cond_eff_brms7$dscaled$Env),
                                geno = as_factor(cond_eff_brms7$dscaled$geno)
                                )

data_small_subsampled7 <- data_small_subsampled %>%
  rename(trial = Env) %>%
  filter(trial %in% loc_medoids, geno == geno_medoid)

cond_eff_brms_df7 %>%
  ggplot() +
  geom_line(aes(x = dscaled, y = fitted, colour = trial)) +
  geom_point(data = data_small_subsampled7,
             aes(x = dscaled, y = biomass, colour = trial),
             alpha = 0.5,
             show.legend = FALSE) +
  scale_colour_viridis_d(option = "turbo", name = "trial") +
  theme(legend.position = "bottom") +
  xlab("time") +
  ylab("biomass") +
  geom_rug(data = data_small_subsampled7, sides = "b", aes(x = dscaled), colour = "red", length = unit(0.05, "npc")) +
  scale_y_continuous(expand = c(0.12, 0.12)) +
  coord_cartesian(clip = "off")

cond_eff_brms_df7 %>%
  ggplot() +
  geom_line(aes(x = dscaled, y = fitted, linetype = trial)) +
  geom_ribbon(aes(x = dscaled, ymin = lower, ymax = upper, fill = trial),
              alpha = 0.55,
              show.legend = FALSE) +
  geom_point(data = data_small_subsampled7,
             aes(x = dscaled, y = biomass, shape = trial),
             show.legend = FALSE) +
  scale_fill_viridis_d(option = "viridis", name = "trial") +
  theme(legend.position = "bottom") +
  xlab("time") +
  ylab("biomass")

cond_eff_brms_df7 %>%
  ggplot() +
  geom_ribbon(aes(x = dscaled, ymin = lower, ymax = upper),
              alpha = 0.55,
              fill = "grey70") +
  geom_point(data = data_small_subsampled7,
             aes(x = dscaled, y = biomass),
             colour = "red") +
  geom_line(aes(x = dscaled, y = fitted)) +
  xlab("time") +
  ylab("biomass") +
  facet_wrap(~ trial)

ggsave(filename = "fits_4loc_medoids.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 18,
       units = "cm",
       device = "png",
       bg = "white")

data_small %>%
  filter(Env %in% loc_medoids) %>%
  ggplot() +
  geom_line(aes(x = das, y = biomass, colour = geno), show.legend = FALSE) +
  facet_wrap(~ Env, nrow = 2) +
  scale_colour_viridis_d(option = "turbo", name = "genotype") +
  theme(legend.position = "none") +
  xlab("days after sowing") +
  ylab("biomass") +
  geom_rug(data = data_small_subsampled7, sides = "b", aes(x = 10 * dscaled + 14), colour = "red", length = unit(0.05, "npc")) +
  scale_y_continuous(expand = c(0.12, 0.12)) +
  coord_cartesian(clip = "off")

ggsave(filename = "fits_4loc_medoids_genotypes.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 18,
       units = "cm",
       device = "png",
       bg = "white")

## Marginal effects for location medoids

cond_eff_brms_df8 <- fit_brms %>%
  add_epred_draws(newdata = expand.grid(Env = loc_medoids,
                                        geno = NA,
                                        dscaled = seq(from = 0.1, to = 16.9, by = 0.7)
                                        ),
                  re_formula = ~ (1 | Env),
                  allow_new_levels = TRUE,
                  sample_new_levels = "uncertainty"
  ) %>%
  summarize(across(.epred, lst(mean, q_low, q_hi), .names = "{.fn}"))

cond_eff_brms_df8 %>%
  ggplot() +
  geom_ribbon(aes(x = dscaled, ymin = q_low, ymax = q_hi, colour = Env), fill = "grey90", show.legend = FALSE) +
  geom_line(aes(x = dscaled, y = mean, colour = Env)) +
  scale_colour_viridis_d(option = "plasma", name = "trial") +
  theme(legend.position = "bottom") +
  xlab("time") +
  ylab("biomass")

## Location medoids with all fitted genotypes

# Extract conditions

conditions_brms9 <- conditions_brms %>%
  filter(Env %in% loc_medoids)

conditions_brms9

cond_eff_brms9 <- conditional_effects(fit_brms,
                                      conditions = conditions_brms9,
                                      re_formula = NULL,
                                      method = "posterior_epred",
                                      resolution = 25
                                      )

# Plot fitted curves via ggplot: extract fitted values

cond_eff_brms_df9 <- data.frame(fitted = cond_eff_brms9$dscaled$estimate__,
                                dscaled = cond_eff_brms9$dscaled$dscaled,
                                trial = as.character(cond_eff_brms9$dscaled$Env),
                                geno = as_factor(cond_eff_brms9$dscaled$geno)
                                )

data_small_subsampled9 <- data_small_subsampled %>%
  rename(trial = Env) %>%
  filter(trial %in% loc_medoids)

cond_eff_brms_df9 %>%
  ggplot() +
  geom_line(aes(x = dscaled, y = fitted, colour = geno)) +
  geom_point(data = data_small_subsampled9,
             aes(x = dscaled, y = biomass, colour = geno),
             alpha = 0.5,
             show.legend = FALSE) +
  scale_colour_viridis_d(option = "turbo", name = "genotype") +
  facet_wrap(~ trial) +
  theme(legend.position = "none") +
  xlab("time") +
  ylab("biomass")

ggsave(filename = "fits_4loc_medoids_allfits.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 18,
       units = "cm",
       device = "png",
       bg = "white")

data_small %>%
  filter(Env %in% loc_medoids) %>%
  ggplot() +
  geom_line(aes(x = das, y = biomass, colour = geno), show.legend = FALSE) +
  facet_wrap(~ Env, nrow = 2) +
  scale_colour_viridis_d(option = "turbo", name = "genotype") +
  theme(legend.position = "none") +
  xlab("days after sowing") +
  ylab("biomass") +
  geom_rug(data = data_small_subsampled5, sides = "b", aes(x = 10 * dscaled + 14), colour = "red", length = unit(0.05, "npc")) +
  scale_y_continuous(expand = c(0.12, 0.12)) +
  coord_cartesian(clip = "off")

ggsave(filename = "fits_4loc_medoids_genotypes.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 18,
       units = "cm",
       device = "png",
       bg = "white")

## Medoid genotype and all years per location

# Extract conditions

conditions_brms10 <- conditions_brms %>%
  filter(geno %in% geno_medoid)

conditions_brms10

cond_eff_brms10 <- conditional_effects(fit_brms,
                                       conditions = conditions_brms10,
                                       re_formula = NULL,
                                       method = "posterior_epred",
                                       resolution = 25
                                       )

# Plot fitted curves via ggplot: extract fitted values

cond_eff_brms_df10 <- data.frame(fitted = cond_eff_brms10$dscaled$estimate__,
                                 dscaled = cond_eff_brms10$dscaled$dscaled,
                                 trial = as.character(cond_eff_brms10$dscaled$Env)
                                 )

loc_year10 <- str_split(as.character(cond_eff_brms10$dscaled$Env), "_", simplify = TRUE)
head(loc_year10)

cond_eff_brms_df10 <- cond_eff_brms_df10 %>%
  mutate(location = as_factor(loc_year10[ , 1]), year = as_factor(loc_year10[ , 2]))

data_small_subsampled10 <- data_small_subsampled %>%
  rename(trial = Env, location = loc) %>%
  filter(geno %in% geno_medoid)

cond_eff_brms_df10 %>%
  ggplot() +
  geom_line(aes(x = dscaled, y = fitted, colour = year)) +
  geom_point(data = data_small_subsampled10,
             aes(x = dscaled, y = biomass, colour = year),
             alpha = 0.5,
             show.legend = FALSE) +
  scale_colour_viridis_d(option = "turbo", name = "year") +
  facet_wrap(~ location) +
  theme(legend.position = "none") +
  xlab("time") +
  ylab("biomass")

ggsave(filename = "fits_4loc_allyears.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 18,
       units = "cm",
       device = "png",
       bg = "white")

################################################################################

## Diagnostics with residuals

# Compute fitted values

newdata_fitted_values <- data_small_subsampled %>%
  filter(geno %in% geno_names[1:4],
         Env %in% c("Emerald_1985", "Merredin_1992", "Narrabri_1999", "Yanco_2006")
         ) %>%
  mutate(Env = as.character(Env), geno = as.character(geno))

fitted_values <- fitted(fit_brms, newdata = newdata_fitted_values)

fitted_values_residuals <- newdata_fitted_values %>%
  mutate(fitted = fitted_values[, "Estimate"]) %>%
  mutate(residuals = biomass - fitted)

# Plot fitted values against residuals

fitted_values_residuals %>%
  ggplot() +
  geom_point(aes(x = fitted, y = residuals)) +
  facet_wrap(geno ~ Env) +
  scale_x_continuous(n.breaks = 3)

ggsave(filename = "residuals.png",
       path = paste0(wd,"/figures_brms"),
       width = 24,
       height = 18,
       units = "cm",
       device = "png",
       bg = "white")


