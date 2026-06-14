# Packages ----------------------------------
carregar_pacotes <- function() {
  pacotes <- c(
    "readxl",
    "tidyverse",
    "vegan",
    "RInSp",
    "ggthemes",
    "pairwiseAdonis",
    "R2jags",
    "rjags",
    "loo",
    "adespatial",
    "purrr",
    "FSA",
    "mclust",
    "MASS",
    "ordinal",
    "ggeffects",
    "brant",
    "cowplot",
    "performance",
    "glmmTMB",
    "broom.mixed", 
    "patchwork",
    "lme4"
  )
  
  instalados <- rownames(installed.packages())
  faltando <- pacotes[!pacotes %in% instalados]
  
  if (length(faltando) > 0) {
    install.packages(faltando, dependencies = TRUE)
  }
  
  invisible(lapply(pacotes, library, character.only = TRUE))
  
  message("Todos os pacotes foram carregados com sucesso.")
}

carregar_pacotes()

# Getting data ---------------------------------
data = read_excel("Atlantirivulus.xlsx")
str(data)
diet = read_excel("Atlantirivulus.xlsx", sheet = "diet")
diet_count = diet[,2:13]

data <- data %>%
  mutate(
    rainfall_period = factor(rainfall_period, levels = c("Dry", "Wet")),
    coordinates = factor(coordinates)
  )


data_dry <- data %>% filter(rainfall_period == "Dry")
data_wet <- data %>% filter(rainfall_period == "Wet")

# Stomach repletion X Period ----------------------------------------------
data$Stomach_repletion_index = as.factor(data$Stomach_repletion_index )
levels(data$Stomach_repletion_index)
data$stomach_bin <- ifelse(data$Stomach_repletion_index == 1, 0, 1)

tab_bin <- table(data$rainfall_period, data$stomach_bin)
prop.table(tab_bin, margin = 1)


model_bin <- glmmTMB(
  stomach_bin ~ rainfall_period + (1 | coordinates),
  data = data,
  family = binomial
)

summary(model_bin)

pred_bin <- ggpredict(model_bin, terms = "rainfall_period")

## Figure 1A ----------------------------------------
fig1_A = ggplot(pred_bin, aes(x = x, y = predicted, fill = x, color = x)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.08, size = 1.8) +
  geom_point(size = 7, shape = 21, stroke = 1.8) +
  scale_fill_manual(values = c(
    "Dry" = "#C7B299",  
    "Wet" = "#B8D8B8"  
  )) +
  scale_color_manual(values = c(
    "Dry" = "#8C6D4F",  
    "Wet" = "#5E8C61"   
  )) +
  labs(
    x = "Hydrological period",
    y = "Probability of stomach with food"
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "none" 
  )

fig1_A 

# Preparing diet data -------------------------
stopifnot("individual_code" %in% names(data))
stopifnot("month" %in% names(data))

data <- data %>%
  mutate(
    rainfall_period = factor(rainfall_period, levels = c("Dry","Wet")),
    coordinates     = factor(coordinates),
    Sex             = factor(Sex, levels = c("Males","Females")),
    month           = as.character(month)  
  )

diet_df <- diet_count %>%
  mutate(individual_code = diet$individual_code) %>%
  left_join(
    data %>% dplyr::select(individual_code, rainfall_period, coordinates, Sex, month),
    by = "individual_code"
  )

stopifnot("month" %in% names(diet_df))
if(any(is.na(diet_df$month))) {
  warning("Tem NA em month após o join. Veja quais individual_code não bateram:")
  print(diet_df %>% filter(is.na(month)) %>% select(individual_code) %>% distinct())
}

diet_items <- setdiff(names(diet_df),
                      c("individual_code","rainfall_period","coordinates","Sex","month"))

diet_mat_all <- as.matrix(diet_df[, diet_items])
ni_all <- rowSums(diet_mat_all)

keep <- ni_all > 0
diet_df2 <- diet_df[keep, ]
diet_mat <- diet_mat_all[keep, , drop = FALSE]
ni <- ni_all[keep]

stopifnot(nrow(diet_mat) == nrow(diet_df2))
stopifnot(all(!is.na(diet_df2$month)))

diet_df2 <- diet_df2 %>%
  mutate(
    pool_month = interaction(coordinates, month, drop = TRUE)
  )

table(diet_df2$rainfall_period)
length(unique(diet_df2$pool_month))

# PERMANOVA, PERMDISP -----------------------
dist_bray <- vegdist(diet_mat, method = "bray")

# PERMANOVA
permanova_period_sex <- adonis2(
  dist_bray ~ rainfall_period + Sex,
  data = diet_df2,
  permutations = 999,
  strata = diet_df2$coordinates,
  by = "margin"
)

permanova_period_sex

# PERMDISP Period
bd_period <- betadisper(dist_bray, diet_df2$rainfall_period)
anova(bd_period)
permutest(bd_period, permutations = 999)

# PERMDISP Sex
bd_sex <- betadisper(dist_bray, diet_df2$Sex)
anova(bd_sex)
permutest(bd_sex, permutations = 999)

## Figure 1B -----------------------------------
diet_long <- as.data.frame(diet_mat) %>%
  mutate(id = row_number()) %>%
  bind_cols(diet_df2 %>% dplyr::select(rainfall_period, Sex)) %>%
  pivot_longer(
    cols = -c(id, rainfall_period, Sex),
    names_to = "item",
    values_to = "abundance"
  )

diet_period <- diet_long %>%
  group_by(rainfall_period, item) %>%
  summarise(total = sum(abundance), .groups = "drop") %>%
  group_by(rainfall_period) %>%
  mutate(prop = total / sum(total))


fig1_B = ggplot(diet_period, aes(x = rainfall_period, y = prop, fill = item)) +
  geom_col(width = 0.7, color = "black", alpha = 0.8) +
  scale_fill_brewer(palette = "Spectral") +
  labs(
    x = "Hydrological period",
    y = "Proportion of diet",
    fill = "Items"
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic(base_size = 16)

fig1_B

## Figure 1 Complete ----------------------
Fig_1 <- (fig1_A | fig1_B) +
  plot_annotation(tag_levels = "A")

Fig_1

ggsave("Figure_1.jpg", Fig_1, height = 5, width = 10)

# Model -----------------------------
sink("PSi_pool_month.txt")
cat("
model {
  for (i in 1:N) {
    yi[i,1:J] ~ dmulti(pi[i,1:J], ni[i])
    pi[i,1:J] ~ ddirch(alpha_pm[pm[i], 1:J])
  }

  for (p in 1:P) {
    for (j in 1:J) {
      alpha_pm[p, j] <- q[p, j] * w[p] + 0.05
    }
    q[p, 1:J] ~ ddirch(alpha[])
    w[p] ~ dunif(0.1, 30)
  }

  for (j in 1:J) {
    alpha[j] <- 1
  }

  for (i in 1:N) {
    for (j in 1:J) {
      diff.prop[i,j] <- abs(pi[i,j] - q[pm[i], j])
    }
    PS[i] <- 1 - 0.5 * sum(diff.prop[i,1:J])
    log_lik[i] <- logdensity.multi(yi[i,], pi[i,], ni[i])
  }

  mean.PS <- mean(PS[])
}
", fill = TRUE)
sink()


# Safety checks -------------------------
stopifnot(nrow(diet_mat) == nrow(diet_df2))
stopifnot("individual_code" %in% names(diet_df2))
stopifnot("coordinates" %in% names(diet_df2))
stopifnot("month" %in% names(diet_df2))
stopifnot("rainfall_period" %in% names(diet_df2))

rownames(diet_mat) <- as.character(diet_df2$individual_code)

fit_ps_model_period <- function(period_label, n_iter = 4000, n_burnin = 2000, n_thin = 5){
  
  df_p <- diet_df2 %>%
    filter(rainfall_period == period_label) %>%
    arrange(individual_code)
  
  mat_p <- diet_mat[as.character(df_p$individual_code), , drop = FALSE]
  
  stopifnot(nrow(df_p) == nrow(mat_p))
  
  ni_p <- rowSums(mat_p)
  
  keep <- ni_p > 0
  df_p  <- df_p[keep, , drop = FALSE]
  mat_p <- mat_p[keep, , drop = FALSE]
  ni_p  <- ni_p[keep]
  
  stopifnot(all(ni_p > 0))
  stopifnot(nrow(df_p) == nrow(mat_p))
  
  pm_p <- as.numeric(factor(df_p$pool_month))
  P_p  <- length(unique(pm_p))
  
  jags_data <- list(
    yi = mat_p,
    N  = nrow(mat_p),
    J  = ncol(mat_p),
    ni = ni_p,
    pm = pm_p,
    P  = P_p
  )
  
  params <- c("pi","q","w","PS","log_lik","mean.PS")
  
  inits_fun <- function() list(w = rep(1, P_p))
  
  mod <- jags(
    data = jags_data,
    inits = inits_fun,
    parameters.to.save = params,
    model.file = "PSi_pool_month.txt",
    n.chains = 3,
    n.iter = n_iter,
    n.burnin = n_burnin,
    n.thin = n_thin,
    DIC = TRUE
  )
  
  list(
    period = period_label,
    df = df_p,
    mat = mat_p,
    ni = ni_p,
    pm = pm_p,
    P = P_p,
    model = mod
  )
}


set.seed(1234)
res_dry <- fit_ps_model_period("Dry")
set.seed(1234)
res_wet <- fit_ps_model_period("Wet")

# TNW Values ----------------------------
extract_q_mean <- function(jags_mod, P, J){
  sumj <- jags_mod$BUGSoutput$summary
  qtab <- sumj[grep("^q\\[", rownames(sumj)), , drop = FALSE]
  rn   <- gsub("\\s+", "", rownames(qtab))        # "q[1, 1]" -> "q[1,1]"
  
  m     <- regexec("^q\\[(\\d+),(\\d+)\\]$", rn)
  parts <- regmatches(rn, m)
  ok    <- lengths(parts) == 3
  if(!any(ok)) stop("Não achei parâmetros q[p,j] no summary do JAGS.")
  
  p <- as.integer(sapply(parts[ok], `[`, 2))
  j <- as.integer(sapply(parts[ok], `[`, 3))
  q <- qtab[ok, "mean"]
  
  q_mat <- matrix(NA_real_, nrow = P, ncol = J)
  q_mat[cbind(p, j)] <- q
  q_mat
}

H_shannon <- function(p){
  p <- p[p > 0 & is.finite(p)]
  -sum(p * log(p))
}

tnw_from_res <- function(res, period_label){
  q_mat <- extract_q_mean(res$model, P = res$P, J = ncol(res$mat))
  
  pm_map <- tibble(
    pm = seq_along(levels(factor(res$df$pool_month))),
    pool_month = levels(factor(res$df$pool_month))
  )
  
  tibble(
    pm = seq_len(res$P),
    H_shannon = apply(q_mat, 1, H_shannon),
    period = period_label
  ) %>%
    left_join(pm_map, by = "pm")
}


TNW_df <- bind_rows(
  tnw_from_res(res_dry, "Dry"),
  tnw_from_res(res_wet, "Wet")
) %>%
  mutate(
    period = factor(period, levels = c("Dry", "Wet")),
    pool_month = as.character(pool_month)
  )

# TNW ~ Period ---------------------------------
mod_tnw <- glm(
  H_shannon ~ period,
  data = TNW_df,
  family = Gamma(link = "log")
)

summary(mod_tnw)

## Figure 2A ------------------------------------
fig2_A = TNW_df %>% 
  ggplot(aes(x = period, y = H_shannon, fill = period)) +
  geom_violin(alpha = 0.7, color = "black", width = 0.9, trim = FALSE) +
  geom_boxplot(width = 0.15,  fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.05, size = 2, alpha = 0.6) +
  scale_fill_manual(values = c(
    "Dry" = "#C7B299",
    "Wet" = "#B8D8B8"   
  )) +
  labs(
    x = "Hydrological period",
    y = "TNW"
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "none"
  )

fig2_A

# Null Model TNW ---------------------------------
simulate_null_period_matrix <- function(df_period, mat_period) {
  ni_p <- rowSums(mat_period)
  stopifnot(all(ni_p > 0))
  
  pooled_counts <- colSums(mat_period)
  prob_items <- pooled_counts / sum(pooled_counts)
  
  null_mat <- t(sapply(ni_p, function(n) {
    as.vector(rmultinom(1, size = n, prob = prob_items))
  }))
  
  colnames(null_mat) <- colnames(mat_period)
  rownames(null_mat) <- rownames(mat_period)
  
  null_mat
}

fit_ps_model_from_matrix <- function(df_p, mat_p,
                                     n_iter = 2000,
                                     n_burnin = 1000,
                                     n_thin = 4) {
  
  stopifnot(nrow(df_p) == nrow(mat_p))
  
  ni_p <- rowSums(mat_p)
  keep <- ni_p > 0
  
  df_p  <- df_p[keep, , drop = FALSE]
  mat_p <- mat_p[keep, , drop = FALSE]
  ni_p  <- ni_p[keep]
  
  stopifnot(all(ni_p > 0))
  stopifnot(nrow(df_p) == nrow(mat_p))
  
  pm_p <- as.numeric(factor(df_p$pool_month))
  P_p  <- length(unique(pm_p))
  
  jags_data <- list(
    yi = mat_p,
    N  = nrow(mat_p),
    J  = ncol(mat_p),
    ni = ni_p,
    pm = pm_p,
    P  = P_p
  )
  
  params <- c("q", "w", "PS", "mean.PS")
  
  inits_fun <- function() list(w = rep(1, P_p))
  
  mod <- jags(
    data = jags_data,
    inits = inits_fun,
    parameters.to.save = params,
    model.file = "PSi_pool_month.txt",
    n.chains = 3,
    n.iter = n_iter,
    n.burnin = n_burnin,
    n.thin = n_thin,
    DIC = TRUE
  )
  
  list(
    df = df_p,
    mat = mat_p,
    ni = ni_p,
    pm = pm_p,
    P = P_p,
    model = mod
  )
}

extract_q_mean <- function(jags_mod, P, J) {
  sumj <- jags_mod$BUGSoutput$summary
  qtab <- sumj[grep("^q\\[", rownames(sumj)), , drop = FALSE]
  rn   <- gsub("\\s+", "", rownames(qtab))
  
  m     <- regexec("^q\\[(\\d+),(\\d+)\\]$", rn)
  parts <- regmatches(rn, m)
  ok    <- lengths(parts) == 3
  if (!any(ok)) stop("Nao achei parametros q[p,j] no summary do JAGS.")
  
  p <- as.integer(sapply(parts[ok], `[`, 2))
  j <- as.integer(sapply(parts[ok], `[`, 3))
  q <- qtab[ok, "mean"]
  
  q_mat <- matrix(NA_real_, nrow = P, ncol = J)
  q_mat[cbind(p, j)] <- q
  q_mat
}


H_shannon <- function(p) {
  p <- p[p > 0 & is.finite(p)]
  -sum(p * log(p))
}

tnw_from_fit <- function(fit_obj, period_label) {
  q_mat <- extract_q_mean(fit_obj$model, P = fit_obj$P, J = ncol(fit_obj$mat))
  
  pm_levels <- levels(factor(fit_obj$df$pool_month))
  pm_map <- tibble(
    pm = seq_along(pm_levels),
    pool_month = pm_levels
  )
  
  tibble(
    pm = seq_len(fit_obj$P),
    H_shannon = apply(q_mat, 1, H_shannon),
    period = period_label
  ) %>%
    left_join(pm_map, by = "pm")
}


df_dry <- diet_df2 %>%
  filter(rainfall_period == "Dry") %>%
  arrange(individual_code)

mat_dry <- diet_mat[as.character(df_dry$individual_code), , drop = FALSE]

df_wet <- diet_df2 %>%
  filter(rainfall_period == "Wet") %>%
  arrange(individual_code)

mat_wet <- diet_mat[as.character(df_wet$individual_code), , drop = FALSE]

stopifnot(nrow(df_dry) == nrow(mat_dry))
stopifnot(nrow(df_wet) == nrow(mat_wet))

set.seed(1234)
fit_obs_dry <- fit_ps_model_from_matrix(df_dry, mat_dry)

set.seed(1234)
fit_obs_wet <- fit_ps_model_from_matrix(df_wet, mat_wet)

TNW_obs <- bind_rows(
  tnw_from_fit(fit_obs_dry, "Dry"),
  tnw_from_fit(fit_obs_wet, "Wet")
) %>%
  mutate(
    period = factor(period, levels = c("Dry", "Wet")),
    type = "Observed"
  )

run_null_rep_period <- function(df_p, mat_p, period_label, rep_id,
                                n_iter = 2000,
                                n_burnin = 1000,
                                n_thin = 4) {
  
  null_mat <- simulate_null_period_matrix(df_p, mat_p)
  
  fit_null <- fit_ps_model_from_matrix(
    df_p = df_p,
    mat_p = null_mat,
    n_iter = n_iter,
    n_burnin = n_burnin,
    n_thin = n_thin
  )
  
  tnw_from_fit(fit_null, period_label) %>%
    mutate(
      rep = rep_id,
      type = "Null"
    )
}



n_reps <- 100

set.seed(2026)

TNW_null_dry <- map_dfr(1:n_reps, function(i) {
  message("Dry null replicate: ", i)
  run_null_rep_period(
    df_p = df_dry,
    mat_p = mat_dry,
    period_label = "Dry",
    rep_id = i,
    n_iter = 2000,
    n_burnin = 1000,
    n_thin = 4
  )
})


set.seed(2027)

TNW_null_wet <- map_dfr(1:n_reps, function(i) {
  message("Wet null replicate: ", i)
  run_null_rep_period(
    df_p = df_wet,
    mat_p = mat_wet,
    period_label = "Wet",
    rep_id = i,
    n_iter = 2000,
    n_burnin = 1000,
    n_thin = 4
  )
})


TNW_null <- bind_rows(TNW_null_dry, TNW_null_wet) %>%
  mutate(
    period = factor(period, levels = c("Dry", "Wet"))
  )


TNW_obs_period <- TNW_obs %>%
  group_by(period) %>%
  summarise(mean_H = mean(H_shannon), .groups = "drop")

TNW_null_period <- TNW_null %>%
  group_by(period, rep) %>%
  summarise(mean_H = mean(H_shannon), .groups = "drop")


TNW_null_summary <- TNW_null_period %>%
  group_by(period) %>%
  summarise(
    null_mean = mean(mean_H),
    null_low  = quantile(mean_H, 0.025),
    null_high = quantile(mean_H, 0.975),
    .groups = "drop"
  )

TNW_obs_summary <- TNW_obs %>%
  group_by(period) %>%
  summarise(
    obs_mean = mean(H_shannon),
    obs_se   = sd(H_shannon) / sqrt(n()),
    .groups = "drop"
  )

TNW_null_summary <- TNW_null_summary %>%
  mutate(xpos = c(1, 2))

TNW_obs_summary <- TNW_obs_summary %>%
  mutate(xpos = c(1, 2))

n <- n_reps  

TNW_test <- TNW_null_period %>%
  left_join(TNW_obs_period, by = "period", suffix = c("_null", "_obs")) %>%
  group_by(period) %>%
  summarise(
    p_upper = (sum(mean_H_null >= mean_H_obs) + 1) / (n + 1),
    p_lower = (sum(mean_H_null <= mean_H_obs) + 1) / (n + 1),
    p_two_tailed = min(1, 2 * pmin(p_upper, p_lower)),
    .groups = "drop"
  )

TNW_test

## Figure 2B -----------------------------------------------
fig2_B <- ggplot() +
  geom_errorbar(
    data = TNW_null_summary,
    aes(x = xpos, ymin = null_low, ymax = null_high),
    width = 0.08,
    linewidth = 1.4,
    color = "grey50"
  ) +
  geom_point(
    data = TNW_null_summary,
    aes(x = xpos, y = null_mean),
    shape = 21,
    size = 6.5,
    fill = "grey70",
    color = "grey50",
    stroke = 1
  ) +
  geom_errorbar(
    data = TNW_obs_summary,
    aes(
      x = xpos + 0.12,
      ymin = obs_mean - obs_se,
      ymax = obs_mean + obs_se,
      color = period
    ),
    width = 0.06,
    linewidth = 1.3
  ) +
  geom_point(
    data = TNW_obs_summary,
    aes(
      x = xpos + 0.12,
      y = obs_mean,
      fill = period,
      color = period
    ),
    shape = 21,
    size = 6.5,
    stroke = 1
  ) +
  
  scale_color_manual(values = c(
    "Dry" = "#8C6D4F",
    "Wet" = "#5E8C61"
  )) +
  
  scale_fill_manual(values = c(
    "Dry" = "#C7B299",
    "Wet" = "#B8D8B8"
  )) +
  
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Dry", "Wet"),
    expand = c(0.1, 0.1)
  ) +
  
  labs(
    x = "Hydrological period",
    y = "Mean TNW"
  ) +
  
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none"
  )

fig2_B

## Figure 2 Complete ---------------------------
Fig_2 <- (fig2_A | fig2_B) +
  plot_annotation(tag_levels = "A")

Fig_2

ggsave("Figure_2.jpg", Fig_2, height = 5, width = 10)

# PSi Values --------------------
get_PSi_means <- function(mod){
  sumj <- mod$BUGSoutput$summary
  idx  <- grep("^PS\\[", rownames(sumj))
  sumj[idx, "mean"]
}

PSi_dry <- get_PSi_means(res_dry$model)
PSi_wet <- get_PSi_means(res_wet$model)

stopifnot(length(PSi_dry) == nrow(res_dry$df))
stopifnot(length(PSi_wet) == nrow(res_wet$df))

PSi_df <- bind_rows(
  tibble(individual_code = res_dry$df$individual_code, period = "Dry", PSi = PSi_dry),
  tibble(individual_code = res_wet$df$individual_code, period = "Wet", PSi = PSi_wet)
) %>%
  mutate(period = factor(period, levels = c("Dry","Wet")))

IS_by_period <- PSi_df %>%
  group_by(period) %>%
  summarise(IS = mean(PSi, na.rm = TRUE), n = n(), .groups = "drop")

IS_by_period

PSi_df2 <- PSi_df %>%
  mutate(individual_code = as.integer(individual_code))

data2 <- data %>%
  mutate(individual_code = as.integer(individual_code))

vars_to_get <- c("pool_type","coordinates","month",
                 "Sex")

PSi_df_enriched <- PSi_df2 %>%
  left_join(
    data2 %>% dplyr::select(individual_code, any_of(vars_to_get)),
    by = "individual_code"
  )

PSi_df_enriched <- PSi_df_enriched %>%
  mutate(
    pool_month = interaction(coordinates, month, drop = TRUE),
    period = factor(period, levels = c("Dry","Wet"))
  )


PSi_df_enriched

## PSi ~ Period -----------------------
mod_psi <- glmmTMB(
  PSi ~ period + Sex + (1 | pool_month),
  data = PSi_df_enriched,
  family = beta_family(link = "logit")
)

summary(mod_psi)

## Figure 3A -------------------------
fig3_A = PSi_df_enriched %>% 
  ggplot(aes(x = period, y = PSi, fill = period)) +
  geom_violin(alpha = 0.7, color = "black", width = 0.9, trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.05, size = 2, alpha = 0.6) +
  scale_fill_manual(values = c(
    "Dry" = "#C7B299",
    "Wet" = "#B8D8B8"   
  )) +
  labs(
    x = "Hydrological period",
    y = expression(PS[i])
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "none"
  )

fig3_A

# Null Model PSi --------------------------------------
get_PSi_means_from_fit <- function(fit_obj, period_label) {
  sumj <- fit_obj$model$BUGSoutput$summary
  idx  <- grep("^PS\\[", rownames(sumj))
  psi_means <- sumj[idx, "mean"]
  
  stopifnot(length(psi_means) == nrow(fit_obj$df))
  
  tibble(
    individual_code = fit_obj$df$individual_code,
    period = period_label,
    PSi = psi_means
  )
}


PSi_obs_period <- PSi_df %>%
  group_by(period) %>%
  summarise(
    mean_PSi = mean(PSi),
    se_PSi   = sd(PSi) / sqrt(n()),
    .groups = "drop"
  )

run_null_rep_period_PSi <- function(df_p, mat_p, period_label, rep_id,
                                    n_iter = 2000,
                                    n_burnin = 1000,
                                    n_thin = 4) {
  
  null_mat <- simulate_null_period_matrix(df_p, mat_p)
  
  fit_null <- fit_ps_model_from_matrix(
    df_p = df_p,
    mat_p = null_mat,
    n_iter = n_iter,
    n_burnin = n_burnin,
    n_thin = n_thin
  )
  
  get_PSi_means_from_fit(fit_null, period_label) %>%
    mutate(
      rep = rep_id,
      type = "Null"
    )
}


n_reps <- 100

set.seed(3030)
PSi_null_dry <- purrr::map_dfr(1:n_reps, function(i) {
  message("Dry PSi null replicate: ", i)
  run_null_rep_period_PSi(
    df_p = df_dry,
    mat_p = mat_dry,
    period_label = "Dry",
    rep_id = i,
    n_iter = 2000,
    n_burnin = 1000,
    n_thin = 4
  )
})


set.seed(4040)
PSi_null_wet <- purrr::map_dfr(1:n_reps, function(i) {
  message("Wet PSi null replicate: ", i)
  run_null_rep_period_PSi(
    df_p = df_wet,
    mat_p = mat_wet,
    period_label = "Wet",
    rep_id = i,
    n_iter = 2000,
    n_burnin = 1000,
    n_thin = 4
  )
})

PSi_null <- bind_rows(PSi_null_dry, PSi_null_wet) %>%
  mutate(
    period = factor(period, levels = c("Dry", "Wet"))
  )

PSi_null_period <- PSi_null %>%
  group_by(period, rep) %>%
  summarise(
    mean_PSi = mean(PSi),
    .groups = "drop"
  )

PSi_null_summary <- PSi_null_period %>%
  group_by(period) %>%
  summarise(
    null_mean = mean(mean_PSi),
    null_low  = quantile(mean_PSi, 0.025),
    null_high = quantile(mean_PSi, 0.975),
    .groups = "drop"
  ) %>%
  mutate(xpos = c(1, 2))


PSi_obs_summary <- PSi_obs_period %>%
  rename(
    obs_mean = mean_PSi,
    obs_se   = se_PSi
  ) %>%
  mutate(xpos = c(1, 2))


n <- n_reps

PSi_test <- PSi_null_period %>%
  left_join(PSi_obs_period, by = "period", suffix = c("_null", "_obs")) %>%
  group_by(period) %>%
  summarise(
    p_upper = (sum(mean_PSi_null >= mean_PSi_obs) + 1) / (n + 1),
    p_lower = (sum(mean_PSi_null <= mean_PSi_obs) + 1) / (n + 1),
    p_two_tailed = min(1, 2 * pmin(p_upper, p_lower)),
    .groups = "drop"
  )

PSi_test

## Figure 3B ---------------------------------
fig3_B = ggplot() +
  geom_errorbar(
    data = PSi_null_summary,
    aes(x = xpos, ymin = null_low, ymax = null_high),
    width = 0.08,
    linewidth = 1.4,
    color = "grey50"
  ) +
  geom_point(
    data = PSi_null_summary,
    aes(x = xpos, y = null_mean),
    shape = 21,
    size = 6.5,
    fill = "grey70",
    color = "grey50",
    stroke = 1
  ) +
  geom_errorbar(
    data = PSi_obs_summary,
    aes(
      x = xpos + 0.12,
      ymin = obs_mean - obs_se,
      ymax = obs_mean + obs_se,
      color = period
    ),
    width = 0.06,
    linewidth = 1.3
  ) +
  geom_point(
    data = PSi_obs_summary,
    aes(
      x = xpos + 0.12,
      y = obs_mean,
      fill = period,
      color = period
    ),
    shape = 21,
    size = 6.5,
    stroke = 1
  ) +
  scale_color_manual(values = c(
    "Dry" = "#8C6D4F",
    "Wet" = "#5E8C61"
  )) +
  scale_fill_manual(values = c(
    "Dry" = "#C7B299",
    "Wet" = "#B8D8B8"
  )) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Dry", "Wet"),
    expand = c(0.1, 0.1)
  ) +
  scale_y_continuous(limits = c(0.55, 0.75))+
  labs(
    x = "Hydrological period",
    y = expression("Mean " * PS[i])
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none"
  )

fig3_B

## Figure 3 Complete ------------------------------------
Fig_3 <- (fig3_A | fig3_B) +
  plot_annotation(tag_levels = "A")

Fig_3

ggsave("Figure_3.jpg", Fig_3, height = 5, width = 10)

# GSI ~ PSi -------------------------------
data2 <- data %>%
  mutate(
    Weight_g         = readr::parse_number(gsub(",", ".", Weight_g)),
    Gonadal_weight_g = readr::parse_number(gsub(",", ".", Gonadal_weight_g))
  )


data2 <- data %>%
  mutate(
    individual_code = as.integer(individual_code),
    Weight_g = readr::parse_number(gsub(",", ".", Weight_g)),
    Gonadal_weight_g = readr::parse_number(gsub(",", ".", Gonadal_weight_g))
  )

PSi_df_gsi <- PSi_df_enriched %>%
  mutate(individual_code = as.integer(individual_code)) %>%
  left_join(
    data2 %>%
      dplyr::select(individual_code, Weight_g, Gonadal_weight_g, Sex, pool_type),
    by = "individual_code",
    suffix = c("", "_new")
  ) %>%
  mutate(
    Sex = factor(dplyr::coalesce(Sex, Sex_new)),
    pool_type = factor(pool_type),
    IGS = Gonadal_weight_g / Weight_g
  ) %>%
  dplyr::select(-Sex_new)


PSi_df_gsi %>% dplyr::select(individual_code, period, 
                             PSi, Weight_g, Gonadal_weight_g, IGS) %>% head(12)


mod_igs = glmmTMB(IGS ~ PSi + Sex +  (1 | period) + (1 | pool_month),
                  data = PSi_df_gsi,
                  family = Gamma(link = "log"))

summary(mod_igs)
