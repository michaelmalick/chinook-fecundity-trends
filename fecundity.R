## Chinook fecundity analysis

dir.create("./figures/fecundity/", showWarnings = FALSE)


## GAM -----------------------------------------------------
dat <- fecund
groups <- levels(dat$group)
mgcv_fits <- vector("list", length(groups))
names(mgcv_fits) <- groups
for(i in seq_along(groups)) {
    grp <- groups[i]
    fit <- mgcv::gam(fecundity_stnd ~ s(year, k = 10),
               data = fecund[group == groups[i], ],
               method = "REML")
    mgcv_fits[[i]] <- fit
}
save(mgcv_fits, file = "./outputs/mgcv_fits.RData")
summary(mgcv_fits[[1]])
summary(mgcv_fits[[2]])

g <- ggplot(fecund) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    aes(x = year, y = fecundity_stnd, color = region) +
    geom_point(na.rm = TRUE, size = 0.7) +
    geom_smooth(na.rm = TRUE, method = "gam", formula = y ~ s(x, k = 10)) +
    labs(y = "Standardized fecundity") +
    scale_color_manual(values = M1) +
    theme_simple() +
    theme(legend.position = "none") +
    facet_wrap( ~ group, ncol = 3)
print(g)
ggsave("./figures/fecundity/gam_fecundity_group.png", units = "in",
       width = 6, height = 6)



## Fit Kalman filter models using MARSS --------------------
fecund_kf <- fecund
## univariate stochastic level model: https://bit.ly/2Pc1zyA
mod  =  list(Z = matrix(1), A = matrix(0), R = matrix("r"),
             B = matrix(1), U = matrix(0), Q = matrix("q"),
             x0 = matrix("pi"), tinitx = 0)
kf <- function(x, model, method = "kem") {
    fit <- MARSS(as.vector(x), model, method = method,
                 control = list(maxit = 2000, conv.test.slope.tol = 1000))
    s <- as.vector(fit$states)
    se <- as.vector(fit$states.se)
    up <- s + 1.96 * se
    lo <- s - 1.96 * se
    m <- rep(fit$convergence, length(s))
    list(s, lo, up, m)
}
fecund_kf[ , c("kalman", "lower", "upper", "converge") := kf(fecundity_stnd, mod),
          by = "id_label"]
table(fecund_kf$converge)

save(fecund_kf, file = "./outputs/fecund_kf.RData")


## Kalman smooths by stock
g <- ggplot(fecund_kf) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    aes(x = year, y = kalman, color = run) +
    geom_ribbon(aes(ymin = lower, ymax = upper, x = year),
                color = NA, fill = "grey80") +
    geom_point() +
    geom_line() +
    labs(y = "Smoothed fecundity") +
    scale_color_manual(values = M1) +
    theme_simple() +
    facet_wrap( ~ id_label)
print(g)
ggsave("./figures/fecundity/kalman_smooth_stock.png", units = "in",
       width = 18, height = 11)


## Kalman smooths by group
g <- ggplot(fecund_kf) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    aes(x = year, y = kalman, group = id_label, color = region) +
    geom_line(na.rm = TRUE) +
    labs(y = "Smoothed fecundity") +
    scale_color_manual(values = M1) +
    theme_simple() +
    theme(legend.position = "none") +
    facet_wrap( ~ group, ncol = 2)
print(g)
ggsave("./figures/fecundity/kalman_smooth_group.png", units = "in",
       width = 6, height = 6)


## Kalman fitted by stock
g <- ggplot(fecund_kf) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    geom_point(aes(x = year, y = fecundity_stnd, color = run), na.rm = TRUE) +
    geom_line(aes(x = year, y = kalman), color = "grey30") +
    labs(y = "Standardized fecundity") +
    scale_color_manual(values = M1) +
    theme_simple() +
    facet_wrap( ~ id_label)
print(g)
ggsave("./figures/fecundity/kalman_fitted_stock.png", units = "in",
       width = 18, height = 11)


## SS RW models --------------------------------------------

## Setup data
fecund_rw <- fecund
fecund_rw_lst <- split(fecund_rw, by = "id_label")
save(fecund_rw, file = "./outputs/fecund_rw.RData")
save(fecund_rw_lst, file = "./outputs/fecund_rw_lst.RData")


## Run JAGS models
ssrw_fit <- vector("list", length(fecund_rw_lst))
jags_params <- c("sigma_pro", "sigma_obs", "X")
n_chains <- 4
n_burnin <- 10000
n_thin   <- 25
n_iter   <- 150000
((n_iter - n_burnin) / n_thin) * n_chains
for(i in seq_along(ssrw_fit)) {
    d <- fecund_rw_lst[[i]]
    y <- as.vector(d$fecundity_stnd)
    jags.data <- list(Y = y, N = length(y))
    set.seed(123)
    mod_ss <- jags(jags.data, parameters.to.save = jags_params,
                   model.file = "./stan/ss_rw.bug",
                   n.chains = n_chains,
                   n.burnin = n_burnin,
                   n.thin = n_thin,
                   n.iter = n_iter,
                   DIC = FALSE)
    ssrw_fit[[i]] <- mod_ss
}
save(ssrw_fit, file = "./outputs/ssrw_fit.RData")
# coda::varnames(as.mcmc(ssrw_fit[[1]]))

## MCMC diagnostics
# for(i in seq_along(ssrw_fit)) {
#     coda_diag(as.mcmc(ssrw_fit[[i]]))
# }

## Yrep
# for(i in seq_along(ssrw_fit)) {
#     yrep <- ss_yrep(as.mcmc(ssrw_fit[[i]]))
#     d <- fecund_rw_lst[[i]]
#     y <- as.vector(d$fecundity_stnd)
#     yrep_plot(y, yrep, main = paste(i, as.character(unique(d$hatchery))))
# }


## Check Neff and Rhat
ssrw_neff <- vector("list", length(ssrw_fit))
for(i in seq_along(ssrw_neff)) {
    s <- ssrw_fit[[i]]$BUGSoutput
    Neff <- s$summary[,"n.eff"]
    Rhat <- s$summary[,"Rhat"]
    ssrw_neff[[i]] <- data.table(model = i,
                                 par = names(Neff),
                                 neff = Neff,
                                 rhat = Rhat)
}
ssrw_neff <- rbindlist(ssrw_neff)
save(ssrw_neff, file = "./outputs/ssrw_neff.RData")
min(ssrw_neff$neff)
max(ssrw_neff$rhat)
ssrw_neff[rhat > 1.05]
ssrw_neff[neff < 500]


## Summarize model parameters
ssrw_summary <- vector("list", length(ssrw_fit))
for(i in seq_along(ssrw_fit)) {
    d <- fecund_rw_lst[[i]]
    y <- as.vector(d$fecundity_stnd)
    mc <- as.mcmc(ssrw_fit[[i]])
    df <- coda_df(mc, parameters = grep("X", coda::varnames(mc), value = TRUE))
    x <- rep(NA, length(y))
    up <- rep(NA, length(y))
    lo <- rep(NA, length(y))
    si <- rep(NA, length(y))
    for(j in 1:length(y)) {
        nam <- paste0("X[", j, "]")
        x[j] <- mean(df[[nam]])
        lo[j] <- quantile(df[[nam]], probs = 0.05)
        up[j] <- quantile(df[[nam]], probs = 0.90)
        si[j] <- ifelse(lo[j] > 0 | up[j] < 0, 1, 0)
    }
    ssrw_summary[[i]] <- data.table(state = x,
                                    lower_90 = lo,
                                    upper_90 = up,
                                    sig_90 = si)
}
ssrw_summary <- rbindlist(ssrw_summary)
fecund_rw_fit <- cbind(fecund_rw, ssrw_summary)
save(fecund_rw_fit, file = "./outputs/fecund_rw_fit.RData")


## Plot state time series
dat <- copy(fecund_rw_fit)
dat[ , fecund_na := ifelse(is.na(fecundity_stnd), 2, 16)]
datl <- unique(dat[ , "id_label"])
g <- ggplot(dat) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    geom_ribbon(aes(ymin = lower_90, ymax = upper_90, x = year, fill = region),
                color = NA, alpha = 0.2) +
    geom_line(aes(x = year, y = state, color = region)) +
    geom_point(data = dat[is.na(fecundity_stnd), ],
               aes(x = year, y = state), color = "grey25", size = 0.8) +
    geom_point(data = dat[sig_90 == 1, ],
               aes(x = year, y = state), color = "red3", size = 0.8) +
    geom_text(data = datl, aes(x = 1995, y = 2.05, label = id_label),
              color = "grey40", size = 2.7, hjust = 0) +
    scale_x_continuous(limits = c(1995, 2020),
                       breaks = c(1995, 2005, 2015)) +
    labs(y = "Smoothed fecundity", x = "Year", color = "", fill = "") +
    scale_color_manual(values = M1) +
    scale_fill_manual(values = M1) +
    theme_simple() +
    theme(legend.position = "top",
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.spacing.y = unit(-0.5, "pt"),
          panel.spacing.x = unit(-0.5, "pt")) +
    facet_wrap( ~ id_label, ncol = 5)
print(g)
ggsave("./figures/fecundity/ssrw_stock.png", units = "in",
       width = 18, height = 11)



## SS AR1 models -------------------------------------------

## Run JAGS models
ssar1_fit <- vector("list", length(fecund_rw_lst))
jags_params <- c("sigma_pro", "sigma_obs", "X", "phi")
n_chains <- 4
n_burnin <- 10000
n_thin   <- 25
n_iter   <- 150000
((n_iter - n_burnin) / n_thin) * n_chains
for(i in seq_along(ssar1_fit)) {
    d <- fecund_rw_lst[[i]]
    y <- as.vector(d$fecundity_stnd)
    jags.data <- list(Y = y, N = length(y))
    set.seed(123)
    mod_ss <- jags(jags.data, parameters.to.save = jags_params,
                   model.file = "./stan/ss_ar1.bug",
                   n.chains = n_chains,
                   n.burnin = n_burnin,
                   n.thin = n_thin,
                   n.iter = n_iter,
                   DIC = FALSE)
    ssar1_fit[[i]] <- mod_ss
}
save(ssar1_fit, file = "./outputs/ssar1_fit.RData")

## Check Neff and Rhat
ssar1_neff <- vector("list", length(ssar1_fit))
for(i in seq_along(ssar1_neff)) {
    s <- ssar1_fit[[i]]$BUGSoutput
    Neff <- s$summary[,"n.eff"]
    Rhat <- s$summary[,"Rhat"]
    ssar1_neff[[i]] <- data.table(model = i,
                                 par = names(Neff),
                                 neff = Neff,
                                 rhat = Rhat)
}
ssar1_neff <- rbindlist(ssar1_neff)
save(ssar1_neff, file = "./outputs/ssar1_neff.RData")
min(ssar1_neff$neff)
max(ssar1_neff$rhat)
ssar1_neff[rhat > 1.05]
ssar1_neff[neff < 500]


## Summarize phi
n <- length(ssar1_fit)
phi <- rep(NA, n)
for(i in 1:n) {
    mc  <- as.mcmc(ssar1_fit[[i]])
    phi_i <- coda_df(mc, parameters = grep ("phi", coda::varnames(mc), value = TRUE))
    phi[i] <- mean(phi_i$phi)
}
range(phi)
plot(phi, ylim = c(-1, 1)); abline(h = c(-1, 1), col = 2)


## Summarize model parameters
ssar1_summary <- vector("list", length(ssar1_fit))
for(i in seq_along(ssar1_fit)) {
    d <- fecund_rw_lst[[i]]
    y <- as.vector(d$fecundity_stnd)
    mc <- as.mcmc(ssar1_fit[[i]])
    df <- coda_df(mc, parameters = grep("X", coda::varnames(mc), value = TRUE))
    x <- rep(NA, length(y))
    up <- rep(NA, length(y))
    lo <- rep(NA, length(y))
    si <- rep(NA, length(y))
    for(j in 1:length(y)) {
        nam <- paste0("X[", j, "]")
        x[j] <- mean(df[[nam]])
        lo[j] <- quantile(df[[nam]], probs = 0.05)
        up[j] <- quantile(df[[nam]], probs = 0.90)
        si[j] <- ifelse(lo[j] > 0 | up[j] < 0, 1, 0)
    }
    ssar1_summary[[i]] <- data.table(state = x,
                                     lower_90 = lo,
                                     upper_90 = up,
                                     sig_90 = si)
}
ssar1_summary <- rbindlist(ssar1_summary)
fecund_ar1_fit <- cbind(fecund_rw, ssar1_summary)
save(fecund_ar1_fit, file = "./outputs/fecund_ar1_fit.RData")


## Plot state time series
dat <- copy(fecund_ar1_fit)
dat[ , fecund_na := ifelse(is.na(fecundity_stnd), 2, 16)]
datl <- unique(dat[ , "id_label"])
g <- ggplot(dat) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    geom_ribbon(aes(ymin = lower_90, ymax = upper_90, x = year, fill = region),
                color = NA, alpha = 0.2) +
    geom_line(aes(x = year, y = state, color = region)) +
    geom_point(data = dat[is.na(fecundity_stnd), ],
               aes(x = year, y = state), color = "grey25", size = 0.8) +
    geom_point(data = dat[sig_90 == 1, ],
               aes(x = year, y = state), color = "red3", size = 0.8) +
    geom_text(data = datl, aes(x = 1995, y = 2.05, label = id_label),
              color = "grey40", size = 2.7, hjust = 0) +
    scale_x_continuous(limits = c(1995, 2020),
                       breaks = c(1995, 2005, 2015)) +
    labs(y = "Smoothed fecundity", x = "Year", color = "", fill = "") +
    scale_color_manual(values = M1) +
    scale_fill_manual(values = M1) +
    theme_simple() +
    theme(legend.position = "top",
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.spacing.y = unit(-0.5, "pt"),
          panel.spacing.x = unit(-0.5, "pt")) +
    facet_wrap( ~ id_label, ncol = 5)
print(g)
ggsave("./figures/fecundity/ssar1_stock.png", units = "in",
       width = 18, height = 11)



## Compare RW + AR1
ss_comp <- data.table(id_label = fecund_rw_fit$id_label,
                      state_rw = fecund_rw_fit$state,
                      state_ar1 = fecund_ar1_fit$state)
ss_cor <- ss_comp[ , .(cor = cor(state_rw, state_ar1)), by = id_label]
mean(ss_cor$cor)


## Correlation matrix --------------------------------------
fecund_cor <- copy(fecund)


cor_all <- dd_cor_group(fecund_cor,
                        id.var = "id_label",
                        grp.var = "group",
                        cor.var = "fecundity", time.var = "year",
                        min.overlap = 5,
                        plot = TRUE)
cor_early <- dd_cor_group(fecund_cor[year <= 2007],
                          id.var = "id_label",
                          grp.var = "group",
                          cor.var = "fecundity", time.var = "year",
                          min.overlap = 5,
                          plot = TRUE)
cor_late <- dd_cor_group(fecund_cor[year > 2007],
                         id.var = "id_label",
                         grp.var = "group",
                         cor.var = "fecundity", time.var = "year",
                         min.overlap = 5,
                         plot = TRUE)

save(cor_all, file = "./outputs/cor_all.RData")
save(cor_early, file = "./outputs/cor_early.RData")
save(cor_late, file = "./outputs/cor_late.RData")

## Global mean correlations
mean(cor_all$cor, na.rm = TRUE)
mean(cor_early$cor, na.rm = TRUE)
mean(cor_late$cor, na.rm = TRUE)

## Mean correlation within groups
cor_all[group1 == group2, .(mean = mean(cor, na.rm = TRUE))]
cor_early[group1 == group2, .(mean = mean(cor, na.rm = TRUE))]
cor_late[group1 == group2, .(mean = mean(cor, na.rm = TRUE))]

## Mean correlation among groups
cor_all[group1 != group2, .(mean = mean(cor, na.rm = TRUE))]
cor_early[group1 != group2, .(mean = mean(cor, na.rm = TRUE))]
cor_late[group1 != group2, .(mean = mean(cor, na.rm = TRUE))]

# Group mean correlations
cor_all[group1 == group2, .(mean = mean(cor, na.rm = TRUE)), by = group1]
cor_early[group1 == group2, .(mean = mean(cor, na.rm = TRUE)), by = group1]
cor_late[group1 == group2, .(mean = mean(cor, na.rm = TRUE)), by = group1]

g0 <- cor_plot_region(cor_all, title = "All years")
ggsave("./figures/fecundity/cor_all.png", units = "in", width = 9, height = 7)

g1 <- cor_plot_region(cor_early, title = "1995-2007")
ggsave("./figures/fecundity/cor_early.png", units = "in", width = 9, height = 7)

g2 <- cor_plot_region(cor_late, title = "2008-2019")
ggsave("./figures/fecundity/cor_late.png", units = "in", width = 9, height = 7)



## DFA: fit ------------------------------------------------
## see: https://cran.r-project.org/web/packages/bayesdfa/vignettes/bayesdfa.html
fecund_dfa <- copy(fecund)
size_id <- as.character(unique(fecund_size$id_label))
fecund_dfa$size <- ifelse(fecund_dfa$id_label %in% size_id, "yes", "no")

mat_dfa <- dcast(fecund_dfa, year ~ id_label, value.var = "fecundity")
mat_dfa[ , year := NULL]
y_dfa <- t(as.matrix(mat_dfa))

save(fecund_dfa, file = "./outputs/fecund_dfa.RData")
save(y_dfa, file = "./outputs/y_dfa.RData")

matplot(t(y_dfa),
        type = "l",
        ylab = "Response", xlab = "Time")


## 1 Trend
set.seed(123)
dfa1 <- fit_dfa(y = y_dfa, num_trends = 1, scale = "zscore",
                iter = 4000, chains = 4, cores = 4, thin = 1,
                par_list = NULL, seed = 29292929)
save(dfa1, file = "./outputs/dfa1.RData")
is_converged(dfa1, 1.04, parameters = c("sigma", "x\\[", "Z\\["))
dfa1rot <- rotate_trends(dfa1, conf_level = 0.90)
plot_trends(dfa1rot)
plot_fitted(dfa1)
plot_loadings(dfa1rot)
dfa1_loo <- loo(dfa1)
save(dfa1_loo, file = "./outputs/dfa1_loo.RData")

## SENSITIVITY: Use 'varIndx' to set unequal observation errors
set.seed(123)
dfa1_unequal = fit_dfa(y = y_dfa, num_trends = 1, scale = "zscore",
                       varIndx = 1:nrow(y_dfa),
                       est_correlation = FALSE,
                       init = 0,
                       iter = 4000, chains = 4, cores = 4, thin = 1,
                       par_list = NULL)
save(dfa1_unequal, file = "./outputs/dfa1_unequal.RData")
is_converged(dfa1_unequal, 1.04, parameters = c("sigma", "x\\[", "Z\\["))
dfa1rot_unequal <- rotate_trends(dfa1_unequal, conf_level = 0.95)
plot_trends(dfa1rot_unequal)

## observation error SDs
sigmas <- rstan::extract(dfa1_unequal$model)$sigma
ncol(sigmas) == nrow(y_dfa)
head(sigmas)


## 2 Trend
set.seed(123)
dfa2 <- fit_dfa(y = y_dfa, num_trends = 2, scale = "zscore",
                iter = 6000, chains = 4, cores = 4, thin = 1,
                par_list = NULL, seed = 1234)
save(dfa2, file = "./outputs/dfa2.RData")
is_converged(dfa2, 1.04, parameters = c("sigma", "x\\[", "Z\\["))
dfa2rot <- rotate_trends(dfa2, conf_level = 0.90)
plot_trends(dfa2rot)
plot_fitted(dfa2)
plot_loadings(dfa2rot)
dfa2_loo <- loo(dfa2)
save(dfa2_loo, file = "./outputs/dfa2_loo.RData")


## 3 Trend
set.seed(123)
dfa3 <- fit_dfa(y = y_dfa, num_trends = 3, scale = "zscore",
                iter = 6000, chains = 4, cores = 4, thin = 1,
                par_list = NULL, seed = 1234)
save(dfa3, file = "./outputs/dfa3.RData")
is_converged(dfa3, 1.04, parameters = c("sigma", "x\\[", "Z\\["))
dfa3rot <- rotate_trends(dfa3, conf_level = 0.90)
plot_trends(dfa3rot)
plot_fitted(dfa3)
plot_loadings(dfa3rot)
dfa3_loo <- loo(dfa3)
save(dfa3_loo, file = "./outputs/dfa3_loo.RData")


## 4 Trend
set.seed(123)
dfa4 <- fit_dfa(y = y_dfa, num_trends = 4, scale = "zscore",
                iter = 8000, chains = 4, cores = 4, thin = 1,
                par_list = NULL, seed = 1234)
save(dfa4, file = "./outputs/dfa4.RData")
is_converged(dfa4, 1.04, parameters = c("sigma", "x\\[", "Z\\["))
dfa4rot <- rotate_trends(dfa4, conf_level = 0.90)
plot_trends(dfa4rot)
plot_fitted(dfa4)
plot_loadings(dfa4rot, names = rownames(y_dfa))
dfa4_loo <- loo(dfa4)
save(dfa4_loo, file = "./outputs/dfa4_loo.RData")


check_divergences(dfa1$model)
check_divergences(dfa2$model)
check_divergences(dfa3$model)
check_divergences(dfa4$model)

check_treedepth(dfa1$model)
check_treedepth(dfa2$model)
check_treedepth(dfa3$model)
check_treedepth(dfa4$model)



## DFA: model comparison -----------------------------------
loo::loo_compare(dfa1_loo, dfa2_loo, dfa3_loo, dfa4_loo)

looic <- c(dfa1_loo$estimates["looic", "Estimate"],
           dfa2_loo$estimates["looic", "Estimate"],
           dfa3_loo$estimates["looic", "Estimate"],
           dfa4_loo$estimates["looic", "Estimate"])
looic_se <- c(dfa1_loo$estimates["looic", "SE"],
              dfa2_loo$estimates["looic", "SE"],
              dfa3_loo$estimates["looic", "SE"],
              dfa4_loo$estimates["looic", "SE"])
dfa_looic <- data.frame(model = c("dfa1", "dfa2", "dfa3", "dfa4"),
                        looic = looic, se = looic_se, trends = 1:4)
dfa_looic <- dfa_looic[order(dfa_looic$looic), ]
dfa_looic$lower <- dfa_looic$looic - (dfa_looic$se * 2)
dfa_looic$upper <- dfa_looic$looic + (dfa_looic$se * 2)
save(dfa_looic, file = "./outputs/dfa_looic.RData")
print(dfa_looic)

g <- ggplot(dfa_looic) +
    geom_point(aes(x = as.factor(trends), y = looic), size = 2) +
    geom_segment(aes(x = as.factor(trends), xend = as.factor(trends),
                     y = lower, yend = upper)) +
    labs(x = "Number of trends", y = "LOOIC") +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/fecundity/dfa_looic.png", units = "in", width = 4, height = 3)



## DFA: best model -----------------------------------------
dfa1rot <- rotate_trends(dfa1, conf_level = 0.95)
dfa2rot <- rotate_trends(dfa2, conf_level = 0.95)
save(dfa1rot, file = "./outputs/dfa1rot.RData")
save(dfa2rot, file = "./outputs/dfa2rot.RData")

## 1. MCMC diagnostics
check_hmc_diagnostics(dfa1$model)
rhat_highest_dfa(dfa1)
neff_lowest_dfa(dfa1)
pdf("./figures/fecundity/dfa1_diag.pdf", width = 6, height  = 4)
    post <- dfa1$samples
    p <- mcmc_trace(post,  regex_pars = "x\\[")
    print(p)
    p <- mcmc_trace(post,  regex_pars = "Z\\[")
    print(p)
    p <- mcmc_trace(post,  regex_pars = "sigma")
    print(p)
dev.off()


## 2. DFA trend
dfa1_trends <- dfa_trends(dfa1rot, years = 1995:2019)
g <- ggplot(dfa1_trends) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    aes(x = time, y = estimate) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4) +
    geom_line() +
    labs(x = "Year", y = "") +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/fecundity/dfa1_trends.png", units = "in",
       width = 4, height = 3)

dfa2_trends <- dfa_trends(dfa2rot, years = 1995:2019)
g <- ggplot(dfa2_trends) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    aes(x = time, y = estimate) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4) +
    geom_line() +
    labs(x = "Year", y = "") +
    facet_wrap( ~ trend_number) +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/fecundity/dfa2_trends.png", units = "in",
       width = 4, height = 3)


## 2. Fitted trends
dfa1_fitted <- dfa_fitted(dfa1, names = rownames(y_dfa))
dfa1_fitted <- as.data.table(dfa1_fitted)
dfa1_fitted[ , year := 1995:2019, by = ID]
grps <- fecund_dfa[ , .(run, region, group, id_label, fecundity_stnd, year)]
names(grps) <- c("run", "region", "group", "ID", "fecundity_stnd", "year")
dfa1_fitted <- merge(grps, dfa1_fitted, by = c("ID", "year"),
                     sort = FALSE)
dfa1_fitted[ , group := factor(group, levels = levels(fecund_dfa$group))]
dfa1_fitted[ , ID := factor(ID, levels = levels(fecund_dfa$id_label))]
dfa1_fitted[ , region := factor(region, levels = levels(fecund_dfa$region))]
dfa1_fitted <- dfa1_fitted[order(year), .SD, by = ID]

g <- ggplot(dfa1_fitted) +
    geom_point(aes(x = year, y = fecundity_stnd, color = run),
               na.rm = TRUE) +
    geom_line(aes(x = year, y = estimate), color = "grey30", na.rm = TRUE) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    scale_color_manual(values = M1) +
    facet_wrap( ~ ID) +
    theme_simple(grid = TRUE) +
    theme(legend.position = "none")
print(g)
ggsave("./figures/fecundity/dfa1_fitted_stock.png", units = "in",
       width = 18, height = 11)

g <- ggplot(dfa1_fitted) +
    geom_line(aes(x = year, y = estimate, group = ID, color = region)) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    scale_color_manual(values = M1) +
    facet_wrap( ~ group, ncol = 2) +
    theme_simple(grid = TRUE) +
    theme(legend.position = "none")
print(g)
ggsave("./figures/fecundity/dfa1_fitted_group.png", units = "in",
       width = 6, height = 6)


## 3. DFA loadings
dfa1_loadings <- dfa_loadings(dfa1rot, names = rownames(y_dfa), conf_level = 0.90)
dfa1_loadings <- as.data.table(dfa1_loadings)
grps <- unique(fecund_dfa[ , .(run, region, group, id_label)])
names(grps) <- c("run", "region", "group", "name")
dfa1_loadings <- merge(dfa1_loadings, grps, by = "name")
dfa1_loadings[ , group := factor(group, levels = levels(fecund_dfa$group))]
dfa1_loadings[ , run := factor(run, levels = levels(fecund_dfa$run))]

setorder(dfa1_loadings, "group")
dfa1_loadings[ , name := factor(name, levels = unique(name))]
dfa1_loadings[ , color := .GRP, by = group]
dfa1_loadings[ , diff := c(0, diff(color))]
dfa1_loadings[ , gr_line := ifelse(diff == 1, .I - 0.5, NA) , by = .I]

g <- ggplot(dfa1_loadings) +
    geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
    geom_point(aes(x = median, y = name, color = group)) +
    geom_segment(aes(x = lower, xend = upper,
                     y = name, yend = name, color = group)) +
    geom_hline(yintercept = dfa1_loadings$gr_line,
               color = "grey70", na.rm = TRUE, size = 0.3) +
    labs(x = "Loading", y = "", color = "") +
    theme_simple()
print(g)
ggsave("./figures/fecundity/dfa1_loadings.png", units = "in",
       width = 8, height = 8)


## DFA: size data stocks only ------------------------------
mat_dfa_s <- dcast(fecund_dfa[size == "yes", ], year ~ id_label, value.var = "fecundity")
mat_dfa_s[ , year := NULL]
y_dfa_s <- t(as.matrix(mat_dfa_s))

matplot(t(y_dfa_s),
        type = "l",
        ylab = "Response", xlab = "Time")


## 1 Trend: size only
set.seed(123)
dfa1_s <- fit_dfa(y = y_dfa_s, num_trends = 1, scale = "zscore",
                  iter = 4000, chains = 4, cores = 4, thin = 1,
                  par_list = NULL, seed = 29292929)
save(dfa1_s, file = "./outputs/dfa1_s.RData")
is_converged(dfa1_s, 1.04, parameters = c("sigma", "x\\[", "Z\\["))
dfa1rot_s <- rotate_trends(dfa1_s, conf_level = 0.95)
save(dfa1rot_s, file = "./outputs/dfa1rot_s.RData")
plot_trends(dfa1rot_s)
plot_fitted(dfa1_s)
plot_loadings(dfa1rot_s)



## DFA: alt var-covar matrices -----------------------------
## base model: diagonal and equal
## model c1: diagonal and unequal
## setting est_correlation = TRUE gives an unconstrained var-covar matrix

dfa1_c1 <- fit_dfa(y = y_dfa, num_trends = 1, scale = "zscore",
                   varIndx = 1:nrow(y_dfa),
                   est_correlation = FALSE,
                   init = 0,
                   iter = 4000, chains = 4, cores = 4, thin = 1,
                   par_list = NULL, seed = 29292929)

save(dfa1_c1, file = "./outputs/dfa1_c1.RData")
is_converged(dfa1_c1, 1.01, parameters = c("sigma", "x\\[", "Z\\["))
dfa1rot_c1 <- rotate_trends(dfa1_c1, conf_level = 0.95, invert = TRUE)
plot_trends(dfa1rot_c1)
plot_fitted(dfa1_c1)
plot_loadings(dfa1rot_c1)
dfa1_c1_loo <- loo(dfa1_c1)
save(dfa1_c1_loo, file = "./outputs/dfa1_c1_loo.RData")
loo_compare(dfa1_loo, dfa1_c1_loo)


## model c2: equal variance, unconstrained covar
## model is too big to estimate -- didn't get past 1% sampling in 12 hours
# dfa1_c2 <- fit_dfa(y = y_dfa, num_trends = 1, scale = "zscore",
#                    varIndx = rep(1, nrow(y_dfa)),
#                    est_correlation = TRUE,
#                    init = 0,
#                    iter = 4000, chains = 4, cores = 4, thin = 1,
#                    par_list = NULL, seed = 29292929)

# save(dfa1_c2, file = "./outputs/dfa1_c2.RData")
# is_converged(dfa1_c2, 1.01, parameters = c("sigma", "x\\[", "Z\\["))
# dfa1rot_c2 <- rotate_trends(dfa1_c2, conf_level = 0.95)
# plot_trends(dfa1rot_c2)
# plot_fitted(dfa1_c2)
# plot_loadings(dfa1rot_c2)
# dfa1_c2_loo <- loo(dfa1_c2)
# save(dfa1_c2_loo, file = "./outputs/dfa1_c2_loo.RData")



## Extract sigmas
s0 <- rstan::extract(dfa1$model)$sigma
head(s0)
s1 <- rstan::extract(dfa1_c1$model)$sigma
head(s1)


