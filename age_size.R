## Analysis of age/size data

dir.create("./figures/age_size/", showWarnings = FALSE)


## Size: year effects --------------------------------------
fit_age3 <- brm(length ~ as.factor(year),
                data = age_size[age == 3, ],
                cores = 4, chains = 4,
                seed = 123, iter = 4000)

fit_age4 <- brm(length ~ as.factor(year),
                data = age_size[age == 4, ],
                cores = 4, chains = 4,
                seed = 123, iter = 4000)

fit_age5 <- brm(length ~ as.factor(year),
                data = age_size[age == 5, ],
                cores = 4, chains = 4,
                seed = 123, iter = 4000)

save(fit_age3, file = "./outputs/fit_age3.RData")
save(fit_age4, file = "./outputs/fit_age4.RData")
save(fit_age5, file = "./outputs/fit_age5.RData")


pp_check(fit_age3, type = "dens_overlay", nsamples = 50)
pp_check(fit_age4, type = "dens_overlay", nsamples = 50)
pp_check(fit_age5, type = "dens_overlay", nsamples = 50)


fits <- list(fit_age3, fit_age4, fit_age5)
lst <- vector("list", 3)
for(i in 1:3) {
    int <- as.data.frame(fits[[i]], variable = "b_Intercept")
    yrs <- as.data.frame(fits[[i]], variable = "b_.*year", regex = TRUE)
    eff <- as.data.frame(lapply(yrs, function(x) x + int[[1]]))
    df <- data.table(age = i + 2,
                     year = 1995:2019,
                     median = c(median(int[[1]]), apply(eff, 2, median)),
                     lower = c(quantile(int[[1]], probs = 0.025),
                               apply(eff, 2, function(x) quantile(x, probs = 0.025))),
                     upper = c(quantile(int[[1]], probs = 0.975),
                               apply(eff, 2, function(x) quantile(x, probs = 0.975))))
    lst[[i]] <- df
}
size_trend <- rbindlist(lst)
save(size_trend, file = "./outputs/size_trend.RData")

g <- ggplot(size_trend) +
    aes(x = year, y = median) +
    geom_point() +
    geom_segment(aes(x = year, xend = year, y = upper, yend = lower)) +
    geom_smooth(method = "loess", formula = y ~ x, span = 0.5, se = FALSE,
                color = "grey60", size = 0.5) +
    labs(x = "Year", y = "Length (mm)") +
    theme_simple() +
    facet_wrap( ~ age, ncol = 1, scale = "free_y")
print(g)
ggsave("./figures/age_size/age_size_year_effects.png",
       width = 5, height = 8)


## Size DFA ------------------------------------------------
set.seed(12345)

size_dfa <- copy(age_size)
save(size_dfa, file = "./outputs/size_dfa.RData")

size_dfa3 <- size_dfa[age == 3, ]
size_dfa4 <- size_dfa[age == 4, ]
size_dfa5 <- size_dfa[age == 5, ]

mat_dfa3 <- dcast(size_dfa3, year ~ id_label, value.var = "length")
mat_dfa4 <- dcast(size_dfa4, year ~ id_label, value.var = "length")
mat_dfa5 <- dcast(size_dfa5, year ~ id_label, value.var = "length")

mat_dfa3[ , year := NULL]
mat_dfa4[ , year := NULL]
mat_dfa5[ , year := NULL]

y_dfa_age3 <- t(as.matrix(mat_dfa3))
y_dfa_age4 <- t(as.matrix(mat_dfa4))
y_dfa_age5 <- t(as.matrix(mat_dfa5))
save(y_dfa_age3, file = "./outputs/y_dfa_age3.RData")
save(y_dfa_age4, file = "./outputs/y_dfa_age4.RData")
save(y_dfa_age5, file = "./outputs/y_dfa_age5.RData")

matplot(t(y_dfa_age3), type = "l", ylab = "Response", xlab = "Time")
matplot(t(y_dfa_age4), type = "l", ylab = "Response", xlab = "Time")
matplot(t(y_dfa_age5), type = "l", ylab = "Response", xlab = "Time")

## Age-3
set.seed(123)
dfa1_age3 <- fit_dfa(y = y_dfa_age3, num_trends = 1, scale = "zscore",
                     iter = 4000, chains = 4, cores = 4, thin = 1,
                     seed = 98765432)
save(dfa1_age3, file = "./outputs/dfa1_age3.RData")
is_converged(dfa1_age3, 1.01, parameters = c("sigma", "x\\[", "Z\\["))
dfa1_age3rot <- rotate_trends(dfa1_age3, conf_level = 0.95)
save(dfa1_age3rot, file = "./outputs/dfa1_age3rot.RData")
plot_trends(dfa1_age3rot)
plot_fitted(dfa1_age3)
plot_loadings(dfa1_age3rot)

## Age-4
set.seed(123)
dfa1_age4 <- fit_dfa(y = y_dfa_age4, num_trends = 1, scale = "zscore",
                     iter = 4000, chains = 4, cores = 4, thin = 1,
                     seed = 918273)
save(dfa1_age4, file = "./outputs/dfa1_age4.RData")
is_converged(dfa1_age4, 1.01, parameters = c("sigma", "x\\[", "Z\\["))
dfa1_age4rot <- rotate_trends(dfa1_age4, conf_level = 0.95)
save(dfa1_age4rot, file = "./outputs/dfa1_age4rot.RData")
plot_trends(dfa1_age4rot)
plot_fitted(dfa1_age4)
plot_loadings(dfa1_age4rot)

## Age-5
set.seed(123)
dfa1_age5 <- fit_dfa(y = y_dfa_age5, num_trends = 1, scale = "zscore",
                     iter = 4000, chains = 4, cores = 4, thin = 1,
                     seed = 191919)
save(dfa1_age5, file = "./outputs/dfa1_age5.RData")
is_converged(dfa1_age5, 1.01, parameters = c("sigma", "x\\[", "Z\\["))
dfa1_age5rot <- rotate_trends(dfa1_age5, conf_level = 0.95)
save(dfa1_age5rot, file = "./outputs/dfa1_age5rot.RData")
plot_trends(dfa1_age5rot)
plot_fitted(dfa1_age5)
plot_loadings(dfa1_age5rot)

check_hmc_diagnostics(dfa1_age3$model)
check_hmc_diagnostics(dfa1_age4$model)
check_hmc_diagnostics(dfa1_age5$model)

rhat_highest_dfa(dfa1_age3)
rhat_highest_dfa(dfa1_age4)
rhat_highest_dfa(dfa1_age5)

neff_lowest_dfa(dfa1_age3)
neff_lowest_dfa(dfa1_age4)
neff_lowest_dfa(dfa1_age5)


pdf("./figures/age_size/dfa1_age3_diag.pdf", width = 6, height  = 4)
    post <- dfa1_age3$samples
    p <- mcmc_trace(post,  regex_pars = "x\\[")
    print(p)
    p <- mcmc_trace(post,  regex_pars = "Z\\[")
    print(p)
    p <- mcmc_trace(post,  regex_pars = "sigma")
    print(p)
dev.off()
pdf("./figures/age_size/dfa1_age4_diag.pdf", width = 6, height  = 4)
    post <- dfa1_age4$samples
    p <- mcmc_trace(post,  regex_pars = "x\\[")
    print(p)
    p <- mcmc_trace(post,  regex_pars = "Z\\[")
    print(p)
    p <- mcmc_trace(post,  regex_pars = "sigma")
    print(p)
dev.off()
pdf("./figures/age_size/dfa1_age5_diag.pdf", width = 6, height  = 4)
    post <- dfa1_age5$samples
    p <- mcmc_trace(post,  regex_pars = "x\\[")
    print(p)
    p <- mcmc_trace(post,  regex_pars = "Z\\[")
    print(p)
    p <- mcmc_trace(post,  regex_pars = "sigma")
    print(p)
dev.off()



##  Age-3 trend
dfa1_age3trends <- dfa_trends(dfa1_age3rot, years = 1995:2019)
g <- ggplot(dfa1_age3trends) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    aes(x = time, y = estimate) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4) +
    geom_line() +
    labs(x = "Year", y = "") +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/age_size/dfa1_age3_trends.png", units = "in",
       width = 4, height = 3)

##  Age-4 trend
dfa1_age4trends <- dfa_trends(dfa1_age4rot, years = 1995:2019)
g <- ggplot(dfa1_age4trends) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    aes(x = time, y = estimate) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4) +
    geom_line() +
    labs(x = "Year", y = "") +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/age_size/dfa1_age4_trends.png", units = "in",
       width = 4, height = 3)

##  Age-5 trend
dfa1_age5trends <- dfa_trends(dfa1_age5rot, years = 1995:2019)
g <- ggplot(dfa1_age5trends) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    aes(x = time, y = estimate) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4) +
    geom_line() +
    labs(x = "Year", y = "") +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/age_size/dfa1_age5_trends.png", units = "in",
       width = 4, height = 3)



## Mean age ------------------------------------------------
mage <- age_size[ , .(mean = weighted.mean(age, w = nage)), by = .(year, id_label)]
mage[ , year_fac := as.factor(year)]

g <- ggplot(mage) +
    # geom_line(aes(x = year, y = mean, group = id_label)) +
    geom_point(aes(x = year, y = mean, group = id_label)) +
    geom_smooth(aes(x = year, y = mean), method = "gam", formula = y ~ s(x)) +
    theme_simple()
print(g)

g <- ggplot(mage) +
    aes(x = year, y = mean, group = id_label) +
    geom_point() +
    # geom_line() +
    geom_smooth(method = "gam", formula = y ~ s(x)) +
    facet_wrap( ~ id_label) +
    theme_simple()
print(g)

fit_mage <- brm(mean ~ year_fac,
                data = mage,
                cores = 4, chains = 4,
                seed = 123, iter = 4000)
save(fit_mage, file = "./outputs/fit_mage.RData")

ce <-  conditional_effects(fit_mage)[[1]]
ce$year <- as.numeric(as.character(ce$year_fac))
g <- ggplot(ce) +
    aes(x = year, y = estimate__) +
    geom_point() +
    # geom_point(data = mage, aes(x = year, y = mean), color = "grey70", size = 0.8) +
    geom_segment(aes(x = year, xend = year, y = upper__, yend = lower__)) +
    geom_smooth(method = "loess", formula = y ~ x, span = 0.5, se = FALSE,
                color = "grey60", size = 0.5) +
    labs(x = "Year", y = "Mean age") +
    scale_x_continuous(limits = c(1995, 2019),
                       breaks = c(1995, 2000, 2005, 2010, 2015)) +
    theme_simple()
print(g)
ggsave("./figures/age_size/mean_size_yr.png", g,
       units = "in", width = 6, height = 4)



## Single: fecundity ~ size --------------------------------
g <- ggplot(fecund_size) +
    aes(x = length_avg_anom2, y = fecundity_anom2) +
    geom_point(na.rm = TRUE) +
    geom_smooth(method = "lm", formula = y ~ x, na.rm = TRUE) +
    facet_wrap( ~ id_label) +
    theme_simple()
print(g)


ids <- unique(fecund_size$id_number)
lm_fits <- vector("list", length(ids))

priors_lm <- c(set_prior("student_t(3, 0, 50)", class = "b"),
               set_prior("student_t(3, 0, 700)", class = "Intercept"),
               set_prior("student_t(3, 0, 500)", class = "sigma"))
form_lm <- fecundity_anom2 ~ length_avg_anom2

for(i in seq_along(lm_fits)) {
    dat_i <- fecund_size[id_number == ids[i], ]
    fit_i <- brm(form_lm,
                 data = dat_i,
                 cores = 4, chains = 4,
                 prior = priors_lm,
                 save_pars = save_pars(all = TRUE),
                 seed = 123, iter = 4000)
    lm_fits[[i]] <- fit_i
}


save(lm_fits, file = "./outputs/lm_fits.RData")



## HBM: fecundity ~ size -----------------------------------

## FS0: independent w/ shared variance
priors_fs0 <- c(set_prior("student_t(3, 0, 50)", class = "b"),
                set_prior("student_t(3, 0, 700)", class = "Intercept"),
                set_prior("student_t(3, 0, 500)", class = "sigma"))
form_fs0 <- fecundity_anom2 ~ id_label + id_label:length_avg_anom2

fs0_prior <- brm(form_fs0,
                 data = fecund_size,
                 prior = priors_fs0,
                 sample_prior = "only",
                 seed = 123, iter = 4000)
pp_check(fs0_prior, type = "dens_overlay", nsamples = 100)
pp_check(fs0_prior, type = "scatter_avg", nsamples = 10)

fit_fs0 <- brm(form_fs0,
               data = fecund_size,
               cores = 4, chains = 4,
               prior = priors_fs0,
               save_pars = save_pars(all = TRUE),
               seed = 123, iter = 4000)
fit_fs0 <- add_criterion(fit_fs0, c("loo", "bayes_R2"), moment_match = TRUE)
save(fit_fs0, file = "./outputs/fit_fs0.RData")
pp_check(fit_fs0, type = "dens_overlay", nsamples = 50)
pp_check(fit_fs0, type = "scatter_avg", nsamples = 10)


## FS1: common across all stocks
priors_fs1 <- c(set_prior("student_t(3, 0, 50)", class = "b"),
                set_prior("student_t(3, 0, 700)", class = "Intercept"),
                set_prior("student_t(3, 0, 500)", class = "sigma"))
form_fs1 <- fecundity_anom2 ~ length_avg_anom2

fs1_prior <- brm(form_fs1,
                 data = fecund_size,
                 prior = priors_fs1,
                 sample_prior = "only",
                 seed = 123, iter = 4000)
pp_check(fs1_prior, type = "dens_overlay", nsamples = 100)
pp_check(fs1_prior, type = "scatter_avg", nsamples = 10)

fit_fs1 <- brm(form_fs1,
               data = fecund_size,
               cores = 4, chains = 4,
               prior = priors_fs1,
               save_pars = save_pars(all = TRUE),
               seed = 123, iter = 4000)
fit_fs1 <- add_criterion(fit_fs1, c("loo", "bayes_R2"))
save(fit_fs1, file = "./outputs/fit_fs1.RData")
pp_check(fit_fs1, type = "dens_overlay", nsamples = 50)
pp_check(fit_fs1, type = "scatter_avg", nsamples = 10)



## FS2: varying slope + intercept (no cor in random effects)
priors_fs2 <- c(set_prior("student_t(3, 0, 50)", class = "b"),
                set_prior("student_t(3, 0, 700)", class = "Intercept"),
                set_prior("student_t(3, 0, 50)", class = "sd"),
                set_prior("student_t(3, 0, 500)", class = "sigma"))
form_fs2 <- fecundity_anom2 ~ length_avg_anom2 + (length_avg_anom2 || id_label)

fs2_prior <- brm(form_fs2,
                 data = fecund_size,
                 prior = priors_fs2,
                 sample_prior = "only",
                 seed = 123, iter = 4000)
save(fs2_prior, file = "./outputs/fs2_prior.RData")
pp_check(fs2_prior, type = "dens_overlay", nsamples = 100)
pp_check(fs2_prior, type = "scatter_avg", nsamples = 10)
pp_check(fs2_prior, type = "dens_overlay", nsamples = 100)

pp = predict(fs2_prior, summary = FALSE)
y = fs2_prior$data$fecundity_anom2
plot(pp[1, ], y)
plot(density(pp[1, ]))
lines(density(y), col = 2)

plot(density(y), xlim = c(-1e4, 1e4))
for(i in 1:50) lines(density(pp[i, ]), col = "grey70")
lines(density(y), xlim = c(-1e4, 1e4))

fit_fs2 <- brm(form_fs2,
               data = fecund_size,
               cores = 4, chains = 4,
               prior = priors_fs2,
               control = list(adapt_delta = 0.9),
               save_pars = save_pars(all = TRUE),
               seed = 123, iter = 4000)
fit_fs2 <- add_criterion(fit_fs2, c("loo", "bayes_R2"), moment_match = TRUE)
save(fit_fs2, file = "./outputs/fit_fs2.RData")
pp_check(fit_fs2, type = "dens_overlay", nsamples = 50)
pp_check(fit_fs2, type = "scatter_avg", nsamples = 10)
plot(conditional_effects(fit_fs2), points = TRUE)


# Exclude Kalama Falls Spring 1996 outlier fecundity (low value)
fit_fs2_noout <- brm(form_fs2,
                     data = fecund_size[fecundity > 2500, ],
                     cores = 4, chains = 4,
                     prior = priors_fs2,
                     save_pars = save_pars(all = TRUE),
                     seed = 123, iter = 4000)
fit_fs2_noout <- add_criterion(fit_fs2_noout, c("loo", "bayes_R2"), moment_match = TRUE)
save(fit_fs2_noout, file = "./outputs/fit_fs2_noout.RData")
pp_check(fit_fs2_noout, type = "dens_overlay", nsamples = 50)
pp_check(fit_fs2_noout, type = "scatter_avg", nsamples = 10)
plot(conditional_effects(fit_fs2_noout), points = TRUE)


## FS3: varying slope + intercept (cor in random effects)
priors_fs3 <- c(set_prior("student_t(3, 0, 50)", class = "b"),
                set_prior("student_t(3, 0, 700)", class = "Intercept"),
                set_prior("student_t(3, 0, 500)", class = "sigma"),
                set_prior("lkj(1)", class = "cor"))
form_fs3 <- fecundity_anom2 ~ length_avg_anom2 + (length_avg_anom2 | id_label)
get_prior(form_fs3, data = fecund_size)

fs3_prior <- brm(form_fs3,
                 data = fecund_size,
                 prior = priors_fs3,
                 sample_prior = "only",
                 seed = 123, iter = 4000)
pp_check(fs3_prior, type = "dens_overlay", nsamples = 100)
pp_check(fs3_prior, type = "scatter_avg", nsamples = 10)

fit_fs3 <- brm(form_fs3,
               data = fecund_size,
               cores = 4, chains = 4,
               prior = priors_fs3,
               save_pars = save_pars(all = TRUE),
               seed = 123, iter = 4000)
fit_fs3 <- add_criterion(fit_fs3, c("loo", "bayes_R2"), moment_match = TRUE)
save(fit_fs3, file = "./outputs/fit_fs3.RData")
pp_check(fit_fs3, type = "dens_overlay", nsamples = 50)
pp_check(fit_fs3, type = "scatter_avg", nsamples = 10)


loo_compare(fit_fs0, fit_fs1, fit_fs2, fit_fs3)
bayes_R2(fit_fs0)
bayes_R2(fit_fs1)
bayes_R2(fit_fs2)
bayes_R2(fit_fs3)

summary(fit_fs2$fit)
fixef(fit_fs2)
ranef(fit_fs2)
neff_lowest(fit_fs2$fit)
rhat_highest(fit_fs2$fit)

slope <- as.data.frame(fit_fs2, variable = "b_length_avg_anom2")
hist(slope[ , 1], xlim = c(0, max(slope[ , 1])))
g <- ggplot(slope) +
    geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
    geom_histogram(aes(x = b_length_avg_anom2), bins = 20, color = "white",
                   fill = "red3", alpha = 0.5) +
    labs(x = "Common size effect (eggs / mm)", y = "Count",
         title = "Posterior distribution for size effect") +
    theme_simple()
print(g)
ggsave("./figures/age_size/posterior_size_effect.png", g,
       units = "in", width = 6, height = 4)



## Extract stock-specific slopes
slope_b  <- as.matrix(fit_fs2, variable = c("^b_.*length_avg_anom2"), regex = TRUE)
slope_r  <- as.matrix(fit_fs2, variable = c("^r_.*length_avg_anom2"), regex = TRUE)
for(i in 1:ncol(slope_r)) slope_r[ , i] <- slope_r[ , i] + slope_b
nam <- unique(fecund_size$id_label[!is.na(fecund_size$length_avg)])
colnames(slope_r) <- as.character(nam)
slope_post <- as.data.table(reshape2::melt(slope_r))

# g <- ggplot(slope_post) +
#     geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
#     aes(x = value, y = variable) +
#     labs(x = "Size effect (eggs / mm)", y = "") +
#     geom_violin() +
#     theme_simple()
# print(g)

slope_sum <- slope_post[ , .(mean = mean(value),
                             upper = quantile(value, probs = 0.975),
                             lower = quantile(value, probs = 0.025)), .(variable)]
g <- ggplot(slope_sum) +
    geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
    geom_vline(xintercept = mean(slope_b[ , 1]), color = "red3", linetype = 2) +
    geom_point(aes(x = mean, y = variable), size = 2) +
    geom_segment(aes(x = lower, xend = upper, y = variable, yend = variable)) +
    labs(x = "Size effect (eggs / mm)", y = "") +
    theme_simple()
print(g)
ggsave("./figures/age_size/posterior_size_effect_stock.png", g,
       units = "in", width = 6, height = 6)
