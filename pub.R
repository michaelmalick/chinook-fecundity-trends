# Publication numbers, tables, figures

dir.create("./figures/pub/", showWarnings = FALSE)


## Numbers -------------------------------------------------

## Number of stocks w/ fecundity data
length(unique(fecund$id_label))

## Number of stocks w/ size data
length(unique(fecund_size$id_label))

## Number of stocks in each region
fecund[ , .(N = length(unique(id_label))) , by = region]

## Number of stocks by run_timing
fecund[ , .(N = length(unique(id_label))) , by = run]

## Number of USFWS hatcheries
length(grep("NFH", levels(fecund$id_label)))

## Percent of stocks only release subyearlings
d <- fecund[ , .(release = unique(release)), by = id_number]
table(d$release) / nrow(d)


## Number of significant SSRW trends
sig_count <- fecund_rw_fit[ , .(count_sig90 = sum(sig_90 == 1)), by = id_label]
sum(sig_count$count_sig90 > 0)
sum(sig_count$count_sig90 == 0)


## Mean correlation
mean(cor_all$cor, na.rm = TRUE)
## Mean correlation w/o Upper Columbia Spring
mean(cor_all$cor[cor_all$group1 != "Upper Columbia Spring" &
                 cor_all$group2 != "Upper Columbia Spring"], na.rm = TRUE)
## Mean correlation within groups
cor_all[group1 == group2, .(mean = mean(cor, na.rm = TRUE))]
## Mean correlation among groups
cor_all[group1 != group2, .(mean = mean(cor, na.rm = TRUE))]
# Mean correlations for each group
cor_all[group1 == group2, .(mean = mean(cor, na.rm = TRUE)), by = group1]


## Fecundity DFA peak 1
dfa1_trends <- dfa_trends(dfa1rot, years = 1995:2019)
dfa1_trends$time[which.max(dfa1_trends$estimate)]
## Fecundity DFA peak 2
dfa1_trends2 <- dfa1_trends[dfa1_trends$time %in% 2005:2015, ]
dfa1_trends2$time[which.max(dfa1_trends2$estimate)]


## DFA size correlations
dfa1_age3trends <- dfa_trends(dfa1_age3rot, years = 1995:2019)
dfa1_age4trends <- dfa_trends(dfa1_age4rot, years = 1995:2019)
dfa1_age5trends <- dfa_trends(dfa1_age5rot, years = 1995:2019)
cor(dfa1_age3trends$estimate, dfa1_age4trends$estimate)
cor(dfa1_age3trends$estimate, dfa1_age5trends$estimate)
cor(dfa1_age4trends$estimate, dfa1_age5trends$estimate)

## DFA size peak fecundity
dfa1_age3trends$time[which.max(dfa1_age3trends$estimate)]
dfa1_age4trends$time[which.max(dfa1_age4trends$estimate)]
dfa1_age5trends$time[which.max(dfa1_age5trends$estimate)]


## Correlate fecundity and size DFA models
cor(dfa1_trends$estimate, dfa1_age3trends$estimate)
cor(dfa1_trends$estimate, dfa1_age4trends$estimate)
cor(dfa1_trends$estimate, dfa1_age5trends$estimate)

## Common size effect
slope <- as.data.frame(fit_fs2, variable = "b_length_avg_anom2")
mean(slope$b_length_avg_anom2)
quantile(slope$b_length_avg_anom2, probs = c(0.025, 0.975))

## HBM R^2
r2 <- as.data.frame(bayes_R2(fit_fs2, summary = FALSE))
mean(r2$R2)
quantile(r2$R2, probs = c(0.025, 0.975))


## % decline in fecundity
m <- fecund_rw_fit[ , .(fecund_2009 = fecundity[year == 2009],
                        fecund_2017 = fecundity[year == 2017]),
                   by = id_label]
m$diff <- m$fecund_2009 - m$fecund_2017
m$pc <- ((m$fecund_2017 - m$fecund_2009) / m$fecund_2009) * 100
mean(m$diff, na.rm = TRUE)
mean(m$pc, na.rm = TRUE)
range(m$pc, na.rm = TRUE)
median(m$pc, na.rm = TRUE)
m[id_label == "Spring Creek NFH Fall", ]

size_id <- unique(fecund_size$id_label)
m$size <- ifelse(m$id_label %in% size_id, 1, 0)
m_size <- m[size == 1]
mean(m_size$diff, na.rm = TRUE)


## % decline in size
s <- age_size[age %in% 3:5 , .(length_2009 = length[year == 2009],
                      length_2017 = length[year == 2017]),
              by = .(age, id_label)]
s$diff <- s$length_2009 - s$length_2017
s$pc <- ((s$length_2017 - s$length_2009) / s$length_2009) * 100
s[id_label == "Spring Creek NFH Fall", ]
s[ , .(mean = mean(pc, na.rm = TRUE),
       min = min(pc, na.rm = TRUE),
       max = max(pc, na.rm = TRUE),
       mean_diff = mean(diff, na.rm = TRUE)), by = age]


## Females spawned
m <- fecund[ , .(min = min(females_spawned, na.rm = TRUE),
                          max = max(females_spawned, na.rm = TRUE),
                          median = median(females_spawned, na.rm = TRUE),
                          mean = mean(females_spawned, na.rm = TRUE)), by = id_label]
range(m$mean)
mean(fecund$females_spawned, na.rm = TRUE)
range(m$median)
median(fecund$females_spawned, na.rm = TRUE)

plot(fecund$females_spawned)
hist(fecund$females_spawned, breaks = 100)



## MCMC diag -----------------------------------------------

## 1. Univariate SS-RW
min(ssrw_neff$neff)
max(ssrw_neff$rhat)


## 2. Fecundity DFA: 1-trend
is_converged(dfa1, 1.01, parameters = c("sigma", "x\\[", "Z\\["))
neff_lowest_dfa(dfa1)
rhat_highest_dfa(dfa1)
check_hmc_diagnostics(dfa1$model)


## 3. Size DFA: 1-trend
is_converged(dfa1_age3, 1.01, parameters = c("sigma", "x\\[", "Z\\["))
is_converged(dfa1_age4, 1.01, parameters = c("sigma", "x\\[", "Z\\["))
is_converged(dfa1_age5, 1.01, parameters = c("sigma", "x\\[", "Z\\["))

neff_lowest_dfa(dfa1_age3)
neff_lowest_dfa(dfa1_age4)
neff_lowest_dfa(dfa1_age5)

rhat_highest_dfa(dfa1_age3)
rhat_highest_dfa(dfa1_age4)
rhat_highest_dfa(dfa1_age5)

check_hmc_diagnostics(dfa1_age3$model)
check_hmc_diagnostics(dfa1_age4$model)
check_hmc_diagnostics(dfa1_age5$model)


## 4. Fecundity ~ size
neff_lowest(fit_fs2$fit)
rhat_highest(fit_fs2$fit)
check_hmc_diagnostics(fit_fs2$fit)



## Table: hatchery info ------------------------------------
tab <- fecund[ , .(id_label  = unique(na.omit(id_label)),
                   hatchery  = unique(na.omit(hatchery)),
                   region    = unique(na.omit(region)),
                   run       = unique(na.omit(run)),
                   release   = unique(na.omit(release)),
                   yr_min    = min(year[!is.na(fecundity)]),
                   yr_max    = max(year[!is.na(fecundity)]),
                   years     = 0,
                   N         = sum(!is.na(fecundity))),
              by = "id_number"]

sz <- as.character(unique(fecund_size$id_label))
sz <- data.table(id_label = sz, size_data = "y")
tab <- merge(tab, sz, all.x = TRUE)
tab$size_data <- ifelse(is.na(tab$size_data), "", "y")

tab[ , id_label := NULL]
tab <- tab[order(region, run), ]
tab$id_number <- 1:nrow(tab)

tab$years <- paste0(tab$yr_min, "--", tab$yr_max)
tab[ , yr_min := NULL]
tab[ , yr_max := NULL]

tab$release <- ifelse(tab$release == "subyearling", "Sub-yearling", tab$release)
tab$release <- ifelse(tab$release == "subyearling and yearling", "Mixed", tab$release)
tab$release <- ifelse(tab$release == "yearling", "Yearling", tab$release)

names(tab) <- c("Number", "Hatchery", "Region", "Run", "Release",
                "Years", "N", "Length")

data.table::fwrite(tab, "./figures/pub/table_data_info.csv")


## Table: females spawned info -----------------------------
tab <- fecund[ , .(hatchery  = unique(na.omit(hatchery)),
                   run       = unique(na.omit(run)),
                   min_females   = min(females_spawned),
                   med_females   = as.double(median(females_spawned)),
                   max_females   = max(females_spawned)
                   ),
              by = "id_number"]

names(tab) <- c("Number", "Hatchery", "Run",
                "Min. N", "Median N", "Max. N")

data.table::fwrite(tab, "./figures/pub/table_females_info.csv")


## Table: size info ----------------------------------------
ss <- age_size[age %in% 3:5, .(n_years = .N,
                               min = min(nsize),
                               median = as.double(median(nsize)),
                               max = max(nsize)
                               # p5 = (sum(nsize < 5) / .N) * 100
                               ) , by = .(id_number, hatchery, run, age)]
ss[order(id_number, age), ]

names(ss) <- c("Number", "Hatchery", "Run", "Age",
               "N years", "Min. N", "Median N", "Max. N")

data.table::fwrite(ss, "./figures/pub/table_size_info.csv")


## Table: DFA LOOIC ----------------------------------------
tab_dfa <- data.table(Trends = dfa_looic$trends,
                      LOOIC = dfa_looic$looic,
                      SE = dfa_looic$se)
data.table::fwrite(tab_dfa, "./figures/pub/table_dfa_looic.csv")



## Table: correlation means --------------------------------
gr <- as.vector(unique(cor_all$group1))
within <- rep(NA, length(gr))
between  <- rep(NA, length(gr))
nn  <- rep(NA, length(gr))
for(i in seq_along(gr)) {
    sub <- cor_all[group1 == gr[i] | group2 == gr[i]]
    wit <- sub[group1 == group2 | group2 == group1, ]
    amo <- sub[group1 != group2 | group2 != group1, ]
    # plot(sub$cor, main = gr[i])
    # abline(h = 0)
    # abline(h = mean(sub$cor, na.rm = TRUE), col = 2)
    within[i] <- mean(wit$cor, na.rm = TRUE)
    between[i]  <- mean(amo$cor, na.rm = TRUE)
    nn[i] <- fecund[group == gr[i], .(N = length(unique(id_label))),
                    by = .(group)]$N
}
tcor <- data.table(group = gr, N = nn,
                   within_group = within,
                   between_group = between)
tcor$within_group <- formatC(tcor$within_group, digits = 2, format = "f")
tcor$between_group <- formatC(tcor$between_group, digits = 2, format = "f")
tcor$within_group[tcor$within_group == "NaN"] <- "-"
names(tcor) <- c("Stock group", "N stocks",
                 "$\\bar{r}$ within group",
                 "$\\bar{r}$ across groups")
data.table::fwrite(tcor, "./figures/pub/table_cor_mean.csv")




## Fig: kalman region --------------------------------------
datl <- data.table(group = levels(fecund_kf$group))
datl$group <- factor(datl$group, levels = levels(fecund_kf$group))
for(i in 1:nrow(datl)) datl$lab[i] <- paste0(LETTERS[i], ": ", datl$group[i])
g <- ggplot(fecund_kf) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    geom_line(aes(x = year, y = kalman, group = id_label, color = region),
              na.rm = TRUE) +
    geom_text(data = datl, aes(x = 1995, y = 2.05, label = lab),
              color = "grey40", size = 2.7, hjust = 0) +
    labs(y = "Smoothed fecundity", x = "Year", color = "") +
    scale_color_manual(values = M1) +
    scale_x_continuous(limits = c(1995, 2020),
                       breaks = c(1995, 2005, 2015)) +
    theme_simple(grid = TRUE) +
    theme(legend.position = "top",
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.spacing.y = unit(-0.5, "pt"),
          panel.spacing.x = unit(-0.5, "pt")) +
    facet_wrap( ~ group, ncol = 3)
print(g)
ggsave("./figures/pub/kalman_group.png", units = "in",
       width = 6, height = 5, bg = "white")
ggsave("./figures/pub/kalman_group.pdf", units = "in",
       width = 6, height = 5, bg = "white")



## Fig: kalman stock ---------------------------------------
dat <- copy(fecund_kf)
dat[ , fecund_na := ifelse(is.na(fecundity_stnd), 2, 16)]
datl <- unique(dat[ , "id_label"])
g <- ggplot(dat) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    geom_ribbon(aes(ymin = lower, ymax = upper, x = year, fill = region),
                color = NA, alpha = 0.2) +
    geom_line(aes(x = year, y = kalman, color = region)) +
    geom_point(data = dat[is.na(fecundity_stnd), ],
               aes(x = year, y = kalman), color = "grey25", size = 0.8) +
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
ggsave("./figures/pub/kalman_stock.png", units = "in",
       width = 8, height = 7.5, bg = "white")
ggsave("./figures/pub/kalman_stock.pdf", units = "in",
       width = 8, height = 7.5, bg = "white")



## Fig: SS RW region ---------------------------------------
datl <- data.table(group = levels(fecund_rw_fit$group))
datl$group <- factor(datl$group, levels = levels(fecund_rw_fit$group))
for(i in 1:nrow(datl)) datl$lab[i] <- paste0(LETTERS[i], ": ", datl$group[i])
g <- ggplot(fecund_rw_fit) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    geom_line(aes(x = year, y = state, group = id_label, color = region),
              na.rm = TRUE) +
    geom_text(data = datl, aes(x = 1995, y = 2.05, label = lab),
              color = "grey40", size = 2.7, hjust = 0) +
    labs(y = "Fecundity (standardized)", x = "Year", color = "") +
    scale_color_manual(values = M1) +
    scale_x_continuous(limits = c(1995, 2020),
                       breaks = c(1995, 2005, 2015)) +
    theme_simple(grid = TRUE) +
    theme(legend.position = "top",
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.spacing.y = unit(-0.5, "pt"),
          panel.spacing.x = unit(-0.5, "pt")) +
    facet_wrap( ~ group, ncol = 3)
print(g)
ggsave("./figures/pub/ssrw_group.png", units = "in",
       width = 6, height = 5, bg = "white")
ggsave("./figures/pub/ssrw_group.pdf", units = "in",
       width = 6, height = 5, bg = "white")



## Fig: SS RW stock ----------------------------------------
dat <- copy(fecund_rw_fit)
dat[ , fecund_na := ifelse(is.na(fecundity_stnd), 2, 16)]
size_id <- as.character(unique(fecund_size$id_label))
dat$size <- ifelse(dat$id_label %in% size_id, "yes", "no")
datl <- unique(dat[ , "id_label"])

g <- ggplot(dat) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    geom_ribbon(aes(ymin = lower_90, ymax = upper_90, x = year, fill = region),
                color = NA, alpha = 0.2) +
    geom_line(aes(x = year, y = state, color = region, linetype = size)) +
    geom_point(data = dat[is.na(fecundity_stnd), ],
               aes(x = year, y = state), color = "grey25", size = 0.6) +
    geom_point(data = dat,
               aes(x = year, y = fecundity_stnd), color = "grey60", size = 0.6, na.rm = TRUE) +
    geom_text(data = datl, aes(x = 1995, y = 2.05, label = id_label),
              color = "grey40", size = 2.7, hjust = 0) +
    scale_x_continuous(limits = c(1995, 2020),
                       breaks = c(1995, 2005, 2015)) +
    labs(y = "Fecundity (standardized)", x = "Year", color = "", fill = "") +
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

ggsave("./figures/pub/ssrw_stock.png", units = "in",
       width = 8, height = 7.5, bg = "white")
ggsave("./figures/pub/ssrw_stock.pdf", units = "in",
       width = 8, height = 7.5, bg = "white")



## Fig: Fecundity DFA --------------------------------------

## Trends
rot1 <- rotate_trends(dfa1, conf_level = 0.95, invert = TRUE)
dfa1_trends <- dfa_trends(rot1, years = 1995:2019)
g1 <- ggplot(dfa1_trends) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    aes(x = time, y = estimate) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line() +
    labs(x = "Year", y = "Fecundity trend", title = "A") +
    scale_x_continuous(limits = c(1995, 2019), expand = c(0, 0),
                       breaks = c(1995, 2000, 2005, 2010, 2015)) +
    theme_simple(grid = TRUE)
print(g1)

f_max <- dfa1_trends[dfa1_trends$time == 2009, ]
f_min <- dfa1_trends[dfa1_trends$time == 2017, ]
(f_max$estimate - f_min$estimate)

## Loadings
dfa1_loadings <- dfa_loadings(rot1, names = rownames(y_dfa), conf_level = 0.90)
dfa1_loadings <- as.data.table(dfa1_loadings)
grps <- unique(fecund_dfa[ , .(run, region, group, id_label)])
names(grps) <- c("run", "region", "group", "name")
dfa1_loadings <- merge(dfa1_loadings, grps, by = "name")
dfa1_loadings[ , group := factor(group, levels = levels(fecund_dfa$group))]
dfa1_loadings[ , run := factor(run, levels = levels(fecund_dfa$run))]
dfa1_loadings[ , N := paste("n =", as.character(.N)), by = group]
dfa1_loadings[ , sig := ifelse(upper < 0 | lower > 0, 1, 0)]
dfa1_loadings[run == "Spring" & region == "Puget Sound", ]

g2 <- ggplot() +
    geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
    geom_point(data = dfa1_loadings,
               aes(x = median, y = group, fill = region,
                   color = region, shape = as.factor(sig)),
               size = 4, alpha = 0.5) +
    labs(x = "Loadings",
         y = NULL,
         fill = "",
         color = "", title = "B") +
    scale_color_manual(values = M1) +
    scale_fill_manual(values = M1) +
    scale_shape_manual(values = c(1, 16)) +
    theme_simple(grid = TRUE) +
    theme(legend.position = "none")
print(g2)

(g1 + theme(axis.title.y = element_text(margin = margin(r = -170, unit = "pt")))) / g2
ggsave("./figures/pub/fecundity_dfa.png", units = "in",
       width = 6, height = 7.5, bg = "white")
ggsave("./figures/pub/fecundity_dfa.pdf", units = "in",
       width = 6, height = 7.5, bg = "white")



## Stock-specific loadings
setorder(dfa1_loadings, "group")
dfa1_loadings[ , name := factor(name, levels = unique(name))]
dfa1_loadings[ , color := .GRP, by = group]
dfa1_loadings[ , diff := c(0, diff(color))]
dfa1_loadings[ , gr_line := ifelse(diff == 1, .I - 0.5, NA) , by = .I]

N <- nrow(dfa1_loadings)
lab <- rep(NA, nrow(dfa1_loadings))
lab[N] <- as.character(dfa1_loadings$group[N])
num <- rep(NA, nrow(dfa1_loadings))
num[N] <- N
for(i in (N-1):1) {
    if(dfa1_loadings$group[i] != dfa1_loadings$group[i+1]) {
        lab[i] <- as.character(dfa1_loadings$group[i])
        num[i] <- i
    }
}
dfa1_loadings$ylabel <- lab
dfa1_loadings$ynumber <- num

 g <- ggplot(dfa1_loadings) +
     geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
     geom_point(aes(x = median, y = name, color = region)) +
     geom_segment(aes(x = lower, xend = upper,
                      y = name, yend = name, color = region)) +
     geom_hline(yintercept = dfa1_loadings$gr_line,
                color = "grey70", na.rm = TRUE, size = 0.3) +
    geom_text(data = dfa1_loadings, aes(x = -1.25, y = ynumber, label = ylabel),
              color = "grey40", size = 2.7, hjust = 0, vjust = 0.3, na.rm = TRUE) +
     labs(x = "Loading", y = "", color = "") +
     scale_x_continuous(expand = c(0, 0)) +
     scale_color_manual(values = M1) +
     theme_simple() +
     theme(legend.position="none")
print(g)
ggsave("./figures/pub/fecundity_dfa_loadings.png", units = "in",
       width = 6, height = 7.5, bg = "white")
ggsave("./figures/pub/fecundity_dfa_loadings.pdf", units = "in",
       width = 6, height = 7.5, bg = "white")




## Fig: Size DFA -------------------------------------------
dfa1_age3trends <- dfa_trends(dfa1_age3rot, years = 1995:2019)
dfa1_age4trends <- dfa_trends(dfa1_age4rot, years = 1995:2019)
dfa1_age5trends <- dfa_trends(dfa1_age5rot, years = 1995:2019)
dfa1_age3trends$lab <- "Age 3"
dfa1_age4trends$lab <- "Age 4"
dfa1_age5trends$lab <- "Age 5"
dat <- as.data.table(rbind(dfa1_age3trends, dfa1_age4trends, dfa1_age5trends))
dat$lab <- ifelse(dat$lab == "Age 5", "A: age 5", dat$lab)
dat$lab <- ifelse(dat$lab == "Age 4", "B: age 4", dat$lab)
dat$lab <- ifelse(dat$lab == "Age 3", "C: age 3", dat$lab)
datl <- dat[ , .(max = max(upper)), by = lab]

g1 <- ggplot(dat) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    aes(x = time, y = estimate) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line() +
    geom_text(data = datl, aes(x = 2017, y = max, label = lab),
              color = "grey40", size = 2.7, hjust = 0.4) +
    labs(x = "Year", y = "Size trend") +
    scale_x_continuous(limits = c(1995, 2019), expand = c(0, 0),
                       breaks = c(1995, 2000, 2005, 2010, 2015)) +
    facet_wrap( ~ lab, ncol = 1) +
    theme_simple(grid = TRUE) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.spacing.y = unit(-0.5, "pt"),
          panel.spacing.x = unit(-0.5, "pt"))
print(g1)
ggsave("./figures/pub/size_dfa.png",
       width = 4, height = 6, bg = "white")
ggsave("./figures/pub/size_dfa.pdf",
       width = 4, height = 6, bg = "white")



## Loadings
dfa1_age3loadings <- dfa_loadings(dfa1_age3rot, names = rownames(y_dfa_age3), conf_level = 0.90)
dfa1_age4loadings <- dfa_loadings(dfa1_age4rot, names = rownames(y_dfa_age4), conf_level = 0.90)
dfa1_age5loadings <- dfa_loadings(dfa1_age5rot, names = rownames(y_dfa_age5), conf_level = 0.90)

dfa1_age3loadings$age <- "Age 3"
dfa1_age4loadings$age <- "Age 4"
dfa1_age5loadings$age <- "Age 5"

dfa1_age_loadings <- rbind(dfa1_age3loadings, dfa1_age4loadings, dfa1_age5loadings)
grps <- unique(size_dfa[ , .(run, region, group, id_label)])
names(grps) <- c("run", "region", "group", "name")
dfa1_age_loadings <- merge(dfa1_age_loadings, grps, by = "name")
dfa1_age_loadings$name <- factor(dfa1_age_loadings$name, levels = levels(age_size$id_label))

g <- ggplot(dfa1_age_loadings) +
    geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
    geom_point(aes(x = median, y = name)) +
    geom_segment(aes(x = lower, xend = upper, y = name, yend = name)) +
    labs(x = "Loading", y = "") +
    scale_x_continuous(expand = c(0, 0)) +
    theme_simple() +
    facet_wrap( ~ age) +
    theme(legend.position="none")
print(g)
ggsave("./figures/pub/size_dfa_loadings.png", units = "in",
       width = 7, height = 5, bg = "white")
ggsave("./figures/pub/size_dfa_loadings.pdf", units = "in",
       width = 7, height = 5, bg = "white")



## Fig: Correlations ---------------------------------------
g0 <- cor_plot_region(cor_all, title = NULL, legend_title = NULL, max_size = 2.5)
ggsave("./figures/pub/cor.png", units = "in",
       width = 6, height = 3.8, bg = "white")
ggsave("./figures/pub/cor.pdf", units = "in",
       width = 6, height = 3.8, bg = "white")


# g0 <- cor_plot_region(cor_all, title = "All years", max_size = 2)
# g1 <- cor_plot_region(cor_early, title = "1995-2007", group_labels = FALSE, max_size = 1)
# g2 <- cor_plot_region(cor_late, title = "2008-2019", group_labels = FALSE, max_size = 1)
#
# g <- (g0 | (g1 / g2)) + plot_layout(widths = c(5, 3), guides = 'collect')
# print(g)
# ggsave("./figures/pub/cor.png", units = "in",
#        width = 7, height = 4, bg = "white")



## Fig: Size yr effect -------------------------------------
dat <- copy(size_trend)
dat[ , lab := paste("Age", age)]
dat$lab <- ifelse(dat$lab == "Age 5", "A: age 5", dat$lab)
dat$lab <- ifelse(dat$lab == "Age 4", "B: age 4", dat$lab)
dat$lab <- ifelse(dat$lab == "Age 3", "C: age 3", dat$lab)
datl <- dat[ , .(max = max(upper)), by = lab]
g <- ggplot(dat) +
    aes(x = year, y = median) +
    geom_point() +
    geom_segment(aes(x = year, xend = year, y = upper, yend = lower)) +
    geom_smooth(method = "loess", formula = y ~ x, span = 0.5, se = FALSE,
                color = "grey60", size = 0.5) +
    geom_text(data = datl, aes(x = 2017, y = max, label = lab),
              color = "grey40", size = 2.7, hjust = 0) +
    labs(x = "Year", y = "Length (mm)") +
    scale_x_continuous(limits = c(1995, 2019),
                       breaks = c(1995, 2000, 2005, 2010, 2015)) +
    facet_wrap( ~ lab, ncol = 1, scale = "free_y") +
    theme_simple() +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.spacing.y = unit(-0.5, "pt"),
          panel.spacing.x = unit(-0.5, "pt"))
print(g)
ggsave("./figures/pub/size_year_effects.png",
       width = 4, height = 6, bg = "white")
ggsave("./figures/pub/size_year_effects.pdf",
       width = 4, height = 6, bg = "white")


## Fig: fecundity ~ size -----------------------------------

## Dotplot: Stock-specific slopes
slope_b  <- as.matrix(fit_fs2, variable = c("^b_.*length_avg_anom2"), regex = TRUE)
slope_r  <- as.matrix(fit_fs2, variable = c("^r_.*length_avg_anom2"), regex = TRUE)
for(i in 1:ncol(slope_r)) slope_r[ , i] <- slope_r[ , i] + slope_b
nam <- unique(fecund_size$id_label[!is.na(fecund_size$length_avg)])
colnames(slope_r) <- as.character(nam)
slope_post <- as.data.table(reshape2::melt(slope_r))
slope_sum <- slope_post[ , .(mean = mean(value),
                             upper = quantile(value, probs = 0.975),
                             lower = quantile(value, probs = 0.025)), .(variable)]
g1 <- ggplot(slope_sum) +
    # geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
    geom_vline(xintercept = mean(slope_b[ , 1]), color = "red3", linetype = 2) +
    geom_point(aes(x = mean, y = variable), size = 2) +
    scale_x_continuous(limits = c(-0.3, 14),
                       breaks = seq(0, 14, by = 2)) +
    geom_segment(aes(x = lower, xend = upper, y = variable, yend = variable)) +
    labs(x = "Size effect (eggs / mm)", y = NULL, title = "A") +
    theme_simple()
print(g1)


## Histogram: Common size effect
slope <- as.data.frame(fit_fs2, variable = "b_length_avg_anom2")
range(unlist(slope))
g2 <- ggplot(slope) +
    # geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
    geom_histogram(aes(x = b_length_avg_anom2), bins = 20, color = "white",
                   fill = "red3", alpha = 0.5, na.rm = TRUE) +
    # scale_x_continuous(limits = c(-0.3, 14),
    #                    breaks = seq(0, 14, by = 2)) +
    labs(x = "Size effect (eggs / mm)", y = "Posterior samples", title = "C") +
    theme_simple()
print(g2)


## R2
r2 <- as.data.frame(bayes_R2(fit_fs2, summary = FALSE))
g3 <- ggplot(r2) +
    # geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
    geom_histogram(aes(x = R2), bins = 20, color = "white",
                   fill = "grey30", alpha = 0.5, na.rm = TRUE) +
    # scale_x_continuous(limits = c(-0.3, 14),
    #                    breaks = seq(0, 14, by = 2)) +
    labs(x = bquote(R^2), y = "Posterior density") +
    theme_simple()
print(g3)
ggsave("./figures/pub/fecundity_size_r2.png", g3,
       width = 5, height = 4, bg = "white")
ggsave("./figures/pub/fecundity_size_r2.pdf", g3,
       width = 5, height = 4, bg = "white")

## Scatterplot

## get common regression line
fmean <- mean(fecund_size$fecundity, na.rm = TRUE)
lmean <- mean(fecund_size$length_avg, na.rm = TRUE)
r <- range(fecund_size$length_avg_anom2, na.rm = TRUE)
l <- seq(r[1], r[2], length.out = 100)
nd <- data.table(length_avg_anom2 = l)
pe <- posterior_epred(fit_fs2, newdata = nd, re_formula = NA)
pr <- data.table(y = apply(pe, 2, median) + fmean,
                 y_lower = apply(pe, 2, function(x) quantile(x, probs = 0.025)) + fmean,
                 y_upper = apply(pe, 2, function(x) quantile(x, probs = 0.975)) + fmean,
                 x = l + lmean)

g4 <- ggplot(fecund_size) +
    geom_point(aes(x = length_avg, y = fecundity), color = "grey70", size = 0.6,
               na.rm = TRUE) +
    # geom_path(data = df, aes(x = x, y = y, group = id_label), color = "grey70",
    #           size = 0.7) +
    geom_path(data = pr, aes(x = x, y = y), color = "red3", size = 1.2) +
    geom_ribbon(data = pr, aes(x = x, ymin = y_lower, ymax = y_upper),
                fill = "grey30", color = NA, alpha = 0.2) +
    labs(x = "Length (mm)", y = "Fecundity", title = "B") +
    theme_simple()
print(g4)



g1 + (g4 / g2) + plot_layout(widths = c(4, 5))
ggsave("./figures/pub/fecundity_size.png",
       width = 7, height = 5, bg = "white")
ggsave("./figures/pub/fecundity_size.pdf",
       width = 7, height = 5, bg = "white")


## Fig: fecundity ~ size by stock --------------------------

## get stock-specific regression lines
fmean <- mean(fecund_size$fecundity, na.rm = TRUE)
lmean <- mean(fecund_size$length_avg, na.rm = TRUE)
sp <- split(fecund_size, by = "id_number")
lst <- lapply(sp, function(x) {
    r <- range(x$length_avg_anom2, na.rm = TRUE)
    l <- seq(r[1], r[2], length.out = 100)
    id <- as.character(unique(x$id_label))
    nd <- data.table(length_avg_anom2 = l, id_label = id)
    pe <- posterior_epred(fit_fs2, newdata = nd)
    pr <- data.table(y = apply(pe, 2, median) + fmean,
                     y_lower = apply(pe, 2, function(x)
                                     quantile(x, probs = 0.025)) + fmean,
                     y_upper = apply(pe, 2, function(x)
                                     quantile(x, probs = 0.975)) + fmean,
                     x = l + lmean,
                     id_label = id)
    pr
})
df <- rbindlist(lst)
df$id_label <- factor(df$id_label, levels = levels(fecund_size$id_label))


g <- ggplot(fecund_size) +
    geom_point(aes(x = length_avg, y = fecundity), color = "grey50", na.rm = TRUE) +
    geom_path(data = df, aes(x = x, y = y), color = "red3", size = 1) +
    geom_ribbon(data = df, aes(x = x, ymin = y_lower, ymax = y_upper),
                fill = "grey30", color = NA, alpha = 0.2) +
    facet_wrap( ~ id_label, scale = "free", ncol = 4) +
    labs(x = "Length (mm)", y = "Fecundity") +
    theme_simple(base_size = 10)
print(g)

ggsave("./figures/pub/fecundity_size_stock.png",
       width = 8, height = 7.5, bg = "white")
ggsave("./figures/pub/fecundity_size_stock.pdf",
       width = 8, height = 7.5, bg = "white")



## Fig: fecundity ~ size single stock ----------------------

## get stock-specific regression lines
fmean <- mean(fecund_size$fecundity, na.rm = TRUE)
lmean <- mean(fecund_size$length_avg, na.rm = TRUE)
sp <- split(fecund_size, by = "id_number")
lst <- vector("list", length(sp))
for(i in seq_along(sp)) {
    x <- sp[[i]]
    r <- range(x$length_avg_anom2, na.rm = TRUE)
    l <- seq(r[1], r[2], length.out = 100)
    id <- as.character(unique(x$id_label))
    nd <- data.table(length_avg_anom2 = l, id_label = id)
    pe <- posterior_epred(lm_fits[[i]], newdata = nd)
    pr <- data.table(y = apply(pe, 2, median) + fmean,
                     y_lower = apply(pe, 2, function(x)
                                     quantile(x, probs = 0.025)) + fmean,
                     y_upper = apply(pe, 2, function(x)
                                     quantile(x, probs = 0.975)) + fmean,
                     x = l + lmean,
                     id_label = id)
    lst[[i]] <- pr
}
df <- rbindlist(lst)
df$id_label <- factor(df$id_label, levels = levels(fecund_size$id_label))


g <- ggplot(fecund_size) +
    geom_point(aes(x = length_avg, y = fecundity), color = "grey50", na.rm = TRUE) +
    geom_path(data = df, aes(x = x, y = y), color = "red3", size = 1) +
    geom_ribbon(data = df, aes(x = x, ymin = y_lower, ymax = y_upper),
                fill = "grey30", color = NA, alpha = 0.2) +
    facet_wrap( ~ id_label, scale = "free", ncol = 4) +
    labs(x = "Length (mm)", y = "Fecundity") +
    theme_simple(base_size = 10)
print(g)

ggsave("./figures/pub/fecundity_size_single_stock.png",
       width = 8, height = 7.5, bg = "white")
ggsave("./figures/pub/fecundity_size_single_stock.pdf",
       width = 8, height = 7.5, bg = "white")



## Fig: DFA trends: fecundity ~ size -----------------------
dfa1_trends <- dfa_trends(dfa1rot, years = 1995:2019)
dfa1_fecund <- data.frame(time = dfa1_trends$time,
                          fecund_mean = dfa1_trends$estimate,
                          fecund_lo = dfa1_trends$lower,
                          fecund_hi = dfa1_trends$upper)
dfa1_age3trends <- dfa_trends(dfa1_age3rot, years = 1995:2019)
dfa1_age4trends <- dfa_trends(dfa1_age4rot, years = 1995:2019)
dfa1_age5trends <- dfa_trends(dfa1_age5rot, years = 1995:2019)
dfa1_age3trends$age <- "Age 3"
dfa1_age4trends$age <- "Age 4"
dfa1_age5trends$age <- "Age 5"

dfa1_age <- rbind(dfa1_age3trends, dfa1_age4trends, dfa1_age5trends)
dfa_merge <- merge(dfa1_age, dfa1_fecund, sort = FALSE)

g1 <- ggplot(dfa_merge) +
    geom_segment(aes(x = estimate, xend = estimate, y = fecund_hi, yend = fecund_lo),
                 color = "grey70", size = 0.3) +
    geom_segment(aes(x = lower, xend = upper, y = fecund_mean, yend = fecund_mean),
                 color = "grey70", size = 0.3) +
    geom_point(aes(x = estimate, y = fecund_mean)) +
    labs(x = "Body length trend", y = "Fecundity trend") +
    facet_wrap( ~ age, ncol = 1, as.table = FALSE) +
    theme_simple()
print(g1)


## Correlation posteriors
trend_post <- function(rotated_trends) {
    ft <- rotated_trends$trends
    nyrs <- dim(ft)[3]
    niter <- dim(ft)[1]
    ft_mat <- matrix(NA, nrow = niter, ncol = nyrs)
    for(i in 1:nyrs) ft_mat[ , i] <- ft[ , 1, i]
    return(ft_mat)
}

fecund_mat <- trend_post(dfa1rot)
age3_mat <- trend_post(dfa1_age3rot)
age4_mat <- trend_post(dfa1_age4rot)
age5_mat <- trend_post(dfa1_age5rot)
plot(fecund_mat[1, ])

age3_cor <- rep(NA, nrow(fecund_mat))
age4_cor <- rep(NA, nrow(fecund_mat))
age5_cor <- rep(NA, nrow(fecund_mat))

for(i in 1:nrow(fecund_mat)) {
    age3_cor[i] <- cor(fecund_mat[i, ], age3_mat[i, ])
    age4_cor[i] <- cor(fecund_mat[i, ], age4_mat[i, ])
    age5_cor[i] <- cor(fecund_mat[i, ], age5_mat[i, ])
}

## Some iterations don't have properly inverted trends, fix these
plot(age3_cor, ylim = c(-1, 1)); abline(h = 0, col = 2)
plot(age4_cor, ylim = c(-1, 1)); abline(h = 0, col = 2)
plot(age5_cor, ylim = c(-1, 1)); abline(h = 0, col = 2)

if(all(age3_cor < 0)) age3_cor <- age3_cor * -1
if(all(age4_cor < 0)) age4_cor <- age4_cor * -1
if(all(age5_cor < 0)) age5_cor <- age5_cor * -1

age3_cor_df <- data.table(cor = age3_cor, age = "Age 3")
age4_cor_df <- data.table(cor = age4_cor, age = "Age 4")
age5_cor_df <- data.table(cor = age5_cor, age = "Age 5")
age_cor <- rbind(age3_cor_df, age4_cor_df, age5_cor_df)
age_cor_mean <- age_cor[cor > 0, .(mean_cor = mean(cor)), by = age]

g2 <- ggplot(age_cor) +
    geom_histogram(aes(x = cor), bins = 50, color = "white",
                   fill = "grey50", size = 0.2) +
    geom_vline(data = age_cor_mean, aes(xintercept = mean_cor),
               color = "red3", linetype = 1) +
    scale_y_continuous(expand = c(0, 25)) +
    scale_x_continuous(breaks = seq(0.4, 1.0 , 0.1)) +
    labs(x = "Correlation coefficient", y = "Posterior density") +
    facet_wrap( ~ age, ncol = 1, as.table = FALSE) +
    theme_simple()
print(g2)


g1 + g2
ggsave("./figures/pub/dfa_fecundity_size_cor.png",
       width = 6, height = 7, bg = "white")
ggsave("./figures/pub/dfa_fecundity_size_cor.pdf",
       width = 6, height = 7, bg = "white")


## Fig: DFA trends: fecundity size only --------------------

## Trends
dfa1_trends_s <- dfa_trends(dfa1rot_s, years = 1995:2019)
g1 <- ggplot(dfa1_trends_s) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    aes(x = time, y = estimate) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line() +
    labs(x = "Year", y = "Fecundity trend", title = "A") +
    scale_x_continuous(limits = c(1995, 2019), expand = c(0, 0),
                       breaks = c(1995, 2000, 2005, 2010, 2015)) +
    theme_simple(grid = TRUE)
print(g1)


## Loadings
dfa1_loadings_s <- dfa_loadings(dfa1rot_s, names = rownames(y_dfa_s), conf_level = 0.90)
dfa1_loadings_s <- as.data.table(dfa1_loadings_s)
grps <- unique(fecund_dfa[ , .(run, region, group, id_label)])
names(grps) <- c("run", "region", "group", "name")
dfa1_loadings_s <- merge(dfa1_loadings_s, grps, by = "name")
dfa1_loadings_s[ , group := factor(group, levels = levels(fecund_dfa[size == "yes"]$group))]
dfa1_loadings_s[ , run := factor(run, levels = levels(fecund_dfa[size == "yes"]$run))]
dfa1_loadings_s[ , N := paste("n =", as.character(.N)), by = group]
dfa1_loadings_s[ , sig := ifelse(upper < 0 | lower > 0, 1, 0)]
dfa1_loadings_s[run == "Spring" & region == "Puget Sound", ]

g2 <- ggplot() +
    geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
    geom_point(data = dfa1_loadings_s,
               aes(x = median, y = group, fill = region,
                   color = region, shape = as.factor(sig)),
               size = 4, alpha = 0.5) +
    labs(x = "Loadings",
         y = NULL,
         fill = "",
         color = "", title = "B") +
    scale_color_manual(values = M1) +
    scale_fill_manual(values = M1) +
    scale_shape_manual(values = c(1, 16)) +
    theme_simple(grid = TRUE) +
    theme(legend.position = "none")
print(g2)

(g1 + theme(axis.title.y = element_text(margin = margin(r = -170, unit = "pt")))) / g2
ggsave("./figures/pub/fecundity_dfa_s.png", units = "in",
       width = 6, height = 7.5, bg = "white")
ggsave("./figures/pub/fecundity_dfa_s.pdf", units = "in",
       width = 6, height = 7.5, bg = "white")


## Fig: DFA trends: unequal trend --------------------------

## Trends
dfa1rot_u <- rotate_trends(dfa1_unequal, conf_level = 0.95, invert = FALSE)
dfa1_trends_u <- dfa_trends(dfa1rot_s, years = 1995:2019)
g1 <- ggplot(dfa1_trends_u) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    aes(x = time, y = estimate) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line() +
    labs(x = "Year", y = "Fecundity trend", title = "A") +
    scale_x_continuous(limits = c(1995, 2019), expand = c(0, 0),
                       breaks = c(1995, 2000, 2005, 2010, 2015)) +
    theme_simple(grid = TRUE)
print(g1)


## Loadings
dfa1_loadings_u <- dfa_loadings(dfa1rot_u, names = rownames(y_dfa), conf_level = 0.90)
dfa1_loadings_u <- as.data.table(dfa1_loadings_u)
grps <- unique(fecund_dfa[ , .(run, region, group, id_label)])
names(grps) <- c("run", "region", "group", "name")
dfa1_loadings_u <- merge(dfa1_loadings_u, grps, by = "name")
dfa1_loadings_u[ , group := factor(group, levels = levels(fecund_dfa[size == "yes"]$group))]
dfa1_loadings_u[ , run := factor(run, levels = levels(fecund_dfa[size == "yes"]$run))]
dfa1_loadings_u[ , N := paste("n =", as.character(.N)), by = group]
dfa1_loadings_u[ , sig := ifelse(upper < 0 | lower > 0, 1, 0)]
dfa1_loadings_u[run == "Spring" & region == "Puget Sound", ]

g2 <- ggplot() +
    geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
    geom_point(data = dfa1_loadings_u,
               aes(x = median, y = group, fill = region,
                   color = region, shape = as.factor(sig)),
               size = 4, alpha = 0.5) +
    labs(x = "Loadings",
         y = NULL,
         fill = "",
         color = "", title = "B") +
    scale_color_manual(values = M1) +
    scale_fill_manual(values = M1) +
    scale_shape_manual(values = c(1, 16)) +
    theme_simple(grid = TRUE) +
    theme(legend.position = "none")
print(g2)

(g1 + theme(axis.title.y = element_text(margin = margin(r = -170, unit = "pt")))) / g2
ggsave("./figures/pub/fecundity_dfa_unequal.png", units = "in",
       width = 6, height = 7.5, bg = "white")
ggsave("./figures/pub/fecundity_dfa_unequal.pdf", units = "in",
       width = 6, height = 7.5, bg = "white")


## Fig: Correlations -- size only --------------------------
cor_all_s <- copy(cor_all)
size_id <- as.character(unique(fecund_size$id_label))
cor_all_s$size1 <- ifelse(cor_all_s$id_label1 %in% size_id, "yes", "no")
cor_all_s$size2 <- ifelse(cor_all_s$id_label2 %in% size_id, "yes", "no")
cor_all_s <- cor_all_s[size1 == "yes" & size2 == "yes", ]

g0 <- cor_plot_region(cor_all_s, title = NULL, legend_title = NULL, max_size = 2.5)
ggsave("./figures/pub/cor_s.png", units = "in",
       width = 6, height = 3.8, bg = "white")
ggsave("./figures/pub/cor_s.pdf", units = "in",
       width = 6, height = 3.8, bg = "white")


# g0 <- cor_plot_region(cor_all, title = "All years", max_size = 2)
# g1 <- cor_plot_region(cor_early, title = "1995-2007", group_labels = FALSE, max_size = 1)
# g2 <- cor_plot_region(cor_late, title = "2008-2019", group_labels = FALSE, max_size = 1)
#
# g <- (g0 | (g1 / g2)) + plot_layout(widths = c(5, 3), guides = 'collect')
# print(g)
# ggsave("./figures/pub/cor.png", units = "in",
#        width = 7, height = 4, bg = "white")



## Fig: predictive checks ----------------------------------

## prior checks
pp = predict(fs2_prior, summary = FALSE)
y = fs2_prior$data$fecundity_anom2

y_df = data.table(source = 'y', draw = 0, y = y)
lst = vector("list", 50)
for(i in 51:100) lst[[i]] = data.table(source = 'y-rep', draw = i, y = pp[i, ])
pp_df = rbindlist(lst)

g1 <- ggplot(pp_df) +
    geom_density(aes(x = y, group = draw), color = "grey75", size = 0.25) +
    geom_density(data = y_df, aes(x = y), size = 0.75) +
    xlim(-1e4, 1e4) +
    labs(title = "A: prior predictive check",
         y = "Density") +
    theme_simple()
print(g1)


## posterior checks
pp = predict(fit_fs2, summary = FALSE)
y = fit_fs2$data$fecundity_anom2

y_df = data.table(source = 'y', draw = 0, y = y)
lst = vector("list", 50)
for(i in 51:100) lst[[i]] = data.table(source = 'y-rep', draw = i, y = pp[i, ])
pp_df = rbindlist(lst)

g2 <- ggplot(pp_df) +
    geom_density(aes(x = y, group = draw), color = "grey75", size = 0.25) +
    geom_density(data = y_df, aes(x = y), size = 0.75) +
    labs(title = "B: posterior predictive check",
         y = "") +
    theme_simple()
print(g2)


g = g1 + g2
print(g)

ggsave("./figures/pub/ppc.png", units = "in",
       width = 6, height = 4, bg = "white")
ggsave("./figures/pub/ppc.pdf", units = "in",
       width = 6, height = 4, bg = "white")


