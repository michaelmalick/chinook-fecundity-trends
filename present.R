## Presentation graphics

dir.create("./figures/present/", showWarnings = FALSE)


## Numbers -------------------------------------------------

x <- fecund[ , .(N = length(year)), by = id_label]
mean(x$N)

x <- fecund_size[ , .(N = length(year)), by = id_label]
mean(x$N)

## Avg decline 2009--2017
m <- fecund[ , .(fecund_2009 = fecundity[year == 2009],
                        fecund_2017 = fecundity[year == 2017]),
                   by = id_label]
m$diff <- m$fecund_2009 - m$fecund_2017
m$pc <- ((m$fecund_2017 - m$fecund_2009) / m$fecund_2009) * 100
mean(m$diff, na.rm = TRUE)
mean(m$pc, na.rm = TRUE)



## Fecundity violin ----------------------------------------
gr <-fecund[ , .N, by = .(group, region)]
gr$fecundity <- 1800
gd <- ggplot(fecund[fecundity > 2000, ]) +
    geom_vline(xintercept = mean(fecund$fecundity, na.rm = TRUE),
               color = "grey50", linetype = 2) +
    aes(x = fecundity, y = group, fill = region, color = region) +
    # geom_jitter(height = 0.2, na.rm = TRUE, size = 0.3) +
    geom_violin(na.rm = TRUE, alpha = 0.5, adjust = 1.3, size = 0.4) +
    geom_text(data = gr, aes(x = fecundity, y = group, label = N), size = 3.5) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 1,
                 color = "black", na.rm = TRUE) +
    labs(x = "Fecundity",
         y = NULL) +
    xlim(1800, 7100) +
    scale_color_manual(values = M1) +
    scale_fill_manual(values = M1) +
    theme_simple(base_size = 14) +
    theme(legend.position = "none")
print(gd)

ggsave("./figures/present/fecund_violin.png", width = 5, height = 5)



## Sim: Random walk ----------------------------------------
set.seed(124)

x <- cumsum(rnorm(100)) + rnorm(100)
plot(x)
f <- StructTS(x, type = "level")
s <- as.vector(tsSmooth(f))
lines(s)
df <- data.frame(time = 1:100, y = x, x = s)


g <- ggplot(df) +
    geom_point(aes(x = time, y = y), color = "steelblue") +
    geom_line(aes(x = time, y = x), color = "red3", size = 1) +
    labs(y = "Fecundity", x = "Time") +
    theme_simple(grid = TRUE) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank())
print(g)
ggsave("./figures/present/sim_rw.png", width = 7, height = 3)



## Sim: DFA ------------------------------------------------
set.seed(124)

time <- 1:100

y1 <- cumsum(rnorm(100))
y2 <- cumsum(rnorm(100))
y3 <- cumsum(rnorm(100))
y4 <- cumsum(rnorm(100))

g1 <- rep("Stock 1", 100)
g2 <- rep("Stock 2", 100)
g3 <- rep("Stock 3", 100)
g4 <- rep("Stock 4", 100)

df <- data.frame(y = c(y1, y2, y3, y4),
                 time = c(time, time, time, time),
                 g = c(g1, g2, g3, g4))

g <- ggplot(df) +
    geom_line(aes(x = time, y = y), size = 1, color = "steelblue") +
    labs(y = "Fecundity", x = "Time") +
    facet_wrap( ~ g, ncol = 1) +
    theme_simple() +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.spacing.y = unit(-0.5, "pt"),
          panel.spacing.x = unit(-0.5, "pt"))
print(g)
ggsave("./figures/present/sim_dfa_y.png", width = 4, height = 4)




f1 <- StructTS(y1, type = "level")
f2 <- StructTS(y2, type = "level")

s1 <- as.vector(tsSmooth(f1))
s2 <- as.vector(tsSmooth(f2))

g1 <- rep("Trend 1", 100)
g2 <- rep("Trend 2", 100)

df <- data.frame(y = c(s1, s2),
                 time = c(time, time),
                 g = c(g1, g2))

g <- ggplot(df) +
    geom_line(aes(x = time, y = y), size = 1, color = "red3") +
    labs(y = "", x = "Time") +
    facet_wrap( ~ g, ncol = 1) +
    theme_simple() +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank())
print(g)
ggsave("./figures/present/sim_dfa_x.png", width = 4, height = 4)



## Normal dist ---------------------------------------------
set.seed(100)

plot(density(rnorm(1e6)))

png("./figures/present/normal_curve.png", res = 300, width = 3, height = 3,
    units = "in")
    plot(dnorm(seq(-4, 4, 0.01)), type = "l", col = "red3", lwd = 3,
         xlab = "", ylab = "", xaxt="n", yaxt = "n", ann = FALSE)
    box(col = "white")
dev.off()



## Fecundity ~ size: methods -------------------------------

id <- unique(fecund_size$id_number)[4:7]
d <- fecund_size[id_number %in% id, ]

g <- ggplot(d) +
    aes(x = length_avg, y = fecundity) +
    geom_point(na.rm = TRUE, color = "steelblue") +
    geom_smooth(method = "lm", se = FALSE, color = "red3", formula = y ~ x,
                na.rm = TRUE) +
    facet_wrap( ~ id_number, nrow = 1) +
    labs(x = "Length", y = "Fecundity") +
    theme_simple() +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.spacing.y = unit(-0.5, "pt"),
          panel.spacing.x = unit(-0.5, "pt"))
print(g)
ggsave("./figures/present/hbm_methods.png", width = 9, height = 3)


## Skagit N ------------------------------------------------
## From @WDFW2017 pg. 141

## number of spawners
yr <- 1992:2016
n_skagit <- c(7348,
              5801,
              5561,
              6892,
              10613,
              4872,
              14609,
              4924,
              16930,
              13793,
              19591,
              9777,
              23553,
              20803,
              20768,
              11281,
              11664,
              6979,
              8017,
              5510,
              13817,
              10882,
              10480,
              13076,
              19388)
df <- data.frame(year = yr, spawners = n_skagit)

g <- ggplot(df) +
    aes(x = year, y = spawners) +
    geom_point(color = "steelblue", size = 2, na.rm = TRUE) +
    geom_line(color = "steelblue", size = 1, na.rm = TRUE) +
    labs(x = "Year", y = "Spawner abundance") +
    theme_simple(base_size = 14) +
    theme(plot.margin = margin(0, 10, 0, 0))
print(g)
ggsave("./figures/present/skagit_abundance.png", width = 7, height = 4)


## Fecundity time series -----------------------------------

unique(fecund$id_label)

g <- ggplot(fecund[id_label == "Soos Creek Fall", ]) +
    aes(x = year, y = fecundity) +
    geom_point(color = "steelblue", size = 2) +
    geom_line(color = "steelblue", size = 1) +
    labs(x = "Year", y = "Fecundity") +
    theme_simple(base_size = 14) +
    theme(plot.margin = margin(0, 10, 0, 0))
print(g)
ggsave("./figures/present/fecundity_ts_green.png", width = 7, height = 4)

g <- ggplot(fecund[id_label == "Marblemount Summer", ]) +
    aes(x = year, y = fecundity) +
    geom_point(color = "steelblue", size = 2) +
    geom_line(color = "steelblue", size = 1) +
    labs(x = "Year", y = "Fecundity") +
    theme_simple(base_size = 14) +
    theme(plot.margin = margin(0, 10, 0, 0))
print(g)
ggsave("./figures/present/fecundity_ts_skagit.png", width = 7, height = 4)

n_mean <- mean(n_skagit)
f <- fecund[id_label == "Marblemount Summer", ]
eggs_max <- f$fecundity[f$year == 2009] * (n_mean / 2)
eggs_min <- f$fecundity[f$year == 2018] * (n_mean / 2)

f$fecundity[f$year == 2009] - f$fecundity[f$year == 2018]

(eggs_max - eggs_min) / 1e6
eggs_max / eggs_min


((eggs_max - eggs_min) / eggs_max) * 100

## Fig: Size DFA -------------------------------------------
dfa1_age3trends <- dfa_trends(dfa1_age3rot, years = 1995:2019)
dfa1_age4trends <- dfa_trends(dfa1_age4rot, years = 1995:2019)
dfa1_age5trends <- dfa_trends(dfa1_age5rot, years = 1995:2019)
dfa1_age3trends$lab <- "Age 3"
dfa1_age4trends$lab <- "Age 4"
dfa1_age5trends$lab <- "Age 5"
dat <- as.data.table(rbind(dfa1_age3trends, dfa1_age4trends, dfa1_age5trends))
dat$lab <- ifelse(dat$lab == "Age 5", "Age 5", dat$lab)
dat$lab <- ifelse(dat$lab == "Age 4", "Age 4", dat$lab)
dat$lab <- ifelse(dat$lab == "Age 3", "Age 3", dat$lab)
datl <- dat[ , .(max = max(upper)), by = lab]

g1 <- ggplot(dat) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    aes(x = time, y = estimate, color = lab) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = lab), color = NA, alpha = 0.2) +
    geom_line() +
    geom_text(data = datl, aes(x = 2017, y = max, label = lab),
              color = "grey40", size = 3.5, hjust = 0.4) +
    labs(x = "Year", y = "Size trend") +
    scale_x_continuous(limits = c(1995, 2019), expand = c(0, 0),
                       breaks = c(1995, 2000, 2005, 2010, 2015)) +
    scale_color_manual(values = c("#239DD7", "#CD2343", "#C28500")) +
    scale_fill_manual(values = c("#239DD7", "#CD2343", "#C28500")) +
    facet_wrap( ~ lab, ncol = 1) +
    theme_simple(grid = TRUE) +
    theme(legend.position="none") +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.spacing.y = unit(-0.5, "pt"),
          panel.spacing.x = unit(-0.5, "pt"))
print(g1)
ggsave("./figures/present/size_dfa.png",
       width = 4, height = 6)



## Fig: fecundity ~ size -----------------------------------

## Dotplot: Stock-specific slopes
slope_b  <- as.matrix(fit_fs2, pars = c("^b_.*length_avg_anom2"))
slope_r  <- as.matrix(fit_fs2, pars = c("^r_.*length_avg_anom2"))
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
    labs(x = "Size effect (eggs / mm)", y = NULL) +
    theme_simple()
print(g1)
ggsave("./figures/present/fecundity_size_dot.png", width = 4, height = 4)


## Histogram: Common size effect
slope <- as.data.frame(fit_fs2, pars = "length_avg_anom2")
range(unlist(slope))
g2 <- ggplot(slope) +
    # geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
    geom_histogram(aes(x = b_length_avg_anom2), bins = 20, color = "white",
                   fill = "red3", alpha = 0.5, na.rm = TRUE) +
    # scale_x_continuous(limits = c(-0.3, 14),
    #                    breaks = seq(0, 14, by = 2)) +
    labs(x = "Size effect (eggs / mm)", y = "Posterior samples") +
    theme_simple(base_size = 14)
print(g2)
ggsave("./figures/present/fecundity_size_slope.png", width = 5, height = 4)


## R2
r2 <- as.data.frame(bayes_R2(fit_fs2, summary = FALSE))
g3 <- ggplot(r2) +
    # geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
    geom_histogram(aes(x = R2), bins = 20, color = "white",
                   fill = "grey30", alpha = 0.5, na.rm = TRUE) +
    # scale_x_continuous(limits = c(-0.3, 14),
    #                    breaks = seq(0, 14, by = 2)) +
    labs(x = bquote(R^2), y = "Posterior density") +
    theme_simple(base_size = 14)
print(g3)
ggsave("./figures/present/fecundity_size_r2.png", g3,
       width = 5, height = 4)


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
    geom_point(aes(x = length_avg, y = fecundity), color = "grey70", size = 1,
               na.rm = TRUE) +
    # geom_path(data = df, aes(x = x, y = y, group = id_label), color = "grey70",
    #           size = 0.7) +
    geom_path(data = pr, aes(x = x, y = y), color = "red3", size = 1.2) +
    geom_ribbon(data = pr, aes(x = x, ymin = y_lower, ymax = y_upper),
                fill = "grey30", color = NA, alpha = 0.2) +
    labs(x = "Length (mm)", y = "Fecundity") +
    theme_simple(base_size = 14)
print(g4)
ggsave("./figures/present/fecundity_size_scatter.png", width = 5, height = 4)


