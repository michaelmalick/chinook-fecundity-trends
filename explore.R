## Exploratory analysis

dir.create("./figures/explore/", showWarnings = FALSE)

## Fecundity: Map facilities -------------------------------
world_sf <- ne_countries(scale = "large", returnclass = "sf")

d <- fecund[ , .(lon = unique(lon),
                          lat = unique(lat),
                          stock = unique(hatchery),
                          run = unique(run),
                          region = unique(region),
                          subregion = unique(subregion)), by = "id_label"]

g <- ggplot(d) +
    geom_sf(data = world_sf) +
    geom_point(aes(x = lon, y = lat, color = region), size = 2) +
    xlim(-126, -118) +
    ylim(45, 50) +
    labs(x = "Longitude",
         y = "Latitude") +
    scale_color_manual(values = M1) +
    facet_wrap( ~ run, nrow = 2) +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/explore/map_regions.png", units = "in",
       width = 8, height = 6)


## Fecundity: Plot raw data --------------------------------

## Fecundity by stock
d <- fecund[ , .(mean = mean(fecundity, na.rm = TRUE)), by = "id_label"]
g <- ggplot(fecund) +
    geom_hline(data = d, aes(yintercept = mean), color = "grey50", linetype = 2) +
    aes(x = year, y = fecundity, group = id_label, color = run) +
    geom_line(na.rm = TRUE) +
    labs(y = "Fecundity") +
    scale_color_manual(values = M1) +
    theme_simple() +
    facet_wrap( ~ id_label)
print(g)
ggsave("./figures/explore/fecundity_raw_stock.png", units = "in",
       width = 18, height = 11)

## Fecundity by group
d <- fecund[ , .(mean = mean(fecundity, na.rm = TRUE)), by = "group"]
g <- ggplot(fecund) +
    geom_hline(data = d, aes(yintercept = mean), color = "grey50", linetype = 2) +
    aes(x = year, y = fecundity, group = id_label, color = region) +
    geom_line(na.rm = TRUE) +
    labs(y = "Fecundity") +
    scale_color_manual(values = M1) +
    theme_simple() +
    theme(legend.position = "none") +
    facet_wrap( ~ group, ncol = 2)
print(g)
ggsave("./figures/explore/fecundity_raw_group.png", units = "in",
       width = 6, height = 6)


## Fecundity by run
m <- mean(fecund$fecundity, na.rm = TRUE)
k <- fecund[ , .(mean = mean(fecundity, na.rm = TRUE)) , by = .(year, run)]
g <- ggplot(fecund) +
    geom_hline(yintercept = m, color = "grey50", linetype = 2) +
    aes(x = year, y = fecundity, group = id_label, color = run) +
    geom_line(na.rm = TRUE) +
    geom_line(data = k, lwd = 2, aes(x = year, y = mean,
                                     color = run, group = run)) +
    labs(y = "Fecundity") +
    scale_color_manual(values = M1) +
    theme_simple()
print(g)
ggsave("./figures/explore/fecundity_raw_run.png", units = "in",
       width = 6, height = 4)


## Histogram
g <- ggplot(fecund) +
    aes(x = fecundity) +
    geom_histogram(bins = 30, color = "white", na.rm = TRUE) +
    labs(x = "Fecundity",
         y = "Count") +
    theme_simple()
print(g)
ggsave("./figures/explore/fecundity_raw_hist.png", units = "in",
       width = 6, height = 4)


## Violin by group
g <- ggplot(fecund) +
    geom_vline(xintercept = mean(fecund$fecundity, na.rm = TRUE),
               color = "grey50", linetype = 2) +
    aes(x = fecundity, y = group, fill = region, color = region) +
    geom_jitter(height = 0.2, na.rm = TRUE, size = 0.5) +
    geom_violin(na.rm = TRUE, alpha = 0.5, adjust = 1.3) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 2,
                 color = "black", na.rm = TRUE) +
    labs(x = "Fecundity",
         y = "") +
    scale_color_manual(values = M1) +
    scale_fill_manual(values = M1) +
    theme_simple() +
    theme(legend.position = "none")
print(g)
ggsave("./figures/explore/fecundity_raw_violin_group.png", units = "in",
       width = 5, height = 6)



## Fecundity: Plot standardized data -----------------------

## Fecundity by group
g <- ggplot(fecund) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    aes(x = year, y = fecundity_stnd, group = id_label, color = region) +
    geom_line(na.rm = TRUE) +
    labs(y = "Standardized fecundity") +
    scale_color_manual(values = M1) +
    theme_simple() +
    theme(legend.position = "none") +
    facet_wrap( ~ group, ncol = 2)
print(g)
ggsave("./figures/explore/fecundity_stnd_group.png", units = "in",
       width = 6, height = 6)


## Fecundity by group GAM
g <- ggplot(fecund) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    aes(x = year, y = fecundity_stnd, color = region) +
    geom_point(na.rm = TRUE, size = 0.7) +
    geom_smooth(na.rm = TRUE, method = "gam", formula = y ~ s(x)) +
    labs(y = "Standardized fecundity") +
    scale_color_manual(values = M1) +
    theme_simple() +
    theme(legend.position = "none") +
    facet_wrap( ~ group, ncol = 2)
print(g)
ggsave("./figures/explore/fecundity_stnd_group_gam.png", units = "in",
       width = 6, height = 6)


## Fecundity all stocks
g <- ggplot(fecund) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    aes(x = year, y = fecundity_stnd, group = id_label) +
    geom_line(color = "grey70", na.rm = TRUE) +
    stat_summary(fun = mean, geom = "line", lwd = 1.2,
                 aes(group = 1), na.rm = TRUE) +
    labs(y = "Standardized fecundity") +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/explore/fecundity_stnd_all.png", units = "in",
       width = 4, height = 3)


## Fecundity by run
g <- ggplot(fecund) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    aes(x = year, y = fecundity_stnd, group = id_label, color = run) +
    geom_line(na.rm = TRUE, alpha = 0.5) +
    stat_summary(fun = mean, geom = "line", lwd = 1.5,
                 aes(group = 1), na.rm = TRUE) +
    labs(y = "Standardized fecundity") +
    scale_color_manual(values = M1) +
    theme_simple() +
    theme(legend.position = "none") +
    facet_wrap( ~ run, ncol = 1)
print(g)
ggsave("./figures/explore/fecundity_stnd_run.png", units = "in",
       width = 6, height = 8)


## Fecundity by region
g <- ggplot(fecund) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    aes(x = year, y = fecundity_stnd, group = id_label, color = region) +
    geom_line(na.rm = TRUE, alpha = 0.5) +
    stat_summary(fun = mean, geom = "line", lwd = 1.5,
                 aes(group = 1), na.rm = TRUE) +
    labs(y = "Standardized fecundity") +
    scale_color_manual(values = M1) +
    theme_simple() +
    theme(legend.position = "none") +
    facet_wrap( ~ region, ncol = 1)
print(g)
ggsave("./figures/explore/fecundity_stnd_region.png", units = "in",
       width = 6, height = 8)


## Fecundity fall variance
d <- fecund[run == "Fall", ]
v <- d[ , .(var = sd(fecundity, na.rm = TRUE)), by = .(year)]
perc <- ecdf(v$var)
v$perc <- perc(v$var)
v$col <- ifelse(v$perc <= 0.25, "blue3", "grey40")
v$lab <- ifelse(v$perc <= 0.25, v$year, "")
v$col <- ifelse(v$perc >= 0.75, "red3", v$col)
v$lab <- ifelse(v$perc >= 0.75, v$year, v$lab)
g <- ggplot(v) +
    geom_hline(yintercept = mean(v$var), color = "grey50", linetype = 2) +
    aes(x = year, y = var) +
    geom_line(color = "grey40", na.rm = TRUE) +
    geom_text(aes(label = lab), nudge_x = 0.75, size = 2.5) +
    geom_point(color = v$col, na.rm = TRUE) +
    labs(y = "Standard deviation") +
    theme_simple()
print(g)
ggsave("./figures/explore/fecundity_fall_var.png", units = "in",
       width = 6, height = 4)


## Fecundity: Multiple run timings -------------------------
x <- unique(fecund[ , c("hatchery", "run")])
dup <- as.character(x$hatchery[duplicated(x$hatchery)])
fecund_run <- fecund[hatchery %in% dup, ]

g <- ggplot(fecund_run) +
    aes(x = year, y = fecundity, color = run) +
    geom_point(na.rm = TRUE) +
    geom_line(na.rm = TRUE) +
    labs(x = "Year", y = "Fecundity", color = "Run") +
    scale_color_manual(values = M1) +
    facet_wrap( ~ hatchery) +
    theme_simple()
print(g)

ggsave("./figures/explore/fecundity_run_timing.png", units = "in",
       width = 8, height = 6)
## Size ----------------------------------------------------

## size ~ year: time series
g <- ggplot(age_size) +
    aes(x = year, y = length, color = as.factor(age)) +
    geom_point(na.rm = TRUE) +
    geom_line(na.rm = TRUE) +
    scale_colour_manual(values = M1) +
    labs(x = "Year", y = "Length",
         title = "Female Chinook",
         color = "Age") +
    facet_wrap(~ id_label) +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/explore/size_trends_line.png",
       units = "in", width = 11, height = 6)


## size ~ year: smooth
g <- ggplot(age_size) +
    aes(x = year, y = length, color = as.factor(age)) +
    geom_point(na.rm = TRUE) +
    geom_smooth(method = "lm", formula = y ~ x, na.rm = TRUE, se = FALSE) +
    scale_colour_manual(values = M1) +
    labs(x = "Year", y = "Length",
         title = "Female Chinook",
         color = "Age") +
    facet_wrap(~ id_label) +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/explore/size_trends_smooth.png",
       units = "in", width = 11, height = 6)


## size ~ year: slopes
slope <- age_size[ , {
    f <- lm(length ~ year)
    .(slope = coef(f)[2])}, by = .(id_label, group, age)]
g <- ggplot(slope) +
    geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
    geom_point(aes(x = slope, y = id_label, color = as.factor(age)), size = 3) +
    scale_colour_manual(values = M1) +
    labs(x = "Slope (mm / year)", y = "Stock",
         title = "Slope: length ~ year",
         color = "Age") +
    facet_wrap( ~ age, nrow = 1) +
    theme_simple() +
    theme(legend.position="none")
print(g)
ggsave("./figures/explore/size_trends_slopes.png",
       units = "in", width = 11, height = 6)


## size combined
g <- ggplot(age_size) +
    aes(x = year, y = length, color = as.factor(age)) +
    geom_point(na.rm = TRUE) +
    geom_smooth(method = "lm", formula = y ~ x, na.rm = TRUE, se = TRUE) +
    scale_colour_manual(values = M1) +
    labs(x = "Year", y = "Length",
         title = "Female Chinook",
         color = "Age") +
    # facet_wrap(~ age, scale = "free_y") +
    facet_wrap(~ age) +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/explore/size_trends_combined.png",
       units = "in", width = 8, height = 6)


## Age proportion ------------------------------------------

## age prop ~ year: time series
g <- ggplot(age_size) +
    aes(x = year, y = age_prop, color = as.factor(age)) +
    geom_point(na.rm = TRUE) +
    geom_line(na.rm = TRUE) +
    scale_colour_manual(values = M1) +
    labs(x = "Year", y = "Proportion",
         title = "Female Chinook",
         color = "Age") +
    facet_wrap(~ id_label) +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/explore/ageprop_trends_line.png",
       units = "in", width = 12, height = 6)

## age prop ~ year: smooth
g <- ggplot(age_size) +
    aes(x = year, y = age_prop, color = as.factor(age)) +
    geom_point(na.rm = TRUE) +
    geom_smooth(method = "lm", formula = y ~ x, na.rm = TRUE, se = FALSE) +
    scale_colour_manual(values = M1) +
    labs(x = "Year", y = "Proportion",
         title = "Female Chinook",
         color = "Age") +
    facet_wrap(~ id_label) +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/explore/ageprop_trends_smooth.png",
       units = "in", width = 11, height = 6)

## age prop ~ year: slopes
slope <- age_size[ , {
    f <- lm(age_prop ~ year)
    .(slope = coef(f)[2])}, by = .(id_label, group, age)]
g <- ggplot(slope) +
    geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
    geom_point(aes(x = slope, y = id_label, color = as.factor(age)), size = 3) +
    scale_colour_manual(values = M1) +
    labs(x = "Slope", y = "Stock",
         title = "Slope: age-prop ~ year",
         color = "Age") +
    facet_wrap( ~ age, nrow = 1) +
    theme_simple() +
    theme(legend.position="none")
print(g)
ggsave("./figures/explore/ageprop_trends_slopes.png",
       units = "in", width = 11, height = 6)



## Fecundity + Size ----------------------------------------

## fecundity ~ size: all
g <- ggplot(fecund_size) +
    aes(x = length_avg, y = fecundity, color = id_label) +
    geom_point(na.rm = TRUE, alpha = 0.3) +
    geom_smooth(method = "lm", formula = y ~ x, na.rm = TRUE, se = FALSE) +
    labs(x = "Length", y = "Fecundity", color = "Stock") +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/explore/fecundity_vs_size_all.png",
       units = "in", width = 8, height = 6)

g <- ggplot(fecund_size) +
    aes(x = length_avg_stnd, y = fecundity_stnd, color = id_label) +
    geom_point(na.rm = TRUE, alpha = 0.3) +
    geom_smooth(method = "lm", formula = y ~ x, na.rm = TRUE, se = FALSE) +
    labs(x = "Length", y = "Fecundity", color = "Stock") +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/explore/fecundity_vs_size_all_stnd.png",
       units = "in", width = 8, height = 6)

g <- ggplot(fecund_size) +
    aes(x = length_avg_anom, y = fecundity_anom, color = id_label) +
    geom_point(na.rm = TRUE, alpha = 0.3) +
    geom_smooth(method = "lm", formula = y ~ x, na.rm = TRUE, se = FALSE) +
    labs(x = "Length", y = "Fecundity", color = "Stock") +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/explore/fecundity_vs_size_all_anom.png",
       units = "in", width = 8, height = 6)

## fecundity ~ size: release
g <- ggplot(fecund_size) +
    aes(x = length_avg, y = fecundity, color = release) +
    geom_point(na.rm = TRUE, alpha = 0.3) +
    geom_smooth(method = "lm", formula = y ~ x, na.rm = TRUE, se = FALSE) +
    labs(x = "Length", y = "Fecundity", color = "Stock") +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/explore/fecundity_vs_size_release.png",
       units = "in", width = 8, height = 6)

## fecundity ~ size: group
g <- ggplot(fecund_size) +
    aes(x = length_avg, y = fecundity, color = group) +
    geom_point(na.rm = TRUE, alpha = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x, na.rm = TRUE, se = FALSE) +
    labs(x = "Length", y = "Fecundity") +
    facet_wrap( ~ group) +
    theme_simple(grid = TRUE) +
    theme(legend.position="none")
print(g)
ggsave("./figures/explore/fecundity_vs_size_group.png",
       units = "in", width = 8, height = 6)

## fecundity ~ size: stock
g <- ggplot(fecund_size) +
    aes(x = length_avg, y = fecundity, color = group) +
    geom_point(na.rm = TRUE, alpha = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x, na.rm = TRUE, se = FALSE) +
    labs(x = "Length", y = "Fecundity") +
    facet_wrap( ~ id_label) +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/explore/fecundity_vs_size_stock.png",
       units = "in", width = 9, height = 6)


