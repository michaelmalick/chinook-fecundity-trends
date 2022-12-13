## Publication map of hatcheries

dir.create("./figures/pub/", showWarnings = FALSE)
dir.create("./figures/present/", showWarnings = FALSE)


## Data prep -----------------------------------------------

## Set x and y limits for map
xlim <- c(-125, -117)
ylim <- c(45, 49)


## Load map data
bmap        <- readRDS("./data/maps/gshhg_pnw_map.rds")
bmap_nation <- readRDS("./data/maps/gshhg_pnw_border_national.rds")
bmap_state  <- readRDS("./data/maps/gshhg_pnw_border_state.rds")
nhd         <- readRDS("./data/maps/nhds17_raw.rds")
ras         <- raster("./data/maps/srtm_5arcsecond_clip_imhof4.tif")
world       <- ne_states(returnclass = "sf",
                         country = c("united states of america", "canada", "mexico"))


## Subsample raster for faster plotting
ras_sub <- sampleRegular(ras, 1e6, asRaster = TRUE)
ras_sub_df <- as.data.frame(ras_sub, xy = TRUE)
names(ras_sub_df) <- c("x", "y", "z")


## Get flow lines
nhd_flow <- lapply(nhd, function(x) {
    k <- x$NHDFlowline
    m <- k[k$StreamOrde > 2, ]
    sf::st_transform(m, proj_wgs)
})


## Get HUC polygons
nhd_huc4 <- lapply(nhd, function(x) {
    k <- x$WBDHU4
    m <- sf::st_geometry(k)
    sf::st_transform(m, proj_wgs)
})


## Get Columbia + Snake Rivers
nhd_cr <- lapply(nhd_flow, function(x) {
    k <- x[x$GNIS_Name %in% c("Columbia River", "Snake River"), ]
    m <- sf::st_geometry(k)
    sf::st_transform(m, proj_wgs)
})


## Clip HUC polgons to xy lims
h02 <- sf::st_crop(nhd_huc4[["hu02"]],
    xmin = xlim[1], xmax = xlim[2], ymin = ylim[1], ymax = ylim[2])
h03 <- sf::st_crop(nhd_huc4[["hu03"]],
    xmin = xlim[1], xmax = xlim[2], ymin = ylim[1], ymax = ylim[2])
h06 <- sf::st_crop(nhd_huc4[["hu06"]],
    xmin = xlim[1], xmax = xlim[2], ymin = ylim[1], ymax = ylim[2])
h07 <- sf::st_crop(nhd_huc4[["hu07"]],
    xmin = xlim[1], xmax = xlim[2], ymin = ylim[1], ymax = ylim[2])
h08 <- sf::st_crop(nhd_huc4[["hu08"]],
    xmin = xlim[1], xmax = xlim[2], ymin = ylim[1], ymax = ylim[2])
h09 <- sf::st_crop(nhd_huc4[["hu09"]],
    xmin = xlim[1], xmax = xlim[2], ymin = ylim[1], ymax = ylim[2])
h10 <- sf::st_crop(nhd_huc4[["hu10"]],
    xmin = xlim[1], xmax = xlim[2], ymin = ylim[1], ymax = ylim[2])
h11 <- sf::st_crop(nhd_huc4[["hu11"]],
    xmin = xlim[1], xmax = xlim[2], ymin = ylim[1], ymax = ylim[2])



## Hatchery data
ff <- copy(fecund_summary)
ff$source <- ifelse(grepl("NFH ", ff$id_label), "USFWS", "WDFW")
ff$source <- factor(ff$source, levels = c("WDFW", "USFWS"))
for(i in 1:nrow(ff)) {
    d <- ff[i, ]
    if(d$duplicated) {
        ind <- which(ff$hatchery == d$hatchery)
        ff$number[ind[1]] <- paste0(ff$number[ind[1]], ",", d$number)
    }
}
ff <- ff[duplicated == FALSE, ]
ff$nudge_x <- 0.06
ff$nudge_y <- 0.06
ff$nudge_y[ff$number == "26,32"] <- -0.05
ff$nudge_x[ff$number == "26,32"] <- -0.025
ff$nudge_y[ff$number == "25,31"] <- -0.05
ff$nudge_y[ff$number == "40"] <- -0.05
ff$nudge_y[ff$number == "7"]  <- -0.05
ff$nudge_x[ff$number == "39"]  <- -0.09
ff$nudge_y[ff$number == "39"]  <- 0.03

size_id <- as.character(unique(fecund_size$id_label))
ff$size <- ifelse(ff$id_label %in% size_id, "yes", "no")



## Create map ----------------------------------------------
g <- ggplot(ras_sub_df) +
    geom_sf(data = bmap, color = "white", size = 0.3, fill = NA) +
    geom_raster(aes(x = x, y = y, fill = z), na.rm = TRUE) +
    geom_sf(data = bmap_nation, color = "grey50", size = 0.3, fill = NA) +
    ## Columbia River
    geom_sf(data = nhd_cr[[1]], color = "steelblue", size = 0.3, fill = NA) +
    geom_sf(data = nhd_cr[[2]], color = "steelblue", size = 0.3, fill = NA) +
    geom_sf(data = nhd_cr[[3]], color = "steelblue", size = 0.3, fill = NA) +
    geom_sf(data = nhd_cr[[4]], color = "steelblue", size = 0.3, fill = NA) +
    geom_sf(data = nhd_cr[[5]], color = "steelblue", size = 0.3, fill = NA) +
    geom_sf(data = nhd_cr[[6]], color = "steelblue", size = 0.3, fill = NA) +
    geom_sf(data = nhd_cr[[7]], color = "steelblue", size = 0.3, fill = NA) +
    geom_sf(data = nhd_cr[[8]], color = "steelblue", size = 0.3, fill = NA) +
    ## Coast
    geom_sf(data = h10, color = NA, fill = M1[1], alpha = 0.15) +
    ## Puget Sound
    geom_sf(data = h11, color = NA, fill = M1[2], alpha = 0.15) +
    ## Lower Columbia
    geom_sf(data = h07, color = NA, fill = M1[3], alpha = 0.15) +
    geom_sf(data = h08, color = NA, fill = M1[3], alpha = 0.15) +
    geom_sf(data = h09, color = NA, fill = M1[3], alpha = 0.15) +
    ## Upper Columbia
    geom_sf(data = h02, color = NA, fill = M1[4], alpha = 0.15) +
    geom_sf(data = h03, color = NA, fill = M1[4], alpha = 0.15) +
    geom_sf(data = h06, color = NA, fill = M1[4], alpha = 0.15) + ## Snake
    #
    geom_point(data = ff,
               aes(x = lon, y = lat, color = region, shape = source), size = 1.4) +
    geom_point(data = ff[size == "yes" & source == "WDFW", ],
               aes(x = lon, y = lat), color = "black", shape = 1, size = 1.4) +
    geom_point(data = ff[size == "yes" & source == "USFWS", ],
               aes(x = lon, y = lat), color = "black", shape = 2, size = 1.4) +
    geom_text(data = ff, aes(x = lon, y = lat, label = number),
              size = 2.4, nudge_x = ff$nudge_x, nudge_y = ff$nudge_y) +
#
    annotate("text", label = "Coast", x = -123.4, y = 46.8,
             size = 5, colour = M1d[1]) +
    annotate("text", label = "Puget Sound", x = -122.1, y = 48.2,
             size = 5, colour = M1d[2]) +
    annotate("text", label = "Lower Columbia", x = -122, y = 45.4,
             size = 5, colour = M1d[3]) +
    annotate("text", label = "Upper Columbia", x = -120.1, y = 47.3,
             size = 5, colour = M1d[4]) +
    #
    scale_x_continuous(limits = c(-125, -117), expand = c(0, 0),
                       breaks = -124:-118) +
    scale_y_continuous(limits = c(45.1, 49.1), expand = c(0, 0),
                       breaks = 46:49) +
    labs(x = NULL, y = NULL) +
    scale_color_manual(values = M1d) +
    # scale_shape_manual(values = c(21, 24)) +
    scale_fill_gradient(low = "grey60",
                        high = "grey98",
                        na.value = "white") +
    theme_simple() +
    theme(legend.position = "none")
print(g)
# ggsave("./figures/pub/map.png", width = 6.9, height = 5, units = "in")

## Inset map
gi <- ggplot() +
    geom_sf(data = world, fill = "grey40", color = "grey40", size = 0.1) +
    geom_rect(data = mtcars[1,],
              aes(xmin = -125, xmax = -117, ymin = 45, ymax = 49.1),
              color = "black", fill = NA, alpha = 1, size = 0.5,
              inherit.aes = FALSE, linetype = 1) +
    ylim(16, 75) +
    xlim(-170, -60) +
    theme_simple() +
    theme(# panel.border = element_blank(),
          panel.border = element_rect(fill = NA, colour = "grey50", size = rel(0.5)),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent", color = NA),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
# print(gi)

## Inset distributions
gr <-fecund[ , .N, by = .(group, region)]
gr$fecundity <- 1800
gd <- ggplot(fecund[fecundity > 2000, ]) +
    geom_vline(xintercept = mean(fecund$fecundity, na.rm = TRUE),
               color = "grey50", linetype = 2) +
    aes(x = fecundity, y = group, fill = region, color = region) +
    # geom_jitter(height = 0.2, na.rm = TRUE, size = 0.3) +
    geom_violin(na.rm = TRUE, alpha = 0.5, adjust = 1.3, size = 0.4) +
    geom_text(data = gr, aes(x = fecundity, y = group, label = N), size = 1.5) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 1,
                 color = "black", na.rm = TRUE) +
    labs(x = "Fecundity",
         y = NULL) +
    xlim(1800, 7100) +
    scale_color_manual(values = M1) +
    scale_fill_manual(values = M1) +
    theme_simple(base_size = 7) +
    theme(legend.position = "none") +
    theme(plot.background = element_rect(fill = alpha("white", 0.85), colour = "grey50",
                                         size = rel(0.5)))
print(gd)


png("./figures/pub/map.png", width = 6.9, height = 5, units = "in", res = 400)
    vp_map <- viewport(width = 0.40, height = 0.30, x = 0.85, y = 0.8)
    vp_dist <- viewport(width = 0.35, height = 0.42, x = 0.80, y = 0.28)
    print(g)
    print(gi, vp = vp_map)
    print(gd, vp = vp_dist)
dev.off()

pdf("./figures/pub/map.pdf", width = 6.9, height = 5)
    vp_map <- viewport(width = 0.40, height = 0.30, x = 0.85, y = 0.8)
    vp_dist <- viewport(width = 0.35, height = 0.42, x = 0.80, y = 0.28)
    print(g)
    print(gi, vp = vp_map)
    print(gd, vp = vp_dist)
dev.off()



## Presentation map ----------------------------------------

g <- ggplot(ras_sub_df) +
    geom_sf(data = bmap, color = "white", size = 0.3, fill = NA) +
    geom_raster(aes(x = x, y = y, fill = z), na.rm = TRUE) +
    geom_sf(data = bmap_nation, color = "grey50", size = 0.3, fill = NA) +
    ## Columbia River
    geom_sf(data = nhd_cr[[1]], color = "steelblue", size = 0.3, fill = NA) +
    geom_sf(data = nhd_cr[[2]], color = "steelblue", size = 0.3, fill = NA) +
    geom_sf(data = nhd_cr[[3]], color = "steelblue", size = 0.3, fill = NA) +
    geom_sf(data = nhd_cr[[4]], color = "steelblue", size = 0.3, fill = NA) +
    geom_sf(data = nhd_cr[[5]], color = "steelblue", size = 0.3, fill = NA) +
    geom_sf(data = nhd_cr[[6]], color = "steelblue", size = 0.3, fill = NA) +
    geom_sf(data = nhd_cr[[7]], color = "steelblue", size = 0.3, fill = NA) +
    geom_sf(data = nhd_cr[[8]], color = "steelblue", size = 0.3, fill = NA) +
    ## Coast
    geom_sf(data = h10, color = NA, fill = M1[1], alpha = 0.15) +
    ## Puget Sound
    geom_sf(data = h11, color = NA, fill = M1[2], alpha = 0.15) +
    ## Lower Columbia
    geom_sf(data = h07, color = NA, fill = M1[3], alpha = 0.15) +
    geom_sf(data = h08, color = NA, fill = M1[3], alpha = 0.15) +
    geom_sf(data = h09, color = NA, fill = M1[3], alpha = 0.15) +
    ## Upper Columbia
    geom_sf(data = h02, color = NA, fill = M1[4], alpha = 0.15) +
    geom_sf(data = h03, color = NA, fill = M1[4], alpha = 0.15) +
    geom_sf(data = h06, color = NA, fill = M1[4], alpha = 0.15) + ## Snake
    #
    geom_point(data = ff,
               aes(x = lon, y = lat, color = region, shape = size), size = 2) +
    # geom_text(data = ff, aes(x = lon, y = lat, label = number),
    #           size = 2.4, nudge_x = ff$nudge_x, nudge_y = ff$nudge_y) +
    #
    annotate("text", label = "Coast", x = -123.4, y = 46.8,
             size = 5, colour = M1d[1]) +
    annotate("text", label = "Puget Sound", x = -122.1, y = 48.2,
             size = 5, colour = M1d[2]) +
    annotate("text", label = "Lower Columbia", x = -122, y = 45.4,
             size = 5, colour = M1d[3]) +
    annotate("text", label = "Upper Columbia", x = -120.1, y = 47.3,
             size = 5, colour = M1d[4]) +
    #
    scale_x_continuous(limits = c(-125, -117), expand = c(0, 0),
                       breaks = -124:-118) +
    scale_y_continuous(limits = c(45.1, 49.1), expand = c(0, 0),
                       breaks = 46:49) +
    labs(x = NULL, y = NULL) +
    scale_color_manual(values = M1d) +
    scale_fill_gradient(low = "grey60",
                        high = "grey98",
                        na.value = "white") +
    theme_simple() +
    theme(legend.position = "none",
          panel.border=element_blank(),
          axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank())
print(g)
ggsave("./figures/present/map.png", width = 6.9, height = 5, units = "in")
