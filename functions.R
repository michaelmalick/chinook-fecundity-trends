## Functions for the analysis


## yrep_plot -----------------------------------------------
yrep_plot <- function(y, yrep, nsamples = 100, main = "Yrep") {
    ## Posterior predictive checks
    ## Density overlay of y and yrep
    ##
    ## y = vector of y values
    ## yrep = matrix of yrep datasets (rows = mcmc draws, columns = y-datapoints)
    ## nsamples = number of yrep samples to plot
    ## main = title for plot
    ysub <- yrep[sample(1:nrow(yrep), nsamples), ]
    ysub <- rbind(ysub, y)
    ysub[nsamples + 1, ] <- y

    dens <- apply(ysub, 1, stats::density, adjust = 1.5, na.rm = TRUE)
    dens.x <- lapply(dens, function(x) x$x)
    dens.x <- do.call("cbind", dens.x)
    dens.y <- lapply(dens, function(x) x$y)
    dens.y <- do.call("cbind", dens.y)

    graphics::matplot(dens.x, dens.y,
                      main = main,
                      type = "l",
                      lty  = 1,
                      col  = "grey50",
                      ylab = "Density",
                      xlab = "Value",
                      axes = FALSE)
    graphics::box(col = "grey50")
    graphics::axis(1, lwd = 0, col = "grey50", lwd.ticks = 1)
    graphics::axis(2, lwd = 0, col = "grey50", las = 1, lwd.ticks = 1)
    graphics::lines(stats::density(y, adjust = 1.5, na.rm = TRUE), lwd = 2)
}


## ss_yrep -------------------------------------------------
ss_yrep <- function(mcmc) {
    ## Get yrep matrix from ss_rw.jags model
    ## mcmc = an mcmc.list object
    mat <- as.matrix(mcmc)
    X <- mat[ , grep("X", coda::varnames(mcmc), value = TRUE)]
    xorder <- paste0("X[", 1:ncol(X), "]")
    X <- X[ , xorder]
    sigma_obs <- mat[ , "sigma_obs", drop = FALSE]

    yrep <- matrix(NA, ncol = ncol(X), nrow = nrow(X))
    for(i in 1:nrow(X)) {
        for(j in 1:ncol(X)) {
            yrep[i,j] <- rnorm(1, X[i, j], sigma_obs[i, 1])
        }
    }
    colnames(yrep) <- paste0("yrep[", 1:ncol(X), "]")
    return(yrep)
}


## cor_plot_region -----------------------------------------
cor_plot_region <- function(dcor, max_size = 5,
                            title = NULL,
                            legend_title = "r",
                            group_labels = TRUE) {
    ## Create a grouped based correlation matrix plot
    ##
    ## dcor = data.frame ouput from dd_cor_group()
    ## max_size = max size of bubbles
    ## group_labels = logical, should grouped text labels be used on y-axis

    dcor <- as.data.table(dcor)
    setorder(dcor, "group1")
    dcor[ , color := .GRP, by = group1]
    labs <- dcor[ , .(group1 = unique(group1), color = unique(color)),
                 by = "id_label1"]
    # lines to demarcate regions in fig
    labs[ , diff := c(0, diff(color))]
    labs[ , gr_line := ifelse(diff == 1, .I - 0.5, NA) , by = .I]

    if(group_labels) {
        labs[ , c("label1", "label2") := {
            med <- median(seq(1, .N))
            a <- rep("", .N)
            b <- rep("", .N)
            a[med] <- paste(unique(as.character(group1)), " ", unique(color))
            b[med] <- as.character(unique(color))
            list(a, b)
        }, by = .(group1)]
    } else {
        labs[ , c("label1", "label2") := {
            med <- median(seq(1, .N))
            a <- rep("", .N)
            b <- rep("", .N)
            a[med] <- as.character(unique(color))
            b[med] <- as.character(unique(color))
            list(a, b)
        }, by = .(group1)]
    }

    g <- ggplot(dcor) +
        geom_point(aes(x = id_label1, y = id_label2, color = cor, size = cor),
                   na.rm = TRUE) +
        geom_hline(data = labs, aes(yintercept = gr_line),
                   color = "grey70", na.rm = TRUE, size = 0.3) +
        geom_vline(data = labs, aes(xintercept = gr_line),
                   color = "grey70", na.rm = TRUE, size = 0.3) +
        scale_size_area(max_size = max_size, guide = "none") +
        scale_color_gradient2(low = "steelblue",
                             mid = "grey90",
                             high = "red3", midpoint = 0,
                             space = "Lab",
                             limits = c(-1, 1),
                             na.value = NA) +
        labs(x = NULL, y = NULL, color = legend_title, title = title, size = legend_title) +
        scale_y_discrete(labels = labs$label1) +
        scale_x_discrete(labels = labs$label2) +
        theme_simple(grid = FALSE)
    print(g)
    return(g)
}



## dd_cor_group --------------------------------------------
dd_cor_group <- function(data,
                         id.var,
                         grp.var,
                         cor.var,
                         time.var,
                         min.overlap, ...) {
    ## Run pairwise correlations and add a grouping variable
    ## see dd_cor for arg help

    dcor <- dd_cor(data,
                   id.var = id.var,
                   cor.var = cor.var,
                   time.var = time.var,
                   min.overlap = min.overlap, ...)
    dcor <- as.data.table(dcor)

    dcor[ , group1 := NA]
    dcor[ , group2 := NA]
    grp <- as.character(data[[grp.var]])
    id  <- as.character(data[[id.var]])
    for(i in 1:nrow(dcor)) {
        gr1 <- unique(grp[id == dcor$id_label1[i]])
        gr2 <- unique(grp[id == dcor$id_label2[i]])
        dcor$group1[i] <- as.character(gr1)
        dcor$group2[i] <- as.character(gr2)
    }

    grp_lev <- levels(data[[grp.var]])
    lev <- grp_lev[grp_lev %in% unique(dcor$group1)]
    dcor[ , group1 := factor(group1, levels = lev)]
    dcor[ , group2 := factor(group2, levels = lev)]
    return(dcor)
}



## dd_cor --------------------------------------------------
# This function takes as input a dataframe giving the variable to be
# correlated, an identifyer for the variables (usually a stock name
# or id number), and a time variable. It breaks apart the data
# according to the levels of specified id variable and computes
# pearson correlation coefficients among all combinations of the
# levels of the factor.
#
# The function outputs a dataframe and heatmap of the cross-
# correlations of all unique combinations between the variables.
#
# There is an option (min.overlap) to specify how many years of
# overlap are required between two data series in order to compute a
# correlation. If the there is not enough overlap an NA is returned
# instead.

dd_cor <- function(
    data,
    id.var,
    cor.var,
    time.var,
    min.overlap     = 10,
    sort            = FALSE,
    only.upper      = TRUE,
    no.diagonal     = TRUE,
    plot            = FALSE,
    plot.title      = "",
    return.matrix   = FALSE,
    return.n.common = FALSE) {

    # data = dataframe in 'long' format
    # id.var = column name in data giving unique identifier for each
    #          variable to correlate
    # cor.var = column name in data identifying variable to correlate
    # time.var = column name in data giving time variable
    # min.overlap = minimum number of years between to data series
    #               required to compute a correlation
    # sort = should the unique identifiers be sorted
    # only.upper = should only upper triangle of the correlation
    #              matrix be returned and plotted
    # no.diagonal = if only upper is true should the diagonal be
    #               included
    # plot = should a plot be created
    # plot.title = title for the resulting plot
    # return.matrix = logical, if FALSE a melted matrix is returned, if TRUE
    #                 a correlation matrix is returned
    # return.n.common = if TRUE a column for giving the number of years of
    #                   overlap is added to the output

    data <- as.data.frame(data)

    # Subset dataset and reshape data
    df.sub <- subset(data, select = c(id.var, time.var, cor.var))

    df.cast <- reshape2::dcast(df.sub, get(time.var) ~ get(id.var),
        value.var = cor.var)

    # Find correlation combinations
    if(sort)
        nn <- sort(unique(data[ , id.var]))
    else
        nn <- unique(data[ , id.var])


    xx <- expand.grid(nn, nn)


    # Create matrix to store correlations
    cor.mat <- matrix(NA, ncol = length(nn),
        nrow = length(nn), dimnames = list(nn, nn))

    # Create vector to store distances
    n.common <- rep(NA, length(nn))

    # Compute correlations looping over correlation combinations
    for(i in 1:length(xx[ , 1])) {

        # Create indices for correlations
        ind.x <- as.character(xx[i, 1])
        ind.y <- as.character(xx[i, 2])

        # Subset timeseries to correlate
        x <- df.cast[ , ind.x]
        y <- df.cast[ , ind.y]

        # Determine no. of overlapping years
        # If less than min.overlap set it to 0 otherwise
        # compute and save the correlation
        test.x <- !is.na(x)
        test.y <- !is.na(y)
        test.x[test.x == FALSE] <- NA
        test.y[test.y == FALSE] <- NA
        n.common[i] <- sum(test.x == test.y, na.rm = TRUE)

        if(n.common[i] >= min.overlap)
            cor.mat[ind.x, ind.y] <- cor(x, y,
                use = "pairwise.complete.obs")

        if(n.common[i] < min.overlap)
            cor.mat[ind.x, ind.y] <- NA

    }


    if(only.upper) cor.mat[lower.tri(cor.mat, no.diagonal)] <- NA

    # Reshape matrix for plotting
    cor.melt <- reshape2::melt(cor.mat)
    if(return.n.common) cor.melt$n.common <- n.common


    if(plot == TRUE) {
        g <- ggplot(cor.melt) +
            aes(as.factor(Var1), as.factor(Var2), fill = value) +
            geom_tile() +
            scale_fill_gradient2(low = "steelblue",
                                 mid = "grey90",
                                 high = "red3", midpoint = 0,
                                 space = "Lab",
                                 limits = c(-1, 1),
                                 na.value = NA) +
            labs(x = "", y = "", fill = "r", title = "") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
            labs(title = plot.title)
        print(g)
    }


    name1 <- paste(id.var, 1, sep = "")
    name2 <- paste(id.var, 2, sep = "")

    names(cor.melt)[1:3] <- c(name1, name2, "cor")

    if(return.matrix)
        return(cor.mat)
    else
        return(cor.melt)

}

if(FALSE) {

    x1 <- rnorm(20)
    x2 <- rnorm(20)
    g1 <- rep("a", 20)
    g2 <- rep("b", 20)
    df <- data.frame(x = c(x1, x2),
                     t = c(1:20, 1:20),
                     g = c(g1, g2), stringsAsFactors = FALSE)
    cor1 <- cor(x1, x2)
    cor2 <- dd.cor(df, id.var = "g", cor.var = "x", time.var = "t", plot = FALSE)
    all.equal(cor1, cor2$cor[!is.na(cor2$cor)])

}


## enviro_avg_months ---------------------------------------
enviro_avg_months <- function(data,
                              first.month,
                              last.month,
                              avg.var,
                              month.var = "month",
                              year.var = "year",
                              grid.id.var = NULL,
                              lat.var = "lat",
                              lon.var = "lon",
                              type = avg.var) {
    ## Compute annual multi-month averages for environmental variables
    ##
    ## Input data should be in a "long" format with a column for year, month,
    ## and the environmental variable to average. `first.month` and
    ## `last.month` give the first and last months of a continuous range to
    ## average the environmental variable over.
    ##
    ## If `first.month` is less than `last.month` the environmental variable is
    ## averaged over the months first.month:last.month within each year. For
    ## example, if `first.month` = 3 and `last.month` = 4, the environmental
    ## variable will be averaged over Mar and Apr for each year.
    ##
    ## If `first.month` equals `last.month`, that month is returned with no
    ## averaging.
    ##
    ## If `first.month` is greater than `last.month` the environmental variable
    ## is averaged for year t starting in `first.month` of year t - 1 and ending
    ## in `last.month` of year t. The output year corresponds to the year
    ## January occurs within the average. For example if `first.month` = 12 and
    ## `last.month` = 3, then the average for the environmental variable will
    ## occur over Dec, Jan, Feb, March and the year is specified by the year
    ## for Jan occurs in.
    ##
    ## The function outputs a data.frame with a `year` column and an `index`
    ## column. When `first.month` is greater than `last.month`, the output
    ## data.frame will have one less year than the input data.frame with no
    ## value for the minimum year in the input data frame.
    ##
    ## If `grid.id.var` is non-null, the averaging is done on a per-grid-cell
    ## basis within a year.
    ##
    ## data = a data.frame
    ## first.month = numeric giving the month to start annual average
    ## last.month = numeric giving the month to stop annual average
    ## avg.var = string, column name in `data` of the variable to average
    ## month.var = string, column name in `data` of the month variable
    ## year.var = string, column name in `data` of the year variable
    ## grid.id.var = string, column name in `data` of the grid cell id column
    ## lon.var = string, column name in `data` of for longitude, only used if
    ##           grid.id.var is non-null
    ## lat.var = string, column name in `data` of for latitude, only used if
    ##           grid.id.var is non-null
    ## type = string, value to set a `type` column in the output data.frame.
    ##        Useful if you want to rbind multiple indices together

    if(!is.data.frame(data))
        stop("Input data is not a data.frame", call. = FALSE)

    if(!first.month %in% 1:12 | !last.month %in% 1:12)
        stop("Months not between 1 and 12", call. = FALSE)

    if(!is.numeric(data[ , month.var]) | !is.numeric(data[ , year.var]))
        stop("Month variable must be numeric", call. = FALSE)

    if(first.month < last.month | first.month == last.month) {
        months <- first.month:last.month
        df <- data[data[ , month.var] %in% months, ]
    }

    if(first.month > last.month) {

        ## Remove months prior to `first.month` in first year and
        ## months after `last.month` in last year
        min.yr <- min(data[ , year.var])
        max.yr <- max(data[ , year.var])
        min.rm <- which(data[ , year.var] == min.yr &
                        data[ , month.var] < first.month)
        max.rm <- which(data[ , year.var] == max.yr &
                        data[ , month.var] > last.month)
        sub <- data[-c(min.rm, max.rm), ]

        ## Remove months not being averaged over
        months <- c(first.month:12, 1:last.month)
        sub2 <- sub[sub[ , month.var] %in% months, ]

        ## Create new year index to average over
        sp <- split(sub2, sub2[ , year.var])
        lst <- lapply(sp, function(x) {
                   x$yr.avg <- ifelse(x[ , month.var] %in% first.month:12,
                                      x[ , year.var] + 1, x[ , year.var])
                   return(x)
               })
        df <- do.call("rbind", c(lst, make.row.names = FALSE))
        df[ , year.var] <- df$yr.avg

    }

    ## Calculate averages
    if(is.null(grid.id.var)) {
        sp.avg <- split(df, df[ , year.var])
        lst.avg <- lapply(sp.avg, function(x) {
                    data.frame(year = unique(x[ , year.var]),
                               index = mean(x[ , avg.var]))
                })
        enviro <- do.call("rbind", c(lst.avg, make.row.names = FALSE))
    } else {
        sp.avg <- split(df, list(df[ , grid.id.var], df[ , year.var]))
        lst.avg <- lapply(sp.avg, function(x) {
                    data.frame(year = unique(x[ , year.var]),
                               id = unique(x[ , grid.id.var]),
                               lon = unique(x[ , lon.var]),
                               lat = unique(x[ , lat.var]),
                               index = mean(x[ , avg.var]))
                })
        enviro <- do.call("rbind", c(lst.avg, make.row.names = FALSE))
    }

    enviro$type <- type

    char.months <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug",
                     "sep", "oct", "nov", "dec")
    enviro$months <- paste(char.months[first.month],
                           char.months[last.month], sep = "-")
    enviro$months <- as.factor(enviro$months)

    return(enviro)
}

if(FALSE) {

    set.seed(101)
    yr <- c(rep(1950, 24), rep(1951, 24))
    id <- c(rep(1, 12), rep(2, 12), rep(1, 12), rep(2, 12))
    mt <- rep(1:12, 4)
    lt <- c(rep(48, 12), rep(50, 12), rep(48, 12), rep(50, 12))
    ln <- c(rep(120, 12), rep(122, 12), rep(120, 12), rep(122, 12))
    vl <- rnorm(48)
    df <- data.frame(year = yr, month = mt, lat = lt, lon = ln, id = id, value = vl)

    enviro_avg_months(df, 1, 12, "value", grid.id.var = NULL)
    enviro_avg_months(df, 1, 12, "value", grid.id.var = "id")

    enviro_avg_months(df, 1, 3, "value", grid.id.var = NULL)
    enviro_avg_months(df, 1, 3, "value", grid.id.var = "id")

    enviro_avg_months(df, 12, 3, "value", grid.id.var = NULL)
    enviro_avg_months(df, 12, 3, "value", grid.id.var = "id")

    enviro_avg_months(df, 1, 3, "value")
    enviro_avg_months(df, 3, 3, "value")
    enviro_avg_months(df, 12, 3, "value")
    enviro_avg_months(df, 6, 5, "value", group = "test")
}


## get_npgo ------------------------------------------------
get_npgo <- function(years) {
    ## This function takes as input a range of years and downloads and processes
    ## the NPGO index. The output of the function is a dataframe in 'long'
    ## format with a column for year, month, and the NPGO index.
    ##
    ## years = vector of years

    if(min(years) < 1950)
        stop("Earliest NPGO year is 1950")

    npgo    <- read.table("http://www.o3d.org/npgo/npgo.php", sep = "\t",
                          strip.white = TRUE)
    n.npgo  <- length(npgo[ , 1])
    npgo    <- npgo[4:n.npgo, ]
    n.npgo  <- length(npgo)
    rm.tail <- n.npgo - 3
    npgo    <- npgo[1:rm.tail]
    npgo    <- as.character(npgo)
    npgo    <- strsplit(npgo, "  ")
    npgo    <- do.call("rbind", npgo)
    npgo    <- data.frame(year  = as.numeric(npgo[ , 1]),
                          month = as.numeric(npgo[ , 2]),
                          npgo  = as.numeric(npgo[ , 3]))
    npgo    <- npgo[npgo$year >= min(years) & npgo$year <= max(years), ]

    return(npgo)
}

if(FALSE) {

    get_npgo(1950:1950)
    get_npgo(1950:2019)

    npgo <- get_npgo(1950:2016)
    head(npgo)
    tail(npgo)
    sapply(npgo, class)
    summary(npgo)

}



## get_pdo -------------------------------------------------
get_pdo <- function(years) {
    ## This function takes as input a range of years and downloads and processes
    ## the PDO index. The output of the function is a dataframe in 'long' format
    ## with a column for year, month, and the PDO index.
    ##
    ## years = vector of years

    if(min(years) < 1854)
        stop("Earliest PDO year is 1854")

    pdo <- read.csv("https://www.ncdc.noaa.gov/teleconnections/pdo/data.csv", skip = 1)
    pdo <- as.data.table(pdo)
    pdo[ , date_string := as.character(Date)]
    pdo[ , year := as.numeric(substr(date_string, 1, 4))]
    pdo[ , month := as.numeric(substr(date_string, 5, 6))]
    pdo <- data.table(year = pdo$year, month = pdo$month, pdo = pdo$Value)
    as.data.frame(pdo[year %in% years])

}

if(FALSE) {

    get_pdo(1950:1950)
    get_pdo(1950:2013)

    pdo <- get_pdo(1900:2016)
    head(pdo)
    tail(pdo)
    sapply(pdo, class)
    summary(pdo)

}




