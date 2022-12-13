## Data wrangling


## Fecundity -----------------------------------------------
fecund <- data.table::fread("./data/fecundity.csv")


## At least 15 females need to be spawned to include
fecund <- fecund[!is.na(females_spawned), ]
fecund <- fecund[females_spawned >= 15, ]
fecund[ , N := sum(!is.na(fecundity)), by = "id_number"]


## Only keep stocks with 10+ years of data
fecund[ , use10 := ifelse(N < 10, 0, 1), by = "id_number"]
length(unique(fecund$id_number[fecund$use10 == 1]))
fecund <- fecund[use10 == 1, ]


## Add standardized fecundity
fecund[ , fecundity_stnd := scale(fecundity), by = "id_number"]
fecund[ , fecundity_anom := scale(fecundity, scale = FALSE), by = "id_number"]


## Set factor levels + order data
fecund$region   <- factor(fecund$region, levels = unique(fecund$region))
fecund$run      <- factor(fecund$run, levels = unique(fecund$run))
fecund$id_label <- factor(fecund$id_label, levels = unique(fecund$id_label))
fecund$group    <- factor(fecund$group, levels = unique(fecund$group))
fecund$hatchery <- factor(fecund$hatchery, levels = unique(fecund$hatchery))
fecund$id_label <- factor(fecund$id_label, levels = unique(fecund$id_label))

unique(fecund$region)
unique(fecund$run)
unique(fecund$id_label)
unique(fecund$group)
unique(fecund$id_label)


## Summarize dataset
fecund_summary <- fecund[ , .(id_label  = unique(na.omit(id_label)),
                              hatchery  = unique(na.omit(hatchery)),
                              region    = unique(na.omit(region)),
                              run       = unique(na.omit(run)),
                              lon       = unique(na.omit(lon)),
                              lat       = unique(na.omit(lat)),
                              N         = sum(!is.na(fecundity)),
                              N_missing = sum(is.na(fecundity))),
                         by = "id_number"]

fecund_summary$number <- 1:nrow(fecund_summary)
fecund_summary$duplicated <- duplicated(fecund_summary$hatchery)
nrow(fecund_summary)
median(fecund_summary$N)


## Save output
save(fecund, file = "./outputs/fecund.RData")
save(fecund_summary, file = "./outputs/fecund_summary.RData")



## Age + Size ----------------------------------------------

## Read data
cwt_size <- data.table::fread("./data/cwt_size.csv")
cwt_age <- data.table::fread("./data/cwt_age.csv")
usfw_age_size <- data.table::fread("./data/usfw_age_size.csv")


## Remove age-1 females -- implausible
cwt_age <- cwt_age[!(age == 1 & sex == "Female"), ]
cwt_size <- cwt_size[!(age == 1 & sex == "Female"), ]
usfw_age_size <- usfw_age_size[!(age == 1 & sex == "Female"), ]


## Subset for females only
cwt_age <- cwt_age[sex == "Female", ]
cwt_size <- cwt_size[sex == "Female", ]
usfw_age_size <- usfw_age_size[sex == "Female", ]


## Remove really small females
range(cwt_size$length)
range(usfw_age_size$length)
cwt_size <- cwt_size[length >= 300, ]


## Summarize cwt size by year
cwt_size_yr <- cwt_size[ , .(nsize = .N,
                             length_mean = mean(length, na.rm = TRUE),
                             length_sd = sd(length, na.rm = TRUE)),
                        by = .(id_number, id_label, hatchery, run, year, age, sex)]


## Combine cwt age and size
## Only keep rows in cwt_size, which has flagged years removed
cwt_age_size <- merge(cwt_age, cwt_size_yr, all.y = TRUE)


## Set common names
setnames(cwt_age_size, "estimated_number", "nage")
setnames(usfw_age_size, "N", "nsize")
usfw_age_size[ , nage := nsize]


## Check we have sample sizes for all size/age data (both should be 0)
sum(is.na(cwt_age_size$nsize) & !is.na(cwt_age_size$length_mean))
sum(is.na(usfw_age_size$nsize) & !is.na(usfw_age_size$length_mean))


## Need at least 3 fish per year and age group
cwt_age_size <- cwt_age_size[nsize >= 3, ]
usfw_age_size <- usfw_age_size[nsize >= 3, ]


## Add age proportions -- sex specific
cwt_age_size[ , nage_total := sum(nage, na.rm = TRUE),
             by = .(id_number, id_label, hatchery, run, sex, year)]
cwt_age_size[ , nsize_total := sum(nsize, na.rm = TRUE),
             by = .(id_number, id_label, hatchery, run, sex, year)]
cwt_age_size[ , age_prop := nage / nage_total]
cwt_age_size[ , age_prop2 := nsize / nsize_total]

usfw_age_size[ , nage_total := sum(nage, na.rm = TRUE),
             by = .(id_number, id_label, hatchery, run, sex, year)]
usfw_age_size[ , nsize_total := sum(nsize, na.rm = TRUE),
             by = .(id_number, id_label, hatchery, run, sex, year)]
usfw_age_size[ , age_prop := nage / nage_total]
usfw_age_size[ , age_prop2 := nsize / nsize_total]


## Combine cwt + usfw
age_size <- rbindlist(list(cwt_age_size, usfw_age_size), fill = TRUE)
length(unique(age_size$id_label))
age_size_raw <- copy(age_size)

## Add modal age (age most frequently observed)
age_size[ , modal_age := {
    tab <- .SD[ , .(avg = mean(age_prop, na.rm = TRUE)), by = age]
    modal_age <- tab$age[which.max(tab$avg)]
    .(modal_age)
}, by = .(id_number)]


## rename length_mean -> length
setnames(age_size, "length_mean", "length")

## Add release info
rel <- unique(fecund[ , .(id_number, release)])
age_size <- merge(age_size, rel, by = "id_number")

## Add standardized length
age_size[ , length_stnd := scale(length), by = .(id_number, age)]


## Order data
age_size <- age_size[order(id_number, year)]
unique(age_size$id_number)

## Only keep stocks with 10+ years of data
age_size[ , N := sum(!is.na(length)), by = .(id_number, age)]
age_size[ , use10 := ifelse(N < 10, 0, 1), by = .(id_number, age)]
length(unique(age_size$id_number[age_size$use10 == 1]))
age_size <- age_size[use10 == 1, ]


## Set factor levels + order data
age_size$region   <- factor(age_size$region, levels = unique(age_size$region))
age_size$run      <- factor(age_size$run, levels = unique(age_size$run))
age_size$id_label <- factor(age_size$id_label, levels = unique(age_size$id_label))
age_size$group    <- factor(age_size$group, levels = unique(age_size$group))
age_size$hatchery <- factor(age_size$hatchery, levels = unique(age_size$hatchery))
age_size$id_label <- factor(age_size$id_label, levels = unique(age_size$id_label))
levels(age_size$id_label)

## Save outputs
save(cwt_age, file = "./outputs/cwt_age.RData")
save(cwt_size, file = "./outputs/cwt_size.RData")
save(cwt_age_size, file = "./outputs/cwt_age_size.RData")
save(usfw_age_size, file = "./outputs/usfw_age_size.RData")
save(age_size_raw, file = "./outputs/age_size_raw.RData")
save(age_size, file = "./outputs/age_size.RData")



## Combine size and fecundity ------------------------------

## We don't have age-specific fecundity -> need one size value per year
## 1. Calc weighted average of length: weights = age proportion
## 2. Use length of most common age (modal)
size_yr <- age_size[ , .(length_avg = ifelse(all(is.na(length)), NA_real_,
                                                 weighted.mean(length, age_prop,
                                                               na.rm = TRUE)),
                         length_avg2 = ifelse(all(is.na(length)), NA_real_,
                                                 weighted.mean(length, age_prop2,
                                                               na.rm = TRUE)),
                         length_avg3 = ifelse(all(is.na(length)), NA_real_,
                                              mean(length, na.rm = TRUE)),
                         length_modal = length[age == modal_age],
                         modal_age = unique(modal_age)),
                    by = .(id_number, id_label, year)]

## Add size data to fecund data frame, retaining all rows for each id_number in fecund
fecund_size <- fecund[id_number %in% unique(size_yr$id_number), ]
fecund_size <- merge(fecund_size, size_yr, sort = FALSE, all.x = TRUE)

## Update N and use10 columns (don't filter on these though)
fecund_size[ , N := sum(!is.na(length_avg) & !is.na(fecundity)), by = .(id_number)]
fecund_size[ , use10 := ifelse(N < 10, 0, 1), by = .(id_number)]
length(unique(fecund_size$id_number[fecund_size$use10 == 1]))

## Scale data by stock
fecund_size[ , fecundity_stnd := scale(fecundity), by = .(id_number)]
fecund_size[ , fecundity_anom := scale(fecundity, scale = FALSE), by = .(id_number)]
fecund_size[ , length_avg_stnd := scale(length_avg), by = .(id_number)]
fecund_size[ , length_avg_anom := scale(length_avg, scale = FALSE), by = .(id_number)]

## Scale across all stocks
fecund_size$fecundity_anom2  <- scale(fecund_size$fecundity, scale = FALSE)
fecund_size$length_avg_anom2 <- scale(fecund_size$length_avg, scale = FALSE)
fecund_size$fecundity_stnd2  <- scale(fecund_size$fecundity)
fecund_size$length_avg_stnd2 <- scale(fecund_size$length_avg)

## Set factor levels
fecund_size$id_label <- factor(fecund_size$id_label, levels = unique(fecund_size$id_label))
fecund_size$hatchery <- factor(fecund_size$hatchery, levels = unique(fecund_size$hatchery))

## Sanity checks
unique(fecund_size$id_label)
length(unique(fecund_size$id_number))

## Save outputs
save(fecund_size, file = "./outputs/fecund_size.RData")
