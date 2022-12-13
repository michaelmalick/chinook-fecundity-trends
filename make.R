## Reproduce project
##
## ** make.R last run on 2022-08-24 **

unlink("./figures", recursive = TRUE)
unlink("./outputs", recursive = TRUE)

rm(list = ls())

cat("Sourcing load.R ...", "\n");       source("./load.R")
cat("Sourcing data.R ...", "\n");       source("./data.R")
cat("Sourcing explore.R ...", "\n");    source("./explore.R")
cat("Sourcing fecundity.R ...", "\n");  source("./fecundity.R")
cat("Sourcing age_size.R ...", "\n");   source("./age_size.R")
cat("Sourcing map.R ...", "\n");        source("./map.R")
cat("Sourcing pub.R ...", "\n");        source("./pub.R")

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
