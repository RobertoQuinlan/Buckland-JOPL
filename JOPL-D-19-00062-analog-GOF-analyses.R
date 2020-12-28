# setwd("your/working/directory/")

##################################################################################
### this script is to run goodness-of-fit and analogue analyses of T1111
### (Buckland Lake) core from Luszczek et al. JOPL-D-19-00062
### using Fortin et al. 2015 as training set data
### Fortin et al. 2015 full training set data (chironomid & environmental data)
### is available at Polar Data Catalogue www.polardata.ca CCIN Ref No. 12504
# user needs to modify spp file "...species-version2-2015.csv"
# by removing site "E511" (outlier in Fortin et al. 2015)
# assuming full chironomid training set (235 lakes) is "TS.csv",
# modified TS is denoted "TS2.csv"; TS2.csv has chironomid data as % relative abundances

TS2<-read.csv(file="TS2.csv", row.names=1)
spp<-TS2[,-1]

# "T1111_Fortin.csv" will be the T1111 core spp file used for analogue analysis
# this core spp file has raw fossil chironomid head capsule counts
# with 78 taxa so that taxonomy is harmonized with Fortin et al. 2015
### this file is available at https://github.com/RobertoQuinlan/Buckland-JOPL

coreinput_analog<-read.csv(file="T1111_Fortin.csv", row.names=1)

# this transforms raw counts to % relative abundances
BCI<-coreinput_analog[,-cbind(1:2)]
core_analog <- BCI / rowSums(BCI) * 100

# Fortin et al. 2015 used filtered >=2%-in-2-lakes spp criteria in their analyses
# "_red" files are those with '2/2' spp filter
# this code is to retain only taxa that are >=2% in 2 intervals in a core
N <- 2
M <- 2
i <- colSums(spp >= M) >= N
spp_red <- spp[, i, drop = FALSE]
# 61 chironomid taxa remain in training set

# "T1111_Fortin_red.csv" will be the core spp file used for GOF analysis,
# matched to training set file "spp_red.csv"
# this file represents "T1111_Fortin.csv" with taxa deleted to match 2/2 filtered TS2
# i.e. "reduced" number of taxa
### user must generate this csv file by applying 2/2 filter to TS2.csv data,
### and matching taxa names with T1111_Fortin.csv file
# "T1111_Fortin_red.csv" has data as % relative abundances

coreinput_gof <- read.csv(file="T1111_Fortin_red.csv", row.names=1)
core_gof <- coreinput_gof[,-cbind(1:2)]

### environmental data also from Polar Data Calalogue in file "...env-version3-2015.csv"
### env data file with E511 removed (outlier in Fortin et al. 2015) here named "env2.csv"
### user must generate this csv file

env2 <- read.csv(file="env2.csv", row.names=1)

# chronological axis for plotting analogue results
chron <- cbind(coreinput_analog[,1] , coreinput_analog[,2])
colnames(chron) <- c("Depth","Year")
chron=as.data.frame(chron)

# chronological axis for plotting goodness-of-fit results
chron_red <- cbind(coreinput_gof[,1] , coreinput_gof[,2])
colnames(chron_red) <- c("Depth","Year")
chron_red=as.data.frame(chron_red)

require("ggpalaeo")
require("analogue")
require("ggplot2")

# goodness-of-fit residuals
rlens <- residLen(spp_red, env2$July, core_gof)
autoplot(rlens, df = data.frame(age = as.numeric(chron_red$Year)), x_axis = "age", fill = c("azure4", "grey","white")) +
  labs(x = "Year", y = "Squared residual distance", fill = "Goodness of fit", categories = c("Good", "Fair", "Poor"))

# analogue distance
# squared chord lengths for Core 
AD <- analogue_distances(spp_red, core_analog)
autoplot(AD, df = data.frame(age = as.numeric(chron$Year)), x_axis = "age", fill = c("azure4", "grey","white")) +
  labs(x = "Year", y = "Squared-chord distance")
