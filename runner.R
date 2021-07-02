require(ALEPlot)
require(cowplot)
require(data.table)
require(doMC)
require(doSNOW)
require(dynatopmodel)
require(foreach)
require(gbm)
require(ggplot2)
require(ggpubr)
require(ggrepel)
require(ggspatial)
require(grid)
require(gridExtra)
require(gstat)
require(hues)
require(iml)
require(lightgbm)
require(lubridate)
require(MASS)
require(mlr)
require(nlme)
require(patchwork) # https://gotellilab.github.io/GotelliLabMeetingHacks/NickGotelli/ggplotPatchwork.html devtools::install_github("thomasp85/patchwork")
require(pbapply)
require(pdist)
require(pdp) # https://bgreenwell.github.io/pdp/articles/pdp-example-xgboost.html
require(plyr)
require(png)
require(psych)
require(purrr)
require(raster)
require(RColorBrewer)
require(readxl)
require(reshape2)
require(rgdal)
require(scales)
require(sf)
require(sp)
require(SpaDES)
require(stringr)
require(stringi)
require(tidyr)
require(tidyverse)
require(unixtools) # devtools::install_github("s-u/unixtools")
require(vegan)
require(velox)
require(R.devices) # https://www.jottr.org/2018/07/21/suppressgraphics/

N_CORE_LARGE = 24 # max. number of cores to use
N_CORE_SMALL = 6 # number of cores for memmory intensive tasks

setwd("/mnt/shared/analysis/Andrews2") # set to project directory

set.tempdir("temp")
setPaths(cachePath="temp", inputPath="temp", modulePath="temp", outputPath="temp", silent = FALSE)

dir.create("data_processed")
dir.create("output")

registerDoMC(N_CORE_LARGE)

start = Sys.time()
source("scripts - Andrews2/utility_functions.R")
source("scripts - Andrews2/000 - conceptual figure.R")
source("scripts - Andrews2/001 - raster setup.R") # 20 min
source("scripts - Andrews2/002 - gridMET comparison.R", encoding = "Latin1") # 7 sec
source("scripts - Andrews2/003 - set up temperature.R") # 27 sec
source("scripts - Andrews2/004 - main models.R") # 45 sec
source("scripts - Andrews2/005 - sd models.R") # 1 sec
source("scripts - Andrews2/006 - spatial prediction.R") # 2 min
source("scripts - Andrews2/007 - ALE plot.R") # 4 sec
Sys.time() - start
