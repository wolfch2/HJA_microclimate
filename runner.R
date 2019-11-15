# assumes data_input is present and has the appropriate files

LD_LIBRARY_PATH=/home/stats/wolfch/lib:/cm/local/apps/gcc/5.1.0/lib64:/cm/shared/apps/gdal/2.1.1/lib:/cm/shared/apps/gdal/2.1.1/bin:/home/stats/wolfch/my_libs:/home/stats/wolfch/R_libs/3.5.1:/cm/shared/apps/R/3.5.1/lib64/R:/home/stats/wolfch/wood/curl/wolf_curl/lib:/cm/shared/apps/R/3.2.2/lib64/R:/cm/shared/apps/slurm/14.11.11/lib64/slurm:/cm/shared/apps/slurm/14.11.11/lib64:/cm/shared/apps/sge/2011.11p1/lib/linux-x64:/cm/local/apps/gcc/5.1.0/lib:/cm/local/apps/gcc/5.1.0/lib64
export LD_RUN_PATH=/cm/shared/apps/gdal/2.1.1/bin
export PKG_CONFIG_PATH=/home/stats/wolfch/lib/pkgconfig
export R_LIBS=/home/stats/wolfch/R_libs/3.5.1
export TMPDIR=/home/stats/wolfch/temp

/cm/shared/apps/R/3.5.1/bin/R

# to install caret if /tmp full (https://stackoverflow.com/questions/17841332/using-install-packages-with-custom-temp-dir)
#        Sys.setenv(TMPDIR="/home/stats/wolfch/temp2")
#        configure.vars="/home/stats/wolfch/temp2"
#        install.packages("mlr")

require(ALEPlot)
require(cowplot)
require(data.table)
require(doMC)
require(doSNOW)
require(dynatopmodel)
require(foreach)
require(gbm)
require(ggplot2)
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
require(StreamThermal) # devtools::install_github('tsangyp/StreamThermal')
require(stringr)
require(tidyr)
require(tidyverse)
require(unixtools) # devtools::install_github("s-u/unixtools")
require(vegan)
require(velox)
require(xgboost) # https://xgboost.readthedocs.io/en/latest/R-package/xgboostPresentation.html
require(xgboostExplainer) # devtools::install_github("AppliedDataSciencePartners/xgboostExplainer")

#setwd("/home/stats/wolfch/Andrews2")
#set.tempdir("/home/stats/wolfch/temp2")
setwd("/home/chrisgraywolf/analysis_desktop/Andrews2")
set.tempdir("/home/chrisgraywolf/analysis_desktop/temp")
setPaths(cachePath="temp", inputPath="temp", modulePath="temp", outputPath="temp", silent = FALSE)

dir.create("data_output")
dir.create("output")

source("scripts - Andrews2/utility_functions.R")

registerDoMC(24)

system("scripts - Andrews2/001 - raster setup.R")
system("scripts - Andrews2/002 - site reduce.R")




