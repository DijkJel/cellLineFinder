library(devtools)

se = readRDS('extData/depmap_cosmic.rds')

se1 = se[1:10000,]
se2 = se[10001:19021,]

#create_package('../cellLineFinder')
use_git()
use_mit_license()
use_build_ignore('setup.R')
use_build_ignore('extData/*')
use_data(se, compress = 'xz', overwrite = T)
use_data(se1, compress = 'xz', overwrite = F)
use_data(se2, compress = 'xz', overwrite = F)

use_package('ComplexHeatmap')
use_package('SummarizedExperiment')

use_r('plotHeatmap')
use_r('data')
use_r('findCellLine')
