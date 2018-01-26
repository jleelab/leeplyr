# leeplyr

## How to install

Open up R and run the following code:
```
if (!require(devtools)) {
  install.packages("devtools")
}
devtools::install_github("jleelab/leeplyr")
```

Then load the package:
```
library(leeplyr)
```
Try out some commands:
```
cytosol.folder<-system.file('data/roi/RoiSet_cytosol', package='leeplyr')
nuclei.folder<-system.file('data/roi/RoiSet_nuclei', package='leeplyr')
transcripts<-system.file('data/fisseq/res_001_001FISSEQ.out', package='leeplyr')
map.to.roi(transcripts, roi.folder = c(cytosol.folder, nuclei.folder), roi.labels =  c('cytoplasm', 'nucleus'), OR.cutoff = 5, p.value.cutoff = 20)
```