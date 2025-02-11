# visiumR

R source file which contains functions

- [x] to run spatial clustering and markers finding for each sample. Thus easy to be applied on any samples
- [ ] to integrate samples and annotate the cell types
- [ ] to score any gene set spatially

Please clone the directory

testdata were found in testdata/

For utils 1.:  

```{r}
source("visiumR.r")
sample_ids <- list.dirs("testdata/", full.names = F, recursive = F)
lapply(sample_ids, run_single)
```

sessionInfo

```{r}
$platform
$version
'R version 4.1.0 (2021-05-18)'
$os
'Ubuntu 20.04.6 LTS'
$system
'x86_64, linux-gnu'
$ui
'X11'
$language
'(EN)'
$collate
'en_US.UTF-8'
$ctype
'en_US.UTF-8'
$date
'2024-08-12'
$pandoc
'2.19.2'
$packages
A packages_info: 113 × 11
package	ondiskversion	loadedversion	path	is_base	date	source	md5ok
<chr>	<chr>	<chr>	<chr>	<lgl>	<chr>	<chr>	<lgl>
abind	abind	1.4.5	1.4-5	FALSE	FALSE	2016-07-21	CRAN (R 4.1.3)
base64enc	base64enc	0.1.3	0.1-3	FALSE	FALSE	2015-07-28	CRAN (R 4.1.3)
bit	bit	4.0.5	4.0.5	FALSE	FALSE	2022-11-15	CRAN (R 4.1.3)
bit64	bit64	4.0.5	4.0.5	FALSE	FALSE	2020-08-30	CRAN (R 4.1.3)
cli	cli	3.6.1	3.6.1	FALSE	FALSE	2023-03-23	CRAN (R 4.1.3)
cluster	cluster	2.1.4	2.1.4	FALSE	FALSE	2022-08-22	CRAN (R 4.1.3)
codetools	codetools	0.2.19	0.2-19	FALSE	FALSE	2023-02-01	CRAN (R 4.1.3)
colorspace	colorspace	2.1.0	2.1-0	FALSE	FALSE	2023-01-23	CRAN (R 4.1.3)
cowplot	cowplot	1.1.1	1.1.1	TRUE	FALSE	2020-12-30	CRAN (R 4.1.3)
crayon	crayon	1.5.2	1.5.2	FALSE	FALSE	2022-09-29	CRAN (R 4.1.3)
data.table	data.table	1.14.8	1.14.8	FALSE	FALSE	2023-02-17	CRAN (R 4.1.3)
deldir	deldir	1.0.9	1.0-9	FALSE	FALSE	2023-05-17	CRAN (R 4.1.3)
digest	digest	0.6.31	0.6.31	FALSE	FALSE	2022-12-11	CRAN (R 4.1.3)
dplyr	dplyr	1.1.2	1.1.2	TRUE	FALSE	2023-04-20	CRAN (R 4.1.3)
ellipsis	ellipsis	0.3.2	0.3.2	FALSE	FALSE	2021-04-29	CRAN (R 4.1.3)
evaluate	evaluate	0.21	0.21	FALSE	FALSE	2023-05-05	CRAN (R 4.1.3)
fansi	fansi	1.0.4	1.0.4	FALSE	FALSE	2023-01-22	CRAN (R 4.1.3)
fastmap	fastmap	1.1.1	1.1.1	FALSE	FALSE	2023-02-24	CRAN (R 4.1.3)
fitdistrplus	fitdistrplus	1.1.11	1.1-11	FALSE	FALSE	2023-04-25	CRAN (R 4.1.3)
future	future	1.32.0	1.32.0	FALSE	FALSE	2023-03-07	CRAN (R 4.1.3)
future.apply	future.apply	1.11.0	1.11.0	FALSE	FALSE	2023-05-21	CRAN (R 4.1.3)
generics	generics	0.1.3	0.1.3	FALSE	FALSE	2022-07-05	CRAN (R 4.1.3)
ggplot2	ggplot2	3.4.2	3.4.2	TRUE	FALSE	2023-04-03	CRAN (R 4.1.3)
ggrepel	ggrepel	0.9.3	0.9.3	FALSE	FALSE	2023-02-03	CRAN (R 4.1.3)
ggridges	ggridges	0.5.4	0.5.4	FALSE	FALSE	2022-09-26	CRAN (R 4.1.3)
globals	globals	0.16.2	0.16.2	FALSE	FALSE	2022-11-21	CRAN (R 4.1.3)
glue	glue	1.6.2	1.6.2	FALSE	FALSE	2022-02-24	CRAN (R 4.1.3)
goftest	goftest	1.2.3	1.2-3	FALSE	FALSE	2021-10-07	CRAN (R 4.1.3)
gridExtra	gridExtra	2.3	2.3	FALSE	FALSE	2017-09-09	CRAN (R 4.1.3)
gtable	gtable	0.3.3	0.3.3	FALSE	FALSE	2023-03-21	CRAN (R 4.1.3)
scales	scales	1.2.1	1.2.1	FALSE	FALSE	2022-08-20	CRAN (R 4.1.3)
scattermore	scattermore	1.1	1.1	FALSE	FALSE	2023-05-17	CRAN (R 4.1.3)
sctransform	sctransform	0.3.5	0.3.5	TRUE	FALSE	2022-09-21	CRAN (R 4.1.3)
sessioninfo	sessioninfo	1.2.2	1.2.2	FALSE	FALSE	2021-12-06	CRAN (R 4.3.1)
Seurat	Seurat	4.0.3	4.0.3	TRUE	FALSE	2021-06-10	CRAN (R 4.1.0)
SeuratDisk	SeuratDisk	0.0.0.9021	0.0.0.9021	TRUE	FALSE	2024-08-12	Github (mojaveazure/seurat-disk@877d4e18ab38c686f5db54f8cd290274ccdbe295)
SeuratObject	SeuratObject	4.1.3	4.1.3	TRUE	FALSE	2022-11-07	CRAN (R 4.1.3)
shiny	shiny	1.7.4	1.7.4	FALSE	FALSE	2022-12-15	CRAN (R 4.1.3)
sp	sp	1.6.1	1.6-1	FALSE	FALSE	2023-05-31	CRAN (R 4.1.3)
spatstat.core	spatstat.core	2.4.4	2.4-4	FALSE	FALSE	2022-05-18	CRAN (R 4.1.3)
spatstat.data	spatstat.data	3.0.1	3.0-1	FALSE	FALSE	2023-03-12	CRAN (R 4.1.3)
spatstat.geom	spatstat.geom	3.2.1	3.2-1	FALSE	FALSE	2023-05-09	CRAN (R 4.1.3)
spatstat.random	spatstat.random	3.1.5	3.1-5	FALSE	FALSE	2023-05-11	CRAN (R 4.1.3)
spatstat.sparse	spatstat.sparse	3.0.1	3.0-1	FALSE	FALSE	2023-03-12	CRAN (R 4.1.3)
spatstat.utils	spatstat.utils	3.0.3	3.0-3	FALSE	FALSE	2023-05-09	CRAN (R 4.1.3)
stringi	stringi	1.7.6	1.7.6	FALSE	FALSE	2021-11-29	CRAN (R 4.1.1)
stringr	stringr	1.5.0	1.5.0	FALSE	FALSE	2022-12-02	CRAN (R 4.1.3)
survival	survival	3.5.5	3.5-5	FALSE	FALSE	2023-03-12	CRAN (R 4.1.3)
tensor	tensor	1.5	1.5	FALSE	FALSE	2012-05-05	CRAN (R 4.1.3)
tibble	tibble	3.2.1	3.2.1	FALSE	FALSE	2023-03-20	CRAN (R 4.1.3)
tidyr	tidyr	1.3.0	1.3.0	FALSE	FALSE	2023-01-24	CRAN (R 4.1.3)
tidyselect	tidyselect	1.2.0	1.2.0	FALSE	FALSE	2022-10-10	CRAN (R 4.1.3)
utf8	utf8	1.2.3	1.2.3	FALSE	FALSE	2023-01-31	CRAN (R 4.1.3)
uuid	uuid	1.1.0	1.1-0	FALSE	FALSE	2022-04-19	CRAN (R 4.1.3)
uwot	uwot	0.1.14	0.1.14	FALSE	FALSE	2022-08-22	CRAN (R 4.1.3)
vctrs	vctrs	0.6.2	0.6.2	FALSE	FALSE	2023-04-19	CRAN (R 4.1.3)
viridisLite	viridisLite	0.4.1	0.4.1	FALSE	FALSE	2022-08-22	CRAN (R 4.1.3)
withr	withr	2.5.0	2.5.0	FALSE	FALSE	2022-03-03	CRAN (R 4.1.3)
xtable	xtable	1.8.4	1.8-4	FALSE	FALSE	2019-04-21	CRAN (R 4.1.3)
zoo	zoo	1.8.12	1.8-12	FALSE	FALSE	2023-04-13	CRAN (R 4.1.3)
```

## Ackowledgement

Codes for spatial clustering and markers finding are copied from [this github repo](https://github.com/Aagaardlab/placenta-spatial-transcriptomics) along with [this peer reviewed paper](https://pubmed.ncbi.nlm.nih.gov/32681824/)

The original author was invited as the collaborator of this repo.
