
25th July 2012
Variable description for 
"Structural Change and Cross-Country Growth Empirics"
published in the World Bank Economic Review

Paper available at:
https://sites.google.com/site/medevecon/publications-and-working-papers

This zip-folder contains one datasets and this readme file.

Note that the data file contains 41 countries, thus including CRI which was dropped from the analysis:
drop if wbcode=="CRI"



(1) Helper variables

clist		country list counter (1 to 41); note: "tsset clist year"
nwbcode		numerical isocode
wbcode		string isocode
aggdata		dummy variable: 1 if observation is in the sample across all 4 datasets



(2) Input/Output variables

Suffix -a	refers to aggriculture data
suffix -m	refers to manufacturing data
suffix -PWT	aggregate data from Penn World Table
no suffix	aggregated data: Y=Ya+Ym, K=Ka+Km, L=La+Lm (mimicking an aggregate economy)

Y-		value added (real, 1990 values)
L-		labour
K-		capital stock (real, 1990 values; see appendix for details)
Na		land (arable)

y-		value added per worker
k-		capital stock per worker
na		land per worker

prefix l-	logarithm



(3) Example regression

Agricultural production function
xi: reg lya lLa lka lna i.year, robust (POLS with year dummies)
 