
```r
# install q value package, and load in ash code
source("ash.R")
source("http://bioconductor.org/biocLite.R")
```

```
## Bioconductor version 2.11 (BiocInstaller 1.8.3), ?biocLite for help
```

```r
biocLite("qvalue")
```

```
## BioC_mirror: http://bioconductor.org
```

```
## Using Bioconductor version 2.11 (BiocInstaller 1.8.3), R version 2.15.
```

```
## Installing package(s) 'qvalue'
```

```
## 
## The downloaded binary packages are in
## 	/var/folders/0t/0teHHybKEZSU9PtJ-X37NU+++TQ/-Tmp-//RtmpnE6LHT/downloaded_packages
```

```
## Old packages: 'foreign', 'lattice', 'MASS', 'nlme', 'rpart', 'survival',
## 'wavethresh'
```

```r
library("qvalue")
```

```
## Loading Tcl/Tk interface ...
```

```
## done
```

```r

# set up some data with mixture of values of sigmaa
set.seed(100)
# sebetahat = sample(c(1,0.1),size=1500,replace=T)
sebetahat = 0.01 * rgamma(1500, 0.1, 0.1)
beta = c(rnorm(500, 0, 1), rnorm(500, 0, 0.5), rnorm(500, 0, 1e-06))
betahat = beta + rnorm(1500, 0, sebetahat)

# compute the usual zscore and corresponding p value
zscore = betahat/sebetahat
pval = pchisq(zscore^2, df = 1, lower.tail = F)

# apply the ash method to do adaptive shrinkage of betahat and the qvalue
# method to get
beta.ash = ash(betahat, sebetahat)
qq = qvalue(pval)
```


Now we want to see how the results compare

```r

attach(beta.ash)
```

```
## The following object(s) are masked from 'package:base':
## 
##     sample
```

```r
conf = ifelse(PositiveProb > 0.5, PositiveProb, 1 - PositiveProb)
sum(conf > 0.95)
```

```
## [1] 1152
```

```r
sum(qq$qvalues < 0.05)
```

```
## [1] 1200
```

So there are more significant at FDR=0.05, than at local bayes FDR 0.05. But this is not suprising because of the difference between local FDR and FDR.

Indeed, let's look at the errors

```r
err = (sign(betahat) != sign(beta))
table(err[conf > 0.95])
```

```
## 
## FALSE  TRUE 
##  1149     3
```

```r
table(err[qq$qvalues < 0.05])
```

```
## 
## FALSE  TRUE 
##  1187    13
```

```r


# check whether ordering by q values does better or worse job than
# ordering by confidence, in terms of identifying betas with the right
# sign


plot(cumsum(err[order(qq$qvalues)]), type = "l")
lines(cumsum(err[order(conf, decreasing = TRUE)]), col = 2)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 

```r

# note: I edited nejm_brca_release.txt by removing columns 1-3
hh = read.table("../data/nejm_brca_release_edit.csv", sep = ",", skip = 3)
subset = apply(hh, 1, max) < 20
hh = hh[subset, ]

labs = read.table("../data/nejm_brca_release_edit.csv", sep = ",", skip = 1, 
    nrows = 1)
labs = 1 * (labs == "BRCA1") + 2 * (labs == "BRCA2")

hh.betahat = apply(hh[, labs == 1], 1, mean) - apply(hh[, labs == 2], 1, mean)
n1 = sum(labs == 1)
n2 = sum(labs == 2)
hh.sebetahat = sqrt(apply(hh[, labs == 1], 1, var)/n1 + apply(hh[, labs == 2], 
    1, var)/n2)

hh.zscore = hh.betahat/hh.sebetahat
hh.pval = pchisq(hh.zscore^2, df = 1, lower.tail = F)
hh.q = qvalue(hh.pval)
sum(hh.q$q < 0.05)
```

```
## [1] 278
```
