Suppose that we measure, with error, a series of ``effects", $\beta_1,\dots,\beta_J$. To take just one concrete 
example, $\beta_j$ could be 
the difference in the mean (log) expression levels of gene $j$ ($j=1,\dots,J$) between 2 conditions. In this case, a measurement of
the effect might be the difference in sample means obtained in the two conditions. We will let $\hat\beta_j$ denote the measured value of $\beta_j$, and assume that
each measurement comes with an associated standard error, 
$s_j$. A key aim here will be to take proper account of the fact that some measurements may be more precise than others: that is, to take proper account of variation in $s_j$ across $j$.

A common goal, particularly in genomic studies, is to identify which $\beta_j$ differ from zero. This is commonly tackled
by first computing a $p$ value for each $j$, testing $H_j:\beta_j=0$, and then applying methods based on False Discovery Rates, such as the qvalue package. For example, we could compute $Z$ scores $Z_j = \beta_j/s_j$, and translate $Z_j$ to a $p$ value, 
$p_j$, and then use now-standard methods to select a $p$ value threshold to control the FDR. 

There are two issues with this approach that I would like to address here. The first is that it really does not take proper account of the measurement errors. To see this, consider an example where half the measurements are quite precise, and the other half are really, really, poor. Intuitively, the poor measurements tell us nothing, and any sane analysis should effectively ignore them. However, in a standard FDR-type analysis, these poor measurements add ``noise" and affect estimated FDRs. This is because the $p$ values from the poor measurements will be effectively uniformly distributed, and some will be significant at any given threshold. 

To see this, take an example. The following R code
simulates some data where the effects are normally distributed,
and the first 500 observations have good precision, and the next 500 have poor precision.


```r
# install q value package
source("http://bioconductor.org/biocLite.R")
```

```
## BiocInstaller version 1.4.9, ?biocLite for help
```

```
## A newer version of Bioconductor is available for this version of R,
## ?BiocUpgrade for help
```

```r
biocLite("qvalue")
```

```
## BioC_mirror: http://bioconductor.org
```

```
## Using R version 2.15, BiocInstaller version 1.4.9.
```

```
## Installing package(s) 'qvalue'
```

```
## 
## The downloaded binary packages are in
## 	/var/folders/1f/d96lz9ts15g81dqs_hcc_9cr0000gq/T//Rtmp2NOcrf/downloaded_packages
```

```
## Old packages: 'boot', 'class', 'cluster', 'foreign', 'knitr', 'lattice',
## 'MASS', 'Matrix', 'mgcv', 'nlme', 'nnet', 'rpart', 'spatial', 'survival',
## 'wavethresh'
```

```r
library("qvalue")
```

```
## Loading Tcl/Tk interface ...
```

```
## Warning: Can't find a usable tk.tcl in the following directories:
## /System/Library/Frameworks/Tcl.framework/Versions/8.5/Resources/Scripts/tk8.5
## /System/Library/Frameworks/Tcl.framework/Versions/8.5/Resources/Scripts/tk8.5/Resources/Scripts
## /System/Library/Frameworks/Tcl.framework/Versions/8.5/Resources/tk8.5
## /System/Library/Frameworks/Tcl.framework/Versions/8.5/Resources/tk8.5/Resources/Scripts
## ./lib/tk8.5 ./lib/tk8.5/Resources/Scripts ~/Library/Tcl/tk8.5
## ~/Library/Tcl/tk8.5/Resources/Scripts /Library/Tcl/tk8.5
## /Library/Tcl/tk8.5/Resources/Scripts /System/Library/Tcl/tk8.5
## /System/Library/Tcl/tk8.5/Resources/Scripts /System/Library/Tcl/8.5/tk8.5
## /System/Library/Tcl/8.5/tk8.5/Resources/Scripts ~/Library/Frameworks/tk8.5
## ~/Library/Frameworks/tk8.5/Resources/Scripts /Library/Frameworks/tk8.5
## /Library/Frameworks/tk8.5/Resources/Scripts
## /System/Library/Frameworks/tk8.5
## /System/Library/Frameworks/tk8.5/Resources/Scripts ./library
## 
## This probably means that tk wasn't installed properly.
```

```
## done
```

```r

# set up some data with mixture of values of s
set.seed(100)
s.good = 0.5
s.poor = 10
J.good = 500
J.poor = 500
J = J.good + J.poor
beta = c(rnorm(J, 0, 1))
s = c(rep(s.good, J.good), rep(s.poor, J.poor))
betahat = beta + rnorm(J, 0, s)
# compute the usual zscore and corresponding p value
zscore = betahat/s
p = pchisq(zscore^2, df = 1, lower.tail = F)
```

As expected the $p$ values from the poor observations are approximately uniform, whereas those from the good observations  are enriched for small $p$ values:


```r
p.good = p[1:J.good]
p.poor = p[(J.good + 1):J]
hist(p, main = "p values from all observations")
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-21.png) 

```r
hist(p.good, main = "p values from good observations")
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-22.png) 

```r
hist(p.poor, main = "p values from poor observations")
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-23.png) 


Now what happens if we apply FDR methods to the all the data is
that the uniform $p$ values from the poor observations
add noise relative to looking at the good observations only. This impacts both the estimate of $pi_0$ from qvalue, and the number of findings that are significant at a given FDR.

```r
qq.all = qvalue(p)
qq.good = qvalue(p.good)
print(c(qq.good$pi0, qq.all$pi0))
```

```
## [1] 0.3368 0.6932
```










The two main assumptions we will make here
are that 
* these effects are exchangeable: that is (by de Finetti's theorem) they can be
thought of as being independent and identically distributed from some unknown distribution $g(\beta)$. 
* these effects are symetrically distributed about 0.

Specifically we will assume a parameteric form for $g$
as a mixture of 0-centered normal distributions

A key issue I want to address here is that if these measurements are made with different precisions, then we want to
take this into account in our analysis.

That is, we assume that $\hat\beta_p \sim






```r
# install q value package, and load in ash code
source("ash.R")
source("http://bioconductor.org/biocLite.R")
biocLite("qvalue")
```

```
## BioC_mirror: http://bioconductor.org
```

```
## Using R version 2.15, BiocInstaller version 1.4.9.
```

```
## Installing package(s) 'qvalue'
```

```
## 
## The downloaded binary packages are in
## 	/var/folders/1f/d96lz9ts15g81dqs_hcc_9cr0000gq/T//Rtmp2NOcrf/downloaded_packages
```

```
## Old packages: 'boot', 'class', 'cluster', 'foreign', 'knitr', 'lattice',
## 'MASS', 'Matrix', 'mgcv', 'nlme', 'nnet', 'rpart', 'spatial', 'survival',
## 'wavethresh'
```

```r
library("qvalue")

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

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 

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

