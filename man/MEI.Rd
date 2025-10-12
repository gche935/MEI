\name{MEI}
\alias{MEI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Measurement Equivalence/Invariance Test
}
\description{
Measurement Equivalence/Invariance Test based on "Cheung, G. W. & Lau, R. S. (2012).  A direct comparison approach for testing measurement invariance.  Organizational Research Methods, 15, 167-198."
Comparison of Defined Parameters based on "Lau, R. S. & Cheung, G. W. (2012).  Estimating and comparing specific mediation effects in complex latent variable models.  Organizational Research Methods, 15, 3-16."
}
\usage{
Full_MEI(model, data.source, Groups, Cluster = "NULL")
CompareLoadings(model, data.source, Groups, Cluster = "NULL", Bootstrap = 0, alpha = 0.01)
CompareMeans(model.PMI, data.source, Groups, Cluster = "NULL", Bootstrap = 0, alpha = 0.01)
CompareParametes(model.PMI, model.PATH, data.source, Groups, Cluster = "NULL", Bootstrap = 0)
MLCompareLoadings(model, data.source, Cluster = "NULL", alpha = 0.05)
LGCompareLoadings(model, data.source, Cluster = "NULL", no.waves = 3, alpha = 0.05)
LGCompareMeans(model.PMI, data.source, Cluster = "NULL", no.waves=3, alpha = 0.01)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
User specified CFA model
}
  \item{model.PMI}{
User specified Partial (or Full) Metric Invariance model
}
  \item{model.PATH}{
User specified Path Model and Defined Parameters for comparisons
}
  \item{data.source}{
A data frame containing the observed variables used in the model.
}
  \item{Groups}{
Grouping variable for cross-group comparisons.
}
  \item{Cluster}{
Cluster  variable for nested data. Monte Carlo simulation method should be used for nested data.
}
  \item{Bootstrap}{
Number of bootstrap samples, must be between 500 and 5,000. If not specify, Monte Carlo simulation (Default) will be used instead of Bootstrapping.
}
  \item{alpha}{
Type I error rate for identifying non-invariant items in List and Delete method. Default is 0.01 (0.05 for MLCFA).  Can also use Bonferroni adjustment (Type I error /No. of comparisons)
}
}
\details{
1.	The significance tests are based on 1 million Monte Carlo simulated results or bootstrap percentile confidence intervals.
2.	Missing values are handled with full-information maximum likelihood.
3.  All non-converged and non-admissible bootstrapped samples are removed in the calculations.
4.  The identification of invariant sets is correct up to a maximum of 25 groups. If there are more than 12 sets of invariant groups, the remaining groups will be considered as non-invariant.
5.  Data files are "Example.A", "Example.B", "Example.C", and "Example.D".
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Cheung, G. W. & Lau, R. S. (2012).  A direct comparison approach for testing measurement invariance.  Organizational Research Methods, 15, 167-198.
Lau, R. S. & Cheung, G. W. (2012). Estimating and comparing specific mediation effects in complex latent variable models. Organizational Research Methods, 15, 3-16.
}
\author{
XXX
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## == Example A - Multiple-Group Comparison == ##

## load data file (Example A.csv)  ##
## Example.A <- read.csv(file = "Example A.csv")  ##
## Groups = "Region" ##


## Specify the measurement model - Model.A
Model.A <- '
        WorkLifeConflict =~ R45a + R45b + R45c + R45d + R45e
        Engagement =~ R90a + R90b + R90c
        Wellbeing =~ R87a + R87b + R87c + R87d + R87e
'

## ===== Full Measurement Invariance Test ===== ##
Full_MEI(Model.A, Example.A, Groups = "Region")

## ===== Compare Loadings ===== ##
CompareLoadings(Model.A, Example.A, Groups = "Region", alpha = 0.001)

## ===== Compare Means ===== ##
CompareMeans(PMI.Model.R, Example.A, Groups = "Region", alpha = 0.001)

## ===== Compare Paths ===== ##
#-- Specify Path model - model.PATH (model.DP) [OrgSize and Tenure are control variables] -- #
model.DP <- '
  Wellbeing ~ Xb1*Engagement + Xc1*WorkLifeConflict + Xd1*OrgSize + Xe1*Tenure
  Engagement ~ Xa1*WorkLifeConflict

  # Defined Parameters #
  IndirectP := Xa1*Xb1
  DirectP := Xc1
  Total := Xa1*Xb1 + Xc1
'

# -- Run function CompareParmeters using  Monte Carlo simulation -- #
CompareParameters(PMI.Model.R, model.DP, Example.A, Groups = "Region")




## == Example B - Panel Data in Longitudinal Studies == ##

## load data file (Example B.csv)  ##
## Example.B <- read.csv("Example B.csv") ##
## All indicators should be ended with _T1, _T2, _T3, etc to indicate the timeframe ##


## Specify the measurement model - Model.B ##
Model.B <- 'SWLC =~ x1 + x2 + x3 + x4 + x5'

## ===== Compare Factor Loadings ===== ##
LGCompareLoadings(Model.B, Example.B, no.waves = 3, alpha =  0.01)

## ===== Compare Means ===== ##
LGCompareMeans(PMI.Model.R, Example.B, no.waves = 3, alpha = 0.01)




## == Example C - Non-independent Data from two sources == ##

## load data file (Example C.csv)  ##
## Example.C <- read.csv("Example C.csv") ##
## All indicators should be ended with _T1 and _T2 to indicate the sources ##


## Specify the measurement model - Model.C ##
Model.C <- 'External =~ x1 + x2 + x3 + x4 + x5 + x6
            Internal =~ x7 + x8 + x9 + x10'

## ===== Compare Factor Loadings ===== ##
LGCompareLoadings(Model.C, Example.C, no.waves = 2, alpha = 0.01)

## ===== Compare Means ===== ##
LGCompareMeans(PMI.Model.R, Example.C, no.waves = 2, alpha = 0.01)




## == Example D - Multilevel Confirmatory Factor Analysis == ##

## load data file (Example D.csv)  ##
## Example.D <- read.csv("Example D.csv") ##
## Cluster = "ID" ##

## Specify the measurement model - Model.D ##
Model.D <- 'OLBI =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8'

## ===== Compare Loadings ===== ##
MLCompareLoadings(Model.D, Example.D, Cluster = "ID", alpha = 0.05)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
