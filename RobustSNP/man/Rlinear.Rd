\name{Rlinear}
\alias{Rlinear}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Robust SNP association tests for quantitative outcomes allowing for covariates
}
\description{
This function performs robust genetic association tests for quantitative  outcomes. SNP associations are tested under 
additive, dominant and recessive models. Valid p-values are returned for each SNP, taking into account the multiple
testing under different inheritance models. Covariates may also be included. 
}
\usage{
Rlinear(data = NA, outcome = NA, SNP = NA, InputFile = NA, covPres = TRUE, 
COVAR = NA, CovarFile = NA, inputheader = TRUE, covarheader = TRUE, 
maxpts = 25000, abseps = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A matrix or dataframe specifying the genotypes. Each row represents an individual and each column 
represents a SNP. SNPs should be coded as 0,1 or 2 according to the allelic counts.  
}
  \item{outcome}{
The quantitative outcome.
}
  \item{SNP}{
SNP labels. If not specified, SNPs will be labelled by numbers. 
}
  \item{InputFile}{
An optional txt file specifying the outcome in the first column and the genotypes in the other columns.
Each row represents an individual. SNPs should be coded as 0,1 or 2 according to the allelic counts. Users may either specify
both \code{data} and \code{outcome} or just specify \code{InputFile}. 
}
  \item{covPres}{
Whether covariates are present. If \code{covPres} is \code{TRUE}, the user must specify either 
\code{COVAR} OR \code{CovarFile}. 
}
  \item{COVAR}{
(Optional) A matrix or dataframe specifying the covariates. Each row represents an individual and each column 
represents a covariate. 
}
  \item{CovarFile}{
(Optional) A file location for covariates. The file should be readable by the \code{read.table} command in R. 
}
  \item{inputheader}{
(Optional) Whether \code{InputFile} contains headers. Only applies when \code{InputFile} is given.
}
  \item{covarheader}{
(Optional) Whether \code{CovarFile} contains headers. Only applies when \code{CovarFile} is given. 
}
  \item{maxpts}{
(Optional) Arguments for \code{pmvnorm} in the \code{mvtnorm} package. Specifies the maximum number of function values (as integer).
}
  \item{abseps}{
(Optional) Arguments for \code{pmvnorm} in the \code{mvtnorm} package. Specifies the absolute error tolerance (as double).
}
}
\details{
See \code{\link{Rbin}}.
}
\value{
Same as the ouput from  \code{\link{Rbin}}.
}
\references{
Hon-Cheong So and Pak C. Sham (2010).Robust association tests under different genetic models, allowing for 
binary or quantitative traits and covariates. Submitted. 
}
\author{
Hon-Cheong So 
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Rlinear.block}}
}
\examples{
## Simuate a case-control study with 300 cases/300 controls
## and 10 SNPs
set.seed(3)
outcome1 = rnorm(600,mean=30,sd=10) 
data1 = matrix(nrow=600, ncol=10) 
for (k in 1:10) {
data1[,k]=sample ( rep(0:2,200) ) }
COVAR1= rnorm(600,mean=100,sd=30) 

## with covariates
Rlinear( data= data1, outcome=outcome1, SNP=1:10,covPres=TRUE,
COVAR= COVAR1)

##without covariates
Rlinear( data= data1, outcome=outcome1, SNP=1:10,covPres=FALSE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ robust }
\keyword{ models }% __ONLY ONE__ keyword per line