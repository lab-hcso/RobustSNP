\name{Rlinear.block}
\alias{Rlinear.block}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Robust SNP association tests for quantitative outcomes allowing for covariates, with genotypes read in blocks
}
\description{
This function is similar to \code{Rlinear}, but the genotypes are read in blocks (instead of reading all of them simultaneously).
The aim to to reduce the memory consumption. This function is particularly suited for large-scale association studies 
like genome-wide association studies (GWAS). 
}
\usage{
Rlinear.block(GenoFile = NA, covPres = TRUE, COVAR = NA, CovarFile = NA, 
SNPrange = NA, covarheader = TRUE, blocksize = 5000, maxpts = 25000, 
abseps = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GenoFile}{
A file location for genotype file in PLINK binary PED format. A collection of three files with the same filename 
but different suffices(bed, bim and fam) should be present. Only the filename but NOT the suffix should be given. 
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
  \item{SNPrange}{
(Optional) A numeric vector specifying the range of SNPs to be tested. 
}
  \item{covarheader}{
(Optional) Whether \code{CovarFile} contains headers. Only applies when \code{CovarFile} is given. 
}
  \item{blocksize}{
The number of SNPs to be read at one time.
}
  \item{maxpts}{
(Optional) Arguments for \code{pmvnorm} in the \code{mvtnorm} package. Specifies the maximum number of function values (as integer).
}
  \item{abseps}{
(Optional) Arguments for \code{pmvnorm} in the \code{mvtnorm} package. Specifies the absolute error tolerance (as double).
}
}
\details{
See \code{\link{Rbin.block}} .
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
 \code{\link{Rlinear}}
}
\examples{
## Assume the binary PED file is named "gwas" and the covariate file is "covar.txt"
## Result = Rlinear.block (GenoFile="gwas" , covPres=TRUE, CovarFile="covar.txt")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ robust }
\keyword{ models }% __ONLY ONE__ keyword per line
