\name{MCIBD-package}
\Rdversion{1.1}
\alias{MCIBD-package}
\docType{package}
\title{
Monte Carlo Identity-By-Descent Matrix Estimation
}
\description{
\code{MCIBD} is a package implemented for flexible identity-by-descent (IBD) matrix estimation in F2 pedigrees. The IBD matrix estimate is approached by a Monte Carlo strategy, where the segregation of the founder alleles is easy to set up. Epistatic IBD matrices can also be calculated for two linked loci. Such IBD matrix estimates are typically used in variance component quantitative trait loci (QTL) analysis. Parallelization using snowfall is implemented to accelerate the Monte Carlo sampling program.
}
\details{
\tabular{ll}{
Package: \tab MCIBD\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2009-10-01\cr
License: \tab GPL-2\cr
Depends: \tab \code{sfsmisc}\cr
Suggests: \tab \code{snow, snowfall}\cr
}
}
\author{
Xia Shen

Maintainer: Xia Shen <xia.shen@lcb.uu.se>
}
\references{
Andersson, L., Haley, C.S., Ellegren, H., Knott, S.A., Johansson, M., Andersson, K., Andersson-Eklund L., Edfors-Lilja, I., Fredholm, M., Hansson, I., Hakansson, J. and Lundstrom, K. (1994). \emph{Genetic Mapping of Quantitative Trait Loci for Growth and Fatness in Pigs}. Science, 263: 1771-1774. \cr \cr
Nettelblad, C., Holmgren, S., Crooks, L. and Carlborg, O. (2009). \emph{cnF2freq: Efficient Determination of Genotype and Haplotype Probabilities in Outbred Populations Using Markov Models}. Lecture Notes in Bioinformatics (LNBI), 5462: 307-319, Springer-Verlag Berlin Heidelberg. \cr \cr
Shen, X., Ronnegard, L. and Carlborg, O. (2009). \emph{A Monte Carlo Full Likelihood Approach in Variance Component Quantitative Trait Loci Analysis}. (Submitted.)
}
