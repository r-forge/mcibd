\name{MCIBD}
\Rdversion{1.1}
\alias{MCIBD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating Identity-By-Descent (IBD) Matrices using Monte Carlo Sampling
}
\description{
\code{MCIBD} estimates the IBD matrix at a given locus or the epistatic IBD matrix for two linked loci. Monte Carlo sampling is used to approach to the matrix estimator, where the output file from \code{\link{cnF2freq}} is required. Segregation of the founder alleles can be set up in the estimation. Parallelzation is available using \code{snowfall} package if executing on a multi-core computer or cluster.
}
\usage{
MCIBD(loci, n.F2, pedigree = NULL, cnF2freq.out = NULL, 
      output.Z = FALSE, read.file = FALSE, segregation = NULL, 
      mc.size = 99, hpc = FALSE, n.cpus = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{loci}{typically an integer specifying a test locus, or a vector containing two linked loci.}
  \item{n.F2}{an integer telling the number of F2 individuals in the pedigree.}
  \item{pedigree}{a matrix or data.frame or file containing pedigree information. For file (when \code{read.file} is \code{TRUE}), see \code{pedi.ric} in \code{other.zip} or \code{other.tar} for example.}
  \item{cnF2freq.out}{a matrix or data.frame containing the output probabilities from cnF2freq. For file (when \code{read.file} is \code{TRUE}), see \code{chro.ric} in \code{other.zip} or \code{other.tar} for example.}
  \item{output.Z}{logical. If \code{TRUE}, average incidence matrix Z is estimated instead of the corresponding IBD matrix.}
  \item{read.file}{logical. If \code{TRUE}, data is loaded from files, where a pedigree information file is specified by \code{pedigree} and a probabilities output file \code{cnF2freq.out} in the working directory.}
  \item{segregation}{a vector specifying the segregation of founder alleles. (See Details.)}
  \item{mc.size}{an integer setting the sample size of the Monte Carlo sampling for one CPU.}
  \item{hpc}{logical. If \code{TRUE}, high performance computing is carried out by parallelization.}
  \item{n.cpus}{an integer telling the number of cores that parallelization is executed on. Only useful when \code{hpc} is \code{TRUE}. Not recommended when mc.size is small.}
}
\details{
To set up the \code{segregation} of the founder alleles, suppose that we have \emph{m} founder(s) from one line and \emph{n} from another line. Then in total there are \emph{2 * (m + n)} alleles needed to be set up. According to the order in the pedigree information file, these alleles for the founders should be filled in \code{segregation} as a vector with length \emph{2 * (m + n)}. For instance, if there is one male in one line and three females in another line, the vector might be created like \code{c(1,1,2,2,3,3,2,3)}, which means that the male has two identical alleles whereas the genotypes of the females are (2,2), (3,3) and (2,3), respectively. If \code{NULL}, all the founder alleles are assumed to be different from each other.
}
\value{
An output IBD matrix is saved in the working directory as a file named by the loci with extension \code{.ibd}.
}
\references{
Shen, X., Ronnegard, L. and Carlborg, O. (2009). \emph{A Monte Carlo Full Likelihood Approach in Variance Component Quantitative Trait Loci Analysis}. (Submitted.)
}
\author{
Xia Shen
}
\seealso{
\code{\link{cnF2freq}}, \code{\link{MCIBD.chro}}, \code{\link{MCIBD.epi2chro}}, \code{\link{MCIBD-package}}
}
\examples{
data(pedigree)
data(probabilities)

## Calculation on one CPU
MCIBD(loci = 80, n.F2 = 191, pedigree = pedigree, 
      cnF2freq.out = probabilities, mc.size = 5)
## IBD matrix of dimension 191 x 191 at locus 80 is accomplished, 
## where 5 imputes were sampled.

## Calculation on 2 CPUs with segregation of the founder alleles
MCIBD(loci = 90, n.F2 = 191, pedigree = pedigree,
      segregation = c(rep(1,4),rep(2,16)),
      cnF2freq.out = probabilities, mc.size = 5, hpc = TRUE, n.cpus = 2)
## IBD matrix of dimension 191 x 191 at locus 90 is accomplished, 
## where 10 imputes were sampled.

## Calculation of epistatic IBD matrix
MCIBD(loci = c(88, 99), n.F2 = 191, pedigree = pedigree,
      segregation = c(rep(1,4),rep(2,16)),
      cnF2freq.out = probabilities, mc.size = 5, hpc = TRUE, n.cpus = 2)
## IBD matrix of dimension 191 x 191 for linked lici 88 and 99 is accomplished, 
## where 10 imputes were sampled.
}