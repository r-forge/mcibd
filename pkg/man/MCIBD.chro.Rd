\name{MCIBD.chro}
\Rdversion{1.1}
\alias{MCIBD.chro}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating Identity-By-Descent (IBD) Matrices along A Whole Chromosome
}
\description{
\code{MCIBD.chro} estimates the IBD matrices along a given chromosome. Monte Carlo sampling is used to approach to the matrix estimator, where the output file from \code{\link{cnF2freq}} is required. Segregation of the founder alleles can be set up in the estimation. Parallelzation is available using \code{snowfall} package if executing on a multi-core computer or cluster.
}
\usage{
MCIBD.chro(dis = 5, n.F2, pedigree = NULL, cnF2freq.out = NULL, 
           IBD.type = "genotypic", output.Z = "none", read.file = FALSE, 
		   segregation = NULL, mc.size = 99, hpc = FALSE, n.cpus = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dis}{the distance between two nearest test loci in centimorgan.}
  \item{n.F2}{an integer telling the number of F2 individuals in the pedigree.}
  \item{pedigree}{a matrix or data.frame or file containing pedigree information. For file (when \code{read.file} is \code{TRUE}), see \code{pedi.ric} in \code{other.zip} or \code{other.tar} for example.}
  \item{cnF2freq.out}{a matrix or data.frame containing the output probabilities from cnF2freq. For file (when \code{read.file} is \code{TRUE}), see \code{chro.ric} in \code{other.zip} or \code{other.tar} for example.}
  \item{IBD.type} {a string specifying the output type of IBD matrix, \code{"genotypic"} or \code{"gametic"}.}
  \item{output.Z}{a string indicates the type of output. If \code{"av"}, average incidence matrix Z is estimated instead of the corresponding IBD matrix; If \code{"all"}, the IBD matrix is still estimated as default but with all the incidence matrix imputes output in folder "Zall"; If \code{"pc"}, an incidence matrix containing the principle components of the corresponding IBD matrix is output.}
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
Output IBD matrices along the whole chromosome are saved in the working directory as files named by the loci with extension \code{.ibd}.
}
\references{
Shen, X., Ronnegard, L. and Carlborg, O. (2009). \emph{Flexible Identity-By-Descent Matrix Estimation with Given Base Generation Structures in F2 Intercross Designs}. (Manuscript)
}
\author{
Xia Shen
}
\seealso{
\code{\link{MCIBD}}, \code{\link{MCIBD.epi2chro}}, \code{\link{MCIBD-package}}
}
\examples{
data(pedigree)
data(probabilities)

MCIBD.chro(80, 191, pedigree, probabilities, output.Z = "pc", mc.size = 10)
}