\name{probabilities}
\Rdversion{1.1}
\alias{probabilities}
\docType{data}
\title{A Probailities Dataset of Pig Chromosome 6 from cnF2freq}
\description{
(Andersson et al. 1994) An F2 cross has been bred from 2 European wild boars mated to 8 large white sows. Four F1 boars were mated to 22 F1 sows to produce 191 recorded F2 offspring in 26 families. The genetic information on chromosome 6 came from 22 genotyped microsatellite markers at: 0.0, 8.6, 36.6, 49.7, 50.5, 62.9, 79.2, 80.4, 83.7, 84.1, 84.8, 90.6, 95.4, 100.7, 101.9, 115.9, 116.7, 119.0, 120.2, 124.0, 127.0 and 170.9 cM.
}
\usage{data(probabilities)}
\format{
  A data frame with 38817 observations on 64 variables, where 38817 is the size of the pedigree 227 times the number of test loci 171 (every 1 cM per locus), and the 64 variables are the probabilities calculated from \code{cnF2freq}.
}
\references{
Andersson, L., Haley, C.S., Ellegren, H., Knott, S.A., Johansson, M., Andersson, K., Andersson-Eklund L., Edfors-Lilja, I., Fredholm, M., Hansson, I., Hakansson, J. and Lundstrom, K. (1994). \emph{Genetic Mapping of Quantitative Trait Loci for Growth and Fatness in Pigs}. Science, 263: 1771-1774. \cr \cr
}
\seealso{
\code{\link{pedigree}}
}
\examples{
data(probabilities)
}
