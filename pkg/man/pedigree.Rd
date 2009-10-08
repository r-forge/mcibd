\name{pedigree}
\Rdversion{1.1}
\alias{pedigree}
\docType{data}
\title{A Pig Pedigree Information Dataset}
\description{
(Andersson et al. 1994) An F2 cross has been bred from 2 European wild boars mated to 8 large white sows. Four F1 boars were mated to 22 F1 sows to produce 191 recorded F2 offspring in 26 families.
}
\usage{data(pedigree)}
\format{
  A data frame with 227 observations on the following 2 variables.
  \describe{
    \item{\code{father}}{a numeric vector indicating fathers' ID}
    \item{\code{mother}}{a numeric vector indicating mothers' ID}
  }
}
\references{
Andersson, L., Haley, C.S., Ellegren, H., Knott, S.A., Johansson, M., Andersson, K., Andersson-Eklund L., Edfors-Lilja, I., Fredholm, M., Hansson, I., Hakansson, J. and Lundstrom, K. (1994). \emph{Genetic Mapping of Quantitative Trait Loci for Growth and Fatness in Pigs}. Science, 263: 1771-1774. \cr \cr
}
\seealso{
\code{\link{probabilities}}
}
\examples{
data(pedigree)
}
