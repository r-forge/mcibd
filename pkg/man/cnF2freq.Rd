\name{cnF2freq}
\Rdversion{1.1}
\alias{cnF2freq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculating Inheritance Probabilities using cnF2freq}
\description{
  \code{cnF2freq} is used to run the \code{cnF2freq} program or even compile \code{cnF2freq} in different operating systems. Usually on a shared cluster, compiling is required for the first time usage.
}
\usage{
cnF2freq(ped.file, chr.file, os = "unix", compile = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ped.file}{a file containing pedigree information. (See \code{pedi.ric} in \code{other.zip} or \code{other.tar} for example.)}
  \item{chr.file}{a file containing chromosome marker information. (See \code{chro.ric} in \code{other.zip} or \code{other.tar} for example.)}
  \item{os}{a string specifying the operating system, typically "unix" or "linux" or "windows".}
  \item{compile}{logical. If \code{TRUE}, \code{cnF2freq} is compiled before running.}
}
\details{
For computers with the user authorities, \code{cnF2freq} may not need to be compiled.
}
\value{
An output file \code{cnF2freq.out} is produced in the working directory.
}
\references{
Nettelblad, C., Holmgren, S., Crooks, L. and Carlborg, O. (2009). \emph{cnF2freq: Efficient Determination of Genotype and Haplotype Probabilities in Outbred Populations Using Markov Models}. Lecture Notes in Bioinformatics (LNBI), 5462: 307-319, Springer-Verlag Berlin Heidelberg.
}
\author{
Xia Shen
}
\seealso{
\code{\link{MCIBD-package}}
}
\examples{
## Unzip other.zip or other.tar in the package library folder. 
## Make sure you are using an operating system with g++ compiler installed.
## Copy the file cnF2freq.cpp in the working directory.
## Execute this function, setting compile = TRUE.

## You may also copy the file cnF2freq_unix or cnF2freq_linux 
## or cnF2freq_windows.exe, then with user authority, 
## this function can be executed without compiling.
}