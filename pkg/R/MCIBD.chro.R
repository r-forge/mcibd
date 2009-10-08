MCIBD.chro <-
function(dis = 5, n.F2, pedigree = NULL, cnF2freq.out = NULL, output.Z = FALSE, read.file = FALSE, segregation = NULL, mc.size = 99, hpc = FALSE, n.cpus = 2) {
## --------------------------------------------- ##
##       Monte Carlo IBD Matrix Calculator       ##
##      cnF2freq has to be run before this.      ##
##      Xia.Shen@lcb.uu.se ---2009-10-01---      ##
## --------------------------------------------- ##
	
## !! PACKAGE REQUIREMENTS: sfsmisc, snow (for HPC), snowfall (for HPC)

## Check some inputs
if (prod(c(dis > 0, n.F2 > 0, mc.size > 0, n.cpus > 0)) == 0) {
	stop("Non-positive input detected!")
}	

## Data & Pre-calculation ##
pb <- txtProgressBar(style = 3)
if (read.file) {
	cat("Data Loading ...", "\n")
	pedi <- read.table(pedigree, header = TRUE)
	setTxtProgressBar(pb, .5)
	carl <- read.table(cnF2freq.out, skip = 2, fill = TRUE)
	setTxtProgressBar(pb, 1)
	cat("\n")
	cat("OKAY.", "\n")
	cat("\n")
} else {
	pedi <- pedigree
	carl <- cnF2freq.out
}
pedi <- as.matrix(pedi)
carl <- as.matrix(carl)
ped.size <- nrow(pedi)
n.loci <- nrow(carl)/ped.size
loci <- seq(0, n.loci - 1, dis)
f2id <- (ped.size - n.F2 + 1):ped.size

cat("Data Transforming ...", "\n")
carlout <- NULL
for (k in 1:ped.size) {
	carlout[[k]] <- carl[((k-1)*n.loci + 1):(k*n.loci),]
	setTxtProgressBar(pb, k/ped.size)
}
cat("\n")
cat("OKAY.", "\n")
cat("\n")

if (is.null(segregation)) {
	n.founder <- 0
	is.founder <- sum(pedi[n.founder + 1,]) == 0
	while (is.founder) {
		n.founder <- n.founder + 1
		is.founder <- sum(pedi[n.founder + 1,]) == 0
	}
	segregation <- 1:(2*n.founder)
}

## Necessary Functions ##
cat("Functions Loading ...", "\n")
## Sampling using cnF2freq Output
require(sfsmisc, quietly = TRUE)
binary <- NULL
for (i in 1:64) {binary <- rbind(binary, as.numeric(digitsBase(i-1,,6)))}
samplecarl <- function(carlout, position, pedigree, f2id) {
	here <- Z <- NULL
	for (id in f2id) {here <- rbind(here, carlout[[id]][position,])}
	for (i in 1:length(f2id)) {
		case <- binary[sample(1:64,1,prob=here[i,]),]
		if (case[6]==0) {
			if (case[5]==1) {
				a1 <- pedigree[pedigree[f2id[i],1],1]*2 - 1
			}
			else {
				a1 <- pedigree[pedigree[f2id[i],1],1]*2
			}
		}
		else {
			if (case[4]==1) {
				a1 <- pedigree[pedigree[f2id[i],1],2]*2 - 1
			}
			else {
				a1 <- pedigree[pedigree[f2id[i],1],2]*2
			}
		}
		if (case[3]==0) {
			if (case[2]==1) {
				a2 <- pedigree[pedigree[f2id[i],2],1]*2 - 1
			}
			else {
				a2 <- pedigree[pedigree[f2id[i],2],1]*2
			}
		}
		else {
			if (case[1]==1) {
				a2 <- pedigree[pedigree[f2id[i],2],2]*2 - 1
			}
			else {
				a2 <- pedigree[pedigree[f2id[i],2],2]*2
			}
		}
		Z.line <- rep(0,length(segregation))
		Z.line[c(a1,a2)] <- 1
		Z <- rbind(Z,Z.line)
	}
	dimnames(Z) <- list(NULL,NULL)
	return(Z)
}
setTxtProgressBar(pb, .5)
samplehere <- function(here, pedigree, f2id) {
	Z <- NULL
	for (i in 1:length(f2id)) {
	case <- binary[sample(1:64,1,prob=here[i,]),]
		if (case[6]==0) {
			if (case[5]==1) {
				a1 <- pedigree[pedigree[f2id[i],1],1]*2 - 1
			}
			else {
				a1 <- pedigree[pedigree[f2id[i],1],1]*2
			}
		}
		else {
			if (case[4]==1) {
				a1 <- pedigree[pedigree[f2id[i],1],2]*2 - 1
		}
			else {
				a1 <- pedigree[pedigree[f2id[i],1],2]*2
			}
		}
		if (case[3]==0) {
			if (case[2]==1) {
				a2 <- pedigree[pedigree[f2id[i],2],1]*2 - 1
			}
			else {
				a2 <- pedigree[pedigree[f2id[i],2],1]*2
			}
		}
		else {
			if (case[1]==1) {
				a2 <- pedigree[pedigree[f2id[i],2],2]*2 - 1
			}
			else {
				a2 <- pedigree[pedigree[f2id[i],2],2]*2
			}
		}
		Z.line <- rep(0,length(segregation))
		Z.line[c(a1,a2)] <- 1
		Z <- rbind(Z,Z.line)
	}
	dimnames(Z) <- list(NULL,NULL)
	return(Z)
}
	
## Adding Segregation Vector for Founder Alleles
sgg <- function(Z, seg.alleles) {
	if (length(seg.alleles) != ncol(Z)) {
		stop("Incorrect number of alleles!")
	}	
	alleles <- seg.alleles[1]
	for (i in 2:length(seg.alleles)) {
		if (seg.alleles[i] != seg.alleles[i - 1]) {
			alleles <- c(alleles, seg.alleles[i])
		}	
	}
	seg.Z <- matrix(0, nrow(Z), length(alleles))
	for (i in 1:length(alleles)) {
		index <- 1:ncol(Z)*(seg.alleles == alleles[i])
		if (sum(index != 0) > 1) {
			seg.Z[,i] <- rowSums(Z[,index])
		}
		else {
			seg.Z[,i] <- Z[,index]	
		}
	}
	return(seg.Z)
}
setTxtProgressBar(pb, 1)
cat("\n")
cat("OKAY.", "\n")
cat("\n")

## Monte Carlo Sampling ##
MIsize <- mc.size
nf2 <- n.F2
if (!hpc) {
	## Single Core
	t0 <- proc.time()[3]
	cat("One master is doing its jobs ...", "\n")
	for (p in loci) {
		sumPi <- matrix(0, nf2, nf2)
		for(i in 1:MIsize) {
			Z <- samplecarl(carlout = carlout, position = p + 1, pedigree = pedi, f2id = f2id)
			Z <- sgg(Z, segregation)
			if (!output.Z) {
				Pi <- .5*Z%*%t(Z)
			}
			else {
				Pi <- Z
			}
			sumPi <- sumPi + Pi
			setTxtProgressBar(pb, i/MIsize)
		}
		cat("\n")
		meanPi <- sumPi/MIsize
		filename <- paste(p, ".ibd", sep = "")
		write.table(meanPi, filename, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
		cat(paste("IBD matrix of dimension", nf2, "x", nf2, "at locus", p, "is accomplished, where", MIsize, "imputes were sampled.", sep = " "), "\n")
	}
	t1 <- proc.time()[3] - t0
}
else {
	## Multiple Cores
	## NOTE! Make sure the multicore environment is set up and snowfall package is installed.
	require(snow, quietly = TRUE)
	require(snowfall, quietly = TRUE)
	sfInit(parallel = TRUE, cpus = n.cpus, type = "SOCK")
	slavejob <- function(idx) {
		sumPi <- matrix(0, nf2, nf2)
		for(i in 1:MIsize) {
			Z <- samplehere(here = here, pedigree = pedi, f2id = f2id)
			Z <- sgg(Z, segregation)
			if (!output.Z) {
				Pi <- .5*Z%*%t(Z)
			}
			else {
				Pi <- Z
			}
			sumPi <- sumPi + Pi
		}
		meanPi <- sumPi/MIsize
		return(meanPi)
	}
	t0 <- proc.time()[3]
	for (p in loci) {
		here <- NULL
		for (id in f2id) {here <- rbind(here, carlout[[id]][p + 1,])}
		cat("Broadcasting objects to slaves ...", "\n")
		sfExport("MIsize")
		setTxtProgressBar(pb, 1/8)
		sfExport("binary")
		setTxtProgressBar(pb, 2/8)
		sfExport("f2id")
		setTxtProgressBar(pb, 3/8)
		sfExport("nf2")
		setTxtProgressBar(pb, 4/8)
		sfExport("pedi")
		setTxtProgressBar(pb, 5/8)
		sfExport("here")
		setTxtProgressBar(pb, 6/8)
		sfExport("samplehere")
		setTxtProgressBar(pb, 7/8)
		sfExport("sgg")
		setTxtProgressBar(pb, 8/8)
		cat("\n")
		cat("OKAY.", "\n")
		cat("\n")
		cat("Slaves are doing their jobs ...", "\n")
		result <- sfLapply(1:n.cpus, slavejob)
		cat("\n")
		cat("Final calculation for the IBD matrix ...", "\n")
		sumPi <- matrix(0, nf2, nf2)
		for(sn in 1:n.cpus) {
			sumPi <- sumPi + result[[sn]]
			setTxtProgressBar(pb, sn/n.cpus)
		}
		cat("\n")
		meanPi <- sumPi/n.cpus
		filename <- paste(p, ".ibd", sep = "")
		write.table(meanPi, filename, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
		cat(paste("IBD matrix of dimension", nf2, "x", nf2, "at locus", p, "is accomplished, where", MIsize*n.cpus, "imputes were sampled.", sep = " "), "\n")
		cat("\n")
		rm(result)
	}
	t1 <- proc.time()[3] - t0
	sfStop()
}
cat("\n")
cat("MC Sampling Execution Time: ", t1, "sec", "\n")
cat("\n")
## End ##
}

