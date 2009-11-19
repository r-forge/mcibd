MCIBD.epi2chro <-
function(dis = 5, n.F2, pedigree = NULL, cnF2freq.out1 = NULL, cnF2freq.out2 = NULL, IBD.type = "genotypic", output.Z = "none", read.file = FALSE, segregation = NULL, mc.size = 99, hpc = FALSE, n.cpus = 2) {
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
	carl1 <- read.table(cnF2freq.out1, skip = 2, fill = TRUE)
	carl2 <- read.table(cnF2freq.out2, skip = 2, fill = TRUE)
	setTxtProgressBar(pb, 1)
	cat("\n")
	cat("OKAY.", "\n")
	cat("\n")
} else {
	pedi <- pedigree
	carl1 <- cnF2freq.out1
	carl2 <- cnF2freq.out2
}
pedi <- as.matrix(pedi)
carl1 <- as.matrix(carl1)
carl2 <- as.matrix(carl2)
ped.size <- nrow(pedi)
n1.loci <- nrow(carl1)/ped.size
n2.loci <- nrow(carl2)/ped.size
loci1 <- seq(0, n1.loci - 1, dis)
loci2 <- seq(0, n2.loci - 1, dis)
f2id <- (ped.size - n.F2 + 1):ped.size

cat("Data Transforming ...", "\n")
carlout1 <- carlout2 <- NULL
for (k in 1:ped.size) {
	carlout1[[k]] <- carl1[((k-1)*n1.loci + 1):(k*n1.loci),]
	carlout2[[k]] <- carl2[((k-1)*n2.loci + 1):(k*n2.loci),]
	setTxtProgressBar(pb, k/ped.size)
}
cat("\n")
cat("OKAY.", "\n")
cat("\n")


n.founder <- 0
is.founder <- sum(pedi[n.founder + 1,]) == 0
while (is.founder) {
	n.founder <- n.founder + 1
	is.founder <- sum(pedi[n.founder + 1,]) == 0
}
if (is.null(segregation)) {
	segregation <- 1:(2*n.founder)
}

## Necessary Functions ##
cat("Functions Loading ...", "\n")
## Sampling using cnF2freq Output
require(sfsmisc, quietly = TRUE)
binary <- NULL
for (i in 1:64) {binary <- rbind(binary, as.numeric(digitsBase(i-1,,6)))}
samplecarl <- function(carlout, position, pedigree, f2id, type) {
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
		if (type == "genotypic") {
			Z.line <- rep(0,length(segregation))
			Z.line[c(a1, a2)] <- 1
			Z <- rbind(Z, Z.line)
		}
		if (type == "gametic") {
			Z.line1 <- Z.line2 <- rep(0,length(segregation))
			Z.line1[a1] <- Z.line2[a2] <- 1
			Z <- rbind(Z, Z.line1, Z.line2)
		}
	}
	dimnames(Z) <- list(NULL,NULL)
	return(Z)
}
setTxtProgressBar(pb, .5)
samplehere <- function(here, pedigree, f2id, type) {
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
		if (type == "genotypic") {
			Z.line <- rep(0,length(segregation))
			Z.line[c(a1, a2)] <- 1
			Z <- rbind(Z, Z.line)
		}
		if (type == "gametic") {
			Z.line1 <- Z.line2 <- rep(0,length(segregation))
			Z.line1[a1] <- Z.line2[a2] <- 1
			Z <- rbind(Z, Z.line1, Z.line2)
		}
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
if (IBD.type == "genotypic") nf2 <- n.F2 else nf2 <- 2*n.F2
if (output.Z == "none") {
	exname <- ".ibd"
	type <- "IBD"
	dim2 <- nf2
}
if (output.Z == "av" | output.Z == "pc") {
	exname <- ".z"
	type <- "Incidence"
	dim2 <- 2*n.founder	
}
	
if (!hpc) {
	## Single Core
	t0 <- proc.time()[3]
	cat("One master is doing its jobs ...", "\n")
	for (p1 in loci1) {
		for (p2 in loci2) {
			if (output.Z == "av") sumPi <- matrix(0, nf2, dim2) else sumPi <- matrix(0, nf2, nf2)
			for(i in 1:MIsize) {
				Z1 <- samplecarl(carlout = carlout1, position = p1 + 1, pedigree = pedi, f2id = f2id, type = IBD.type)
				Z2 <- samplecarl(carlout = carlout2, position = p2 + 1, pedigree = pedi, f2id = f2id, type = IBD.type)
				Z1 <- sgg(Z1, segregation)
				Z2 <- sgg(Z2, segregation)
				if (output.Z == "none" | output.Z == "pc") {
					Pi1 <- .5*Z1%*%t(Z1)
					Pi2 <- .5*Z2%*%t(Z2)
				}
				if (output.Z == "av") {
					Pi1 <- Z1
					Pi2 <- Z2
				}
				## Hadamard product applied (Shen et al. 2009)
				Pi <- Pi1*Pi2
				sumPi <- sumPi + Pi
				setTxtProgressBar(pb, i/MIsize)
			}
			cat("\n")
			meanPi <- sumPi/MIsize
			filename <- paste(paste(p1, p2, sep = "_x_"), exname, sep = "")
			write.table(meanPi, filename, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
			cat(paste(type, "matrix of dimension", nf2, "x", dim2, "for epistatic loci", p1, "and", p2, "is accomplished, where", MIsize, "imputes were sampled.", sep = " "), "\n")
		}
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
		if (output.Z == "av") sumPi <- matrix(0, nf2, dim2) else sumPi <- matrix(0, nf2, nf2)
		for(i in 1:MIsize) {
			Z1 <- samplehere(here = here1, pedigree = pedi, f2id = f2id, type = IBD.type)
			Z2 <- samplehere(here = here2, pedigree = pedi, f2id = f2id, type = IBD.type)
			Z1 <- sgg(Z1, segregation)
			Z2 <- sgg(Z2, segregation)
			if (output.Z == "none" | output.Z == "pc") {
				Pi1 <- .5*Z1%*%t(Z1)
				Pi2 <- .5*Z2%*%t(Z2)
			}
			if (output.Z == "av") {
				Pi1 <- Z1
				Pi2 <- Z2
			}
			## Hadamard product applied (Shen et al. 2009)
			Pi <- Pi1*Pi2
			sumPi <- sumPi + Pi
		}
		meanPi <- sumPi/MIsize
		return(meanPi)
	}
	t0 <- proc.time()[3]
	for (p1 in loci1) {
		for (p2 in loci2) {
			here1 <- here2 <- NULL
			for (id in f2id) {
				here1 <- rbind(here1, carlout1[[id]][p1 + 1,])
				here2 <- rbind(here2, carlout2[[id]][p2 + 1,])
			}
			cat("Broadcasting objects to slaves ...", "\n")
			sfExport("MIsize")
			setTxtProgressBar(pb, 1/10)
			sfExport("binary")
			setTxtProgressBar(pb, 2/10)
			sfExport("f2id")
			setTxtProgressBar(pb, 3/10)
			sfExport("nf2")
			setTxtProgressBar(pb, 4/10)
			sfExport("dim2")
			setTxtProgressBar(pb, 5/10)
			sfExport("pedi")
			setTxtProgressBar(pb, 6/10)
			sfExport("here1")
			setTxtProgressBar(pb, 7/10)
			sfExport("here2")
			setTxtProgressBar(pb, 8/10)
			sfExport("samplehere")
			setTxtProgressBar(pb, 9/10)
			sfExport("sgg")
			setTxtProgressBar(pb, 10/10)
			cat("\n")
			cat("OKAY.", "\n")
			cat("\n")
			cat("Slaves are doing their jobs ...", "\n")
			result <- sfLapply(1:n.cpus, slavejob)
			cat("\n")
			cat("Final calculation for the", type, "matrix ...", "\n")
			if (output.Z == "av") sumPi <- matrix(0, nf2, dim2) else sumPi <- matrix(0, nf2, nf2)
			for(sn in 1:n.cpus) {
				sumPi <- sumPi + result[[sn]]
				setTxtProgressBar(pb, sn/n.cpus)
			}
			cat("\n")
			meanPi <- sumPi/n.cpus
			if (output.Z == "pc") {
				A <- eigen(meanPi)
				v <- A$values[1:(2*n.founder)]
				pc <- A$vectors[,1:(2*n.founder)]
				meanPi <- pc%*%diag(sqrt(v))
			}
			filename <- paste(paste(p1, p2, sep = "_x_"), exname, sep = "")
			write.table(meanPi, filename, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
			cat(paste(type, "matrix of dimension", nf2, "x", dim2, "for epistatic loci", p1, "and", p2, "is accomplished, where", MIsize*n.cpus, "imputes were sampled.", sep = " "), "\n")
			cat("\n")
			rm(result)
		}
	}
	t1 <- proc.time()[3] - t0
	sfStop()
}
cat("\n")
cat("MC Sampling Execution Time: ", t1, "sec", "\n")
cat("\n")
## End ##
}

