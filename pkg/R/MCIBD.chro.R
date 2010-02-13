MCIBD.chro <-
function(dis = 5, n.F2, pedigree = NULL, cnF2freq.out = NULL, IBD.type = "genotypic", output.Z = "none", read.file = FALSE, segregation = NULL, mc.size = 99, hpc = FALSE, n.cpus = 2) {
## --------------------------------------------- ##
##       Monte Carlo IBD Matrix Calculator       ##
##      cnF2freq has to be run before this.      ##
##      Xia.Shen@lcb.uu.se ---2009-10-17---      ##
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
		for (id in f2id) {here <- rbind(here, rbind(carlout[[id]])[position,])}
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
		dimnames(Z) <- list(NULL, NULL)
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
if (output.Z == "all") dir.create("Zall")
if (output.Z == "none" | output.Z == "all") {
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
	for (p in loci) {
		if (output.Z == "av") sumPi <- matrix(0, nf2, dim2) else sumPi <- matrix(0, nf2, nf2)
		if (output.Z == "all") dir.create(paste("Zall/", p, sep = ""))
		for(i in 1:MIsize) {
			Z <- samplecarl(carlout = carlout, position = p + 1, pedigree = pedi, f2id = f2id, type = IBD.type)
			Z <- sgg(Z, segregation)
			if (output.Z == "none" | output.Z == "pc") {
				Pi <- .5*Z%*%t(Z)
			}
			if (output.Z == "av") {
				Pi <- Z
			}
			if (output.Z == "all") {
				Pi <- .5*Z%*%t(Z)
				filename <- paste("Zall/", p, "/", i, ".z", sep = "")
				write.table(Z, filename, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
			}
			sumPi <- sumPi + Pi
			setTxtProgressBar(pb, i/MIsize)
		}
		cat("\n")
		meanPi <- sumPi/MIsize
		if (output.Z == "pc") {
			A <- eigen(meanPi)
			v <- A$values[1:(2*n.founder)]
			pc <- A$vectors[,1:(2*n.founder)]
			meanPi <- pc%*%diag(sqrt(v))
		}
		filename <- paste(p, exname, sep = "")
		write.table(meanPi, filename, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
		cat(paste(type, "matrix of dimension", nf2, "x", dim2, "at locus", p, "is accomplished, where", MIsize, "imputes were sampled.", sep = " "), "\n")
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
			Z <- samplehere(here = here, pedigree = pedi, f2id = f2id, type = IBD.type)
			Z <- sgg(Z, segregation)
			if (output.Z == "none" | output.Z == "pc") {
				Pi <- .5*Z%*%t(Z)
			}
			if (output.Z == "av") {
				Pi <- Z
			}
			if (output.Z == "all") {
				Pi <- .5*Z%*%t(Z)
				dir.create("Zall")
				filename <- paste("Zall/", p, "/", i, ".z", sep = "")
				write.table(Z, filename, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
			}
			sumPi <- sumPi + Pi
		}
		meanPi <- sumPi/MIsize
		return(meanPi)
	}
	t0 <- proc.time()[3]
	for (p in loci) {
		if (output.Z == "all") dir.create(paste("Zall/", p, sep = ""))
		here <- NULL
		for (id in f2id) {here <- rbind(here, carlout[[id]][p + 1,])}
		cat("Broadcasting objects to slaves ...", "\n")
		sfExport("MIsize")
		setTxtProgressBar(pb, 1/10)
		sfExport("p")
		setTxtProgressBar(pb, 2/10)
		sfExport("binary")
		setTxtProgressBar(pb, 3/10)
		sfExport("f2id")
		setTxtProgressBar(pb, 4/10)
		sfExport("nf2")
		setTxtProgressBar(pb, 5/10)
		sfExport("dim2")
		setTxtProgressBar(pb, 6/10)
		sfExport("pedi")
		setTxtProgressBar(pb, 7/10)
		sfExport("here")
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
		filename <- paste(p, exname, sep = "")
		write.table(meanPi, filename, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
		cat(paste(type, "matrix of dimension", nf2, "x", dim2, "at locus", p, "is accomplished, where", MIsize*n.cpus, "imputes were sampled.", sep = " "), "\n")
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

