.onAttach <- 
function(...)
{
 cat("qtl.outbred: QTL mapping for outbred line crosses\n")
 cat('Version 2011.03.17 installed\n')
}




.onLoad <- 
function(...)
{
 if (require(qtl)) {
 	scanone.rqtl <<- scanone
 	scanone <<- scanone.raw
 }
}




`calc.prob` <-
function(os = 'unix', stepsize = 1, marker.info.file = 'mi.txt', 
	pedigree.file = 'mp.txt', genotype.file = 'mg.txt', output.file = 'cnout') {
		
#if (compile) {
	# ----- compile cnF2freq
#		system(paste("g++ -I . -O2 cnF2freq.cpp -o cnF2freq_", os, sep = ""))
#	}
	# ----- run cnF2freq

	if (os == "windows") {
		if(file.exists("cnF2freq_windows.exe")==FALSE){
			file.copy(paste(R.home(), '/library/qtl.outbred/cnF2freq_windows.zip', sep = ''), getwd())
			zip.file.extract('cnF2freq_windows.exe',zipname="cnF2freq_windows.zip",dir = getwd())
			zip.file.extract('libgcc_s_dw2-1.dll',zipname="cnF2freq_windows.zip",dir = getwd())
		}
        system(paste("cnF2freq_", os, sep = ""), 
        input = c(stepsize, '1', marker.info.file, pedigree.file, genotype.file, output.file))    
    
	} else {#unix here
		if(file.exists("cnF2freq_unix")==FALSE){
			file.copy(paste(R.home(), '/library/qtl.outbred/cnF2freq_unix.zip', sep = ''), getwd())
			zip.file.extract('cnF2freq_unix',zipname="cnF2freq_unix.zip",dir = getwd())
		}		
        system(paste("./cnF2freq_", os, sep = ""),
        input = c(stepsize, '1', marker.info.file, pedigree.file, genotype.file, output.file)) 
    }
    
# perl script automatically called still, this means the main file is still called 'cnout'
cat('Transforming genotype probabilities ...\n')
file.copy(paste(R.home(), '/library/qtl.outbred/lp.pl', sep = ''), getwd())
system('perl lp.pl')
cat('\n')
}





`impo.prob` <-
function(cross.data, path = NULL, stepsize = 1, Grid = FALSE) {

#transforming grid probabilities	
if (Grid){
	cat('Transforming genotype probabilities from GridQTL...\n')
	file.copy(paste(R.home(), '/library/qtl.outbred/lp2.pl', sep = ''), getwd())
	system('perl lp2.pl')
	cat('\n')
}

real.chrno <- nchr(cross.data)
real.indno <- nind(cross.data)
real.crosstype <- class(cross.data)[1]
	
data.list <- list()
probno.per.chr <- marker.per.chr <- numeric(real.chrno)
cat('Importing genotype probabilities from files ...\n')
pb <- txtProgressBar(style = 3)
for (i in 1:real.chrno) {
	loading.name <- paste(path, 'p_output_chrom_', i, '.txt', sep='')
	data.list[[i]] <- read.table(loading.name)
	probno.per.chr[i] <- length(data.list[[i]][,1])
	marker.per.chr[i] <- probno.per.chr[i]/real.indno
	setTxtProgressBar(pb, i/real.chrno)
}
cat('\n\n')
	
# ----- creating a fake dataset with prob spaces
fake.map <- sim.map(len = c(marker.per.chr), n.mar = marker.per.chr, eq.spacing = TRUE, include.x = FALSE)
fake.data <- sim.cross(fake.map, type = real.crosstype, n.ind = real.indno)
fake.data$pheno <- cross.data$pheno
fake.data <- calc.genoprob(fake.data, step = 0, error.prob = 0.01)

# ----- filling in the genotype probs
cat('Swapping genotype probabilities ...\n')
for (i.chrom in 1:real.chrno) { #for each chromosome
	cat('Chromosome ', i.chrom, '\n')
	newvec <- c(1:length(fake.data$geno[[i.chrom]]$prob))
	counter <- 1
	for (i.old in 1:marker.per.chr[i.chrom]){
		counter2 <- i.old
		for(ii.old in 1:real.indno){
			newvec[counter] <- data.list[[i.chrom]]$V2[counter2] #aa
			newvec[counter + probno.per.chr[i.chrom]]<- data.list[[i.chrom]]$V3[counter2] #ab
			newvec[counter + probno.per.chr[i.chrom] + probno.per.chr[i.chrom]]<- data.list[[i.chrom]]$V4[counter2] #bb
			counter  <- counter + 1
			counter2  <- counter2 + marker.per.chr[i.chrom]
		}
	setTxtProgressBar(pb, i.old/marker.per.chr[i.chrom])
	}
	cat('\n')
	# ----- refilling the genoprobs in fakedata (1 to chrno)
	fake.data$geno[[i.chrom]]$prob[1:probno.per.chr[i.chrom]] <- newvec[1:probno.per.chr[i.chrom]]
	fake.data$geno[[i.chrom]]$prob[(probno.per.chr[i.chrom] + 1):((probno.per.chr[i.chrom])*2)] <- newvec[(probno.per.chr[i.chrom] + 1):((probno.per.chr[i.chrom])*2)]
	fake.data$geno[[i.chrom]]$prob[((probno.per.chr[i.chrom]*2) + 1):((probno.per.chr[i.chrom])*3)] <- newvec[((probno.per.chr[i.chrom]*2) + 1):((probno.per.chr[i.chrom])*3)]
}
cat('\n\n')
	
# ----- compressing the imported genotype probs to fit
compression <- 1/stepsize	
cat('Rescaling the data ...\n')
for (i.chrom in 1:real.chrno){
	fake.map[[i.chrom]] <- fake.map[[i.chrom]]/compression
	fake.data$geno[[i.chrom]]$prob <- fake.data$geno[[i.chrom]]$prob/compression
	fake.data$geno[[i.chrom]]$map <- fake.data$geno[[i.chrom]]$map/compression
	attr(fake.data$geno[[i.chrom]]$prob, 'map') <- attr(fake.data$geno[[i.chrom]]$prob, 'map')/compression
	setTxtProgressBar(pb, i.chrom/real.chrno)
}
cat('\n')

summary(fake.data)

return(fake.data)

}





`scanone.raw` <- function(cross, chr, pheno.col=1, model=c("normal","binary","2part","np"),
         method=c("em","imp","hk","ehk","mr","mr-imp","mr-argmax"),
         addcovar=NULL, intcovar=NULL, weights=NULL,
         use=c("all.obs", "complete.obs"), upper=FALSE,
         ties.random=FALSE, start=NULL, maxit=4000, tol=1e-4,
         n.perm, perm.Xsp=FALSE, perm.strata=NULL, verbose, batchsize=250,
         n.cluster=1, ind.noqtl) 
{
res <- scanone.rqtl(cross, chr, pheno.col, model, method, addcovar, intcovar, weights, use, upper,
         ties.random, start, maxit, tol,
         n.perm, perm.Xsp, perm.strata, verbose, batchsize,
         n.cluster, ind.noqtl)
class(res) <- c('outbred', class(res))
return(res)
}




`plot.outbred` <- function(x, old.data = data, chr = -1, gap = 25, ...){

if (inherits(x, 'scanone')) {

if (chr[1] == -1) { # plot all the chromosomes   
	spacer <- gap
	plot.scanone(x, incl.markers = FALSE, gap = spacer) 
	incer <- 0
	scale2 <- 0
	for (i in 1:nchr(old.data)) {
		points((old.data$geno[[i]]$map+incer), rep(-1, length(old.data$geno[[i]]$map)), pch = "|", cex = 0.5, ylim = c(0, max(x$lod)))
		incer <- incer + ceiling(old.data$geno[[i]]$map[length(old.data$geno[[i]]$map)]) + spacer
	}
} else { # only do some chromosomes
	spacer <- gap
	plot.scanone(x, incl.markers = FALSE, gap = spacer, chr = chr) 
	incer <- 0
	chr <- sort(chr)
	for (i in 1:length(chr)) {
		points((old.data$geno[[chr[i]]]$map+incer), rep(-1, length(old.data$geno[[chr[i]]]$map)), pch = "|", cex = 0.5, ylim = c(0, max(x$lod)))
 		incer <- incer + ceiling(old.data$geno[[chr[i]]]$map[length(old.data$geno[[chr[i]]]$map)]) + spacer
	}
}

}
else {
	plot.scantwo(x, ...)
}

}






