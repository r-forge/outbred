impo.prob <-
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

