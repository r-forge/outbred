plot_outbred <- function(model, old.data, dimension = '1D', chr = -1, gap = 25, ...){

if (dimension == '1D') {

if (chr[1] == -1) { # plot all the chromosomes   
	spacer <- gap
	plot.scanone(model, incl.markers = FALSE, gap = spacer) 
	incer <- 0
	scale2 <- 0
	for (i in 1:nchr(old.data)) {
		points((old.data$geno[[i]]$map+incer), rep(-1, length(old.data$geno[[i]]$map)), pch = "|", cex = 0.5, ylim = c(0, max(model$lod)))
		incer <- incer + ceiling(old.data$geno[[i]]$map[length(old.data$geno[[i]]$map)]) + spacer
	}
} else { # only do some chromosomes
	spacer <- gap
	plot.scanone(model, incl.markers = FALSE, gap = spacer, chr = chr) 
	incer <- 0
	chr <- sort(chr)
	for (i in 1:length(chr)) {
		points((old.data$geno[[chr[i]]]$map+incer), rep(-1, length(old.data$geno[[chr[i]]]$map)), pch = "|", cex = 0.5, ylim = c(0, max(model$lod)))
 		incer <- incer + ceiling(old.data$geno[[chr[i]]]$map[length(old.data$geno[[chr[i]]]$map)]) + spacer
	}
}

}
else {
	scantwo(model, ...)
}

}
