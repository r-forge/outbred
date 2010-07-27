calc.prob <-
function(os = 'unix', compile = FALSE, stepsize = 1, marker.info.file = 'mi.txt', 
	pedigree.file = 'mp.txt', genotype.file = 'mg.txt', output.file = 'cnout') {
		
if (compile) {
	# ----- compile cnF2freq
		system(paste("g++ -I . -O2 cnF2freq.cpp -o cnF2freq_", os, sep = ""))
	}
	# ----- run cnF2freq
	if (os == "windows") {
        system(paste("cnF2freq_", os, sep = ""), 
        input = c(stepsize, '1', marker.info.file, pedigree.file, genotype.file, output.file))    
    } else {
        system(paste("./cnF2freq_", os, sep = ""),
        input = c(stepsize, '1', marker.info.file, pedigree.file, genotype.file, output.file)) 
    }
    
# perl script automatically called still, this means the main file is still called 'cnout'
cat('Transforming genotype probabilities ...\n')
file.copy(paste(R.home(), '/library/qtl.outbred/lp.pl', sep = ''), getwd())
system('perl lp.pl')
cat('\n')

}

