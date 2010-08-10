calc.prob <-
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
