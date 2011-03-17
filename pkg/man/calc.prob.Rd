\name{calc.prob}
\alias{calc.prob}
\title{
Calculate genotype probabilities for outbred intercrosses using cnF2freq 
}
\description{
The function calculates genotype probabilities from outbred intercross data using published software 'cnF2freq'. The software is compiled for Windows and Unix operating systems. However the source code is included in the package installation folder or can be downloaded from: http://user.it.uu.se/~carln/cnF2freq.html
}
\usage{
calc.prob(os = 'unix', stepsize = 1, 
	marker.info.file = 'mi.txt', pedigree.file = 'mp.txt', 
	genotype.file = 'mg.txt', output.file = 'cnout')
}
\arguments{
  \item{os}{A string specifying the operating system, typically 'unix' or 'windows'.}
  \item{stepsize}{An integer giving the distance that the two consecutive test loci are spaced in centi-morgan.}
  \item{marker.info.file}{A marker information file required by \code{cnF2freq}, see Details.}
  \item{pedigree.file}{A pedigree information file required by \code{cnF2freq}, see Details.}
  \item{genotype.file}{A genotype data file required by \code{cnF2freq}, see Details.}
  \item{output.file}{An output file name required by \code{cnF2freq}, see Details.}
}
\details{
Three files are required for the calculating of genotype probabilities. All the files are tab-delimited. The marker information files provides information on the number of chromosomes, the total number of markers, followed by the number of per chromosomes and the name of each. Finally the space between each marker for each chromosome is provided. Note that the first column of each row should be filled with a '1' followed by the position data (each chromosome in a single row). The genotype file provide the genotypic information for each individual. The individuals IDs are indicated in the first column of each row. The following columns are filled with integer values indicating the genotype at each locus in sequence. Every 2 columns corresponds to one marker (one allele in each column) and the columns are arranged sequentially in the order described in the marker file. The marker pedigree file are arranged in fullsib families. For each family the the number of F2 individuals within the family are provided. This is followed wit the individual ID's in the first column (starting with the F0 generation, then the F1s and then the F2s). For each individual his/her parent's ID's are provided in the next two columns followed by it's sex. Note for the F0 generation their parents are indicated with a '0' while an additional column with their line origins are provided.
 
Three example files with the required format described above are provided in the package installation folder. 'mp.txt', 'mi.txt', and 'mg.txt' are the pedigree information file, marker information file, and the marker genotype file, respectively.

Note that Perl is required for some of the functions in the package \code{qtl.outbred}.  and can be downloaded from:
http://www.perl.org/
 
}
\value{
The function generates a file named by \code{output.file} for each chromosome containing genotype probabilities that can be thereafter imported using \code{impo.prob}.
}
\author{
Ronnie Nelson and Xia Shen
}
%\note{
%Further notes...
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{impo.prob}}, \code{\link{qtl.outbred-package}}
}
\examples{
\dontrun{
# ----- calc.prob example (calculating the genotype probabilities from the files 'mi.txt', 'mp.txt' and 'mg.txt' at 5 cM intervals and produce the files 'p_chrom_A.txt',where A indicates the chromosome number for each chromosome).  
setwd(paste(R.home(), '/library/qtl.outbred/', sep = ''))
 calc.prob(stepsize = 5)
}
}