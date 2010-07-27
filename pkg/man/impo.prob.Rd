\name{impo.prob}
\alias{impo.prob}
\title{
Import genotype probabilities into an existing cross object
}
\description{
The function imports genotype probabilities into an existing qtl object. 
}
\usage{
impo.prob(cross.data, stepsize = 1)
}
\arguments{
  \item{cross.data}{A cross object typically read in via \code{read.cross} in package \code{qtl}.}
  \item{stepsize}{An integer giving the distance that the two consecutive test loci are spaced in centi-morgan.}
}
\details{
The data in the object of class \code{cross} is used to create a new object of the same class, while importing the genotype probabilities from a set of files. It is necessary to import create the original object prior to using the impo.prob function (see \code{\link{read.cross}} for more information). 

The genotype probabilities are imported from a file for each chromosome separately. In this format first column indicates the individual id. The next three columns are the genotype probabilities of being homozygous, heterozygous and homozygous for the other alleles respectively (i.e. AA, AB and BB). A separate row for each centi-morgan interval per individual is provided. The rows are arranged to have the probabilities per individual together (i.e. from position 1 to end of chromosome for individual 1 followed by position 1 to end of chromosome for individual 2 etc.). It is important to note that the filenames should be for each chromosome as follows (where A indicates the chromosome number): 'p_output_chrom_A.txt'. This is also the file output created by \code{calc.prob}.

A set of example files with the required format described above are provided in the package installation folder. 'p_output_chrom_1.txt' and 'p_output_chrom_2.txt' contains the genotype probabilities for the two chromosomes to be imported. 'example_data.csv' is an F2 intercross file that can be read by \code{\link{read.cross}} function in package \code{qtl}.
}
\value{
The new cross object \code{cross.data} is returned, with the genotype probabilities imported.
}
\author{
Ronnie Nelson and Xia Shen
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{calc.prob}}, \code{\link{qtl.outbred-package}}, \code{\link{read.cross}}
}
\examples{
# ----- impo.prob example
file <- paste(R.home(), '/library/qtl.outbred/example_data.csv', sep = '')
data <- read.cross('csv', file = file, genotypes = c('AA', 'AB', 'BB'))
outbred.data <- impo.prob(cross.data = data, stepsize = 5)

# ----- QTL scan using Haley-Knott regression and plotting
hk.model <- scanone(outbred.data, method = 'hk', pheno.col = 2) 
plot_outbred(hk.model, data)
}