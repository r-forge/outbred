\name{qtl.outbred-package}
\alias{qtl.outbred-package}
\alias{qtl.outbred}
\docType{package}
\title{
Interface for Genotype Probabilities of Outbred Intercross
}
\description{
The package calculates and imports genotype probabilities from outbred intercross designs into an R object. This allows fast scans for main effect QTL and two way interactions. Using this interface, most functions in package \code{qtl} can be used for outbred designs.
}
\details{
\tabular{ll}{
Package: \tab qtl.outbred\cr
Type: \tab Package\cr
Version: \tab 2010.11.03\cr
Date: \tab 2010-11-03\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
Depends: \tab qtl\cr
}
}
\author{
Ronnie Nelson and Xia Shen

Maintainer: Xia Shen <xia.shen@lcb.uu.se>
}
\keyword{QTL, outbred lines, genotype probabilities}
\seealso{
\code{\link{impo.prob}}, \code{\link{calc.prob}}, \code{\link{plot_outbred}}
}
\examples{
\dontrun{
# ----- impo.prob example
file <- paste(R.home(), '/library/qtl.outbred/example_data.csv', sep = '')
data <- read.cross('csv', file = file, genotypes = c('AA', 'AB', 'BB'))
outbred.data <- impo.prob(cross.data = data, stepsize = 5, 
                          path= paste(R.home(), '/library/qtl.outbred/', sep = ''))

# ----- QTL scan using Haley-Knott regression and plotting
hk.model <- scanone(outbred.data, method = 'hk', pheno.col = 2) 
plot_outbred(hk.model, data)
}
}