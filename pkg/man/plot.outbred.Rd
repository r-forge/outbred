\name{plot.outbred}
\alias{plot.outbred}
\alias{outbred}
\title{
Plot LOD curves for outbred data
}
\description{
Plot the LOD curve including the markers for a genome scan with the single-QTL model. 
}
\usage{
\method{plot}{outbred}(x, old.data = data, chr = -1, gap = 25, ...)
}
\arguments{
  \item{x}{An object of class \code{outbred}, as output by \code{\link{outbred}}.}
  \item{old.data}{An object of class \code{cross}. See \code{\link{read.cross}} for details.}
  \item{chr}{Optional vector indicating the chromosomes to plot.}
  \item{gap}{Gap separating chromosomes (in cM).}
  \item{...}{Other arguments. See package \code{qtl} for details.}
}
\details{
This function allow the user to plot the results produced by \code{scanone}. The LOD scores for single or multiple chromosomes are plotted with the markers indicated. 

See also \code{\link{plot.scanone}} for plotting the scan results but with the probabilities on the marker positions.

}
\value{
None
}
\author{
Ronnie Nelson and Xia Shen
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{impo.prob}}, \code{\link{qtl.outbred-package}}, \code{\link{plot.scanone}}
}
\examples{
\dontrun{
# ----- impo.prob example
file <- paste(R.home(), '/library/qtl.outbred/example_data.csv', sep = '')
data <- read.cross('csv', file = file, genotypes = c('AA', 'AB', 'BB'))
path <- paste(R.home(), '/library/qtl.outbred/', sep = '')
outbred.data <- impo.prob(cross.data = data, path = path, stepsize = 5)

# ----- QTL scan using Haley-Knott regression and plotting
hk.model <- Scanone(outbred.data, method = 'hk', pheno.col = 2) 
plot(hk.model)
}
}