\name{Scanone}
\alias{Scanone}

\title{Genome scan with a single QTL model} 

\description{
  Genome scan with a single QTL model, with possible allowance for
  covariates, using any of several possible models for the phenotype and
  any of several possible numerical methods.
}

\usage{
Scanone(cross, chr, pheno.col=1, model=c("normal","binary","2part","np"),
        method=c("em","imp","hk","ehk","mr","mr-imp","mr-argmax"),
        addcovar=NULL, intcovar=NULL, weights=NULL,
        use=c("all.obs", "complete.obs"), upper=FALSE,
        ties.random=FALSE, start=NULL, maxit=4000,
        tol=1e-4, n.perm, perm.Xsp=FALSE, perm.strata=NULL, verbose,
        batchsize=250, n.cluster=1, ind.noqtl)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{chr}{Optional vector indicating the chromosomes for which LOD
   scores should be calculated.  This should be a vector of character
   strings referring to chromosomes by name; numeric values are
   converted to strings.  Refer to chromosomes with a preceding \code{-}
   to have all chromosomes but those considered.  A logical (TRUE/FALSE)
   vector may also be used.}
 \item{pheno.col}{Column number in the phenotype matrix which should be
   used as the phenotype.  This can be a vector of integers; for methods
   \code{"hk"} and \code{"imp"} this can be considerably faster than doing
   them one at a time.  One may also give a character strings matching
   the phenotype names.  Finally, one may give a numeric vector of
   phenotypes, in which case it must have the length equal to the number
   of individuals in the cross, and there must be either non-integers or
   values < 1 or > no. phenotypes; this last case may be useful for studying
   transformations.}
 \item{model}{The phenotype model: the usual normal model, a model for
   binary traits, a two-part model or non-parametric analysis}
 \item{method}{Indicates whether to use the EM algorithm, 
   imputation, Haley-Knott regression, the extended Haley-Knott method,
   or marker regression.  Not all methods are available for all models.
   Marker regression is performed either by dropping individuals with
   missing genotypes (\code{"mr"}), or by first filling in missing data
   using a single imputation (\code{"mr-imp"}) or by the Viterbi
   algorithm (\code{"mr-argmax"}).}
 \item{addcovar}{Additive covariates;
   allowed only for the normal and binary models.}
 \item{intcovar}{Interactive covariates (interact with QTL genotype);
   allowed only for the normal and binary models.}
 \item{weights}{Optional weights of individuals.  Should be either NULL
   or a vector of length n.ind containing positive weights.  Used only
   in the case \code{model="normal"}.}
 \item{use}{In the case that multiple phenotypes are selected to be
   scanned, this argument indicates whether to use all individuals, 
   including those missing some phenotypes, or just those individuals
   that have data on all selected phenotypes.}
 \item{upper}{Used only for the two-part model; if true, the
   "undefined" phenotype is the maximum observed phenotype; otherwise,
   it is the smallest observed phenotype.} 
 \item{ties.random}{Used only for the non-parametric "model"; if TRUE,
   ties in the phenotypes are ranked at random.  If FALSE, average ranks
   are used and a corrected LOD score is calculated.}
 \item{start}{Used only for the EM algorithm with the normal model and
   no covariates.  If \code{NULL}, use the usual starting values; if
   length 1, use random initial weights for EM; otherwise, this should
   be a vector of length n+1 (where n is the number of possible
   genotypes for the cross), giving the initial values for EM.}
 \item{maxit}{Maximum number of iterations for methods \code{"em"} and
   \code{"ehk"}.}
 \item{tol}{Tolerance value for determining convergence for methods
   \code{"em"} and \code{"ehk"}.} 
 \item{n.perm}{If specified, a permutation test is performed rather than
   an analysis of the observed data.  This argument defines the number
   of permutation replicates.}
 \item{perm.Xsp}{If \code{n.perm} > 0, so that a permutation test will
   be performed, this indicates whether separate permutations should be
   performed for the autosomes and the X chromosome, in order to get an
   X-chromosome-specific LOD threshold.  In this case, additional
   permutations are performed for the X chromosome.}
 \item{perm.strata}{If \code{n.perm} > 0, this may be used to perform a
   stratified permutation test.  This should be a vector with the same
   number of individuals as in the cross data.  Unique values indicate
   the individual strata, and permutations will be performed within the
   strata.}
 \item{verbose}{In the case \code{n.perm} is specified, display
   information about the progress of the permutation tests.}
 \item{batchsize}{The number of phenotypes (or permutations) to be run
   as a batch; used only for methods \code{"hk"} and \code{"imp"}.}
 \item{n.cluster}{If the package \code{snow} is available and
   \code{n.perm} > 0, permutations are run in parallel using this number
   of nodes.}
 \item{ind.noqtl}{Indicates individuals who should not be allowed a QTL
   effect (used rarely, if at all); this is a logical vector of same
   length as there are individuals in the cross.}
}


\details{
  See \code{"scanone"} in \code{"qtl"} package for details.
}

\value{
  If \code{n.perm} is missing, the function returns a data.frame whose
  first two columns contain the chromosome IDs and cM positions.
  Subsequent columns contain the LOD scores for each phenotype.
  In the case of the two-part model, there are three LOD score columns
  for each phenotype: LOD(\eqn{p,\mu}{p,mu}), LOD(\eqn{p}) and
  LOD(\eqn{\mu}{mu}).   The result is given class \code{"scanone"} and
  has attributes  \code{"model"}, \code{"method"}, \code{"df"} and
  \code{"type"} (the latter is the type of cross analyzed). 

  If \code{n.perm} is specified, the function returns the results of a
  permutation test and the output has class \code{"scanoneperm"}.  If
  \code{perm.Xsp=FALSE}, the function returns a matrix with
  \code{n.perm} rows, each row containing the genome-wide 
  maximum LOD score for each of the phenotypes.  In the case of the
  two-part model, there are three columns for each phenotype,
  corresponding to the three different LOD scores. If
  \code{perm.Xsp=TRUE}, the result contains separate permutation results
  for the autosomes and the X chromosome respectively, and an attribute
  indicates the lengths of the chromosomes and an indicator of which
  chromosome is X.
}

