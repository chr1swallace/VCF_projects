\name{snpSummary}
\alias{snpSummary}


\title{Summarize SNPs in a CollapsedVCF object}

\description{Calculate genotype and allele counts, and Hardy Weinburg
  chisquare statistics for SNPs or any variant with exactly two alleles
  in a CollapsedVCF object}

\usage{
  \S4method{snpSummary}{CollapsedVCF}(x, ...)
}

\arguments{
  \item{x}{A CollapsedVCF object}.
}

\details{Genotype counts, allele counts and Hardy Weinberg equilibrium
  (HWE) statistics are calculated.  HWE has been established as a useful
  quality filter on genotype data.  Hardy Weinberg equilibrium should be
  attained in a single generation of random mating, and departures,
  indicated by small p values, are almost invariably indicative of a
  problem with genotype calls.}

\value{
  The object returned is a list with four elements:
  \describe{
    \item{genotype.counts:}{
      Three column matrix with the counts for each genotype
    }
    \item{allele.counts:}{
      Two column matrix with the counts for each allele
    }
    \item{hwe:}{
      Two column matrix giving a Z score and p value for departure from
      a null hypothesis of Hardy Weinberg equilibrium.  NB, this is only
      calculated when there are exactly two alleles and both are observed.
    }
    \item{criteria:}{
      Two column data.frame.  Summary stats are only calculated variants
      with \code{diploid==TRUE} and \code{n.alt.alleles<=1} (ie at most
      one alternative allele).
    }
}


\author{
  Chris Wallace <cew54@cam.ac.uk>
}

\examples{
  fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
  x <- readVcf(fl, "hg19")
  snpSummary(x)

}

\keyword{manip}
