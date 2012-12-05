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
  The object returned is a data.frame with nine columns
  \describe{
    \item{genotype.counts.*:}{
      Three columns with the counts for each genotype, g00, g01, g11
    }
    \item{allele.freq.*:}{
      Two columns with the frequencies for each allele, a1, a2
    }
    \item{hwe.*:}{
      Two columns giving a Z score (hwe.Z) and p value (hwe.p.value) for departure from
      a null hypothesis of Hardy Weinberg equilibrium.  NB, this is only
      calculated when there are exactly two alleles and both are observed.
    }
    \item{criteria.*:}{
      Two columns.  Summary stats are only calculated variants
      with \code{criteria.diploid==TRUE} and \code{criteria.n.alt.alleles<=1} (ie at most
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
