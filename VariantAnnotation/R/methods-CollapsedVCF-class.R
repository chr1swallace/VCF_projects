### =========================================================================
### CollapsedVCF class methods
### =========================================================================

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###
### See VCF() on VCF-class page

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and Setters
###

### alt() setter
setReplaceMethod("alt", c("CollapsedVCF", "DNAStringSetList"),
    function(x, value)
{
    if (length(value) != length(rowData(x)))
        stop("length(value) must equal length(rowData(x))")
    slot(x, "fixed")$ALT <- value
    x
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### expand
###

setMethod("expand", "CollapsedVCF",
    function(x, ...)
    {
        elt <- elementLengths(alt(x))
        if (all(elt == 1L) || is(alt(x), "CharacterList")) {
            fxd <- mcols(fixed(x))[-1]
            fxd$ALT <- unlist(alt(x), use.names=FALSE)
            return(VCF(rowData=rowData(x), colData=colData(x),
                       exptData=exptData(x), fixed=fxd,
                       info=mcols(info(x))[-1], geno=geno(x),
                       ..., collapsed=FALSE))
        }
        idx <- rep.int(seq_len(nrow(x)), elt)

        ## info and fixed:
        iexp <- .expandInfo(x, elt, idx)
        fexp <- mcols(fixed(x))[idx, -1]
        fexp$ALT <- unlist(alt(x), use.names=FALSE)

        ## geno:
        gexp <- endoapply(geno(x), function(i) {
                    if (is(i, "matrix"))
                        i[idx, ]
                    else
                        i[idx, , ]
                 })

        ## rowData:
        rdexp <- rowData(x)[idx, ]

        ## exptData, colData untouched
        VCF(rowData=rdexp, colData=colData(x), exptData=exptData(x),
            fixed=fexp, info=iexp, geno=gexp, ..., collapsed=FALSE)
    }
)

.expandInfo <- function(x, elt, idx)
{
    hd <- exptData(x)$header
    icol <- mcols(info(x))[-1]
    ## FIXME: handle inconsistent elementLengths in readVcf()?
    inms <- rownames(info(hd))[info(hd)$Number == "A"]
    ielt <- sapply(inms, function(i)
                elt - elementLengths(mcols(info(x))[[i]]))
    if (any(colSums(ielt) == 0L))
        res <- expand(icol, inms[colSums(ielt) == 0L], TRUE)
    if (any(colSums(ielt) != 0)) {
        res <- icol[idx, ]
        for (i in inms[colSums(ielt) != 0L]) {
            irep <- rep(seq_len(nrow(icol)), ielt[colnames(ielt) == i] + 1L)
            res[[i]] <- unlist(icol[[i]][irep], use.names=FALSE)
        }
    }
    res
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show
###

setMethod(show, "CollapsedVCF",
    function(object)
{
    .showVCFSubclass(object)
})

### --------------------------------------------------------------------------
### snpSummary
###

## snpSummary
setMethod("snpSummary", "CollapsedVCF", function(x, ...) {
  ## alleles
  g.ref <- ref(x)
  g.alt <- alt(x)

  ## called genotypes - ignore phasing and ignore, for now, probabilistic calls
  g<-sub("/|\\|","",geno(x)$GT)

  ## identify those that meet criteria
  diploid<-apply(nchar(g)==2,1,all) # diploid
  n.alt <- sapply(g.alt,length) # exactly one alt allele
  use<-n.alt<=1 & diploid

  ## genotype counts of allele 1
  g.ct<-matrix(NA,nrow(g),3)
  dimnames(g.ct) <- list(Variant=rownames(g),Genotype=c("g00","g01","g11"))
  g.ct[use,] <- cbind(rowSums(g[use,]=="00"),
                      rowSums(g[use,]=="10"|g[use,]=="01"),
                      rowSums(g[use,]=="11"))
  ## matrix explaining whether variants got dropped
  g.use <- data.frame(diploid=diploid,
                      n.alt.alleles=n.alt)

  ## allele counts
  a.ct <- g.ct %*% matrix(c(2,1,0,0,1,2),nrow=3)
  a.n=rowSums(a.ct)
  a.fq<-a.ct/a.n
  colnames(a.fq) <- sprintf("a%s",0:1)

  ## hwe
  g.n=rowSums(g.ct)
  hwe.expt <- a.n * cbind(a.fq[,1]^2/2, a.fq[,1]*a.fq[,2], a.fq[,2]^2/2)
  hwe.chisq<-rowSums((g.ct-hwe.expt)^2/hwe.expt)
  Z <- rep(NA,nrow(g))
  use1 <- diploid & n.alt==1 & apply(a.ct>0,1,all)
  Z[use1] <- sign(g.ct[use1,2] - hwe.expt[use1,2]*g.n[use1]) * sqrt(hwe.chisq[use1])

  ## a list may not be the best thing to return
  ## alternatives: a matrix with lots of columns or a new class
  ## with methods to show all/some elements
  return(list(genotype.counts=g.ct,
              allele.freq=a.fq,
              hwe=cbind(Z=Z,
                p.value=pchisq(Z^2,1,lower=FALSE)),
              criteria=g.use))

})
