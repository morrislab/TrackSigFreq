# TrackSig.R
# Defines main functions for user to interact with package TrackSig.
# Author: Cait Harrigan


#' Get a list of signatures active in a sample
#'
#' @description
#' \code{detectActiveSignatures} determines what signatures are active above an
#' activity threshold, from a list of signature definintions. To do this, a
#' multinomial mixture model of signature activities is fit to the collection of
#' all mutations in the sample.
#'
#' @param vcfFile path to variant calling format (vcf) file
#' @param cnaFile path to copy number abberation (cna) file. If not provided,
#'   all copy numbers default to 2.
#' @param purity number between 0 and 1 of the percentage of cells in the sample
#'   that are cancerous
#' @param threshold minimum activity level that signature must have to be
#'   detected
#' @param prior prior on the likelihood of observing a given signature (must
#'   match signatures present in referenceSignatures)
#' @param binSize number of mutations per bin (default 100)
#' @param referenceSignatures dataframe containing definitions of mutational
#'   signatures.
#' @param refGenome BSgenome to use as reference
#'
#' @return Names of signatures active in sample.
#' @export
detectActiveSignatures <- function(vcfFile, cnaFile = NULL, purity = 1,
                                   threshold = 0.05, prior = NULL, binSize = 100,
                                   referenceSignatures = alex_merged,
                                   refGenome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19){

  # return list of active signatures in sample, whether by matching per-cancer-type to provided data,
  # or fitting all counts by EM. If not using this function, must provide active signatures per sample

  context = generateContext(c("CG", "TA"))

  binCounts <- vcfToCounts(vcfFile = vcfFile, cnaFile = cnaFile,
                           purity = purity, binSize = binSize,
                           context = context, refGenome = refGenome)[[2]]

  counts <-  rowSums(binCounts)
  mixtures <- fitMixturesEM(counts, referenceSignatures, prior=prior)
  mixtures <- mixtures[mixtures >= threshold]
  mixtures <- sort(mixtures, decreasing = T)

  return(names(mixtures))
}


#' Determine an evolutionary trajectory.
#'
#'
#' @description
#' \code{TrackSig} will take an input VCF file and corresponding annotation, and determine an eveolutionary trajectory for the sample, based on changepoints found using the PELT segmentation algorithm.
#'
#' @param vcfFile path to variant calling format (vcf) file
#' @param activeInSample list of signatures names to fit the exposures of. All listed signatures must be present in the referenceSignatures dataframe.
#' @param cnaFile path to copy number abberation (cna) file. If not provided, all copy numbers default to 2.
#' @param sampleID name to call sample. If none provided, name will be automatically drawn from the provided vcf file name.
#' @param referenceSignatures dataframe containing definitions of mutational signatures.
#' @param purity number between 0 and 1 of the percentage of cells in the sample that are cancerous
#' @param scoreMethod string to indicate what scoring method to apply when finding changepoints. Options are "Signature", "SigFreq" and "Frequency"
#' @param binSize number of mutations per bin (default 100)
#' @param nCutoff maximum number of total mutations to consider (samples with more than nCutoff muations will be down-sampled)
#' @param desiredMinSegLen minimum number of mutations to include in a PELT segment (the desiredMinSegLen will be overridden if there are too few for accurate scoring)
#' @param refGenome BSgenome to use as reference
#'
#' @export

TrackSig <- function(vcfFile,
                     activeInSample,
                     cnaFile = NULL,
                     purity = 1,
                     sampleID = NULL,
                     referenceSignatures = TrackSig:::alex_merged,
                     scoreMethod = "SigFreq",
                     binSize = 100,
                     nCutoff = 10000,
                     desiredMinSegLen = 1,
                     refGenome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19) {

  # input checking

  assertthat::assert_that(grepl(".vcf$", vcfFile) | grepl(".txt$", vcfFile), msg = "Unsupported VCF file extension. Expected file type .vcf or .txt\n")

  if(missing(activeInSample)){
    assertthat::assert_that(scoreMethod == "Frequency", msg = "When scoreMethod is not equal to \"Frequency\", activeInSample must be provided.")
  }else{
    assertthat::assert_that(all(activeInSample %in% colnames(referenceSignatures)))
  }

  assertthat::assert_that(is.numeric(purity) & (0 < purity) & (purity <= 1),
                          msg = "Purity should be a proportion between 0 and 1\n")


  # TODO: activeSignatures %in% colnames(referenceSignatures) must be TRUE
  # TODO: length(activeInSample) >1 should be true, else no mixture to fit
  # TODO: binSize has to make sense; positive, not larger than nMut, maybe throw warning if it's some ratio too large for low-resolution.
  # TODO: generateContext and mut types in referenceSignatures should make sense together.

  # take sampleID from file name if not provided
  if (is.null(sampleID)){

    if (grepl(".txt$", vcfFile)){
      sampleID <- strsplit( unlist(strsplit(vcfFile, "/"))[ length( strsplit(vcfFile, "/")[[1]] ) ] , ".txt")[[1]]
    }

    if (grepl(".vcf$", vcfFile)){
      sampleID <- strsplit( unlist(strsplit(vcfFile, "/"))[ length( strsplit(vcfFile, "/")[[1]] ) ] , ".vcf")[[1]]
    }else{stop("Failed setting sampleID. Please check input vcf file.")}

  }

  # TODO: get context from supplied referenceSignatures
  context <- generateContext(c("CG", "TA"))

  # TODO: other parameters non-default options
  vcaf <- countsPerBin <- NULL
  list[vcaf, countsPerBin] <- vcfToCounts(vcfFile = vcfFile, cnaFile = cnaFile,
                                          purity = purity, binSize = binSize,
                                          nCutoff = nCutoff, context = context,
                                          refGenome = refGenome)


  assertthat::assert_that(all(rownames(countsPerBin) %in%
                                rownames(referenceSignatures)),
                          msg = "Mutation type counts failed.")

  countsPerBin <- countsPerBin[rownames(referenceSignatures),,drop = FALSE]

  # subset referenceSignatures with activeInSample
  referenceSignatures <- referenceSignatures[activeInSample]

  if ( any(rowSums(countsPerBin)[rowSums(referenceSignatures) == 0] != 0) ) {
    print(sprintf("Error in sample %s: Some mutation types have probability 0 under the model, but their count is non-zero. This count vector is impossible under the model.", sampleID))
  }

  # compute results
  mixtures <- changepoints <- NULL
  list[mixtures, changepoints] <- getChangepointsPELT(vcaf = vcaf,
                                                      countsPerBin = countsPerBin,
                                                      referenceSignatures = referenceSignatures,
                                                      scoreMethod = scoreMethod,
                                                      binSize = binSize,
                                                      desiredMinSegLen = desiredMinSegLen)
  # assign mutations to clusters
  if(is.null(changepoints)){
    vcaf$clust = 1
  }else{
    clustIdx = rep(1:(length(changepoints) + 1),
                   times = c(changepoints, max(vcaf$bin)) - c(0, changepoints))
    vcaf$clust = clustIdx[vcaf$bin]
  }

  return (list(mixtures = mixtures, changepoints = changepoints, sampleID = sampleID, binData = vcaf))
}


# list unpacker util: used internally in package TrackSig
# source: https://stat.ethz.ch/pipermail/r-help/2004-June/053343.html

list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}

# TODO: function recommendBinSize()


# TODO: see http://adv-r.had.co.nz/S3.html for best practices
# constructor function for tracksig results. Returned by TrackSig when class=T
TS.trajectory <- function(sampleID = NULL, scoreMethod = NULL,
                          mixtures = NULL, changepoints = NULL,
                          binData = NULL, binSize = NULL){

  return(structure(list(), class = "TS.trajectory",
                      sampleID = sampleID, scoreMethod = scoreMethod,
                      mixtures = mixtures, changepoints = changepoints,
                      binData = binData, binSize = binSize
                    )
        )
}

is.TS.trajectory <- function(x) inherits(x, "TS.trajectory")

# [END]
