# TrackSig.R
# Defines main functions for user to interact with package TrackSig.
# Author: Cait Harrigan


detectActiveSignatures <- function(sample, referenceSignatures){

  # return list of active signatures in sample, whether by matching per-cancer-type to provided data,
  # or fitting all counts by EM. If not using this function, must provide active signatures per sample

  NULL
}

#' \code{TrackSig} Take an input vcf file and annotation and generate the counts data.
#' Create all plotting output that compute_signatures_for_all_examples does.
#'
#' @rdname TrackSig
#' @name TrackSig
#'
#' @param vcfFile path to variant calling format (vcf) file
#' @param cnaFile path to copy number abberation (cna) file
#' @param purityFile path to sample purity file
#' @param sampleID name to call sample. If none provided, name will be automatically drawn from the provided vcf file name.
#' @param saveIntermediate boolean whether to save intermediate results (mutation types)
#'
#'
#' activeInSample is list used to subset referenceSignatures
#'
#' @export

TrackSig <- function(vcfFile,
                     activeInSample,
                     cnaFile = NULL,
                     purity = NULL,
                     sampleID = NULL,
                     referenceSignatures = alex,
                     scoreMethod = "SigFreq",
                     binSize = 100,
                     desiredMinSegLen = NULL,
                     refGenome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19) {

  # input checking

  assertthat::assert_that(grepl(".vcf$", vcfFile) | grepl(".txt$", vcfFile), msg = "Unsupported VCF file extension. Expected file type .vcf or .txt")

  assertthat::assert_that(scoreMethod %in% c("SigFreq", "Signature", "Frequency"),
  msg = "scoreMethod should be one of \"SigFreq\", \"Signature\", \"Frequency\". \n Please see documentation for more information on selecting a scoreMethod)")

  # TODO: activeSignatures %in% rownames(referenceSignatures) must be TRUE
  # TODO: length(activeInSample) >1 should be true, else no mixture to fit
  # TODO: binSize has to make sense; positive, not larger than nMut, maybe throw warning if it's some ratio too large for low-resolution.

  # take sampleID from file name if not provided
  if (is.null(sampleID)){

    if (grepl(".txt$", vcfFile)){
      sampleID <- strsplit( unlist(strsplit(vcfFile, "/"))[ length( strsplit(vcfFile, "/")[[1]] ) ] , ".txt")[[1]]
    }

    if (grepl(".vcf$", vcfFile)){
      sampleID <- strsplit( unlist(strsplit(vcfFile, "/"))[ length( strsplit(vcfFile, "/")[[1]] ) ] , ".vcf")[[1]]
    }else{stop("Failed setting sampleID. Please check input vcf file.")}

  }

  # TODO: other parameters non-default options
  list[vcaf, countsPerBin] <- vcfToCounts(vcfFile, cnaFile, purity, binSize, context = generateContext(c("CG", "TA")), refGenome)

  assertthat::assert_that(all(rownames(countsPerBin) == rownames(referenceSignatures)), msg = "Mutation type counts failed.")

  # subset referenceSignatures with activeInSample
  referenceSignatures <- referenceSignatures[activeInSample]

  if ( any(rowSums(countsPerBin)[rowSums(referenceSignatures) == 0] != 0) ) {
    print(sprintf("Error in sample %s: Some mutation types have probability 0 under the model, but their count is non-zero. This count vector is impossible under the model.", sampleID))
  }

  # compute results
  list[changepoints, mixtures] <- getChangepointsPELT(countsPerBin, referenceSignatures, vcaf, scoreMethod, binSize, desiredMinSegLen)

  plot <- NULL

  # side effect: plot
  tryCatch({

            binned_phis <- aggregate(vcaf$phi, by = list(vcaf$bin), FUN = mean)$x

            plot <- ( plotTrajectory(mixtures * 100, phis = binned_phis, changepoints, linearX = T, anmac = T)
                      + ggtitle(paste0(sampleID, " Signature Trajectory"))
                    )

            print(plot)

           },
           warning = function(w){w},
           error = function(e){print("Error: failed to plot signature trajectory")}
          )


  return (list(mixtures = mixtures, changepoints = changepoints, plot = plot, vcaf = vcaf))
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



# [END]
