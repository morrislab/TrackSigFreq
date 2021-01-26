# computeTrajectories.R
# Authors: Yulia Rubanova, Cait Harrigan

## \code{generateContext} Generate a trinucleotide context from an alphabet. Note: this involves finding all three-member
## permutations of the alphabet, which can be inconveinent for large alphabets. Nucleotides are assumed to be provided as complementary pairs,
## where the first of each pair is used as the reference to build the context.
##
## @param alphabet list of pairs of characters to create combinations of as a mutation context type
## @return data.frame containing all the possible trinucleotide contextes for a mutation in the supplied alphabet
##
## @examples
## context <- TrackSig:::generateContext(c("CG", "TA"))
## dim(context) == c(96, 3)
## head(context)
##
## @name generateContext

generateContext <- function(alphabet){

  if (any(nchar(alphabet) != 2)){
    stop("Alphabet is malformed. Please provide alphabet as a list of complementary pairs")
  }

  allpha <- unlist(strsplit(alphabet, split=NULL))
  nTypes <- (length(allpha) - 1) * length(allpha)^3 * 1/2

  context <- data.frame()

  for (i in seq(1, length(allpha), by = 2)){

    midRef <- allpha[i]
    rest <- base::setdiff(allpha, midRef)
    repSize <- length(allpha)^2 - length(allpha)

    midSet <- base::cbind(rep(midRef, length.out = repSize), rep(rest, length.out=repSize),
                    paste0(sort(rep(allpha, repSize)), rep(midRef, length.out = repSize), rep(allpha, repSize)))
    context <- base::rbind(context, midSet)
  }

  stopifnot( dim(context)[1] == nTypes )

  return (context)
}

makeBinaryTable <- function(multinomial_vector)
{

  nMut <- sum(multinomial_vector)
  nMutTypes <- length(multinomial_vector)

  # one-hot encoded column vector of mutation type for each mutation
  oneHotTypes <- diag(1, nMutTypes, nMutTypes)

  sel <- rep(1:nMutTypes, times = multinomial_vector)
  binaryTable <- oneHotTypes[,sel]

  assertthat::assert_that(all(dim(binaryTable) == c(nMutTypes, nMut)),
                          msg = "Binary matrix construction failed: dimensions don't match input\n")

  return(binaryTable)

}



# fit mixture of multinomials to the vector
fitMixturesEM <- function(counts, composing_multinomials, prior=NULL)
{

  multinomial_vector <- rowSums(as.matrix(counts))


  # Number of mutations to fit
  nMut = sum(multinomial_vector)

  # Number of mutation types / categories of mutinomial
  nMutTypes <- length(multinomial_vector)

  # Number of multinomials/signatures to fit and to make mixture of
  nSig <- ncol(composing_multinomials)


  assertthat::assert_that(length(multinomial_vector) == nrow(composing_multinomials),
                          msg = "Length of data vector is not equal to nrow of matrix to fit. Did you forget to transpose the matrix?\n")

  mutType <- makeBinaryTable(multinomial_vector)

  applyMutTypeMask <- function(sigMultinom, mutType){
    return(exp(colSums(log(sigMultinom^mutType))))
  }

  # pDataGivenClass[i,n] corresponds to class/signature i and sample/mutation n
  pDataGivenClass <- t(apply(composing_multinomials, MARGIN = 2,
                             mutType = mutType, FUN = applyMutTypeMask))

  # Mixtures of multinomials. Use uniform prior unless the prior is specified
  pi <- rep(1/nSig, nSig)

  if (!is.null(prior)){

    assertthat::assert_that(length(prior) == nSig,
                            msg = sprintf("Length of prior should be equal to %s\n", nSig))
    pi <- prior
  }

  pi_diff <- Inf
  iteration <- 1

  while (pi_diff > 0.001 & iteration < 1000)
  {
    # E-step: update posterior.
    p_x <- colSums(pDataGivenClass * pi)

    # class_given_data[i,n] corresponds to class/signature i and sample/mutation n
    class_given_data <- t(t(pDataGivenClass * pi) / p_x)

    # S-step: update mixtures
    pi_new <- 1/nMut * rowSums(class_given_data)

    if (sum(pi_new > 1) != 0) {
      stop("Mixture ratio is greater than 1")
    }

    if (sum(pi_new < 0) != 0)
      stop("Mixture ratio is less than 0")

    if (sum(pi_new) > 1.5)
      stop("Sum of mixture ratios is greater than 1")

    pi_diff <- sum(abs(pi_new - pi))
    pi <- pi_new
    iteration <- iteration + 1
  }

  return(pi)
}


distributeBinCounts <- function(binCounts){
  # split the counts of a changepoint bin and add them to the chunks that surround it

  # distribute evenly what can be
  leftChunk <- rightChunk <- floor(binCounts / 2)

  # sample masks to assign remaining counts
  binCounts <- binCounts - (2 * floor(binCounts/2))
  leftMask <- sample( which(binCounts > 0), (length(which(binCounts > 0)) / 2) )
  rightMask <- base::setdiff(which(binCounts >0), leftMask)

  # distribute remaining counts
  leftChunk[leftMask] <- leftChunk[leftMask] + binCounts[leftMask]
  rightChunk[rightMask] <- rightChunk[rightMask] + binCounts[rightMask]

  return(list(leftAdd = leftChunk, rightAdd = rightChunk))

}


# fit mixture of mutinomials in each time slice specified by change_points
fitMixturesInTimeline <- function(data, changepoints, alex.t, split_data_at_change_point = T)
{

  # cast to matrix if possible
  data <- as.matrix(data)

  # allocate fitted values matrix
  fitted_values <- matrix(NA, ncol=ncol(data), nrow=ncol(alex.t),
                          dimnames = list(colnames(alex.t), colnames(data)))

  sumSlice <- function(slice, data){return(rowSums(subset(data, select = slice)))}
  repChunk <- function(chunkFit, times, nSig){return(matrix(rep(chunkFit, times), nrow = nSig))}

  # if no changepoints, use all data
  if (length(changepoints) == 0) {

    fitted_for_time_slice <- fitMixturesEM(rowSums(data), alex.t)
    fitted_values <- matrix(rep(fitted_for_time_slice, ncol(fitted_values)), nrow=nrow(fitted_values),
                            dimnames = list(colnames(alex.t), colnames(data)))
    return(fitted_values)

  }

  # ensure changepoints are valid
  assertthat::assert_that((1 %in% changepoints) == FALSE, msg = "Impossible changepoint, cannot segment before first timepoint\n")
  assertthat::assert_that((dim(data)[2] %in% changepoints) == FALSE, msg = "Impossible changepoint, cannot segment after last timepoint\n")
  changepoints <- sort(changepoints)

  # TODO: distributing counts between bins with two immediately adjacent cahngepoints
  # is undefined behaviour. Minseglen >2 should avoid this, but may not be guarenteed.
  if (any((c(1, changepoints + 1) - c(changepoints - 1, dim(data)[2])) == 1 )){
    split_data_at_change_point <- F
  }

  # if changepoints, get changepoints as data indices
  if (split_data_at_change_point){

    slices <- mapply(c(1, changepoints + 1), c(changepoints - 1, dim(data)[2]), FUN = `:`, SIMPLIFY = F)
    chunkSums <- lapply(slices, data, FUN = sumSlice)

    # split change point bins over chunks
    for (cp_i in 1:length(changepoints)){
      leftAdd <- rightAdd <- NULL
      list[leftAdd, rightAdd] <- distributeBinCounts(data[,changepoints[cp_i]])
      chunkSums[[cp_i]] <- chunkSums[[cp_i]] + leftAdd
      chunkSums[[(cp_i + 1)]] <- chunkSums[[(cp_i + 1)]] + rightAdd
    }

    # all counts should be present
    assertthat::assert_that(all(rowSums(data) == rowSums(do.call(cbind,chunkSums))),
                            msg = "Timepoints lost in chunking\n")


  } else {
    slices <- mapply(c(1, changepoints + 1), c(changepoints, dim(data)[2]), FUN = `:`, SIMPLIFY = F)
    chunkSums <- lapply(slices, data, FUN = sumSlice)

    # all counts should be present
    assertthat::assert_that(all(base::rowSums(data) == rowSums(do.call(cbind,chunkSums))),
                            msg = "Timepoints lost in chunking\n")
  }


  chunkFits <- lapply(chunkSums, composing_multinomials = alex.t, FUN = fitMixturesEM)
  chunkFits <- mapply(chunkFits, times = c(changepoints, dim(data)[2]) - c(0, changepoints),
                      nSig = dim(alex.t)[2], FUN = repChunk, SIMPLIFY = F)

  fitted_values <- do.call(cbind, chunkFits)
  dimnames(fitted_values) <- list(colnames(alex.t), colnames(data))

  return(fitted_values)
}



mixtureLL <- function(counts, composing_multinomials, mixtures, ...) {
  # replaces log_likelihood_mixture_multinomials
  multinomial_vector <- rowSums(as.matrix(counts))

  mutation_binary_table <- makeBinaryTable(multinomial_vector)

  # mutation_probabilities_under_signature_mixture[i,n] corresponds to class/signature i and sample/mutation n
  mutation_probabilities_under_signature_mixture <- matrix(0, nrow=ncol(composing_multinomials), ncol=ncol(mutation_binary_table))

  for (sig in 1:ncol(composing_multinomials)) {
    mutation_probabilities_under_signature_mixture[sig,] <- apply(composing_multinomials[,sig]^mutation_binary_table,2,prod)
  }

  mutation_probabilities_under_mixture <-  log(t(mutation_probabilities_under_signature_mixture) %*% as.matrix(mixtures))
  stopifnot(length(mutation_probabilities_under_mixture) == sum(multinomial_vector))

  return(sum(mutation_probabilities_under_mixture))
}

# beta likelihood maximization
betaLL <- function(qis, ...){

  #qis are the VAFs for the subproblem

  n <- length(qis)

  #assertthat::assert_that(length(qis) == length(vis), length(qis) == length(ris), msg = "problem subsetting is not good!")

  alpha <- sum(qis) + 1
  beta <- sum(1-qis) + 1

  LL <- lbeta(alpha, beta) + log(stats::pbeta(max(qis), alpha, beta) - stats::pbeta(min(qis), alpha, beta))

  #print(c(alpha, beta, LL))

  return(LL)

}

sumBetaMixtureLL <- function(qis, counts,
                                composing_multinomials, mixtures, ...){


  score <- sum(
               mixtureLL(counts, composing_multinomials, mixtures),
               betaLL(qis)
              )

  return(score)

}


parseScoreMethod <- function(scoreMethod){
  # return the penalty and score function to use when computing partitions

  #assertthat::assert_that(scoreMethod %in% c("SigFreq", "Signature", "Frequency"),
  #msg = "scoreMethod should be one of \"SigFreq\", \"Signature\", \"Frequency\". \n Please see documentation for more information on selecting a scoreMethod)")

  if(scoreMethod == "SigFreq"){
    return(list(penalty = expression(-log(0.1) + (n_sigs + 1) * log(n_mut)),
                score_fxn = sumBetaMixtureLL))
  }

  if(scoreMethod == "Signature"){
    return(list(penalty = expression((n_sigs - 1) * log(n_mut)),
                score_fxn = mixtureLL))
  }

  if(scoreMethod == "Frequency"){
    return(list(penalty = expression((n_sigs + 2) * log(n_mut)),
                score_fxn = betaLL))
  }

  stop("scoreMethod should be one of \"SigFreq\", \"Signature\", \"Frequency\". \n Please see documentation for more information on selecting a scoreMethod\n)")

}

getActualMinSegLen <- function(desiredMinSegLen, binSize, n_mut){
  # return the minimum segment length (in bins) to use

  assertthat::assert_that(is.numeric(desiredMinSegLen) & (0 < desiredMinSegLen),
                          msg = "desiredMinSegLen should be an integer greater than 0\n")

  # for high resolution segment scoring, use at least 400 mutations per segment.
  if( (desiredMinSegLen == 1) & (n_mut > 400)){
    return ( ceiling(400/binSize) )
  }

  # for accurate segment scoring, reqire at least 100 mutations per segment.
  actualMinSegLen <- max(as.integer(desiredMinSegLen), ceiling(100/binSize) )

  if (actualMinSegLen != desiredMinSegLen){
    warning(sprintf("Could not use desiredMinSegLen, too few mutations for accurate segment scoring. minSegLen set to: %s\n", actualMinSegLen))
  }

  if (n_mut < binSize * actualMinSegLen){
    warning("Not enough total mutations to segment at selected binSize.\n")
  }

  return(actualMinSegLen)
}

# Find optimal changepoint and mixtures using PELT method.
# if desiredMinSegLen is NULL, the value will be selected by default based off binSize to try to give good performance

getChangepointsPELT <- function(vcaf, countsPerBin, referenceSignatures, scoreMethod, binSize = 100, desiredMinSegLen = NULL)

{

  minSegLen <- getActualMinSegLen(desiredMinSegLen, binSize, dim(vcaf)[1])
  score_matrix <- scorePartitionsPELT(countsPerBin, referenceSignatures, vcaf, scoreMethod, binSize, minSegLen)

  #print(score_matrix[1:15, 1:15])

  changepoints <- recoverChangepoints(score_matrix)
  mixtures <- fitMixturesInTimeline(countsPerBin, changepoints, referenceSignatures)

  # mixtures should also contain binned phi
  binned_phis <- stats::aggregate(vcaf$phi, by = list(vcaf$bin), FUN = mean)$x
  colnames(mixtures) <- binned_phis

  return(list(mixtures = mixtures, changepoints = changepoints))
}

# Calculate penalized BIC score for all partitions using PELT method.
scorePartitionsPELT <- function(countsPerBin, referenceSignatures, vcaf, scoreMethod, binSize, minSegLen)
{
  n_bins <- dim(countsPerBin)[2]
  n_sigs <- dim(referenceSignatures)[2]
  n_mut <- dim(vcaf)[1]


  penalty <- score_fxn <- NULL
  list[penalty, score_fxn] <- parseScoreMethod(scoreMethod)
  penalty <- eval(penalty)

  # Store score for all partitions of all sub-problems
  # Rows are length of sub-problem. Columns correspond to last changepoint
  sp_scores <- matrix(nrow=n_bins, ncol=n_bins)

  max_sp_scores <- numeric(n_bins)
  prune_set <- c()

  # Replace print msg with progress bar
  pb <- progress::progress_bar$new(format = "Scoring subpartitions: [:bar] :percent",
                         total = n_bins, clear = FALSE, width = 60)

  # Score all subproblems of length sp_len using last_cp as last changepoint
  for (sp_len in 1:n_bins)
  {
    valid_cps <- setdiff(0:(sp_len - 1), prune_set)
    pb$tick()

    for (last_cp in valid_cps){

      # check segment length
      if (sp_len - last_cp < minSegLen){
        sp_scores[sp_len, last_cp + 1] <- -Inf
        next
      }

      # score segment
      sp_slice <- c((last_cp + 1), sp_len)
      r_seg_qis <- vcaf$qi[vcaf$bin %in% (sp_slice[1] : sp_slice[2])]
      r_seg_counts <- (countsPerBin[, sp_slice[1] : sp_slice[2]])
      r_seg_mix <- fitMixturesEM(r_seg_counts, referenceSignatures)


      r_seg_score <- 2 * score_fxn(counts = r_seg_counts, composing_multinomials = referenceSignatures,
                                   mixtures = r_seg_mix, qis = r_seg_qis, sp_len = sp_len)
      l_seg_score <- ifelse(last_cp == 0, penalty, max_sp_scores[last_cp])

      sp_scores[sp_len, last_cp + 1] <- l_seg_score + r_seg_score - penalty

    }

    max_sp_scores[sp_len] <- max(sp_scores[sp_len, ][!is.na(sp_scores[sp_len, ])])

    # Evaluate all changepoints for pruning condition
    for (cp in valid_cps){

      # check segment length
      if (sp_len - cp < minSegLen){
        next
      }

      # prune
      if (sp_scores[sp_len, cp + 1] + penalty < max_sp_scores[sp_len]){
        prune_set <- c(prune_set, cp)
      }
    }

  }

  return(sp_scores)
}

# Recover optimal changepoints by from subproblem matrix
recoverChangepoints <- function(sp_score_matrix)
{
  changepoints <- c()

  continue <- TRUE
  current <- dim(sp_score_matrix)[1]
  while (continue){

    prev <- which.max(sp_score_matrix[current, ])

    if (prev - 1 <= 1){
      continue <- FALSE
    }else{
      changepoints <- c(prev - 1, changepoints)
    }

    current <- prev - 1
  }

  return(changepoints)
}


# [END]
