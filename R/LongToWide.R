#' Function to format the dataset (transposing the datset from Long to Wide structure)
#' @export
#' 
LongToWide <- function(subject.id, time.sequence, outcome, verbose = FALSE) {
  # Transpose data from Long to Wide format
  #
  # Args:
  #   subject.id    -->   a column of subject identifiers
  #   time.sequence -->   a column of time points
  #   outcome       -->   a column of measures to be transposed
  #
  # Optional:
  #   verbose       -->   TRUE or FALSE (default = FALSE)
  #                       prints output to screen. Usefull for ensuring accuracy
  #
  # Returns:
  #   mat.wide      -->  a matrix of data in short format
  #                     columns = timeseq, rows = outcome for each subject
  #
  #
  ##

  if (missing(verbose)) {
    verbose == FALSE
  } else if (verbose != TRUE && verbose != FALSE) {
    stop("Input 'verbose' should be TRUE or FALSE")
  }

  testM <- cbind(subject.id, time.sequence, outcome)
  if (is.unsorted(testM[, 1])) {
    warning("Data was not sorted by subject.id.\n  Sorting now.")
    testM <- testM[order(testM[,1]), ]
  }
  N.sub <- unique(testM[, 1])
  N.time <- sort(unique(testM[, 2]))
  mat.wide <- matrix(data = NA, nrow = length(N.sub), ncol = length(N.time))


  counter.sub <- 1
  counter.time <- 1
  while (counter.time <= dim(testM)[1]) {
    sub.id <- testM[counter.time]
    while ((sub.id == N.sub[counter.sub]) & (counter.sub <= length(N.sub))) {
      j <- 1
      while ((j <= length(N.time)) & (counter.time <= dim(testM)[1]) & sub.id == testM[counter.time]) {
        if (testM[counter.time, 2] == N.time[j]) {
          if (verbose == TRUE) {
            print(sprintf("MATCH    | sub counter = %04d | sub ID = %d | time counter = %02d | j = %02d | current time = %04.1f | sub time = %04.1f", counter.sub, sub.id, counter.time, j, N.time[j], testM[counter.time, 2]))
          }
          mat.wide[counter.sub, j] <- testM[counter.time,3]
          counter.time <- counter.time + 1
        } else {
          if (verbose == TRUE) {
            print(sprintf("NO MATCH | sub counter = %04d | sub ID = %d | time counter = %02d | j = %02d | current time = %04.1f | sub time = %04.1f", counter.sub, sub.id, counter.time, j, N.time[j], testM[counter.time, 2]))
          }

        }
        j <- j + 1
      }
      counter.sub <- counter.sub + 1
    }
  }
  dimnames(mat.wide) <- list(N.sub, N.time)
  mat.wide <- t(mat.wide)
  return(mat.wide)
}


