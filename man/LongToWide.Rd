\name{LongToWide}
\alias{LongToWide}
\title{Long To Wide}

\description{
Reshape data from long to wide format
} 

\usage{
mat.wide <- LongToWide(subject.id, time.sequence, outcome, verbose = FALSE)
}

\arguments{
\item{subject.id}{A column of subject identifiers}
\item{time.sequence}{A column of time points}
\item{outcome}{A column of measures to be transposed}
\item{verbose}{TRUE or FALSE (DEFAULT = FALSE) - prints output to screen. Usefull for ensuring accuracy}
}

\value{
\item{mat.wide}{a matrix of data in wide format. columns = timeseq, rows = outcome for each subject}
}

\details{
If data is not sorted by subject.id, a warning message
will appear. The function will then sort the data by subject.id.
It is reccomended that the user sorts the data prior to using this function.
}

\examples{
set.seed(4)
df <- data.frame("subject.id" = rep(seq(from = 101, to = 105, by = 1), 4),
        "time.sequence" = sort(rep(seq(from = 0.1, to = 0.4, by = .1), 5)),
        "outcome" = sample(1:20, 20, replace=F))
df
df <- arrange(df, subject.id)
df
mat.wide <- LongToWide(df$subject.id, df$time.sequence, df$outcome)
mat.wide
}

\seealso{
\code{\link{tvmediation}}, \code{\link{tvma}}, \code{\link{smoker}}
}