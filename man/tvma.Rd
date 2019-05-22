\name{tvma}
\title{Time Varying Mediation Analysis}
\alias{tvma}

\description{
Calculates Time-Varying Mediation Effect and SD
} 

\usage{
results <- tvma(treatment, t.seq, mediator, outcome)
}

\arguments{
\item{treatment}{A vector with treatment values.}
\item{t.seq}{A vector of time points of observations.}
\item{mediator}{A matrix of mediator values.}
\item{outcome}{A matrix of outcome values.}
\item{t.est}{A vector of time points to make the estimation. DEFAULT is t.seq}
\item{plot}{A logical TRUE/FALSE argument that plots the mediation effect. DEFAULT is TRUE.}
\item{CI}{A string to determine either "asymptotic" or "boot" derived confidence intervals. DEFAULT is "boot".}
\item{replicates}{Number of replicates for bootstrapping confidence intervals. DEFAULT is 1000.}
}

\value{
\item{hat.alpha.1}{Estimated treatment effect on mediator.}
\item{hat.beta.2}{Estimated mediator effect on outcome.}
\item{mediation.effect}{Time varying mediation effect.}
\item{sd}{Estimated standard deviation of the mediation effect. *Only returned for asymptotic estimation.}
\item{CI.upper}{Upper confidence intervals.}
\item{CI.lower}{Lower confidence intervals.}
}

\note{
*Currently supports 2 treatment options. Future realeases may expand number of treatment options
}

\examples{
data(smoker)
# REDUCE DATA SET TO ONLY 2 TREATMENT CONDITIONS
smoker.parsed <- smoker[smoker$ConditionID != "comboNRT", ]
# GENERATE WIDE FORMATTED MEDIATORS
mediator <- LongToWide(smoker.parsed$SubjectID, smoker.parsed$timeseq, smoker.parsed$NegMoodLst15min)
# GENERATE WIDE FORMATTED OUTCOMES
outcome <- LongToWide(smoker.parsed$SubjectID, smoker.parsed$timeseq, smoker.parsed$cessFatig)
# GENERATE A BINARY TREATMENT VARIABLE
treatment <- as.numeric(unique(smoker.parsed[ , c("SubjectID","varenicline")])[ , 2])-1
# GENERATE A VECTOR OF UNIQUE TIME POINTS
t.seq <- sort(unique(smoker.parsed$timeseq))

# COMPUTE TIME VARYING MEDIATION ANALYSIS USING BOOTSTRAPPED CONFIDENCE INTERVALS
results <- tvm(treatment, t.seq, mediator, outcome)
# COMPUTE TIME VARYING MEDIATION ANALYSIS FOR SPECIFIED POINTS IN TIME USING 500 REPLICATES
results <- tvm(treatment, t.seq, mediator, outcome, t.est = c(0.2, 0.4, 0.6, 0.8), replicates = 500)
# COMPUTE TIME VARYING MEDIATION ANALYSIS WITH ASYMPTOTIC CONFIDENCE INTERVALS
results <- tvm(treatment, t.seq, mediator, outcome, CI = "asymptotic")
}

\references{
Fan, J. and Gijbels, I. (1996). Local polynomial modelling and its applications: monographs on statistics and applied probability 66 66. CRC Press.
Fan, J. and Zhang, W. (1999). Statistical estimation in varying coefficient models. The annals of Statistics 27 1491-1518.
Fan, J. and Zhang, W. (2000). Two-step estimation of functional linear models with applications to longitudinal data. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 62 303-322.
}

\seealso{
\code{\link{tvmediation}}, \code{\link{LongToWide}}, \code{\link{smoker}}
}