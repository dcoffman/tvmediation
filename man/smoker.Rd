\name{smoker}
\alias{smoker}
\non_function{}
\title{Wisconsin Smokers' Health Study 2}
\description{
The dataset is \strong{simulated} based on the Wisconsin Smokers' Health Study 2.
}
\usage{data(smoker)}
\format{
	A data frame with 40,130 observations on 12 variables.
	\enumerate{
		\item \strong{SubjectID:}	Subject ID
		\item \strong{treatment:}	Treatment group (2 = patch, 3 = varenicline, 4 = combination nicotine replacement therapy)
		\item \strong{patch:}	Received patch (0 = No, 1 = Yes)
		\item \strong{varenicline:}	Received varenicline (0 = No, 1 = Yes)
		\item \strong{comboNRT:}	Received combination nicotine replacement therapy (0 = No, 1 = Yes)
		\item \strong{DaysFromTQD:}	Number of days from quit date
		\item \strong{time.of.day:}	Time of day (0 = am, 1 = pm)
		\item \strong{timeseq:}	Number of days from quit date (.5 indicates pm)
		\item \strong{WantToSmokeLst15min:}	How did you feel in the last 15 min: wanting to smoke (1 = not at all, 7 = extremely)
		\item \strong{NegMoodLst15min:}	How did you feel in the last 15 min: Negative mood (1 = not at all, 7 = extremely)
		\item \strong{cessFatig:}	Cessation fatigue - I am tired of trying to quit smoking (1 = strongly disagree, 7 = strongly agree)
		\item \strong{CigCount:}	Cigarettes smoked over entire day
		\item \strong{smoke_status:} Did you smoke (0 = No, 1 = Yes)
	}
}
\examples{
data(smoker)
}
\seealso{
\code{\link{tvmb}}, \code{\link{tvma_3trt}} \code{\link{tvma}}, \code{\link{LongToWide}}
}
\keyword{datasets}
