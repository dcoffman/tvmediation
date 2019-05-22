## tidyverse MUST BE CALLED BEFORE Hmisc OR THERE WILL BE ISSUES AFTERUSING UpData() LATER
library(dplyr)
library(tidyverse)
library(Hmisc)
library(ggplot2)
library(gridExtra)

# library(MASS)
# library(np)
# library(locpol)
# library(locfit)
# library(rqPen)
# library(KernSmooth)

## LOAD RAW DATA
load("C:/Users/hzjr/Box Sync/R Package/Smoking data/Donna's initial exploration R code/wisconsin.Rdata")


## DEFINE VARIABLES OF INTEREST, CREATE smoker DATAFRAME, STORE A COPY
my.vars <- c("SubjectID","DaysFromTQD","time.of.day","WantToSmokeLst15min","NegMoodLst15min","cessFatig","ConditionID","CigCount_Today")
smoker <- dplyr::select(dat, one_of(my.vars))  # CREATE smoker WITH SELECTED VARIABLES
rawsmoker <- smoker  # STORE A COPY


## CLEAN UP SOME DATA
rownames(smoker) <- NULL  # REMOVE ROWNAMES
smoker <- smoker %>%
  filter(!DaysFromTQD < 0) %>%  # ONLY USE + VALUES FROM DaysFrom TQD
  filter(!time.of.day == 2) %>%  # GET RID OF WHEN time.of.day == 2 (random), THEN RECODE 1 TO 0 AND 3 TO 1 (AM & PM RESPECTIVELY)
  filter(!ConditionID < 2 & !ConditionID >4)  # GET RID OF CONDITIONS NOT ASSOCIATED WITH TREATMENT

smoker$time.of.day[smoker$time.of.day == 1] <- 0
smoker$time.of.day[smoker$time.of.day == 3] <- 1


## CREATE NEW VARIABLES
smoker <- add_column(smoker, timeseq = as.numeric(smoker$DaysFromTQD) + (as.numeric(smoker$time.of.day) / 2),
                     .before = "WantToSmokeLst15min")  # CREATE timeseq VARIABLE IN dataframe
# CREATE DUMMY TREATMENT VARIABLES
smoker <- add_column(smoker, patch = as.numeric(smoker$ConditionID == 2),
                     .before = "CigCount_Today")
smoker <- add_column(smoker, varenicline = as.numeric(smoker$ConditionID == 3),
                     .before = "CigCount_Today")
smoker <- add_column(smoker, comboNRT = as.numeric(smoker$ConditionID == 4),
                     .before = "CigCount_Today")

## CREATE FACTORS & LEVELS
my.factors <- c("time.of.day", "WantToSmokeLst15min", "NegMoodLst15min", "cessFatig", "ConditionID", "patch", "varenicline", "comboNRT")
smoker[my.factors] <- lapply(smoker[my.factors], factor)
levels(smoker$time.of.day) <- c("AM","PM")
levels(smoker$WantToSmokeLst15min) <- c("Not at all - 1","2", "3", "4", "5", "6", "Extremely - 7")
levels(smoker$NegMoodLst15min) <- c("Not at all - 1","2", "3", "4", "5", "6", "Extremely - 7")
levels(smoker$cessFatig) <- c("Strongly Disagree - 1", "2", "3", "4", "5", "6", "Strongly Agree - 7")
levels(smoker$ConditionID) <- c("patch", "varenicline", "comboNRT")
levels(smoker$patch) <- c("No", "Yes")
levels(smoker$varenicline) <- c("No", "Yes")
levels(smoker$comboNRT) <- c("No", "Yes")


## DEFINE & APPLY LABELS TO VARIABLES
my.labels <- c(SubjectID = "subject ID",
               DaysFromTQD = "Number of days from quit date",
               time.of.day = "Time of day (0 = am, 1 = pm)",
               timeseq = "Number of days from quit data (.5 indicates pm)",
               WantToSmokeLst15min = "How did you feel in the last 15 min: wanting to smoke (1 = not at all, 7 = extremely)",
               NegMoodLst15min = "How did you feel in the last 15 min: Negative mood (1 = not at all, 7 = extremely)",
               cessFatig = "Cessation Fatigue - I am tired of trying to quit smoking (1 = strongly disagree, 7 = strongly agree)",
               ConditionID = "Treatment group (2 = patch, 3 = varenicline, 4 = combination nicotine replacement therapy",
               patch = "Received patch (0 = No, 1 = Yes)",
               varenicline = "Received varenicline (0 = No, 1 = Yes)",
               comboNRT = "Received combination nicotine replacement therapy (0 = No, 1 = Yes)",
               CigCount_Today = "Cigarettes smoked over entire day")

smoker <- Hmisc::upData(smoker, labels = my.labels)


## ARRANGE BY ID, THEN DAYS FROM QUIT, THEN TIME OF DAY
smoker <- arrange(smoker, SubjectID, DaysFromTQD, time.of.day)


## LOOK AT SOME TABLES...
Hmisc::contents(smoker)
Hmisc::describe(smoker)
psych::describe(smoker)



## AGGREGATE THE MEAN AND SD OF THE DATA BY TIME AND CONDITION
smoker.agg <- aggregate(cbind(NegMoodLst15min = as.numeric(smoker[, 6]), cessFatig = as.numeric(smoker[, 7])),
          list(timeseq = smoker$timeseq, ConditionID = smoker$ConditionID),
          function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE)))




## PLOT LINE GRAPHS OF MEAN AND SD OF NegMood & cessFatig BY TREATMENT GROUP OVER TIME
g1 <- ggplot(data = smoker.agg, aes(x = timeseq, y = NegMoodLst15min[, 1], color = ConditionID)) + #group=ConditionID, color = ConditionID)) +
  #geom_bar(stat = "identity", color = "black",
  #         position = position_dodge()) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = smoker.agg$NegMoodLst15min[, 1] - smoker.agg$NegMoodLst15min[, 2],
                    ymax = smoker.agg$NegMoodLst15min[, 1] + smoker.agg$NegMoodLst15min[, 2]), width = .2,
                position = position_dodge(.9)) +
  facet_wrap(~ConditionID, nrow = 3) +
  ylim(0, 7) +
  labs(title = "",
       y = "Negative Mood in last 15 min",
       x = "Number of days since quit date",
       color = "Treatment Condition") +
  theme(legend.position = "none")

g2 <- ggplot(data = smoker.agg, aes(x = timeseq, y = cessFatig[, 1], color = ConditionID)) + #group=ConditionID, color = ConditionID)) +
  #geom_bar(stat = "identity", color = "black",
  #         position = position_dodge()) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = smoker.agg$cessFatig[, 1] - smoker.agg$cessFatig[, 2],
                    ymax = smoker.agg$cessFatig[, 1] + smoker.agg$cessFatig[, 2]), width = .2,
                position = position_dodge(.9)) +
  facet_wrap(~ConditionID, nrow =3) +
  ylim(0, 7) +
  labs(title = "",
       y = "Cessation Fatigue",
       x = "Number of days since quit date",
       color = "Treatment Condition") +
  theme(legend.position = "none")

grid.arrange(g1, g2, nrow = 1, top = "Changes in mood over time by treatment")





############## WORK ##################
# Y = smoker$cessFatig
# M = smoker$NegMoodLst15min
# t.seq = smoker$timeseq
# trt = smoker$varenicline
# t.est = t.seq

# for (i in 1:length(colnames(smoker))) {
#   b <- colnames(smoker)[i]
#   a <- attributes(eval(parse(text=paste("smoker$",b,sep=""))))
#   print(a)
# }




