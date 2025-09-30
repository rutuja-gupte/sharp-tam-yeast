library(tidyverse)
library(MASS)

print("loading data")
d <- read.delim("/outs/the_greater_big_data_file.txt", header=TRUE, sep=" ")
print("loaded data. running models. breathe.")

d$rep.fitted <- (d$rep.fitted.dip + d$rep.fitted.hap)/2
d$collision <- (d$collision.dip + d$collision.hap)/2

d <- d %>% filter(transcripts < 15)

head(d)

## stepping for SNM

lm1 <- glm(SNM ~ 
             base.type +
             collision +
             transcripts +
             rep.fitted +
             protein, d, family="binomial")
print(summary(lm1))

lm0 <- glm(SNM ~ base.type + rep.fitted, d, family="binomial")
print(summary(lm0))

lm2 <- glm(SNM ~ base.type +
             collision +
             transcripts +
             rep.fitted +
             protein +
             base.type * collision +
             transcripts * collision +
             collision * rep.fitted +
             rep.fitted * protein +
             base.type * transcripts +
             base.type * rep.fitted +
             base.type * protein +
             collision * protein +
             transcripts * rep.fitted +
             transcripts * protein, 
           d, family="binomial")
print(summary(lm2))

mod <- stepAIC(lm1, direction="both", scope=list(upper=lm2,lower=lm0))

print(summary(mod))

## Stepping for INDEL

lm1 <- glm(INDEL ~ 
             collision +
             transcripts +
             rep.fitted +
             protein, d, family="binomial")
print(summary(lm1))

lm0 <- glm(INDEL ~ rep.fitted, d, family="binomial")
print(summary(lm0))

lm2 <- glm(INDEL ~ base.type +
             collision +
             transcripts +
             rep.fitted +
             protein +
             base.type * collision +
             transcripts * collision +
             collision * rep.fitted +
             rep.fitted * protein +
             base.type * transcripts +
             base.type * rep.fitted +
             base.type * protein +
             collision * protein +
             transcripts * rep.fitted +
             transcripts * protein, 
           d, family="binomial")
print(summary(lm2))

mod <- stepAIC(lm1, direction="both", scope=list(upper=lm2,lower=lm0))

print(summary(mod))