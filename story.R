library(tidyverse)

print("loading data")
d <- read.delim("/outs/the_greater_big_data_file.txt", header=TRUE, sep=" ")
print("loaded data. breathe.")

d$rep.fitted <- (d$rep.fitted.dip + d$rep.fitted.hap)/2
d$collision <- (d$collision.dip + d$collision.hap)/2

head(d)

d <- d %>% filter(transcripts < 15)

head(d)

mod <- glm(SNM ~ 
             base.type +
             stops, d,
           family="binomial")

print(summary(mod))

mod <- glm(SNM ~ 
             base.type +
             transcripts, d,
           family="binomial")

print(summary(mod))

mod <- glm(SNM ~ 
             base.type +
             transcripts +
             rep.fitted, d,
           family="binomial")

print(summary(mod))

mod <- glm(SNM ~ 
             base.type +
             transcripts*rep.fitted, d,
           family="binomial")

print(summary(mod))

mod <- glm(SNM ~ 
             base.type +
             protein, d,
           family="binomial")

print(summary(mod))

mod <- glm(SNM ~ 
             base.type +
             protein +
             rep.fitted, d,
           family="binomial")

print(summary(mod))

mod <- glm(SNM ~ 
             base.type +
             protein*rep.fitted, d,
           family="binomial")

print(summary(mod))

mod <- glm(SNM ~ 
             base.type +
             protein*rep.fitted +
             transcripts*protein +
             transcripts*rep.fitted, d,
           family="binomial")

print(summary(mod))

mod <- glm(SNM ~ 
             base.type +
             collision + 
             protein*rep.fitted +
             transcripts*protein +
             transcripts*rep.fitted, d,
           family="binomial")

print(summary(mod))

mod <- glm(SNM ~ 
             base.type*collision +
             protein*rep.fitted +
             transcripts*protein +
             transcripts*rep.fitted, d,
           family="binomial")

print(summary(mod))


mod <- glm(SNM ~ 
             base.type*collision +
             base.type*transcripts +
             protein*rep.fitted +
             transcripts*protein +
             transcripts*rep.fitted, d,
           family="binomial")

print(summary(mod))

mod <- glm(SNM ~ 
             base.type*collision +
             base.type*transcripts +
             transcripts*collision +
             protein*rep.fitted +
             transcripts*protein +
             transcripts*rep.fitted, d,
           family="binomial")

mod <- glm(SNM ~ 
             base.type*collision +
             base.type*transcripts +
             transcripts*collision +
             rep.fitted*collision +
             protein*rep.fitted +
             transcripts*protein +
             transcripts*rep.fitted, d,
           family="binomial")

print(summary(mod))

### INDEL time

mod <- glm(INDEL ~ 
             base.type +
             stops, d,
           family="binomial")

print(summary(mod))

mod <- glm(INDEL ~ 
             base.type +
             transcripts, d,
           family="binomial")

print(summary(mod))

mod <- glm(INDEL ~ 
             base.type +
             transcripts +
             rep.fitted, d,
           family="binomial")

print(summary(mod))

mod <- glm(INDEL ~ 
             base.type +
             transcripts*rep.fitted, d,
           family="binomial")

print(summary(mod))

mod <- glm(INDEL ~ 
             base.type +
             protein, d,
           family="binomial")

print(summary(mod))

mod <- glm(INDEL ~ 
             protein +
             rep.fitted, d,
           family="binomial")

print(summary(mod))

mod <- glm(INDEL ~ 
             protein*rep.fitted, d,
           family="binomial")

print(summary(mod))

mod <- glm(INDEL ~ 
             protein*rep.fitted +
             transcripts*protein +
             transcripts*rep.fitted, d,
           family="binomial")

print(summary(mod))

mod <- glm(INDEL ~ 
             collision + 
             protein*rep.fitted +
             transcripts*protein +
             transcripts*rep.fitted, d,
           family="binomial")

print(summary(mod))

mod <- glm(INDEL ~ 
             protein*rep.fitted +
             transcripts*protein +
             transcripts*rep.fitted, d,
           family="binomial")

print(summary(mod))


mod <- glm(INDEL ~ 
             protein*rep.fitted +
             transcripts*protein +
             transcripts*rep.fitted, d,
           family="binomial")

print(summary(mod))

mod <- glm(INDEL ~ 
             transcripts*collision +
             protein*rep.fitted +
             transcripts*protein +
             transcripts*rep.fitted, d,
           family="binomial")

mod <- glm(INDEL ~ 
             transcripts*collision +
             rep.fitted*collision +
             protein*rep.fitted +
             transcripts*protein +
             transcripts*rep.fitted, d,
           family="binomial")

print(summary(mod))




