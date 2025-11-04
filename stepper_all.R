library(tidyverse)
library(MASS)

print("loading data")
d <- read.delim("/outs/the_greater_big_data_file_bigger.txt", header=TRUE, sep=" ")
print("loaded data. running models. breathe.")

d$rep.fitted <- (d$rep.fitted.dip + d$rep.fitted.hap)/2
d$collision <- (d$collision.dip + d$collision.hap)/2

# d <- d %>% filter(transcripts < 15)
d <- d %>% filter(transcripts < 2.5e4)

# d <- d %>% filter(!is.na(length))

head(d)

## Stepping again for SNM

lm1 <- glm(SNM ~ 
             PLACEHOLDER, d, family="binomial")
print(summary(lm1))
print(1 - (summary(lm1)$deviance/summary(lm1)$null.deviance))

lm0 <- glm(SNM ~ base.type + rep.fitted, d, family="binomial")
print(summary(lm0))
print(1 - (summary(lm0)$deviance/summary(lm0)$null.deviance))

lm2 <- glm(SNM ~ (base.type + rep.fitted + protein +
                    transcripts + collision)^2, 
           d, family="binomial")
print(summary(lm2))
print(1 - (summary(lm2)$deviance/summary(lm2)$null.deviance))

mod <- stepAIC(lm1, direction="both", scope=list(upper=lm2,lower=lm0))

print(summary(mod))