library(tidyverse)

print("loading data")
d <- read.delim("/outs/the_greater_big_data_file.txt", header=TRUE, sep=" ")
print("loaded data. running models. breathe.")

d$rep.fitted <- (d$rep.fitted.dip + d$rep.fitted.hap)/2
d$collision <- (d$collision.dip + d$collision.hap)/2

d <- d %>% filter(transcripts < 15)

head(d)

## This is the model with all two way options

mod <- glm(SNM ~
             base.type +
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
             transcripts * protein, d,
           family="binomial")

print(summary(mod))

## This the model after stepping

mod <- glm(SNM ~ base.type + collision + transcripts + rep.fitted +
             protein + base.type * collision + transcripts * collision +
             collision * rep.fitted + rep.fitted * protein, d,
           family="binomial")

print(summary(mod))

## INDEL model with all two-way interactions

mod <- glm(INDEL ~
             collision +
             transcripts +
             rep.fitted +
             protein +
             transcripts * collision +
             collision * rep.fitted +
             rep.fitted * protein +
             collision * protein +
             transcripts * rep.fitted +
             transcripts * protein, d,
           family="binomial")

print(summary(mod))


library(MASS)
library(glmnet)

d$transcripts.scaled <- scale(d$transcripts, scale=FALSE)

## LASSO for SNM

x = model.matrix(SNM ~
                   base.type +
                   collision +
                   transcripts.scaled +
                   rep.fitted +
                   protein +
                   base.type * collision +
                   transcripts.scaled * collision +
                   collision * rep.fitted +
                   rep.fitted * protein +
                   base.type * transcripts.scaled +
                   base.type * rep.fitted +
                   base.type * protein +
                   collision * protein +
                   transcripts.scaled * rep.fitted +
                   transcripts.scaled * protein, d)[,-1]
y = d$SNM
print("data ready")

lm.lasso = glmnet(x,y,alpha=1, family="binomial")
print("cv started")
cv.out = cv.glmnet(x,y,alpha=1, family="binomial")
lambda.best = cv.out$lambda.min
lambda.best

print(coef(lm.lasso, s=lambda.best))

## LASSO for INDEL

x = model.matrix(INDEL ~
                   base.type +
                   collision +
                   transcripts.scaled +
                   rep.fitted +
                   protein +
                   base.type * collision +
                   transcripts.scaled * collision +
                   collision * rep.fitted +
                   rep.fitted * protein +
                   base.type * transcripts.scaled +
                   base.type * rep.fitted +
                   base.type * protein +
                   collision * protein +
                   transcripts.scaled * rep.fitted +
                   transcripts.scaled * protein, d)[,-1]
y = d$INDEL
print("data ready")

lm.lasso = glmnet(x,y,alpha=1, family="binomial")
print("cv started")
cv.out = cv.glmnet(x,y,alpha=1, family="binomial")
lambda.best = cv.out$lambda.min
lambda.best

print(coef(lm.lasso, s=lambda.best))

## Stepping again for SNM

lm1 <- glm(SNM ~ 
             base.type +
             collision +
             transcripts +
             rep.fitted +
             protein +
             base.type * collision + 
             transcripts * collision + 
             collision * rep.fitted + 
             rep.fitted * protein, d, family="binomial")
print(summary(lm1))
print(1 - (summary(lm1)$deviance/summary(lm1)$null.deviance))

lm0 <- glm(SNM ~ base.type + rep.fitted, d, family="binomial")
print(summary(lm0))
print(1 - (summary(lm0)$deviance/summary(lm0)$null.deviance))

lm2 <- glm(SNM ~ base.type + rep.fitted + protein +
             transcripts + collision +
             rep.fitted*transcripts*collision*base.type*protein, 
           d, family="binomial")
print(summary(lm2))
print(1 - (summary(lm2)$deviance/summary(lm2)$null.deviance))

mod <- stepAIC(lm1, direction="both", scope=list(upper=lm2,lower=lm0))

print(summary(mod))

## Stepping for INDELs

lm1 <- glm(INDEL ~ 
             collision +
             transcripts +
             rep.fitted +
             protein +
             transcripts * collision + 
             collision * rep.fitted + 
             rep.fitted * protein, d, family="binomial")
print(summary(lm1))
print(1 - (summary(lm1)$deviance/summary(lm1)$null.deviance))

lm0 <- glm(INDEL ~ rep.fitted, d, family="binomial")
print(summary(lm0))
print(1 - (summary(lm0)$deviance/summary(lm0)$null.deviance))

lm2 <- glm(INDEL ~ base.type + rep.fitted + protein +
             transcripts + collision +
             rep.fitted*transcripts*collision*base.type*protein, 
           d, family="binomial")
print(summary(lm2))
print(1 - (summary(lm2)$deviance/summary(lm2)$null.deviance))

mod <- stepAIC(lm1, direction="both", scope=list(upper=lm2,lower=lm0))

print(summary(mod))
