library(tidyverse)

print("loading data")
d <- read.delim("/outs/the_greater_big_data_file.txt", header=TRUE, sep=" ")
print("loaded data. running models. breathe.")

d$rep.fitted <- (d$rep.fitted.dip + d$rep.fitted.hap)/2
d$collision <- (d$collision.dip + d$collision.hap)/2

d <- d %>% filter(transcripts < 15)

head(d)

mod <- glm(SNM ~ 
             base.type +
             collision +
             transcripts +
             rep.fitted +
             protein +
             base.type * collision + 
             transcripts * collision + 
             collision * rep.fitted + 
             rep.fitted * protein, 
           d,
           family="binomial")

print(summary(mod))

mod <- glm(SNM ~ 
             triplet +
             collision +
             transcripts +
             rep.fitted +
             protein +
             base.type * collision + 
             transcripts * collision + 
             collision * rep.fitted + 
             rep.fitted * protein, 
           d,
           family="binomial")

print(summary(mod))



library(MASS)
library(glmnet)

d$transcripts <- scale(d$transcripts, scale=FALSE)

x = model.matrix(SNM ~ 
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
                   transcripts * protein, d)[,-1]
y = d$SNM
print("data ready")

lm.lasso = glmnet(x,y,alpha=1, family="binomial")
print("cv started")
cv.out = cv.glmnet(x,y,alpha=1, family="binomial")
lambda.best = cv.out$lambda.min
lambda.best

print(coef(lm.lasso, s=lambda.best))
