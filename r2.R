library(tidyverse)
library(MASS)

print("loading data")
d <- read.delim("/outs/the_greater_big_data_file_bigger.txt", header=TRUE, sep=" ")
print("loaded data. running models. breathe.")

d$rep.fitted <- (d$rep.fitted.dip + d$rep.fitted.hap)/2
d$collision <- (d$collision.dip + d$collision.hap)/2

# d <- d %>% filter(transcripts < 15)
d <- d %>% filter(transcripts < 2.5e4)

head(d)

mod <- glm(formula = SNM ~ base.type + rep.fitted + protein + distance + 
             collision + I(distance^2) + distance:collision + transcripts:length + 
             protein:distance + base.type:distance + base.type:collision, 
           family = "binomial", data = d)
summary(mod)
print(1 - (summary(mod)$deviance/summary(mod)$null.deviance))
total_r2 <- 1 - (summary(mod)$deviance/summary(mod)$null.deviance)

# distance^2
mod <- glm(formula = SNM ~ base.type + rep.fitted + protein + distance + 
             collision + distance:collision + transcripts:length + 
             protein:distance + base.type:distance + base.type:collision, 
           family = "binomial", data = d)
summary(mod)
print(1 - (1 - (summary(mod)$deviance/summary(mod)$null.deviance))/total_r2)

# distance:collision
mod <- glm(formula = SNM ~ base.type + rep.fitted + protein + distance + 
             collision + I(distance^2) + transcripts:length + 
             protein:distance + base.type:distance + base.type:collision, 
           family = "binomial", data = d)
summary(mod)
print(1 - (1 - (summary(mod)$deviance/summary(mod)$null.deviance))/total_r2)

# protein:distance
mod <- glm(formula = SNM ~ base.type + rep.fitted + protein + distance + 
             collision + I(distance^2) + distance:collision + transcripts:length + 
             base.type:distance + base.type:collision, 
           family = "binomial", data = d)
summary(mod)
print(1 - (1 - (summary(mod)$deviance/summary(mod)$null.deviance))/total_r2)

# base.type:distance
mod <- glm(formula = SNM ~ base.type + rep.fitted + protein + distance + 
             collision + I(distance^2) + distance:collision + transcripts:length + 
             protein:distance + base.type:collision, 
           family = "binomial", data = d)
summary(mod)
print(1 - (1 - (summary(mod)$deviance/summary(mod)$null.deviance))/total_r2)

# base.type:collision 
mod <- glm(formula = SNM ~ base.type + rep.fitted + protein + distance + 
             collision + I(distance^2) + distance:collision + transcripts:length + 
             protein:distance + base.type:distance, 
           family = "binomial", data = d)
summary(mod)
print(1 - (1 - (summary(mod)$deviance/summary(mod)$null.deviance))/total_r2)
