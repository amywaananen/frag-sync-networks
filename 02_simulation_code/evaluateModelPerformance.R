##  Evaluate model performance ####

nboot <- 100
levels <- c(1, 5, 20, 40, 50, 80, 100, 200)
out <- matrix(nrow = nboot*length(levels), ncol = 7)
row <- 1

for(m in 1:nboot){

  for(j in 1:length(levels)){

    SD_ST <- levels[j]  # standard deviation of distribution for sampling timing by space
    SD_SG <- levels[j] # standard deviation of distribution for sampling S-alleles by space
    SD_TG <- levels[j] # standard deviation of distribution for sampling S-alleles by time
    SD_STG <- levels[j] # standard deviation of distribution for sampling S-alleles by time

    source('makePop.R')

    out[row, 1] <- m
    out[row, 2] <- levels[j]
    out[row, 3] <- cor(pop$x, pop$start2)
    out[row, 4] <- cor(pop$x, pop$S1SG)
    out[row, 5] <- cor(pop$start2, pop$S1TG)
    out[row, 6] <- cor(pop$start2, pop$S1STG)
    out[row, 7] <- cor(pop$x, pop$S1STG)

    row <- row + 1
  }
}

out <- as.data.frame(out)
colnames(out) <- c('iteration','sd','xStart2','xS1SG','start2S1TG','start2S1STG','cS1STG')
# install.packages('tidyverse')
summary(out)

library(tidyverse)
ggplot(out, aes(x = log(sd), y = xStart2)) +
  geom_point() +
  geom_boxplot(aes(group =sd)) +
  theme_bw()