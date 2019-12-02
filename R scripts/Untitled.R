library(devtools)
install_github('samstroupe/SimSNPData')

library(SimSNPData)
library(sequoia)
library(Rcpp)
library(pedantics)
library(kinship2)
library(stringr)

df <- SimSnpData(5, 10)
df <- AssignSex(df, 3, 2, 0)

pop2 <- MakeBabies(df, 2)

pop <- rbind(df, pop2$pop)

pop3 <- MakeBabies(pop, 3)
pop <- rbind(pop, pop3$pop)

Ped <- rbind(pop2$pedigree, pop3$pedigree)


Ped <- fixPedigree(Ped)
Ped = Ped[order(as.numeric(str_sub(Ped$id, 4,))), ]

ped <- pedigree(id = Ped$id, dadid = Ped$sire, momid = Ped$dam, sex = pop$sex)
plot(ped)
