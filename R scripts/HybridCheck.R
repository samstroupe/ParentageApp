library(bit64)
library("HybridCheck")
library("ape")

####################### 
Seq <- read.dna("Chr29.fasta", format= "fasta")
Test <- HC$new(Seq)

# Set the parameters for analysis #
Test$setParameters("TripletGeneration", Method = 2L, DistanceThreshold = 0.1)
Test$setParameters("SSAnalysis", WindowSize = 1000L, StepSize = 1L)

triplets <- list(c("29 BosTau_Ref", "29 RRYNP2", "29 CCSP662"))
Test$analyzeSS(triplets)
Test$findBlocks(triplets)
 
Test$setParameters("Plotting", What = "Lines", LegendFontSize = 10)
CCSP662 <- Test$plotTriplets(c("29 BosTau_Ref", "29 RRYNP2", "29 CCSP662"))
                 
rm(Seq)
rm(Test)
                 
        
########################
## CCSP662 ##
df <- CCSP662[[1]][["plot_env"]][["data"]]
df$Chromosome <- xyz
data = df[c("Chromosome","ActualStart", "ActualEnd", "AB", "BC")]
names(data)[names(data) == "AB"] <- "Bos_Bis"
names(data)[names(data) == "BC"] <- "Bis_Bis"
data$ActualStart <- as.integer64(data$ActualStart)
data$ActualEnd <- as.integer64(data$ActualEnd)
write.table(data, file = "CCSP662_xyz.txt", sep = "\t", col.names = TRUE, row.names = FALSE)


## Wrie in R code
# do awk '$4 > $5' ${i}.txt > ${i}_introgression_blocks.txt

# Sort and merge within 15000 base pairs

# Sum base pair and get average seqsim for Bis-Bos and Bis-Bis for introgress regions