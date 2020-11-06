library(ggplot2)
library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
library(reshape2)
library(sequoia)

input.VCF <- read.vcfR("test1.vcf")
pop.data <- read.csv("test1pops.csv")


# We can now check that all the samples in the VCF and the population data frame are included:
all(colnames(input.VCF@gt)[-1] == pop.data$AccessID)
# TRUE

# Converting the dataset to a genlight object
gl.object <- vcfR2genlight(input.VCF)

# Specify ploidy
ploidy(gl.object) <- 2

# Add population data to genlight object
pop(gl.object) <- pop.data$Population

gl.object

########### Subset VCF 1000 random varients ###############

subset.1 <- sample(size = 1000, x= c(1:nrow(input.VCF)))
input.VCF.sub1 <- input.VCF[subset.1,]

write.vcf(input.VCF.sub1, "VCFsub1")

VCF <- read.vcfR("VCFsub1")


############ Tree #################

tree <- aboot(gl.object, tree = "upgma",
              distance = bitwise.dist, sample = 100, showtree = F,
              cutoff = 50, quiet = T)

cols <- brewer.pal(n = nPop(gl.object), name = "Dark2")
plot.phylo(tree, cex = 0.8, font = 2, adj = 0, tip.color =  cols[pop(gl.object)])
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
legend('topleft', legend = c("B. bubalis", "Cattle", "Plains", "Wood"), fill = cols,
       border = FALSE, bty = "n", cex = 2)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")


########### PCA Plot
pca <- glPca(gl.object)
barplot(100*pca$eig/sum(pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)


pca.scores <- as.data.frame(pca$scores)
pca.scores$pop <- pop(gl.object)


  # determine the number of pc's to retain

w <- which(cumsum(100*pca$eig/sum(pca$eig)) >= 80)
w[1]

optim.num <- optim.a.score(pnw.dapc)
optim.num$best




set.seed(9)
p <- ggplot(pca.scores, aes(x=PC1, y=PC2, colour=pop))
p <- p + geom_point(size=2)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = cols)
p <- p + geom_hline(yintercept = 0)
p <- p + geom_vline(xintercept = 0)
p <- p + theme_bw()

p

################ Structure-like plot ##############
pnw.dapc <- dapc(gl.object, n.pca = 39, n.da = 1)

dapc.results <- as.data.frame(pnw.dapc$posterior)
dapc.results$pop <- pop(gl.test2)
dapc.results$indNames <- rownames(dapc.results)

dapc.results <- melt(dapc.results)

colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity')
p <- p + scale_fill_manual(values = cols)
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p


###############

dp <- extract.gt(input.VCF, element = "DP", as.numeric=TRUE)
sum(is.na(dp[,1]))

myMiss <- apply(dp, MARGIN = 2, function(x){ sum(is.na(x)) })
myMiss <- myMiss/nrow(input.VCF)


par(mar = c(12,4,4,2))
barplot(myMiss, las = 2, col = 1:12)
title(ylab = "Missingness (%)")


##################################

system.time( my_loci <- vcfR2loci(input.VCF) )
class(my_loci)


geno <- read.csv("Geno.csv")
class(geno)

ID5885 <- my_loci[2,1:100]
