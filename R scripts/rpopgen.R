library(ggplot2)
library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
library(reshape2)

input.VCF <- read.vcfR("test2.vcf")

VCFFile <- input$VCF
input.VCF <- read.vcfR(VCFFile$datapath)

output$vcfsum <- renderPrint(input.VCF)
# ***** Object of Class vcfR *****
#   62 samples
# 1 CHROMs
# 7,750 variants

pop.data <- read.table("vcfpops.txt", sep = "\t", header = TRUE)

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




tree <- aboot(gl.test2, tree = "upgma",
              distance = bitwise.dist, sample = 100, showtree = F,
              cutoff = 50, quiet = T)

cols <- brewer.pal(n = nPop(gl.test2), name = "Dark2")
plot.phylo(tree, cex = 0.8, font = 2, adj = 0, tip.color =  cols[pop(gl.test2)])
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
legend('topleft', legend = c("B. bubalis", "Cattle", "Plains", "Wood"), fill = cols,
       border = FALSE, bty = "n", cex = 2)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")


########### PCA Plot
pca <- glPca(gl.test2, nf = 4)
barplot(100*pca$eig/sum(pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)


pca.scores <- as.data.frame(pca$scores)
pca.scores$pop <- pop(gl.test2)


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
pnw.dapc <- dapc(gl.test2, n.pca = 4, n.da = 2)

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
