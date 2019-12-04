library(shiny)
library(sequoia)
library(Rcpp)
library(pedantics)
library(kinship2)
library(ggplot2)
library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
library(reshape2)

###### UI #######
ui <- fluidPage(
  HTML('<p><img class="one" src="BisonPic.png" width="100%" height="100%"/></p>'),
  
      # <p><img src="BisonPic.png"/></p> class="one"
   # Application title
   titlePanel("Parentage Testing"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        fileInput("Geno", label = h3("Genotype File input")),
        
        fileInput("LHis", label = h3("Life History File input (.csv)")),
        actionButton("run", "Run Parentage"),
        
         radioButtons("radio", label = h3("Plot Pedigree"),
                     choices = list("Pedigree 1" = 1, "Pedigree 2" = 2), 
                     selected = 1),
        
        fileInput("VCF", label = h3("VCF File input")),
        fileInput("Pop", label = h3("Population File input (.csv)")),
        actionButton("runvcf", "Run VCF"),
        radioButtons("PCradio", label = h3("Radio buttons"),
                     choices = list("Optimal number of PCs retained" = 1, 
                                    "PCs that account for 80% of variance" = 2), 
                     selected = 2)
        
       
      ),
    
      # Show the pedigree
      mainPanel(
          
        # Output: Tabset
        tabsetPanel(type = "tabs",
                    tabPanel("Parentage Assignment",
                             downloadButton("downloadData", "Download"),
                             tableOutput("assignment")
                             ),
                    tabPanel("Pedigree Plot", plotOutput("pedigree")),
                    tabPanel("VCF Summary", verbatimTextOutput("vcfsum"),
                             textOutput("allIncludetext"),
                             verbatimTextOutput("allInclude"),
                             verbatimTextOutput("gl"),
                             plotOutput("eigenval")),
                    tabPanel("PCA Plot", plotOutput("PCA")),
                    tabPanel("Structure Plot", plotOutput("Structure"))
                    
        )
      )
   )
)



########### SERVER #############
server <- function(input, output) {
  
  
 observeEvent(input$run, {
    # Read in the Life History File
   LHisFile <- input$LHis
   if (is.null(LHisFile))
     return(NULL)
   LHis <- read.csv(LHisFile$datapath)
    # Remove population column for sequoia
    LHis1 <- LHis[ , c("ID", "Sex", "BY")]
    # Life History for pedigree 2
    LHis2 <- LHis

    # Read in the Genotype File
   GenoFile <- input$Geno
   if (is.null(GenoFile))
     return(NULL)
   Geno <- read.csv(GenoFile$datapath)
   row.names(Geno) <- Geno$X
   Geno <- Geno[,-1]
   Geno <- as.matrix(Geno)

   levels(LHis1$Sex) <- sub("female", "1", levels(LHis1$Sex))
   levels(LHis1$Sex) <- sub("male", "2", levels(LHis1$Sex))
    ParOUT <- sequoia(GenoM = Geno,
                    LifeHistData = LHis1,
                    MaxSibIter = 0)

    SeqOUT <- sequoia(GenoM = Geno,
                    SeqList = ParOUT,
                    MaxSibIter = 5)

    Ped <- SeqOUT$Pedigree
    Ped = Ped[order(as.numeric(Ped$id)), ]

  output$assignment <- renderTable({Ped[complete.cases(Ped), ]
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$Geno, "Parentage", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(Ped[complete.cases(Ped), ], file, row.names = F)
    }
  )
  
  output$pedigree <- renderPlot({
  
  if (input$radio == 1){
    plot <- drawPedigree(Ped, dots = "y", dotSize = .005)
  }
  if (input$radio == 2){
    # Plot the pedegree
    ped2 <- fixPedigree(Ped)
    # Order the id columns the same
    ped2 = ped2[order(ped2$id),]
    LHis2 = LHis2[order(LHis2$ID),]
    # Assign the Life History sex to the pedigree 
    ped2$Sex = LHis2$Sex
    
    ped2 <- pedigree(id = ped2$id, dadid = ped2$sire, momid = ped2$dam, sex = LHis2$Sex)
    plot <- plot(ped2)
  }
  plot
  })
  })
 
  
      ### VCF ###
 observeEvent(input$runvcf, {
   # read in VCF File
   VCFFile <- input$VCF
   input.VCF <- read.vcfR(VCFFile$datapath)
   # read in Population txt File
   PopFile <- input$Pop
   if (is.null(PopFile))
     return(NULL)
   pop.data <- read.csv(PopFile$datapath)

   # VCF summary output
   output$vcfsum <- renderPrint(input.VCF)
   
   # We can now check that all the samples in the VCF and the population data frame are included:
   output$allIncludetext <- renderText(" All samples are included in the VCF File and 
                                       the population data frame.")
   output$allInclude <- renderPrint(all(colnames(input.VCF@gt)[-1] == pop.data$AccessID))
   
   # Converting the dataset to a genlight object
   gl.object <- vcfR2genlight(input.VCF)
   
   # Specify ploidy
   ploidy(gl.object) <- 2
   
   # Add population data to genlight object
   pop(gl.object) <- pop.data$Population
   
   # genlight object summary output
   output$gl <- renderPrint(gl.object)
  
   ## Number of PCs to use in Analysis  
   # Optimal number of PCs based off of a-score
   optim.num <- optim.a.score(pnw.dapc)
   optim.num$best
   # Number of eigenvalues to retain to account for 80% of the variation
   w <- which(cumsum(100*pca$eig/sum(pca$eig)) >= 80)
   
   pca <- glPca(gl.object, nf = w[1])
   pnw.dapc <- dapc(gl.object, n.pca = w[1], n.da = nlevels(gl.object$pop)-1)
   
   if (input$PCradio == 1){
     pc.num <- optim.num$best
    }
   if (input$PCradio == 2){
     pc.num <- w[1]
   }
   
   # Make the PCA
   pca <- glPca(gl.object, nf = pc.num)
   output$eigenval <- renderPlot(barplot(100*pca$eig/sum(pca$eig), col = heat.colors(50), 
                                         main="PCA Eigenvalues", 
                                         ylab="Percent of variance explained",  
                                         xlab="Eigenvalues"))
   
   # Create a PCA plot
   pca.scores <- as.data.frame(pca$scores)
   pca.scores$pop <- pop(gl.object)
      
   set.seed(9) # random seed
   output$PCA <- renderPlot(ggplot(pca.scores, aes(x=PC1, y=PC2, colour=pop)) + 
     geom_point(size=2) + 
     stat_ellipse(level = 0.95, size = 1) + 
     geom_hline(yintercept = 0) + 
     geom_vline(xintercept = 0))
   
   # Create a compoplot or Structure Plot
   pnw.dapc <- dapc(gl.object, n.pca = pc.num, n.da = nlevels(gl.object$pop)-1)
   
   dapc.results <- as.data.frame(pnw.dapc$posterior)
   dapc.results$pop <- pop(gl.object)
   dapc.results$indNames <- rownames(dapc.results)
   
   dapc.results <- melt(dapc.results)
   
   colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")
   
   output$Structure <- renderPlot(ggplot(dapc.results, aes(x=Sample, 
                                                y=Posterior_membership_probability,
                                                fill=Assigned_Pop)) + 
     geom_bar(stat='identity') + 
     facet_grid(~Original_Pop, scales = "free") + 
     theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)))
 })
 
}

# Run the application 
shinyApp(ui = ui, server = server)

