#install.packages(c("shiny", "sequoia", "Rcpp", "kinship2", "ggplot2", "vcfR", "poppr", "ape", "RColorBrewer",
    #               "reshape2", "adegenet", "cowplot", "Cairo", "shinyWidgets", "grDevices", "shinyjs", "shinythemes",
    #               "markdown", "DT"))
library(shiny)
library(sequoia)
library(Rcpp)
#library(pedantics)
library(kinship2)
library(ggplot2)
library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
library(reshape2)
library(adegenet)
library(cowplot)
library(Cairo)
library(shinyWidgets)
library(grDevices)
library(shinyjs)
library(shinythemes)
library(markdown)
library(DT)

###### UI #######
ui <- fluidPage(
   useShinyjs(),
      ## Set the theme of the APP
     theme = shinytheme("sandstone"),
    # theme = "sandstone.css",
         # Bison Picture as background image    
      setBackgroundImage(src="BisonPic.png"),
         # Bison picture as image fixed to top
   # HTML('<p><img class="one" src="BisonPic.png" width="100%" height="100%"/></p>'),
 tagList(
    
    
    navbarPage( position = "fixed-top", #tags$style(type="text/css", "body {padding-top: 70px;}"),
                
       title = div(img(src='bison_dna_white.png',
                               style="margin-top: -14px; padding-right:10px;padding-bottom:10px", 
                               height = 60)),
        tabPanel("Home",
                #includeHTML(rmarkdown::render("about.Rmd")))
                includeMarkdown("about.Rmd")),
       
            ## Panel for Sequoia Parentage Testing ##
   tabPanel("Parentage Assignment",
      sidebarLayout(
         sidebarPanel(width = 3,
               fileInput("Geno", h3("Genotype File:")),
                  helpText("*Must be in Plink.raw format"),
        
        
               fileInput("LHis", label = h3("Life History File:"),
                  accept = ".csv"),
                  helpText("*CSV file with ID, Sex, and BY (Birth Year)"
            ),
 
               tags$div(dropdownButton(
                  tags$h3("Parameters for Sequoia"),
                  checkboxInput("checkbox", label = "Print only complete parentage matches", value = FALSE),
                  numericInput("Err", "Genotype Error Rate:", min = 0, max = 1, value = 0.01),
                  numericInput("MM", "Max Mismatches:", min = 0, max = 10000, value = 20),
                  circle = FALSE, status = "warning", icon = icon("gear"), width = "300px",
                  tooltip = tooltipOptions(title = "Adjust Parameters for Parentage Assignment")),
                  style="display:inline-block"
            ),
               tags$div(dropdownButton(
                  tags$h3("Pedigree Options"),
                  radioButtons("radio", label = NULL,
                        choices = list("Pedigree 1" = 1, "Pedigree 2" = 2), 
                        selected = 1), circle = FALSE, status = "warning", icon = icon("gear"), width = "300px",
                        tooltip = tooltipOptions(title = "Pedigree Options")),
                        style="display:inline-block"
            ),
                           
         actionButton("run", "Run Parentage", class = "btn-success"
                      )
      ), ##End of parentage side Bar
         
        mainPanel(
           tabsetPanel(type = "tabs",
               tabPanel("Parentage Assignment",
                             downloadButton("downloadData", "Download"),
##                             tableOutput("assignment")
                              DTOutput("assignment")
                             ),
               tabPanel("Pedigree Plot", plotOutput("pedigree"))
            )
         ) ## End of parentage Main panel
      ) ##End Parentage sidebar Layout
   ), ## End parentage Tab Panel
   
      
      
      ## Panel for Population Genetic Testing using various methods ##
            tabPanel("Population Evaluation",
               sidebarLayout(
                  sidebarPanel(width = 3,  
                     # input files #
                     fileInput("VCF", label = h3("VCF File input")),
                     fileInput("Pop", label = h3("Population File input (.csv)")),
                           tags$div(dropdownButton(
                              tags$h3("Discriminant Analysis of Principal Components"),
                              radioButtons("PCradio", label = h4("Number of Principal Components to be retained:"),
                                 choices = list("4" = 1,
                                 "Optimal number determined by a-score" = 2,
                                 "80% of variance" = 3), 
                                 selected = 1),
                                 circle = FALSE, status = "warning", icon = icon("gear"), width = "300px",
                                 tooltip = tooltipOptions(title = "Adjust Parameters for DAPC")),
                                 style="display:inline-block"),
                           tags$div(dropdownButton(
                              tags$h3("Parameters for Cluster Plots"),
                                 sliderInput("n", "Max K:", min = 2, max = 10, value = 8),
                                 numericInput("pca", "Number of PCAs:", min = 2, max = 1000, value = 40),
                                 circle = FALSE, status = "warning", icon = icon("gear"), width = "300px",
                                 tooltip = tooltipOptions(title = "Adjust Parameters for Cluster Plots")),
                                 style="display:inline-block"),
        
                           actionButton("runvcf", "Run VCF", class = "btn-success")
         
      ), # End of PopGen sidebar Panel
      
          mainPanel(   
             tabsetPanel(type = "tabs",
                  tabPanel("VCF Summary", verbatimTextOutput("vcfsum"),
                             textOutput("allIncludetext"),
                             verbatimTextOutput("allInclude"),
                             verbatimTextOutput("gl"),
                             plotOutput("eigenval")),
                  tabPanel("PCA Plot", 
                             downloadButton("downloadPCA", "Download PCA"),
                             plotOutput("PCA")),
                  tabPanel("Discriminant Analysis of Principal Components", 
                             downloadButton("downloadStruct", "Download Plot"),
                             plotOutput("Structure")),
                  tabPanel("Clustering Plots", 
                             downloadButton("downloadClust", "Download Plot"),
                             plotOutput("PopGenPlots", height = 1000
                        ))
         )           
          ) # End of PopGen main Panel
       
      ) #End PopGen Sidebar Layout
   ) #End PopGen tab Panel
   
      
 ), #End NavbarPage  
 
 ## Style Tags
   tags$style(type="text/css", "body {padding-top: 70px;}"),
   tags$style(type = "text/css", "#map {height: calc(100vh - 53px) !important;}"), 
   tags$style(type = "text/css", ".container-fluid {padding-left:0px;padding-right:0px;}"),
   tags$style(type = "text/css", ".navbar {margin-bottom: .5px;}"),
   tags$style(type = "text/css", ".container-fluid .navbar-header 
               .navbar-brand {margin-left: 0px;}"

              
         )
      )
)






########### SERVER #############

server <- function(input, output) {
  
   # observeEvent(input$toggleSidebar, {
   #    shinyjs::toggle(id = "Sidebar")
   # })

   ############### SEQUOIA PARENTAGE ASSIGNMENT ###################  
 observeEvent(input$run, {
   
    withProgress(message = "Running Sequoia", value = 0, {
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
 
    incProgress(amount = 1/10, detail = "Converting Genotype File")
     Geno <- GenoConvert(InFile = GenoFile$datapath, InFormat = "raw")

   levels(LHis1$Sex) <- sub("female|f", "1", levels(LHis1$Sex), ignore.case = TRUE)
   levels(LHis1$Sex) <- sub("male|m", "2", levels(LHis1$Sex), ignore.case = TRUE)
   levels(LHis1$Sex) <- sub("unknown|u", "3", levels(LHis1$Sex), ignore.case = TRUE)
   
   incProgress(amount = 3/10, detail = "Assigning Potential Parents")
            Err <- input$Err  # adjust the parameters for sequoia using the widgets
            MaxMismatch <- input$MM
        
      ParOUT <- sequoia(GenoM = Geno,
                    LifeHistData = LHis1,
                    Err = Err, 
                    MaxMismatch = MaxMismatch, 
                    FindMaybeRel = T, CalcLLR = T)
      
   incProgress(amount = 6/10, detail = "Constructing Full Pedigree")
      

    Ped <- ParOUT$PedigreePar
    Ped <- Ped[order(as.numeric(Ped$id)), ]
                 })  


    if (input$checkbox == TRUE){
       Ped <- Ped[complete.cases(Ped),]}
    else { # Pedigree remains the same
    }
    
    incProgress(amount = 8/10, detail = "Done")
## output$assignment <- renderTable({Ped[complete.cases(Ped), ] })
    
    
    output$assignment <- DT::renderDT({
      datatable(Ped)   %>%
        DT::formatStyle(
          columns = 1:ncol(Ped),  # Apply style to all columns
          backgroundColor = "rgba(245, 245, 220, 1)"
        )
       })
  
    
    incProgress(amount = 10/10, detail = "Done") #Progress Meter
  
  
  
  #Download Button for Sequoia Package
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$Geno, "Parentage", ".csv", sep = "")
         },
    content = function(file) {
       write.csv(Ped, file)
       }
          )
  
  incProgress(amount = 10/10, detail = "Done")
    })

    ## Plot Pedigree- **Needs to befixed ##  
  #                output$pedigree <- renderPlot({
  # if (input$radio == 1){
  #   plot <- drawPedigree(Ped, dots = "y", dotSize = .005)
  # }
  # if (input$radio == 2){
  #   # Plot the pedegree
  #   ped2 <- fixPedigree(Ped)
  #   # Order the id columns the same
  #   ped2 = ped2[order(ped2$id),]
  #   LHis2 = LHis2[order(LHis2$ID),]
  #   # Assign the Life History sex to the pedigree 
  #   ped2$Sex = LHis2$Sex
  #   
  #   ped2 <- pedigree(id = ped2$id, dadid = ped2$sire, momid = ped2$dam, sex = LHis2$Sex)
  #   plot <- plot(ped2)
  # }
  # plot
  # })
 
  ############ POPULATION GENETIC EVALUATION ################
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
                                       the population data frame:")
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
    pca <- glPca(gl.object, nf = 4)
   pnw.dapc <- dapc(gl.object, n.pca = 4, n.da = nlevels(gl.object$pop)-1)
   
   # Optimal number of PCs based off of a-score
   optim.num <- optim.a.score(pnw.dapc)
   optim.num$best
   # Number of eigenvalues to retain to account for 80% of the variation
   w <- which(cumsum(100*pca$eig/sum(pca$eig)) >= 80)
   
  # PC Radio buttons
   if (input$PCradio == 1){
     pc.num <- 4
   }
   if (input$PCradio == 2){
     pc.num <- optim.num$best
    }
   if (input$PCradio == 3){
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
  PCAPlot <- ggplot(pca.scores, aes(x=PC1, y=PC2, colour=pop)) + 
     geom_point(size=2) + 
     stat_ellipse(level = 0.95, size = 1) + 
     geom_hline(yintercept = 0) + 
     geom_vline(xintercept = 0)
  output$PCA <- renderPlot(PCAPlot)
  
   ## Download the PCA Plot
   output$DownloadPCA <- downloadHandler(
     filename = function() { paste("PCA_plot") },
     content = function(file) {
       pdf(file) #open pdf device
       PCAPlot()
       dev.off() # turn the device off
     }
   )
   
   # Create a compoplot or Structure Plot
   pnw.dapc <- dapc(gl.object, n.pca = pc.num, n.da = nlevels(gl.object$pop)-1)
   
   dapc.results <- as.data.frame(pnw.dapc$posterior)
   dapc.results$pop <- pop(gl.object)
   dapc.results$indNames <- rownames(dapc.results)
   
   dapc.results <- melt(dapc.results)
   
   colnames(dapc.results) <- c("Original_Population","Sample","Assigned_Population","Posterior_membership_probability")
   
   output$Structure <- renderPlot(ggplot(dapc.results, aes(x=Sample,
                                    y=Posterior_membership_probability,
                                    fill=Assigned_Population)) +
                                    geom_bar(stat='identity') +
                                    facet_grid(~Original_Population, scales = "free") +
                                    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)))

    # Download Structure-Like Plot 
      output$downloadStruct <- downloadHandler(
         filename = function(){paste("DAPC_plot", '.pdf', sep = '')},
            # content is a function with argument file. content writes the plot to the device
         content = function(file) {
            cairo_pdf(filename = file, width = 18, height = 6) # open the pdf device
               print(ggplot(dapc.results, aes(x=Sample,
                                           y=Posterior_membership_probability,
                                           fill=Assigned_Population)) +
                     geom_bar(stat='identity') +
                     facet_grid(~Original_Population, scales = "free") +
                     theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))) # for GGPLOT
            dev.off()  # turn the device off
         }, 
         contentType = "application/pdf"
      )
      
      
         ##### POP-GEN Clustering Plots ######   
   # Set the Max & Min K and number of PCAs to be retained
   
   maxK <- input$n
   minK <- 2
   n.pca <- input$pca
   myMat <- matrix(nrow=10, ncol=maxK)
   colnames(myMat) <- 1:ncol(myMat)
   for(i in 1:nrow(myMat)){
      grp <- find.clusters(gl.object, n.pca = n.pca, choose.n.clust = FALSE,  max.n.clust = maxK)
      myMat[i,] <- grp$Kstat
   }

   my_df <- melt(myMat)
   colnames(my_df)[1:3] <- c("Group", "K", "BIC")
   my_df$K <- as.factor(my_df$K)

   # BIC vs K plot
   p1 <- ggplot(my_df, aes(x = K, y = BIC)) +
      geom_boxplot() +
      theme_bw() +
      xlab("Number of groups (K)") +
      ylab("Bayesian Information Criterion (BIC)")

   # Plot a scatter plot for the Max K value
   my_k <- minK:maxK
   grp_l <- vector(mode = "list", length = length(my_k))
   dapc_l <- vector(mode = "list", length = length(my_k))

   for(i in 1:length(dapc_l)){
      set.seed(9)
      grp_l[[i]] <- find.clusters(gl.object, n.pca = n.pca, n.clust = my_k[i])
      dapc_l[[i]] <- dapc(gl.object, pop = grp_l[[i]]$grp, n.pca = n.pca, n.da = my_k[i])
   }

   my_df <- as.data.frame(dapc_l[[ length(dapc_l) ]]$ind.coord)
   my_df$Group <- dapc_l[[ length(dapc_l) ]]$grp
   my_pal <- RColorBrewer::brewer.pal(n=8, name = "Dark2")

   p2 <- ggplot(my_df, aes(x = LD1, y = LD2, color = Group, fill = Group)) +
      geom_point(size = 4, shape = 21) +
      theme_bw() +
      scale_color_manual(values=c(my_pal))+
      scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))

   # Set up the K-means plot
   tmp <- as.data.frame(dapc_l[[1]]$posterior)
   tmp$K <- my_k[1]
   tmp$Sample <- rownames(tmp)
   tmp <- melt(tmp, id = c("Sample", "K"))
   names(tmp)[3:4] <- c("Group", "Posterior")
   tmp$Population <- pop.data$Population
   my_df <- tmp
   for(i in 2:length(dapc_l)){
      tmp <- as.data.frame(dapc_l[[i]]$posterior)
      tmp$K <- my_k[i]
      tmp$Sample <- rownames(tmp)
      tmp <- melt(tmp, id = c("Sample", "K"))
      names(tmp)[3:4] <- c("Group", "Posterior")
      tmp$Population <- pop.data$Population

      my_df <- rbind(my_df, tmp)
   }

   # Plot the K-means Clustering plot #
   grp.labs <- paste("K =", my_k)
   names(grp.labs) <- my_k

   p3 <- ggplot(my_df, aes(x = Sample, y = Posterior, fill = Group)) +
      geom_bar(stat = "identity") +
      facet_grid(K ~ Population , scales = "free_x", space = "free",
                 labeller = labeller(K = grp.labs)) +
      theme_bw() +
      ylab("Posterior membership probability") +
      theme(legend.position='none') +
      #scale_color_brewer(palette="Dark2") +
      scale_fill_manual(values=c(my_pal)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))

   ## Multi-panel Plot
   # p1= K vs BIC plot, p2 = scatter plot, p3= K-means/clustering plot
   output$PopGenPlots <- renderPlot(plot_grid(plot_grid(p1, p2,labels = c("A", "B"), ncol = 2),
             p3, labels = c("", "C"), ncol = 1))

   ## Download the Cluster Plot
   output$DownloadClust <- downloadHandler(
      filename = function() { paste("Cluster_plot") },
      content = function(file) {
         ggsave(file, plot = plotInput(), device = "png")
      }
   )


 })
 
}

# Run the application 
shinyApp(ui = ui, server = server)

