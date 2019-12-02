library(shiny)
library(sequoia)
library(Rcpp)
library(pedantics)
library(kinship2)
library(shinycssloaders)

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
        
        fileInput("LHis", label = h3("Life History File input")),
        
        actionButton("run", "Run"),
        
        radioButtons("radio", label = h3("Plot Pedigree"),
                     choices = list("Pedigree 1" = 1, "Pedigree 2" = 2), 
                     selected = 1)
      ),
    
      # Show the pedigree
      mainPanel(
          
        # Output: Tabset
        tabsetPanel(type = "tabs",
                    tabPanel("Parentage Assignment",# Button
                             withSpinner(downloadButton("downloadData", "Download")),
                             tableOutput("assignment")
                             ),
                    tabPanel("Pedigree Plot", plotOutput("pedigree"))
                    
        )
      )
   )
)



########### SERVER #############
server <- function(input, output) {
  
  output$myImage <- renderImage({
    # A temp file to save the output.
    # This file will be removed later by renderImage
    outfile <- tempfile(fileext = '.png')
    
    # Generate the PNG
    png(outfile, width = 400, height = 300)
    hist(rnorm(input$obs), main = "Generated in renderImage()")
    dev.off()
  })
  
 observeEvent(input$run, {
    # Read in the Life History File
   LHisFile <- input$LHis
   if (is.null(LHisFile))
     return(NULL)
   LHis <- read.csv(LHisFile$datapath)
    LHis2 <- LHis  

    # Read in the Genotype File
   GenoFile <- input$Geno
   if (is.null(GenoFile))
     return(NULL)
   Geno <- read.csv(GenoFile$datapath)
   row.names(Geno) <- Geno$X
   Geno <- Geno[,-1]
   Geno <- as.matrix(Geno)

   levels(LHis$Sex) <- sub("female", "1", levels(LHis$Sex))
   levels(LHis$Sex) <- sub("male", "2", levels(LHis$Sex))
    ParOUT <- sequoia(GenoM = Geno,
                    LifeHistData = LHis,
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
}

# Run the application 
shinyApp(ui = ui, server = server)

