#Experiment: Gene Regulatory Networks of Developmental Plasticity
#Last updated: June 27, 2024
#https://julie-jung.shinyapps.io/grn-shiny/

library(shiny)
library(ggplot2)
library(data.table)

# My functions

# Function to produce summary statistics (mean and +/- se)
se <- function(x) sqrt(var(x)/length(x))

data_summary <- function(x) {
  m <- mean(x, na.rm=T)
  ymin <- m-se(x)
  ymax <- m+se(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

#Read in data

#read in our csv file with the expression values from DEseq
DEseq <- read.csv("values_All.csv", stringsAsFactors=TRUE)
str(DEseq)
#we see that we have a dataframe of 24,678 genes with 173 variables (expression)
#ultimately we want to be able to search (use input text widget and search for an exact string match) among the 24678 genes
#The 173 variables are a mix of condition_development_replicate
#The main challenge here will be to tidy the data from the input in the right long format. 

#First I want to convert the values in the first column into row names in the existing dataframe
DEseq2<-DEseq[,-1]
rownames(DEseq2) <- DEseq[,1]

#Let's start by transposing the dataframe
t_DEseq<-as.data.frame(t(DEseq2))

#now, split that string into condition_timepoint_replicate
suppressMessages(library(tidyverse))

rn <- rownames(t_DEseq) #rn lists the conditions, dev hours, and replicates
#use strsplit to get the output into a 'list', rbind the output to get a matrix "m1"
m1 <- do.call(rbind, strsplit(rn, "_"))
m1 #now they're split
t_DEseq$condition <- m1[,1] #assigns condition as new column
t_DEseq$development <- m1[,2]
t_DEseq$replicate <- m1[,3]

# Define UI -----
# Peruse widgets from here: https://shiny.posit.co/r/gallery/widgets/widget-gallery/
# More/extension of widgets: https://dreamrs.github.io/shinyWidgets/

ui <- fluidPage(
  
  # App title -----
  titlePanel("Search gene expression values from DEseq"),
  
  # Sidebar layout with input and output definitions -----
  sidebarLayout(position="left",
    
    # Sidebar panel for inputs -----
    sidebarPanel( #"sidebar panel", #insert name here if you want
      
      # Text input to subset by gene (must be an exact match)
      textInput(inputId = "geneID", #the name to use to look up the value of the widget (as a character string)
                label = "Gene:", #a label to display above the text field
                value = "Contig0-snapTAU.703"), #placeholder
      
      # Choose gene from checkboxes of possible genes (instead of providing an exact match)
      checkboxGroupInput(inputId = "Question", 
                         label = "Conditions: ", 
                         choices = c("Agar WT" = "AG", #a list of values. The widget will include a menu option for each value of the list. If the list has names, these will be displayed in the drop down menu. Otherwise the values themselves will be displayed. 
                                     "eud1 k/o" = "eud1", 
                                     "eud1 t/g" = "eud1TG", 
                                     "Expectatus WT" = "Pex", 
                                     "Liquid Culture WT" = "LC", 
                                     "lsy12 k/o" = "lsy12", 
                                     "mbd2 k/o" = "mbd2", 
                                     "nag1 and nag2 k/o" = "nag12", 
                                     "nhr40 k/o" = "nhr40", 
                                     "RS5410 St bias" = "RS5410", 
                                     "RS5427 Eu bias" = "RS5427",
                                     "RSA100 Eu bias" = "RSA100",
                                     "RSC017 St bias" = "RSC017", 
                                     "sult1 k/o" = "sult1"), 
                         selected="AG")
    ), # This closes out the sidebar panel. 
    
    # Main panel for displaying outputs -----
    
    mainPanel(
      uiOutput('ui_plot'), 
      downloadButton("downloadPlot", "Download")
    ) # This closes out the main panel. 
  
    )) # This closes out the UI

# Define server logic

server <- function(input, output) {

  output$ui_plot <- renderUI({
    
    out <- list()
    if (length(input$Question)==0){return(NULL)}
    for (i in 1:length(input$Question)){
      out[[i]] <-  plotOutput(outputId = paste0("plot",i))
    }  
    return(out) 
    
  })
  
  #render plots
  CountPlotFunction <- function(MyData) {
    MyPlot <- ggplot(data = MyData, aes(x=development, y=geneID)) +
      geom_jitter(pch=16, size=2.5, width=0.1, na.rm=T) +
      #stat_summary(fun.data=data_summary, color='gray')+
      labs(x = "development", y = "gene expression") +
      theme_bw(base_size=12, base_family="Palatino")
    return(MyPlot)
  }
  
  
  observe({  
    for (i in 1:14){  
      local({  #because expressions are evaluated at app init
        ii <- i 
        output[[paste0('plot',ii)]] <- renderPlot({ 
          
          if ( length(input$Question) > ii-1 ){
            #I want MyData to be a subset of t_DEseq where the rows include input$Question[[ii]]
            #do the rownames include the condition in question?? if so keep those
            subset_conditions<-subset(t_DEseq, t_DEseq$condition==input$Question[[ii]])
            #subset_conditions$development shows each column as a geneID! and each row as a condition
            #let's make sure the development column is there/included! 
            rn <- rownames(subset_conditions) #rn lists the conditions, dev hours, and replicates
            #use strsplit to get the output into a 'list', rbind the output to get a matrix "m1"
            m1 <- do.call(rbind, strsplit(rn, "_"))
            m1 #now they're split
            subset_conditions$condition <- m1[,1] #assigns condition as new column
            subset_conditions$development <- m1[,2]
            subset_conditions$replicate <- m1[,3]
            #str(subset_conditions) to check
            subset_temp2<-subset(subset_conditions, select=c(input$geneID, "development"))
            MyPlot <- ggplot(data = subset_temp2, aes(x=development, y=subset_temp2[,1])) +
              geom_jitter(pch=16, size=2.5, width=0.1, na.rm=T) +
              stat_summary(fun.data=data_summary, color='gray')+
              labs(x = "development", y = "gene expression") +
              theme_bw(base_size=12, base_family="Palatino")
            return(MyPlot)
          } 
          NULL
        })
      })
    }                                  
  })
  
  # Download button for plot
 output$downloadPlot <- downloadHandler(
   filename = function() { paste(input$dataset, '.png', sep='') },
   content = function(file) {
     ggsave(file, plot = plotInput(), device = "png")
   }
 )
}

# Create Shiny app
shinyApp(ui = ui, server = server)
