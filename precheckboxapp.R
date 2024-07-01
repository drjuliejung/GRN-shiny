#Experiment: Gene Regulatory Networks of Developmental Plasticity
#Last updated: June 24, 2024
#https://julie-jung.shinyapps.io/grn-shiny/

library(shiny)
library(ggplot2)
library(data.table)

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
      
      #checkboxInput("do2", "Make 2 plots", value = T), 
      
      # Drop down selection for conditions
      selectInput(inputId = "conditionInput", #name to use to look up the value of the widget (as a character string)
                  label = "Condition:", #label to display above the drop-down box
                  choices = c("Agar WT" = "AG", #a list of values. The widget will include a menu option for each value of the list. If the list has names, these will be displayed in the drop down menu. Otherwise the values themselves will be displayed. 
                              "eud1 k/o" = "eud1", 
                              "eud1 t/g" = "eud1TG", 
                              "Expectatus WT" = "Exp", 
                              "Liquid Culture WT" = "LC", 
                              "lsy12 k/o" = "lsy12", 
                              "mbd2 k/o" = "mbd2", 
                              "nag1 and nag2 k/o" = "nag12", 
                              "nhr40 k/o" = "nhr40", 
                              "RS5410 St bias" = "RS5410", 
                              "RS5427 Eu bias" = "RS5427",
                              "RSA100 Eu bias" = "RSA100",
                              "RSC017 St bias" = "RSC017", 
                              "sult1 k/o" = "sult1")), 
      
      # Text input to subset by gene (must be an exact match)
      textInput(inputId = "geneID", #the name to use to look up the value of the widget (as a character string)
                label = "Gene:", #a label to display above the text field
                value = "Contig0-snapTAU.703"), #placeholder
      
      # ALTERNATIVE: Choose gene from checkboxes of possible genes (instead of providing an exact match)
      
      ## Text input to subset by a second gene (must be an exact match)
      #textInput(inputId = "geneID2", #the name to use to look up the value of the widget (as a character string)
      #          label = "Another Gene:", #a label to display above the text field
      #          value = "Contig13-snapTAU.388"), #placeholder
      
      # Drop down selection for a second condition
      selectInput(inputId = "conditionInput2", #name to use to look up the value of the widget (as a character string)
                  label = "Another Condition:", #label to display above the drop-down box
                  choices = c("Agar WT" = "AG", #a list of values. The widget will include a menu option for each value of the list. If the list has names, these will be displayed in the drop down menu. Otherwise the values themselves will be displayed.
                              "eud1 k/o" = "eud1",
                              "eud1 t/g" = "eud1TG",
                              "Expectatus WT" = "Exp",
                              "Liquid Culture WT" = "LC",
                              "lsy12 k/o" = "lsy12",
                              "mbd2 k/o" = "mbd2",
                              "nag1 and nag2 k/o" = "nag12",
                              "nhr40 k/o" = "nhr40",
                              "RS5410 St bias" = "RS5410",
                              "RS5427 Eu bias" = "RS5427",
                              "RSA100 Eu bias" = "RSA100",
                              "RSC017 St bias" = "RSC017",
                              "sult1 k/o" = "sult1"))
      
    ), # This closes out the sidebar panel. 
    
    # Main panel for displaying outputs -----
    
    mainPanel( #"main panel", 
              
              #fluidRow(
              #  splitLayout(cellWidths = c("50%", "50%"), plotOutput("boxplot1", plotOutput("boxplot2"))
              #  )), #close fluidRow
      
     plotOutput("boxplot"), #delete this comma and on
    
    # I think this is also where we'd want to put our download button
    downloadButton("downloadPlot", "Download")
    #fluidPage(downloadButton('downloadPlot'))
    
    ) # This closes out the main panel. 
  
    )) # This closes out the UI

# Define server logic

# Function to produce summary statistics (mean and +/- se)
se <- function(x) sqrt(var(x)/length(x))

data_summary <- function(x) {
  m <- mean(x, na.rm=T)
  ymin <- m-se(x)
  ymax <- m+se(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

#server

server <- function(input, output) {
    
  
  output$boxplot <- renderPlot({
    
    
    #you can access the value of the select widget with input$select (or in our case, input$conditionInput)  
    subset.DEseq <- subset(t_DEseq, t_DEseq$condition == input$conditionInput) 
    # still have all columns
    
    #separate out the development column so that we don't lose it. 
    development<-subset.DEseq$development
    
    #subset the gene too! currently column name. but without losing development column
    keeps <- c(colnames(subset.DEseq) %in% input$geneID) #specify which columns to keep in a vector
    
    subset.DEseq<- cbind(subset.DEseq[keeps], development)
    
    colnames(subset.DEseq)<-c("geneID", "development")
    
##second conditions input
#    subset.DEseq <- subset(t_DEseq, t_DEseq$condition == input$conditionInput | t_DEseq$condition == input$conditionInput2)
#    development_cond2<-subset.DEseq_cond2$development
#    keeps_cond2 <- c(colnames(subset.DEseq_cond2) %in% input$geneID)
#    subset.DEseq_cond2<- cbind(subset.DEseq_cond2[keeps_cond2], development)
#    colnames(subset.DEseq_cond2)<-c("geneID", "development")

#    str(subset.DEseq_cond2)
    
    # draw the plot of expression over development
    boxplot <- ggplot(data=subset.DEseq, aes(x=development, y=geneID))+
      geom_jitter(pch=16, size= 2.5, width=0.1, na.rm=T) +
      stat_summary(fun.data=data_summary, color='gray')+
      
      #geom_jitter(data= subset.DEseq_cond2, pch=16, size= 2.5, width=0.1, na.rm=T, color="red") +
      #stat_summary(fun.data=data_summary, color='pink')+
      
      labs(x = "development",
           y = "gene expression") +
      theme_bw(base_size=12, base_family="Palatino")
    boxplot
    #plotly::ggplotly(expression_plot)
    

    # 
    # # draw the plot of expression over development
    # boxplot2 <- ggplot(data=subset.DEseq_cond2, aes(x=development, y=geneID))+
    #   geom_jitter(pch=16, size= 2.5, width=0.1, na.rm=T) +
    #   stat_summary(fun.data=data_summary, color='gray')+
    #   labs(x = "development",
    #        y = "gene expression") +
    #   theme_bw(base_size=12, base_family="Palatino")
    # boxplot2
    # #plotly::ggplotly(expression_plot)
      
    }) 
    #output$boxplot2 <- renderPlot({boxplot2()})
  
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
