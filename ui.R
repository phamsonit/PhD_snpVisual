library(shiny)
library(shinythemes)
library(cluster)
library(shape)
library(ggplot2)
options(shiny.maxRequestSize = 50*1024^2)

# global variables ... 
max_cluster_height <- 0

#hierarchical clustering methods
cmeths = c("ward.D", "ward.D2",
           "single", "complete", "average", "mcquitty",
           "median", "centroid")
#distance measures
dmeths = c("euclidean", "maximum", "manhattan", "canberra")

#   main shiny components: ui and server
#   ui: defines page layout and components
#   server: defines operations
fluidPage(
  #shinythemes::themeSelector(),
  #theme = shinytheme("cerulean"),
  #theme = shinytheme("lumen"),
  
  titlePanel("Genetic variant visualization"),
  
  #define input parameters
  sidebarPanel(
    #Input data
    #pattern file
    fileInput('datafile', 'Choose pattern file:',
              accept=c('text/csv', 'text/comma-separated-values,text/plain')),
    #index file
    fileInput('indexfile', 'Choose index file:',
              accept=c('text/csv', 'text/comma-separated-values,text/plain')),
    #chromosome file
    fileInput('chrfile', 'Choose chromosome file:',
              accept=c('text/csv', 'text/comma-separated-values,text/plain')),
    #gene file
    fileInput('genefile', 'Choose genes file:',
              accept=c('text/csv', 'text/comma-separated-values,text/plain')),
    #QTL file
    fileInput('QTLfile', 'Choose QTL file:',
              accept=c('text/csv', 'text/comma-separated-values,text/plain')),
    
    #Clustering methods and distance measures input
    fluidRow(
      selectInput("meth", "Select clustering method:", choices=cmeths,
                  selected=cmeths[2])),
    fluidRow(
      selectInput("dmeth", "Select distance measure:", choices=dmeths,
                  selected=dmeths[1])),
    
    #Cut threshold value and method to create pattern groups
    fluidRow(
       uiOutput("cutval")),
    
    fluidRow(
      radioButtons("mp_method","Choose create pattern groups method:",
      c("Union" = "union", "Intersection"="intersection","Majority"="majority")
    )),
    
    fluidRow(
      sliderInput("fre", "Frequency", value = 0.5, min = 0.1, max = 0.9 )
    ),
    
    #Select pattern groups to visualize
    fluidRow(
      uiOutput("select_all")),
      uiOutput("mp_checkgroup"),

    width=2),
  
  
  #
  # main panel (graphical presentation)
  #
  mainPanel(
    
    tabsetPanel(
      tabPanel("Main", 
               htmlOutput("main")),
      #show hierarchical clustering tree
      tabPanel("Clustering", 
               plotOutput("cluster", height = 800)),
      #show visualization of SNP in chromosome
      tabPanel("Visualization",
              fluidRow(
                #pattern group plot: overview SNP on chromosomes
                column(width = 6, align="center",
                  plotOutput("meta_plot", height = 800,#width = 800,
                             click = "meta_plot_click", #click to show the position
                             dblclick = "meta_plot_bdlclick", #double click to zoom
                             brush = brushOpts(
                               id = "meta_brush",
                               resetOnNew = FALSE)
                             ),
                  verbatimTextOutput("meta_pos_info") ),
                #zoom plot: focus on a region of chromosome
                column(width = 6,
                       plotOutput("zoom_meta_plot", height = 800,
                                  click = "zoom_meta_plot_click"),
                       verbatimTextOutput("zoom_meta_pos_info"))
                )),
      #show detail of pattern groups
      tabPanel("Detail", 
               dataTableOutput("meta_pattern_detail"))
      
    ),width = 10) #end mainPanel
)  # end fluidPage

