#install.packages("shiny")
#install.packages('visNetwork')

library(shiny)
library(shinyalert)
library(shinyjs)
library(reticulate)
library(uuid)
library(ggplot2)
library(RSQLite)

source("functions_lib.R")

require(visNetwork)

isValidEmail <- function(x) {
    grepl("\\<[A-Z0-9._%+-]+@[A-Z0-9.-]+\\.[A-Z]{2,}\\>", as.character(x), ignore.case=TRUE)
}

# Define UI for dataset viewer app ----
ui <- fluidPage(
    useShinyjs(),
  # App title ----
  titlePanel("Relationship among Biological Entities"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      fileInput("file1", "Choose txt Uniprot Proteins File",
                multiple = FALSE,
                accept = c(".txt")),
                                      
      sliderInput("tanimoto", "Tanimoto score:",
                  min = 0.5, max = 1,
                  value = 0.5, step = 0.1),
      
      checkboxInput("liginfer", "Enable Ligands Inference", TRUE),
      
      numericInput("num_contacts", "Number of Ligand Contacts:", 1, min = 1),
      
      numericInput("perc_contacts", "Percentage of conserved AAs in domain:", 10, min = 1, max = 100),
      
      textInput("email", label = "E-mail", value = "user@site.com"),
      hidden(p(id = "notice_error", "Use a valid e-mail")),
      
      # Input: Select a dataset ----
      #selectInput("dataset", "Choose Node Type:",
      #            choices = c("Protein", "Structure", "Family", "Ligand")),
      
      # Input: Specify the number of observations to view ----
      #numericInput("obs", "Number of observations to view:", 10),

      # Include clarifying text ----
      helpText("Note: a job will be created and a token will be given to you to retrieve the results."),
      hidden(p(id = "notice", "Processing...")),

      # Input: actionButton() to defer the rendering of output ----
      # until the user explicitly clicks the button (rather than
      # doing it immediately when inputs change). This is useful if
      # the computations required to render output are inordinately
      # time-consuming.
      actionButton("update", "Build Network")

    ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Header + summary of distribution ----
      #h4("Summary"),
      #verbatimTextOutput("summary"),

      # Output: Header + table of distribution ----
      #h4("Observations"),
      #tableOutput("view"),
      
      hidden( selectInput("protein_selection",  "Select Protein", choices = c('All'), selected="All" ) ),
      
      h3("Network"),
      visNetworkOutput("bionet"),
      textOutput("log"),
      
      h3("Statistics"),
      plotOutput("plot", click = "plot_click"),
    )

  )
)

# Define server logic to summarize and view selected dataset ----
server <- function(input, output, session) {

  # Return the requested dataset ----
  # Note that we use eventReactive() here, which depends on
  # input$update (the action button), so that the output is only
  # updated when the user clicks the button
  dataset <- rock

  # Generate a summary of the dataset ----
  output$summary <- renderPrint({
    summary(dataset)
  })
  
  observe( {
      toggle("num_contacts", condition=input$liginfer)
      toggle("perc_contacts", condition=input$liginfer)
      toggle("email", condition=input$liginfer)
  })
  
  # Show the first "n" observations ----
  # The use of isolate() is necessary because we don't want the table
  # to update whenever input$obs changes (only when the user clicks
  # the action button)
  output$view <- renderTable({
    head(dataset, n = isolate(10))
  })
  
  msg=""
  output$log <- renderText({msg})
  
  #file.remove("bionetwork.db")
  conn =  dbConnect(RSQLite::SQLite(), "bionetwork.db")
  init_db(conn)
  
  graph=NULL
  nameFile=NULL
  observeEvent(input$update, {
      inFile <- input$file1
      msg="Reading file..."
      output$log <- renderText({msg})
      if (is.null(inFile)){
          shinyalert("Error", "Upload proteins file", type = "error")
          #return(NULL)
      }
      else{
          nameFile=inFile$datapath
          proteins = read_lines(nameFile)
          if(length(proteins)>0 && length(proteins)<=100 ){
              print(input$email)
              flag=TRUE
              if(input$liginfer){
                  if( !isValidEmail(input$email) || input$email=='' || input$email=='user@site.com' ){
                      shinyalert(title="Error", text="Email Address: Please Input a valid E-mail address", type = "error")
                      shinyjs::show("notice_error")
                      flag=FALSE
                  }
              }
              
              if(flag){
                  idu<-UUIDgenerate(use.time = NA, n = 1L)
                  file.copy(nameFile, paste0(idu,".txt") )
                  nameFile=paste0(idu,".txt")
                  
                  shinyjs::disable("update")
                  shinyjs::show("notice")
                  shinyjs::hide("notice_error")
                  
                  if( is.null(conn)){
                      conn =  dbConnect(RSQLite::SQLite(), "bionetwork.db") 
                      init_db(conn)
                  }
                  
                  msg=paste0(msg, "<br>Loading proteins information...")
                  output$log <- renderText({msg})
                  process_protein_list(nameFile, conn)
                  
                  msg=paste0(msg, "<br>Calculating edges and nodes...")
                  output$log <- renderText({msg})
                  
                  #separate function to generate node sand edges from the function that renders the plot
                  #dir.create( idu, showWarnings = FALSE)
                  
                  #idu="44aff6c3-0f96-45e6-9492-dfb2a2dc3430"
                  #nameFile="../proteins.txt"
                  #process_protein_list(nameFile, conn)
                  #data=get_network(conn, nameFile, 0.6, idu)
                  
                  data=get_network(conn, nameFile, input$tanimoto, idu, 'All')
                  nodes=data[1][[1]]
                  edges=data[2][[1]]
                  
                  print( paste(length(nodes), length(edges)) )
                  msg=paste0(msg, "<br>Rendering network")
                  output$log <- renderText({msg})
                  
                  index=length(proteins)
                  proteins[index+1]='All'
                  updateSelectInput(session, input="protein_selection",  "Select Protein", choices = proteins, selected="All" )
                    
                  shinyjs::show("protein_selection")
                  
                  data=get_network(conn, nameFile, input$tanimoto, idu, input$protein_selection)
                  nodes=data[1][[1]]
                  edges=data[2][[1]]
                  output$bionet <- renderVisNetwork({
                      
                      render_network(nodes, edges)
                      
                  })
                  
                  if(input$liginfer){
                      system(paste0("python3 back.py ", nameFile, " ", input$num_contacts, " ", input$perc_contacts, " ", idu, " ", input$email))
                  }
                  
                  netproxy <- visNetworkProxy("bionet")
                  
                  x=c()
                  y=c()
                  i=1
                  for(t in unique(nodes$type) ){
                    x[i]=""
                    y[i]=0
                    i=i+1
                  }
                  #x=unique(nodes$type)
                  i=1
                  for(t in nodes$type_id ){
                      y[t]=y[t]+1
                      x[t]=as.character(nodes$type[i])
                      i=i+1
                  }
                  print(x)
                  output$plot <- renderPlot({
                    if(input$protein_selection=="All"){
                            
                          df=data.frame(type=x, occurrences=y)
                          
                          ggplot(data=df, aes(x=type, y=occurrences, fill=type) ) + geom_bar(stat="identity") + geom_text(aes(label=y), vjust=0.1) 
                    }
                    else{
                        
                    }
                  }, res = 96)
                  
                  observeEvent(input$protein_selection, {
                      data=get_network(conn, nameFile, input$tanimoto, idu, input$protein_selection)
                      nodes=data[1][[1]]
                      edges=data[2][[1]]
                      output$bionet <- renderVisNetwork({
                          
                          render_network(nodes, edges)
                          
                      })
                      netproxy <- visNetworkProxy("bionet")
                      
                      x=c()
                      y=c()
                      i=1
                      for(t in unique(nodes$type) ){
                        x[i]=""
                        y[i]=0
                        i=i+1
                      }
                      #x=unique(nodes$type)
                      i=1
                      for(t in nodes$type_id ){
                          y[t]=y[t]+1
                          x[t]=as.character(nodes$type[i])
                          i=i+1
                      }
                      print(x)
                      output$plot <- renderPlot({
                        df=data.frame(type=x, occurrences=y)
                        ggplot(data=df, aes(x=type, y=occurrences, fill=type) ) + geom_bar(stat="identity") + geom_text(aes(label=y), vjust=0.1) 
                      }, res = 96)
                  })
                  
                  shinyjs::enable("update")
                  shinyjs::hide("notice")
              }
          }
          else{
              shinyalert("Error", "The protein seeds must have more than 1 and maximum of 100", type = "error")
          }
      }
  })
  

  
  observe ({
      output$log <- renderText({msg})
      if(! is.null(nameFile)){
          get_network(conn, nameFile, input$tanimoto)
      }
  })
  
  output$contents <- renderTable({
  
  })
  
  #process_protein_list(inFile$datapath, conn)
  
  # deal with file: https://shiny.rstudio.com/reference/shiny/0.14/fileInput.html
  
  

}

# Create Shiny app ----
shinyApp(ui, server)

# runApp()
