library(shiny)
source("LabTestScript.R")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Virtual Metal Salts and Acid Molarities for Chemistry 1 Lab Test"),

    # Sidebar with...
    sidebarLayout(
        sidebarPanel(
            column(6, 
                   numericInput("s1m", 
                                h6("Salt 1 Mass (g):"), 
                                value = 1.1)),  
            column(6, 
                   numericInput("s1v", 
                                h6("Salt 1 Acid Titration Volume (mL):"), 
                                value = 1.8)),  
            column(6, 
                   numericInput("s2m", 
                                h6("Salt 2 Mass (g):"), 
                                value = 24.57)),  
            column(6, 
                   numericInput("s2v", 
                                h6("Salt 2 Acid Titration Volume (mL):"), 
                                value = 27.11))  
        ),

        # Show a plot of the generated distribution
        mainPanel(
           tableOutput("SaltPairs")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$SaltPairs <- renderTable(
        fullcalculation(input$s1m,input$s1v,input$s2m,input$s2v), striped = FALSE, hover = FALSE, bordered = FALSE,
        spacing = c("s"), width = "auto", align = NULL,
        rownames = FALSE, colnames = TRUE, digits = 4, na = "NA",
        env = parent.frame(), quoted = FALSE, outputArgs = list())
    
}

# Run the application 
shinyApp(ui = ui, server = server)
