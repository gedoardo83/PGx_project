#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(ggrepel)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Prediction score for SERPINA1"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        selectInput(inputId = "myscore", 
                    label = "Choose a score to display",
                    choices = c(
                      "REVEL" = "REVEL_score", 
                      "PolyPhen2 HVAR" = "Polyphen2_HVAR_score", 
                      "SIFT" = "SIFT_score",
                      "MetaSVM" = "MetaSVM_score"),
                    selected = "REVEL"),
        checkboxInput("VarLabels", "Show known vars labels", value = FALSE),
        checkboxInput("usemyvar", "Use my variant", value = FALSE),
        textInput("user_var", "Input your var as genomic (chr_pos_ref_alt)\nor protein (A25P):", value = "Enter your var here...")
        ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("distPlot"),
         plotOutput("jitterPlot"),
         plotOutput("AFPlot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  output$distPlot <- renderPlot({
    if (input$VarLabels == FALSE) {
      # generate bins based on input$bins from ui.R
      
      ggplot(data=SERPINA1, aes(x=Known_vars,y=SERPINA1[,input$myscore])) + geom_violin(scale = "width") + geom_boxplot(width=0.2)
    } else {
      mysubset <- SERPINA1[SERPINA1$Known_vars=="P" | SERPINA1$Known_vars=="P*",]
      ggplot() + geom_violin(data=SERPINA1, aes(x=Known_vars,y=SERPINA1[,input$myscore]), scale = "width") + geom_boxplot(width=0.2) + geom_point(data=mysubset, aes(x=Known_vars, y=mysubset[,input$myscore])) + geom_label_repel(data=mysubset, aes(x=Known_vars, y=mysubset[,input$myscore], label=var_name))
    }
  })
  
  output$jitterPlot <- renderPlot({
    if (input$usemyvar == TRUE) {
      myvartab <- SERPINA1[SERPINA1$var_id == input$user_var | SERPINA1$aachange == input$user_var,]
      ggplot(data=SERPINA1, aes(x=U_notU,y=SERPINA1[,input$myscore],color=Known_vars)) + geom_jitter() + geom_point(data=myvartab,aes(x=U_notU,y=myvartab[,input$myscore]),color="black", size=5) + geom_label_repel(data=myvartab,aes(x=U_notU,y=myvartab[,input$myscore], label=paste0(var_id,"\n",aachange)),point.padding = 0.5, show.legend=FALSE, color="black")
    } else {
      ggplot(data=SERPINA1, aes(x=U_notU,y=SERPINA1[,input$myscore],color=Known_vars)) + geom_jitter()
    }
  })
  
  output$AFPlot <- renderPlot({
    ggplot(data=SERPINA1, aes(x=-log10(ExAC_Adj_AF),y=SERPINA1[,input$myscore], color=Known_vars)) + geom_point() 
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

