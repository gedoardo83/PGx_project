#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("PGx Report"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         textInput("ClinCode", h3("Clinician code"),
          value = "Enter clinician code..."
          ),
         uiOutput("Sample")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        tabsetPanel(id = "mainarea",
        tabPanel("Summary", 
          h3(textOutput("RepTitle"), align="center"),
          h4(textOutput("firstDrug")),
          DT::dataTableOutput("mygenos") ),
        tabPanel("Drug details",
          uiOutput("Drugs"),
          h3(textOutput("DrugCategory")),
          h4("Recommendations based on CYP"),
          DT::dataTableOutput("DrugCYPDetails"),
          h4("Recommendations based on other SNPs"),
          DT::dataTableOutput("DrugSNPDetails")),
        tabPanel("Genotypes",
          h4("Recommendations based on other SNPs"),
          DT::dataTableOutput("GenoDetails") )
        )
      )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  all_drugs <- scan("AD_ITA.list", what="", sep="\n")
  drugs_category <- reactiveValues(red = character(), yellow = character(), green = character())
  samples_info <- read.table("samples_info.csv", header=T, as.is=T)
  samples_genos <- read.table("samples_genos.csv", header=T, as.is=T)
  guidelines_CYP <- read.table("GuideLines_table_CYP.csv", header=T, sep="\t", as.is=T)
  guidelines_SNP <- read.table("GuideLines_table_SNPs.csv", header=T, sep="\t", as.is=T)
  clinicians <- unique(samples_info$Clinician)
  mysample_df <- data.frame(Green=character(), Yellow=character(), Red=character())
  
  output$Sample <- renderUI({
    selectInput("SampleList", h3("Subject code:"), choices = unique(samples_info$Sample[samples_info$Clinician==input$ClinCode & samples_info$Group_PGx==1]) )
  })
  
  output$RepTitle <- renderText({
    selected_sample <- input$SampleList
    paste("PGx report for sample", selected_sample)
  })
  
  output$firstDrug <- renderText({
    selected_sample <- input$SampleList
    first_drug <- samples_genos$first_drug[samples_genos$Sample==selected_sample]
    paste("1st drug: ", first_drug)
  })
  
  output$Drugs <- renderUI({
    selected_sample <- input$SampleList
    first_drug <- samples_genos$first_drug[samples_genos$Sample==selected_sample]
    mydrugs <- all_drugs[which(all_drugs != first_drug)]
    selectInput("DrugList", h4("Select drug to obtain details:"), choices = mydrugs)
  })
  
  output$DrugCYPDetails <- DT::renderDataTable({
    selected_sample <- input$SampleList
    selected_drug <- input$DrugList
    
    mydetail_CYP2D6<-guidelines_CYP[guidelines_CYP$Drug == selected_drug & guidelines_CYP$Gene == "CYP2D6" & guidelines_CYP$Phenotype == samples_genos$CYP2D6_pheno[samples_genos$Sample == selected_sample],] 
    if (nrow(mydetail_CYP2D6)>0) {
      mydetail_CYP2D6$Genotype <- samples_genos$CYP2D6_alleles[samples_genos$Sample == selected_sample]
    }
    mydetail_CYP2C19<-guidelines_CYP[guidelines_CYP$Drug == selected_drug & guidelines_CYP$Gene == "CYP2C19" & guidelines_CYP$Phenotype == samples_genos$CYP2C19_pheno[samples_genos$Sample == selected_sample],]
    if (nrow(mydetail_CYP2C19) >0) {
      mydetail_CYP2C19$Genotype <- samples_genos$CYP2C19_alleles[samples_genos$Sample == selected_sample]
    }
    mydetail_CYP <- rbind(mydetail_CYP2D6,mydetail_CYP2C19)
    if (nrow(mydetail_CYP)>0) {
      mydetail_CYP[,c(2,8,3:6)]
    } else {mydetail_CYP}
  })
  
  output$DrugSNPDetails <- DT::renderDataTable({
    selected_sample <- input$SampleList
    selected_drug <- input$DrugList
    
    mydetail_SNP <- guidelines_SNP[guidelines_SNP$Drug == selected_drug & guidelines_SNP$Gene == "SLC6A4" & guidelines_SNP$Genotype == samples_genos$SLC6A4[samples_genos$Sample == selected_sample],]
    mydetail_SNP <- rbind(mydetail_SNP, guidelines_SNP[guidelines_SNP$Drug == selected_drug & guidelines_SNP$Variant == "rs489693" & guidelines_SNP$Genotype == samples_genos$rs489693[samples_genos$Sample == selected_sample],])
    mydetail_SNP <- rbind(mydetail_SNP, guidelines_SNP[guidelines_SNP$Drug == selected_drug & guidelines_SNP$Variant == "rs4713916" & guidelines_SNP$Genotype == samples_genos$rs4713916[samples_genos$Sample == selected_sample],])
    mydetail_SNP <- rbind(mydetail_SNP, guidelines_SNP[guidelines_SNP$Drug == selected_drug & guidelines_SNP$Variant == "rs7997012" & guidelines_SNP$Genotype == samples_genos$rs7997012[samples_genos$Sample == selected_sample],])
    mydetail_SNP <- rbind(mydetail_SNP, guidelines_SNP[guidelines_SNP$Drug == selected_drug & guidelines_SNP$Variant == "rs6295" & guidelines_SNP$Genotype == samples_genos$rs6295[samples_genos$Sample == selected_sample],])    
    mydetail_SNP[,2:5]
  })
  
  output$mygenos = DT::renderDataTable({
    selected_sample <- input$SampleList
    first_drug <- samples_genos$first_drug[samples_genos$Sample==selected_sample]
    mydrugs <- all_drugs[which(all_drugs != first_drug)]
    red_drugs <- character()
    yellow_drugs <- character()
    green_drugs <- character()
    for (d in mydrugs) {
      suggestions <- character()
      suggestions <- guidelines_CYP$Category[guidelines_CYP$Drug == d & guidelines_CYP$Gene == "CYP2D6" & guidelines_CYP$Phenotype == samples_genos$CYP2D6_pheno[samples_genos$Sample == selected_sample]]
      suggestions <- c(suggestions, guidelines_CYP$Category[guidelines_CYP$Drug == d & guidelines_CYP$Gene == "CYP2C19" & guidelines_CYP$Phenotype == samples_genos$CYP2C19_pheno[samples_genos$Sample == selected_sample]])
      suggestions <- c(suggestions, guidelines_SNP$Category[guidelines_SNP$Drug == d & guidelines_SNP$Gene == "SLC6A4" & guidelines_SNP$Genotype == samples_genos$SLC6A4[samples_genos$Sample == selected_sample]])
      suggestions <- c(suggestions, guidelines_SNP$Category[guidelines_SNP$Drug == d & guidelines_SNP$Variant == "rs489693" & guidelines_SNP$Genotype == samples_genos$rs489693[samples_genos$Sample == selected_sample]])
      suggestions <- c(suggestions, guidelines_SNP$Category[guidelines_SNP$Drug == d & guidelines_SNP$Variant == "rs4713916" & guidelines_SNP$Genotype == samples_genos$rs4713916[samples_genos$Sample == selected_sample]])
      suggestions <- c(suggestions, guidelines_SNP$Category[guidelines_SNP$Drug == d & guidelines_SNP$Variant == "rs7997012" & guidelines_SNP$Genotype == samples_genos$rs7997012[samples_genos$Sample == selected_sample]])
      suggestions <- c(suggestions, guidelines_SNP$Category[guidelines_SNP$Drug == d & guidelines_SNP$Variant == "rs6295" & guidelines_SNP$Genotype == samples_genos$rs6295[samples_genos$Sample == selected_sample]])
      if ("red" %in% suggestions) {
        red_drugs <- c(red_drugs, d)
      } else if ("yellow" %in% suggestions) {
        yellow_drugs <- c(yellow_drugs, d)
      } else {
        green_drugs <- c(green_drugs, d)
      }
      
    }
    drugs_category$red <- red_drugs
    drugs_category$yellow <- yellow_drugs
    drugs_category$green <- green_drugs
    max.len <- max(length(red_drugs), length(yellow_drugs), length(green_drugs))
    red_drugs <- c(red_drugs, rep(NA, max.len - length(red_drugs)))
    yellow_drugs <- c(yellow_drugs, rep(NA, max.len - length(yellow_drugs)))
    green_drugs <- c(green_drugs, rep(NA, max.len - length(green_drugs)))
    
    mysample_df <- data.frame(Green=green_drugs, Yellow=yellow_drugs, Red=red_drugs)
    mysample_df
    
  }, options=list(pageLength = 25))

  output$DrugCategory <- renderText({
    if (input$DrugList %in% drugs_category$red) {
      mycategory <- "red"
    } else if (input$DrugList %in% drugs_category$yellow) {
      mycategory <- "yellow"
    } else {
      mycategory <- "green"
    }
    paste("This drug is marked as", mycategory)
  })
  
  output$GenoDetails <- DT::renderDataTable({
    selected_sample <- input$SampleList
    mygenos <- as.data.frame(t(samples_genos[samples_genos$Sample==selected_sample,c(3,5:(ncol(samples_genos)-1))]))
    colnames(mygenos) <- "Genotype"
    mygenos
  }, options=list(pageLength = 25))
}

# Run the application 
shinyApp(ui = ui, server = server)

