#library(R6)
#library(shiny)
#library(shinyjs)
library(ProteomicsFilterWidget)
#----------------------------------------------------------------------------------------------------
f <- "tbl.all-dia-rna.38662-14.RData"
stopifnot(file.exists(f))
tbl.all <- get(load(f))
# f <- "tbl.complexes.RData"
f <- "tbl.complexes-12jan2023.RData"
stopifnot(file.exists(f))
tbl.complexes <- get(load(f))
#----------------------------------------------------------------------------------------------------
buttonStyle <- "margin: 5px; margin-right: 0px; font-size: 14px;"

DiaVizApp = R6Class("DiaVizApp",

    #--------------------------------------------------------------------------------
    private = list(proteomicsFilter1 = NULL,
                   proteomicsFilter2 = NULL
                   ),
    #--------------------------------------------------------------------------------
    public = list(

        initialize = function(){
            message(sprintf("initializing demo"))
            private$proteomicsFilter1 <- ProteomicsFilteringWidget$new(id="d3.dia", tbl.all, tbl.complexes)
            },

        #------------------------------------------------------------
        ui = function(){

           fluidPage(
              tabsetPanel(
                     tabPanel("Filtering", private$proteomicsFilter1$ui())
                     ) # tabsetPanel
              ) # fluidPage
            },

        #------------------------------------------------------------
        server = function(input, output, session){
            message(sprintf("entering ProteomicsFilterWidgetDemo::server"))
            private$proteomicsFilter1$server(input, output, session)
            message(sprintf("leaving ProteomicsFilterWidgetDemo::server"))
            }
       ) # public
    ) # class
#--------------------------------------------------------------------------------
deploy <-function()
{
   repos <- options("repos")[[1]]
   stopifnot(sort(names(repos)) == c("BioCann", "BioCsoft", "CRAN"))
   stopifnot(repos$BioCann=="https://bioconductor.org/packages/3.16/data/annotation")
   stopifnot(repos$BioCsoft=="https://bioconductor.org/packages/3.16/bioc")
   stopifnot(repos$CRAN=="https://cran.microsoft.com")

   require(devtools)

   install_github("PriceLab/BioShiny/ProteomicsFilteringWidget", force=TRUE)
   install_github("daattali/shinyjs", force=TRUE)

      # jim hester suggests, with reference
      # Setting R_REMOTES_NO_ERRORS_FROM_WARNINGS="false" will cause warning
      # messages during calls to install.packages() to become errors. Often warning
      # messages are caused by dependencies failing to install.
   Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")

   require(rsconnect)

   deployApp(account="hoodlab",
             appName="diaViz2",
             appTitle="DiaViz2",
             appFiles=c("DiaViz.R",
                        "tbl.all-dia-rna.38662-14.RData",
                        "tbl.complexes-12jan2023.RData"),
              appPrimaryDoc="DiaViz.R"
              )

} # deploy
#----------------------------------------------------------------------------------------------------
app <- DiaVizApp$new()
# shinyApp(app$ui(), app$server)
#if(grepl("hagfish", Sys.info()[["nodename"]])){ #  & interactive()){
#   printf("--- on hagfish")
#   runApp(shinyApp(app$ui(), app$server), port=1112)
#   } else {
   shinyApp(app$ui(), app$server)
#   }



