library(shiny)
library(R6)
#----------------------------------------------------------------------------------------------------
buttonStyle <- "margin: 5px; font-size: 20px;"
textOutputStyle <- paste0("margin:10px; margin-left: 50px;",
		          " padding:5px; width: 200px; height: 60px; color:red; ",
		          "border: 1px solid black; font-size: 20px;")

DiaVizApp = R6Class("DiaVizApp",

    #--------------------------------------------------------------------------------
    private = list(latestText=NULL
                   ),

    #--------------------------------------------------------------------------------
    public = list(

        initialize = function(){
            print(noquote(sprintf("initializing DiaViz")))
            },

        #------------------------------------------------------------
        ui = function(){
           fluidPage(
              wellPanel(style="width: 1000px;",
                actionButton("randomTextButton", label= "Generate random Text",
                             style=buttonStyle)),
                div(style=textOutputStyle, textOutput("textDisplay"))

            )},

        #------------------------------------------------------------
        server = function(input, output, session){

            print(noquote(sprintf("entering server")))

            observeEvent(input$randomTextButton, ignoreInit=TRUE, {
              randomText <- paste(sample(c(LETTERS, letters), 10, replace=TRUE), collapse="")
              private$latestText <- randomText;  # though not currently reused
	      output$textDisplay <- renderText(randomText);
              })
            } # server

       ) # public
    ) # class
#--------------------------------------------------------------------------------
deploy <-function()
{
   repos <- options("repos")[[1]]
   stopifnot(sort(names(repos)) == c("BioCann", "BioCsoft", "CRAN"))
   stopifnot(repos$BioCann=="https://bioconductor.org/packages/3.13/data/annotation")
   stopifnot(repos$BioCsoft=="https://bioconductor.org/packages/3.13/bioc")
   stopifnot(repos$CRAN=="https://cran.microsoft.com")
   require(devtools)

      # jim hester suggests, with reference
      # Setting R_REMOTES_NO_ERRORS_FROM_WARNINGS="false" will cause warning
      # messages during calls to install.packages() to become errors. Often warning
      # messages are caused by dependencies failing to install.
   Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")

   require(rsconnect)

   deployApp(account="paulshannon",
              appName="diaViz",
              appTitle="DiaViz",
              appFiles=c("DiaViz.R"),
              appPrimaryDoc="DiaViz.R"
              )

} # deploy
#------------------------------------------------------------------------------------------------------------------------
shinyApp(ui = ui, server = server)

app <- DiaVizApp$new()
x <- shinyApp(app$ui, app$server)
runApp(x, port=1156)

