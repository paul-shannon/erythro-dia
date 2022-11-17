library(shiny)
library(R6)
library(r2d3)
#----------------------------------------------------------------------------------------------------
buttonStyle <- "margin: 5px; font-size: 20px;"
textOutputStyle <- paste0("margin:10px; margin-left: 50px;",
		          " padding:5px; width: 200px; height: 60px; color:red; ",
		          "border: 1px solid black; font-size: 20px;")

#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
print(load("tbl.all-11492x14.RData"))  # tbl.all
max.time.points <- 9
fraction.names <- sort(unique(tbl.all$fraction))

tbl.complexes <- get(load("tbl.complexes.RData"))
complexes <- sort(unique(tbl.complexes$complex))
#----------------------------------------------------------------------------------------------------
DiaVizApp = R6Class("DiaVizApp",

    #--------------------------------------------------------------------------------
    private = list(proteinCount=NULL,
                   tbl.all=data.frame(),
                   tbl.current=data.frame(),
                   tbl.selected=data.frame(),
                   currentProteins=NULL,
                   transform="None"
                   ),

    #--------------------------------------------------------------------------------
    public = list(

        initialize = function(){
            print(noquote(sprintf("initializing DiaViz")))
            private$tbl.all <- tail(get(load("tbl.all-11492x14.RData")), n=-1)
            private$tbl.current <- private$tbl.all
            },

        findCorrelated = function(targetTF, threshold, direction){
            if(direction == "None")
                return(targetTF)

            suppressWarnings({
                mtx <- as.matrix(private$tbl.current[, grep("^D", colnames(private$tbl.current))])
                rownames(mtx) <- private$tbl.current$gene
                correlations <- apply(mtx, 1,
                                      function(row) cor(mtx[targetTF,], row,  use="complete.obs"))
                })

            if(direction == "-")
                result <- names(which(correlations <= (-1 * threshold)))
            else  # must be "+"
                result <- names(which(correlations >= threshold))
            return(unique(c(targetTF, result)))
            }, # findCorrelated

        #------------------------------------------------------------
        plotCorrelatedProteins=function(input, output){
            proteins <- input$proteinSelector[1]
            correlationThreshold <- input$correlationThresholdSlider;
            correlationDirection <- isolate(input$correlationDirectionChooser)
            proteins.all <- self$findCorrelated(protein, correlationThreshold, correlationDirection)
            transform <- input$srm.transformChoice
            proteinCategories <- input$fractionSelector

            output$srm.d3 <- renderD3({
                self$plotProteins(proteins.all, input, output, transform, proteinCategories)
                })
            }, # plotCorrelatedProteins

        #------------------------------------------------------------
        plotProteins=function(input, output){
            printf("plotProteins, nrow(tbl.selected): %d", nrow(private$tbl.selected))
            printf("plotProteins, currentProteins:  %d", length(private$currentProteins))
            print("--- tbl.selected")
            #print(head(private$tbl.selected))
            tbl.sub <- subset(private$tbl.selected, gene %in% private$currentProteins)
            printf("plotProteins, nrow(tbl.sub): %d", nrow(tbl.sub))
            protein.fraction.names <- sprintf("%s-%s", tbl.sub$gene, tbl.sub$fraction) # eg, "A2M-cyto" & "A2M-ne1"
            mtx <- tbl.sub[, grep("^D", colnames(private$tbl.current), ignore.case=TRUE)] # just the time columns

            if(nrow(mtx) > 0 & private$transform == "Normalized"){
                mtx <- t(apply(mtx, 1, function(row) row/max(row)))
                }

            rownames(mtx) <- protein.fraction.names
            timePoints <- as.numeric(sub("D", "", colnames(mtx)))
            protein.count.vectors <- lapply(seq_len(nrow(mtx)), function(row) as.numeric(mtx[row,]))
            names(protein.count.vectors) <- protein.fraction.names

            xMin <- min(timePoints)
            xMax <- max(timePoints)
            yMin <- 0
            yMax <- max(mtx)

            vectorsWithTimes <- vector(mode="list", length(rownames(mtx)))
            names(vectorsWithTimes) <- rownames(mtx)

            for(protein.fraction.name in protein.fraction.names){
                vector <- protein.count.vectors[[protein.fraction.name]]
                vectorsWithTimes[[protein.fraction.name]] <- lapply(seq_len(length(timePoints)),
                                                                    function(i)
                                                                        return(list(x=timePoints[i], y=vector[i])))
                } # for protein.fraction.name

            lineSmoothing <- input$srm.lineTypeSelector
            data <- list(vectors=vectorsWithTimes, xMax=xMax, yMax=yMax, cmd="plot", smoothing=lineSmoothing)
            output$srm.d3 <- renderD3({
                r2d3(data, script = "multiPlot.js")
                })
            }, # plotProteins

        #------------------------------------------------------------
        maxOfVectors = function(vectorList){
            max <- 0
            for(vector in vectorList){
                vector.max <- max(vector, na.rm=TRUE)
                                        #if(is.na(vector.max)) browser()
                if(vector.max > max)
                    max <- vector.max
                } # for vector
            return(max)
            }, # maxOfVectors

        #------------------------------------------------------------
        ui = function(){
           fluidPage(
               br(), br(),
               sidebarLayout(
                   sidebarPanel(
                       checkboxGroupInput("fractionSelector",
                                          label="Choose Cellular Fractions:",
                                          choices = fraction.names,
                                          selected = fraction.names[1:3]),
                       selectizeInput("complexSelector", "Choose Complexes:", complexes, selected=NULL,
                                      multiple=TRUE, options=list(maxOptions=length(complexes))),
                       h6("Filtered Set Size:"),
                       verbatimTextOutput(outputId="currentCurveCountDisplay"),
                       hr(),
                       h6("Choose Proteins From Filtered Set:"),
                       selectizeInput("proteinSelector",
                                      label=NULL,
                                      choices=NULL, #sort(unique(private$tbl.current$gene)),
                                      selected=NULL,
                                      multiple=TRUE,
                                      options=list(maxOptions=nrow(private$tbl.current))),
                       verbatimTextOutput(outputId="currentSubsetCountDisplay"),
                       actionButton("plotCurrentSelectionButton", "Plot Current Selection"),
                       br(),
                       #br(),
                       #sliderInput("correlationThresholdSlider", label = "Pearson", min = 0, max = 1, value = 0.9, step = 0.01),
                       #radioButtons("correlationDirectionChooser", "Find Correlations", c("None", "+", "-")),
                       br(),
                       radioButtons("srm.transformChoice", "Data Transform", c("None", "Normalized")), # , "Arcsinh")),
                       #radioButtons("srm.lineTypeSelector", "Smoothing", c("No", "Yes")),
                       verbatimTextOutput(outputId="currentVectorDisplay"),
                       width=3
                       ),
                   mainPanel(
                       d3Output("srm.d3", height="80vh"),
                       width=9
                       )
                   ) # sidebarLayout
           )},

        #------------------------------------------------------------
        server = function(input, output, session){

            print(noquote(sprintf("entering server")))


            observeEvent(input$proteinSelector, ignoreInit=FALSE, {
                proteins <- input$proteinSelector
                private$currentProteins <- proteins
                private$tbl.selected <- subset(private$tbl.current, gene %in% proteins)
                row.count <- nrow(private$tbl.selected)
                printf("tbl.selected has %d rows", row.count)
                text <- sprintf("%d rows, %d proteins", row.count, length(proteins))
                output$currentSubsetCountDisplay <- renderText(text)
                })


            currentTable <- reactive({
                printf("------------------------------------------- entering currentTable()")
                tbl.tmp <- private$tbl.all
                complexes <- input$complexSelector
                fractions <- input$fractionSelector
                #proteins <- isolate(input$proteinSelector)
                #proteins <- NULL
                #    proteins <- input$proteinSelector
                #    }

                printf("==== nrow(tbl.tmp) 1: %d", nrow(tbl.tmp))

                printf("complexes: %s", paste(complexes, collapse=", "))
                if(!is.null(complexes)){
                    printf("filtering on chosen complexes")
                    tbl.complexes.sub <- subset(tbl.complexes, complex %in% complexes)
                    printf(" nrow(tbl.complexes.sub): %d", nrow(tbl.complexes.sub))
                    tbl.tmp <- subset(tbl.tmp, gene %in% tbl.complexes.sub$gene)
                    printf("tbl.tmp, filtered for complexes, now has %d rows", nrow(tbl.tmp))
                    }
                printf("==== nrow(tbl.tmp) 2: %d", nrow(tbl.tmp))
                printf("fractions: %s", paste(fractions, collapse=", "))
                if(length(fractions) == 0)
                    tbl.tmp <- subset(tbl.tmp, fraction == "none specified")
                printf("==== nrow(tbl.tmp) 3: %d", nrow(tbl.tmp))
                if(length(fractions) > 0)
                    tbl.tmp <- subset(tbl.tmp, fraction %in% fractions)

                printf("==== nrow(tbl.tmp) 4: %d", nrow(tbl.tmp))

                printf("==== nrow(tbl.tmp) 5: %d", nrow(tbl.tmp))
                printf("==== nrow(tbl.tmp) 6: %d", nrow(tbl.tmp))
                private$tbl.current <- tbl.tmp
                return(nrow(private$tbl.current))
                }) # currentTable

            observe({
                row.count <- currentTable()
                unique.proteins <- sort(unique(private$tbl.current$gene))
                printf("    unique.proteins: %d", length(unique.proteins))
                updateSelectizeInput(session = session,
                                     inputId = "proteinSelector",
                                     choices = unique.proteins,
                                     server=TRUE
                                     )
                printf("--- observe, new row count: %d", row.count)
                output$currentCurveCountDisplay <-
                    renderText(sprintf("%d rows, %d proteins", row.count, length(unique.proteins)))
                })

            observeEvent(input$plotCurrentSelectionButton, ignoreInit=FALSE, {
                printf("--- plotCurrentSelectionButton")
                printf("plot %d proteins", length(private$currentProteins))
                self$plotProteins(input, output)
                #self$plotCorrelatedProteins(input, output)
                })

            # output$currentCurveCountDisplay <- renderText(currentTable())

            #observeEvent(input$fractionSelector, ignoreInit=FALSE, {
            #    newChoices <- input$fractionSelector
            #    printf("new fractions: %s", paste(newChoices, collapse=", "))
                #output$currentCurveCountDisplay <- renderText({nrow(private$tbl.current)})
                                        #plotCorrelatedProteins(input, output)
             #   })

            #observeEvent(input$proteinSelector, ignoreInit=FALSE, {
            #    #self$plotCorrelatedProteins(input, output)
            #    })

            #observeEvent(input$complexSelector, ignoreInit=FALSE, {
            #    printf("--- complex selected: %s", input$complexSelector)
            #    })

            #observeEvent(input$geneSelector, ignoreInit=TRUE, {
            #    tf <- input$geneSelector
            #    })

            observeEvent(input$srm.transformChoice, ignoreInit=TRUE, {
                new.choice <- input$srm.transformChoice
                printf("--- setting private$transform: %s", new.choice)
                private$transform <- new.choice
                })

            observeEvent(input$correlationThresholdSlider, ignoreInit=TRUE, {
                #self$plotCorrelatedProteins(input, output)
                })

            observeEvent(input$correlationDirectionChooser, ignoreInit=TRUE, {
                #self$plotCorrelatedProteins(input, output)
                })


            observeEvent(input$currentlySelectedVector, ignoreInit=FALSE, {
                newValue <- input$currentlySelectedVector
                printf("newValue: %s", newValue)
                if(nchar(newValue) == 0) newValue <- "   "
                output$currentVectorDisplay <- renderText({newValue})
                })

            } # server

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

      # jim hester suggests, with reference
      # Setting R_REMOTES_NO_ERRORS_FROM_WARNINGS="false" will cause warning
      # messages during calls to install.packages() to become errors. Often warning
      # messages are caused by dependencies failing to install.
   Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")

   require(rsconnect)

   deployApp(account="hoodlab",
             appName="diaViz",
             appTitle="DiaViz",
             appFiles=c("DiaViz.R",
                        "multiPlot.js", "linePlot.js", "linePlot.css", "multiPlot.css",
                        "tbl.all-11492x14.RData",
                        "tbl.complexes.RData"),
              appPrimaryDoc="DiaViz.R"
              )

} # deploy
#------------------------------------------------------------------------------------------------------------------------
#shinyApp(ui = ui, server = server)

app <- DiaVizApp$new()
shinyApp(app$ui, app$server)
#runApp(x, port=1156)

