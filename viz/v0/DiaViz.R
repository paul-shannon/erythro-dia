library(shiny)
library(R6)
library(r2d3)
#----------------------------------------------------------------------------------------------------
buttonStyle <- "margin: 5px; font-size: 20px;"
textOutputStyle <- paste0("margin:10px; margin-left: 50px;",
		          " padding:5px; width: 200px; height: 60px; color:red; ",
		          "border: 1px solid black; font-size: 20px;")

#----------------------------------------------------------------------------------------------------
print(load("tbl.all-11492x14.RData"))  # tbl.all
max.time.points <- 9
fraction.names <- sort(unique(tbl.all$fraction))
#goi <- sort(unique(tbl.all$gene))

print(load("tbl.complexes.RData"))
complexes <- sort(unique(tbl.complexes$complex))
#----------------------------------------------------------------------------------------------------
DiaVizApp = R6Class("DiaVizApp",

    #--------------------------------------------------------------------------------
    private = list(proteinCount=NULL,
                   tbl.all=data.frame(),
                   tbl.current=data.frame()
                   ),

    #--------------------------------------------------------------------------------
    public = list(

        initialize = function(){
            print(noquote(sprintf("initializing DiaViz")))
            private$tbl.all <- tail(get(load("tbl.all-11492x14.RData")), n=30)
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
            protein <- input$proteinSelector[1]
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
        plotProteins= function(proteins, input, output, transform, proteinCategories){
            printf("plotProteins (%s): %s", transform, paste(proteins, collapse=", "))
            printf("  categories to include: %s", paste(proteinCategories, collapse=", "))

            tbl.sub <- subset(private$tbl.current, gene %in% proteins & fraction %in% proteinCategories)
            protein.fraction.names <- sprintf("%s-%s", tbl.sub$gene, tbl.sub$fraction) # eg, "A2M-cyto" & "A2M-ne1"
            mtx <- tbl.sub[, grep("^D", colnames(private$tbl.current), ignore.case=TRUE)] # just the time columns

            if(nrow(mtx) > 0 & transform == "Normalized"){
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
            r2d3(data, script = "multiPlot.js")

            },

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
               sidebarLayout(
                   sidebarPanel(
                       verbatimTextOutput(outputId="currentCurveCountDisplay"),
                       verbatimTextOutput(outputId="currentVectorDisplay"),
                       checkboxGroupInput("fractionSelector",
                                          label="category",
                                          choices = fraction.names,
                                          selected = fraction.names[1:3]),
                       selectizeInput("proteinSelector",
                                      label="Protein/s:",
                                      choices=NULL, #sort(unique(private$tbl.current$gene)),
                                      selected=NULL,
                                      multiple=TRUE,
                                      options=list(maxOptions=nrow(private$tbl.current))),
                       selectizeInput("complexSelector", "Draw from Complex:", complexes, selected=NULL,
                                      multiple=TRUE, options=list(maxOptions=length(complexes))),
                       actionButton("plotCurrentSelelctionButton", "Plot Current Selection"),
                       br(),
                       br(),
                       sliderInput("correlationThresholdSlider", label = "Pearson", min = 0, max = 1, value = 0.9, step = 0.01),
                       radioButtons("correlationDirectionChooser", "Find Correlations", c("None", "+", "-")),
                       br(),
                       radioButtons("srm.transformChoice", "Data Transform", c("None", "Normalized")), # , "Arcsinh")),
                       radioButtons("srm.lineTypeSelector", "Smoothing", c("No", "Yes")),
                       width=2
                       ),
                   mainPanel(
                       d3Output("srm.d3", height="80vh"),
                       width=10
                       )
                   ) # sidebarLayout
           )},

        #------------------------------------------------------------
        server = function(input, output, session){

            print(noquote(sprintf("entering server")))

            currentCount <- reactive({
                printf("--- entering currentCount()")
                tbl.tmp <- private$tbl.all
                printf("==== nrow(tbl.tmp) 1: %d", nrow(tbl.tmp))

                complexes <- input$complexSelector
                printf("complexes: %s", paste(complexes, collapse=", "))
                if(!is.null(complexes)){
                    printf("filtering on chosen complexes")
                    tbl.complexes.sub <- subset(tbl.complexes, complex %in% complexes)
                    printf(" nrow(tbl.complexes.sub): %d", nrow(tbl.complexes.sub))
                    tbl.tmp <- subset(tbl.tmp, gene %in% tbl.complexes.sub$gene)
                    printf("tbl.tmp, filtered for complexes, now has %d rows", nrow(tbl.tmp))
                    }
                printf("==== nrow(tbl.tmp) 2: %d", nrow(tbl.tmp))
                fractions <- input$fractionSelector
                printf("fractions: %s", paste(fractions, collapse=", "))
                if(length(fractions) == 0)
                    tbl.tmp <- subset(tbl.tmp, fraction == "none specified")
                printf("==== nrow(tbl.tmp) 3: %d", nrow(tbl.tmp))
                if(length(fractions) > 0)
                    tbl.tmp <- subset(tbl.tmp, fraction %in% fractions)

                printf("==== nrow(tbl.tmp) 4: %d", nrow(tbl.tmp))
                proteins <- input$proteinSelector
                printf("--- current proteins: ")
                print(proteins)
                if(!is.null(proteins))
                   tbl.tmp <- subset(tbl.tmp, gene %in% proteins)
                printf("==== nrow(tbl.tmp) 5: %d", nrow(tbl.tmp))
                #if(is.null(proteins)){
                #    surviving.proteins <- sort(unique(tbl.tmp$gene))
                #    printf("head(surviving.proteins): %s", paste(head(surviving.proteins), collapse=", "))
                #    printf("--- updating protein selectizor with %d proteins", length(surviving.proteins))
                #    updateSelectizeInput(session, input$proteinSelector,
                #                         choices=surviving.proteins, server=TRUE)
                #    }
                printf("==== nrow(tbl.tmp) 6: %d", nrow(tbl.tmp))
                private$tbl.current <- tbl.tmp
                result <- nrow(tbl.tmp)
                printf("  currentCount: %d", result)
                return(result)
                })

            observeEvent(input$fractionSelector, ignoreInit=FALSE, {
                newChoices <- input$fractionSelector
                printf("new fractions: %s", paste(newChoices, collapse=", "))
                output$currentCurveCountDisplay <- renderText({currentCount()})
                                        #plotCorrelatedProteins(input, output)
                })

            observeEvent(input$proteinSelector, ignoreInit=FALSE, {
                self$plotCorrelatedProteins(input, output)
                })

            observeEvent(input$complexSelector, ignoreInit=FALSE, {
                printf("--- complex selected: %s", input$complexSelector)
                })

            observeEvent(input$geneSelector, ignoreInit=TRUE, {
                tf <- input$geneSelector
                })

            observeEvent(input$srm.transformChoice, ignoreInit=TRUE, {
                self$plotCorrelatedProteins(input, output)
                })

            observeEvent(input$correlationThresholdSlider, ignoreInit=TRUE, {
                self$plotCorrelatedProteins(input, output)
                })

            observeEvent(input$correlationDirectionChooser, ignoreInit=TRUE, {
                self$plotCorrelatedProteins(input, output)
                })


            observeEvent(input$currentlySelectedVector, ignoreInit=FALSE, {
                newValue <- input$currentlySelectedVector
                                        # printf("newValue: %s", newValue)
                if(nchar(newValue) == 0)
                    newValue <- "   "
                output$currentVectorDisplay <- renderText({newValue})
                                        #output$currentVectorDisplay <- renderText({newValue})
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
#shinyApp(ui = ui, server = server)

app <- DiaVizApp$new()
x <- shinyApp(app$ui, app$server)
runApp(x, port=1156)

