library(shiny)
library(r2d3)
#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
print(load("tbl.all-11492x14.RData"))  # tbl.all
max.time.points <- 9
fraction.names <- sort(unique(tbl.all$fraction))
goi <- c("ALL proteins", sort(unique(tbl.all$gene)))

print(load("tbl.complexes.RData"))
complexes <- c("ALL", sort(unique(tbl.complexes$complex)))
#------------------------------------------------------------------------------------------------------------------------
srm.coexpression.tab <- function()
{
   sidebarLayout(
      sidebarPanel(
         verbatimTextOutput(outputId="currentCurveCountDisplay"),
         checkboxGroupInput("fractionSelector",
                            label="category",
                            choices = fraction.names,
                            selected = fraction.names[1:3]),
         radioButtons("srm.transformChoice", "Data Transform", c("None", "Normalized")), # , "Arcsinh")),
         radioButtons("srm.lineTypeSelector", "Smoothing", c("No", "Yes")),
         selectizeInput("proteinSelector", "Protein/s:", goi, selected=goi[2],  multiple=TRUE,
                        options=list(maxOptions=nrow(tbl.all))),
         selectizeInput("complexSelector", "Draw from Complex:", complexes, selected=complexes[2],
                        multiple=TRUE, options=list(maxOptions=length(complexes))),
         sliderInput("correlationThresholdSlider", label = "Pearson", min = 0, max = 1, value = 0.9, step = 0.01),
         radioButtons("correlationDirectionChooser", "Find Correlations", c("None", "+", "-")),
         br(),
         verbatimTextOutput(outputId="currentVectorDisplay"),
         width=2
         ),
      mainPanel(
         d3Output("srm.d3", height="80vh"),
         width=10
         )
      )

} # srm.coexpression.tab
#------------------------------------------------------------------------------------------------------------------------
ui <- fluidPage(

   tags$head(
        tags$style(".tab-pane {margin-top: 20px;}"),
        tags$link(rel = "stylesheet", type = "text/css", href = "app.css")
        ),

   # titlePanel("Transcription Factor Protein and RNA Expression Profiles During Erythropoiesis"),

   tabsetPanel(
       tabPanel("Temporal Protein Abundance", srm.coexpression.tab())
       ) # tabsetPanel

   ) # fluidPage

#------------------------------------------------------------------------------------------------------------------------
server <- function(input, output, session) {

   reactiveState <- reactiveValues(selectedTF=NULL, correlatedTFs=list())

    currentCount <- reactive({
      printf("--- entering currentCount()")
      proteins <- input$proteinSelector
      if("ALL proteins" %in% proteins)
          proteins <- unique(tbl.all$gene)
      printf("proteins: %s", paste(head(proteins), collapse=", "))
      fractions <- input$fractionSelector
      complexes <- input$complexSelector
      tbl.complexes.sub <- subset(tbl.complexes, complex %in% complexes)
      proteins <- intersect(proteins, tbl.complexes.sub$gene)
      printf("fractions: %s", paste(fractions, collapse=", "))
      printf("complexes: %s", paste(complexes, collapse=", "))
      tbl.sub <- subset(tbl.all, gene %in% proteins & fraction %in% fractions)
      result <- nrow(tbl.sub)
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
      plotCorrelatedProteins(input, output)
      })

   observeEvent(input$complexSelector, ignoreInit=FALSE, {
       printf("--- complex selected: %s", input$complexSelector)
       })

   observeEvent(input$geneSelector, ignoreInit=TRUE, {
      tf <- input$geneSelector
      })

   observeEvent(input$srm.transformChoice, ignoreInit=TRUE, {
     plotCorrelatedProteins(input, output)
     })

   observeEvent(input$correlationThresholdSlider, ignoreInit=TRUE, {
      plotCorrelatedProteins(input, output)
      })

   observeEvent(input$correlationDirectionChooser, ignoreInit=TRUE, {
      plotCorrelatedProteins(input, output)
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
#------------------------------------------------------------------------------------------------------------------------
transformData.rna.srm <- function(rna, srm, transformName)
{
   printf("--- transform by %s", transformName)

   if(transformName == "None"){
      rna.out <- rna;
      srm.out <- srm;
      }

   if(transformName == "Normalized"){
      rna.out <- rna/max(rna)
      srm.out <- srm/max(srm);
      }

   if(transformName == "Arcsinh"){
      rna.out <- asinh(rna)
      srm.out <- asinh(srm)
      }

   return(list(rna=rna.out, srm=srm.out))

} # transformData.rna.srm
#------------------------------------------------------------------------------------------------------------------------
findCorrelated <- function(targetTF, threshold, direction)
{
   if(direction == "None")
      return(targetTF)

   suppressWarnings({
      mtx <- as.matrix(tbl.all[, grep("^D", colnames(tbl.all))])
      rownames(mtx) <- tbl.all$gene
      correlations <- apply(mtx, 1,
                           function(row) cor(mtx[targetTF,], row,  use="complete.obs"))
      })

   if(direction == "-")
      result <- names(which(correlations <= (-1 * threshold)))
   else  # must be "+"
      result <- names(which(correlations >= threshold))

   return(unique(c(targetTF, result)))

} # findCorrelated
#------------------------------------------------------------------------------------------------------------------------
plotCorrelatedProteins <- function(input, output)
{
   protein <- input$proteinSelector[1]
   correlationThreshold <- input$correlationThresholdSlider;
   correlationDirection <- isolate(input$correlationDirectionChooser)
   proteins.all <- findCorrelated(protein, correlationThreshold, correlationDirection)
   transform <- input$srm.transformChoice
   proteinCategories <- input$fractionSelector

   output$srm.d3 <- renderD3({
     plotProteins(proteins.all, input, output, transform, proteinCategories)
     })

} # plotCorrelatedProteins
#------------------------------------------------------------------------------------------------------------------------
plotProteins <- function(proteins, input, output, transform, proteinCategories)
{
   printf("plotProteins (%s): %s", transform, paste(proteins, collapse=", "))
   printf("  categories to include: %s", paste(proteinCategories, collapse=", "))

   tbl.sub <- subset(tbl.all, gene %in% proteins & fraction %in% proteinCategories)
   protein.fraction.names <- sprintf("%s-%s", tbl.sub$gene, tbl.sub$fraction) # eg, "A2M-cyto" & "A2M-ne1"
   mtx <- tbl.sub[, grep("^D", colnames(tbl.all), ignore.case=TRUE)] # just the time columns

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

} # plotProteins
#------------------------------------------------------------------------------------------------------------------------
maxOfVectors <- function(vectorList)
{
   max <- 0
   for(vector in vectorList){
      vector.max <- max(vector, na.rm=TRUE)
      #if(is.na(vector.max)) browser()
      if(vector.max > max)
         max <- vector.max
      } # for vector

   return(max)

} # maxOfVectors
#------------------------------------------------------------------------------------------------------------------------
transformData.srm <- function(srm, transformName)
{
   # printf("--- transform by %s", transformName)

   if(transformName == "None"){
      srm.out <- srm;
      }

   if(transformName == "Normalized"){
      srm.out <- lapply(srm, function(vec) vec/max(vec, na.rm=TRUE))
      }

   if(transformName == "Arcsinh"){
      srm.out <- lapply(srm, asinh)
      }

   return(srm.out)

} # transformData.srm
#------------------------------------------------------------------------------------------------------------------------
app <- shinyApp(ui = ui, server = server)

