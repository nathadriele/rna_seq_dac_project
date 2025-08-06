library(shiny)
library(ggplot2)
library(DT)
library(readr)

res <- read_csv("../results/deseq2_results.csv")
go <- read_csv("../results/go_enrichment.csv")

ui <- fluidPage(
  titlePanel("RNA-Seq Dashboard - DAC"),
  sidebarLayout(
    sidebarPanel(
      selectInput("selected_gene", "Select a gene to inspect:",
                  choices = res$gene[1:100], selected = res$gene[1])
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Volcano Plot", plotOutput("volcanoPlot")),
        tabPanel("Top Genes Table", DTOutput("geneTable")),
        tabPanel("GO Enrichment", DTOutput("goTable"))
      )
    )
  )
)

server <- function(input, output) {
  output$volcanoPlot <- renderPlot({
    res$significant <- res$padj < 0.05
    ggplot(res, aes(log2FoldChange, -log10(padj), color = significant)) +
      geom_point(alpha = 0.6) +
      geom_text(data = subset(res, gene == input$selected_gene),
                aes(label = gene), vjust = -1, hjust = 1.1, color = "black") +
      theme_minimal()
  })

  output$geneTable <- renderDT({
    datatable(res[, c("gene", "log2FoldChange", "pvalue", "padj")], options = list(pageLength = 10))
  })

  output$goTable <- renderDT({
    datatable(go[, c("Description", "GeneRatio", "p.adjust", "Count")], options = list(pageLength = 10))
  })
}

shinyApp(ui = ui, server = server)
