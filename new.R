# Shiny app for miRNA qPCR analysis

library(shiny)
library(readxl)
library(dplyr)
library(ggplot2)

ui <- fluidPage(
  titlePanel("miRNA qPCR Analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload Excel File", accept = c(".xlsx")),
      helpText("Excel file must have columns: Well, Fluor, Target, Sample, Cq, Cq Mean, Type, Sd"),
      uiOutput("targetSelector")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Dot Plots", uiOutput("dotPlotsUI")),
        tabPanel("Bar Plot", plotOutput("barPlot")),
        tabPanel("Table", tableOutput("summaryTable"))
      )
    )
  )
)

server <- function(input, output) {

  # Reactive expression to load and process raw data, removing NA
  rawData <- reactive({
    req(input$file)
    df <- read_excel(input$file$datapath)
    df <- df %>% filter(!is.na(`Cq Mean`) & !is.na(Target) & !is.na(Type))
    return(df)
  })

  # Summary means for table and bar plot
  dataInput <- reactive({
    df <- rawData()
    df_means <- df %>%
      group_by(Target, Type) %>%
      summarise(mean_Cq = mean(`Cq Mean`, na.rm = TRUE),
                sd_Cq = sd(`Cq Mean`, na.rm = TRUE), .groups = "drop")
    return(df_means)
  })

  # Dynamic UI for target selection
  output$targetSelector <- renderUI({
    df <- rawData()
    targets <- unique(df$Target)
    checkboxGroupInput("selectedTargets", "Select miRNAs to display:",
                       choices = targets, selected = targets)
  })

  # Generate dot plots per target showing individual samples
  output$dotPlotsUI <- renderUI({
    df <- rawData()
    targets <- input$selectedTargets
    if (is.null(targets)) return(NULL)

    plot_output_list <- lapply(targets, function(tg) {
      plotname <- paste0("plot_", tg)
      plotOutput(plotname, height = 300)
    })

    do.call(tagList, plot_output_list)
  })

  observe({
    df <- rawData()
    targets <- input$selectedTargets
    if (is.null(targets)) return(NULL)

    for (tg in targets) {
      local({
        tg_local <- tg
        plotname <- paste0("plot_", tg_local)
        output[[plotname]] <- renderPlot({
          df_subset <- df %>% filter(Target == tg_local & !is.na(`Cq Mean`) & !is.na(Type))
          ggplot(df_subset, aes(x = Type, y = `Cq Mean`, color = Type)) +
            geom_jitter(width = 0.2, size = 3, alpha = 0.7) +  # individual sample dots
            stat_summary(fun = mean, geom = "point", shape = 18, size = 5, color = "black",
                         position = position_dodge(width = 0.5)) +
            scale_color_manual(values = c("Normal" = "grey", "Tumor" = "darkred")) +
            labs(title = paste("Cq values for", tg_local),
                 y = "Cq",
                 x = "") +
            theme_minimal()
        })
      })
    }
  })

  # Bar Plot of overall means
  output$barPlot <- renderPlot({
    df_means <- dataInput()
    df_means <- df_means %>% filter(!is.na(mean_Cq) & !is.na(Type))
    overall_means <- df_means %>%
      group_by(Type) %>%
      summarise(overall_Cq = mean(mean_Cq, na.rm = TRUE))

    ggplot(overall_means, aes(x = Type, y = overall_Cq, fill = Type)) +
      geom_bar(stat = "identity", width = 0.6) +
      labs(title = "Overall Average Cq (All miRNAs)",
           y = "Mean Cq",
           x = "") +
      theme_minimal()
  })

  # Summary Table
  output$summaryTable <- renderTable({
    dataInput()
  })
}

shinyApp(ui = ui, server = server)
