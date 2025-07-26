required_packages <- c(
  "shiny", "readxl", "dplyr", "tidyr", "ggplot2",
  "writexl", "shinyWidgets", "DT", "rmarkdown",
  "knitr", "openxlsx"
)

installed_packages <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Load them all after installing
lapply(required_packages, library, character.only = TRUE)

ui <- fluidPage(
  titlePanel("ΔΔCt Analysis Tool"),

  sidebarLayout(
    sidebarPanel(
      fileInput("normal_file", "Upload Normal Sample Excel", accept = ".xlsx"),
      fileInput("tumor_file", "Upload Tumor Sample Excel", accept = ".xlsx"),
      pickerInput("ref_gene", "Select Reference Gene(s)", choices = NULL, multiple = TRUE),
      pickerInput("internal_controls", "Select Internal Controls", choices = NULL, multiple = TRUE),
      actionButton("analyze", "Analyze"),
      downloadButton("download_results", "Download Results"),
      downloadButton("download_plot_excel", "Download Internal Controls Plot as Excel"),
      downloadButton("download_mirna_plots_excel", "Download miRNA Plots as Excel")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("ΔΔCt Tables",
                 uiOutput("normal_tables"),
                 uiOutput("tumor_tables")
        ),
        tabPanel("Internal Control Plot", plotOutput("internal_control_plot")),
        tabPanel("miRNA Ct Plots", uiOutput("mirna_plots"))
      )
    )
  )
)

server <- function(input, output, session) {
  dataInput <- reactive({
    req(input$normal_file, input$tumor_file)
    normal_df <- read_excel(input$normal_file$datapath)
    tumor_df <- read_excel(input$tumor_file$datapath)
    bind_rows(normal_df %>% mutate(Type = "Normal"),
              tumor_df %>% mutate(Type = "Tumor"))
  })

  observeEvent(dataInput(), {
    gene_choices <- unique(dataInput()$Target)
    updatePickerInput(session, "ref_gene", choices = gene_choices)
    updatePickerInput(session, "internal_controls", choices = gene_choices)
  })

  analysis <- eventReactive(input$analyze, {
    df <- dataInput()
    req(input$ref_gene)

    df <- df %>%
      group_by(Sample, Target, Type) %>%
      summarize(Cq = mean(Cq, na.rm = TRUE), .groups = 'drop') %>%
      pivot_wider(names_from = Target, values_from = Cq)

    df$Ct_ref <- rowMeans(df[, input$ref_gene], na.rm = TRUE)

    targets <- setdiff(names(df), c("Sample", "Type", input$ref_gene, "Ct_ref"))

    for (target in targets) {
      df[[paste0("dCt_", target)]] <- df[[target]] - df$Ct_ref
    }

    normal_df <- df %>% filter(Type == "Normal")
    tumor_df <- df %>% filter(Type == "Tumor")

    means_normal <- sapply(targets, function(t) mean(normal_df[[paste0("dCt_", t)]], na.rm = TRUE))

    for (target in targets) {
      normal_df[[paste0("ddCt_", target)]] <- normal_df[[paste0("dCt_", target)]] - means_normal[[target]]
      normal_df[[paste0("FC_", target)]] <- 2^(-normal_df[[paste0("ddCt_", target)]])

      tumor_df[[paste0("ddCt_", target)]] <- tumor_df[[paste0("dCt_", target)]] - means_normal[[target]]
      tumor_df[[paste0("FC_", target)]] <- 2^(-tumor_df[[paste0("ddCt_", target)]])
    }

    list(normal = normal_df, tumor = tumor_df, targets = targets)
  })

  output$normal_tables <- renderUI({
    req(analysis())
    targets <- analysis()$targets
    normal_df <- analysis()$normal
    tables <- lapply(targets, function(target) {
      datatable(normal_df[, c("Sample", input$ref_gene, target,
                              paste0("dCt_", target),
                              paste0("ddCt_", target),
                              paste0("FC_", target))],
                options = list(scrollX = TRUE),
                caption = htmltools::tags$caption(
                  style = 'caption-side: top; text-align: left;',
                  paste("Normal Tissue -", target)))
    })
    do.call(tagList, tables)
  })

  output$tumor_tables <- renderUI({
    req(analysis())
    targets <- analysis()$targets
    tumor_df <- analysis()$tumor
    tables <- lapply(targets, function(target) {
      datatable(tumor_df[, c("Sample", input$ref_gene, target,
                            paste0("dCt_", target),
                            paste0("ddCt_", target),
                            paste0("FC_", target))],
                options = list(scrollX = TRUE),
                caption = htmltools::tags$caption(
                  style = 'caption-side: top; text-align: left;',
                  paste("Tumor Tissue -", target)))
    })
    do.call(tagList, tables)
  })

  output$internal_control_plot <- renderPlot({
    req(dataInput(), input$internal_controls)
    df <- dataInput()
    internal_controls <- df %>% filter(Target %in% input$internal_controls)

    ggplot(internal_controls, aes(x = Sample, y = Cq, fill = Target)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      theme_minimal() +
      labs(title = "Ct of Internal Controls across Samples",
           x = "Sample",
           y = "Ct Value") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })

  output$mirna_plots <- renderUI({
    req(dataInput())
    df <- dataInput()
    mirna_targets <- setdiff(unique(df$Target), input$internal_controls)

    plots <- lapply(mirna_targets, function(target) {
      plot_data <- df %>% filter(Target == target)

      plot_output <- renderPlot({
        ggplot(plot_data, aes(x = Sample, y = Cq, fill = Type)) +
          geom_bar(stat = "identity", position = position_dodge()) +
          theme_minimal() +
          labs(title = paste("Ct Values for", target), x = "Sample", y = "Ct") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      })

      plotOutput(outputId = paste0("mirna_plot_", target))
    })
    do.call(tagList, plots)
  })

  observe({
    req(dataInput())
    df <- dataInput()
    mirna_targets <- setdiff(unique(df$Target), input$internal_controls)

    for (target in mirna_targets) {
      local({
        local_target <- target
        output[[paste0("mirna_plot_", local_target)]] <- renderPlot({
          plot_data <- df %>% filter(Target == local_target)
          ggplot(plot_data, aes(x = Sample, y = Cq, fill = Type)) +
            geom_bar(stat = "identity", position = position_dodge()) +
            theme_minimal() +
            labs(title = paste("Ct Values for", local_target), x = "Sample", y = "Ct") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        })
      })
    }
  })

  output$download_results <- downloadHandler(
    filename = function() {
      paste0("ddct_results_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      result <- analysis()
      wb <- createWorkbook()
      addWorksheet(wb, "Normal")
      writeData(wb, "Normal", result$normal)
      addWorksheet(wb, "Tumor")
      writeData(wb, "Tumor", result$tumor)
      saveWorkbook(wb, file)
    }
  )

  output$download_plot_excel <- downloadHandler(
    filename = function() {
      paste0("Internal_Controls_Plot_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      req(dataInput(), input$internal_controls)
      df <- dataInput()
      internal_controls <- df %>% filter(Target %in% input$internal_controls)

      p <- ggplot(internal_controls, aes(x = Sample, y = Cq, fill = Target)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        theme_minimal() +
        labs(title = "Ct of Internal Controls across Samples",
             x = "Sample",
             y = "Ct Value") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

      tmpfile <- tempfile(fileext = ".png")
      ggsave(tmpfile, plot = p, width = 8, height = 5)

      wb <- createWorkbook()
      addWorksheet(wb, "Internal Controls Plot")
      insertImage(wb, sheet = 1, tmpfile, startRow = 1, startCol = 1, width = 15, height = 6)
      saveWorkbook(wb, file)
    }
  )

  output$download_mirna_plots_excel <- downloadHandler(
    filename = function() {
      paste0("miRNA_Ct_Plots_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      df <- dataInput()
      mirna_targets <- setdiff(unique(df$Target), input$internal_controls)

      wb <- createWorkbook()

      for (target in mirna_targets) {
        plot_data <- df %>% filter(Target == target)

        p <- ggplot(plot_data, aes(x = Sample, y = Cq, fill = Type)) +
          geom_bar(stat = "identity", position = position_dodge()) +
          theme_minimal() +
          labs(title = paste("Ct Values for", target), x = "Sample", y = "Ct") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))

        tmpfile <- tempfile(fileext = ".png")
        ggsave(tmpfile, plot = p, width = 8, height = 5)

        addWorksheet(wb, target)
        insertImage(wb, sheet = target, file = tmpfile, startRow = 1, startCol = 1, width = 15, height = 6)
      }

      saveWorkbook(wb, file)
    }
  )
}

shinyApp(ui, server)
