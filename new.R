library(shiny)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)
library(shinyWidgets)
library(DT)
library(rmarkdown)
library(knitr)

ui <- fluidPage(
  titlePanel("ΔΔCt Analysis Tool"),

  sidebarLayout(
    sidebarPanel(
      fileInput("normal_file", "Upload Normal Sample Excel", accept = ".xlsx"),
      fileInput("tumor_file", "Upload Tumor Sample Excel", accept = ".xlsx"),
      pickerInput("ref_genes", "Select Reference Genes", choices = NULL, multiple = TRUE),
      actionButton("analyze", "Analyze"),
      downloadButton("download_report", "Download HTML Report")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("ΔCt Tables",
                 h4("Normal Samples"), DTOutput("normal_table"),
                 h4("Tumor Samples"), DTOutput("tumor_table")),
        tabPanel("ΔΔCt Table", DTOutput("comparison_table")),
        tabPanel("Fold Change Plot", plotOutput("fold_change_plot")),
        tabPanel("ΔΔCt SD Plot", plotOutput("ddct_sd_plot")),
        tabPanel("ΔCt Plot", plotOutput("dct_plot")),
        tabPanel("Sample-level Dysregulation", 
                 DTOutput("sample_ddCt_table"))
      )
    )
  )
)

server <- function(input, output, session) {
  dataList <- reactiveValues()

  observeEvent(c(input$normal_file, input$tumor_file), {
    req(input$normal_file, input$tumor_file)
    normal <- read_excel(input$normal_file$datapath)
    tumor <- read_excel(input$tumor_file$datapath)

    required_cols <- c("Target", "Sample", "Cq")
    if (!all(required_cols %in% names(normal)) || !all(required_cols %in% names(tumor))) {
      showModal(modalDialog("Each file must contain columns: Target, Sample, Cq", easyClose = TRUE))
      return()
    }

    updatePickerInput(session, "ref_genes", choices = unique(normal$Target))

    dataList$normal_raw <- normal
    dataList$tumor_raw <- tumor
  })

  observeEvent(input$analyze, {
    req(input$ref_genes)

    process_data <- function(df) {
      df <- df %>%
        filter(!is.na(Cq)) %>%
        group_by(Target, Sample) %>%
        summarise(Cq = mean(Cq), .groups = "drop") %>%
        pivot_wider(names_from = Target, values_from = Cq)

      ref_mean <- df %>%
        select(all_of(input$ref_genes)) %>%
        rowMeans(na.rm = TRUE)

      df$Ct_ref <- ref_mean
      targets <- setdiff(names(df), c("Sample", input$ref_genes, "Ct_ref"))

      for (target in targets) {
        df[[paste0("dCt_", target)]] <- df[[target]] - df$Ct_ref
      }
      df
    }

    normal_dCt <- process_data(dataList$normal_raw)
    tumor_dCt <- process_data(dataList$tumor_raw)



    comparison <- data.frame(Target = setdiff(names(normal_dCt), c("Sample", input$ref_genes, "Ct_ref")))
    comparison$Target <- gsub("dCt_", "", comparison$Target)

    comparison <- comparison %>%
      rowwise() %>%
      mutate(
        mean_dCt_normal = mean(normal_dCt[[paste0("dCt_", Target)]], na.rm = TRUE),
        mean_dCt_tumor = mean(tumor_dCt[[paste0("dCt_", Target)]], na.rm = TRUE),
        ddCt = mean_dCt_tumor - mean_dCt_normal,
        ddCt_sd = sd(tumor_dCt[[paste0("dCt_", Target)]] - normal_dCt[[paste0("dCt_", Target)]], na.rm = TRUE),
        fold_change = 2^(-ddCt),
        regulation = case_when(
          fold_change > 1.5 ~ "Upregulated",
          fold_change < 0.66 ~ "Downregulated",
          TRUE ~ "No Change"
        )
      ) %>%
      ungroup()

    normal_dCt_long <- normal_dCt %>%
      select(Sample, starts_with("dCt_")) %>%
      pivot_longer(-Sample, names_to = "Target", values_to = "dCt") %>%
      mutate(Target = gsub("dCt_", "", Target), Group = "Normal")

    tumor_dCt_long <- tumor_dCt %>%
      select(Sample, starts_with("dCt_")) %>%
      pivot_longer(-Sample, names_to = "Target", values_to = "dCt") %>%
      mutate(Target = gsub("dCt_", "", Target), Group = "Tumor")

baseline_normal <- normal_dCt %>%
  select(starts_with("dCt_")) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  pivot_longer(everything(), names_to = "Target", values_to = "baseline_dCt") %>%
  mutate(Target = gsub("dCt_", "", Target))

all_samples <- bind_rows(
  normal_dCt %>% mutate(Group = "Normal"),
  tumor_dCt %>% mutate(Group = "Tumor")
)

all_samples_long <- all_samples %>%
  select(Sample, Group, starts_with("dCt_")) %>%
  pivot_longer(cols = starts_with("dCt_"), names_to = "Target", values_to = "dCt") %>%
  mutate(Target = gsub("dCt_", "", Target))

sample_ddCt <- all_samples_long %>%
  left_join(baseline_normal, by = "Target") %>%
  mutate(
    ddCt_sample = dCt - baseline_dCt,
    fold_change_sample = 2^(-ddCt_sample),
    regulation_sample = case_when(
      fold_change_sample > 1.5 ~ "Upregulated",
      fold_change_sample < 0.66 ~ "Downregulated",
      TRUE ~ "No Change"
    )
  )

# Save all data to reactiveValues
dataList$normal <- normal_dCt
dataList$tumor <- tumor_dCt
dataList$comparison <- comparison
dataList$plot_data <- bind_rows(normal_dCt_long, tumor_dCt_long)
dataList$sample_ddCt <- sample_ddCt
  })

output$sample_ddCt_table <- renderDT({
  req(dataList$sample_ddCt)
  datatable(dataList$sample_ddCt, filter = "top",
            options = list(pageLength = 15, scrollX = TRUE))
})

  output$normal_table <- renderDT({
    req(dataList$normal)
    datatable(dataList$normal)
  })

  output$tumor_table <- renderDT({
    req(dataList$tumor)
    datatable(dataList$tumor)
  })

  output$comparison_table <- renderDT({
    req(dataList$comparison)
    datatable(dataList$comparison)
  })

  output$fold_change_plot <- renderPlot({
    req(dataList$comparison)
    ggplot(dataList$comparison, aes(x = reorder(Target, fold_change), y = fold_change, fill = regulation)) +
      geom_col() +
      coord_flip() +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
      labs(title = "Fold Change (Tumor vs Normal)", y = "Fold Change (2^-ΔΔCt)", x = "miRNA") +
      theme_minimal() +
      scale_fill_manual(values = c("Upregulated" = "#E74C3C", "Downregulated" = "#3498DB", "No Change" = "gray"))
  })

  output$ddct_sd_plot <- renderPlot({
    req(dataList$comparison)
    ggplot(dataList$comparison, aes(x = reorder(Target, ddCt_sd), y = ddCt_sd)) +
      geom_col(fill = "purple") +
      coord_flip() +
      labs(title = "ΔΔCt Standard Deviation", y = "SD of ΔΔCt", x = "miRNA") +
      theme_minimal()
  })

  output$dct_plot <- renderPlot({
    req(dataList$plot_data)
    ggplot(dataList$plot_data, aes(x = Target, y = dCt, fill = Group)) +
      geom_boxplot() +
      labs(title = "ΔCt Comparison Between Groups", x = "miRNA", y = "ΔCt") +
      theme_minimal()
  })

  output$download_report <- downloadHandler(
    filename = function() { "ddct_report.html" },
    content = function(file) {
      tempReport <- tempfile(fileext = ".Rmd")
      reportContent <- "---\ntitle: 'ΔΔCt Expression Report'\noutput: html_document\n---\n\n```{r setup, include=FALSE}\nlibrary(dplyr)\nlibrary(ggplot2)\nnormal <- dataList$normal\ntumor <- dataList$tumor\ncomparison <- dataList$comparison\nplot_data <- dataList$plot_data\n```\n\n# ΔCt Tables\n\n## Normal Samples\n```{r}\nknitr::kable(normal)\n```\n\n## Tumor Samples\n```{r}\nknitr::kable(tumor)\n```\n\n# Comparison Table\n```{r}\nknitr::kable(comparison)\n```\n\n# Fold Change Plot\n```{r}\nggplot(comparison, aes(x = reorder(Target, fold_change), y = fold_change, fill = regulation)) +\n  geom_col() +\n  coord_flip() +\n  geom_hline(yintercept = 1, linetype = 'dashed', color = 'gray40') +\n  labs(title = 'Fold Change (Tumor vs Normal)', y = 'Fold Change (2^-ΔΔCt)', x = 'miRNA') +\n  theme_minimal() +\n  scale_fill_manual(values = c('Upregulated' = '#E74C3C', 'Downregulated' = '#3498DB', 'No Change' = 'gray'))\n```\n\n# ΔΔCt SD Plot\n```{r}\nggplot(comparison, aes(x = reorder(Target, ddCt_sd), y = ddCt_sd)) +\n  geom_col(fill = 'purple') +\n  coord_flip() +\n  labs(title = 'ΔΔCt Standard Deviation', y = 'SD of ΔΔCt', x = 'miRNA') +\n  theme_minimal()\n```\n\n# ΔCt Plot\n```{r}\nggplot(plot_data, aes(x = Target, y = dCt, fill = Group)) +\n  geom_boxplot() +\n  labs(title = 'ΔCt Comparison Between Groups', x = 'miRNA', y = 'ΔCt') +\n  theme_minimal()\n```"
      writeLines(reportContent, tempReport)
      rmarkdown::render(tempReport, output_file = file, envir = new.env(parent = globalenv()))
    }
  )
}

shinyApp(ui = ui, server = server)
