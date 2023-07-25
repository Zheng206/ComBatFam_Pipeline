require(shiny)
require(DT)
require(ggplot2)
require(bslib)
require(bsicons)
require(shinydashboard)

#' Batch Effect Interactive Visualization
#'
#' Provide interactive batch/site effect visualization through shiny app.
#'
#' @param result A list derived from `visual_prep()` that contains datasets for shiny visualization.
#'
#' @import ggplot2
#' @import shiny
#' @import bsicons
#' @import shinydashboard
#' @import bslib 
#' @importFrom DT datatable formatStyle styleEqual DTOutput renderDT
#' @importFrom stats reorder
#' 
#' @export

comfam_shiny = function(result){
  info = result$info
  type = info$type
  df = info$df
  batch = info$batch
  features = info$features
  covariates = info$cov_shiny
  char_var = info$char_var
  num_var = setdiff(covariates, char_var)
  ui = function(request) {
    fluidPage(
      theme = bslib::bs_theme(version = 4, bootswatch = "minty"),
      titlePanel("Batch Effect Visualization"),
      sidebarLayout(
        sidebarPanel(
          conditionalPanel(condition="input.tabselected==2",
                           radioButtons("type", "Select output type", choices = c("Plot", "Table"), selected = "Plot"),
                           #radioButtons("cov_type", "Select covariates type", choices = c("numerical", "categorical"), selected = "numerical"),
                           selectInput("cov", "Select covariate", choices = covariates, selected = covariates[1])
          ),
          conditionalPanel(condition="input.tabselected==3",
                           selectInput("feature", "Select Feature", choices = features, selected = features[1])
          ),
          conditionalPanel(condition="input.tabselected==4",
                           selectInput("PC1", "Select the first PC", choices = colnames(result$pr.feature$x), selected = colnames(result$pr.feature$x)[1]),
                           selectInput("PC2", "Select the second PC", choices = colnames(result$pr.feature$x), selected = colnames(result$pr.feature$x)[2])
          ),
          conditionalPanel(condition="input.tabselected==5",
                           selectInput("com_type", "Select ComBatFam type", choices = c("comfam", "covfam"), selected = "comfam"),
                           selectInput("batch_selection", "Select the batches to be shown on the graph", choices = c("All", levels(result$info$df[[batch]])), selected = "All"),
                           #uiOutput("eb")
          ),
          conditionalPanel(condition="input.tabselected==6",
                           radioButtons("test_batch", "Batch Effect Test", choices = c("MDMR", "Kenward-Roger (liner mixed model)", "ANOVA", "Kruskal-Wallis"), selected = "MDMR"),
                           radioButtons("test_variance", "Equality of Variance Test", choices = c("Fligner-Killeen", "Levene's Test", "Bartlett's Test"), selected = "Fligner-Killeen")
          ),
          bookmarkButton()),
        mainPanel(
          tabsetPanel(
            tabPanel("Summary", value = 2, 
                     fluidRow(
                       shinydashboard::box(
                        width = 12,
                        title = "Batch Sample Size Summary",
                       shiny::uiOutput("output"))),
                     fluidRow(
                       shinydashboard::box(
                         width = 12,
                         title = "Covariate Distribution",
                         shiny::uiOutput("cov_output")))),
            tabPanel("Residual Plot", value = 3, 
                     fluidRow(
                       shinydashboard::box(
                         width = 12,
                         title = "Additive Batch Effect",
                         shiny::plotOutput("res_add"))),
                     fluidRow(
                       shinydashboard::box(
                         width = 12,
                         title = "Multiplicative Batch Effect",
                         shiny::plotOutput("res_ml")))),
            tabPanel("Residual Dimensionality Reduction", value = 4, 
                     fluidRow(
                       shinydashboard::box(
                         width = 12,
                         title = "PCA",
                         shiny::plotOutput("pca"))),
                     fluidRow(
                       shinydashboard::box(
                         width = 12,
                         title = "T-SNE",
                         shiny::plotOutput("tsne")))),
            tabPanel("Empirical Distribution", value = 5, 
                     fluidRow(
                       shinydashboard::box(
                         width = 12,
                         title = "Location Parameters",
                         shiny::plotOutput("eb_location"))),
                     fluidRow(
                       shinydashboard::box(
                         width = 12,
                         title = "Scale Paramaters",
                         shiny::plotOutput("eb_scale")))),
            tabPanel("Statistical Test", value = 6,
                     fluidRow(
                       shinydashboard::box(
                         width = 12,
                         title = "Batch Effect Test",
                         DT::DTOutput("test_batch"),
                         shiny::textOutput("sig_pct_batch"))),
                     fluidRow(
                       shinydashboard::box(
                         width = 12,
                         title = "Equality of Variance Test",
                         DT::DTOutput("test_variance"),
                         shiny::textOutput("sig_pct_variance")))),
            id = "tabselected"
          )   
        )
      )
    )
  }
  
  server = function(input, output, session) {
    output$output = shiny::renderUI({
      if (input$type == "Plot") {
        plotOutput("plot")
      } else if (input$type == "Table") {
        DT::DTOutput("table")
      }
    })
    output$plot = shiny::renderPlot({
      ggplot(result$summary_df, aes(x = count, y = eval(parse(text = batch)))) +
        geom_bar(stat = "identity", fill = "aquamarine") +
        geom_text(aes(label = count), hjust = 1.5, position = position_dodge(0.9), size = 3, colour = "black") +
        labs(x = "Count", y = "Batch") +
        theme(plot.title = element_text(hjust = 0.5))
    })
    output$table = DT::renderDT({
      result$summary_df %>% mutate(`percentage (%)` = sprintf("%.3f", `percentage (%)`)) %>% DT::datatable()
    })
    output$cov_output = shiny::renderUI({
      if (input$type == "Plot") {
        plotOutput("cov_plot")
      } else if (input$type == "Table") {
        DT::DTOutput("cov_table")
      }
    })
    output$cov_plot = shiny::renderPlot({
      if(input$cov %in% num_var){
        ggplot(df, aes(x = eval(parse(text = input$cov)), y = reorder(as.factor(eval(parse(text = batch))), eval(parse(text = input$cov)), Fun = median), fill = eval(parse(text = batch))))+
          geom_boxplot() +
          labs(x = input$cov, y = "Batch", fill = "Covariate") +
          theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
      }else if(input$cov %in% char_var){
        df_c = df %>% group_by(eval(parse(text = batch)), eval(parse(text = input$cov))) %>% dplyr::tally() %>% mutate(percentage = n/sum(n))
        colnames(df_c) = c(batch, input$cov, "n", "percentage")
        ggplot(df_c, aes(y = as.factor(eval(parse(text = batch))), x = n, fill = eval(parse(text = input$cov)))) +
          geom_bar(stat="identity", position ="fill") +
          geom_text(aes(label = paste0(sprintf("%1.1f", percentage*100),"%")), position = position_fill(vjust=0.5), colour="black", size = 3) +
          scale_fill_brewer(palette = "Pastel1") +
          labs(x = "Percentage", y = "Batch", title = "Covariate Distribution", fill = input$cov) +
          theme(plot.title = element_text(hjust = 0.5)) 
      }
    })
    output$cov_table =  DT::renderDT({
      if(input$cov %in% num_var){
        cov_summary_table = df %>% group_by(eval(parse(text = batch))) %>% summarize(min = min(eval(parse(text = input$cov))), mean = mean(eval(parse(text = input$cov))), max = max(eval(parse(text = input$cov))))
        colnames(cov_summary_table) = c(batch, "min", "mean", "max")
        cov_summary_table = cov_summary_table %>% mutate(mean = round(mean, 3))
        cov_summary_table %>% DT::datatable()
      }else if(input$cov %in% char_var){
        cov_summary_table = df %>% group_by(eval(parse(text = batch)), eval(parse(text = input$cov))) %>% dplyr::tally() %>% mutate(percentage = 100*n/sum(n))
        colnames(cov_summary_table) = c(batch, input$cov, "n", "percentage (%)")
        cov_summary_table %>% mutate(`percentage (%)` = sprintf("%.3f", `percentage (%)`)) %>% DT::datatable()
      }
    })
    output$res_add = shiny::renderPlot({
      add_mean = result$residual_add_df %>% group_by(result$residual_add_df[[batch]]) %>% summarize(across(features, median, .names = "mean_{.col}")) %>% ungroup()
      colnames(add_mean) = c(batch, colnames(add_mean)[-1])
      result$residual_add_df = result$residual_add_df %>% left_join(add_mean, by = c(batch))
      ggplot(result$residual_add_df, aes(x = reorder(as.factor(eval(parse(text = batch))), result$residual_add_df[[paste0("mean_",input$feature)]]), y = result$residual_add_df[[input$feature]])) +
        geom_boxplot() +
        geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
        labs(x = "Batch", y = "Residual")
    })
    output$res_ml = shiny::renderPlot({
      ml_mean = result$residual_ml_df %>% group_by(result$residual_ml_df[[batch]]) %>% summarize(across(features, median, .names = "mean_{.col}")) %>% ungroup()
      colnames(ml_mean) = c(batch, colnames(ml_mean)[-1])
      result$residual_ml_df = result$residual_ml_df %>% left_join(ml_mean, by = c(batch))
      ggplot(result$residual_ml_df, aes(x = reorder(as.factor(eval(parse(text = batch))), result$residual_ml_df[[paste0("mean_",input$feature)]]), y = result$residual_ml_df[[input$feature]])) +
        geom_boxplot() +
        geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
        labs(x = "Batch", y = "Residual")
    })
    output$pca = shiny::renderPlot({
      ggplot(result$pca_df, aes(x = result$pca_df[[input$PC1]], y = result$pca_df[[input$PC2]], color = result$pca_df[[batch]])) +
        geom_point() +
        labs(x = input$PC1, y = input$PC2, color = "Batch")
    })
    output$tsne = shiny::renderPlot({
      ggplot(result$tsne_df, aes(x = cor_1, y = cor_2, color = result$tsne_df[[batch]])) +
        geom_point() +
        labs(x = "Dim 1", y = "Dim 2", color = "Batch")
    })
    #output$eb = shiny::renderUI({
    #  if(input$com_type == "covfam"){
    #    radioButtons("eb_type", "Select the type of covfam empirical parameters", choices = c("com", "com.scores"), selected = "com")
    #  }
    #})
    output$eb_location = shiny::renderPlot({
      if(input$com_type == "comfam"){
        min_x = result$eb_df %>% filter(grepl("^gamma_*", type)) %>% pull(eb_values) %>% min()
        max_x = result$eb_df %>% filter(grepl("^gamma_*", type)) %>% pull(eb_values) %>% max()
        if(input$batch_selection == "All"){
          ggplot(result$eb_df %>% filter(grepl("^gamma_*", type)), aes(x = eb_values, color = batch, linetype = type)) +
            geom_density() +
            xlim(min_x, max_x) +
            labs(x = "Gamma", y = "Density", color = "Batch", linetype = "Estimate Type")
        }else{
          ggplot(result$eb_df %>% filter(grepl("^gamma_*", type), batch == input$batch_selection), aes(x = eb_values, linetype = type)) +
            geom_density() +
            xlim(min_x, max_x) +
            labs(x = "Gamma", y = "Density", linetype = "Estimate Type")
          }
      }else if(input$com_type == "covfam"){
        if(input$batch_selection == "All"){
            min_x = result$eb_df %>% filter(grepl("^gamma_*", type)) %>% pull(eb_values) %>% min()
            max_x = result$eb_df %>% filter(grepl("^gamma_*", type)) %>% pull(eb_values) %>% max()
            ggplot(result$eb_df %>% filter(grepl("^gamma_*", type)), aes(x = eb_values, color = batch, linetype = type)) +
              geom_density() +
              xlim(min_x, max_x) +
              labs(x = "Gamma", y = "Density", color = "Batch", linetype = "Estimate Type")
        }else{
            min_x = result$eb_df %>% filter(grepl("^gamma_*", type)) %>% pull(eb_values) %>% min()
            max_x = result$eb_df %>% filter(grepl("^gamma_*", type)) %>% pull(eb_values) %>% max()
            ggplot(result$eb_df %>% filter(grepl("^gamma_*", type), batch == input$batch_selection), aes(x = eb_values, linetype = type)) +
              geom_density() +
              xlim(min_x, max_x) +
              labs(x = "Gamma", y = "Density", linetype = "Estimate Type")
        }
      }
    })
    output$eb_scale = shiny::renderPlot({
      min_x = result$eb_df %>% filter(grepl("^delta_*", type)) %>% pull(eb_values) %>% min()
      max_x = result$eb_df %>% filter(grepl("^delta_*", type)) %>% pull(eb_values) %>% max()
      if(input$com_type == "comfam"){
        if(input$batch_selection == "All"){
          ggplot(result$eb_df %>% filter(grepl("^delta_*", type)), aes(x = eb_values, color = batch, linetype = type)) +
            geom_density() +
            xlim(min_x, max_x) +
            labs(x = "Delta", y = "Density", color = "Batch", linetype = "Estimate Type")
        }else{
          ggplot(result$eb_df %>% filter(grepl("^delta_*", type), batch == input$batch_selection), aes(x = eb_values, linetype = type)) +
            geom_density() +
            xlim(min_x, max_x) +
            labs(x = "Delta", y = "Density", linetype = "Estimate Type")
        }
      }else if(input$com_type == "covfam"){
        if(input$batch_selection == "All"){
            min_x = result$eb_df %>% filter(grepl("^delta_*", type)) %>% pull(eb_values) %>% min()
            max_x = result$eb_df %>% filter(grepl("^delta_*", type)) %>% pull(eb_values) %>% max()
            ggplot(result$eb_df %>% filter(grepl("^delta_*", type)), aes(x = eb_values, color = batch, linetype = type)) +
              geom_density() +
              xlim(min_x, max_x) +
              labs(x = "Delta", y = "Density", color = "Batch", linetype = "Estimate Type")
        }else{
            min_x = result$eb_df %>% filter(grepl("^delta_*", type)) %>% pull(eb_values) %>% min()
            max_x = result$eb_df %>% filter(grepl("^delta_*", type)) %>% pull(eb_values) %>% max()
            ggplot(result$eb_df %>% filter(grepl("^delta_*", type), batch == input$batch_selection), aes(x = eb_values, linetype = type)) +
              geom_density() +
              xlim(min_x, max_x) +
              labs(x = "Delta", y = "Density", linetype = "Estimate Type")
        }
      }
    })
    output$test_batch = DT::renderDT({
      if(input$test_batch == "Kenward-Roger (liner mixed model)"){
        if(type == "lmer"){
          result$kr_test_df %>% datatable(options = list(columnDefs = list(list(className = 'dt-center', 
                                                                              targets = "_all")))) %>% formatStyle(columns = c("p.value"),color = styleEqual(result$red, "red")) %>% formatStyle(
                                                                                'sig',
                                                                                target = 'row',
                                                                                backgroundColor = styleEqual(c("*", "**", "***"), "lightyellow")
                                                                              )}else{
                                                                                result$kr_test_df %>% DT::datatable()
                                                                              }
        }else if(input$test_batch == "MDMR"){
                                                                                result$mdmr.summary %>% DT::datatable(options = list(columnDefs = list(list(className = 'dt-center', 
                                                                                                                                                            targets = "_all")))) %>% formatStyle(columns = c("p.value"),color = styleEqual(result$red, "red")) %>% formatStyle(
                                                                                                                                                              'sig',
                                                                                                                                                              target = 'row',
                                                                                                                                                              backgroundColor = styleEqual(c("*", "**", "***"), "lightyellow"))
                                                                              }else if(input$test_batch== "ANOVA"){
                                                                                result$anova_test_df %>% datatable(options = list(columnDefs = list(list(className = 'dt-center', 
                                                                                                                                                         targets = "_all")))) %>% formatStyle(columns = c("p.value"),color = styleEqual(result$red, "red")) %>% formatStyle(
                                                                                                                                                           'sig',
                                                                                                                                                           target = 'row',
                                                                                                                                                           backgroundColor = styleEqual(c("*", "**", "***"), "lightyellow"))
                                                                              }else if(input$test_batch == "Kruskal-Wallis"){
                                                                                result$kw_test_df %>% datatable(options = list(columnDefs = list(list(className = 'dt-center', 
                                                                                                                                                         targets = "_all")))) %>% formatStyle(columns = c("p.value"),color = styleEqual(result$red, "red")) %>% formatStyle(
                                                                                                                                                           'sig',
                                                                                                                                                           target = 'row',
                                                                                                                                                           backgroundColor = styleEqual(c("*", "**", "***"), "lightyellow"))
                                                                              }
    })
    output$test_variance = DT::renderDT({
      if(input$test_variance == "Fligner-Killeen"){
        result$fk_test_df %>% datatable(options = list(columnDefs = list(list(className = 'dt-center', 
                                                                              targets = "_all")))) %>% formatStyle(columns = c("p.value"),color = styleEqual(result$red, "red")) %>% formatStyle(
                                                                                'sig',
                                                                                target = 'row',
                                                                                backgroundColor = styleEqual(c("*", "**", "***"), "lightyellow"))
      }else if(input$test_variance == "Levene's Test"){
        result$lv_test_df %>% datatable(options = list(columnDefs = list(list(className = 'dt-center', 
                                                                              targets = "_all")))) %>% formatStyle(columns = c("p.value"),color = styleEqual(result$red, "red")) %>% formatStyle(
                                                                                'sig',
                                                                                target = 'row',
                                                                                backgroundColor = styleEqual(c("*", "**", "***"), "lightyellow"))
      }else if(input$test_variance == "Bartlett's Test"){
        result$bl_test_df %>% datatable(options = list(columnDefs = list(list(className = 'dt-center', 
                                                                              targets = "_all")))) %>% formatStyle(columns = c("p.value"),color = styleEqual(result$red, "red")) %>% formatStyle(
                                                                                'sig',
                                                                                target = 'row',
                                                                                backgroundColor = styleEqual(c("*", "**", "***"), "lightyellow"))
      }
    })
    output$sig_pct_batch = shiny::renderText({
      if(input$test_batch == "Kenward-Roger (liner mixed model)"){
        if(type == "lmer"){
          n = nrow(result$kr_test_df)
          pct = 100 * (n - sum(is.na(result$kr_test_df$sig)))/n
          print(paste0("The percentage of significant features is: ", round(pct,2), "%."))}else{
            print("The Kenward-Roger test is a modification of the degrees of freedom in linear mixed models. Not appropriate for the current type of model.")
          }
      }else if(input$test_batch == "MDMR"){
        if(is.na(result$mdmr.summary$sig[2])){
          print("The batch effect seems not significant.")
        }else{print("The batch effect seems significant.")}
      }else if(input$test_batch == "ANOVA"){
        n = nrow(result$anova_test_df)
        pct = 100 * (n - sum(is.na(result$anova_test_df$sig)))/n
        print(paste0("The percentage of significant features is: ", round(pct,2), "%."))
      }
    })
    output$sig_pct_variance = shiny::renderText({
      if(input$test_variance == "Fligner-Killeen"){
        n = nrow(result$fk_test_df)
        pct = 100 * (n - sum(is.na(result$fk_test_df$sig)))/n
        print(paste0("The percentage of significant features is: ", round(pct,2), "%."))
      }else if(input$test_variance == "Levene's Test"){
        n = nrow(result$lv_test_df)
        pct = 100 * (n - sum(is.na(result$lv_test_df$sig)))/n
        print(paste0("The percentage of significant features is: ", round(pct,2), "%."))
      }
    })
  }
  
  shinyApp(ui = ui, server = server, enableBookmarking = "url")
}

utils::globalVariables(c("count", "cor_1", "cor_2", "percentage", "batch", "features", "eb_values", "type", "percentage (%)"))

