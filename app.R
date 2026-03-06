# ==============================================================================
# Limit of Detection (LoD) Analysis – ddPCR (Shiny App)
# ==============================================================================

library(shiny)
library(bslib)
library(tidyverse)
library(patchwork)
library(cowplot)
library(colourpicker)

# ================================ UI ==========================================

# 1. Define the documentation text (using markdown)

ui <- page_navbar(
  title = "ddPCR Limit of Detection (LoD) Calculator",
  header = tagList(
    withMathJax(),
    tags$style(HTML("
    h1 { font-size: 22px; }
    h2 { font-size: 18px; }
    h3 { font-size: 16px; }
    h4 { font-size: 14px; }
  "))
  ),
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  fillable = FALSE,
  
  # Custom CSS
  tags$head(
    tags$style(HTML("
      .table { width: 100% !important; }
      .table th, .table td { white-space: nowrap; padding-left: 15px !important; padding-right: 15px !important; }
      .table td:first-child { font-weight: bold; color: #2c3e50; }
    "))
  ),
  
  # ================= TAB 1: ANALYSIS =================
  nav_panel(
    title = "Analysis Dashboard",
    layout_sidebar(
      sidebar = sidebar(
        fileInput("file", "1. Upload Bio-Rad CSV", accept = c(".csv")),
        helpText("Sample Description 1 should contain the theoretical concentration in whichever unit you choose (gc/uL or reaction, or uL reaction). Used only for grouping."),
        h5("2. Analysis Parameters"),
        numericInput("pos_threshold", "Positive Droplet Threshold", value = 3, min = 0, step = 1),
        helpText("Concentrations with fewer droplets will be set to 0."),
        numericInput("acc_threshold", "Min Accepted Droplets", value = 10000, min = 0),
        numericInput("reaction_vol", "Reaction Volume (µL)", value = 20, min = 1),
        numericInput("lod_prob", "Detection Probability", value = 0.95, min = 0.1, max = 0.99, step = 0.01),
        
        h5("3. Model Selection"),
        radioButtons("method", "LoD Calculation Method", 
                     choices = c("Both", "Probit", "Weibull"), 
                     selected = "Both"),
        
        hr(),
        
        h5("4. Customize Plot Colors"),
        colourInput("probit_line_color", "Probit Line Color", value = "#1F78B4"),
        colourInput("probit_ribbon_color", "Probit CI Color", value = "#A6CEE3"),
        colourInput("weibull_line_color", "Weibull Line Color", value = "#33A02C"),
        colourInput("weibull_ribbon_color", "Weibull CI Color", value = "#B2DF8A"),
        checkboxInput("show_axis_titles", "Show Axis Titles", value = TRUE),
        
        h5("5. Download Outputs"),
        fluidRow(
          column(6, numericInput("plot_width", "Width (in)", value = 18)),
          column(6, numericInput("plot_height", "Height (in)", value = 12))
        ),
        numericInput("plot_dpi", "Plot DPI", value = 300),
        downloadButton("download_plot", "Download Figures", class = "btn-primary w-100 mb-2"),
        downloadButton("download_csv", "Download Results CSV", class = "btn-success w-100")
      ),
      
      # Main Content Area for Tab 1
      card(
        card_header("LoD Plots"),
        plotOutput("lod_plot", height = "500px"),
        height = "700px"
      ),
      layout_column_wrap(
        width = 1/2,
        card(card_header("Method Recommendation"), verbatimTextOutput("fit_recommendation"), height = "350px"),
        card(card_header("LoD Results Table"), tableOutput("lod_table"), height = "350px")
      ),
      card(
        card_header(
          div(class = "d-flex justify-content-between align-items-center",
              "Replicate Summary & Data Quality",
              downloadButton("download_replicate_summary", "Export CSV", class = "btn-sm btn-success"))
        ),
        card_body(tableOutput("replicate_summary_table")),
        full_screen = TRUE
      )
    )
  ),
  
  # ================= TAB 2: DOCUMENTATION =================
  nav_panel(
    title = "Documentation & Methodology",
    card(
      card_header("Statistical Methods and Assumptions"),
      card_body(
        includeMarkdown("docs.Rmd")
      )
    )
  )
)

# ============================== SERVER ========================================
server <- function(input, output, session) {
  
  # 0. REPLICATE SUMMARY TABLE
  
  # A. Make the summary data reactive so it can be shared by the table and the download handler
  replicate_summary_data <- reactive({
    req(input$file)
    
    # Read the raw data
    df_raw <- read_csv(input$file$datapath, show_col_types = FALSE)
    
    df_raw %>%
      filter(`Sample description 1` != 'NTC') %>%
      group_by(Target, `Sample description 1`) %>%
      summarise(
        `Total Reps` = n(),
        
        # Count wells that failed the droplet threshold
        `Omit (Low Drplts)` = sum(`Accepted Droplets` < input$acc_threshold),
        
        # Count wells that passed droplet threshold
        `Valid Reps` = sum(`Accepted Droplets` >= input$acc_threshold),
        
        # Measured Concentration (mean of valid replicates only)
        `Meas. Conc. (gc/rxn)` = mean(`Conc(copies/µL)`[`Accepted Droplets` >= input$acc_threshold] * input$reaction_vol, na.rm = TRUE),
        
        # Count 'Detected'
        `Detect` = sum(Positives >= input$pos_threshold & 
                                 `Accepted Droplets` >= input$acc_threshold),
        
        # Count 'Not Detected'
        `Not Detect` = sum(Positives < input$pos_threshold & 
                               `Accepted Droplets` >= input$acc_threshold),
        .groups = "drop"
      ) %>%
      mutate(
        # Calculate % Detection
        `Detect Rate (%)` = if_else(`Valid Reps` > 0, 
                                       (`Detect` / `Valid Reps`) * 100, 
                                       0)
      ) %>%
      # Format the numeric columns to 2 decimal places for a cleaner look
      mutate(
        `Meas. Conc. (gc/rxn)` = sprintf("%.2f", `Meas. Conc. (gc/rxn)`),
        `Detect Rate (%)` = sprintf("%.1f%%", `Detect Rate (%)`)
      ) %>%
      rename(`Theoretical Conc.` = `Sample description 1`)
  })
  
  # B. Render the table using the reactive data
  output$replicate_summary_table <- renderTable({
    replicate_summary_data()
  }, striped = TRUE, hover = TRUE, bordered = TRUE, align = 'c')
  
  
  # 1. Reactive Data Processing
  processed_data <- reactive({
    req(input$file) # Ensure file is uploaded
    
    # Read raw CSV
    df_raw <- read_csv(input$file$datapath, show_col_types = FALSE)
    
    # Filter NTCs and poorly accepted wells, apply zeroing logic
    df_clean <- df_raw %>%
      filter(`Sample description 1` != 'NTC') %>%
      filter(`Accepted Droplets` >= input$acc_threshold) %>%
      mutate(
        # Zero out concentrations if droplets are below user-defined threshold
        Conc_filtered = if_else(Positives < input$pos_threshold, 0, `Conc(copies/µL)`),
        detected_rep = if_else(Conc_filtered > 0, 1, 0),
        Conc_reaction = Conc_filtered * input$reaction_vol
      )
    
    # Aggregate data by target and concentration
    df_avg <- df_clean %>%
      group_by(Target, `Sample description 1`) %>%
      summarise(
        mean_conc = mean(Conc_reaction, na.rm = TRUE),
        sd_conc   = sd(Conc_reaction, na.rm = TRUE),
        n_reps    = n(),
        n_detected = sum(detected_rep),
        detection_probability = n_detected / n_reps,
        .groups = "drop"
      )
    
    return(df_avg)
  })
  
  # 2. Main Analysis Logic (Models, Plotting, Lack-of-Fit)
  analysis_results <- reactive({
    df_avg <- processed_data()
    targets <- unique(df_avg$Target)
    
    # Output containers
    results_table <- list()
    plot_list <- list()
    fit_texts <- c()
    
    # Progress bar for bootstrapping
    withProgress(message = 'Calculating LoD models...', value = 0, {
      
      for(i in seq_along(targets)) {
        tgt <- targets[i]
        df_tgt <- df_avg %>% filter(Target == tgt)
        
        incProgress(1/length(targets), detail = paste("Target:", tgt))
        
        # Adjust probabilities slightly for Weibull NLS
        epsilon <- 0.001
        df_tgt <- df_tgt %>%
          mutate(detection_probability_adj = pmin(pmax(detection_probability, epsilon), 1 - epsilon))
        
        # Sequence for smooth curve plotting
        pred_seq <- seq(0, max(df_tgt$mean_conc, na.rm=T) * 1.1, length.out = 200)
        pred_df <- tibble(mean_conc = pred_seq)
        
        # Initialize variables
        lod_probit <- NA; lod_probit_lwr <- NA; lod_probit_upr <- NA
        lod_weibull <- NA; lod_weibull_lwr <- NA; lod_weibull_upr <- NA
        rss_probit <- NA; rss_weibull <- NA
        
        # ================= PROBIT MODEL =================
        if (input$method %in% c("Probit", "Both")) {
          probit_model <- glm(detection_probability ~ mean_conc, 
                              data = df_tgt, 
                              family = binomial(link = "probit"), 
                              weights = n_reps)
          
          # Predictions
          pred_probit <- predict(probit_model, newdata = tibble(mean_conc = pred_seq), se.fit = TRUE, type = "link")
          pred_df <- pred_df %>% mutate(
            Probit = pnorm(pred_probit$fit),
            probit_lwr = pnorm(pred_probit$fit - 1.96 * pred_probit$se.fit),
            probit_upr = pnorm(pred_probit$fit + 1.96 * pred_probit$se.fit)
          )
          
          # LoD & Confidence Intervals (Delta Method)
          lod_probit <- (qnorm(input$lod_prob) - coef(probit_model)[1]) / coef(probit_model)[2]
          coef_se <- summary(probit_model)$coefficients[, "Std. Error"]
          se_lod <- sqrt((coef_se[1]/coef(probit_model)[2])^2 + 
                           ((qnorm(input$lod_prob)-coef(probit_model)[1])/coef(probit_model)[2]^2 * coef_se[2])^2)
          lod_probit_lwr <- lod_probit - 1.96 * se_lod
          lod_probit_upr <- lod_probit + 1.96 * se_lod
          
          # Lack of Fit: Residual Sum of Squares (RSS)
          obs_prob <- df_tgt$detection_probability
          pred_obs_probit <- pnorm(predict(probit_model, type = "link"))
          rss_probit <- sum((obs_prob - pred_obs_probit)^2)
        }
        
        # ================= WEIBULL MODEL =================
        if (input$method %in% c("Weibull", "Both")) {
          wb_start <- list(lambda = median(df_tgt$mean_conc[df_tgt$detection_probability_adj > 0]), k = 1)
          wb <- tryCatch({
            nls(detection_probability_adj ~ 1 - exp(- (mean_conc / lambda)^k),
                data = df_tgt, start = wb_start,
                algorithm = "port", lower = c(lambda=1e-6, k=0.1), upper = c(lambda=Inf, k=10))
          }, error = function(e) NULL)
          
          if(!is.null(wb)) {
            wb_coef <- coef(wb)
            lod_weibull <- wb_coef["lambda"] * (-log(1 - input$lod_prob))^(1/wb_coef["k"])
            
            # Predictions
            pred_df <- pred_df %>% mutate(
              Weibull = 1 - exp(- (mean_conc / wb_coef["lambda"])^wb_coef["k"])
            )
            
            # Bootstrapping for CI (condensed size to 100 for Shiny responsiveness)
            n_boot <- 100 
            set.seed(123)
            wb_boot_preds <- replicate(n_boot, {
              df_b <- df_tgt[sample(nrow(df_tgt), replace = TRUE), ]
              tryCatch({
                fit <- nls(detection_probability_adj ~ 1 - exp(- (mean_conc / lambda)^k),
                           data = df_b, start = wb_start, algorithm = "port",
                           lower = c(1e-6, 0.1), upper = c(Inf, 10))
                l <- coef(fit)["lambda"]; k <- coef(fit)["k"]
                # Return list of curve predictions and the LOD
                list(curve = 1 - exp(- (pred_seq / l)^k), lod = l * (-log(1 - input$lod_prob))^(1 / k))
              }, error = function(e) list(curve = rep(NA, length(pred_seq)), lod = NA))
            }, simplify = FALSE)
            
            # Extract CI bounds
            curves <- do.call(rbind, lapply(wb_boot_preds, `[[`, "curve"))
            lods <- sapply(wb_boot_preds, `[[`, "lod")
            
            pred_df <- pred_df %>% mutate(
              weibull_lwr = pmin(pmax(apply(curves, 2, quantile, probs=0.025, na.rm=T), 0), 1),
              weibull_upr = pmin(pmax(apply(curves, 2, quantile, probs=0.975, na.rm=T), 0), 1)
            )
            lod_weibull_lwr <- quantile(lods, 0.025, na.rm = TRUE)
            lod_weibull_upr <- quantile(lods, 0.975, na.rm = TRUE)
            
            # Lack of Fit (RSS)
            obs_prob <- df_tgt$detection_probability
            pred_obs_weibull <- 1 - exp(- (df_tgt$mean_conc / wb_coef["lambda"])^wb_coef["k"])
            rss_weibull <- sum((obs_prob - pred_obs_weibull)^2)
          }
        }
        
        # ================= LACK OF FIT RECOMMENDATION =================
        if (input$method == "Both") {
          if (!is.na(rss_probit) && !is.na(rss_weibull)) {
            diff <- rss_probit - rss_weibull
            best <- if(abs(diff) < 0.05) "Either model is fine (similar fit)." 
            else if (diff > 0) "Weibull is recommended (lower error)." 
            else "Probit is recommended (lower error)."
            
            fit_texts <- c(fit_texts, sprintf("%s: Probit RSS = %.3f | Weibull RSS = %.3f \n    -> %s", 
                                              tgt, rss_probit, rss_weibull, best))
          }
        }
        
        # ================= STORE RESULTS (Formatted for Table) =================
        results_table[[tgt]] <- tibble(
          Target = tgt,
          `Probit LoD [95% CI]` = if(!is.na(lod_probit)) {
            sprintf("%.2f [%.2f - %.2f]", lod_probit, lod_probit_lwr, lod_probit_upr)
          } else { "N/A" },
          
          `Weibull LoD [95% CI]` = if(!is.na(lod_weibull)) {
            sprintf("%.2f [%.2f - %.2f]", lod_weibull, lod_weibull_lwr, lod_weibull_upr)
          } else { "N/A" }
        )
        
        # ================= PLOTTING =================
        p <- ggplot() +
          geom_errorbarh(data = df_tgt, height = 0.02, alpha = 0.6,
                         aes(y = detection_probability, 
                             xmin = pmax(mean_conc - sd_conc, 0), 
                             xmax = mean_conc + sd_conc)) +
          geom_point(data = df_tgt, aes(x = mean_conc, y = detection_probability), size = 3, color = "grey40") +
          geom_hline(yintercept = input$lod_prob, linetype = "dashed", color = "gray50") +
          scale_x_continuous(limits = c(0, max(df_tgt$mean_conc, na.rm=TRUE)*1.2), breaks = scales::pretty_breaks(n=5)) +
          labs(
            title = tgt,
            x = if(input$show_axis_titles) "Average Concentration (gc/reaction)" else NULL,
            y = if(input$show_axis_titles) "Detection Probability" else NULL
          ) +
          theme_minimal(base_size = 14) +
          theme(panel.grid.minor = element_blank(),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
                plot.title = element_text(hjust = 0.5, face = "bold"),
                axis.title = element_text(size = 15), # Slightly smaller so they don't overlap
                axis.text.x = element_text(size = 13),       # X-axis tick labels
                axis.text.y = element_text(size = 13),       # Y-axis tick labels
                legend.position = "none"
          )
        
        # Add Probit specifics
        if (input$method %in% c("Both", "Probit") && !is.na(lod_probit)) {
          p <- p + 
            geom_ribbon(data = pred_df, aes(x = mean_conc, ymin = probit_lwr, ymax = probit_upr),
                        fill = input$probit_ribbon_color, alpha = 0.3) +
            geom_line(data = pred_df, aes(x = mean_conc, y = Probit),
                      color = input$probit_line_color, linewidth = 1.2) +
            geom_vline(xintercept = lod_probit, linetype = "dotted",
                       color = input$probit_line_color, linewidth = 1) +
            annotate("text", x = lod_probit*1.05, y = 0.3,
                     label = paste0("LoD_pro = ", round(lod_probit, 2)),
                     color = input$probit_line_color, angle = 0, hjust = 0)
        }
        
        # Add Weibull specifics
        if (input$method %in% c("Both", "Weibull") && !is.na(lod_weibull)) {
          p <- p + 
            geom_ribbon(data = pred_df, aes(x = mean_conc, ymin = weibull_lwr, ymax = weibull_upr),
                        fill = input$weibull_ribbon_color, alpha = 0.3) +
            geom_line(data = pred_df, aes(x = mean_conc, y = Weibull),
                      color = input$weibull_line_color, linewidth = 1.2) +
            geom_vline(xintercept = lod_weibull, linetype = "dotted",
                       color = input$weibull_line_color, linewidth = 1) +
            annotate("text", x = lod_weibull*1.05, y = 0.2,
                     label = paste0("LoD_wei = ", round(lod_weibull, 2)),
                     color = input$weibull_line_color, angle = 0, hjust = 0)
        }
        
        p <- p + scale_color_manual(name = "Model", values = c("Probit" = input$probit_line_color, "Weibull" = input$weibull_line_color))
        
        plot_list[[tgt]] <- p
      }
    })
    
    # Return everything needed for the UI and Downloads
    list(
      table = bind_rows(results_table),
      plots = plot_list,
      fit_recommendation = if(length(fit_texts) > 0) paste(fit_texts, collapse = "\n") else "Run 'Both' methods to see lack-of-fit comparison."
    )
  })
  
  # ================= OUTPUT RENDERERS =================
  output$docs_content <- renderUI({
    HTML(markdown::markdownToHTML(
      text = docs_text,
      fragment.only = TRUE
    ))
  })
  
  output$fit_recommendation <- renderText({
    req(analysis_results())
    analysis_results()$fit_recommendation
  })
  
  output$lod_table <- renderTable({
    req(analysis_results())
    analysis_results()$table
  })
  
  final_dashboard <- reactive({
    req(analysis_results())
    plots <- analysis_results()$plots
    
    # Use patchwork to wrap all generated plots
    wrapped <- wrap_plots(plots, ncol = min(3, length(plots)))
    
    # Extract shared legend to put on the side if "Both" is selected
    if (input$method == "Both") {
      dummy_plot <- ggplot(data.frame(x=1, y=1, m=c("Probit", "Weibull")), aes(x,y, color=m)) + 
        geom_line(linewidth=1.2) + 
        scale_color_manual(name="Model", values=c("Probit"=input$probit_line_color, "Weibull"=input$weibull_line_color)) +
        theme_minimal(base_size = 14)
      legend <- get_legend(dummy_plot)
      wrapped <- (wrapped | legend) + plot_layout(widths = c(8, 1))
    }
    
    wrapped + plot_annotation(title = "Limit of Detection", theme = theme(plot.title = element_text(size = 18, face = "bold")))
  })
  
  output$lod_plot <- renderPlot({
    final_dashboard()
  })
  
  # ================= DOWNLOAD HANDLERS =================
  output$download_plot <- downloadHandler(
    filename = function() { paste0("LoD_Dashboard_", Sys.Date(), ".png") },
    content = function(file) {
      ggsave(file, plot = final_dashboard(), 
             width = input$plot_width, 
             height = input$plot_height, 
             dpi = input$plot_dpi)
    }
  )
  
  output$download_csv <- downloadHandler(
    filename = function() { paste0("LoD_Results_", Sys.Date(), ".csv") },
    content = function(file) {
      write_csv(analysis_results()$table, file)
    }
  )
  # ================= DOWNLOAD HANDLERS =================
  
  output$download_plot <- downloadHandler(
    filename = function() { paste0("LoD_Dashboard_", Sys.Date(), ".png") },
    content = function(file) {
      ggsave(file, plot = final_dashboard(), 
             width = input$plot_width, 
             height = input$plot_height, 
             dpi = input$plot_dpi)
    }
  )
  
  output$download_csv <- downloadHandler(
    filename = function() { paste0("LoD_Results_", Sys.Date(), ".csv") },
    content = function(file) {
      write_csv(analysis_results()$table, file)
    }
  )
  
  # Replicate Summary Download Handler
  output$download_replicate_summary <- downloadHandler(
    filename = function() { paste0("Replicate_Summary_Quality_", Sys.Date(), ".csv") },
    content = function(file) {
      write_csv(replicate_summary_data(), file)
    }
  )
  
}


shinyApp(ui, server)