# libraries
library(shiny)
library(deSolve)
library(ggplot2)
library(tidyr)
library(patchwork)
library(bslib)

# model ---- 
model <- function(time, state, parms, host_type) {
  with(as.list(c(state, parms)), {
    
    # host parameters
    r_h <- parms[[paste0("r_h_", host_type)]]
    K_h <- parms[[paste0("K_h_", host_type)]]
    delta_hp <- parms[[paste0("delta_hp_", host_type)]]
    delta_ph <- parms[[paste0("delta_ph_", host_type)]]
    alpha <- parms[[paste0("alpha_", host_type)]]
    
    # rates & functional responses
    frac <- (Z / H) / (c_impact + (Z / H))
    
    # ODEs (updated to match formulation with beta)
    d_H <- (r_h - r_h * beta * frac) * H * (1 - (H + delta_hp * P) / K_h)
    d_P <- r_p * P * (1 - (P + delta_ph * H) / K_p) + alpha * beta * frac
    d_Z <- gamma * P - mu * Z
    
    # prevent negative values
    if (H <= 0) {
      d_H <- 0
      H <- 0
    }
    if (P <= 0 && d_P < 0) d_P <- 0
    if (Z <= 0 && d_Z < 0) d_Z <- 0
    
    return(list(c(d_H, d_P, d_Z)))
  })
}

# ui ----
ui <- page_sidebar(
  title = "Hemiparasite Model (Pairwise Comparisons)",
  sidebar = sidebar(
    width = 400,
    fileInput("upload_csv", "Upload parameters from csv:", accept = ".csv"),
    numericInput("t_max", "Simulation time steps:", value = 100, min = 10, 
                 max = 500, step = 10),
    
    accordion(
      accordion_panel("initial state",
                      numericInput("init_H", "initial host biomass (H):", 
                                   value = 10, min = 1, step = 1),
                      numericInput("init_P", "initial parasite biomass (P):",
                                   value = 2, min = 0, step = 1),
                      numericInput("init_Z", "initial haustoria (Z):", 
                                   value = 0, min = 0, step = 1)
      ),
      accordion_panel("host growth rates (r_h)",
                      sliderInput("r_h_forb", "forb growth rate:", min = 0.001,
                                  max = 0.2, value = 0.08, step = 0.01),
                      sliderInput("r_h_grass", "grass growth rate:", min = 0.001,
                                  max = 0.2, value = 0.15, step = 0.01),
                      sliderInput("r_h_legume", "legume growth rate:", 
                                  min = 0.001, max = 0.2, value = 0.04, 
                                  step = 0.01)
      ),
      accordion_panel("host carrying capacity (K_h)",
                      sliderInput("K_h_forb", "forb carrying capacity:", 
                                  min = 10, max = 200, value = 40, step = 5),
                      sliderInput("K_h_grass", "grass carrying capacity:", 
                                  min = 10, max = 200, value = 60, step = 5),
                      sliderInput("K_h_legume", "legume carrying capacity:", 
                                  min = 10, max = 200, value = 35, step = 5)
      ),
      
      accordion_panel("parasite efficiency (alpha)",
                      sliderInput("alpha_forb", "forb attachment rate:", 
                                  min = 0.05, max = 1.0, value = 0.30, 
                                  step = 0.05),
                      sliderInput("alpha_grass", "grass attachment rate:", 
                                  min = 0.05, max = 1.0, value = 0.15, 
                                  step = 0.05),
                      sliderInput("alpha_legume", "legume attachment rate:", 
                                  min = 0.05, max = 1.0, value = 0.50, 
                                  step = 0.05)
      ),
      accordion_panel("competition (delta_hp & delta_ph)",
                      sliderInput("delta_hp_forb", 
                                  "forb on parasite (delta_hp):", 
                                  min = 0.0, max = 1, value = 0.001, 
                                  step = 0.01),
                      sliderInput("delta_ph_forb", 
                                  "parasite on forb (delta_ph):", 
                                  min = 0.0, max = 1, value = 0.001, 
                                  step = 0.01),
                      sliderInput("delta_hp_grass", 
                                  "grass on parasite (delta_hp):", 
                                  min = 0.0, max = 1, value = 0.001, 
                                  step = 0.01),
                      sliderInput("delta_ph_grass", 
                                  "parasite on grass (delta_ph):", 
                                  min = 0.0, max = 1, value = 0.001, 
                                  step = 0.01),
                      sliderInput("delta_hp_legume", 
                                  "legume on parasite (delta_hp):", 
                                  min = 0.0, max = 1, value = 0.001, 
                                  step = 0.01),
                      sliderInput("delta_ph_legume", 
                                  "parasite on legume (delta_ph):", 
                                  min = 0.0, max = 1, value = 0.001, 
                                  step = 0.01)),
      accordion_panel("parasite & haustoria parameters",
                      sliderInput("r_p", "parasite base growth (r_p):", 
                                  min = 0.01, max = 1.0, value = 0.25, 
                                  step = 0.05),
                      sliderInput("K_p", "parasite carrying capacity (K_p):",
                                  min = 5, max = 100, value = 20, step = 5),
                      sliderInput("gamma", "haustoria creation (gamma):", 
                                  min = 0.1, max = 5.0, value = 1.2, 
                                  step = 0.1),
                      sliderInput("c_impact", "host growth impact (c_impact):",
                                  min = 0, max = 2, value = 1, 
                                  step = 0.01),
                      sliderInput("beta", "resource conversion efficiency (beta):",
                                  min = 0.01, max = 2.0, value = 1.0,
                                  step = 0.05),
                      sliderInput("mu", "haustoria death rate (mu):", 
                                  min = 0.01, max = 0.5, value = 0.1, 
                                  step = 0.01)
      ))
  ),
  navset_card_tab(
    id = "tabs",
    
    # tab 1: about
    nav_panel(title = "About",
              div(style = "padding: 15px;",
                  h3("Project Questions & Main Objectives"),
                  p("Broad Objective Placeholder"),
                  h4("Main Questions:"),
                  tags$ol(
                    tags$li("Question 1 placeholder"),
                    tags$li("Question 2 placeholder"),
                    tags$li("Question 3 placeholder")))),
    
    # tab 2: model structure & parameters
    nav_panel(title = "Model Structure & Parameters",
              div(style = "padding: 15px;",
                  h3("Model Formulation"),
                  withMathJax(
                    p("$$\\frac{dH}{dt} = \\left(r_h - r_h \\beta \\left(\\frac{\\frac{Z}{H}}{\\text{c_impact} + \\frac{Z}{H}}\\right)\\right) H \\left(1 - \\frac{H + \\delta_{hp} P}{K_h}\\right)$$"),
                    p("$$\\frac{dP}{dt} = r_p P \\left(1 - \\frac{P + \\delta_{ph} H}{K_p}\\right) + \\alpha \\beta \\left(\\frac{\\frac{Z}{H}}{\\text{c_impact} + \\frac{Z}{H}}\\right)$$"),
                    p("$$\\frac{dZ}{dt} = \\gamma P - \\mu Z$$")
                  ),
                  hr(),
                  tableOutput("param_table"))),
    
    # tab 3: simulation outputs
    nav_panel(title = "Simulation Outputs",
              div(style = "padding: 15px;",
                  card(
                    card_header("Model Simulation Outputs"),
                    plotOutput("six_panel_plot", height = "750px"),
                    card_footer(
                      downloadButton("download_results", 
                                     "download simulation results (csv)")))
              ))
  )
)

# server ----
server <- function(input, output, session) {
  
  # reupload parameters from csv
  observeEvent(input$upload_csv, {
    req(input$upload_csv)
    lines <- readLines(input$upload_csv$datapath)
    comment_lines <- lines[grep("^#", lines)]
    
    kv_lines <- comment_lines[grep(":", comment_lines)]
    for (line in kv_lines) {
      clean_line <- sub("^#\\s*", "", line)
      parts <- strsplit(clean_line, ":")[[1]]
      if (length(parts) == 2) {
        param_name <- trimws(parts[1])
        param_val <- as.numeric(trimws(parts[2]))
        if (!is.na(param_val)) {
          input_id <- if (param_name %in% c("H", "P", "Z")) {
            paste0("init_", param_name)
          } else {
            param_name
          }
          if (input_id %in% c("init_H", "init_P", "init_Z", "t_max")) {
            updateNumericInput(session, input_id, value = param_val)
          } else {
            updateSliderInput(session, input_id, value = param_val)
          }
        }
      }
    }
  })
  
  parameters <- reactive({
    list(
      r_p = input$r_p, K_p = input$K_p, gamma = input$gamma,
      c_impact = input$c_impact, beta = input$beta, mu = input$mu,
      
      r_h_forb = input$r_h_forb, K_h_forb = input$K_h_forb, 
      alpha_forb = input$alpha_forb,
      delta_hp_forb = input$delta_hp_forb, delta_ph_forb = input$delta_ph_forb,
      
      r_h_grass = input$r_h_grass, K_h_grass = input$K_h_grass, 
      alpha_grass = input$alpha_grass,
      delta_hp_grass = input$delta_hp_grass, delta_ph_grass = input$delta_ph_grass,
      
      r_h_legume = input$r_h_legume, K_h_legume = input$K_h_legume, 
      alpha_legume = input$alpha_legume,
      delta_hp_legume = input$delta_hp_legume, delta_ph_legume = input$delta_ph_legume)
  })
  
  sim_results <- reactive({
    times <- seq(0, input$t_max, by = 1)
    parms <- parameters()
    inits <- c(H = input$init_H, P = input$init_P, Z = input$init_Z)
    
    hosts <- c("forb", "grass", "legume")
    res_list <- list()
    
    for (h_type in hosts) {
      out <- tryCatch({
        withWarnings <- function(expr) {
          myWarnings <- NULL
          wHandler <- function(w) {
            myWarnings <<- c(myWarnings, w$message)
            invokeRestart("muffleWarning")
          }
          val <- withCallingHandlers(expr, warning = wHandler)
          list(value = val, warnings = myWarnings)
        }
        
        res <- withWarnings(ode(y = inits, times = times, func = model, 
                                parms = parms, host_type = h_type))
        
        if (!is.null(res$warnings)) {
          showNotification(
            paste0("lsoda warning (", h_type, "): ", res$warnings[1]),
            type = "warning", duration = 6)
        }
        res$value
      }, error = function(e) {
        showNotification(paste0("ode error (", h_type, "): ", e$message), 
                         type = "error")
        return(matrix(c(times, rep(NA, length(times) * 3)), ncol = 4))
      })
      
      df <- as.data.frame(out)
      colnames(df)[1:4] <- c("time", "H", "P", "Z")
      df$host_type <- h_type
      res_list[[h_type]] <- df
    }
    
    do.call(rbind, res_list)
  })
  
  # parameter metadata table
  output$param_table <- renderTable({
    data.frame(
      Symbol = c("H", "P", "Z", "r_h", "K_h", "r_p", "K_p", "c_impact", 
                 "beta", "alpha", 
                 "gamma", "mu", "delta_hp / delta_ph"),
      Parameter = c(
        "Host Biomass", "Parasite Biomass", "Haustoria Density", 
        "Host Growth Rate", "Host Carrying Capacity", "Parasite Growth Rate", 
        "Parasite Carrying Capacity", "Parasite Impact Coefficient", 
        "Resource Conversion Efficiency", "Attachment Rate Constant", 
        "Haustoria Creation Rate", "Haustoria Loss Rate", 
        "Competition Coefficients"
      ),
      Description = c(
        "Total living biomass of host plant species",
        "Total living biomass of hemiparasite species",
        "Total number of functional haustorial connections on host",
        "Baseline intrinsic growth rate of host plant without parasites or competitors",
        "Maximum host biomass supported by abiotic environmental resources",
        "Baseline intrinsic growth rate of parasite from photosynthesis",
        "Maximum parasite biomass sustainable by baseline structures",
        "Per-haustorium rate constant altering host growth efficiency reduction",
        "Maximum biomass conversion rate from stolen host resources",
        "Baseline probability coefficient of haustorium attachment",
        "Production rate of new haustoria per unit of parasite biomass",
        "Combined rate of haustorial senescence, decay, or host ejection",
        "Interspecific competitive suppression factors between host and parasite"
      )
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  # plotting function
  output$six_panel_plot <- renderPlot({
    df_all <- sim_results()
    
    host_types <- c("forb", "grass", "legume")
    host_colors <- c("forb" = "#e66101", 
                     "grass" = "#4dac26", 
                     "legume" = "#2b83ba")
    
    plots <- list()
    
    for (h_type in host_types) {
      sub_df <- df_all[df_all$host_type == h_type, ]
      
      p_biomass <- ggplot(sub_df, aes(x = time)) +
        geom_line(aes(y = H, color = "Host"), linewidth = 1.1, na.rm = TRUE) +
        geom_line(aes(y = P, color = "Parasite"), linewidth = 1.1, na.rm = TRUE) +
        scale_color_manual(values = c("Host" = host_colors[[h_type]], 
                                      "Parasite" = "#d7191c")) +
        coord_cartesian(xlim = c(0, input$t_max)) +
        labs(title = paste0(h_type), x = "time", y = "biomass", color = "species") +
        theme_classic() +
        theme(legend.position = "bottom")
      
      p_haustoria <- ggplot(sub_df, aes(x = time, y = Z)) +
        geom_line(color = "#7b3294", linewidth = 1.1, na.rm = TRUE) +
        coord_cartesian(xlim = c(0, input$t_max)) +
        labs(title = paste0(h_type), x = "time", y = "count") +
        theme_classic()
      
      plots[[paste0(h_type, "_bio")]] <- p_biomass
      plots[[paste0(h_type, "_hau")]] <- p_haustoria
    }
    
    (plots[["forb_bio"]] | plots[["grass_bio"]] | plots[["legume_bio"]]) /
      (plots[["forb_hau"]] | plots[["grass_hau"]] | plots[["legume_hau"]])
  })
  
  # download handler
  output$download_results <- downloadHandler(
    filename = function() {
      paste0("hemiparasite_model_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      data <- sim_results()
      params_vec <- unlist(parameters())
      inits_vec <- c(H = input$init_H, P = input$init_P, Z = input$init_Z)
      
      meta <- c(
        "# initial states:",
        paste0("#    ", names(inits_vec), ": ", inits_vec),
        # model parameters:
        paste0("#    ", names(params_vec), ": ", params_vec),
        "#"
      )
      
      writeLines(meta, file)
      write.table(data, file, append = TRUE, sep = ",", 
                  row.names = FALSE, col.names = TRUE)
    })
}

shinyApp(ui = ui, server = server)
