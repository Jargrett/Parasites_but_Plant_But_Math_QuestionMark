# shiny app - grouped by parameter

# libraries
library(shiny)
library(deSolve)
library(ggplot2)
library(tidyr)
library(patchwork)
library(bslib)

# model
model <- function(time, state, parms, new_h) {
  with(as.list(c(state, parms)), {
    # forb
    d_hf <- r_hf * hf * (1 - (hf / k_hf)) - 
      ((1 - psi_f * (hf / k_hf)) * (a_f * alpha_f) / 
         (1 + a_f * h_f * alpha_f)) -
      delta_lf * hf * hl - delta_gf * hg * hf 
    
    # legume
    d_hl <- r_hl * hl * (1 - (hl / k_hl)) -
      ((1 - psi_l * (hl / k_hl)) * (a_l * alpha_l) / 
         (1 + a_l * h_l * alpha_l)) -
      delta_fl * hf * hl - delta_gl * hg * hl 
    
    # grass
    d_hg <- r_hg * hg * (1 - (hg / k_hg)) -
      ((1 - psi_g * (hg / k_hg)) * (a_g * alpha_g) / 
         (1 + a_g * h_g * alpha_g)) -
      delta_fg * hf * hg - delta_lg * hg * hl
    
    # parasite
    d_p <- r_p * p * (1 - (p / k_p)) + 
      beta_f * ((a_f * alpha_f^2) / (1 + a_f * h_f * alpha_f)) + 
      beta_l * ((a_l * alpha_l) / (1 + a_l * h_l * alpha_l)) +
      beta_g * ((a_g * alpha_g) / (1 + a_g * h_g * alpha_g)) -
      (delta_fp * hf + delta_gp * hg + delta_lp * hl)
    
    # allocation strategy
    alpha_val <- switch(
      new_h,
      "bm" = {
        out <- c(hf, hl, hg) / sum(hf, hl, hg, na.rm = TRUE)
        ifelse(is.na(out), 0, out)
      },
      "unif" = {
        c(1/3, 1/3, 1/3)
      },
      "pc" = {
        vec <- c(
          ifelse(hf > 0, beta_f * ((a_f * (alpha_f + gamma * p)) / (1 + a_f * h_f * (alpha_f + gamma * p)) - (a_f * alpha_f) / (1 + a_f * h_f * alpha_f)), 0),
          ifelse(hl > 0, beta_l * ((a_l * (alpha_l + gamma * p)) / (1 + a_l * h_l * (alpha_l + gamma * p)) - (a_l * alpha_l) / (1 + a_l * h_l * alpha_l)), 0),
          ifelse(hg > 0, beta_g * ((a_g * (alpha_g + gamma * p)) / (1 + a_g * h_g * (alpha_g + gamma * p)) - (a_g * alpha_g) / (1 + a_g * h_g * alpha_g)), 0)
        )
        s_vec <- sum(vec, na.rm = TRUE)
        if (s_vec == 0 || is.na(s_vec)) c(0, 0, 0) else vec / s_vec
      }
    ) 
    
    d_alpha_f <- gamma * alpha_val[1] * p 
    d_alpha_l <- gamma * alpha_val[2] * p 
    d_alpha_g <- gamma * alpha_val[3] * p 
    
    return(list(c(d_hf, d_hl, d_hg, d_p, d_alpha_f, d_alpha_l, d_alpha_g)))
  })
}

# plotting function
plotting <- function(results, lab) {
  plot_df <- as.data.frame(results)
  
  long <- pivot_longer(
    plot_df,
    cols = c("hf", "hl", "hg", "p"),
    names_to = "variable",
    values_to = "value"
  )
  
  plot <- ggplot(long, aes(x = time, y = value, colour = variable)) +
    geom_line(linewidth = 1.2) +
    labs(
      title = paste0("host and parasite biomass - strategy: ", lab),
      x = "time",
      y = "biomass"
    ) +
    theme_classic() +
    scale_color_manual(values = c(
      hf = "#ff8000",
      hl = "#2c85b2",
      hg = "#6b990f",
      p  = "red"
    ))
  
  alpha_long <- pivot_longer(
    plot_df,
    cols = c("alpha_f", "alpha_l", "alpha_g"),
    names_to = "variable",
    values_to = "value"
  )
  
  alpha_plot <- ggplot(alpha_long, aes(x = time, y = value, colour = variable)) +
    geom_line(linewidth = 1.2) +
    labs(
      title = "haustoria",
      x = "time",
      y = "count / level"
    ) +
    theme_classic() +
    scale_color_manual(values = c(
      alpha_f = "#ffbf80",
      alpha_l = "#7ec3e5",
      alpha_g = "#a3cc51"
    ))
  
  plot / alpha_plot
}

# ui
ui <- page_sidebar(
  title = "hemiparasite model",
  sidebar = sidebar(
    width = 380,
    selectInput("new_h", "haustoria allocation strategy (new_h):", 
                choices = c("uniform (unif)" = "unif", 
                            "biomass proportional (bm)" = "bm", 
                            "parasite choice (pc)" = "pc")),
    numericInput("t_max", "simulation time steps:", value = 100, min = 10, max = 1000, step = 10),
    
    accordion(
      accordion_panel("initial state",
                      numericInput("init_hf", "initial forb biomass (hf):", value = 50, min = 0, max = 500),
                      numericInput("init_hl", "initial legume biomass (hl):", value = 50, min = 0, max = 500),
                      numericInput("init_hg", "initial grass biomass (hg):", value = 50, min = 0, max = 500),
                      numericInput("init_p",  "initial parasite biomass (p):", value = 10, min = 0, max = 200),
                      numericInput("init_alpha_f", "initial alpha_f:", value = 0, min = 0, max = 10, step = 0.1),
                      numericInput("init_alpha_l", "initial alpha_l:", value = 0, min = 0, max = 10, step = 0.1),
                      numericInput("init_alpha_g", "initial alpha_g:", value = 0, min = 0, max = 10, step = 0.1)
      ),
      accordion_panel("forb parameters",
                      sliderInput("r_hf", "growth (r_hf):", min = 0.01, max = 2.0, value = 0.40, step = 0.05),
                      sliderInput("k_hf", "carrying capacity (k_hf):", min = 10, max = 200, value = 50, step = 5),
                      sliderInput("psi_f", "defense (psi_f):", min = 0, max = 1, value = 0.5, step = 0.05),
                      sliderInput("beta_f", "parasite efficiency (beta_f):", min = 0, max = 1, value = 0.2, step = 0.05),
                      sliderInput("h_f", "handling time (h_f):", min = 0.01, max = 2.0, value = 0.5, step = 0.05),
                      sliderInput("a_f", "parasitism rate (a_f):", min = 0.001, max = 0.5, value = 0.005, step = 0.005)
      ),
      accordion_panel("legume parameters",
                      sliderInput("r_hl", "growth (r_hl):", min = 0.01, max = 2.0, value = 0.50, step = 0.05),
                      sliderInput("k_hl", "carrying capacity (k_hl):", min = 10, max = 200, value = 50, step = 5),
                      sliderInput("psi_l", "defense (psi_l):", min = 0, max = 1, value = 0.4, step = 0.05),
                      sliderInput("beta_l", "parasite efficiency (beta_l):", min = 0, max = 1, value = 0.3, step = 0.05),
                      sliderInput("h_l", "handling time (h_l):", min = 0.01, max = 2.0, value = 0.5, step = 0.05),
                      sliderInput("a_l", "parasitism rate (a_l):", min = 0.001, max = 0.5, value = 0.005, step = 0.005)
      ),
      accordion_panel("grass parameters",
                      sliderInput("r_hg", "growth (r_hg):", min = 0.01, max = 2.0, value = 0.60, step = 0.05),
                      sliderInput("k_hg", "carrying capacity (k_hg):", min = 10, max = 200, value = 65, step = 5),
                      sliderInput("psi_g", "defense (psi_g):", min = 0, max = 1, value = 0.2, step = 0.05),
                      sliderInput("beta_g", "parasite efficiency (beta_g):", min = 0, max = 1, value = 0.1, step = 0.05),
                      sliderInput("h_g", "handling time (h_g):", min = 0.01, max = 2.0, value = 0.5, step = 0.05),
                      sliderInput("a_g", "parasitism rate (a_g):", min = 0.001, max = 0.5, value = 0.090, step = 0.005)
      ),
      accordion_panel("parasite & global parameters",
                      sliderInput("r_p", "parasite growth (r_p):", min = 0.01, max = 2.0, value = 0.30, step = 0.05),
                      sliderInput("k_p", "parasite carrying capacity (k_p):", min = 1, max = 100, value = 10, step = 1),
                      sliderInput("gamma", "haustoria creation:", min = 0.1, max = 5, value = 1.0, step = 0.1),
                      sliderInput("l_p", "l_p:", min = 0, max = 0.5, value = 0.06, step = 0.01)
      ),
      accordion_panel("interspecific competition",
                      sliderInput("delta_lf", "delta_lf:", min = 0, max = 0.01, value = 0.0008, step = 0.0001),
                      sliderInput("delta_fl", "delta_fl:", min = 0, max = 0.01, value = 0.0006, step = 0.0001),
                      sliderInput("delta_gf", "delta_gf:", min = 0, max = 0.01, value = 0.0004, step = 0.0001),
                      sliderInput("delta_gl", "delta_gl:", min = 0, max = 0.01, value = 0.0004, step = 0.0001),
                      sliderInput("delta_fg", "delta_fg:", min = 0, max = 0.01, value = 0.0009, step = 0.0001),
                      sliderInput("delta_lg", "delta_lg:", min = 0, max = 0.01, value = 0.0008, step = 0.0001),
                      sliderInput("delta_fp", "delta_fp:", min = 0, max = 0.01, value = 0.0009, step = 0.0001),
                      sliderInput("delta_lp", "delta_lp:", min = 0, max = 0.01, value = 0.0008, step = 0.0001),
                      sliderInput("delta_gp", "delta_gp:", min = 0, max = 0.01, value = 0.0009, step = 0.0001)
      )
    )
  ),
  card(
    card_header("model simulation outputs"),
    plotOutput("model_plots", height = "600px"),
    card_footer(
      downloadButton("download_results", "download simulation results (csv)")
    )
  )
)

# server
server <- function(input, output, session) {
  init_states <- reactive({
    c(
      hf = input$init_hf, 
      hl = input$init_hl, 
      hg = input$init_hg, 
      p  = input$init_p, 
      alpha_f = input$init_alpha_f, 
      alpha_l = input$init_alpha_l, 
      alpha_g = input$init_alpha_g
    )
  })
  
  parameters <- reactive({
    c(
      r_hf = input$r_hf, k_hf = input$k_hf,
      r_hl = input$r_hl, k_hl = input$k_hl,
      r_hg = input$r_hg, k_hg = input$k_hg,
      
      delta_lf = input$delta_lf, delta_fl = input$delta_fl,
      delta_gf = input$delta_gf, delta_gl = input$delta_gl,
      delta_fg = input$delta_fg, delta_lg = input$delta_lg,
      
      delta_fp = input$delta_fp, delta_lp = input$delta_lp, delta_gp = input$delta_gp,
      
      psi_f = input$psi_f, psi_l = input$psi_l, psi_g = input$psi_g,
      
      beta_g = input$beta_g, beta_f = input$beta_f, beta_l = input$beta_l,
      r_p = input$r_p, k_p = input$k_p, gamma = input$gamma, l_p = input$l_p,
      
      h_f = input$h_f, h_g = input$h_g, h_l = input$h_l,
      a_f = input$a_f, a_g = input$a_g, a_l = input$a_l
    )
  })
  
  sim_data <- reactive({
    times <- seq(0, input$t_max, by = 1)
    
    res <- ode(
      y = init_states(),
      times = times,
      func = model,
      parms = parameters(),
      new_h = input$new_h
    )
    
    as.data.frame(res)
  })
  
  output$model_plots <- renderPlot({
    plotting(sim_data(), input$new_h)
  })
  
  output$download_results <- downloadHandler(
    filename = function() {
      paste0("model_results_", input$new_h, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      data <- sim_data()
      params_vec <- parameters()
      inits_vec <- init_states()
      
      meta <- c(
        paste0("# strategy: ", input$new_h),
        "# initial densities:",
        paste0("#   ", names(inits_vec), ": ", inits_vec),
        "# model parameters:",
        paste0("#   ", names(params_vec), ": ", params_vec),
        "#"
      )
      
      writeLines(meta, file)
      write.table(data, file, append = TRUE, sep = ",", row.names = FALSE, col.names = TRUE)
    }
  )
}

shinyApp(ui = ui, server = server)
