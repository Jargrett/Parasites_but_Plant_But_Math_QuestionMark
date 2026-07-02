library(deSolve)
library(ggplot2)
library(tidyr)

# define model
hemiparasite_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    H <- state[1:3]
    P <- state[4]
    
    # competition coefficients, L-V
    competition_effects <- c(
      H[1] + alpha12*H[2] + alpha13*H[3], # legume surrounding pressure
      alpha21*H[1] + H[2] + alpha23*H[3], # grass surrounding pressure
      alpha31*H[1] + alpha32*H[2] + H[3]  # forb surrounding pressure
    )
    
    # shading effects - aboveground
    total_host_biomass <- sum(H)
    autotrophic_growth <- max(0, eps * P * (1 - total_host_biomass / K_canopy))
    
    # resource stealing - belowground
    gross_intake <- sum((1 - R) * q * H * beta)
    heterotrophic_growth <- (M_max * gross_intake) / (K_sat + gross_intake) * P
    
    # ODEs for each host (H1,2,3 and the parasite)
    dH1 <- r[1] * H[1] * (1 - competition_effects[1] / K[1]) - 
      ((1 - R[1]) * q[1] * H[1] * P)
    dH2 <- r[2] * H[2] * (1 - competition_effects[2] / K[2]) - 
      ((1 - R[2]) * q[2] * H[2] * P)
    dH3 <- r[3] * H[3] * (1 - competition_effects[3] / K[3]) - 
      ((1 - R[3]) * q[3] * H[3] * P)
    
    dP  <- autotrophic_growth + heterotrophic_growth - (m * P)
    
    return(list(c(dH1, dH2, dH3, dP)))
  })
}

# define parameters ----
# order is legume, dominant grass, forb with defenses
# host 1: legume w/ moderate growth rate, high value host
# host 2: dominant grass w/ fast growth rate, more moderate value host
# host 3: forb with defenses against haustoria, lowest growth rate, 
#         lower value host

base_params <- list(
  r        = c(0.45, 0.75, 0.35),  # host growth rate
  K        = c(150, 250, 120),     # host carrying capacity
  R        = c(0.0, 0.0, 0.80),    # resistance to haustoria attaching
  q        = c(0.04, 0.05, 0.03),  # attachment affinities ? - modify based on 
                                   # haustoria already attached, parasite 
                                   # preferences, no preference (random) etc. 
  beta     = c(0.60, 0.20, 0.10),  # host quality
  eps      = 0.04,                 # parasite autotrophy  what does it mean? 
                                   # 20-60%
  K_canopy = 450,                  # carrying capacity of total canopy (light)
  M_max    = 0.8,                  # max heterotrophic efficiency
  K_sat    = 25,                   # multi-host saturation constant
  m        = 0.08,                 # metabolic maintenance cost
  
  # incorporate more timing, lags, etc. 1-2 wks from germinated to connected 
  # depending on biomass, attachment
  
  # competition matrix
  alpha12 = 1.5, alpha13 = 0.5, # grass (2) strongly suppresses legume (1)
  alpha21 = 0.4, alpha23 = 0.4, # grass ignores others
  alpha31 = 0.6, alpha32 = 1.8  # grass (2) completely suppresses forb (3)
)

times <- seq(0, 150, by = 0.2)

# biology points:
# haustorial count and rates of stealing

# scenario 1: parasite invasion into an even community ----
out_control_A <- as.data.frame(ode(y = c(H=c(40,40,40), P=0.0), times=times, func=hemiparasite_model, parms=base_params))
out_paras_A   <- as.data.frame(ode(y = c(H=c(40,40,40), P=1.0), times=times, func=hemiparasite_model, parms=base_params))

total_biomass_data_A <- data.frame(
  Time = times,
  Control = rowSums(out_control_A[, 2:4]),
  With_Parasite = rowSums(out_paras_A[, 2:4]) + out_paras_A[, 5]
) %>% pivot_longer(-Time, names_to = "Condition", values_to = "Total_Biomass")

ggplot(total_biomass_data_A, aes(x=Time, y=Total_Biomass, color=Condition)) +
  geom_line(linewidth=1.2) + theme_minimal() +
  labs(y = "Total Ecosystem Biomass (Hosts + Parasite)")
# parasite reduces community biomass, tax on whole community

# scenario 2: parasite with legume compensation ----
params_B <- base_params
params_B$alpha12 <- 1.1 

out_control_B <- as.data.frame(ode(y = c(H=c(30,70,20), P=0.0), times=times, func=hemiparasite_model, parms=params_B))
out_paras_B   <- as.data.frame(ode(y = c(H=c(30,70,20), P=1.0), times=times, func=hemiparasite_model, parms=params_B))

total_biomass_data_B <- data.frame(
  Time = times,
  Control = rowSums(out_control_B[, 2:4]),
  With_Parasite = rowSums(out_paras_B[, 2:4]) + out_paras_B[, 5]
) %>% pivot_longer(-Time, names_to = "Condition", values_to = "Total_Biomass")

ggplot(total_biomass_data_B, aes(x=Time, y=Total_Biomass, color=Condition)) +
  geom_line(size=1.2) + theme_minimal() + scale_color_manual(values=c("Control"="black", "With_Parasite"="purple")) +
  labs(y = "Total Ecosystem Biomass")

# scenario 3 ----
params_C <- base_params
params_C$alpha12 <- 2.2
params_C$alpha32 <- 2.5
params_C$K[1] <- 220 

# run
out_control_C <- as.data.frame(ode(y = c(H=c(20,90,15), P=0.0), times=times, func=hemiparasite_model, parms=params_C))
out_paras_C   <- as.data.frame(ode(y = c(H=c(20,90,15), P=1.0), times=times, func=hemiparasite_model, parms=params_C))

total_biomass_data_C <- data.frame(
  Time = times,
  Control = rowSums(out_control_C[, 2:4]),
  With_Parasite = rowSums(out_paras_C[, 2:4]) + out_paras_C[, 5]
) %>% pivot_longer(-Time, names_to = "Condition", values_to = "Total_Biomass")

plot_biomass_over_time <- ggplot(total_biomass_data_C, aes(x=Time, y=Total_Biomass, color=Condition)) +
  geom_line(size=1.5) + theme_bw() + scale_color_manual(values=c("Control"="red", "With_Parasite"="darkgreen"))

colnames(out_control_C) <- c("Time", "Legume", "Grass", "Forb", "Parasite")
colnames(out_paras_C)   <- c("Time", "Legume", "Grass", "Forb", "Parasite")

long_comp_control <- pivot_longer(out_control_C, cols=2:5, names_to="Species", values_to="Biomass") %>% dplyr::mutate(System="Control")
long_comp_paras   <- pivot_longer(out_paras_C, cols=2:5, names_to="Species", values_to="Biomass") %>% dplyr::mutate(System="With Parasite")

composition_data <- rbind(long_comp_control, long_comp_paras)

ggplot(composition_data, aes(x=Time, y=Biomass, 
                             color=Species, linetype=Species)) +
  geom_line(size=1.2) + 
  facet_wrap(~System) + 
  theme_classic() +
  scale_color_manual(values = c("Legume" = "orange3", "Grass" = "forestgreen", "Forb" = "purple4", "Parasite" = "black")) +
  scale_linetype_manual(values = c("Legume" = 1, "Grass" = 1, "Forb" = 1, "Parasite" = 2)) 

# monocultures???
# does the parasite have a growth rate?