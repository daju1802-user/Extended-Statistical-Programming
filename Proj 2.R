#Group 6 :Xuyi Shi(S2796586), Chi Zhang(S2828433), Jiachen Guang(S2789777)
#Repo Link: https://github.com/daju1802-user/Extended-Statistical-Programming.git

#This code implements an extended SEIR epidemic model that incorporates social 
#network structure, used to investigate the impact of family relationships and 
#individual differences in sociability on the spread of infectious diseases. 
#The model distributes the population into households of random size and 
#constructs a contact network based on sociability parameters. It simulates three 
#infection pathways: intra-household transmission, transmission with fixed contacts, 
#and random mixing. By comparing four different scenarios (a complete model, 
#a model without social structure, a model with homogeneous sociability, and a baseline model), 
#we analyze how social structure influences the dynamics, scale, and peak of the epidemic.


#Contribution:
#Xuyi Shi(34%): Established vector h, sorted out details and modeling ideas, and 
#completed the "get.net" function part, providing ideas for the "nseir" function, 
#and doing final checks and tests.

#Chi Zhang(33%): Completed the part on the "nseir" function and checked the idea and 
#organization of the function. Discuss ideas and final model comparison with teammates,
#and added comments, and final checks.

#Jiachen Guang(33%): Completed the drawing part, discussed the problem of comparing models, 
#completed the code for this part, and provided subsequent comments through images, doing
#the final check.


#get.net(): Constructing social networks based on individual social parameters

#nseir(): Core simulation function, running the SEIR model including social structure

#plot.nseir(): Visualize the epidemic development curve

#compare.scenarios(): Comparing models under different conditions


##============================================================================================================
##------------------Create the net funtion between regular (non-household) contact people---------------------
## Note: 
# Construct a social network where connections form based on individual sociability parameters.
# The probability of connection between individuals i and j is proportional to beta[i] * beta[j].
# No connections are created between people in the same household.

get.net = function(beta, h, nc = 15) {
  # Creates regular contact network based on sociability parameters (OPTIMIZED)
  # beta: vector of sociability parameters for each person
  # h: household assignment vector
  # nc: average number of contacts per person
  # Returns: list where element i contains indices of i's regular contacts
  
  n = length(beta)
  beta_bar = mean(beta)
  
  # Vectorized approach using expand.grid (faster for moderate n)
  # Generate all unique pairs (i, j) where i < j
  pairs = which(lower.tri(matrix(0, n, n)), arr.ind = TRUE)
  i = pairs[, 1]
  j = pairs[, 2]
  
  # Filter out same-household pairs (vectorized)
  diff_household <- h[i] != h[j]
  i = i[diff_household]
  j = j[diff_household]
  
  # Calculate link probabilities (vectorized)
  probs = nc * beta[i] * beta[j] / (beta_bar^2 * (n - 1))
  
  # Determine which links to create (vectorized)
  link <- runif(length(probs)) < probs
  
  # Extract pairs that form links
  linked_i = i[link]
  linked_j = j[link]
  
  # Build contact list efficiently using split
  contacts = vector("list", n)
  
  # Add links in both directions
  if (length(linked_i) > 0) {
    # For each i, collect all j's
    contacts_i = split(linked_j, linked_i)
    # For each j, collect all i's
    contacts_j = split(linked_i, linked_j)
    
    # Merge into final contact list
    for (id in names(contacts_i)) {
      contacts[[as.integer(id)]] = contacts_i[[id]]
    }
    for (id in names(contacts_j)) {
      idx = as.integer(id)
      contacts[[idx]] = c(contacts[[idx]], contacts_j[[id]])
    }
  }
  
  return(contacts)
}
##============================================================================================================
##-----------------------------------Create a function of infections------------------------------------------
## Simulate the SEIR epidemic dynamics with three infection pathways:
## 1. Household transmission: infection probability alpha_h between household members
## 2. Network transmission: infection probability alpha_c between social contacts  
## 3. Random mixing: infection probability proportional to alpha_r and sociability parameters
## Individuals progress through states: Susceptible -> Exposed -> Infectious -> Recovered

nseir = function(beta, h, alink, alpha = c(.1, .01, .01), delta = .2, 
                  gamma = .4, nc = 15, nt = 100, pinf = .005) {
  # Simulates SEIR epidemic with household and network structure
  # beta: vector of sociability parameters for each person
  # h: household assignment vector
  # alink: list of regular contacts for each person (from get.net)
  # alpha: c(alpha_h, alpha_c, alpha_r) - infection probabilities for 
  #        household, contact network, and random mixing
  # delta: daily probability of recovery (I -> R)
  # gamma: daily probability of becoming infectious (E -> I)
  # nc: average number of contacts per person
  # nt: number of days to simulate
  # pinf: initial proportion in I state
  # Returns: list with S, E, I, R (daily counts) and t (days)
  
  n = length(beta)  # population size
  beta_bar = mean(beta)  # mean sociability
  
  # Extract alpha parameters
  alpha_h = alpha[1]  # household infection probability
  alpha_c = alpha[2]  # contact network infection probability
  alpha_r = alpha[3]  # random mixing infection probability
  
  # Initialize states: 1=S, 2=E, 3=I, 4=R
  state = rep(1, n)
  
  # Randomly select initial infected individuals
  initial_infected = sample(n, round(n * pinf))
  state[initial_infected] = 3
  
  # Storage for daily counts
  S = E = I = R = numeric(nt + 1)
  t = 0:nt                        #It starts from the initial state, so time should start from 0
  
  # Record initial state
  S[1] = sum(state == 1)
  E[1] = sum(state == 2)
  I[1] = sum(state == 3)
  R[1] = sum(state == 4)
  
  # Simulate epidemic dynamics
  for (day in 1:nt) {
    new_state = state
    # Find individuals in each state
    susceptible = which(state == 1) 
    infectious = which(state == 3)
    exposed = which(state == 2)
    
    
    # Process transitions E -> I (exposed become infectious)
    if (length(exposed) > 0) {
      become_infectious = exposed[runif(length(exposed)) < gamma]
      new_state[become_infectious] = 3
    }
    
    # Process transitions I -> R (infectious recover)
    if (length(infectious) > 0) {
      recover = infectious[runif(length(infectious)) < delta]
      new_state[recover] = 4
    }
    
    # Process infections S -> E (susceptible become exposed)
    if (length(infectious) > 0 && length(susceptible) > 0) {
      
      # Loop through each infectious person
      for (i in infectious) {
        
        # 1. Household infections
        if (alpha_h > 0) {
          # Find household members who are susceptible
          household_members = which(h == h[i] & state == 1)
          if (length(household_members) > 0) {
            # Each susceptible household member gets infected with probability alpha_h
            infected = household_members[runif(length(household_members)) < alpha_h]
            new_state[infected] = 2
          }
        }
        
        # 2. Regular contact network infections
        if (alpha_c > 0 && length(alink[[i]]) > 0) {
          # Find contacts who are susceptible
          susceptible_contacts = alink[[i]][state[alink[[i]]] == 1]
          if (length(susceptible_contacts) > 0) {
            # Each susceptible contact gets infected with probability alpha_c
            infected = susceptible_contacts[runif(length(susceptible_contacts)) < alpha_c] 
            new_state[infected] = 2
          }
        }
        
        # 3. Random mixing infections
        if (alpha_r > 0) {
          # Calculate infection probability for each susceptible person
          susceptible_now = which(state == 1)
          if (length(susceptible_now) > 0) {
            # Infection probability depends on both individuals' sociability
            prob = alpha_r * nc * beta[i] * beta[susceptible_now] / (beta_bar^2 * (n - 1))
            #Find all susceptible people who can be infected by the current infected person i
            infected = susceptible_now[runif(length(susceptible_now)) < prob]
            new_state[infected] = 2
          }
        }
      }
    }
    
    # Update state for next day
    state = new_state
    
    # Record daily counts
    S[day + 1] = sum(state == 1)
    E[day + 1] = sum(state == 2)
    I[day + 1] = sum(state == 3)
    R[day + 1] = sum(state == 4)
  }
  
  # Return results
  return(list(S = S, E = E, I = I, R = R, t = t,beta = beta))
}


##============================================================================================================
##-----------------------------Create the function to plot the dynamics states----------------------------------
plot.nseir <- function(nseir_results, main = "SEIR Epidemic Dynamics") {
  # More detailed plotting with individual panels and summary statistics
  # nseir_results: list output from nseir function
  # main: main title for the plot
  
  S <- nseir_results$S
  E <- nseir_results$E
  I <- nseir_results$I
  R <- nseir_results$R
  t <- nseir_results$t
  beta <- nseir_results$beta  # Extract beta values if available
  
  n <- S[1] + E[1] + I[1] + R[1]  # total population
  
  # Set up 2x3 plot layout to accommodate beta histogram
  old_par = par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
  on.exit(par(old_par))
  
  # Plot 1: All classes together
  plot(t, S, type = "l", col = "black", lwd = 2.5, main = main,
       xlab = "Day", ylab = "Population", ylim = c(0, n))
  lines(t, E, col = "blue", lwd = 2.5)
  lines(t, I, col = "red", lwd = 2.5)
  lines(t, R, col = "green", lwd = 2.5)
  legend("right", legend = c("S", "E", "I", "R"),
         col = c("black", "blue", "red", "green"), lwd = 2.5, bty = "n")
  grid(col = "gray80", lty = "dotted")
  
  # Plot 2: Susceptible only
  plot(t, S, type = "l", col = "black", lwd = 2.5, 
       main = "Susceptible", xlab = "Day", ylab = "Population")
  grid(col = "gray80", lty = "dotted")
  
  # Plot 3: Infectious only
  plot(t, I, type = "l", col = "red", lwd = 2.5,
       main = "Infectious", xlab = "Day", ylab = "Population")
  grid(col = "gray80", lty = "dotted")
  
  # Plot 4: Cumulative infections (E + I + R)
  cumulative_infected <- E + I + R
  plot(t, cumulative_infected, type = "l", col = "purple", lwd = 2.5,
       main = "Cumulative Infections(E + I + R)", xlab = "Day", ylab = "Population")
  grid(col = "gray80", lty = "dotted")
  
  # Plot 5: Beta distribution histogram
  if (!is.null(beta)) {
    hist(beta, col = "lightblue", main = "Distribution of Beta Values",
         xlab = "Beta (Sociability)", ylab = "Frequency", border = "white")
    abline(v = mean(beta), col = "red", lwd = 2, lty = 2)
    legend("topright", legend = paste("Mean =", round(mean(beta), 4)), 
           col = "red", lwd = 2, bty = "n")
    grid(col = "gray80", lty = "dotted")
  } else {
    # If no beta values available, show empty plot with message
    plot(0, 0, type = "n", xlab = "", ylab = "", main = "No Beta Data", axes = FALSE)
    text(0, 0, "Beta values not available", cex = 1.2)
  }
  
  # Plot 6: Epidemic summary statistics
  peak_infectious <- max(I)
  peak_day <- t[which.max(I)]
  final_size <- R[length(R)]
  attack_rate <- final_size / n * 100
  
  # Create summary plot
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), 
       main = "Epidemic Summary", xlab = "", ylab = "", axes = FALSE)
  
  text(0.5, 0.8, paste("Peak Infectious:", round(peak_infectious)), cex = 1.2)
  text(0.5, 0.6, paste("Peak Day:", peak_day), cex = 1.2)
  text(0.5, 0.4, paste("Final Size:", round(final_size)), cex = 1.2)
  text(0.5, 0.2, paste("Attack Rate:", round(attack_rate, 1), "%"), cex = 1.2)
}


#Compare the different scenarios with the function
compare.scenarios <- function(results_list, scenario_names = NULL, main_title = "", 
                          show_beta = TRUE) {
  # Compares multiple nseir results in panels side-by-side
  # by using plot.nseir for each scenario
  # results_list: list of nseir outputs
  # scenario_names: names for each scenario
  # main_title: overall title for the comparison
  # show_beta: if TRUE, includes beta distributions in the plots
  
  n_scenarios <- length(results_list)
  
  # Set default scenario names if not provided
  if (is.null(scenario_names)) {
    scenario_names <- paste("Scenario", 1:n_scenarios)
  }
  
  # Set up plot layout - if show_beta, each scenario needs 2 rows (SEIR + beta)
  # Otherwise just 1 row per scenario
  n_cols <- n_scenarios  
  n_rows <- ceiling(n_scenarios / n_cols)
  
  if (show_beta && !is.null(results_list[[1]]$beta)) {
    # Each scenario takes 2 rows: one for SEIR, one for beta
    n_plot_rows <- n_rows * 2
    n_plot_cols <- n_cols
  } else {
    n_plot_rows <- n_rows
    n_plot_cols <- n_cols
  }
  
  old_par <- par(mfrow = c(n_plot_rows, n_plot_cols), 
                 mar = c(5, 2.5, 2, 1),
                 oma = c(0, 0, 2, 0))
  on.exit(par(old_par)) # Restore the original par settings
  
  # Find global max for consistent y-axis across all scenarios
  max_pop <- 0
  for (res in results_list) {
    max_pop <- max(max_pop, max(res$S, res$E, res$I, res$R))
  }
  # Plot beta distributions for each scenario if requested
  if (show_beta && !is.null(results_list[[1]]$beta)) {
    for (i in 1:n_scenarios) {
      res <- results_list[[i]]
      hist(res$beta, col = "steelblue", border = "white",
           main = "", xlab = "β", ylab = "Frequency",
           cex.axis = 0.9, cex.lab = 0.9)
      grid(col = "gray80", lty = "dotted", lwd = 0.8)
    }
  }
  
  # Plot SEIR dynamics for each scenario
  for (i in 1:n_scenarios) {
    res <- results_list[[i]]
    
    # Plot SEIR trajectories
    plot(res$t, res$S, type = "l", col = "black", lwd = 2.2,
         main = scenario_names[i], xlab = "Day", ylab = "Population",
         ylim = c(0, max_pop),
         cex.axis = 0.9, cex.lab = 0.9, cex.main = 1.0)
    
    lines(res$t, res$E, col = "blue", lwd = 2.2)
    lines(res$t, res$I, col = "red", lwd = 2.2)
    lines(res$t, res$R, col = "green", lwd = 2.2)
    
    grid(col = "gray80", lty = "dotted", lwd = 0.8)
  }
  # Add overall legend at bottom
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0),
      new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  
  legend("bottom", ncol = 4,
         legend = c("Susceptible", "Exposed", "Infectious", "Recovered"),
         col = c("black", "blue", "red", "green"),
         lwd = 2.2, bty = "n", cex = 1.0,
         xpd = TRUE)
  
  # Add main title if provided
  if (main_title != "") {
    mtext(main_title, outer = TRUE, cex = 1.2, line = 0.5)
  }
}


##============================================================================================================
##---------------Create the vector of integers indicating which household each person belongs to--------------
set.seed(0)
n = 10000
h_max = 5

#generates household IDs repeated by random household sizes
h = rep(1:n,sample(1:h_max,n,replace = TRUE))[1:n] #[1:n] ensures the total length is exactly n by truncating any excess

#check the household size and NAs
tabulate(tabulate(h))
sum(table(table(h)))
table(table(h))
which(is.na(h)==TRUE)


##============================================================================================================
##-------------------------Setting beta, Use function,Compare scenarios and Visualization---------------------
set.seed(0)
beta = runif(n, 0, 1)
alink = get.net(beta, h, nc = 15)

Full_model = nseir(beta, h, alink, alpha = c(.1, .01, .01))  # Full model
Random_mixing_only = nseir(beta, h, alink, alpha = c(0, 0, 0.04))    # Random mixing only
Constant_beta = nseir(rep(mean(beta), n), h, alink, alpha = c(.1, .01, .01))  # Constant beta
Constant_beta_random = nseir(rep(mean(beta), n), h, alink, alpha = c(0, 0, 0.04))    # Constant beta, random only

# Plot individual scenarios
plot.nseir(Full_model, main = "Full model")
plot.nseir(Random_mixing_only, main = "Random mixing only")
plot.nseir(Constant_beta, main = "Constant beta")
plot.nseir(Constant_beta_random, main = "Constant β + random")

# Compare all scenarios side-by-side
compare.scenarios(
  list(Full_model, Random_mixing_only, Constant_beta, Constant_beta_random),
  scenario_names = c("Full model", "Random mixing", 
                     "Constant β", "Constant β + random"),
  main_title = "Effect of Structure on Epidemic Dynamics",
  show_beta = TRUE
)



##============================================================================================================
##-------------------------------------Comment by using the plots---------------------------------------------

#In the Full model and Random mixing graphs, beta is random, that is, everyone has different social abilities. 
#It can be seen that after removing the household and regular network structure, 
#each person comes into contact with different people every day. 
#This led to a higher peak in the number of infections, earlier occurrence of infections, 
#resulting in a wide spread of the epidemic, rapid speed, and a rapid decline in the number of infections. 
#In contrast, the Full model, due to having a fixed social group and more localized contact, 
#has a slower spread of the epidemic, a lower overall peak, but a longer duration of the epidemic.

#In the Constant beta and Constant beta+random models, beta is fixed, meaning that everyone's social skills are the same. 
#In Constant beta+random, without the isolation of families and social networks, the entire population was exposed, 
#with the highest number of infections, almost affecting the entire population.

#Therefore, under the influence of the household and regular network structure, the epidemic spread within a small range, 
#spread slowly among the entire population, and the peak number of infections was relatively low, but the epidemic lasted for a long time.
