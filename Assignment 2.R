#Group 6 :Xuyi Shi(S2796586), Chi Zhang(S2828433), Jiachen Guang(S2789777)
#Contribution:



#This code implements an extended SEIR epidemic model that incorporates social 
#network structure, used to investigate the impact of family relationships and 
#individual differences in sociability on the spread of infectious diseases. 
#The model distributes the population into households of random size and 
#constructs a contact network based on sociability parameters. It simulates three 
#infection pathways: intra-household transmission, transmission with fixed contacts, 
#and random mixing. By comparing four different scenarios (a complete model, 
#a model without social structure, a model with homogeneous sociability, and a baseline model), 
#we analyze how social structure influences the dynamics, scale, and peak of the epidemic.



#create_households(): 创建家庭结构，将人口分配到随机规模的家庭中

#get.net(): 基于个体社交性参数构建社交网络

#nseir(): 核心模拟函数，运行包含社会结构的SEIR模型

#plot_epidemic_curve(): 可视化疫情发展曲线

#情景分析: 比较不同社会结构设置对疫情的影响

##============================================================================================================
##---------------Create the vector of integers indicating which household each person belongs to--------------
set.seed(0)
n = 1000
h_max = 5

h = rep(1:n,sample(1:h_max,n,replace = TRUE))[1:n]
tabulate(tabulate(h))
sum(table(table(h)))
table(table(h))
which(is.na(h)==TRUE)


##============================================================================================================
##------------------Create the net funtion between regular (non-household) contact people---------------------
get.net = function(beta, h, nc = 15) {
  # Creates regular contact network based on sociability parameters
  # beta: vector of sociability parameters for each person
  # h: household assignment vector
  # nc: average number of contacts per person
  # Returns: list where element i contains indices of i's regular contacts
  
  n = length(beta)  # population size
  beta_bar = mean(beta)  # mean sociability
  
  # Initialize empty contact list for each person
  contacts = vector("list", n)
  
  # Loop through all possible pairs (i,j) where i < j to avoid duplicates
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      
      # Skip if i and j are in the same household
      if (h[i] == h[j]) next
      
      # Calculate probability of link between i and j
      prob = nc * beta[i] * beta[j] / (beta_bar^2 * (n - 1))
      
      # Create link with calculated probability
      if (runif(1) < prob) {
        # Record link in both directions (undirected network)
        contacts[[i]] = c(contacts[[i]], j)
        contacts[[j]] = c(contacts[[j]], i)
      }
    }
  }
  
  return(contacts)
}
##============================================================================================================
##-----------------------------------Create a function of infections------------------------------------------

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
  S = numeric(nt + 1)
  E = numeric(nt + 1)
  I = numeric(nt + 1)
  R = numeric(nt + 1)
  t = 0:nt
  
  # Record initial state
  S[1] = sum(state == 1)
  E[1] = sum(state == 2)
  I[1] = sum(state == 3)
  R[1] = sum(state == 4)
  
  # Simulate epidemic dynamics
  for (day in 1:nt) {
    
    new_state = state  # copy current state
    
    # Find individuals in each state
    infectious = which(state == 3)
    exposed = which(state == 2)
    susceptible = which(state == 1)
    
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
            infected = household_members[runif(length(household_members)) < alpha_h]
            new_state[infected] = 2
          }
        }
        
        # 2. Regular contact network infections
        if (alpha_c > 0 && length(alink[[i]]) > 0) {
          # Find contacts who are susceptible
          susceptible_contacts = alink[[i]][state[alink[[i]]] == 1]
          if (length(susceptible_contacts) > 0) {
            infected <- susceptible_contacts[runif(length(susceptible_contacts)) < alpha_c] 
            new_state[infected] = 2
          }
        }
        
        # 3. Random mixing infections
        if (alpha_r > 0) {
          # Calculate infection probability for each susceptible person
          susceptible_now = which(state == 1)
          if (length(susceptible_now) > 0) {
            prob <- alpha_r * nc * beta[i] * beta[susceptible_now] / (beta_bar^2 * (n - 1))
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

#test
beta <- runif(n, 0, 1)
alink <- get.net(beta, h, nc = 15)
results <- nseir(beta, h, alink)

##============================================================================================================
##-----------------------------Create a function to plot the dynamics states----------------------------------




