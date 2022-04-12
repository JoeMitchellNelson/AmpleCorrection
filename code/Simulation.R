require(pacman)
p_load(pacman,tidyverse,MASS,evd,sampleSelection,foreach,data.table,
       doParallel,tictoc,patchwork,matrixcalc,survival,plotly,apollo,haven,readstata13,broom,stringr)


# this code handles data generation and estimation for the sim results

# Running this simulation with 15 threads in parallel (AMD Ryzen 7) takes about 22 to 23 hours
# Recommend reducing n to something smaller than 20k, then go for a walk and get some coffee
n <- 20000 #number of simulated invitees


#load helper functions
source("~/AmpleCorrection/code/selection_helpers.R")


# names of correlated variables in selection equation (alpha, in this case)
# and in the choice equation: statusquo, benefit (called att in this code), and cost
data1 <- c("alpha","statusquo","att","cost")

# adds a function to the global environment called "rescaler"
# that is used to rescale the index in the choice equation
build_rescale_function(corvars=data1)


# set true param values
alpha <- 0
a_iv <- 1
b_statquo <- -1
b_cost <- -5
b_att <- 5

# utility functions
u_respond <- function (iv,eps_alpha,omega) {(alpha + eps_alpha) + a_iv*iv + omega}
u_choice <- function (statquo,eps_statquo,att,eps_att,cost,eps_cost,eta) {
  (b_statquo+eps_statquo)*statquo + (b_att+eps_att)*att + (b_cost+eps_cost)*cost + eta
}


### for the correlation matrix for data generation

sigma_alpha_beta_11 <- 0.5 # alpha constant and statquo
sigma_alpha_beta_12 <- 0.5 # alpha constant and att
sigma_alpha_beta_13 <- 0.5 # alpha constant and cost
sigma_beta_12 <- 0 # statquo and att
sigma_beta_13 <- 0 # statquo and cost
sigma_beta_23 <- 0 # att and cost




set.seed(1234)

# draw the correlated components
Sigma <- matrix(nrow=4,ncol=4,byrow=T,
                data=c(1,sigma_alpha_beta_11,sigma_alpha_beta_12,sigma_alpha_beta_13,
                       sigma_alpha_beta_11,1,sigma_beta_12,sigma_beta_13,
                       sigma_alpha_beta_12,sigma_beta_12,1,sigma_beta_23,
                       sigma_alpha_beta_13,sigma_beta_13,sigma_beta_23,1))

# eta is the policy choice error
# omega is the response/non-response error
# epsilons are random components of utility parameters

epsilons <- mvrnorm(n=n,mu=c(0,0,0,0),Sigma=Sigma)


dataf <- data.frame(id=1:n,
                    eta1=rlogis(n), # independent draw
                    eta2=rlogis(n), # independent draw
                    omega=rlogis(n), # independent draw
                    eps_alpha=epsilons[,1], # correlated with choice parameters
                    eps_statquo=epsilons[,2], # correlated with alpha
                    eps_att=epsilons[,3], # correlated with alpha
                    eps_cost=epsilons[,4], # correlated with alpha
                    iv=rnorm(n)) # independent draw

dataf <- rbind(dataf,dataf) %>% arrange(id) # two rows per person

# generate alternatives (all available alternatives on a single row for apollo)
dataf <- dataf %>% mutate(choicetype=rep(c("R","S"),n),
                          cost1=rep(runif(n,0,3),each=2),
                          cost2=rep(runif(n,0,3),each=2),
                          att1 = rep(runif(n,0,3),each=2),
                          att2 = rep(runif(n,0,3),each=2),
                          statquo=-1)

# calculate utilities for r/nr and policy-choice alternatives-- u(statquo)=0, u(non-response)=0
dataf <- dataf %>% mutate(u_respond1=u_respond(iv,eps_alpha,omega),
                          u_choice1=u_choice(statquo,eps_statquo,att1,eps_att,cost1,eps_cost,eta1),
                          u_choice2=u_choice(statquo,eps_statquo,att2,eps_att,cost2,eps_cost,eta2))

# choice key
# 1 - don't respond
# 2 - respond
# 3 - choose statquo
# 4 - choose alt 1
# 5 - choose alt 2

# make the choices
dataf$choice <- ifelse(dataf$choicetype=="R" & dataf$u_respond1 < 0,1,NA)
dataf$choice <- ifelse(dataf$choicetype=="R" & dataf$u_respond1 >= 0,2,dataf$choice)
dataf$choice <- ifelse(dataf$choicetype=="S" & dataf$u_choice1 <= 0 & dataf$u_choice2 <= 0,3,dataf$choice)
dataf$choice <- ifelse(dataf$choicetype=="S" & dataf$u_choice1 > 0 & dataf$u_choice1 >= dataf$u_choice2,4,dataf$choice)
dataf$choice <- ifelse(dataf$choicetype=="S" & dataf$u_choice2 > 0 & dataf$u_choice2 > dataf$u_choice1,5,dataf$choice)



# label individuals responders/nonresponders
dataf <- dataf %>% group_by(id) %>%
  mutate(response_status = min(choice)) %>%
  ungroup %>%
  mutate(response_status = ifelse(response_status==1,"NONRESPONDER","RESPONDER"))

# apollo requires "availability" dummies
dataf$av1 <- ifelse(dataf$choice %in% 1:2,1,0)
dataf$av2 <- ifelse(dataf$choice %in% 1:2,1,0)
dataf$av3 <- ifelse(dataf$choice %in% 3:5,1,0)
dataf$av4 <- ifelse(dataf$choice %in% 3:5,1,0)
dataf$av5 <- ifelse(dataf$choice %in% 3:5,1,0)

# reshape data for clogit

clogitdata <- dataf %>% filter(choicetype=="S") %>% mutate(cost0=0,att0=0)

clogitdata <- clogitdata %>% pivot_longer(cols=starts_with("cost"),names_to="cost_drop",values_to="cost") %>%
  pivot_longer(cols=starts_with("att"),names_to="att_drop",values_to="att") %>%
  select(id,choice,response_status,att,cost,att_drop,cost_drop) %>%
  mutate(att_drop=str_extract(att_drop,"[[:digit:]]"),
         cost_drop=str_extract(cost_drop,"[[:digit:]]")) %>%
  filter(att_drop==cost_drop) %>%
  mutate(alt=as.numeric(att_drop)+3) %>%
  mutate(choice=ifelse(choice==alt,1,0)) %>%
  mutate(statquo=ifelse(alt %in% 4:5,-1,0))

# clogit everyone
right <- clogit(choice ~ att + cost + statquo + strata(id),data=clogitdata)
# clogit responders only (uncorrected)
wrong <- clogit(choice ~ att + cost + statquo + strata(id),data=clogitdata[clogitdata$response_status=="RESPONDER",])
summary(right)
summary(wrong)


df <- dataf %>% filter(!(response_status=="NONRESPONDER" & choicetype=="S"))

database <- df %>% ungroup() %>% as.data.frame()

### Initialise code
apollo_initialise()

### Set core controls
### NOTE: YOU MAY NEED TO LOWER NCORES BEFORE RUNNING. MOST COMPUTERS MAX OUT AT 8, IN WHICH CASE YOU SHOULD USE NO MORE THAN 7.
### ALSO: "NCORES" SHOULD REALLY BE CALLED "NTHREADS". TYPICAL TO HAVE 2 THREADS PER CORE.

apollo_control <- list(
  modelName ="Selection correction",
  modelDescr ="Mixed logit model",
  indivID   ="id",  
  mixing    = TRUE,
  nCores    = 7,
  panelData = TRUE
)

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #

## Vector of parameters, including any that are kept fixed in estimation
apollo_beta <- c(
  #response variables
  
  alpha_iv = 5,
  mu_alpha = 0,
  
  # choice variables
  
  mu_att = 3,
  mu_statquo = 0,
  mu_cost = -3,
  beta_imr = -1.7,
  
  
  
  # covariances
  # var alpha = a1^2
  # var statquo = b1^2 + b2^2
  # var att = c1^2 + c2^2 + c3^2
  # var cost = d1^2 + d2^2 + d3^2 + d4^2
  
  # cov alpha, statquo = a1*b1
  # cov alpha, att = a1*c1
  # cov alpha, cost = a1*d1
  
  # cov statquo, att = b1*c1 + b2*c2
  # cov statquo, cost = b1*d1 + b2*d2
  
  # cov att, cost = c1*d1 + c2*d2 + c3*d3
  
  
  # VCOV MATRIX, order is alpha, statquo, att, cost
  
  # vcov =  c(a1^2,  a1*b1,       a1*c1,             a1*d1,
  #           a1*b1, b1^2+b2^2,   b1*c1+b2*c2,       b1*d1+b2*d2, 
  #           a1*c1, b1*c1+b2*c2, c1^2+c2^2+c3^2,    c1*d1+c2*d2+c3*d3,
  #           a1*d1, b1*d1+b2*d2, c1*d1+c2*d2+c3*d3, d1^2+d2^2+d3^2+d4^2)
  
  # VCOV for just the k x k matrix of choice parameters
  
  # vcovkxk =  c(b1^2+b2^2,   b1*c1+b2*c2,       b1*d1+b2*d2, 
  #              b1*c1+b2*c2, c1^2+c2^2+c3^2,    c1*d1+c2*d2+c3*d3,
  #              b1*d1+b2*d2, c1*d1+c2*d2+c3*d3, d1^2+d2^2+d3^2+d4^2)
  
  
  
  
  a1 = 1,
  b1 = .7,
  b2 = .3,
  c1 = .7,
  c2 = 1,
  c3 = 0,
  d1 = 1,
  d2 = 0,
  d3 = 0,
  d4 = 1,
  
  # scales
  
  scale_RP =  1,
  scale_SP = 1
)



#apollo_beta <- startpars

### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
apollo_fixed <- c("scale_RP","scale_SP")

# ################################################################# #
#### DEFINE RANDOM COMPONENTS                                    ####
# ################################################################# #

### Set parameters for generating draws
apollo_draws <- list(
  interDrawsType = "halton",
  interNDraws    = 500,
  interUnifDraws = c(),
  interNormDraws = c("draws_alpha","draws_statquo","draws_att","draws_cost"),
  intraDrawsType = "halton",
  intraNDraws    = 0,
  intraUnifDraws = c(),
  intraNormDraws = c()
)




### Create random parameters
apollo_randCoeff <- function(apollo_beta, apollo_inputs){
  
  randcoeff = list()
  
  randcoeff[["random_alpha"]]   = mu_alpha +   a1*draws_alpha
  
  randcoeff[["random_statquo"]] = mu_statquo + b1*draws_alpha + b2*draws_statquo
  
  randcoeff[["random_att"]]     = mu_att +     c1*draws_alpha + c2*draws_statquo + c3*draws_att
  
  randcoeff[["random_cost"]]    = mu_cost +    d1*draws_alpha + d2*draws_statquo + d3*draws_att + d4*draws_cost
  
  
  return(randcoeff)
}





# ################################################################# #
#### GROUP AND VALIDATE INPUTS                                   ####
# ################################################################# #

apollo_inputs <- apollo_validateInputs(silent=F)

# ################################################################# #
#### DEFINE MODEL AND LIKELIHOOD FUNCTION                        ####
# ################################################################# #



apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality="estimate"){
  
  ### Function initialisation: do not change the following three commands
  ### Attach inputs and detach after function exit
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  ### Create list of probabilities P
  P = list()
  
  ### List of utilities: these must use the same names as in mnl_settings, order is irrelevant
  
  # 1 - don't respond
  # 2 - respond
  # 3 - status quo
  # 4 - take the policy
  
  
  
  V = list()
  V[['alt1']] = 0
  V[['alt2']] = random_alpha + alpha_iv * iv
  V[['alt3']] = 0
  

  V[['alt4']] = (random_cost*cost1 + random_att*att1 + random_statquo*statquo +
                   beta_imr*imr((alpha_iv*iv + mu_alpha)/sqrt(a1^2+(pi^2/3))))/
    rescaler(a1,b1,b2,c1,c2,c3,d1,d2,d3,d4,0,statquo,(att1+att2)/2,(cost1+cost2)/2)
  V[['alt5']] = (random_cost*cost2 + random_att*att2 + random_statquo*statquo +
                   beta_imr*imr((alpha_iv*iv + mu_alpha)/sqrt(a1^2+(pi^2/3))))/
    rescaler(a1,b1,b2,c1,c2,c3,d1,d2,d3,d4,0,statquo,(att1+att2)/2,(cost1+cost2)/2)
  

  
  
  
  mnl_settings = list(
    alternatives  = c(alt1=1, alt2=2, alt3=3, alt4=4,alt5=5),
    avail         = list(alt1=av1, alt2=av2,alt3=av3,alt4=av4,alt5=av5),
    choiceVar     = choice,
    V             = lapply(V,"*",scale_RP),
    rows          = (choicetype=="R"),
    componentName = "MNL-RP"
  )
  
  P[['RP']] = apollo_mnl(mnl_settings, functionality)
  
  ### Compute probabilities for the SP part of the data using MNL model
  mnl_settings$V = lapply(V, "*", scale_SP)
  mnl_settings$rows = (choicetype=="S")
  mnl_settings$componentName = "MNL-SP"
  
  P[['SP']] = apollo_mnl(mnl_settings, functionality)
  
  ### Combined model
  P = apollo_combineModels(P, apollo_inputs, functionality)
  
  ### Compute probabilities using MNL model
  #P[["model"]] = apollo_mnl(mnl_settings, functionality)
  
  ### Take product across observation for same individual
  P = apollo_panelProd(P, apollo_inputs, functionality)
  
  ### Average across inter-individual draws
  P = apollo_avgInterDraws(P, apollo_inputs, functionality)
  
  # use selectionpopweight as weights
  #P = apollo_weighting(P, apollo_inputs, functionality)
  
  
  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

# ################################################################# #
#### MODEL ESTIMATION                                            ####
# ################################################################# #


model <- apollo_estimate(apollo_beta, apollo_fixed,
                         apollo_probabilities, apollo_inputs,
                         estimate_settings=list(
                           hessianRoutine="numDeriv",
                           # hessianRoutine="none",
                           silent=F))

startpars <- model$estimate

# ################################################################# #
#### MODEL OUTPUTS                                               ####
# ################################################################# #

# ----------------------------------------------------------------- #
#---- FORMATTED OUTPUT (TO SCREEN)                               ----
# ----------------------------------------------------------------- #

apollo_modelOutput(model)

#save(dataf, database,clogitdata,wrong,right,model, file = "~/AmpleCorrection/results/simresults20k.RData")