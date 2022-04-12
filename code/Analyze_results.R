require(pacman)
p_load(pacman,tidyverse,MASS,evd,sampleSelection,foreach,data.table,patchwork,
       doParallel,tictoc,patchwork,matrixcalc,survival,plotly,apollo,haven,readstata13,broom,stringr,ggridges,usethis)

####################################################
#### load in results and data for 20k simulation ###
####################################################

load("~/AmpleCorrection/results/simresults20k.Rdata") # heterosked scaling using the average attribute levels in a choice set

# You should get the following objects in your environment:
# clogitdata - simulated data, reshaped for R's clogit function (includes "choice" rows for non-responders)
# database - simulated data, shaped for apollo (non-responders' choice rows are dropped)
# dataf - simulated data, shaped for apollo (including choice rows for non-responders)
# model - results of selection corrected mixed logit estimation in apollo
# right - clogit results on the full set of invitees
# wrong - clogit results for responders only (no correction)

# More notes:
# the attribute called "att" is the "benefit" variable in the main paper
# the code uses a status-quo indicator rather than an any-policy indicator
# (which makes no substantive difference)

# Let's have a look

# print the results the apollo model
apollo_modelOutput(model)



# do some checks on estimated VCOV MATRIX
# (matrix order is alpha, statquo, att, cost)

ests <- model$estimate


# eigen values positive?
with(as.list(ests),{
  vcov =  c(a1^2,  a1*b1,       a1*c1,             a1*d1,
            a1*b1, b1^2+b2^2,   b1*c1+b2*c2,       b1*d1+b2*d2,
            a1*c1, b1*c1+b2*c2, c1^2+c2^2+c3^2,    c1*d1+c2*d2+c3*d3,
            a1*d1, b1*d1+b2*d2, c1*d1+c2*d2+c3*d3, d1^2+d2^2+d3^2+d4^2)
  matrix(vcov,nrow=4)
}) %>% eigen()
# yup, looks good

# convert vcov matrix to a correlation matrix
# not as useful as you might think
with(as.list(ests),{
  vcov =  c(a1^2,  a1*b1,       a1*c1,             a1*d1,
            a1*b1, b1^2+b2^2,   b1*c1+b2*c2,       b1*d1+b2*d2,
            a1*c1, b1*c1+b2*c2, c1^2+c2^2+c3^2,    c1*d1+c2*d2+c3*d3,
            a1*d1, b1*d1+b2*d2, c1*d1+c2*d2+c3*d3, d1^2+d2^2+d3^2+d4^2)
  matrix(vcov,nrow=4)
}) %>% cov2cor()

# compare that cor matrix to the cor matrix in the simulated data
getcors <- dataf %>% filter(choicetype=="R") %>% select(starts_with("eps"))
cor(getcors)


# store the vcov matrix as "mat" and use it to draw cost, attribute, statquo parameters
mat <- with(as.list(ests),{
  vcov =  c(a1^2,  a1*b1,       a1*c1,             a1*d1,
            a1*b1, b1^2+b2^2,   b1*c1+b2*c2,       b1*d1+b2*d2,
            a1*c1, b1*c1+b2*c2, c1^2+c2^2+c3^2,    c1*d1+c2*d2+c3*d3,
            a1*d1, b1*d1+b2*d2, c1*d1+c2*d2+c3*d3, d1^2+d2^2+d3^2+d4^2)
  matrix(vcov,nrow=4)
})

drawwtp <- mvrnorm(n=20000,mu=ests[c("mu_alpha","mu_statquo","mu_att","mu_cost")],Sigma=mat)

# print estimated median and 95 CI for att WTP (should be 1)
quantile(-1*drawwtp[,3]/drawwtp[,4],c(0.025,0.5,0.975))

# print estimated median and 95 CI for statquo WTP (should be -0.2)
quantile(-1*drawwtp[,2]/drawwtp[,4],c(0.025,0.5,0.975))


# compare those to the TRUE distributions of WTPs in the simulated data (including non-responders)
singleserve <- dataf %>% filter(choicetype=="R") # one row per person (drop response/non-response rows)
responders <- singleserve %>% filter(response_status=="RESPONDER") # just the responder sample

quantile(-1*(5+singleserve$eps_att)/(-5+singleserve$eps_cost),c(0.025,0.5,0.975)) # WTP for att
quantile(-1*(-1+singleserve$eps_statquo)/(-5+singleserve$eps_cost),c(0.025,0.5,0.975)) # WTP for statquo

# plot those distributions
# red is responders
# green is everyone
# blue is estimates
# vertical lines are medians, dashed line is for estimates

drawwtpdf <- as.data.frame(drawwtp)

estimatedwtp <- drawwtpdf %>% 
  mutate(wtp_ben = -mu_att/mu_cost) %>% 
  mutate(wtp_statquo = -mu_statquo/mu_cost) %>% 
  select(wtp_ben,wtp_statquo) %>% 
  mutate(Type="Estimate")

simulatedwtp <- singleserve %>% 
  mutate(wtp_ben=-(5+eps_att)/(-5+eps_cost)) %>% 
  mutate(wtp_statquo=-(1+eps_statquo)/(-5+eps_cost)) %>% 
  mutate(Type="Responders &\nNon-responders") %>% 
  select(wtp_ben,wtp_statquo,Type)

forplotting <- rbind(estimatedwtp,simulatedwtp)


# plots

dataf$response_status <- ifelse(dataf$response_status=="RESPONDER","Responders",dataf$response_status)
dataf$response_status <- ifelse(dataf$response_status=="NONRESPONDER","Non-responders",dataf$response_status)


p1 <- ggplot(dataf) +
  stat_density_ridges(aes(x=5+eps_att,y=response_status,fill=response_status,color=response_status),
                      alpha=1,show.legend = F,quantile_lines = TRUE, quantiles = 2) +
  labs(x="βbenefit",y="") +
  scale_fill_manual(values=c("#f0f037","#1ABDB7")) +
  scale_color_manual(values=c("#84841e","#0c5552")) +
  lims(x=c(-10,10))+ 
  theme_minimal()

p2 <- ggplot(dataf) +
  stat_density_ridges(aes(x=-5+eps_cost,y=response_status,fill=response_status,color=response_status),
                      alpha=1,show.legend = F,quantile_lines = TRUE, quantiles = 2) +
  labs(x="βcost",y="") +
  lims(x=c(-10,10))+
  scale_fill_manual(values=c("#f0f037","#1ABDB7")) +
  scale_color_manual(values=c("#84841e","#0c5552")) +
  theme_ridges(center_axis_labels = TRUE) +
  theme_minimal()

p3 <- ggplot(dataf) +
  stat_density_ridges(aes(x=-(5+eps_att)/(-5+eps_cost),y=response_status,fill=response_status,color=response_status),
                      alpha=1,show.legend = F,quantile_lines = TRUE, quantiles = 2) +
  labs(x="",y="") +
  scale_fill_manual(values=c("#f0f037","#1ABDB7")) +
  scale_color_manual(values=c("#84841e","#0c5552")) +
  theme_minimal() +
  lims(x=c(0,3))


p4 <- ggplot(forplotting) +
  stat_density_ridges(aes(x=wtp_ben,y=Type,fill=Type,color=Type),
                      quantile_lines = TRUE, quantiles = 2,show.legend = F) +
  scale_fill_manual(values=c("#af1abd","#85d25f")) +
  scale_color_manual(values=c("#460a4c","#355426")) +
  theme_minimal() +
  labs(x="Marginal WTP for benefit",y="") +
  lims(x=c(0,3))

(p1/p2)

(p3/p4)

####### total wtp ############

atta <- 1.5 # pick a level of "benefit" 

dataf <- dataf %>% mutate(twtp_i = -atta*(5+eps_att)/(eps_cost-5) + (eps_statquo-1)/(eps_cost-5))
choices <- dataf %>% filter(choicetype=="S")
choicesr <- choices %>% filter(response_status=="Responders")
choicesnr <- choices %>% filter(response_status=="Non-responders")
summary(choices$twtp_i)
summary(choicesr$twtp_i)
summary(choicesnr$twtp_i)

simchoices <- drawwtpdf %>% mutate(twtp_i = (-mu_att/mu_cost)*atta + mu_statquo/mu_cost)
summary(simchoices$twtp_i)


twtp <- rbind(
  choices %>% mutate(type="Responders &\nNon-responders") %>% select(type,twtp_i),
  simchoices %>% mutate(type="Estimate") %>% select(type,twtp_i),
  choicesr %>% mutate(type="Responders") %>% select(type,twtp_i),
  choicesnr %>% mutate(type="Non-responders") %>% select(type,twtp_i)
)


p5 <- ggplot(twtp[twtp$type %in% c("Responders","Non-responders"),]) +
  stat_density_ridges(aes(x=twtp_i,y=type,fill=type,color=type),
                      alpha=1,show.legend = F,quantile_lines = TRUE, quantiles = 2) +
  labs(x="",y="") +
  scale_fill_manual(values=c("#f0f037","#1ABDB7")) +
  scale_color_manual(values=c("#84841e","#0c5552")) +
  theme_minimal() +
  lims(x=c(0,6))


p6 <- ggplot(twtp[twtp$type %in% c("Responders &\nNon-responders","Estimate"),]) +
  stat_density_ridges(aes(x=twtp_i,y=type,fill=type,color=type),
                      alpha=1,show.legend = F,quantile_lines = TRUE, quantiles = 2) +
  scale_fill_manual(values=c("#af1abd","#85d25f")) +
  scale_color_manual(values=c("#460a4c","#355426")) +
  theme_minimal() +
  labs(x="Total WTP (benefit = 1.5)",y="") +
  lims(x=c(0,6))

p5/p6


twtpstats <- twtp %>% 
  filter(type!="Non-responders") %>%  
  group_by(type) %>% 
  summarise(mean=mean(twtp_i),
            med=median(twtp_i),
            n=n(),
            se=sd(twtp_i)/sqrt(n))

# bootstrap to get standard error of medians because no closed-form formula

m <- 5000
estmed <- c()
allmed <- c()
resmed <- c()
for (i in 1:m) {
  estmed <- c(estmed,median(sample(twtp$twtp_i[twtp$type=="Estimate"],20000,replace = T)))
  allmed <- c(allmed,median(sample(twtp$twtp_i[twtp$type=="Responders &\nNon-responders"],20000,replace = T)))
  resmed <- c(resmed,median(sample(twtp$twtp_i[twtp$type=="Responders"],9947,replace = T)))
  
}
medse <- c(sd(estmed),sd(resmed),sd(allmed))
medse

twtpstats$se_med <- medse

mean(estmed-allmed)/sd(estmed-allmed)




##############################################
############## many small runs ###############
##############################################

res <- readRDS("~/AmpleCorrection/results/100runsof2k.rds")

res <- res[2:nrow(res),] # first row is trash


wtpatt <- NULL
wtpstatquo <- NULL

for (i in 1:nrow(res)) {
  
  mat <- with(as.list(res[i,]),{
    vcov =  c(a1^2,  a1*b1,       a1*c1,             a1*d1,
              a1*b1, b1^2+b2^2,   b1*c1+b2*c2,       b1*d1+b2*d2,
              a1*c1, b1*c1+b2*c2, c1^2+c2^2+c3^2,    c1*d1+c2*d2+c3*d3,
              a1*d1, b1*d1+b2*d2, c1*d1+c2*d2+c3*d3, d1^2+d2^2+d3^2+d4^2)
    matrix(vcov,nrow=4)
  })
  
  drawwtp <- mvrnorm(n=50000,mu=unlist(res[i,c("mu_alpha","mu_statquo","mu_att","mu_cost")]),Sigma=mat)
  
  
  wtpatt <- rbind(wtpatt,quantile(-1*drawwtp[,3]/drawwtp[,4],c(0.025,0.5,0.975)))
  wtpstatquo <- rbind(wtpstatquo,quantile(-1*drawwtp[,2]/drawwtp[,4],c(0.025,0.5,0.975)))
  
}

wtpatt <- wtpatt %>% as.data.frame()
wtpstatquo <- wtpstatquo %>% as.data.frame()

names(wtpatt) <- c("lb","corrected","ub")

mean(wtpatt$corrected)
median(wtpatt$corrected)

quantile(wtpatt$`50%`,c(0.025,.5,0.975))

wtpatt$wrong <- -1*res$wrong_att/res$wrong_cost
wtpatt$right <- -1*res$right_att/res$right_cost

wtpatt2 <- wtpatt %>% pivot_longer(cols=c("right","wrong","corrected"))

p100 <- ggplot(wtpatt2) +
  stat_density_ridges(aes(x=value,y=name,group=name,fill=name,color=name),scale=2.5,
                      quantile_lines = TRUE,quantiles = 2,show.legend = F) +
  scale_fill_manual(values=c("#af1abd","#85d25f","#1ABDB7")) +
  scale_color_manual(values=c("#460a4c","#355426","#0c5552")) +
  scale_y_discrete(labels = c('Corrected','Clogit\n(full sample)','Clogit\n(responders only)')) +
  theme_minimal() +
  labs(x="WTP for benefit",y="")


p100


