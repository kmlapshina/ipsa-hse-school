library(haven)
library(dplyr)
library(ggplot2)
library(knitr)
library(MCMCpack)
library(bayesplot)
library(rstanarm)
library(broom.mixed)
library(coda)

#--------READING THE FILE, PREPARING THE DATA--------#

data_de <- read_sav('data-ipsa-de.sav') %>% dplyr:::select(imwbcnt, trstprl, 
                                                           trstlgl, trstplc, 
                                                           ppltrst, pplfair,
                                                           pplhlp)

data_de <- na.omit(data_de)
summary(data_de)

#--------RStanARM fitting---------#

stan_reg <- stan_glm(imwbcnt ~ trstprl + trstlgl + trstplc + 
                         ppltrst + pplfair + 
                         pplhlp,
                       data = data_de, iter = 200000, 
                       chains = 4, thin=10)

kable(tidyMCMC(stan_reg,
               conf.int = TRUE, conf.method = 'quantile'),
      digits = 3)

mcmc_areas(as.matrix(stan_reg),
           pars = c('trstprl', 'trstlgl', 'trstplc', 
                    'ppltrst', 'pplfair', 'pplhlp'))

###-------------------------------------------------###
###---------  CONVERGENCE DIAGNOSTICS  -------------###
###-------------------------------------------------###


color_scheme_set('mix-orange-darkgray')
bayesplot_theme_update(text = element_text(size = 10))

posterior_stan <- as.array(stan_reg)

mcmc_dens_overlay(posterior_stan,
                  pars = c('trstprl', 'trstlgl', 'trstplc', 
                           'ppltrst', 'pplfair', 'pplhlp'))

mcmc_trace(posterior_stan, 
           pars = c('trstprl', 'trstlgl', 'trstplc', 
                    'ppltrst', 'pplfair', 'pplhlp'))

mcmc_rhat(rhat = rhat(stan_reg)) + 
  yaxis_text(hjust = 0)

kable(rhat(stan_reg)) # preferably less than 1.05

#The Rhat function produces R-hat convergence diagnostic, 
#which compares the between- and within-chain estimates 
#for model parameters and other univariate quantities 
#of interest. If chains have not mixed well 
#(ie, the between- and within-chain estimates don't agree), 
#R-hat is larger than 1. 


heidel.diag(stan_reg) #Stationarity test 
# If passed then the chain has 
# converged to a stationary distribution

#Autocorrelation
plot(stan_reg, "acf", c('trstprl', 'trstlgl', 'trstplc', 
                          'ppltrst', 'pplfair', 'pplhlp')) # plot

mcmc_acf_bar(posterior_stan, 
             pars = c('trstprl', 'trstlgl', 'trstplc', 
                      'ppltrst', 'pplfair', 'pplhlp')) # chart

#Effective samples
mcmc_neff(neff_ratio(stan_reg, size = 8))

#If the ESS of a parameter is small then the estimate of the posterior 
#distribution of that parameter will be poor. In Tracer you can calculate 
#the standard deviation of the estimated mean of a parameter. 
#If the ESS is small then the standard deviation will be large. 



###-------------------------------------------------###
###-------  POSTERIOR PREDICTIVE CHECKS  -----------###
###-----------------  bayesplot  -------------------###
###-------------------------------------------------###

# generate replicated data
yrep <- posterior_predict(stan_reg, draws = 500) 


ppc_dens_overlay(as.numeric(data_de$imwbcnt), yrep[1:50, ]) +
  theme(legend.text = element_text(face="bold", 
                                   size=20)) +
  theme(axis.text.x = element_text(face="bold", 
                                   size=14)) +
  theme(axis.text.y = element_text(face="bold", 
                                   size=14))


# histograms of observed and rep. data
ppc_hist(as.numeric(data_de$imwbcnt), yrep[1:5, ]) +
  theme(legend.text = element_text(face="bold", 
                                   size=20)) +
  theme(axis.text.x = element_text(face="bold", 
                                   size=14)) +
  theme(axis.text.y = element_text(face="bold", 
                                   size=14))


# Observed mean vs. replicated means
ppc_stat(as.numeric(data_de$imwbcnt), yrep, stat = "mean") + 
  theme(legend.title=element_text(size=20),
        legend.text=element_text(size=10),
        axis.text=element_text(size=10))

# Observed mean vs. replicated means
ppc_stat(as.numeric(data_de$imwbcnt), yrep, stat = "median") + 
  theme(legend.title=element_text(size=20),
        legend.text=element_text(size=10),
        axis.text=element_text(size=10))


# Observed variance vs. replicated variances
ppc_stat(as.numeric(data_de$imwbcnt), yrep, stat = "sd") + 
  theme(legend.title=element_text(size=20),
        legend.text=element_text(size=10),
        axis.text=element_text(size=10))



###-------------------------------------------------###
###-------  POSTERIOR PREDICTIVE CHECKS  -----------###
###------------------  rstanarm  -------------------###
###-------------------------------------------------###

#2-D Hist
pp_check(stan_reg, plotfun = "stat_2d", stat = c("mean", "sd")) + 
  theme(legend.title=element_text(size=20),
        legend.text=element_text(size=10),
        axis.text=element_text(size=10))


# Average predictive error
pp_check(stan_reg, plotfun='error_scatter_avg')+ 
  theme(text=element_text(size=18))




###-------------------------------------------------###
###--------------  COMPARING MODELS  ---------------###
###-------------------------------------------------###

library(BayesFactor)

bf_nest1 <- lmBF(imwbcnt ~ trstprl + trstlgl + trstplc + 
                       ppltrst + pplfair,
                     data = data_de, iterations=1000)

bf_nest2 <- lmBF(imwbcnt ~ trstprl + trstlgl + trstplc + 
                       ppltrst + pplhlp,
                     data = data_de)

bf_full <- lmBF(imwbcnt ~ trstprl + trstlgl + trstplc + 
                   ppltrst + pplfair + 
                   pplhlp,
                 data = data_de, iterations=1000)


bf_nest1/bf_full
bf_nest2/bf_full


######################


## FINAL MODEL

stan_reg_1 <- stan_glm(imwbcnt ~ trstprl + trstlgl + trstplc + 
                       ppltrst + pplfair,
                     data = data_de, iter = 200000, 
                     chains = 4, thin=10)

kable(tidyMCMC(stan_reg_1,
               conf.int = TRUE, conf.method = 'quantile'),
      digits = 3)

mcmc_areas(as.matrix(stan_reg_1),
           pars = c('trstprl', 'trstlgl', 'trstplc', 
                    'ppltrst', 'pplfair'))
