
library(ergm.count)
library(ergm.rank)
library(latentnet)

library('tidyverse')
library('readxl')
library('statnet')

source('/Users/ravigoyal/Dropbox/Academic/Research/Projects/Viral_migration/_vici/codes/VMN/VMN_func_git.R', echo=TRUE)

directory_loc = "/Users/ravigoyal/Dropbox/Academic/Research/Projects/Viral_migration/_vici/dta/allintactsequences"

VMN.df = read_csv(file = paste(directory_loc, "_LastGift_DTA_intact.baselinesequences.2021-08-02.csv", sep="/")) 
Tissue.df = read_csv(file = paste(directory_loc, "_allpredictors.2021-05-31.csv", sep="/")) 
Mapping.df = read_csv(file = paste(directory_loc, "_mapppingID_colortable_all_withGroups.csv", sep="/")) 

Mapping.df$compartment = recode(Mapping.df$compartment, 
                                RIGHTCOLON = "COLONRIGHT",
                                PERITRACHLYMPHNODE = "LNPERITRACHEAL",
                                LEFTCOLON = "COLONLEFT"
                                
)

Tissue.df = left_join(Tissue.df, Mapping.df %>% select(compartment, Groups), by = c("location" = "compartment"))

print("########################")
print("########################")
print("########################")

VMN_est.df = data.frame(pid = NULL,
                        param = NULL,
                        exp = NULL,
                        est = NULL,
                        stderr = NULL,
                        pval = NULL,
                        model = NULL)

models_run = c(1)

for (pid in unique(VMN.df$pid)) {
  
  print("########################")
  print(pid)
  print("########################")
  
  VMN_ind.net = VMN_create_missing(VMN.df, 
                                   Tissue.df, 
                                   pid_ind = pid,
                                   missing_edges = TRUE,
                                   VMN_wgt_metric = "log_totalCount",
                                   TCR_wgt_metric = "numbersharedclones") 
  
  VMN_est_ind.df = data.frame(pid = NULL,
                              param = NULL,
                              exp = NULL,
                              est = NULL,
                              stderr = NULL,
                              pval = NULL,
                              model = NULL)
  
  VMN_est_ind.df = bind_rows(VMN_est_ind.df, data.frame(pid = pid,
                                                        param = c("Nodes", "Edges", "Density"),
                                                        exp = rep(0,3),
                                                        est = rep(0,3),
                                                        pval = c(network.size(VMN_ind.net),
                                                                 as.numeric(summary(VMN_ind.net ~ edges)),
                                                                 as.numeric(summary(VMN_ind.net ~ density))
                                                        )
  )
  )
  
  #Model mutual (1)
  if (1 %in% models_run) {
    VMN_ind.net.est <- ergm(VMN_ind.net ~ sum + nonzero + nodeocov("sequence") + mutual("min"),
                            reference = ~Poisson,
                            response = "events") #transitive, triangle, isolates, idegree1.5, gwesp, gwdsp, gwodegree, twopath
    x = summary(VMN_ind.net.est)
    
    VMN_est_ind.df = bind_rows(VMN_est_ind.df, data.frame(pid = pid,
                                                          param = x$coefficients %>% row.names(),
                                                          exp = 0,
                                                          est = x$coefficients[,"Estimate"],
                                                          stderr = x$coefficients[,"Std. Error"],
                                                          pval = x$coefficients[,"Pr(>|z|)"],
                                                          model = 1))
    
  }
  print("########################")
  print("Model 1")
  print("########################")
  
  #Model homophily (2)
  if (2 %in% models_run) {
    VMN_ind.net.vest <- ergm(VMN_ind.net ~ sum + nonzero + nodeocov("sequence") + mutual("min") +
                               nodematch("group", diff = FALSE, form = "sum"),
                             reference = ~Poisson,
                             response = "events")
    
    x = summary(VMN_ind.net.vest)
    
    VMN_est_ind.df = bind_rows(VMN_est_ind.df, data.frame(pid = pid,
                                                          param = x$coefficients %>% row.names(),
                                                          exp = 0,
                                                          est = x$coefficients[,"Estimate"],
                                                          stderr = x$coefficients[,"Std. Error"],
                                                          pval = x$coefficients[,"Pr(>|z|)"],
                                                          model = 2))
  }
  print("########################")
  print("Model 2")
  print("########################")
  
  ##Model triads (3)
  if (3 %in% models_run) {
    VMN_ind.net.vest <- ergm(VMN_ind.net ~ sum + nonzero + mutual("min") +
                               nodeocov("sequence") +
                               transitiveweights("min", "max", "min") +
                               cyclicalweights("min", "max", "min"),
                             reference = ~Poisson,
                             response = "events")
    
    x = summary(VMN_ind.net.vest)
    
    VMN_est_ind.df = bind_rows(VMN_est_ind.df, data.frame(pid = pid,
                                                          param = x$coefficients %>% row.names(),
                                                          exp = 0,
                                                          est = x$coefficients$Estimate,
                                                          pval = x$coefficients$`Pr(>|z|)`,
                                                          model = 3))
  }
  print("########################")
  print("Model 3")
  print("########################")
  
  ##Model DNA (4)
  if (4 %in% models_run) {
    VMN_ind.net.est <- ergm(VMN_ind.net ~ sum + nonzero + nodeocov("sequence") + mutual("min") +
                              nodeocov("dna"),
                            reference = ~Poisson,
                            response = "events") #transitive, triangle, isolates, idegree1.5, gwesp, gwdsp, gwodegree, twopath
    x = summary(VMN_ind.net.est)
    
    VMN_est_ind.df = bind_rows(VMN_est_ind.df, data.frame(pid = pid,
                                                          param = x$coefficients %>% row.names(),
                                                          exp = 0,
                                                          est = x$coefficients$Estimate,
                                                          pval = x$coefficients$`Pr(>|z|)`,
                                                          model = 4))
    
  }
  print("########################")
  print("Model 4")
  print("########################")
  
  #Model homophily - subgroups (5)
  if (5 %in% models_run) {
    VMN_ind.net.vest <- ergm(VMN_ind.net ~ sum + nonzero + nodeocov("sequence") + mutual("min") +
                               nodematch("group", diff = TRUE, form = "sum"),
                             reference = ~Poisson,
                             response = "events")
    
    x = summary(VMN_ind.net.vest)
    
    VMN_est_ind.df = bind_rows(VMN_est_ind.df, data.frame(pid = pid,
                                                          param = x$coefficients %>% row.names(),
                                                          exp = 0,
                                                          est = x$coefficients[,"Estimate"],
                                                          stderr = x$coefficients[,"Std. Error"],
                                                          pval = x$coefficients[,"Pr(>|z|)"],
                                                          model = 5))
  }
  print("########################")
  print("Model 5")
  print("########################")
  
  VMN_est.df = bind_rows(VMN_est.df, VMN_est_ind.df)
  
}
