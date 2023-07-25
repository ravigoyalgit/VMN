

VMN_create_missing <- function(VMN.df, Tissue.df, pid_ind, 
                               missing_edges = TRUE,
                               VMN_wgt_metric = "log_totalCount",
                               TCR_wgt_metric = "numbersharedclones") {
  
  VMN_ind.df = VMN.df %>% filter(pid == pid_ind) %>% 
    mutate(across(where(is.character), toupper)) %>%
    mutate(across(where(is.character), str_remove_all, pattern = fixed(" ")))
  
  VMN_ind.df$from = recode(VMN_ind.df$from, 
                           `LYMPHNODES(PERITRACHEAL)` = "LNPERITRACHEAL",
                           `LYMPHNODES(OTHERS)` = "LNOTHERS",
                           `LYMPHNODES(AXILLARY)` = "LNAXILLARY",
                           `LYMPHNODES(AORTIC)` = "LNAORTIC",
                           `LYMPHNODES(INGUINAL)` = "LNINGUINAL",
                           `LYMPHNODES(MESENTARY)` = "LNMESENTARY",
                           `LYMPHNODES(MEDIASTINAL)` = "LNMEDIASTINAL",
                           `LYMPHNODES(HILAR)` = "LNHILAR"
  )
  
  VMN_ind.df$to = recode(VMN_ind.df$to, 
                         `LYMPHNODES(PERITRACHEAL)` = "LNPERITRACHEAL",
                         `LYMPHNODES(OTHERS)` = "LNOTHERS",
                         `LYMPHNODES(AXILLARY)` = "LNAXILLARY",
                         `LYMPHNODES(AORTIC)` = "LNAORTIC",
                         `LYMPHNODES(INGUINAL)` = "LNINGUINAL",
                         `LYMPHNODES(MESENTARY)` = "LNMESENTARY",
                         `LYMPHNODES(MEDIASTINAL)` = "LNMEDIASTINAL",
                         `LYMPHNODES(HILAR)` = "LNHILAR"
  )
  
  Tissue_mod.df = Tissue.df %>% 
    mutate(across(where(is.character), toupper))
  
  Tissue_mod.df$location = recode(Tissue_mod.df$location,
                                 COLONRIGHT = "COLON"
  )
  
  Tissue_mod.df = Tissue_mod.df %>% mutate(loc_fac = as.factor(location),
                                           loc_int = as.numeric(loc_fac)) 
  
  Tissue_mod.df = Tissue_mod.df %>% mutate_at(vars(dna), ~replace(., is.na(.), 0))
  
  Tissue_ind.df = Tissue_mod.df %>% 
    group_by(location) %>% 
    slice_head() %>% 
    select(location, Groups, loc_fac, loc_int)
  
  Tissue_ind.df = left_join(Tissue_ind.df, Tissue_mod.df %>% 
                              filter(pid == pid_ind) %>%
                              select(location, dna, sequences) %>%
                              mutate(obs = 1),
                            by="location")
  
  VMN_ind.df = left_join(VMN_ind.df, Tissue_ind.df %>% select(c(location, loc_int)), by = c("from" = "location")) %>%
    rename(from_int = loc_int)
  
  VMN_ind.df = left_join(VMN_ind.df, Tissue_ind.df %>% select(c(location, loc_int)), by = c("to" = "location")) %>%
    rename(to_int = loc_int)
  
  if ((sum(is.na(VMN_ind.df$from_int)) + sum(is.na(VMN_ind.df$to_int))) > 0) {
    print("Error: Missing tissue name")
  }
  
  #add.edges(VMN_ind.net, tail = VMN_ind.df$from_int, head = VMN_ind.df$to_int)
  tissue_missing = Tissue_ind.df %>% filter(is.na(obs)) %>% pull(loc_int)
  
  VMN_ind_missing.df = crossing(from_int = Tissue_ind.df %>% pull(loc_int), to_int = Tissue_ind.df %>% pull(loc_int)) %>%
    mutate(from_int_missing = from_int %in% tissue_missing,
           to_int_missing = to_int %in% tissue_missing,
           non_miss_edge = from_int_missing == FALSE & to_int_missing == FALSE) %>%
    filter(non_miss_edge == FALSE) %>%
    filter(from_int != to_int) %>%
    select(from_int, to_int) %>%
    mutate(pid = pid_ind,
           from = NA,
           to = NA,
           totalCount = 2)
  
  if (missing_edges) {
    VMN_ind.df = bind_rows(VMN_ind.df, VMN_ind_missing.df)
  }
  
  VMN_ind.net = network::network.initialize(
    n = max(Tissue_mod.df$loc_int),
    directed = TRUE,
    hyper = FALSE,
    loops = FALSE,
    multiple = FALSE,
    bipartite = FALSE
  )
  
  if (VMN_wgt_metric == "log_totalCount") {
    VMN_ind.net[VMN_ind.df[,5:6], names.eval="events", add.edges=TRUE] <- VMN_ind.df %>% pull(totalCount) %>% log() %>% round()
  }
  if (VMN_wgt_metric == "totalCount") {
    VMN_ind.net[VMN_ind.df[,5:6], names.eval="events", add.edges=TRUE] <- VMN_ind.df %>% pull(totalCount) 
  }
  if (VMN_wgt_metric == "totalCount_d_sequences") {
    
    events_TEMP = left_join(VMN_ind.df, Tissue_ind.df, by = c("from" = "location")) %>%
      mutate(events = totalCount/sequences  %>% round()) %>% pull(events)
    
    VMN_ind.net[VMN_ind.df[,5:6], names.eval="events", add.edges=TRUE] <- events_TEMP
  }
  
  if (missing_edges) {  
    for (i in c(1:nrow(VMN_ind_missing.df))) { 
      VMN_ind.net[VMN_ind_missing.df$from_int[i],VMN_ind_missing.df$to_int[i]] <- NA
    }
    missing_tissue.df = tibble(tissue_id = c(1:max(Tissue_mod.df$loc_int)),
                               missing = 0)
    missing_tissue.df$missing[tissue_missing] = 1
  }
  
  Tissue_ind.df$dna[is.na(Tissue_ind.df$dna)] <- mean(Tissue_ind.df$dna, na.rm = TRUE)
  Tissue_ind.df$sequences[is.na(Tissue_ind.df$sequences)] <- mean(Tissue_ind.df$sequences, na.rm = TRUE)
  
  set.vertex.attribute(VMN_ind.net,"group",Tissue_ind.df$Groups)
  set.vertex.attribute(VMN_ind.net,"dna",Tissue_ind.df$dna)
  set.vertex.attribute(VMN_ind.net,"sequence",Tissue_ind.df$sequences)
  
  if (missing_edges) { 
    set.vertex.attribute(VMN_ind.net,"missing",missing_tissue.df$missing)
  }
  
  return(VMN_ind.net)

}

