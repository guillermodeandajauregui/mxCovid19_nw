library(igraph)
library(tidygraph)

#read network, "gt"
#gt


#name a sequence of numbers for ease of access
my_seeds = 1:100
names(my_seeds) = my_seeds

#attributes to keep 
my_node_attr = names(vertex.attributes(gt))

my_keepers = c("name", "cvegeo", "cveedo", 
               "estado", "cvemun", "nombre", 
               "poblacion_total", "LAT", "LON")

my_deleters = setdiff(my_node_attr, my_keepers)
#model 1: basic rewiring with prob 0.5 (all nodes) ----


NullModell_01 <- 
  lapply(X = my_seeds, FUN = function(i){
    set.seed(i)
    g_null = rewire(graph = gt, with = each_edge(prob = 0.5, loops = F, multiple = F))
    delete_vertex_attr(g_null, name = "my_infomap")
    g_null = 
      g_null %>% 
      as_tbl_graph() %>% 
      activate("nodes") %>% 
      tidygraph::select(all_of(my_keepers)) 
  })

summary_nullmodel_01 <-
  lapply(X = NullModell_01, FUN = function(g){
    set.seed(725)
    data.frame(avg_k = mean(degree(g)),
               clus_coef  = transitivity(g, type = "global"),
               avg.path.length = average.path.length(g),
               no_conn_comp = components(g)$no,
               max_conn_comp_size  = max(components(g)$csize),
               max_louvain_module  = max(sizes(igraph::cluster_louvain(g))),
               no_louvain_module   = length(which(sizes(igraph::cluster_louvain(g))!=1))
    ) %>% as_tibble()
  }) %>% bind_rows() %>% summarize_all(mean, na.rm = T)

#model 2: basic rewiring with prob 0.5 (main component) ----

NullModell_02 <- 
  lapply(X = my_seeds, FUN = function(i){
    set.seed(i)
    
    #keep only main component
    gtb = gt %>% activate("nodes") %>% filter(my_compone ==1)
    
    g_null = rewire(graph = gtb, with = each_edge(prob = 0.5, loops = F, multiple = F))
    delete_vertex_attr(g_null, name = "my_infomap")
    g_null = 
      g_null %>% 
      as_tbl_graph() %>% 
      activate("nodes") %>% 
      tidygraph::select(all_of(my_keepers)) 
  })

summary_nullmodel_02 <-
  lapply(X = NullModell_02, FUN = function(g){
    set.seed(725)
    data.frame(avg_k = mean(degree(g)),
               clus_coef  = transitivity(g, type = "global"),
               avg.path.length = average.path.length(g),
               no_conn_comp = components(g)$no,
               max_conn_comp_size  = max(components(g)$csize),
               max_louvain_module  = max(sizes(igraph::cluster_louvain(g))),
               no_louvain_module   = length(which(sizes(igraph::cluster_louvain(g))!=1))
    ) %>% as_tibble()
  }) %>% bind_rows() %>% summarize_all(mean, na.rm = T)

#model 3: configuration rewiring (all nodes) ----

NullModell_03 <- 
  lapply(X = my_seeds, FUN = function(i){
    set.seed(i)
    g_null = rewire(graph = gt, with = keeping_degseq(loops = F, niter = 100))
    delete_vertex_attr(g_null, name = "my_infomap")
    g_null = 
      g_null %>% 
      as_tbl_graph() %>% 
      activate("nodes") %>% 
      tidygraph::select(all_of(my_keepers)) 
  })

summary_nullmodel_03 <-
  lapply(X = NullModell_03, FUN = function(g){
    set.seed(725)
    data.frame(avg_k = mean(degree(g)),
               clus_coef  = transitivity(g, type = "global"),
               avg.path.length = average.path.length(g),
               no_conn_comp = components(g)$no,
               max_conn_comp_size  = max(components(g)$csize),
               max_louvain_module  = max(sizes(igraph::cluster_louvain(g))),
               no_louvain_module   = length(which(sizes(igraph::cluster_louvain(g))!=1))
    ) %>% as_tibble()
  }) %>% bind_rows() %>% summarize_all(mean, na.rm = T)

#model 4: configuration rewiring (main component) ----

NullModell_04 <- 
  lapply(X = my_seeds, FUN = function(i){
    set.seed(i)
    
    gtb = gt %>% activate("nodes") %>% filter(my_compone ==1)
    
    g_null = rewire(graph = gtb, with = keeping_degseq(loops = F, niter = 100))
    delete_vertex_attr(g_null, name = "my_infomap")
    g_null = 
      g_null %>% 
      as_tbl_graph() %>% 
      activate("nodes") %>% 
      tidygraph::select(all_of(my_keepers)) 
  })

summary_nullmodel_04 <-
  lapply(X = NullModell_04, FUN = function(g){
    set.seed(725)
    data.frame(avg_k = mean(degree(g)),
               clus_coef  = transitivity(g, type = "global"),
               avg.path.length = average.path.length(g),
               no_conn_comp = components(g)$no,
               max_conn_comp_size  = max(components(g)$csize),
               max_louvain_module  = max(sizes(igraph::cluster_louvain(g))),
               no_louvain_module   = length(which(sizes(igraph::cluster_louvain(g))!=1))
               ) %>% as_tibble()
  }) %>% bind_rows() %>% summarize_all(mean, na.rm = T)

#for comparison, results of real network
results_full = 
  data.frame(avg_k = mean(degree(gt)),
             clus_coef  = transitivity(gt, type = "global"),
             avg.path.length = average.path.length(gt),
             no_conn_comp = components(gt)$no,
             max_conn_comp_size  = max(components(gt)$csize),
             max_louvain_module  = max(sizes(igraph::cluster_louvain(gt))),
             no_louvain_module   = length(which(sizes(igraph::cluster_louvain(gt))!=1))
  ) %>% as_tibble()

main_component = gt %>% activate("nodes") %>% filter(my_compone ==1)

results_main = 
  data.frame(avg_k = mean(degree(main_component)),
             clus_coef  = transitivity(main_component, type = "global"),
             avg.path.length = average.path.length(main_component),
             no_conn_comp = components(main_component)$no,
             max_conn_comp_size  = max(components(main_component)$csize),
             max_louvain_module  = max(sizes(igraph::cluster_louvain(main_component))),
             no_louvain_module   = length(which(sizes(igraph::cluster_louvain(main_component))!=1))
  ) %>% as_tibble()


#merge them 
models_with_all_nodes = 
  bind_rows(list(model_01  = summary_nullmodel_01,
                 model_03  = summary_nullmodel_03,
                 empirical = results_full), 
            .id = "model"
  ) %>% 
  pivot_longer(cols = -model) %>% 
  pivot_wider(names_from = name, values_from = value)

models_with_main_component = 
  bind_rows(list(model_02  = summary_nullmodel_02,
                 model_04  = summary_nullmodel_04,
                 empirical = results_main), 
            .id = "model"
  ) %>% 
  pivot_longer(cols = -model) %>% 
  pivot_wider(names_from = name, values_from = value)

models_with_all_nodes
models_with_main_component



