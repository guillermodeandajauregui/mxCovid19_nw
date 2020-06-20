###############################################################################
#Spatial correlations of COVID-19 epidemic curves
#Mexican municipios
#gdeanda@inmegen.edu.mx
###############################################################################

###############################################################################
#libraries
###############################################################################

library(tidyverse)
library(vroom)
library(lubridate)
library(infotheo)
library(parallel)
library(igraph)
library(tidygraph)
library(pheatmap)
###############################################################################
#read data
###############################################################################

############################
###today 
############################
el_path <- "http://187.191.75.115/gobmx/salud/datos_abiertos/historicos/datos_abiertos_covid19_18.06.2020.zip"
temp <- tempfile()
download.file(url = el_path, destfile = temp)
el_file <- unzip(zipfile = temp, list = T)
el_file <- unz(description = temp, filename = el_file[1])
el_data <- vroom::vroom(el_file)
unlink(temp)
rm(temp, el_file)


today_data  = list(el_data) 
names(today_data) <-paste(format((Sys.Date()-1), "%Y%m%d"))

data_gobmx = el_data

###############################################################################
#get case counts per municipio per date
###############################################################################
my_data <-
  data_gobmx %>% 
  select(ID_REGISTRO, ENTIDAD_RES, MUNICIPIO_RES, TIPO_PACIENTE, FECHA_SINTOMAS, FECHA_DEF, RESULTADO)


mis_fechas = seq(from = ymd('2020-03-01'), to = ymd("2020-05-31"), by="days")


#confirmed counts 
daily_cases = 
  my_data %>% 
  filter(RESULTADO == 1) %>%
  #time period between march and june 
  filter(FECHA_SINTOMAS > "2020-03-01") %>%
  filter(FECHA_SINTOMAS < "2020-06-01") %>% 
  group_by(ENTIDAD_RES, MUNICIPIO_RES, FECHA_SINTOMAS) %>% 
  tally()
###############################################################################
#normalize per 100,000 people?
###############################################################################
#https://datos.gob.mx/busca/dataset/proyecciones-de-la-poblacion-de-mexico-y-de-las-entidades-federativas-2016-2050/resource/751728c1-e0cf-4fe8-b0fb-55b17d22bac4?inner_span=True

#intercensus data, parsed by Irving Morales 
pop_data = vroom::vroom("data/INTERCENSAL_pob_mun_2015_edades.csv")
pop_totals = pop_data %>% janitor::clean_names() %>% filter(grupos_quinquenales_de_edad=="Total")
pop_totals.light = pop_totals %>% select(cveedo, cvemun, poblacion_total, cvegeo)

#join and normalize
#MI is scale invariant so any normalization will do
daily_cases.norm <-
  left_join(daily_cases, pop_totals.light, by=c("ENTIDAD_RES" = "cveedo", "MUNICIPIO_RES" = "cvemun")) %>% 
  mutate(casos_fraq = n/poblacion_total)

###############################################################################
#fill missing dates for completeness
###############################################################################

mis_fechas.df = data.frame(date = mis_fechas)

#split and get keys
daily_cases.norm.split <- 
  daily_cases.norm %>% 
  group_by(cvegeo) %>% 
  group_split()


daily_cases.norm.split.nona <- 
  lapply(X = daily_cases.norm.split, function(i){
    
    my_completor = data.frame(FECHA_SINTOMAS = mis_fechas, 
                              ENTIDAD_RES    = i[["ENTIDAD_RES"]][1],
                              MUNICIPIO_RES  = i[["MUNICIPIO_RES"]][1],
                              cvegeo  = i[["cvegeo"]][1],
                              poblacion_total = i[["poblacion_total"]][1]
    )
    
    left_join(my_completor, i) %>% 
      mutate(casos_fraq  = ifelse(is.na(casos_fraq), 0, casos_fraq),
             n = ifelse(is.na(n), 0, n))
  }) %>% bind_rows()

#


###############################################################################
#make list of count vector
###############################################################################

cvegeo_ids = daily_cases.norm.split.nona$cvegeo %>% unique
names(cvegeo_ids) <- cvegeo_ids
## free up some memory
rm(daily_cases.norm.split)
## 

daily_cases.norm.split.nona %>% 
  filter(cvegeo == "01001") %>% 
  pull(casos_fraq)



count_vctrs <- 
  lapply(cvegeo_ids, function(i){
    daily_cases.norm.split.nona %>% 
      filter(cvegeo == i) %>% 
      pull(casos_fraq)
  })

###############################################################################
#discretize
###############################################################################

names(count_vctrs) <- paste0("mun_", names(count_vctrs))
which(names(count_vctrs)=="mun_NA")

count_vctrs = count_vctrs[-1624] #remove an na

count_vctrs.discrete <- lapply(X = count_vctrs, 
                               FUN = function(i){discretize(i)}
)

###############################################################################
#calculate mutual information
###############################################################################

tempus = Sys.time()
mi_list <- 
  mclapply(X = count_vctrs.discrete, 
           mc.cores = 40, 
           FUN = function(i){
             sapply(count_vctrs.discrete, FUN = function(j){
               infotheo::mutinformation(i,j, method = "mm")
             })
           })

tempus = Sys.time() - tempus
print(tempus)

###############################################################################
#work with matrix
###############################################################################

mi_mx   = matrix(unlist(mi_list), nrow = length(mi_list), byrow = T)

colnames(mi_mx) <- names(mi_list)
rownames(mi_mx) <- names(mi_list)


density(mi_mx) %>% plot

#normalize matrix and remove diagonal
mi_mx.clean = mi_mx
diag(mi_mx.clean) <- 0


mi_mx.clean = mi_mx.clean/max(mi_mx.clean)

mi_mx.clean[1:5, 1:5]

#exploratory plots
density(mi_mx.clean) %>% plot
pheatmap::pheatmap(mat = mi_mx.clean, main = "normalizado por población")
###############################################################################
#network construction scans
###############################################################################

my_thresholds = seq(from = 0.05, to = 0.95, by = 0.05)
names(my_thresholds) = my_thresholds

nw_list = 
  lapply(X = my_thresholds, function(i){
    
    mx = ifelse(mi_mx.clean <= i, 0, 1)
    
    g = graph_from_adjacency_matrix(adjmatrix = mx, mode = "undirected", weighted = NULL, diag = F)
    
    # data.frame(th = i,
    #            nodes = length(V(g)),
    #            edges = length(E(g)),
    #            density = edge_density(g),
    #            concomp = components(g)$no
    #            )
    # 
  }) #%>% bind_rows()

###############################################################################
#network construction scans
###############################################################################


gt = tidygraph::as_tbl_graph(nw_list$`0.55`)

gt <- 
  gt %>% 
  activate("nodes") %>% 
  mutate(my_louvain = group_louvain(),
         my_compone = group_components())

gt <-
gt %>%
  activate("nodes") %>%
  mutate(cvegeo = str_remove(name, "mun_")) %>%
  left_join(y = pop_totals) %>%
  mutate(edo_mun = paste0(str_sub(string = estado, start = 1, end = 4),
                          "_",
                          str_sub(string = nombre, start = 1, end = 8)
                          )
         )



write_rds(gt, "nw_055.rds")


