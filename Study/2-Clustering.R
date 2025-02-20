# On to clustering
longCovid_cases <- read.csv(paste0(dir_results,"/longcovid_cases_clustering.csv")) |>
  as_tibble() |>
  dplyr::select(-"X")
longCovid_controls <- read.csv(paste0(dir_results,"/longcovid_controls_clustering.csv")) |>
  as_tibble() |>
  dplyr::select(-"X")
longCovid_random <- read.csv(paste0(dir_results,"/longcovid_random_clustering.csv")) |>
  as_tibble() |>
  dplyr::select(-"X")

names_symptoms <- colnames(longCovid_cases %>%
                             dplyr::select(-c(eid, specdate, questionnaire_started, symptom, length)))

# Load characteristics data ----
baselineCharacteristics <- as_tibble(read.csv(paste0(dir_results,"/CleanData_baselineCharacteristics.csv"))) |> dplyr::select(-c("X"))
biomarkers <- as_tibble(read.csv(paste0(dir_results,"/CleanData_biomarkers.csv")))  |> dplyr::select(-c("X"))
comorbidities <- as_tibble(read.csv(paste0(dir_results,"/CleanData_comorbidities.csv")))  |> dplyr::select(-c("X"))

run_clustering <- function(mydata,numclust, numsymp, counter, namefolder, results, nameclust) {
  output_clustering_w <- file.path(namefolder,paste0(numclust,"_clust_",numsymp,"_symp"))
  if (!file.exists(output_clustering_w)){
    dir.create(output_clustering_w, recursive = TRUE)}
  
  # Get people with more than the required number of clusters
  working_data <- mydata %>%
    dplyr::select(-"specdate") %>%
    dplyr::compute()
  
  x_vars <- c("1")
  cols <- paste0("cbind(", paste(names_symptoms, collapse = ","), ")")
  f <- with(working_data, as.formula(sprintf("%s ~ %s", cols, paste(x_vars, collapse = " + "))))
  
  entropy <- function(p) sum(-p*log(p))
  
  if(working_data %>% tally() > 100) {
    lc <- poLCA(f, working_data, nclass=numclust, maxiter=5000, graphs = FALSE,
                tol=1e-5, na.rm=FALSE,
                nrep=10, verbose=TRUE, calc.se=FALSE)
    
    # nrep, maxiter or the number of classes can be tuned later on if for some databases convergence or IC output does not look right
    
    # Create table with information criteria for the model
    error_prior <- entropy(lc$P)
    error_post <- mean(apply(lc$posterior,1, entropy),na.rm = TRUE)
    
    df <- lc$posterior
    df[] <- t(apply(df, 1, function(x) replace(x, x != max(x, na.rm = TRUE), NA)))
    numerator <- apply(df,2,sum,na.rm=T)
    denominator <- ifelse(is.na(df),0,1)
    denominator <- apply(denominator,2,sum,na.rm=T)
    
    results[[counter]] <- tibble::tibble(Clust = numclust,
                                         Symp = numsymp,
                                         log_likelihood= lc$llik,
                                         df = lc$resid.df,
                                         BIC= lc$bic,
                                         ABIC=  (-2*lc$llik) + ((log((lc$N + 2)/24)) * lc$npar),
                                         CAIC = (-2*lc$llik) + lc$npar * (1 + log(lc$N)),
                                         entropy = round(((error_prior-error_post) / error_prior),3),
                                         mean_posterior = numerator/denominator*100)
    
    
    write.csv(
      lc[["probs"]],
      file = here::here(output_clustering_w, paste0("Clustering_LCA_clust_",numclust,"_probs.csv")),
      row.names = FALSE
    )
    write.csv(
      lc[["P"]],
      file = here::here(output_clustering_w, paste0("Clustering_LCA_clust_",numclust,"_pclust.csv")),
      row.names = FALSE
    )
    write.csv(
      lc[["N"]],
      file = here::here(output_clustering_w, paste0("Clustering_LCA_clust_",numclust,"_N.csv")),
      row.names = FALSE
    )
    
    lcmodel <- reshape2::melt(lc$probs, level=2)
    lcmodel$L2 <- stringr::str_to_title(lcmodel$L2)
    
    # Change cluster assignment number depending on size
    cluster_sizes <- lc$predclass %>%
      table() %>%
      order(decreasing = TRUE)  # Sort by size (largest first)
    
    # Create a mapping from old to new cluster labels
    new_labels <- setNames(cluster_sizes, 1:numclust)
    
    # Reassign clusters based on size ranking
    if(numclust == 2) {
      lcmodel <- lcmodel %>%
        dplyr::mutate(
          Var1 = dplyr::case_when(
            Var1 == "class 1: " ~ paste0("class ", new_labels[[1]],":"),
            Var1 == "class 2: " ~ paste0("class ", new_labels[[2]],":"),
            TRUE ~ Var1
          )
        )
    } else if(numclust == 3){
      lcmodel <- lcmodel %>%
        dplyr::mutate(
          Var1 = dplyr::case_when(
            Var1 == "class 1: " ~ paste0("class ", new_labels[[1]],":"),
            Var1 == "class 2: " ~ paste0("class ", new_labels[[2]],":"),
            Var1 == "class 3: " ~ paste0("class ", new_labels[[3]],":"),
            TRUE ~ Var1
          )
        )
    } else {
      lcmodel <- lcmodel %>%
        dplyr::mutate(
          Var1 = dplyr::case_when(
            Var1 == "class 1: " ~ paste0("class ", new_labels[[1]],":"),
            Var1 == "class 2: " ~ paste0("class ", new_labels[[2]],":"),
            Var1 == "class 3: " ~ paste0("class ", new_labels[[3]],":"),
            Var1 == "class 4: " ~ paste0("class ", new_labels[[4]],":"),
            TRUE ~ Var1
          )
        )
    }
    
    # Get nice plot of class membership
    zp1 <- ggplot(lcmodel,aes(x = L2, y = value, fill = Var2))
    zp1 <- zp1 + geom_bar(stat = "identity", position = "stack")
    zp1 <- zp1 + facet_grid(Var1 ~ .)
    zp1 <- zp1 + scale_fill_brewer(type="seq", palette="Blues",labels = c("No symptom", "Symptom")) +theme_bw()
    zp1 <- zp1 + labs(x = "Symptoms",y=paste0("Prevalence symptoms in UKBB"), fill ="Response categories")
    zp1 <- zp1 + theme( axis.text.y=element_blank(),
                        axis.ticks.y=element_blank(),
                        panel.grid.major.y=element_blank(),
                        axis.text.x=element_text(size=10, angle = 90, vjust = 0.5, hjust=1))
    zp1 <- zp1 + guides(fill = guide_legend(reverse=TRUE))
    
    # Include horizontal lines with Prevalence of each symptom in all patients
    # List prevalence values (range from 0 to 1)
    # Transform back (1,2) to (0,1) binary
    working_data[,2:(ncol(working_data))] <- working_data[,2:(ncol(working_data))] -1
    factors <- (colSums(working_data[,2:(ncol(working_data))], na.rm = TRUE))/(length(working_data[[1]]))
    factors <- as.matrix(factors)
    rownames(factors) <- names_symptoms
    factors <- factors[order(rownames(factors)),]
    
    write.csv(
      factors,
      file = here::here(output_clustering_w, paste0("Clustering_LCA_clust_",numclust,"_average.csv")),
      row.names = TRUE
    )
    
    factors <- factors[factors != 0]
    
    # Include prevalence average lines
    for (i in 1:length(factors)) {
      zp1 <- zp1 + geom_segment(x = (i - 0.5), y = factors[[i]], xend = (i + 0.5), yend = factors[[i]])
    }
    
    ggsave(here::here(output_clustering_w, paste0("Clustering_LCA_clust_",numclust,"_figure.jpg")),
           height = 8, width = 8)
    
    # Characterise clusters
    # Look at characterisation of clusters: age and sex
    working_data <- working_data %>%
      mutate(cluster_assignment = lc$predclass) %>%
      compute()
    
    # Change cluster assignment number depending on size
    cluster_sizes <- working_data %>%
      count(cluster_assignment) %>%
      arrange(desc(n))  # Sort by size (largest first)
    
    # Create a mapping from old to new cluster labels
    new_labels <- setNames(cluster_sizes$cluster_assignment, 1:numclust)
    
    # Reassign clusters based on size ranking
    working_data <- working_data %>%
      mutate(new_cluster = recode(cluster_assignment, 
                                  !!!setNames(1:numclust, new_labels))) |>
      dplyr::select(-"cluster_assignment") |>
      rename("cluster_assignment" = "new_cluster")
    
    # Look at number of people with symptom per cluster
    number_people <- working_data %>% dplyr::select(-c("eid"))
    clust <- list()
    for(i in 1:numclust) {
      clust[[i]] <- apply(number_people %>% dplyr::filter(cluster_assignment == i) %>%
                            dplyr::select(-cluster_assignment), 2, sum, na.rm = TRUE) |>
        as_tibble() |>
        dplyr::mutate(cluster = i)
      clust[[i]] <-  clust[[i]] |>
        dplyr::mutate(symptom = apply(number_people %>% dplyr::filter(cluster_assignment == i) %>%
                                        dplyr::select(-cluster_assignment), 2, sum, na.rm = TRUE) |>
                        names())
    }
    
    clust <- dplyr::bind_rows(clust) |>
      dplyr::mutate(cluster = paste0("cluster_",cluster)) |>
      tidyr::pivot_wider(names_from = "cluster",
                         values_from = "value")
    
    write.csv(
      clust,
      file = here::here(output_clustering_w, paste0("Clustering_LCA_clust_",numclust,"_n_each_symp.csv")),
      row.names = TRUE
    )
    
    num_symp_total <- working_data |>
      rowwise() |>
      mutate("symptom" = across(starts_with("symptom")) |> sum(na.rm = TRUE)) |>
      ungroup() |>
      dplyr::group_by(cluster_assignment, symptom) |>
      dplyr::summarise(n = n()) |>
      dplyr::mutate(cluster_assignment = paste0("n_clust_",cluster_assignment)) |>
      tidyr::pivot_wider(names_from = "cluster_assignment",
                         values_from = "n")
    
    write.csv(
      num_symp_total,
      file = here::here(output_clustering_w, paste0("Clustering_LCA_clust_",numclust,"_n_symp_total.csv")),
      row.names = TRUE
    )
    
    num_people_cluster <- working_data |>
      dplyr::group_by(cluster_assignment) |>
      dplyr::summarise(n = n())
    
    write.csv(
      num_people_cluster,
      file = here::here(output_clustering_w, paste0("Clustering_LCA_clust_",numclust,"_n_per_cluster.csv")),
      row.names = TRUE
    )
    
    working_data <- working_data %>%
      left_join(
        mydata |>
          dplyr::select(c("eid", "specdate")),
        by = "eid"
      ) %>%
      compute()
    
    # Baseline Characteristics and biomarkers if nclust is 3
    
    if(numclust == 3){
      # Long covid clustering
      longCovid_clust <- tableOneStep1(
        working_data |> rename("state" = "cluster_assignment") |>
          mutate(state = state - 1), 
        baselineCharacteristics, biomarkers, "Long COVID cluster", c("Cluster 1", "Cluster 2", "Cluster 3"))
      
      # Merge all tables ----
      y <- longCovid_clust |>
        mutate(order = row_number()) |>
        dplyr::select(-c("order")) |>
        mutate(`Risk factor` = if_else(`Risk factor` == "Sociodemographics factors", "Sociodemographic factors",`Risk factor`))
      
      y <- y |>
        filter(`Risk factor` == "N") |>
        rbind(
          y |>
            filter(`Risk factor` != "N")
        )
      
      merge <- c(1:nrow(y))[y$`Risk factor` %in% c("Sociodemographic factors", "Comorbidities [Cases (%)]", "Biomarkers [Mean (SD)]")]
      row_fill <- c(1:nrow(y))[!y$`Long COVID cluster_Cluster 1` == " "]
      y <- y |>
        flextable() |>
        span_header() |>
        bold(i = 1, part = "header") |>
        bold(i = 2, part = "header") |>
        align(j = c(2,3,4), align = "center", part = "all") |>
        width(j = 1, width = 5, unit = "cm") |>
        bg(bg = "#F2F2F2", i = 1, j = c(2,3), part = "header") |>
        bg(bg = "#F2F2F2", i = 2, j = c(2,4), part = "header") |>
        bg(bg = "#F2F2F2", i = row_fill, j = c(2,4), part = "body")
      
      for(ii in merge){
        y <- y |>
          merge_at(i = ii, j = c(1:4)) |>
          hline(i = ii-1, border = fp_border(color = "black")) |>
          bold(i = ii)
      }
      
      y |>
        save_as_docx(path = paste0(dir_results, "/Tables/Clust_",nameclust,"_BaselineCharacteristics.docx"))
      
      saveRDS(y, paste0(dir_results, "/Tables/Clust_",nameclust,"_BaselineCharacteristics.rds"))
      
      rm(list = c("baselineCharacteristics", "biomarkers","x","name","name_cohort",
                  "longCovid_cohort","pacs_cohort","y", "merge", "row_fill"))
    }
  }
  return(results)
}

do_clustering <- function(cohort, foldername, nameclust) {
  mydata <- cohort %>%
    dplyr::select(-c(questionnaire_started, symptom, length))
  
  mydata[, 3:ncol(mydata)] <- mydata[, 3:ncol(mydata)] + 1  # Adjust for poLCA
  
  results <- list()
  counter <- 1
  for(nc in c(2:4)) { # 2 to 4 clusters
    #for(ns in c(1:3)) { # 1 to 3 number of symptoms
    results <- run_clustering(mydata,nc,1,counter,foldername, results, nameclust)
    counter <- counter + 1
    #}
  }
  results <- bind_rows(results)
  results2 <- tidyr::gather(results,Criteria,Value,5:9)
  
  # Save table with all IC information
  write.csv(
    results,
    file = here::here(foldername, paste0("Information_criteria.csv")),
    row.names = FALSE
  )
  # Plots with information criteria results for all models tried (one per # symptoms)
  result_working <- results2 %>% dplyr::filter(Symp == 1)
  
  if(result_working %>% dplyr::tally() %>% dplyr::pull() != 0) {
    fit.plot <- ggplot(result_working) +
      geom_point(aes(x=Clust,y=Value),size=3) +
      geom_line(aes(Clust, Value, group = 1)) +
      theme_bw()+
      labs(x = "", y="", title = "") +
      facet_grid(Criteria ~. ,scales = "free") +
      theme_bw(base_size = 16, base_family = "") +
      theme(panel.grid.major.x = element_blank() ,
            panel.grid.major.y = element_line(colour="grey", size=0.5),
            legend.title = element_text(size = 16, face = 'bold'),
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 16),
            legend.text=  element_text(size=16),
            axis.line = element_line(colour = "black"))
    
    ggsave(here::here(foldername, paste0("Clustering_LCA_symp_","1","_IC.jpg")))
    
  }
  
  network_data <- cohort %>%
    dplyr::select(-c(eid,specdate, questionnaire_started, symptom, length)) %>%
    mutate(across(starts_with("symptom"), ~ replace(., . == -1, 0))) %>%
    as.matrix()
  
  colnames(network_data) <- gsub("symptom_", "", colnames(network_data))
  network_data[is.na(network_data)] <- 0
  
  graph <- IsingFit(network_data)
  png(here::here(foldername,"plot_isinggraph.png"))
  plot(graph)
  dev.off()
  adj_matrix <- graph$weiadj
  
  # Adjust weights and create igraph object
  if (min(adj_matrix) < 0) adj_matrix <- adj_matrix - min(adj_matrix)
  igraph_obj <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)
  
  # Community detection
  clusters <- cluster_louvain(igraph_obj)
  plot(igraph_obj, vertex.color = clusters$membership)
  membership(clusters)
  write.csv(membership(clusters), file = here::here(foldername, "Louvain.csv"))
  
  clusters <- cluster_walktrap(igraph_obj)
  plot(igraph_obj, vertex.color = clusters$membership)
  plot(as.dendrogram(clusters))
  png(here::here(foldername,"Walktrap.png"), height = 10, width = 10, units = "in", res = 100)
  par(mar=c(3,1,1,19))
  plot(as.dendrogram(clusters), horiz = TRUE)
  dev.off()
  membership(clusters)
  write.csv(membership(clusters), file = here::here(foldername, "Walktrap.csv"))
}

do_clustering(longCovid_cases, paste0(dir_results,"/Clustering_cases"), "cases")
do_clustering(longCovid_controls, paste0(dir_results,"/Clustering_controls"), "controls")
do_clustering(longCovid_random, paste0(dir_results,"/Clustering_random"), "random")

# Table ones for the whole clustering cohorts (cases and controls)

# Long covid clustering
x <- longCovid_cases |>
  mutate(state =1) |>
  union_all(
    longCovid_controls |>
      mutate(state = 0)
  )
name <- "Long COVID clustering"
name_cohort <- c("Controls", "Cases")
longCovid_clust <- tableOneStep1(x, baselineCharacteristics, biomarkers, name, name_cohort)

# Merge all tables ----
y <- longCovid_clust |>
  mutate(order = row_number()) |>
  dplyr::select(-c("order")) |>
  mutate(`Risk factor` = if_else(`Risk factor` == "Sociodemographics factors", "Sociodemographic factors",`Risk factor`))

y <- y |>
  filter(`Risk factor` == "N") |>
  rbind(
    y |>
      filter(`Risk factor` != "N")
  )

merge <- c(1:nrow(y))[y$`Risk factor` %in% c("Sociodemographic factors", "Comorbidities [Cases (%)]", "Biomarkers [Mean (SD)]")]
row_fill <- c(1:nrow(y))[!y$`Long COVID clustering_Controls` == " "]
y <- y |>
  flextable() |>
  span_header() |>
  bold(i = 1, part = "header") |>
  bold(i = 2, part = "header") |>
  align(j = c(2,3), align = "center", part = "all") |>
  width(j = 1, width = 5, unit = "cm") |>
  bg(bg = "#F2F2F2", i = 1, j = c(2,3), part = "header") 

for(ii in merge){
  y <- y |>
    merge_at(i = ii, j = c(1:3)) |>
    hline(i = ii-1, border = fp_border(color = "black")) |>
    bold(i = ii)
}

y |>
  save_as_docx(path = paste0(dir_results, "/Tables/BaselineCharacteristics_clustering.docx"))

saveRDS(y, paste0(dir_results,"/Tables/BaselineCharacteristics_clustering.rds"))

rm(list = c("baselineCharacteristics", "biomarkers","x","name","name_cohort",
            "longCovid_cohort","pacs_cohort","y", "merge", "row_fill"))
