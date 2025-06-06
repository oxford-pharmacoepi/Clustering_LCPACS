library(shiny)
library(shinydashboard)
library(officer)
library(tidyverse)
library(magrittr)
library(fs)
library(flextable)
library(here)
library(DT)

# Load Data in Global Environment
attrition_files <- list.files(here::here("Results"), pattern = "attrition_longcovid_.*\\.csv", full.names = TRUE)
attrition_tables <- setNames(lapply(attrition_files, read_csv), tools::file_path_sans_ext(basename(attrition_files)))

char_files <- list.files(here::here("Results/Tables"), pattern = "Clust_.*_BaselineCharacteristics\\.rds", full.names = TRUE)
clust_cases_BaselineCharacteristics <- readRDS(char_files[1])
clust_cases_matched_BaselineCharacteristics <- readRDS(char_files[2])
clust_controls_BaselineCharacteristics <- readRDS(char_files[3])
clust_controls_matched_BaselineCharacteristics <- readRDS(char_files[4])
clust_random_BaselineCharacteristics <- readRDS(char_files[5])

baseline_file <- here::here("Results/Tables/BaselineCharacteristics_clustering.rds")
baseline_file_matched <- here::here("Results/Tables/BaselineCharacteristics_clustering_matched.rds")
baseline_table <- if (file.exists(baseline_file)) {
  readRDS(baseline_file)
} else {
  "Baseline Characteristics file not found."
}
baseline_table_matched <- if (file.exists(baseline_file_matched)) {
  readRDS(baseline_file_matched)
} else {
  "Baseline Characteristics for matched cohort file not found."
}

# Define UI
ui <- dashboardPage(
  dashboardHeader(title = "Clustering Results"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Introduction", tabName = "intro", icon = icon("home")),
      menuItem("Attrition", tabName = "attrition", icon = icon("arrow-down-wide-short")),
      menuItem("Characteristics", tabName = "characteristics", icon = icon("user-group"),
               menuSubItem("Baseline Characteristics", tabName = "char_BaselineCharacteristics"),
               menuSubItem("Baseline Characteristics Matched", tabName = "char_BaselineCharacteristics_matched"),
               menuSubItem("Clusters Cases", tabName = "charac_cases"),
               menuSubItem("Clusters Controls", tabName = "charac_controls"),
               menuSubItem("Clusters Random", tabName = "charac_random"),
               menuSubItem("Clusters Cases Matched", tabName = "charac_cases_matched"),
               menuSubItem("Clusters Controls Matched", tabName = "charac_controls_matched")),
      menuItem("Clustering", tabName = "clustering", icon = icon("network-wired"),
               menuSubItem("2 clusters", tabName = "clust_2"),
               menuSubItem("3 clusters", tabName = "clust_3"),
               menuSubItem("4 clusters", tabName = "clust_4"),
               menuSubItem("Information criteria", tabName = "ic"),
          #     menuSubItem("Network Ising clustering", tabName = "network_ising"),
               menuSubItem("Network Walktrap clustering", tabName = "network_walktrap"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(
        tabName = "intro",
        fluidRow(
          column(12,
                 h3("Clustering of patients with post-acute COVID-19-related symptoms"),
                 p("This project explores the clustering of patients with post-acute COVID-19 symptoms using self-reported data from the UK Biobank (UKBB). The analysis uses latent class analysis (LCA) to identify patterns of symptom co-occurrence and group patients accordingly. This is carried out for people with a previous SARS-CoV-2 infection ('true cases'), people without it ('controls') and a random sample of both ('random')."),
                 p("The goal is to gain insights into the different phenotypes of post-acute COVID-19, which could inform future interventions and clinical management, and particularly whether they seem to be related to COVID-19 or not."),
                 img(src = "oxford.jpeg", height = "150px", align = "left")  
          )
        )
      ),
      tabItem(
        tabName = "attrition",
        tabsetPanel(
          tabPanel("Attrition cases",  DTOutput("table_attrition_longcovid_cases")),
          tabPanel("Attrition controls",  DTOutput("table_attrition_longcovid_controls_clustering"))
        )
      ),
      tabItem(tabName = "charac_cases", tableOutput("table_charac_cases")),
      tabItem(tabName = "charac_controls", tableOutput("table_charac_controls")),
      tabItem(tabName = "charac_cases_matched", tableOutput("table_charac_cases_matched")),
      tabItem(tabName = "charac_controls_matched", tableOutput("table_charac_controls_matched")),
      tabItem(tabName = "charac_random", tableOutput("table_charac_random")),
      tabItem(tabName = "char_BaselineCharacteristics", tableOutput("table_BaselineCharacteristics")),
      tabItem(tabName = "char_BaselineCharacteristics_matched", tableOutput("table_BaselineCharacteristics_matched")),
      tabItem(
        tabName = "clust_2",
        tabsetPanel(
          tabPanel("Cluster cases", imageOutput("img_clust_2")),
          tabPanel("Cluster controls", imageOutput("img_clust_2_co")),
          tabPanel("Cluster random", imageOutput("img_clust_2_ra")),
          tabPanel("Cluster cases matched", imageOutput("img_clust_2_m")),
          tabPanel("Cluster controls matched", imageOutput("img_clust_2_co_m")),
        )
      ),
      tabItem(
        tabName = "clust_3",
        tabsetPanel(
          tabPanel("Cluster cases", imageOutput("img_clust_3")),
          tabPanel("Cluster controls", imageOutput("img_clust_3_co")),
          tabPanel("Cluster random", imageOutput("img_clust_3_ra")),
          tabPanel("Cluster cases matched", imageOutput("img_clust_3_m")),
          tabPanel("Cluster controls matched", imageOutput("img_clust_3_co_m")),
        )
      ),
      tabItem(
        tabName = "clust_4",
        tabsetPanel(
          tabPanel("Cluster cases", imageOutput("img_clust_4")),
          tabPanel("Cluster controls", imageOutput("img_clust_4_co")),
          tabPanel("Cluster random", imageOutput("img_clust_4_ra")),
          tabPanel("Cluster cases matched", imageOutput("img_clust_4_m")),
          tabPanel("Cluster controls matched", imageOutput("img_clust_4_co_m")),
        )
      ),
      tabItem(
        tabName = "ic",
        tabsetPanel(
          tabPanel("IC cases", imageOutput("img_info_criteria")),
          tabPanel("IC controls", imageOutput("img_info_criteria_co")),
          tabPanel("IC random", imageOutput("img_info_criteria_ra")),
          tabPanel("IC cases matched", imageOutput("img_info_criteria_m")),
          tabPanel("IC controls matched", imageOutput("img_info_criteria_co_m")),
        )
      ),
      # tabItem(
      #   tabName = "network_ising",
      #   tabsetPanel(
      #     tabPanel("Ising cases", imageOutput("img_network_ising")),
      #     tabPanel("Ising controls", imageOutput("img_network_ising_co")),
      #     tabPanel("Ising random", imageOutput("img_network_ising_ra"))
      #   )
      # ),
      tabItem(
        tabName = "network_walktrap",
        tabsetPanel(
          tabPanel("Walktrap cases", imageOutput("img_network_walktrap")),
          tabPanel("Walktrap controls", imageOutput("img_network_walktrap_co")),
          tabPanel("Walktrap random", imageOutput("img_network_walktrap_ra")),
          tabPanel("Walktrap cases matched", imageOutput("img_network_walktrap_m")),
          tabPanel("Walktrap controls matched", imageOutput("img_network_walktrap_co_m"))
        )
      )
    )
  )
)

# Define Server
server <- function(input, output, session) {
  # Render tables for attrition_tables
  output$table_attrition_longcovid_cases <- DT::renderDataTable({
    attrition_tables$attrition_longcovid_cases |>
      dplyr::select(-"...1")
  })
  
  output$table_attrition_longcovid_controls_clustering <- DT::renderDataTable({
    attrition_tables$attrition_longcovid_controls_clustering |>
      dplyr::select(-"...1")
  })
  
  # Render tables for characteristics_tables
  output$table_charac_cases <- renderUI({
    clust_cases_BaselineCharacteristics |>
      htmltools_value()
  })
  
  output$table_charac_controls <- renderUI({
    clust_controls_BaselineCharacteristics |>
      htmltools_value()
  })
  
  output$table_charac_random <- renderUI({
    clust_random_BaselineCharacteristics |>
      htmltools_value()
  })
  
  output$table_charac_cases_matched <- renderUI({
    clust_cases_matched_BaselineCharacteristics |>
      htmltools_value()
  })
  
  output$table_charac_controls_matched <- renderUI({
    clust_controls_matched_BaselineCharacteristics |>
      htmltools_value()
  })

  # Render BaselineCharacteristics table
  output$table_BaselineCharacteristics <- renderUI({
    baseline_table %>%
      htmltools_value() 
  })
  
  output$table_BaselineCharacteristics_matched <- renderUI({
    baseline_table_matched %>%
      htmltools_value() 
  })
  
  # Render images
  render_image <- function(img_path, width = "100%") {
    if (file.exists(img_path)) {
      return(renderImage({ list(src = img_path, contentType = "image/jpeg", width = width) }, deleteFile = FALSE))
    } else {
      return(renderText("Image not found"))
    }
  }
  
  output$img_clust_2 <- render_image(here::here("Results/Clustering_cases/2_clust_1_symp/Clustering_LCA_clust_2_figure.jpg"), width = "50%")
  output$img_clust_2_co <- render_image(here::here("Results/Clustering_controls/2_clust_1_symp/Clustering_LCA_clust_2_figure.jpg"), width = "50%")
  output$img_clust_2_ra <- render_image(here::here("Results/Clustering_random/2_clust_1_symp/Clustering_LCA_clust_2_figure.jpg"), width = "50%")
  output$img_clust_2_m <- render_image(here::here("Results/Clustering_cases_matched/2_clust_1_symp/Clustering_LCA_clust_2_figure.jpg"), width = "50%")
  output$img_clust_2_co_m <- render_image(here::here("Results/Clustering_controls_matched/2_clust_1_symp/Clustering_LCA_clust_2_figure.jpg"), width = "50%")
  
  output$img_clust_3 <- render_image(here::here("Results/Clustering_cases/3_clust_1_symp/Clustering_LCA_clust_3_figure.jpg"), width = "50%")
  output$img_clust_3_co <- render_image(here::here("Results/Clustering_controls/3_clust_1_symp/Clustering_LCA_clust_3_figure.jpg"), width = "50%")
  output$img_clust_3_ra <- render_image(here::here("Results/Clustering_random/3_clust_1_symp/Clustering_LCA_clust_3_figure.jpg"), width = "50%")
  output$img_clust_3_m <- render_image(here::here("Results/Clustering_cases_matched/3_clust_1_symp/Clustering_LCA_clust_3_figure.jpg"), width = "50%")
  output$img_clust_3_co_m <- render_image(here::here("Results/Clustering_controls_matched/3_clust_1_symp/Clustering_LCA_clust_3_figure.jpg"), width = "50%")
  
  output$img_clust_4 <- render_image(here::here("Results/Clustering_cases/4_clust_1_symp/Clustering_LCA_clust_4_figure.jpg"), width = "50%")
  output$img_clust_4_co <- render_image(here::here("Results/Clustering_controls/4_clust_1_symp/Clustering_LCA_clust_4_figure.jpg"), width = "50%")
  output$img_clust_4_ra <- render_image(here::here("Results/Clustering_random/4_clust_1_symp/Clustering_LCA_clust_4_figure.jpg"), width = "50%")
  output$img_clust_4_m <- render_image(here::here("Results/Clustering_cases_matched/4_clust_1_symp/Clustering_LCA_clust_4_figure.jpg"), width = "50%")
  output$img_clust_4_co_m <- render_image(here::here("Results/Clustering_controls_matched/4_clust_1_symp/Clustering_LCA_clust_4_figure.jpg"), width = "50%")
  
#  output$img_network_ising <- render_image(here::here("Results/Clustering_cases/plot_isinggraph.png"), width = "50%")
#  output$img_network_ising_co <- render_image(here::here("Results/Clustering_controls/plot_isinggraph.png"), width = "50%")
#  output$img_network_ising_ra <- render_image(here::here("Results/Clustering_random/plot_isinggraph.png"), width = "50%")
  
  output$img_network_walktrap <- render_image(here::here("Results/Clustering_cases/Walktrap.png"), width = "50%")
  output$img_network_walktrap_co <- render_image(here::here("Results/Clustering_controls/Walktrap.png"), width = "50%")
  output$img_network_walktrap_ra <- render_image(here::here("Results/Clustering_random/Walktrap.png"), width = "50%")
  output$img_network_walktrap_m <- render_image(here::here("Results/Clustering_cases_matched/Walktrap.png"), width = "50%")
  output$img_network_walktrap_co_m <- render_image(here::here("Results/Clustering_controls_matched/Walktrap.png"), width = "50%")
  
  output$img_info_criteria <- render_image(here::here("Results/Clustering_cases/Clustering_LCA_symp_1_IC.jpg"), width = "50%")
  output$img_info_criteria_co <- render_image(here::here("Results/Clustering_controls/Clustering_LCA_symp_1_IC.jpg"), width = "50%")
  output$img_info_criteria_ra <- render_image(here::here("Results/Clustering_random/Clustering_LCA_symp_1_IC.jpg"), width = "50%")
  output$img_info_criteria_m <- render_image(here::here("Results/Clustering_cases_matched/Clustering_LCA_symp_1_IC.jpg"), width = "50%")
  output$img_info_criteria_co_m <- render_image(here::here("Results/Clustering_controls_matched/Clustering_LCA_symp_1_IC.jpg"), width = "50%")
  
}


# Run App
shinyApp(ui, server)
