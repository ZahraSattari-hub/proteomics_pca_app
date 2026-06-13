
getwd()

library(shiny)
library(plotly)
library(readxl)
library(dplyr)
library(DT)



# --------------------
# DATA
# --------------------

df <- read_excel(
  "mer_proteins.xlsx",
  sheet = 4
)

sample_info <- df[,1:2]

protein_data0 <- df[,3:ncol(df)]

# protein_data <- log1p(protein_data)


col2remove <- which(apply(protein_data0, 2, sd) == 0)
View(col2remove)

# remove cols with no variation

protein_data <- protein_data0[, -col2remove]
View(protein_data)


# log transform
protein_data <- log1p(protein_data)
View(protein_data)



pca <- prcomp(
  protein_data,
  center = TRUE,
  scale. = TRUE
)

scores <- as.data.frame(pca$x)

scores$sample_nr <- sample_info$sample_nr

scores$treatment <- factor(
  sample_info$treatment,
  labels = c("Control","Treated")
)

# --------------------
# LOADINGS
# --------------------

loadings <- as.data.frame(pca$rotation)

loadings$protein <- rownames(loadings)


# --------------------
# DIFFERENTIAL ABUNDANCE
# --------------------

protein_df <- protein_data

protein_df$treatment <- sample_info$treatment

diff_results <- data.frame()

for(p in names(protein_data)){
  
  control_mean <-
    mean(protein_df[protein_df$treatment==0,p])
  
  treated_mean <-
    mean(protein_df[protein_df$treatment==1,p])
  
  diff_results <-
    rbind(
      diff_results,
      data.frame(
        protein = p,
        control_mean = control_mean,
        treated_mean = treated_mean,
        difference = treated_mean - control_mean
      )
    )
}

diff_results <-
  diff_results %>%
  arrange(desc(abs(difference)))



DTOutput("loading_table")

tabsetPanel(
  
  tabPanel(
    "PC Loadings",
    DTOutput("loading_table")
  ),
  
  tabPanel(
    "Treatment Differences",
    DTOutput("diff_table")
  )
  
)


# Add new server output

output$diff_table <- renderDT({
  
  diff_results %>%
    select(
      protein,
      control_mean,
      treated_mean,
      difference
    ) %>%
    head(input$nproteins)
  
})



# --------------------
# UI
# --------------------

ui <- fluidPage(
  
  titlePanel("Proteomics PCA Explorer"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      selectInput(
        "xpc",
        "X axis",
        choices = paste0("PC",1:5),
        selected = "PC1"
      ),
      
      selectInput(
        "ypc",
        "Y axis",
        choices = paste0("PC",1:5),
        selected = "PC2"
      ),
      
      sliderInput(
        "nproteins",
        "Number of proteins",
        min = 5,
        max = 50,
        value = 15
      )
      
    ),
    
    mainPanel(
      
      plotlyOutput("pca_plot"),
      
      br(),
      
      DTOutput("loading_table")
      
    )
  )
)


# --------------------
# SERVER
# --------------------

server <- function(input, output){
  
  output$pca_plot <- renderPlotly({
    
    plot_ly(
      data = scores,
      x = ~get(input$xpc),
      y = ~get(input$ypc),
      color = ~treatment,
      text = ~paste(
        "Sample:", sample_nr,
        "<br>Treatment:", treatment
      ),
      type = "scatter",
      mode = "markers+text",
      source = "pca"
    )
    
  })
  
  output$loading_table <- renderDT({
    
    pc <- input$xpc
    
    loadings %>%
      arrange(desc(abs(.data[[pc]]))) %>%
      select(protein, all_of(pc)) %>%
      head(input$nproteins)
    
  })
  
}

shinyApp(ui, server)


# Rewrite 2 ===============================================
# =======================================================


library(shiny)
library(plotly)
library(readxl)
library(dplyr)
library(DT)

# ==================================================
# DATA
# ==================================================

getwd()

df <- read_excel(
  "mer_proteins.xlsx",
  sheet = 4
)

sample_info <- df[,1:2]

protein_data0 <- df[,3:ncol(df)]

# protein_data <- log1p(protein_data)


col2remove <- which(apply(protein_data0, 2, sd) == 0)
View(col2remove)

# remove cols with no variation

protein_data <- protein_data0[, -col2remove]
View(protein_data)


# log transform
protein_data <- log1p(protein_data)
View(protein_data)


# ==================================================
# PCA
# ==================================================

pca <- prcomp(
  protein_data,
  center = TRUE,
  scale. = TRUE
)

scores <- as.data.frame(pca$x)

scores$sample_id <- sample_info$sample_nr

scores$treatment <- factor(
  sample_info$treatment,
  levels = c(0,1),
  labels = c("Control","Treated")
)

# ==================================================
# PCA LOADINGS
# ==================================================

loadings <- as.data.frame(pca$rotation)

loadings$protein <- rownames(loadings)

rownames(loadings) <- NULL

# ==================================================
# TREATMENT COMPARISON
# ==================================================

protein_df <- protein_data

protein_df$treatment <- sample_info$treatment

diff_results <- data.frame()

for(p in names(protein_data)){
  
  control_values <-
    protein_df[
      protein_df$treatment == 0,
      p
    ][[1]]
  
  treated_values <-
    protein_df[
      protein_df$treatment == 1,
      p
    ][[1]]
  
  control_mean <- mean(
    control_values,
    na.rm = TRUE
  )
  
  treated_mean <- mean(
    treated_values,
    na.rm = TRUE
  )
  
  test_result <- tryCatch(
    t.test(
      treated_values,
      control_values
    ),
    error = function(e) NULL
  )
  
  p_value <-
    ifelse(
      is.null(test_result),
      NA,
      test_result$p.value
    )
  
  diff_results <-
    rbind(
      diff_results,
      data.frame(
        protein = p,
        control_mean = control_mean,
        treated_mean = treated_mean,
        difference = treated_mean - control_mean,
        p_value = p_value
      )
    )
}

# FDR correction

diff_results$fdr <-
  p.adjust(
    diff_results$p_value,
    method = "BH"
  )

# Significance labels

diff_results$significance <-
  case_when(
    diff_results$fdr < 0.001 ~ "***",
    diff_results$fdr < 0.01  ~ "**",
    diff_results$fdr < 0.05  ~ "*",
    TRUE ~ "ns"
  )

# Sort by effect size

diff_results <-
  diff_results %>%
  arrange(desc(abs(difference)))

# ==================================================
# UI
# ==================================================

ui <- fluidPage(
  
  titlePanel(
    "Proteomics PCA Explorer"
  ),
  
  wellPanel(
    
    h4("About this application"),
    
    p(
      "This application provides an interactive version of the PCA analysis."
    ),
    
    p(
      "Users can explore principal components, identify proteins contributing most strongly to sample separation, and compare protein abundances between treatment groups."
    )
    
  ),
  
  sidebarLayout(
    
    sidebarPanel(
      
      selectInput(
        "xpc",
        "X-axis Principal Component",
        choices = paste0("PC",1:5),
        selected = "PC1"
      ),
      
      selectInput(
        "ypc",
        "Y-axis Principal Component",
        choices = paste0("PC",1:5),
        selected = "PC2"
      ),
      
      sliderInput(
        "nproteins",
        "Proteins displayed in tables",
        min = 5,
        max = 50,
        value = 15
      )
      
    ),
    
    mainPanel(
      
      plotlyOutput(
        "pca_plot",
        height = "600px"
      ),
      
      br(),
      
      tabsetPanel(
        
        tabPanel(
          "PCA Driver Proteins",
          DTOutput("loading_table")
        ),
        
        tabPanel(
          "Treatment Comparison",
          DTOutput("diff_table")
        )
        
      )
      
    )
  )
)

# ==================================================
# SERVER
# ==================================================

server <- function(input, output){
  
  output$pca_plot <- renderPlotly({
    
    plot_ly(
      data = scores,
      x = ~get(input$xpc),
      y = ~get(input$ypc),
      color = ~treatment,
      type = "scatter",
      mode = "markers+text",
      text = ~paste(
        "Sample ID:", sample_id,
        "<br>Treatment:", treatment
      ),
      hoverinfo = "text"
    )
    
  })
  
  output$loading_table <- renderDT({
    
    pc <- input$xpc
    
    loadings %>%
      arrange(
        desc(
          abs(.data[[pc]])
        )
      ) %>%
      select(
        Protein = protein,
        Loading = all_of(pc)
      ) %>%
      head(input$nproteins)
    
  },
  rownames = FALSE)
  
  output$diff_table <- renderDT({
    
    diff_results %>%
      select(
        Protein = protein,
        `Control Mean` = control_mean,
        `Treated Mean` = treated_mean,
        Difference = difference,
        `P-value` = p_value,
        FDR = fdr,
        Significance = significance
      ) %>%
      head(input$nproteins)
    
  },
  rownames = FALSE)
  
}

# ==================================================
# RUN APP
# ==================================================

shinyApp(ui, server)




# Rewrite 3 ===============================================
# =======================================================




if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
library(limma)



library(shiny)
library(plotly)
library(readxl)
library(dplyr)
library(DT)
library(BiocManager) 
library(limma)


# =====================================================
# LOAD DATA
# =====================================================

getwd()

df <- read_excel(
  "mer_proteins.xlsx",
  sheet = 4
)

sample_info <- df[,1:2]

protein_data0 <- df[,3:ncol(df)]



col2remove <- which(apply(protein_data0, 2, sd) == 0)
View(col2remove)

# remove cols with no variation

protein_data <- protein_data0[, -col2remove]
View(protein_data)


# log transform
protein_data <- log1p(protein_data)
View(protein_data)


sample_info$treatment <- factor(
  sample_info$treatment,
  levels = c(0,1),
  labels = c("Control","Treated")
)

# =====================================================
# PCA
# =====================================================

pca <- prcomp(
  protein_data,
  center = TRUE,
  scale. = TRUE
)

scores <- as.data.frame(pca$x)

scores$sample_id <- sample_info$sample_nr

scores$treatment <- sample_info$treatment

# =====================================================
# PCA LOADINGS
# =====================================================

loadings <- as.data.frame(
  pca$rotation
)

loadings$protein <- rownames(loadings)

rownames(loadings) <- NULL

# =====================================================
# LIMMA DIFFERENTIAL ANALYSIS
# =====================================================

design <- model.matrix(
  ~ treatment,
  data = sample_info
)

protein_matrix <- t(
  as.matrix(protein_data)
)

fit <- lmFit(
  protein_matrix,
  design
)

fit <- eBayes(fit)

limma_results <- topTable(
  fit,
  coef = 2,
  number = Inf,
  adjust.method = "BH"
)

limma_results$protein <- rownames(
  limma_results
)

rownames(limma_results) <- NULL

# =====================================================
# MERGE PCA LOADINGS + LIMMA RESULTS
# =====================================================

combined_results <- loadings %>%
  left_join(
    limma_results,
    by = "protein"
  )

# =====================================================
# UI
# =====================================================

ui <- fluidPage(
  
  titlePanel(
    "Proteomics PCA Explorer"
  ),
  
  wellPanel(
    
    h4("About this application"),
    
    p(
      "Interactive exploration of principal component analysis (PCA) and treatment-associated protein abundance differences."
    ),
    
    p(
      "The PCA plot visualises sample clustering while the accompanying tables identify proteins contributing to principal components and proteins associated with treatment effects."
    )
    
  ),
  
  sidebarLayout(
    
    sidebarPanel(
      
      selectInput(
        "xpc",
        "X-axis Principal Component",
        choices = paste0("PC",1:5),
        selected = "PC1"
      ),
      
      selectInput(
        "ypc",
        "Y-axis Principal Component",
        choices = paste0("PC",1:5),
        selected = "PC2"
      ),
      
      sliderInput(
        "nproteins",
        "Proteins displayed",
        min = 5,
        max = 50,
        value = 15
      )
      
    ),
    
    mainPanel(
      
      plotlyOutput(
        "pca_plot",
        height = "600px"
      ),
      
      br(),
      
      tabsetPanel(
        
        tabPanel(
          "Proteins Driving PCA",
          DTOutput("loading_table")
        ),
        
        tabPanel(
          "Treatment Effect (limma)",
          DTOutput("limma_table")
        )
        
      )
      
    )
    
  )
  
)

# =====================================================
# SERVER
# =====================================================

server <- function(input, output){
  
  output$pca_plot <- renderPlotly({
    
    plot_ly(
      data = scores,
      x = ~get(input$xpc),
      y = ~get(input$ypc),
      color = ~treatment,
      type = "scatter",
      mode = "markers+text",
      text = ~paste(
        "Sample ID:", sample_id,
        "<br>Treatment:", treatment
      ),
      hoverinfo = "text"
    )
    
  })
  
  # ---------------------------------------------
  # PCA DRIVER PROTEINS
  # ---------------------------------------------
  
  output$loading_table <- renderDT({
    
    pc <- input$xpc
    
    combined_results %>%
      arrange(
        desc(
          abs(.data[[pc]])
        )
      ) %>%
      select(
        Protein = protein,
        Loading = all_of(pc),
        logFC,
        adj.P.Val
      ) %>%
      head(input$nproteins)
    
  },
  rownames = FALSE)
  
  # ---------------------------------------------
  # LIMMA RESULTS
  # ---------------------------------------------
  
  output$limma_table <- renderDT({
    
    limma_results %>%
      select(
        Protein = protein,
        logFC,
        AveExpr,
        t,
        P.Value,
        adj.P.Val,
        B
      ) %>%
      head(input$nproteins)
    
  },
  rownames = FALSE)
  
}

# =====================================================
# RUN APP
# =====================================================

shinyApp(ui, server)



# Rewrite 3 ===============================================
# =======================================================




if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
library(limma)



library(shiny)
library(plotly)
library(readxl)
library(dplyr)
library(DT)
library(BiocManager) 
library(limma)


# =====================================================
# LOAD DATA
# =====================================================

getwd()

df <- read_excel(
  "mer_proteins.xlsx",
  sheet = 4
)

sample_info <- df[,1:2]

protein_data0 <- df[,3:ncol(df)]



col2remove <- which(apply(protein_data0, 2, sd) == 0)
View(col2remove)

# remove cols with no variation

protein_data <- protein_data0[, -col2remove]
View(protein_data)


# log transform
protein_data <- log1p(protein_data)
View(protein_data)


sample_info$treatment <- factor(
  sample_info$treatment,
  levels = c(0,1),
  labels = c("Control","Treated")
)

# =====================================================
# PCA
# =====================================================

pca <- prcomp(
  protein_data,
  center = TRUE,
  scale. = TRUE
)

scores <- as.data.frame(pca$x)

scores$sample_id <- sample_info$sample_nr

scores$treatment <- sample_info$treatment

# =====================================================
# PCA LOADINGS
# =====================================================

loadings <- as.data.frame(
  pca$rotation
)

loadings$protein <- rownames(loadings)

rownames(loadings) <- NULL

# =====================================================
# LIMMA DIFFERENTIAL ANALYSIS
# =====================================================

design <- model.matrix(
  ~ treatment,
  data = sample_info
)

protein_matrix <- t(
  as.matrix(protein_data)
)

fit <- lmFit(
  protein_matrix,
  design
)

fit <- eBayes(fit)

limma_results <- topTable(
  fit,
  coef = 2,
  number = Inf,
  adjust.method = "BH"
)

limma_results$protein <- rownames(
  limma_results
)

rownames(limma_results) <- NULL

# =====================================================
# MERGE PCA LOADINGS + LIMMA RESULTS
# =====================================================

combined_results <- loadings %>%
  left_join(
    limma_results,
    by = "protein"
  )

# =====================================================
# UI
# =====================================================

ui <- fluidPage(
  
  titlePanel(
    "Proteomics PCA Explorer"
  ),
  
  wellPanel(
    
    h4("About this application"),
    
    p(
      "Interactive exploration of principal component analysis (PCA) and treatment-associated protein abundance differences."
    ),
    
    p(
      "The PCA plot visualises sample clustering while the accompanying tables identify proteins contributing to principal components and proteins associated with treatment effects."
    )
    
  ),
  
  sidebarLayout(
    
    sidebarPanel(
      
      selectInput(
        "xpc",
        "X-axis Principal Component",
        choices = paste0("PC",1:5),
        selected = "PC1"
      ),
      
      selectInput(
        "ypc",
        "Y-axis Principal Component",
        choices = paste0("PC",1:5),
        selected = "PC2"
      ),
      
      sliderInput(
        "nproteins",
        "Proteins displayed",
        min = 5,
        max = 50,
        value = 15
      )
      
    ),
    
    mainPanel(
      
      plotlyOutput(
        "pca_plot",
        height = "600px"
      ),
      
      br(),
      
      tabsetPanel(
        
        tabPanel(
          "Proteins Driving PCA",
          DTOutput("loading_table")
        ),
        
        tabPanel(
          "Treatment Effect (limma)",
          DTOutput("limma_table")
        )
        
      )
      
    )
    
  )
  
)

# =====================================================
# SERVER
# =====================================================

server <- function(input, output){
  
  output$pca_plot <- renderPlotly({
    
    plot_ly(
      data = scores,
      x = ~get(input$xpc),
      y = ~get(input$ypc),
      color = ~treatment,
      type = "scatter",
      mode = "markers+text",
      text = ~paste(
        "Sample ID:", sample_id,
        "<br>Treatment:", treatment
      ),
      hoverinfo = "text"
    )
    
  })
  
  # ---------------------------------------------
  # PCA DRIVER PROTEINS
  # ---------------------------------------------
  
  output$loading_table <- renderDT({
    
    pc <- input$xpc
    
    combined_results %>%
      arrange(
        desc(
          abs(.data[[pc]])
        )
      ) %>%
      select(
        Protein = protein,
        Loading = all_of(pc),
        logFC,
        adj.P.Val
      ) %>%
      head(input$nproteins)
    
  },
  rownames = FALSE)
  
  # ---------------------------------------------
  # LIMMA RESULTS
  # ---------------------------------------------
  
  output$limma_table <- renderDT({
    
    limma_results %>%
      select(
        Protein = protein,
        logFC,
        AveExpr,
        t,
        P.Value,
        adj.P.Val,
        B
      ) %>%
      head(input$nproteins)
    
  },
  rownames = FALSE)
  
}

# =====================================================
# RUN APP
# =====================================================

shinyApp(ui, server)