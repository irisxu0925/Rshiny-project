library(shiny)
library(bslib)
library(ggplot2)
library(colourpicker)
library(dplyr)
library(DT)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(PCAtools)


# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Final project app"),
  
  tabsetPanel(
    #sample
    tabPanel('Sample',
             sidebarLayout(
               #input matrix and control
               sidebarPanel(
                 fileInput("sample_data", "Please input a sample information file: ", accept = ".csv"),
                 actionButton("button1", "Submit")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel('Summary', DT::dataTableOutput('sum_table')),
                   tabPanel('Table', DT::dataTableOutput('std_table')),
                   tabPanel('Plots', sidebarLayout(
                     sidebarPanel(
                       #Select the column draw the histogram
                       selectInput("columnSelect", "Select Column:", choices = NULL)
                     ),
                     mainPanel(
                       plotOutput("histogram")
                     )
                   )) 
                 )
               )
             )
             
    ),
    #count
    tabPanel('Counts',
             sidebarLayout(
               #input matrix and control
               sidebarPanel(
                 fileInput("count_data", "Please input a normalized counts matrix file: ", accept=c(".csv",".tsv",".txt")),
                 sliderInput("slider1",
                             "genes with at least X percentile of variance",
                             min = 1,max = 100,value = 30),
                 sliderInput("slider2",
                             "genes with at least X samples that are non-zero",
                             min = 1,max = 20,value = 5),
                 actionButton("button2", "Submit")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel('Summary', DT::dataTableOutput("sum_count")),
                   tabPanel('Diagnostic scatter plot', plotOutput("ds_plot")),
                   tabPanel('Heatmap',plotOutput("heatmap")),
                   tabPanel('PCA',sidebarLayout(
                     sidebarPanel(
                       #choose the pc of x/y axis
                       radioButtons('xpca', 'Please select a PC for the x-axis',
                                    c('PC1','PC2','PC3','PC4')),
                       radioButtons('ypca', 'Please select a PC for the y-axis',
                                    c('PC1','PC2','PC3','PC4')),
                     ),
                     mainPanel(plotOutput("pca"))
                   ))
                 )
               )
             )
    ),
    #DE
    tabPanel('DE',
             sidebarLayout(
               #either input DE result matrix or count matrix to do DE analysis
               sidebarPanel(
                 fileInput("DE_data", "Please input a differencial expression matrix file: ", accept = c(".csv",".tsv",".txt")),
                 sliderInput("bins2",
                             "padj threshold",
                             min = 0,max = 1,value = 0.05),
                 actionButton("button3", "Submit")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel('Summary', DT::dataTableOutput('sum_volcano')),
                   tabPanel('Volcano Plot',sidebarLayout(
                     sidebarPanel(
                       #choose the color of volcano plot
                       colourInput('color3', 'Choose a color for up regulated gene', 'red'), 
                       colourInput('color4', 'Choose a color for down regulated gene', 'blue'),
                     ),
                     mainPanel(plotOutput("scatter"))
                   )),
                   tabPanel('Histogram',sidebarLayout(
                     sidebarPanel(
                       #choose the x axis of histogram
                       radioButtons('xaxis', 'Please select a column for the x-axis',
                                    choices = c('log2FoldChange','pvalue'))
                     ),
                     mainPanel(plotOutput("Hist"))
                   ))
                 )
               )
             )
    ),
    #GSEA
    tabPanel('GSEA',
             sidebarLayout(
               #input matrix and control
               sidebarPanel(
                 fileInput("GSEA_data", "Please input the fgsea result file: ", accept = ".csv"),
                 actionButton("button4", "Submit")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel('Barplot',sidebarLayout(
                     sidebarPanel(
                       #filter table by adjusted p-value 
                       sliderInput("padj_thr","set adjusted p-value",min = 0,max = 1,value = 0.05),
                       actionButton("button7", "filter")
                     ),
                     mainPanel(plotOutput("gsea_bar"))
                   )),
                   tabPanel('Sort_table',sidebarLayout(
                     sidebarPanel(
                       #filter table by adjusted p-value 
                       sliderInput("filter1","set adjusted p-value",min = 0,max = 1,value = 0.05),
                       #select all, positive or negative NES pathways
                       radioButtons('pathway', 'Please select NES pathways',
                                    choices = c('all','positive','negative')),
                       actionButton("button5", "filter"),
                       #export current filtered and displayed table results
                       DTOutput("table"),
                       downloadButton("downloadBtn", "Download Data")
                     ),
                     mainPanel(DT::dataTableOutput('gsea_table'))
                   )),
                   tabPanel('Scatter_plot',sidebarLayout(
                     sidebarPanel(
                       #filter table by adjusted p-value 
                       sliderInput("filter2","set adjusted p-value",min = 0,max = 1,value = 0.05),
                       actionButton("button6", "filter")
                     ),
                     mainPanel(plotOutput("gsea_scatter"))
                   ))
                   
                 )
               )
             )
    )
    
  )
)



# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  #define the file size limit
  options(shiny.maxRequestSize = 30*1024^2)
  
  #function to detect the file format eg:csv/tsv
  detect_separator <- function(file_path) {
    #read the first few lines of the file
    sample_lines <- readLines(file_path, n = 5)
    potential_separators <- c(',', '\t')
    for (separator in potential_separators) {
      if (all(grepl(separator, sample_lines))) {
        return(separator)
      }
    }
    #defult
    return(',')
  }
  
  #sample data input
  load_data <- reactive({
    req(input$sample_data)
    #if(detect_separator(input$sample_data)==","){
    #  df <- read.csv(input$sample_data$datapath,header = TRUE)
    #}else{df <- read.table(input$sample_data$datapath,sep = "\t",header = TRUE)}
    df <- read.csv(input$sample_data$datapath,header = TRUE)
    return(df)
  })
  
  
  #DE data input
  load_DE<- reactive({
    req(input$DE_data)
    df <- read.table(input$DE_data$datapath,sep = "\t", header = TRUE)
    return(df)
  })
  
  
  
  ##sample
  #summarizing sample info
  sample_summarize <- function(df) {
    summary_info <- sapply(names(df), function(col_name) {
      # preprocess data
      df$PMI[df$PMI == "unk"] <- NA
      df$PMI <- as.numeric(df$PMI)
      #deal with mRNA.Seq.reads
      if (col_name == "mRNA.Seq.reads" && any(!is.na(df[[col_name]]))) {
        df[[col_name]] <- as.numeric(gsub(",", "", df[[col_name]]))
      }
      #get column and column type
      col <- df[[col_name]]
      col_type <- class(col)[1]
      #deap with PMI 
      if (col_name == "PMI" && col_type == "character" && any(!is.na(col))) {
        col <- as.numeric(col)
      }
      #calculate
      if (col_type %in% c("numeric", "integer")) {
        mean_sd <- paste(round(mean(col, na.rm = TRUE), 2), " ( +/-", round(sd(col, na.rm = TRUE), 2), ")", sep = "")
      } else if (col_type == "factor") {
        mean_sd <- paste(levels(col), collapse = ", ")
      } else if (col_type == "character") {
        mean_sd <- paste(col, collapse = ", ")
      }
      
      c(col_name = col_name, data_type = col_type, "Mean (sd) or Distinct Values" = mean_sd)
    })
    return(as.data.frame(t(summary_info)))
  }
  
  #plot continuous sample info
  observe({
    #dynamically update the selection column
    choices <- names(load_data())
    updateSelectInput(session, "columnSelect", choices = choices, selected = choices[1])
  })
  histogram_plot <- function(df,col){
    # preprocess data
    df$PMI[df$PMI == "unk"] <- NA
    df$PMI <- as.numeric(df$PMI)
    #deal with mRNA.Seq.reads
    col_name = col
    if (col_name == "mRNA.Seq.reads" && any(!is.na(df[[col_name]]))) {
      df[[col_name]] <- as.numeric(gsub(",", "", df[[col_name]]))
    }
    #get column and column type
    col <- df[[col_name]]
    col_type <- class(col)[1]
    #deap with PMI 
    if (col_name == "PMI" && col_type == "character" && any(!is.na(col))) {
      col <- as.numeric(col)
    }
    if( col_type %in% c("numeric", "integer")){
      plot <- ggplot(df, aes_string(x = col)) + geom_histogram(color="black", fill="white")
      return(plot)
    }
  }
  
  #output sample info
  clickCount <- reactiveVal(0)
  observeEvent(input$button1, {
    clickCount(clickCount() + 1)
    sample_df <- load_data()
    output$sum_table<- DT::renderDataTable({
      sample_summarize(sample_df)
    })
    #output table of sample
    output$std_table <- DT::renderDataTable({
      return(sample_df)
    })
    
  })
  #output histogram
  output$histogram <- renderPlot({
    sample_df <- load_data()
    col <- input$columnSelect
    print(col)
    histogram_plot(sample_df,col)
  })
  
  
  
  ##count
  #norm counts data input
  load_norm_counts <- reactive({
    req(input$count_data)
    df <- read.table(input$count_data$datapath,sep = "\t",header = TRUE)
    return(df)
  })
  
  #summarizing count
  count_summarize <- function(df, slider1, slider2){
    #remove rows with na
    df <- df[complete.cases(df),]
    #number of samples
    num_samples <- ncol(df[,-1]) 
    #number of genes
    num_genes <- nrow(df) 
    #calculate variance of each gene
    df$rowVar<- apply(df[,-1], 1, var) 
    #determine the X percentile cutoff for variance
    variance_cutoff <- quantile(df$rowVar, probs = slider1/100)
    df$threshold[df$rowVar >= variance_cutoff | rowSums(df==0) >= slider2] = "P"
    df$threshold[df$rowVar < variance_cutoff | rowSums(df==0) < slider2] = "F"
    #number and % of genes passing current filter
    num_passed <- sum(df$threshold=="P")
    per_passed <- num_passed/num_genes*100
    #number and % of genes not passing current filter
    num_failed <- sum(df$threshold=="F")
    per_failed <- num_failed/num_genes*100
    return(data.frame("number of samples" = num_samples, "number of genes" = num_genes, 
                      "number of genes passing current filter" = num_passed,
                      "% of genes passing current filter" = per_passed,
                      "number of genes not passing current filter" = num_failed,
                      "% of genes not passing current filter" = per_failed
    ) %>% as_tibble())
  }
  #diagnostic scatter plot
  scatterplot <- function(df, variance_threshold, zero_threshold){
    #median count
    df$median <- apply(df[,-1], 1, median) 
    #variance
    df$var <- apply(df[,-1], 1, var) 
    #no. of zeros
    df$zeros <-  rowSums(df == 0)
    #filter
    variance_cutoff <- quantile(df$var, probs = variance_threshold/100)
    df$filter_status <- ifelse(df$var > variance_cutoff/100 & df$zeros < zero_threshold, "P", "F")
    #Plot by median and variance
    plot1 <- #mutate(df,  `log10(var)` = log10(var)) %>% 
      ggplot(df, aes(x = median, y = var, color = filter_status)) +
      geom_point() +
      scale_color_manual(values = c("F" = "lightgray", "P" = "darkblue")) +
      labs(title = paste("Diagnostic Plot: Median vs Variance"),
           x = "Median", y = "log(Variance)", color = "Filter Status") + 
      xlim(0, 25000) +
      scale_x_log10() + scale_y_log10()
    #Plot by median and zeros
    plot2 <- ggplot(df, aes(x = median, y = zeros, color = filter_status)) +
      geom_point() +
      scale_color_manual(values = c("F" = "lightgray", "P" = "darkblue")) + xlim(0, 25000) +
      scale_x_log10() 
    labs(title = paste("Diagnostic Plot: Median vs Zeros"),
         x = "Median", y = "Zeros", color = "Filter Status") 
    #Combine plots horizontally
    return(plot1 / plot2)
  }
  #heatmap
  heatmap_function <- function(df, variance_threshold, zero_threshold) {
    # Calculate row-wise variances
    df$variance <- apply(df[, -1], 1, var)
    # Number of zeros
    df$zeros <- rowSums(df == 0)
    # Filter
    variance_cutoff <- quantile(df$variance, probs = variance_threshold / 100)
    df <- dplyr::filter(df, variance > variance_cutoff & zeros < zero_threshold)
    # Log transforming without the first and last 2 columns (variance and zeros)
    log_df <- log10(as.matrix(df[, -c(1, ncol(df) - 1, ncol(df))] + 1))
    # Set rownames on the log-transformed matrix
    rownames(log_df) <- df[, 1]
    
    # Create a heatmap using pheatmap
    col.pal <- colorRampPalette(brewer.pal(8, 'RdBu'))(256)
    pheatmap(log_df, color = col.pal, show_rownames = FALSE, show_colnames = FALSE)
  }
  
  #PCA
  pca_plot <- function(df, PC1, PC2){
    pca <- pca(df[,-1])
    biplot(pca, x = PC1, y = PC2) %>%
      return()
  }
  
  
  #output summary data of counts file
  clickCount <- reactiveVal(0)
  observeEvent(input$button2, {
    clickCount(clickCount() + 1)
    norm_df <- load_norm_counts()
    output$sum_count <- DT::renderDataTable({
      count_summarize(norm_df, isolate(input$slider1), isolate(input$slider2))
    })
    #output diagnostic scatter plot
    output$ds_plot <- renderPlot({
      scatterplot(norm_df,isolate(input$slider1), isolate(input$slider2))
    })
    #output heatmap
    output$heatmap <- renderPlot({
      heatmap_function(norm_df,isolate(input$slider1), isolate(input$slider2))
    })
    output$pca <- renderPlot(pca_plot(norm_df, input$xpca, input$ypca))
  })

  
  #DE
  #threshold setting
  threshold_function <- function(df, padj_threshold){
    UP <- which(df$padj<padj_threshold & df$log2FoldChange>0)
    DOWN <- which(df$padj<padj_threshold & df$log2FoldChange<0)
    df <- df %>% mutate(volc_plot_status = rep("NS",nrow(df))) %>% relocate(volc_plot_status,.before =  log2FoldChange)
    df$volc_plot_status[UP] = "UP"
    df$volc_plot_status[DOWN] = "DOWN"
    df <- cbind(genes = rownames(df), df)
    return(df)
  }
  #P-value/logFC plot (heatmap)
  plot_pvals <- function(df,xaxis) {
    p <-  ggplot(df,aes(x = !!sym(xaxis))) + 
      geom_histogram(binwidth = 0.01, fill = "lightblue2", color = "black") +
      labs(title = "Histogram of raw p-values from DE analysis",
           x =  xaxis,
           y = "count") +
      theme_minimal()
    return(p)
  }
  #volcano plot
  plot_volcano <- function(df, padj_threshold, color3, color4) {
    labeled_df <- threshold_function(df, padj_threshold)
    p <- mutate(labeled_df, `padj<0.05`=padj<0.05, `-log10(padj)` = -log10(padj)) %>% 
      ggplot(aes(x = log2FoldChange, y = `-log10(padj)`, color=volc_plot_status)) +
      geom_point()+
      geom_hline(yintercept = 0,linetype = "dashed", colour = "black") +
      scale_color_manual(values = c("UP" = color3, "DOWN" = color4, "NS" = "grey")) +
    labs(title = "Volcano plot of differential expression results")
    return(p)
  }
  
  #DE output
  clickCount <- reactiveVal(0)
  observeEvent(input$button3, {
    clickCount(clickCount() + 1)
    de_df <- load_DE()
    #output sum table
    output$sum_volcano <- DT::renderDataTable({
      datatable(de_df, options = list(ordering = TRUE, searching = TRUE))
    })
    #output histogram
    output$Hist <- renderPlot({
      plot_pvals(de_df, input$xaxis)
    })
    #output volcano plot
    output$scatter <- renderPlot({
      plot_volcano(de_df, isolate(input$bins2), isolate(input$color3), isolate(input$color4))
    })
  })
  
  
  ##GSEA
  
  #GSEA data input
  load_GSEA <- reactive({
    req(input$GSEA_data)
    df <- read.csv(input$GSEA_data$datapath, header = TRUE)
    return(df)
  })
  
  #barplot
  barplot <- function(fgsea_results,padj_thr){
    #filter by FDR threshold and subset significant gene sets by the direction of their NES
    top_positive_nes <- fgsea_results %>% filter(padj < padj_thr)
    top_negative_nes <- fgsea_results %>% filter(padj < padj_thr)
    top_nes <- rbind(top_positive_nes,top_negative_nes)
    #plot
    p <- top_nes %>%
      mutate(pathway = forcats::fct_reorder(pathway, NES)) %>%
      ggplot() +
      geom_bar(aes(x=pathway, y=NES, fill = ifelse(NES > 0, 'Positive', 'Negative')), stat='identity', width = 0.7) +
      scale_fill_manual(values = c('Positive' = 'red', 'Negative' = 'blue')) + 
      theme_minimal() +
      ggtitle('fgsea results for Hallmark MSigDB gene sets') +
      ylab('Normalized Enrichment Score (NES)') +
      xlab('') +
      coord_flip() + 
      theme(axis.text.y = element_text(size = 5), axis.title.y = element_text(size = 5),
            plot.title = element_text(size = 8), legend.position = "none")
    
    return(p)
  }
  #filter result table by padj
  filter_result <- function(fgsea_results, padj_thr, pathway){
    #filter by FDR threshold and subset significant gene sets by the direction of their NES
    top_positive_nes <- fgsea_results %>% filter(padj < padj_thr & NES > 0)
    top_negative_nes <- fgsea_results %>% filter(padj < padj_thr & NES < 0)
    if(pathway == "positive"){return(top_positive_nes)}
    else if(pathway == "negative"){return(top_negative_nes)}
    else if(pathway == "all"){return(rbind(top_positive_nes,top_negative_nes))}
  }
  #scatter plot
  scatter_gsea <- function(fgsea_results, padj_thr){
    p <- mutate(fgsea_results, `padj<padj_thr`=padj<padj_thr, `-log10(padj)` = -log10(padj)) %>% 
      ggplot(aes(x = NES, y = `-log10(padj)`)) +
      geom_point(aes(color = ifelse(padj < padj_thr, "True", "False")))+
      scale_color_manual(values = c("False" = "grey", "True" = "red")) +
      labs(title = "Scatter plot of GSEA results")
    return(p)
  }
  
  #output
  observeEvent(input$button4, {
    clickCount(clickCount() + 1)
    fgsea_results <- load_GSEA()
    #output histogram
    output$gsea_bar <- renderPlot({
      input$button7
      barplot(fgsea_results, isolate(input$padj_thr))
    })
    #output sum table
    output$gsea_table <- DT::renderDataTable({
      input$button5
      filter_result(fgsea_results, isolate(input$filter1), isolate(input$pathway))
    })
    #Download table
    output$downloadBtn <- downloadHandler(
      filename = function() {
        paste("fgsea_data", Sys.Date(), ".csv", sep = "_")
      },
      content = function(file) {
        write.csv(filter_result(fgsea_results, isolate(input$filter1), isolate(input$pathway)), file)
      }
    )
    #output scatter plot
    output$gsea_scatter <- renderPlot({
      input$button6
      scatter_gsea(fgsea_results, isolate(input$filter2))
    })
  })
  
  
   
  
}

# Run the application 
shinyApp(ui = ui, server = server)

