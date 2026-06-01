#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
library(visNetwork)
library(shiny)
library(tidyverse)
library(igraph)
library(plotly)
library(RColorBrewer)
library(shinycssloaders)
library(shinythemes)
library(DT)
# library(shinymanager)
set.seed(2021)
options(spinner.type = 1)
getwd()



#########Read Edges table#############################################

lr_pairs_cells <-
  read.csv("LR_interactions.csv")

lr_pairs_cells <-
  lr_pairs_cells %>% dplyr::select("source_nodes", "target_nodes", everything()) # since first two columns are taken as source and target, bringing these two first
lr_pairs_cells <-
  lr_pairs_cells %>% dplyr::mutate_if(is.character, function(y)
    gsub("CM\\-", "CM_Sox10_neg", y)) # changing CM- to sox10nrg, otherwise grep later causes problems
lr_pairs_cells <-
  lr_pairs_cells %>% dplyr::mutate_if(is.character, function(y)
    gsub("CM\\+", "CM_Sox10_pos", y)) # changing CM- to sox10pos, otherwise grep later causes problems



#######Read Vertices table##############################################


lr_pairs_cells_vertices <-
  read.csv("LR_longformat.csv")


lr_pairs_cells_vertices <-
  lr_pairs_cells_vertices[, !(colnames(lr_pairs_cells_vertices) %in% c("heart_lit", "regdev_lit"))] #remove these columns as they contain all the pubmedids, which are huge and not required for this app


lr_pairs_cells_vertices <-
  lr_pairs_cells_vertices %>% dplyr::mutate_if(is.character, function(y)
    gsub("CM\\-", "CM_Sox10_neg", y)) # changing CM- to sox10neg, otherwise grep later causes problems
lr_pairs_cells_vertices <-
  lr_pairs_cells_vertices %>% dplyr::mutate_if(is.character, function(y)
    gsub("CM\\+", "CM_Sox10_pos", y))# changing CM- to sox10pos, otherwise grep later causes problems


lr_pairs_cells_vertices$label <- lr_pairs_cells_vertices$gene_symbol # assign gene symbol as the label for the nodes

# assign the color to the nodes
lr_pairs_cells_vertices$color[(lr_pairs_cells_vertices$celltype) == "FB"] <-
  "lightgreen" 
lr_pairs_cells_vertices$color[(lr_pairs_cells_vertices$celltype) == "EC"] <-
  "lightblue"
lr_pairs_cells_vertices$color[(lr_pairs_cells_vertices$celltype) == "CM_Sox10_pos"] <-
  "orange"
lr_pairs_cells_vertices$color[(lr_pairs_cells_vertices$celltype) == "CM_Sox10_neg"] <-
  "pink"

# assign size to the nodes
lr_pairs_cells_vertices$size <- lr_pairs_cells_vertices$log2FC # assign the log2fc as the size
lr_pairs_cells_vertices$size <-
  abs(lr_pairs_cells_vertices$size) + 10 # increase the size of nodes which are differentially regulated
lr_pairs_cells_vertices$size[(lr_pairs_cells_vertices$log2FC == 0)] <-
  2 # make the non differential genes small


# assign shape to the nodes - circle to the upregulated genes, square to the downregulated genes, and dot to the non differential genes
lr_pairs_cells_vertices$shape[(lr_pairs_cells_vertices$log2FC) > 0] <-
  "circle" 
lr_pairs_cells_vertices$shape[(lr_pairs_cells_vertices$log2FC) < 0] <-
  "square"
lr_pairs_cells_vertices$shape[(lr_pairs_cells_vertices$log2FC) == 0] <-
  "dot"


################Perform the scoring of the differential genes in the nodes#################################

# take only the nodes with Differential genes
deg_df <-
  lr_pairs_cells_vertices[lr_pairs_cells_vertices$log2FC != 0, ]

# copy the dataframe to get scoring
score_df <- deg_df
rownames(score_df) <- score_df$nodes

#select only the numeric columns
score_df <- score_df %>% dplyr::select(where(is.numeric))
score_df <- apply(score_df, 2,function(x)(abs(x)))
score_df <- as.data.frame(score_df)



#Scaling in log integers
cols_logs <- c("heart_lit_len", "regdev_lit_len")
score_df[cols_logs] <-  log(score_df[cols_logs] + 1, base = 2)


#Scaling in log decimal
cols_log_dec = c("pagerank", "pagerank_directed")
score_df[cols_log_dec] = log(score_df[cols_log_dec] + 0.0000001, base = 2) #since page rank is in very small numbers so taking lg2fc in decimal

#Scaling 0 to 1
score_df <-
  apply(score_df, 2, function(x)
    (x - min(x)) / (max(x) - min(x)))


# Changing Weight
cols_half <-
  c("influenced_genes_len",
    "influenced_celltypes_len",
    "influencer_celltypes_len")
cols_double <-  c("pagerank", "pagerank_directed")
score_df[, cols_half] <-
  apply(score_df[, cols_half], 2, function(x)
    (x * 0.5))
score_df[, cols_double] <-
  apply(score_df[, cols_double], 2, function(x)
    (x * 3))

# Inverting the scale
cols_invert <-
  c("heart_lit_len",
    "regdev_lit_len",
    "influencer_celltypes_len")

score_df[, cols_invert] <-
  apply(score_df[, cols_invert], 2, function(x)
    (1 - x))
score_df <- as.data.frame(score_df)


# calculate the sum of the scaled scores
cols_not_score <-
  c("means",
    "influenced_genes_len",
    "pagerank_directed",
    "size",
    "quality")
score_df <-
  score_df %>% mutate(score = rowSums(.[!(colnames(.) %in% cols_not_score)]))
score_df$nodes <- rownames(score_df)


# merge back score column with the score to the deg dataframe
deg_with_score <-
  merge(x = deg_df, y = score_df[, c("nodes", "score")], by = "nodes")
deg_with_score <-
  deg_with_score %>% dplyr::mutate_if(is.character, function(y)
    gsub("CM\\-", "CM_Sox10_neg", y))
deg_with_score <-
  deg_with_score %>% dplyr::mutate_if(is.character, function(y)
    gsub("CM\\+", "CM_Sox10_pos", y))


# merge this dataframe with the vertices dataframe to give score to each node
colnames(lr_pairs_cells_vertices)

lr_pairs_cells_vertices <-
  merge(
    lr_pairs_cells_vertices,
    deg_with_score[, c("ensemblID", "celltype", "score")],
    by = c("ensemblID", "celltype"),
    all = TRUE
  )

# remove any duplicates if present
lr_pairs_cells_vertices <-
  lr_pairs_cells_vertices[!duplicated(lr_pairs_cells_vertices$nodes), ]

lr_pairs_cells_vertices <-
  lr_pairs_cells_vertices %>% dplyr::select(nodes, everything()) %>% dplyr::arrange(desc(score))
lr_pairs_cells_vertices$pagerank <-
  as.numeric(lr_pairs_cells_vertices$pagerank)


###############Legend for network##############################################

legendNodes <- data.frame(
  label = c(
    "Fibroblasts",
    "Endocardial Cells",
    "Cardiomyocytes_Sox10_pos",
    "Cardiomyocytes_Sox10_neg",
    "Upregulated",
    "Downregulated",
    "No Change",
    "Small log2FC",
    "Large log2FC"
  ),
  color = c(
    "lightgreen",
    "lightblue",
    "orange",
    "pink",
    "lightgrey",
    "lightgrey",
    "lightgrey",
    "lightgrey",
    "lightgrey"
  ),
  shape = c("dot", "dot", "dot", "dot", "dot", "square", "dot", "dot", "dot"),
  size = c(20, 20, 20, 20, 20, 20, 2, 10, 30)
  
)

############### UI #################################################
# Define UI for application that draws a histogram

ui <- navbarPage(
  # Name that comes next to the tab
    "Ligand Receptor Analysis",
  # Make tabs 
    tabPanel(
      ########## LR Network tab ################
    "L-R Network",
    fluidPage(
      theme = shinytheme("lumen"),
      # Application title
      titlePanel("Ligand Receptor Network for Zebrafish Heart Regeneration at 7dpi"),
      
      ####### Show Network  ########
      fluidRow(
        visNetworkOutput("lr_network_out", width = "100%", height = "600px") %>% withSpinner(),
      ),
      br(),
      br(),
      
      ###### Show filter buttons ###########
      fluidRow(
       
          column(
            3,
            selectInput(         # button to select genes
              inputId = "gene",
              label = "Gene Symbol",
              choices = c(as.character(sort(
                        unique(lr_pairs_cells_vertices$gene_symbol)
              ))),
              selected = NULL,
              multiple = TRUE
            ),
          ),
          
          column(
            3,
            selectInput(      # button to select cell types
              inputId = "cell_type",
              label = "Cell Type",
              choices = c(as.character(sort(
                unique(lr_pairs_cells_vertices$celltype)
              ))),
              selected = NULL,
              multiple = TRUE
              ),
          ),
          column(3,         # button to download the table with filtered genes and cell types
                 downloadButton("download_df_filtered", "Download"))
        ),
        ######## Show the filtered table
      br(),
      br(),
      fluidRow(
        dataTableOutput('lr_df_filtered', width = "100%") %>% withSpinner(), 
      ),
      br(),
      br(),
      br(),
      ############ Show the scoring graph for the selected nodes  
      fluidRow(
      plotlyOutput(
        outputId = "score_plot",
        width = "100%",
        height = "600px"
      ) %>% withSpinner()
    )
    ),
    
  ),
  ############ Interactive graph Tab ###################
  tabPanel(
    "Interactive graph",
  # Show the interactive graph
  fluidRow(  
    plotlyOutput(
      outputId = "interactive_plot",
      width = "100%",
      height = "600px"
    ) %>% withSpinner(),
  ),
  br(),
  br(),
  br(),
  
  ##### Checkbox to select the Wang et al paper
  fluidRow( column(
    2),
    column(
    3,
    checkboxInput(
      "wangetal",
      "Genes common with Wang et al, 2020, Cell Reports",
      FALSE
    ),
  ),
  
  column(
    3,
    downloadButton("download_df_interactive", "Download")
  ),),
  br(),
  br(),
    dataTableOutput('lr_df_mouse', width = "100%") %>% withSpinner(),
  
  ),
  
  
  ############ Documentation Tab ########################
  
  tabPanel("Documentation",
           tags$html(
             tags$head(tags$title('Read Me')),
             tags$body(
               h1('Read Me'),
               
               div(
                 id = 'Steps',
                 class = 'simpleDiv',
                 tags$br(),
                 tags$br(),
                 tags$ol(
                   h3(
                     'Steps for making a list of relevant ligand and
                              receptors during in zebrafish heart regeneration
                              after cryoinjury at 7days post injury (7dpi) :',
                     br(),
                     br()
                   ),
                   
                   tags$li(
                     "Take genes from the chosen cell types at 7pdi
                               – periostin+ve Fibroblasts",
                     strong("(FB)"),
                     ", Sox10+ve",
                     strong("(CM_Sox10_pos)"),
                     ", Sox10-ve",
                     strong("(CM_Sox10_neg)"),
                     ", Endocardial cells",
                     strong("(EC).")
                   ),
                   tags$br(),
                   tags$li("Subset for Ligands and Recpetors curated from OmniPathDB;"),
                   tags$br(),
                   tags$li(
                     "To check the interactions we look at a human ligand receptor database. We have to convert the zebrafish genes to its human orthologue correspondents (we lose 20-30% of the genes in this step);"
                   ),
                   tags$br(),
                   tags$li(
                     "To consider an interaction, at least one interacting ligand or receptor must be differentially expressed between regenerating and physiological conditions. The other partner of the interaction can be also differentially expressed or just expressed;"
                   ),
                   tags$br(),
                   tags$li(
                     "Add relevant metadata (like positioning in a network of interactions, and literature search) and making comparisons to external databases (like neonatal mouse) to choose the most relevant genes."
                   )
                   
                 )
                 
                 
               ),
               
               div(
                 id = 'Steps',
                 class = 'simpleDiv',
                 tags$br(),
                 
                 tags$ul(
                   h3(
                     'Explanation of the columns in table of genes (alphabetical order): ',
                     br(),
                     br()
                   ),
                   
                   tags$li(
                     strong("celltype:"),
                     "the cell type in zebrafish (where it is differentially expressed)"
                   ),
                   tags$br(),
                   tags$li(
                     strong("common_injury_genes:"),
                     "analysis of different injury types at 7dpi in zebrafish. If a 'yes' is present it means that gene is differentially expressed in all injury types in zebrafish"
                   ),
                   tags$br(),
                   tags$li(
                     strong("Direction:"),
                     " interactions have directionality A ? B. This indicates whether a gene is the 'source' (ligand), the 'target'(receptor) or it could function as 'both'"
                   ),
                   tags$br(),
                   tags$li(strong("ensemblID:"), " ensembID of the gene"),
                   tags$br(),
                   tags$li(
                     strong("gene_symbol:"),
                     " zebrafish gene symbol. (if there was no match an ID appears instead of the symbol)"
                   ),
                   tags$br(),
                   tags$li(
                     strong("Gene_Symbol_human:"),
                     " the human gene orthologue used to check which zebrafish LR interacts with another LR"
                   ),
                   tags$br(),
                   tags$li(
                     strong("heart_lit_len:"),
                     " number of papers found on pubmed matching the gene name and its synonyms to the 'heart' topic"
                   ),
                   tags$br(),
                   tags$li(
                     strong("influenced_cell_types:"),
                     " cell types that express an interacting partner"
                   ),
                   tags$br(),
                   tags$li(
                     strong("influencer_cell_types:"),
                     " other cell types where the gene is present"
                   ),
                   tags$br(),
                   tags$li(
                     strong("influenced_genes:"),
                     " the genes that are expressed and interact with the gene "
                   ),
                   tags$br(),
                   tags$li(
                     strong("log2FC:"),
                     "log2 scale of the fold of the x fold of differential expression"
                   ),
                   tags$br(),
                   tags$li(strong("means:"), " mean of counts of sequenced fragments "),
                   tags$br(),
                   tags$li(
                     strong("pagerank:"),
                     " between 0 and 1. The higher the value the more important is relative position of a gene in the network of interactions between all the genes"
                   ),
                   tags$br(),
                   tags$li(
                     strong("reg_dev_lit_len:"),
                     " same as heart literature search but for regeneration/development topic"
                   ),
                   tags$br(),
                   tags$li(
                     strong("score:"),
                     " final score used to rank the genes (information about score calculation in another file)"
                   ),
                   tags$br(),
                   tags$li(
                     strong("synonyms:"),
                     " all the possible designations of the gene in human and fish"
                   ),
                   tags$br(),
                   tags$li(
                     strong("WangEtAl:"),
                     " neonatal mouse single-cell analysis from Wang et al 2020, Cell Reports. This column informs about the cell type where a gene is present (FB, or EC). It also informs whether this gene is differentially expressed in P1 vs P8 (P1, first week of life which is regenerative vs P8, the second week of life which is non-regenerative). 'dataset' means the gene is differentially expressed. 'interactome' means its expressed and interacting with an LR that’s differentially expressed."
                   ),
                   tags$br(),
                   
                   
                 )
                 
                 
               )
               
             )
           ))
)




########################## Server #####################################

# Define server logic required to draw a histogram

server <- function(input, output) {

 
    
    
  ########### make the network ##########
    
  output$lr_network_out <- renderVisNetwork({
    
    
    
    network_graph <- # convert LR edge dataframe to network using igraph
      graph_from_data_frame(lr_pairs_cells,
                            directed = F,
                            vertices = lr_pairs_cells_vertices)
    
    igraph_to_visnetwork <-
      toVisNetworkData(network_graph) # convert igraph to visnetwork
    
    igraph_to_visnetwork$nodes$label <-
      igraph_to_visnetwork$nodes$gene_symbol # give labels as gene_symbols
    
    
    #### output the network for visualizing
    
    visNetwork( 
      nodes = igraph_to_visnetwork$nodes,
      edges = igraph_to_visnetwork$edges,
      width = "100%",
      height = "100%"
    ) %>%
      visIgraphLayout(layout = "layout_with_graphopt", randomSeed = 2021) %>% # use layout to help deal with ohyics of visnetwork
      visOptions(                       # set options for network visualization
        highlightNearest = list(
          enabled = TRUE,
          degree = 1,
          hover = T
        ), ) %>%
      visLegend(                   # visualize the netowrk
        position = "right",
        width = 0.125,
        addNodes = legendNodes,
        zoom = FALSE,
        useGroups = FALSE
      ) %>%
      visInteraction(            # show navigation buttons
        zoomView = TRUE,
        navigationButtons = TRUE,
        dragNodes = TRUE
      )
  })
  
  ######## focus the network on the selected genes and cell types
  
  observe({
    if ((!isTruthy(input$gene)) & (!isTruthy(input$cell_type))) { # if nothing is selected
      visNetworkProxy("lr_network_out") %>%
        visOptions(
          highlightNearest = list(
            enabled = TRUE,
            degree = 0,
            hover = T
          )) %>%
        visFit(nodes = NULL) %>%
        visInteraction(
          zoomView = TRUE,
          navigationButtons = TRUE,
          dragNodes = TRUE
        )
      
    } else if ((isTruthy(input$gene)) & # if only gene is selected
               (!isTruthy(input$cell_type))) {
      gene_selection <-
        
        lr_pairs_cells_vertices %>% dplyr::filter(gene_symbol %in% input$gene) %>% pull(nodes)
      visNetworkProxy("lr_network_out") %>%
        visOptions(
          highlightNearest = list(
            enabled = TRUE,
            degree = 1,
            hover = T
          )) %>%
        visFit(nodes = gene_selection) %>%
        visSelectNodes(id = gene_selection)  %>%
        visInteraction(
          zoomView = TRUE,
          navigationButtons = TRUE,
          dragNodes = TRUE
        )
      
    } else if (!(isTruthy(input$gene))  & # if only cell type is selected
               (isTruthy(input$cell_type))) {
      cell_selection <-
        lr_pairs_cells_vertices %>% dplyr::filter(celltype %in% input$cell_type) %>% pull(nodes)
      
      
      visNetworkProxy("lr_network_out") %>%
        visOptions(       
          highlightNearest = list(
            enabled = TRUE,
            degree = 0,
            hover = T
          )) %>%
        visFit(nodes = cell_selection) %>%
        visSelectNodes(id = cell_selection) %>%
        visInteraction(
          zoomView = TRUE,
          navigationButtons = TRUE,
          dragNodes = TRUE
        )
      
      
    } else { # if both gene and cell type are selected
      gene_cell_selection <-
        lr_pairs_cells_vertices %>% dplyr::filter(celltype %in% input$cell_type) %>% 
        dplyr::filter(gene_symbol %in% input$gene) %>% 
        pull(nodes)
      
      visNetworkProxy("lr_network_out") %>%
        visOptions(
          highlightNearest = list(
            enabled = TRUE,
            degree = 1,
            hover = T
          )) %>%
        visFit(nodes = gene_cell_selection) %>%
        visSelectNodes(id = gene_cell_selection) %>%
        visInteraction(
          zoomView = TRUE,
          navigationButtons = TRUE,
          dragNodes = TRUE
        )
      
      
    }
    
  })

  
  
  
  ######## Dataframe to plot the score bar plot based on filtered genes ########

  for_score_graph <- reactive({
    
    # split nodes to gene symbol and cell types
    df_score <-
      score_df %>% separate(
        col = nodes,
        into = c("gene_symbol", "celltype"),
        remove = FALSE,
        extra = "merge"
      ) %>%
      select(-all_of(c(cols_not_score, "score")))
    
    # filter the dataframe for the gene and cell types
    
    if ((!isTruthy(input$gene)) & 
        (!isTruthy(input$cell_type))) {       # if nothing is selected
      df_score
      
    } else if ((isTruthy(input$gene)) &
               (!isTruthy(input$cell_type))) {  # if only gene is selected
      df_score <- df_score %>%
        dplyr::filter(gene_symbol %in% input$gene)
      
      
      
    } else if (!(isTruthy(input$gene))  &
               (isTruthy(input$cell_type))) {   # if only cell type is selected
      df_score <- df_score %>%
        dplyr::filter(celltype %in% input$cell_type)
      
    } else {                                    # if both gene and cell types are selected
      df_score <- df_score %>%
        dplyr::filter(gene_symbol %in% input$gene) %>%
        dplyr::filter(celltype %in% input$cell_type)
      
      
    }
    
    
    # convert the data frame into long format to allow plotting by ggplot
    
    df_score <-  df_score %>%
      pivot_longer(cols = (-c("nodes", "gene_symbol",  "celltype")),
                   names_to = "parameters",
                   values_to = "weighted_scores")
    
    df_score
  })

  ##### Make the plotly graph for the score of the nodes #########
  
  output$score_plot <- renderPlotly({
    score_plot <-
      ggplot(data = for_score_graph(), aes(x = nodes, y = weighted_scores, fill = parameters)) +
      geom_bar(stat = "identity") +
      theme_bw() +
      theme(axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      )) +
      labs(x = 'Nodes from Network', y = 'Weighted Scores')
    
    ggplotly(score_plot) %>%
      layout(title = list(
        text = paste0(
          "Score distribution of the Differential Genes that form the network",
          "<br>",
          "<sup>",
          "If the node is not differential no graph is shown",
          "</sup>"
        )
      ),
      margin = list(t = 50))%>% # adding margin to remove overlapping title to the graph
      config(displaylogo = FALSE
      )
  })
  
  
  
  ############ Data frame to show filtered data #####################
  
  
  # make a vector to remove unwanted columns 
  
  remove_cols <- c(
    "nodes",
    "label",
    "color",
    "size",
    "shape",
    "influenced_genes_len",
    "influenced_celltypes_len",
    "influencer_celltypes_len",
    "quality",
    "pagerank_directed"
  )
  
  # subset the nodes dataframe to remove unqanted columns
  df_to_render <-  lr_pairs_cells_vertices %>%
    dplyr::select(-dplyr::all_of(remove_cols)) %>%
    dplyr::select("gene_symbol", "celltype", "score", everything())
  
  
  # Filter the dataframe
  
  df_filtered <- reactive({
    if ((!isTruthy(input$gene)) & (!isTruthy(input$cell_type))) { # if nothing is selected
      df_to_render
      
    } else if ((isTruthy(input$gene)) &                           # if only gene is selected
               (!isTruthy(input$cell_type))) {
      df_to_render %>% dplyr::filter(gene_symbol %in% input$gene)

      
      
    } else if (!(isTruthy(input$gene))  &                        # if only cell type is selected
               (isTruthy(input$cell_type))) {
      df_to_render %>% dplyr::filter(celltype %in% input$cell_type)
      
      
    } else {                                                      # if both gene and cell types are selected
      df_to_render %>% dplyr::filter(gene_symbol %in% input$gene) %>%
        dplyr::filter(celltype %in% input$cell_type)
      
    }
  })
  
  # Output the filtered data frame  
  
  output$lr_df_filtered <-   renderDataTable(datatable(
    df_filtered(),
    filter = list(position = 'top', clear = FALSE),
    options = list(pageLength = 5, scrollX = TRUE)
  ),
  server = TRUE)
  
  
  # Download the filtered data frame
  
  output$download_df_filtered <- downloadHandler(
    filename = function() {
      gene_name <- paste0(input$gene, collapse = "_")
      cell_name <- paste0(input$cell_type, collapse = "_")
      
      if ((!isTruthy(input$gene))) {
        gene_name <- "all_genes"
      }
      
      if (!isTruthy(input$cell_type)) {
        cell_name <- "all_cells"
      }
      paste0(gene_name, "_", cell_name, ".csv")
    },
    content = function(file) {
      write.csv(df_filtered(), file, row.names = FALSE)
    }
  )
  
  
  
  ################# Make the data frame for the differential genes ##################
  
  
  ## check box to allow filtering of the data for wang et al
  
  df_interactive <- reactive({
    df_int <- df_to_render[df_to_render$log2FC != 0, ]
      if (input$wangetal == FALSE) {
      df_int
  
    } else {
      df_int[grepl(pattern = "dataset", x = df_int$WangEtAl), ]
      
    }
  })
  
  # Output wang et al filtered differential gene data
  
  output$lr_df_mouse <- renderDataTable(datatable(
    df_interactive(),
    filter = list(position = 'top', clear = FALSE),
    options = list(pageLength = 5, scrollX = TRUE)
  ),
  server = TRUE)
  
  # download wang et al filtered differential gene data
  
  output$download_df_interactive <- downloadHandler(
    filename = function() {
      if (input$wangetal == FALSE) {
        fname <- "all_nodes.csv"
        
      } else {
        fname <- "zf_wang_common.csv"
        
      }
      
    },
    content = function(file) {
      write.csv(df_interactive(), file, row.names = FALSE)
    }
  )
  
  
  ################# Interactive Graph plot #################
  
  symbols_all <- c("circle", "square", "triangle-up", "diamond")
  names(symbols_all) <- unique(df_to_render$celltype)
  
  symbols_plotly <- reactive({
    symbols_all[names(symbols_all) %in% df_interactive()$celltype]
  })
  
  # Show the interactive plot 
  
  output$interactive_plot <- renderPlotly({
    p <- ggplot(
      data = df_interactive(),
      aes(
        x = log(heart_lit_len),
        y = log(regdev_lit_len),
        size = pagerank,
        color = log2FC,
        shape = celltype,
        text = sprintf(
          "Gene Symbol : %s<br>Human Gene : %s<br>Cell type : %s<br>Influenced Genes : %s<br>Influenced Cell types : %s<br>Present in Wang et al : %s<br>Common in Core Injury genes(Marius)  : %s<br>No. of Heart Papers : %s<br>No. of Regernation Papers : %s",
          gene_symbol,
          Gene_Symbol_human,
          celltype,
          influenced_genes,
          influenced_celltypes,
          WangEtAl,
          common_injury_genes,
          heart_lit_len,
          regdev_lit_len
        )
      )
    ) +
      geom_point() +
      scale_shape_manual(values = symbols_plotly()) +
      theme_bw() +
      scale_color_gradient2(
        midpoint = 0,
        low = "blue",
        high = "red",
        space = "Lab"
      ) + 
      labs(x = 'Heart Literature (in Log Scale)', y = 'Regeneration Literature (in Log Scale)', title = "Ligand Receptor Interactive Graph for Zebrafish Heart Regeneration at 7dpi")
    
    ggplotly(p, tooltip = "text") %>% 
      config(displaylogo = FALSE
            )
    
  })
  
}

####### Run the application ###########

shinyApp(ui = ui, server = server)
