
#This is the shiny application of scHD4E R Github Package
#scHD4E package is developed for the scHD4E ensemble learning based DE analysis
#This package include four individual and scDEA methods
#Simultaneosly run the the methods and calculate p-values and adjusted p-values

library(shiny)
library("scDEA")
library(MAST)
library(data.table)
library(plyr)
library(dplyr)
library(gmodels)
library(SwarnSeq)
library(ROSeq)
library(edgeR)
library(limma)

# Define UI
ui <- fluidPage(
  titlePanel("scHD4E Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Choose CSV File"),
      fileInput("file1", "Choose CSV File1"),
      numericInput("ncluster", "Enter Number of Clusters:", value = 5),
      actionButton("run_analysis", "Run Analysis")
    ),
    
    mainPanel(
      tableOutput("results_table")
    )
  )
)

scHD4E<-function(scData,condition,ncluster){
  #This looping is necessary when a number of genes of a dataset contains all zeros across cell#
  
  index<-c()
  for (i in 1:dim(scData)[1]){
    index[i]=length(which(as.numeric(scData[i,])!=0))/length(as.numeric(scData[i,]))
  }
  keep1<-which(index>=0.001)
  filterData<-scData[keep1,]
  
  ############ROSeq##########
  pred.ROSeq <- Sys.time()
  output<-ROSeq(countData=filterData, condition = condition, numCores=1)
  data_ROSeq<-data.frame(output)
  
  #### The program of the next section is essential when a method produce missing p-values##
  dat1<-data.frame(scData[,1:2])
  library(tibble)
  dat1 = rownames_to_column(dat1, "Plants")
  dat2 = rownames_to_column(data_ROSeq, "Plants")
  library(dplyr)
  dat = full_join(dat1, dat2, )
  dat = dat %>% replace(is.na(.), 0.90)   ####when a gene produce missing value we replace it by 0.90. Thats means we assume it is non-differential. This is nessary to avoid the disruption of algorithm#####
  #####################################
  
  ROSeq<-dat$pVals
  end.ROSeq<- Sys.time()
  time.ROSeq<- difftime(end.ROSeq, pred.ROSeq, units = "mins")
  cat("Run time for ROSeq: ", time.ROSeq, "min","\n")
  
  #########Seurat##############
  pred.Seurat <- Sys.time()
  Pvals_Seurat_full<-scDEA_individual_methods(
    raw.count=scData,
    cell.label=condition,
    is.normalized = FALSE,
    verbose = TRUE,
    BPSC = F,
    DEsingle = F,
    DESeq2 = F,
    edgeR = F,
    MAST = F,
    monocle = F,
    scDD = F,
    Ttest = F,
    Wilcoxon = F,
    limma = F,
    Seurat = TRUE,
    zingeR.edgeR = F,
    Seurat.normalize = "CPM",
    Seurat.method = "bimod")
  Seurat<-Pvals_Seurat_full
  end.Seurat<- Sys.time()
  time.Seurat <- difftime(end.Seurat, pred.Seurat, units = "mins")
  cat("Run time for Seurat: ", time.Seurat, "min","\n")
  
  ###########Limma-Voom##################
  pred.Limma <- Sys.time()
  design_full<- model.matrix(~condition)
  nf_full<- calcNormFactors(scData)
  v_full<- voom(scData, design_full, lib.size=colSums(scData)*nf_full,
                normalize.method="quantile", plot=F)
  fit.voom_full<- lmFit(v_full, design_full)
  fit.voom_full<- eBayes(fit.voom_full)
  Limma<-fit.voom_full$p.value[,2]
  end.Limma<- Sys.time()
  time.Limma <- difftime(end.Limma, pred.Limma, units = "mins")
  cat("Run time for Limma: ", time.Limma, "min","\n")
  
  ##############TPMM##############
  pred.TPMM <- Sys.time()
  fData<- data.frame(primerid=rownames(scData))
  #results <- optimcluster(CountData=as.matrix(scDataB), n = 10, seed = 108, Threshold = 0.3, plot = T)
  set.seed(30)
  km.res<- kmeans(t(scData), ncluster, nstart = 10)
  PoolID<-as.vector(km.res$cluster)
  names(condition)=colnames(scData)
  cData<- data.frame(wellKey=colnames(scData),condition=condition,PoolID)
  count<-as.matrix(scData)
  scaRaw<- FromMatrix(count, cData, fData,check_sanity =F)
  
  lmer.output<- zlm(~condition+(1|PoolID), scaRaw, method='bayesglm', ebayes=T)
  lmer.lr<- lrTest(lmer.output, 'condition')
  mixed= cbind(lmer.lr[,,'Pr(>Chisq)'], lmer.lr[,,'lambda'])
  colnames(mixed) = c('gaussian_mixed_pvalue','binomial_mixed_pvalue','combine_mixed_pvalue','gaussian_mixed_lambda','binomial_mixed_lambda','combine_mixed_lambda')
  
  data_TPMM<-data.frame(mixed)
  TPMM<-data_TPMM$combine_mixed_pvalue
  
  end.TPMM<- Sys.time()
  time.TPMM <- difftime(end.TPMM, pred.TPMM, units = "mins")
  cat("Run time for TPMM: ", time.TPMM, "min","\n")

  #######scDEA##############
  pred.scDEA <- Sys.time()
  #condition <- factor(c(rep("1", 229), rep("2", 161)))
  Pvals_scDEA<-scDEA_individual_methods(
    raw.count=scData,
    cell.label=condition,
    is.normalized = FALSE,
    verbose = TRUE,
    BPSC = TRUE,
    DEsingle = TRUE,
    DESeq2 = TRUE,
    edgeR = TRUE,
    MAST = TRUE,
    monocle = TRUE,
    scDD = TRUE,
    Ttest = TRUE,
    Wilcoxon = TRUE,
    limma = TRUE,
    Seurat = TRUE,
    zingeR.edgeR = TRUE,
    BPSC.coef = 2,
    BPSC.normalize = "CPM",
    BPSC.parallel = TRUE,
    DEsingle.parallel = TRUE,
    DEsingle.normalize = "CPM",
    DESeq2.test = "LRT",
    DESeq2.parallel = TRUE,
    DESeq2.beta.prior = TRUE,
    DESeq2.fitType = "parametric",
    DESeq2.normalize = "CPM",
    edgeR.Test = "QLFT",
    edgeR.normalize = "TMM",
    limma.method.fit = "ls",
    limma.trend = TRUE,
    limma.robust = TRUE,
    limma.normalize = "CPM",
    Seurat.normalize = "CPM",
    Seurat.method = "bimod",
    MAST.method = "bayesglm",
    MAST.normalize = "CPM",
    MAST.parallel = TRUE,
    monocle.cores = 1,
    monocle.normalize = "CPM",
    scDD.alpha1 = 0.01,
    scDD.mu0 = 0,
    scDD.s0 = 0.01,
    scDD.a0 = 0.01,
    scDD.b0 = 0.01,
    scDD.normalize = "CPM",
    scDD.permutation = 0,
    Ttest.normalize = "CPM",
    Wilcoxon.normalize = "CPM",
    zingeR.edgeR.normalize = "CPM",
    zingeR.edgeR.maxit.EM = 100
  )
  
  combination.Pvals<- lancaster.combination(Pvals_scDEA, weight = TRUE, trimmed = 0.2)
  adjusted.Pvals<- scDEA.p.adjust(combination.Pvals, adjusted.method = "BH")
  scDEA<-data.frame(Pvals_scDEA,combination.Pvals,adjusted.Pvals)
  
  end.scDEA<- Sys.time()
  time.scDEA <- difftime(end.scDEA, pred.scDEA, units = "mins")
  cat("Run time for scDEA: ", time.scDEA, "min","\n")
  
  #############scHD4E############
  pred.HD4E<- Sys.time()
  Pvals_HN_full<-as.matrix(data.frame(Seurat,Limma,ROSeq,TPMM))
  Pvals_HN_full= Pvals_HN_full%>% replace(is.na(.), 0.20)      # We impute the missing p-values by 0.20 if produce any missing value(rare). This indicates that the genes are non-differentially expressed##
  
  combination.Pvals_HN_full<- lancaster.combination(Pvals_HN_full, weight = TRUE, trimmed = 0.2)
  adjusted.Pvals_HN_full<- scDEA.p.adjust(combination.Pvals_HN_full, adjusted.method = "BH")
  end.HD4E<- Sys.time()
  time.HD4E1<- difftime(end.HD4E, pred.HD4E, units = "mins")
  time.HD4E<-time.Seurat+time.Limma+time.ROSeq+time.HD4E1
  cat("Run time for HD4E: ", time.HD4E, "min","\n")
  
  data_HN_full<-data.frame(Pvals_HN_full,combination.Pvals_HN_full,adjusted.Pvals_HN_full,scDEA) ####Data frame of p-values including our proposed method#############
  return(data_HN_full)
}

# Define server logic
server <- function(input, output) {
  
  observeEvent(input$run_analysis, {
    req(input$file)
    req(input$file1)
    # Read data from CSV file
    scData <- read.csv(input$file$datapath, header = TRUE, row.names = 1)
    
    # Example condition (replace with actual condition)
    #condition <- factor(c(rep("1", 229), rep("2", 161)))
    
    condition1 = read.csv(file = input$file1$datapath,header = T)
    condition = factor(condition1$V1)
    #names(condition) = colnames(scData)
    
    
    # Run scHD4E analysis
    Pvals_full <- scHD4E(scData,condition,ncluster = input$ncluster)
    
    # Show results table
    output$results_table <- renderTable({
      Pvals_full
    })
  })
}

# Run the application
shinyApp(ui = ui, server = server)