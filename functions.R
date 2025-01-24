
tcga_remove_duplicate = function(df) {
  df[,1] = substring(df[,1], 1, 15)
  df = df[!duplicated(df[,1]),]
  return(df)
}

tcga_to_symbol = function(df) {
  
  library(org.Hs.eg.db)
  
  g2s <- toTable(org.Hs.egSYMBOL)
  g2e <- toTable(org.Hs.egENSEMBL)
  merged_id <- merge(g2s, g2e, by='gene_id')
  merged_id <- subset(merged_id, select=-gene_id)
  
  colnames(df)[1] = 'ensembl_id'
  df$ensembl_id <- substr(df$ensembl_id, 1, 15)
  
  df <- merge(df, merged_id, by='ensembl_id')
  df <- df[!duplicated(df$symbol),]
  rownames(df) <- df$symbol
  df <- subset(df, select=-c(symbol, ensembl_id))
  
  return(df)
}

tcga_pam50 = function(expression = read.table('D:/Data/TCGA/BRCA/HiSeqV2', sep='\t', header=T)) {
  
  library(genefu)
  library(org.Hs.eg.db)
  
  rownames(expression) <- expression[,1]
  
  symbols <- expression[,1]##提取基因symbol
  ##id 转换 to entrz
  
  s2g <- toTable(org.Hs.egSYMBOL)
  ids <- s2g[match(symbols, s2g$symbol),1]
  ##构造注释文件
  #  probe Gene.symbol Gene.ID
  id_df <- data.frame(probe = symbols,
                      "Gene.Symbol" = symbols, 
                      "EntrezGene.ID" = ids)
  
  ##保留已注释的基因及注释文件
  expression <- expression[!is.na(id_df$EntrezGene.ID),]
  id_df <- id_df[!is.na(id_df$EntrezGene.ID),] 
  head(id_df)
  
  expression <- expression[,-1]
  expression <- as.data.frame(t(expression))
  
  ##使用genefu包进行分子分型
  data(pam50.robust)
  subtypes <- molecular.subtyping(sbt.model = "pam50", data=expression,
                                  annot=id_df, do.mapping=TRUE)
  
  return(data.frame(id=names(subtypes[[1]]), subtype=subtypes[[1]]))
}


combat_correction = function(a, b) {
  a = as.data.frame(a)
  b = as.data.frame(b)
  
  # filters so that a and b contains the same rownames
  a = a %>% filter(row.names(a) %in% rownames(b))
  b = b %>% filter(row.names(b) %in% rownames(a))
  
  a = as.data.frame(t(a))
  b = as.data.frame(t(b))
  
  data = as.data.frame(t(rbind(a, b)))
  
  a$label = 'a'
  b$label = 'b'
  
  batch = rbind(a, b)
  data_combat = as.data.frame(ComBat(dat=as.matrix(data), batch=batch$label, par.prior=F))
  
  print(nrow(a))
  print(ncol(data_combat))
  
  new_a = data_combat[,c(1:nrow(a))]
  new_b = data_combat[,c((nrow(a)+1):ncol(data_combat))]
  
  return(list(new_a, new_b))
}

# Calls combat_correction
HIM_predict = function(testData, metb) {
  
  library(sva)
  
  meta_df = as.data.frame(metb[[1]])
  rownames(meta_df) = metb[[3]]
  
  # Normalization
  return_combat = combat_correction(testData, meta_df)
  
  # Reconstructing training dataframe
  testData = return_combat[[1]]
  trainData = list(x=as.matrix(return_combat[[2]]),
                   y=metb[[2]],
                   genenames=rownames(return_combat[[2]]),
                   samplelabels=metb[[5]])
  
  # Training
  trainedData <- pamr.train(trainData)
  
  # Predicting subtypes
  classPredict <- pamr.predict(trainedData, as.matrix(testData), threshold=0, type="class")
  probPredict <- pamr.predict(trainedData, as.matrix(testData), threshold=0, type="posterior")
  
  prediction <- cbind(testData$samplelabels, classPredict, probPredict)
  prediction <- as.data.frame(prediction)
  
  prediction$classPredict <- paste0('CS',prediction$classPredict)
  colnames(prediction)[1] <- 'clust'
  
  return(prediction)
}

# Calls combat_correction
subtype_casting = function(testData, him) {
  
  library(sva)
  
  meta_df = as.data.frame(him[[1]])
  
  # Normalization
  return_combat = combat_correction(testData, meta_df)
  
  # Reconstructing training dataframe
  return_combat[[1]] <- na.omit(return_combat[[1]])
  return_combat[[2]] <- na.omit(return_combat[[2]])
  
  testData = return_combat[[1]]
  trainData = list(x=as.matrix(return_combat[[2]]),
                   y=him[[2]],
                   genenames=rownames(return_combat[[2]]),
                   samplelabels=him[[3]])
  
  # Training
  trainedData <- pamr.train(trainData) # requires pamr ver 1.55
  
  # Predicting subtypes
  classPredict <- pamr.predict(trainedData, as.matrix(testData), threshold=0, type="class")
  probPredict <- pamr.predict(trainedData, as.matrix(testData), threshold=0, type="posterior")
  
  prediction <- cbind(testData$samplelabels, classPredict, probPredict)
  prediction <- as.data.frame(prediction)
  
  prediction$classPredict <- paste0('CS',prediction$classPredict)
  colnames(prediction)[1] <- 'clust'
  
  return(prediction)
}


# General purpose limma de
limma_de = function(df, group) {
  library(limma)
  
  vec = colnames(df)
  for (i in c(1:length(vec))) if (vec[i] %in% group) vec[i]='select' else vec[i]='other'
  
  list <- vec %>% factor(., levels = c('select', 'other'), ordered = F)
  
  list <- model.matrix(~factor(list)+0)  # Set group as model matrix
  
  colnames(list) <- c('select', 'other')
  
  df.fit <- lmFit(df, list)
  df.matrix <- makeContrasts(select - other, levels = list)
  
  fit <- contrasts.fit(df.fit, df.matrix)
  fit <- eBayes(fit)
  tempOutput <- topTable(fit, coef=1, n = Inf)
  
  return(tempOutput)
}

# Limma de for HIM subtyping
limma_de_HIM = function(df, pam, group) {
  library(limma)
  
  df = df[,colnames(df) %in% pam$samplelabels]
  
  df = df[,pam$samplelabels]
  
  colnames(df) == pam$samplelabels
  
  vec = pam$subtype
  for (i in c(1:length(vec))) if (vec[i]!=group) vec[i]='other' else vec[i]='select'
  
  list <- vec %>% factor(., levels = c('select', 'other'), ordered = F)
  
  list <- model.matrix(~factor(list)+0)  #把group设置成一个model matrix
  colnames(list) <- c('select', 'other')
  
  df.fit <- lmFit(df, list)
  df.matrix <- makeContrasts(select - other, levels = list)
  
  fit <- contrasts.fit(df.fit, df.matrix)
  fit <- eBayes(fit)
  tempOutput <- topTable(fit, coef=1, n = Inf)
  
  return(tempOutput)
}

gsea_limma = function(res, pathways) {
  
  res2 <- data.frame(symbol=rownames(res), res$t)
  colnames(res2) <- c('SYMBOL', 'stat')
  res2 <- res2[complete.cases(res2),]
  
  res2 <- res2[order(res2$stat, decreasing = T),]
  
  ranks <- deframe(res2)
  
  fgseaRes <- fgsea(pathways=pathways, stats=ranks, nproc=1)
  
  return(fgseaRes)
}

deseq_de = function(df, group) {
  library(DESeq2)
  
  df1 = df[,colnames(df) %in% group]
  df2 = df[,!(colnames(df) %in% group)]
  
  df = cbind(df1, df2)
  
  condition <- factor(c(rep('select',ncol(df1)),rep('other',ncol(df2))))
  colData <- data.frame(row.names=colnames(df), condition)
  
  dds <- DESeqDataSetFromMatrix(countData = round(df), colData = colData, design = ~ condition)
  dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE) 
  
  res <- results(dds1)
  
  return(as.data.frame(res))
}

gsea_deseq = function(res, pathways) {
  
  res2 <- data.frame(symbol=rownames(res), res$stat)
  colnames(res2) <- c('SYMBOL', 'stat')
  res2 <- res2[complete.cases(res2),]
  
  res2 <- res2[order(res2$stat, decreasing = T),]
  
  ranks <- deframe(res2)
  
  fgseaRes <- fgsea(pathways=pathways, stats=ranks, nproc=1)
  
  return(fgseaRes)
}

