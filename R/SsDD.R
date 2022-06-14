SsDD <- function(KTU,kmer,cores=1){
  require(foreach)
  kmer <- kmer[,match(rownames(KTU),colnames(kmer))]
  #distM <- matrix(data = NA,nrow = ncol(KTU),ncol = ncol(KTU))
  cos <- 1-coop::cosine(kmer)
  
  #for(i in 1:ncol(KTU)){
  #  for(j in 1:ncol(KTU)){
  #    us <- c(1:nrow(KTU))[which((KTU[,i]>0)+(KTU[,j]>0)==1)]
  #    #ss <- intersect(which(KTU[,i]>0),which(KTU[,j]>0))
  #    #distM[i,j] <- sum(cos[us,us][lower.tri(cos[us,us])])/sum(cos[ss,ss][lower.tri(cos[ss,ss])])
  #    distM[i,j] <- sum(cos[us,us][lower.tri(cos[us,us])])/sum(cos[lower.tri(cos)])
  #  }
  #}
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  distM <- foreach::foreach(i = 1:ncol(KTU), .combine = rbind, .inorder = T) %:% 
    foreach::foreach(j = 1:ncol(KTU), .combine = rbind, .inorder = T) %dopar%{
      us <- c(1:nrow(KTU))[which((KTU[,i]>0)+(KTU[,j]>0)==1)]
      sum(cos[us,us][lower.tri(cos[us,us])])/sum(cos[lower.tri(cos)])
    }
  distM <- Matrix::sparseMatrix(x = distM[,1],i=rep(1:ncol(KTU),each=ncol(KTU)),j=rep(1:ncol(KTU),ncol(KTU)))
  parallel::stopCluster(cl)  
  colnames(distM) <- rownames(distM) <- colnames(KTU)
  
  bc.dist <- as.matrix(vegan::vegdist(t(KTU)))
  
  distMu <- as.dist(distM)
  distMw <- as.dist(distM * bc.dist)
  return(list(distMu=distMu,distMw=distMw))
}


SsDi <- function(KTU,kmer){
  kmer <- kmer[,match(rownames(KTU),colnames(kmer))]
  ssd <- c()
  for(i in 1:ncol(KTU)){
    kmer1 <- kmer[,match(rownames(KTU)[which(KTU[,i]>0)],colnames(kmer))]
    cos1 <- (1-coop::cosine(kmer1))
    diag(cos1) <- NA
    ssd[i] <- mean(colMeans(cos1,na.rm = T)*(sum(KTU[,i]>0)))
  }
  names(ssd) <- names(KTU)
  return(ssd)
}

tetra.freq2 <- function (repseq, pscore = FALSE, file = TRUE, cores = 1) 
{
  #require(foreach, quietly = T)
  #cl <- parallel::makeCluster(cores)
  #doParallel::registerDoParallel(cl)
  rev_comp <- function(x) {
    x.r <- paste0(rev(strsplit(x, "")[[1]]), collapse = "")
    x.r <- gsub("A", "1", x.r)
    x.r <- gsub("T", "2", x.r)
    x.r <- gsub("C", "3", x.r)
    x.r <- gsub("G", "4", x.r)
    x.r <- gsub("1", "T", x.r)
    x.r <- gsub("2", "A", x.r)
    x.r <- gsub("3", "G", x.r)
    x.r <- gsub("4", "C", x.r)
    return(x.r)
  }
  readseq <- function(seq) {
    x <- seq
    x <- paste0(x, collapse = "")
    x.r <- rev_comp(x)
    return(list(x, x.r))
  }
  if (isTRUE(file)) {
    scaffolds = Biostrings::readDNAStringSet(filepath = repseq, 
                                             use.names = T)
  }
  else scaffolds = Biostrings::DNAStringSet(repseq)
  asv.id <- scaffolds@ranges@NAMES
  species <- lapply(scaffolds, function(x) readseq(as.character(x)))
  DNA <- c("A", "T", "C", "G")
  tetra.mer <- expand.grid(DNA, DNA, DNA, DNA)
  tetra.mer <- do.call(paste0, tetra.mer)
  tetra.table <- matrix(nrow = 256, ncol = length(species))
  rownames(tetra.table) <- tetra.mer
  if (isTRUE(pscore)) {
    for (i in 1:256) {
      for (j in 1:length(species)) {
        single.forward <- ifelse(length(grep(tetra.mer[i], 
                                             species[[j]][[1]])) > 0, grep(tetra.mer[i], 
                                                                           species[[j]][[1]]), 0)
        single.reverse <- ifelse(length(grep(tetra.mer[i], 
                                             species[[j]][[2]])) > 0, grep(tetra.mer[i], 
                                                                           species[[j]][[2]]), 0)
        tetra.table[i, j] <- (single.forward + single.reverse)
      }
    }
  }
  else if (!isTRUE(pscore)) {
#    tetra.table <- foreach(j = 1:length(species), .combine = cbind) %dopar% 
#      {
#        single.forward <- S4Vectors::elementNROWS(Biostrings::matchPDict(Biostrings::PDict(tetra.mer), 
#                                                                         Biostrings::DNAString(unlist(species[[j]][1]))))
#        single.reverse <- S4Vectors::elementNROWS(Biostrings::matchPDict(Biostrings::PDict(tetra.mer), 
#                                                                         Biostrings::DNAString(unlist(species[[j]][2]))))
#        tetras <- (single.forward + single.reverse)
#        tetras
#      }
#    rownames(tetra.table) <- tetra.mer
    single.forward <- sapply(lapply(species, "[[", 1), function(x){
      S4Vectors::elementNROWS(Biostrings::matchPDict(Biostrings::PDict(tetra.mer), 
                                                     Biostrings::DNAString(x)))})
    single.reverse <- sapply(lapply(species, "[[", 2), function(x){
      S4Vectors::elementNROWS(Biostrings::matchPDict(Biostrings::PDict(tetra.mer), 
                                                     Biostrings::DNAString(x)))})
    tetra.table <- (single.forward + single.reverse)
    rownames(tetra.table) <- tetra.mer
  }
  colnames(tetra.table) <- asv.id
  tetra.table <- prop.table(tetra.table, 2)
  return(tetra.table)
}

ssd.rarecurve <- function(ktu,kmer,steps=20,plot=TRUE){
  ssdi <- list()
  for(i in 1:ncol(ktu)){
    subsampsize <- round(seq(0,colSums(ktu)[i],length.out=(steps+1)))
    subssdi <- c()
    for(j in 1:c(steps+1)){
      subdataset <- as.data.frame(table(sample(rep.int(rownames(ktu),ktu[,i]),size = subsampsize[j])),row.names = 1)
      subssdi[j] <- ifelse(length(subdataset)==0,0,SsDi(subdataset,kmer))
    }
    names(subssdi) <- subsampsize
    ssdi[[i]] <- subssdi
  }
  if(isTRUE(plot)){
    plot(seq(1,max(as.numeric(names(unlist(ssdi)))),length.out=steps),seq(0,max(unlist(ssdi)),length.out=steps),
         xlab="Seq_depth",ylab="SSD",las=1,type="n")
    for(i in 1:ncol(ktu)){
      lo <- loess(ssdi[[i]]~as.numeric(names(ssdi[[i]])))
      lines(as.numeric(names(ssdi[[i]])),predict(lo))
    } 
  }
  invisible(ssdi)
}


#KTU::klustering2

##### hind functions #####
rev_comp <- function(x) {
  x.r <- paste0(rev(strsplit(x, "")[[1]]), collapse = "")
  x.r <- gsub("A", "1", x.r)
  x.r <- gsub("T", "2", x.r)
  x.r <- gsub("C", "3", x.r)
  x.r <- gsub("G", "4", x.r)
  x.r <- gsub("1", "T", x.r)
  x.r <- gsub("2", "A", x.r)
  x.r <- gsub("3", "G", x.r)
  x.r <- gsub("4", "C", x.r)
  return(x.r)
}
readseq <- function(seq) {
  x <- seq
  x <- paste0(x, collapse = "")
  x.r <- rev_comp(x)
  return(list(x, x.r))
}
aggregate2df <- function(data, groups, FUN) {
  agg.data <- aggregate(data, list(groups), FUN)
  rownames(agg.data) <- agg.data[, 1]
  agg.data <- as.data.frame(t(agg.data[, -1]))
  return(agg.data)
}
##### hind functions #####
##### update functions #####
tetra.freq <- function (repseq, pscore = FALSE, file = TRUE, output.seq = FALSE, cores=1) 
{
  require(parallel)
  cl <- makeCluster(cores)
  if (isTRUE(file)) {
    scaffolds = Biostrings::readDNAStringSet(filepath = repseq, 
                                             use.names = T)
  }
  else scaffolds = Biostrings::DNAStringSet(repseq)
  asv.id <- scaffolds@ranges@NAMES
  species <- lapply(scaffolds, function(x) readseq(as.character(x)))
  DNA <- c("A", "T", "C", "G")
  tetra.mer <- expand.grid(DNA, DNA, DNA, DNA)
  tetra.mer <- do.call(paste0, tetra.mer)
  tetra.table <- matrix(nrow = 256, ncol = length(species))
  rownames(tetra.table) <- tetra.mer
  if (isTRUE(pscore)) {
    for (i in 1:256) {
      for (j in 1:length(species)) {
        single.forward <- ifelse(length(grep(tetra.mer[i], 
                                             species[[j]][[1]])) > 0, grep(tetra.mer[i], 
                                                                           species[[j]][[1]]), 0)
        single.reverse <- ifelse(length(grep(tetra.mer[i], 
                                             species[[j]][[2]])) > 0, grep(tetra.mer[i], 
                                                                           species[[j]][[2]]), 0)
        tetra.table[i, j] <- (single.forward + single.reverse)
      }
    }
  }
  else if (!isTRUE(pscore)) {
    invisible(clusterEvalQ(cl, library(S4Vectors)))
    invisible(clusterEvalQ(cl, library(Biostrings)))
    single.forward <- parSapply(cl,lapply(species, "[[", 1), function(x){
      elementNROWS(matchPDict(PDict(tetra.mer),DNAString(x)))})
    single.reverse <- parSapply(cl,lapply(species, "[[", 2), function(x){
      elementNROWS(matchPDict(PDict(tetra.mer),DNAString(x)))})
    tetra.table <- (single.forward + single.reverse)
    rownames(tetra.table) <- tetra.mer
  }
  colnames(tetra.table) <- asv.id
  tetra.table <- prop.table(tetra.table, 2)
  stopCluster(cl)
  if(isTRUE(output.seq)){
    return(list(tetra.table=tetra.table,output.seq = as.character(scaffolds)))
  } else return(tetra.table)
}


klustering2 <- function (repseq, seqfromfile=TRUE, pscore = FALSE, feature.table = NULL, write.fasta = TRUE,
                        step = c(5, 10), search.min = NULL, search.max = NULL, cores = 1)
{
  require(parallel)
  oldw <- getOption("warn")
  options(warn = -1)
  if(class(repseq)=="list" & all(names(repseq)==c("tetra.table","output.seq"))){
    message("loading k-mer frequency table...")
    tetra.table <- repseq[[1]]
    asv.id <- colnames(tetra.table)
    species  <- repseq[[2]]
  } else{
    message("k-mer frequency calling...")
    tetra.table <- tetra.freq(repseq, pscore = as.logical(pscore), file = as.logical(seqfromfile), cores=cores) 
    asv.id <- colnames(tetra.table)
    species  <- as.character(Biostrings::readDNAStringSet(filepath = repseq,use.names = T)) #2022/5/18
  }
  message("cosine dissimilarity measurement...")
  cos <- as.dist(1 - coop::cosine(tetra.table))
  tree <- hclust(cos)
  if (is.null(search.max))
    search.max <- max(cutree(tree, h = 0.015))
  else search.max <- search.max
  if (is.null(search.min))
    search.min <- max(cutree(tree, h = 0.03))
  else search.min <- search.min
  message(paste("Searching range:",search.min,"to",search.max))
  if(search.min==search.max | search.max==length(species)){
    message(paste("Best KTU #:", search.max))
    kms <- list()
    kms$clustering <- 1:search.max
    centroid <- names(species)
    crepseq <- species
  } else if(search.min!=search.max){
    first_run_done <- FALSE 
    steps <- ceiling(seq(search.min, search.max, length.out = step[1]))
    cl <- makeCluster(cores) 
    invisible(clusterEvalQ(cl, library(cluster)))
    clusterExport(cl, c("cos"))
    repeat {
      asw <- parLapply(cl,unique(steps), function(x) pam(cos, x, do.swap = F,pamonce = 5)) 
      silhouette <- unlist(sapply(asw, "[[",7)[3,])
      k.best <- steps[order(silhouette, decreasing = T)][1:2]
      steps.update <- ceiling(seq(k.best[2], k.best[1], length.out = step[1]))
      if (all(steps.update == steps))
        break else if(length(unique(steps.update))==1){
          kms <- asw[which.max(silhouette)][[1]]
          message(paste("Best KTU #:", max(kms$clustering)))
          first_run_done <- TRUE 
          break
        } else {message(paste(k.best[1],"to",k.best[2]))
          steps <- steps.update
        }
    } 
    if(!isTRUE(first_run_done)){
      if(length(unique(steps.update))>=2){
        steps <- unique(ceiling(seq(k.best[2], k.best[1], length.out = step[2])))
        repeat {
          asw <- parLapply(cl,unique(steps), function(x) pam(cos, x, do.swap = F,pamonce = 5)) 
          silhouette <- unlist(sapply(asw, "[[",7)[3,])
          if (length(steps) > 1) {
            k.best <- steps[order(silhouette, decreasing = T)][1:2]
          } else k.best <- c(steps[1], steps[1])
          message(paste(k.best[1],"to",k.best[2]))
          steps.update <- ceiling(seq(k.best[2], k.best[1], length.out = step[2]))
          if (all(steps.update == steps)){
            k.best <- steps[which.max(silhouette)]
            message(paste("Best KTU #:", k.best))
            kms <- asw[which.max(silhouette)][[1]]
            break} else if(length(unique(steps.update))==1){
              kms <- asw[[1]]
              message(paste("Best KTU #:", max(kms$clustering)))
              break
            } else {message(paste(k.best[1],"to",k.best[2]))
              steps <- steps.update
            }
        }
      } else if(abs(diff(unique(steps.update)))==1) {
        kms <- asw[which.max(silhouette)][[1]]
        message(paste("Best KTU #:", max(kms$clustering)))
      }
    }
    parallel::stopCluster(cl)
    
    centroid <- kms$id.med
    crepseq <- sapply(species[centroid], "[[", 1)
    for (i in 1:length(centroid)) names(crepseq)[i] <- digest::digest(crepseq[i],algo = "md5", serialize = F)
  }

  kmer.table <- data.frame(tetra.table[, centroid])
  colnames(kmer.table) <- names(crepseq)
  if (isTRUE(write.fasta)) {
    repseq <- matrix(rbind(paste0(">", names(crepseq)),
                           crepseq), ncol = 1)
    write(repseq, file = "ktu-sequence.fasta")
  } else if(is.character(write.fasta)){
    repseq <- matrix(rbind(paste0(">", names(crepseq)),
                           crepseq), ncol = 1)
    write(repseq, file = write.fasta)
  }
  if (!is.null(feature.table)) {
    feature.table <- feature.table[match(asv.id, feature.table[,1], nomatch = 0), ]
    otu <- feature.table[, -1]
    ktu <- as.data.frame(t(aggregate2df(otu, kms$clustering,sum)))
    rownames(ktu) <- names(crepseq)
    colnames(ktu) <- colnames(otu)
    return(list(KTU.table = ktu, ReqSeq = crepseq, kmer.table = kmer.table,
                clusters = kms$clustering))
  }
  else return(list(ReqSeq = crepseq, kmer.table = kmer.table,
                   clusters = kms$clustering))
  options(warn = oldw) 
}

split_tetra <- function(repseq,split_tree_init=5,split_lwrlim=10000,split_reassemble=1000,cores=1){
  #tmpdir <- tempdir()
  tetra <- tetra.freq(repseq,output.seq = TRUE, cores = cores)
  cos <- as.dist(1-coop::cosine(tetra$tetra.table))
  invisible(gc())
  #saveRDS(cos,paste0(tmpdir,"/cos.RDS")
  invisible(gc())
  tree <- hclust(cos)
  
  k=split_tree_init
  message(paste("Initialized cutree:",k))
  ct <- cutree(tree,k=k)
  ct.count <- table(ct)
  while(sum(ct.count>split_lwrlim)>0){
    k=k+1
    #message(paste("cutree:",k))
    ct <- cutree(tree,k=k)
    ct.count <- table(ct)
  }
  ct[ct %in% which(ct.count<=split_reassemble)] <- min(which(ct.count<=split_reassemble))
  ct.count <- table(ct)
  keepname <- names(ct)
  for(i in 1:length(ct.count)) {
    ct <- gsub(pattern = names(ct.count)[i],replacement = i,x = ct)
    ct.count <- table(ct)
  }
  ct <- as.integer(ct)
  names(ct) <- keepname
  invisible(gc())
  message(paste("Finalized cutree:",length(ct.count),"\n subtree sizes:", paste(ct.count,collapse = ", ")))
  
  rm(cos)
  subtetra <- list()
  subtetra.elements <- list()
  for(i in 1:length(ct.count)) {
    subtetra.elements[[1]] <- tetra$tetra.table[,which(ct==i)]
    subtetra.elements[[2]] <- tetra$output.seq[which(ct==i)]
    names(subtetra.elements) <- c( "tetra.table", "output.seq" )
    subtetra[[i]] <- subtetra.elements
  }
  rm(subtetra.elements)
  invisible(gc())
  return(subtetra)
}

trim.primer <- function (fasta, forseq, revseq, truncFL = NULL, output.name = NULL, progress = TRUE) 
{
  fastasets <- Biostrings::readDNAStringSet(fasta)
  seqNames <- fastasets@ranges@NAMES
  Lpattern = Biostrings::DNAString(forseq)
  Forward <- Biostrings::vmatchPattern(pattern = Lpattern, 
                                       subject = fastasets, fixed = F)
  if(!missing(revseq) & is.null(truncFL)){
    Rpattern = Biostrings::reverseComplement(Biostrings::DNAString(revseq))
    Reverse <- Biostrings::vmatchPattern(pattern = Rpattern, 
                                         subject = fastasets, fixed = F)
    trimmed <- c()
    for (i in 1:length(fastasets)) {
      if (length(Forward[[i]]) > 0 & length(Reverse[[i]]) > 
          0) {
        Fstart <- ifelse(length(Forward[[i]]) > 1, Forward[[i]][1]@start, 
                         Forward[[i]]@start)
        Fwidth <- ifelse(length(Forward[[i]]) > 1, Forward[[i]][1]@width, 
                         Forward[[i]]@width)
        Rend <- ifelse(length(Reverse[[i]]) > 1, Reverse[[i]][length(Reverse[[i]])]@start, 
                       Reverse[[i]]@start)
        if (Rend - (Fstart + Fwidth) > 0) {
          trimmed[i] <- as.character(fastasets[[i]][(Fstart + 
                                                       Fwidth):(Rend - 1)])
        }
        else trimmed[i] <- NA
      }
      else trimmed[i] <- NA
      if (progress == TRUE) {
        print(paste0(i, "/", length(fastasets), " sequences from DB are processed"))
        flush.console()
        Sys.sleep(0.01)
      }
    }
  } else if(missing(revseq) & !is.null(truncFL)){
    trimmed <- c()
    for (i in 1:length(fastasets)) {
      if (length(Forward[[i]]) > 0) {
        Fstart <- ifelse(length(Forward[[i]]) > 1, Forward[[i]][1]@start, 
                         Forward[[i]]@start)
        Fwidth <- ifelse(length(Forward[[i]]) > 1, Forward[[i]][1]@width, 
                         Forward[[i]]@width)
        Rend <- truncFL
        if (!is.na(Rend - (Fstart + Fwidth)) & (nchar(as.character(fastasets[[i]])) - (Fstart + Fwidth)) > truncFL) {
          trimmed[i] <- as.character(fastasets[[i]][(Fstart + Fwidth):(Fstart + Fwidth + Rend-1)])
        }
        else trimmed[i] <- NA
      }
      else trimmed[i] <- NA
      if (progress == TRUE) {
        print(paste0(i, "/", length(fastasets), " sequences from DB are processed"))
        flush.console()
        Sys.sleep(0.01)
      }
    }
  }
  names(trimmed) <- seqNames
  trimmed <- trimmed[!is.na(trimmed)]
  output <- matrix(rbind(paste0(">", names(trimmed)), trimmed), 
                   ncol = 1)
  if (is.null(output.name)) {
    write(output, file = "trimmed-sequence.fasta")
  }
  else if (!is.null(output.name)) 
    write(output, file = paste0(output.name, "_trimmed-sequence.fasta"))
}

ktusp <- function (repseq, feature.table = NULL, write.fasta = TRUE,
                   split_tree_init=5,split_lwrlim=10000,split_reassemble=1000,cores = 1, ...)
{
  subtetra <- split_tetra(repseq=repseq,split_tree_init=split_tree_init,split_lwrlim=split_lwrlim,split_reassemble=split_reassemble,cores = cores)
  klusterX <- lapply(subtetra, function(x) klustering2(x,seqfromfile = F,feature.table = feature.table,write.fasta = F,cores = cores, ...))
  
  message("\nSplitted-KTU has done, merging splitted-KTU starts...")
  
  ReqSeq=do.call(c,lapply(klusterX, "[[","ReqSeq"))
  kmer.table=do.call(cbind,lapply(klusterX, "[[","kmer.table"))
  clusterX <- lapply(klusterX, "[[","clusters")
  clusters <- clusterX[[1]]
  for(i in 2:length(clusterX)) clusters <- append(clusters,c(clusterX[[i]]+max(clusters)))
  message(paste0("\n",max(clusters)," KTUs are generated!"))
  
  if (isTRUE(write.fasta)) {
    repseq <- matrix(rbind(paste0(">", names(ReqSeq)),
                           ReqSeq), ncol = 1)
    write(repseq, file = "ktu-sequence.fasta")
  } else if(is.character(write.fasta)){
    repseq <- matrix(rbind(paste0(">", names(ReqSeq)),
                           ReqSeq), ncol = 1)
    write(repseq, file = write.fasta)
  }
  if (!is.null(feature.table)) {
    KTU.table=do.call(rbind,lapply(klusterX, "[[","KTU.table"))
    return(list(KTU.table=KTU.table,ReqSeq=ReqSeq,kmer.table=kmer.table,clusters=clusters))
  }
  else return(list(ReqSeq=ReqSeq,kmer.table=kmer.table,clusters=clusters))
}

##### update functions #####