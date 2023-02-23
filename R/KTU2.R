#' @keywords internal
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
#' @keywords internal
readseq <- function(seq) {
  x <- seq
  x <- paste0(x, collapse = "")
  x.r <- rev_comp(x)
  return(list(x, x.r))
}
#' @keywords internal
aggregate2df <- function(data, groups, FUN) {
  agg.data <- aggregate(data, list(groups), FUN)
  rownames(agg.data) <- agg.data[, 1]
  agg.data <- as.data.frame(t(agg.data[, -1]))
  return(agg.data)
}
#' Generate a kmer table from representative sequences
#'
#' Convert DNA sequences to tetranucleotide frequency table
#' @param repseq The fasta file path
#' @param file logical, default is TRUE, loading DNA sequences from file. if 'file=FALSE', the 'repseq' is loaded from object.
#' @param output.seq logical, default if FALSE, outputting data of tetranucleotide frequency table and sequences.
#' @param cores Numbers of CPUs
#' @return A 256 tetranucleotide combinations frequency table
#' @export
tetra.freq <- function (repseq, pscore = FALSE, file = TRUE, output.seq = FALSE, cores=1)
{
  require(parallel)
  cl <- makeCluster(cores)
  if (isTRUE(file)) {
    if(is.character(repseq)){
      scaffolds = Biostrings::readDNAStringSet(filepath = repseq,
                                               use.names = T)
    } else if(class(repseq)=="DNAStringSet"){
      scaffolds <- repseq
    }
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


#' KTU clustering
#'
#' K-mer taxonomic unit clustering
#' @param repseq The fasta file path
#' @param seqformfile logical, default is TRUE if representative ASV sequences are loaded from a fasta file. seqformfile=FALSE if repseq using pre-calculated tetranucleotide frequency table (with output.seq=TRUE).
#' @param pscore pscore=TRUE/FLASE: using k-mer present scores (tetranucleotide features 0:absent; 1:present in one direction; 2:present in both direction of a sequence)/k-mer frequency
#' @param feature.table (optional) 'data.frame' formated ASV/OTU table. Note: first column must be ID.
#' @param write.fasta (optional) write out a representative KTU sequence fasta file. Logical (output with ktu-sequences.fasta) or characters+.fasta (user-defined names).
#' @param step Split searching range for optimal K. Two step searching process: large scale searching in the first round and smaller scale searching in the second round. Default is 'step=c(5,10)'.
#' @param search.min The minimum K for searching, default is 'NULL' = tip numbers at 0.03 height of cosine hierarchical clustering tree
#' @param search.max The maximum K for searching, default is 'NULL' = tip numbers at 0.015 height of cosine hierarchical clustering tree
#' @param cores Numbers of CPUs
#' @return KTU.table Aggregated KTU table
#' @return ReqSeq Representative KTU sequences
#' @return kmer.table Tetranucleotide present score table or tetranucleotide frequency table
#' @return clusters K-clusters of input features
#' @export
klustering <- function (repseq, seqfromfile=TRUE, pscore = FALSE, feature.table = NULL, write.fasta = TRUE,
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
    if(class(repseq)=="DNAStringSet"){
      species  <- as.character(repseq)
    } else species  <- as.character(Biostrings::readDNAStringSet(filepath = repseq,use.names = T))
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
    #for (i in 1:length(centroid)) names(crepseq)[i] <- digest::digest(crepseq[i],algo = "md5", serialize = F)
    names(crepseq) <- sapply(as.list(crepseq),FUN = function(x) digest::digest(x,algo = "md5", serialize = F))
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



#' Prepare KTU taxonomy assignment database by trimming primers
#'
#' Trim primer sequences of full length 16S/18S/target gene database by PCR primers
#' @param fasta input database fasta file
#' @param forseq forward primer sequence (5' -> 3')
#' @param revseq reverse primer sequence (5' -> 3')
#' @param truncFL trim n-bp length of the sequences from forward primer (for single-end analysis studies)
#' @param output.name (optional) the name for output file
#' @param progress show the processing progress
#' @export
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

#' Build KTU taxonomy assignment database by tetranucleotide frequency table
#'
#' make KTU db
#' @param input.fasta input database fasta file (trimmed)
#' @param input.taxa input database taxonomy file
#' @param pscore pscore=TRUE/FLASE: using k-mer present scores (tetranucleotide features 0:absent; 1:present in one direction; 2:present in both direction of a sequence)/k-mer frequency
#' @param output.file (optional) the directory for output RDS format file
#' @export
makektudb <- function(input.fasta,input.taxa,pscore=FALSE,output.file=NULL){
  rev_comp <- function(x){
    x.r <- paste0(rev(strsplit(x,"")[[1]]),collapse = "")
    x.r <- gsub("A","1",x.r);x.r <- gsub("T","2",x.r);x.r <- gsub("C","3",x.r);x.r <- gsub("G","4",x.r)
    x.r <- gsub("1","T",x.r);x.r <- gsub("2","A",x.r);x.r <- gsub("3","G",x.r);x.r <- gsub("4","C",x.r)
    return(x.r)
  }
  readseq <- function(seq){
    x <- seq
    x <- paste0(x,collapse = "")
    x.r <- rev_comp(x)
    return(list(x,x.r))
  }
  scaffolds = Biostrings::readDNAStringSet(filepath = input.fasta, use.names = T)

  DNA <- c("A","T","C","G")
  tetra.mer <- expand.grid(DNA,DNA,DNA,DNA)
  tetra.mer <- do.call(paste0,tetra.mer)

  species <- lapply(scaffolds, readseq)

  tetra.table <- matrix(nrow=256,ncol=length(species))
  rownames(tetra.table) <- tetra.mer

  if(isTRUE(pscore)){
    for(i in 1:256){
      for(j in 1:length(species)){
        single.forward <- ifelse(length(grep(tetra.mer[i],species[[j]][[1]])) > 0, grep(tetra.mer[i], species[[j]][[1]]),0)
        single.reverse <- ifelse(length(grep(tetra.mer[i],species[[j]][[2]])) > 0, grep(tetra.mer[i], species[[j]][[2]]),0)
        tetra.table[i,j] <- (single.forward+single.reverse)
      }
    }
  } else if (!isTRUE(pscore)) {
    for (j in 1:length(species)) {
      #print(j)
      single.forward <- S4Vectors::elementNROWS(Biostrings::matchPDict(Biostrings::PDict(tetra.mer),Biostrings::DNAString(unlist(species[[j]][1]))))
      single.reverse <- S4Vectors::elementNROWS(Biostrings::matchPDict(Biostrings::PDict(tetra.mer),Biostrings::DNAString(unlist(species[[j]][2]))))
      tetra.table[, j] <- (single.forward + single.reverse)
    }
  }

  colnames(tetra.table) <- scaffolds@ranges@NAMES
  tetra.table <- prop.table(tetra.table,2)

  taxa <- read.delim(input.taxa,header = F)
  taxa <- taxa[match(scaffolds@ranges@NAMES,taxa[,1]),]

  if(!is.null(output.file)){
    saveRDS(tetra.table,file = paste0(output.file,"_DB.RDS"))
    saveRDS(taxa,file = paste0(output.file,"_TX.RDS"))
  }
  return(list(tetra.table,taxa))
}

#' KTU annotation
#'
#' KTU-based alignment-free taxonomy assignment
#' @param dbRDS Kmer formated database RDS file or choosing a built-in database. see \link[KTU2]{kaxaDB}
#' @param taxaRDS taxonomy ID RDS file, NULL if dbRDS is used by built-in database
#' @param kmer.table tetranucleotide frequency table for annotation
#' @param cos.cutff cosine similarity cutoff, default=0.95
#' @param consensus taxonomy assignment by consensus [0-1]
#' @param candidate how many candidates (≥ cos.cutff) for taxonomy assignment
#' @param cores Numbers of CPUs
#' @export
kaxonomy <- function(dbRDS,taxaRDS=NULL,kmer.table,cos.cutff=0.95,consensus=0.8,candidate=10,cores=1){
  require(foreach,quietly = T)
  machine.core <- function(x){
    xxx <- matrix(ncol = 8)
    xx <- cbind(taxa[which(x>=cos.cutff),],cos.score=x[which(x>=cos.cutff)])
    if(nrow(xx)>0){
      ranks <- rank(1/(xx$cos.score))
      if(min(ranks)>candidate){
        xx <- xx[which(ranks==min(ranks)),]
      } else xx <- xx[which(ranks<=candidate),]
      xx[,-8] <- apply(xx[,-8], 2, as.character)
      xy <- apply(xx,2,function(x) names(table(x))[which(prop.table(table(x))>consensus)])
      xxx[,1:7] <- c(xy[[1]][1],xy[[2]][1],xy[[3]][1],xy[[4]][1],xy[[5]][1],xy[[6]][1],xy[[7]][1])
      #xxx[,8] <- mean(unlist(apply(xx,2,function(x) prop.table(table(x))[which(prop.table(table(x))>consensus)])))
      xxx[,8] <- mean(mean(xx[,8]) * unlist(apply(xx,2,function(x) prop.table(table(x))[which(prop.table(table(x))>consensus)])))
      xxx[,7] <- ifelse(!is.na(xxx[,7]) & is.na(xxx[,5]),NA,xxx[,7])
      for(i in 7:2) xxx[,i] <- ifelse(!is.na(xxx[,i]) & is.na(xxx[,i-1]),NA,xxx[,i])
    } else if(nrow(xx)==0) xxx[,c(1,8)] <- c("Unassigned",1)
    return(xxx)
  }
  cl <- parallel::makeCluster(cores) #not to overload your computer
  doParallel::registerDoParallel(cl)

  if(is.null(taxaRDS) & dbRDS=="NCBI_V3V4" | dbRDS=="NCBI_V4" | dbRDS=="SILVA138_V3V4" | dbRDS=="SILVA138_V4"){
    path <- eval(parse(text = 'system.file("extdata",package = "KTU2")'))
    db.tetra <- readRDS(paste0(path,"/",dbRDS,"_DB.RDS"))
    db.taxa <- readRDS(paste0(path,"/",dbRDS,"_TX.RDS"))
  } else if(!is.null(dbRDS) & !is.null(taxaRDS)){
    db.tetra <- readRDS(dbRDS)
    db.taxa <- readRDS(taxaRDS)
  }
  {
    taxa <- as.character(db.taxa[,2])
    sep <- ifelse(grepl("; ",taxa[1],fixed = T),"; ",";")
    taxa <- strsplit(taxa,sep)
    taxa <- data.frame(do.call(rbind,taxa),row.names = db.taxa[,1])
  }

  i1 <- seq(1,ncol(data.frame(kmer.table)),20)
  i2 <- c((i1[-1]-1),ncol(data.frame(kmer.table)))
  i2 <- i2[i2>0]

  for(i in 1:length(i1)){
    kmer <-  kmer.table[,i1[i]:i2[i]]
    cos.table <- data.frame(matrix(data = NA,nrow = nrow(db.taxa),ncol=ncol(kmer)))
    cos.table <- foreach::foreach(cln=1:ncol(kmer), .combine=cbind) %dopar% {
      temp.cos.table = apply(db.tetra,2,function(x) coop::cosine(kmer[,cln],x)) #calling a function
      temp.cos.table #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
    }
    taxa.consensus <- lapply(X = as.list(as.data.frame(cos.table)),
                             FUN = machine.core)
    taxa.consensus <- data.frame(do.call(rbind,taxa.consensus))
    colnames(taxa.consensus) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Cos.score")
    taxa.consensus <- apply(taxa.consensus, 2, as.character)
    taxa.consensus[is.na(taxa.consensus)] <- "Unassigned"
    taxa.consensus <- as.data.frame(taxa.consensus)
    taxa.consensus[,8] <- as.numeric(as.character(taxa.consensus[,8]))
    rownames(taxa.consensus) <- colnames(kmer)
    invisible(gc())
    if(i==1){
      write.table(taxa.consensus,file = "kaxa-temp.tsv",sep = "\t")
    } else if(i>1){
      write.table(taxa.consensus,file = "kaxa-temp.tsv",append = T,sep = "\t",col.names = FALSE)
    }
  }
  parallel::stopCluster(cl) #stop cluster
  taxa.consensus <- read.delim("kaxa-temp.tsv")
  file.remove("kaxa-temp.tsv")
  return(taxa.consensus)
}



#' KTU evaluation
#'
#' Evaluate sequence similarity within KTUs
#' @param klusterRDS klustering output RDS file
#' @param ASVfasta representative ASV sequence fasta file
#' @export
KTUsim.eval <- function(klusterRDS,ASVfasta){
  kluster <- readRDS(klusterRDS)
  dna <- Biostrings::readDNAStringSet(ASVfasta,use.names = T)
  dna <- lapply(dna, as.character)
  dna <- sapply(dna, "[[", 1)
  dna <- dna[match(names(kluster$clusters),names(dna))]

  mean.simi <- c()
  Ks <- c()
  kn.stats <- table(kluster$clusters)
  for(k in 1:length(kn.stats)){
    if(kn.stats[k]>1){
      dna.cluster.n <- dna[which(kluster$clusters==k)]
      dna.cluster.cb <- combn(x = 1:length(dna.cluster.n),m = 2)
      dna.seq.cb <- apply(dna.cluster.cb,2,function(x) dna.cluster.n[x])
      within.simi <- c()
      for(i in 1:ncol(dna.seq.cb)) within.simi[i] <- Biostrings::pid(Biostrings::pairwiseAlignment(dna.seq.cb[1,i],dna.seq.cb[2,i]))
      mean.simi[k] <- mean(within.simi)

      tetra <- tetra.freq(dna[which(kluster$clusters==k)],file = F)
      Ks[k] <- mean(as.dist(1-coop::cosine(tetra)))
      Ks[k] <- ifelse(Ks[k]>0,Ks[k],0) # check no negative values
    } else if(kn.stats[k]==1){
      mean.simi[k] <- 100
      Ks[k] <- NA
    }
    print(paste("mean similarity within Kluster",k,"=",mean.simi[k],"% from",kn.stats[k],"seq features; divergence =",Ks[k]))
  }
  g.mean.simi <- mean(mean.simi)
  g.Ks <- mean(Ks[!is.na(Ks)])
  print(paste("Mean similarity of 'within KTU' of", length(kn.stats), "KTUs =",g.mean.simi,"; Mean divergence within KTU =",g.Ks))
  #hist(mean.simi)
  return(list(eachmean=data.frame(kluster=1:length(kn.stats),similarity.pct=mean.simi,divergence=Ks,n.feature=as.data.frame(kn.stats)$Freq),globalmean=g.mean.simi,globaldivergence=g.Ks))
}


#' plot KTU evaluations
#'
#' Evaluate sequence similarity within KTUs by histogram
#' @param histdata KTU evaluation variable
#' @param coll color for distribution curve
#' @export
histNdistri <- function(histdata,coll,xlab.text=NULL,legend.posit="topright",legend.text="",...){
  h <- hist(histdata,main="",xlab=xlab.text,...)
  xfit<-seq(min(histdata),max(histdata),length=100)
  yfit<-dnorm(xfit,mean=mean(histdata),sd=sd(histdata))
  yfit <- yfit*diff(h$mids[1:2])*length(histdata)
  lines(xfit,yfit,type="l",lwd=2,col=coll,xlab="",ylab="",xaxt="n",yaxt="n",bty="n",xpd=T)
  legend(x=legend.posit,legend = legend.text,bty="n")
}

#' Subset large tetranucleotide frequency table
#'
#' Subset large tetranucleotide frequency table by splitting kmer tree
#' @param repseq The fasta file path
#' @param split_tree_init Initial numbers for cutting kmer tree, default is 5
#' @param split_lwrlim Split kmer subtrees if tree tips > n, n default is 10000
#' @param split_reassemble Reassemble tree tips when subtrees are too small, default is 1000
#' @param cores Numbers of CPUs
#' @return Splitted tetranucleotide frequency tables
#' @export
split_tetra <- function(repseq,split_tree_init=5,split_lwrlim=10000,split_reassemble=1000,cores=1){
  tetra <- tetra.freq(repseq,output.seq = TRUE, cores = cores)
  cos <- as.dist(1-coop::cosine(tetra$tetra.table))
  invisible(gc())
  invisible(gc())
  tree <- hclust(cos)

  k=split_tree_init
  message(paste("Initialized cutree:",k))
  ct <- cutree(tree,k=k)
  ct.count <- table(ct)
  while(sum(ct.count>split_lwrlim)>0){
    k=k+1
    ct <- cutree(tree,k=k)
    ct.count <- table(ct)
  }
  ct[ct %in% which(ct.count<=split_reassemble)] <- min(which(ct.count<=split_reassemble))
  ct.count <- table(ct)
  keepname <- names(ct)
  for(i in 1:length(ct.count)) {
    ct <- gsub(pattern = names(ct.count)[i],replacement = i,x = ct)
    ct <- as.integer(ct)
    ct.count <- table(ct)
  }
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

#' Splitting KTU algorithm
#'
#' Implement KTU algorithm with subsisting a large tetranucleotide frequency table by splitting kmer tree
#' @param repseq The fasta file path
#' @param feature.table (optional) 'data.frame' formated ASV/OTU table. Note: first column must be ID.
#' @param write.fasta (optional) write out a representative KTU sequence fasta file. Logical (output with ktu-sequences.fasta) or characters+.fasta (user-defined names).
#' @param split_tree_init Initial numbers for cutting kmer tree, default is 5
#' @param split_lwrlim Split kmer subtrees if tree tips > n, n default is 10000
#' @param split_reassemble Reassemble tree tips when subtrees are too small, default is 1000
#' @param cores Numbers of CPUs
#' @param ... other parameters of `klustering` function
#' @return KTU.table Aggregated KTU table
#' @return ReqSeq Representative KTU sequences
#' @return kmer.table Tetranucleotide present score table or tetranucleotide frequency table
#' @return clusters K-clusters of input features
#' @export
ktusp <- function (repseq, feature.table = NULL, write.fasta = TRUE,
                   split_tree_init=5,split_lwrlim=10000,split_reassemble=1000,cores = 1, ...)
{
  subtetra <- split_tetra(repseq=repseq,split_tree_init=split_tree_init,split_lwrlim=split_lwrlim,split_reassemble=split_reassemble,cores = cores)
  klusterX <- lapply(subtetra, function(x) klustering(x,seqfromfile = F,feature.table = feature.table,write.fasta = F,cores = cores, ...))

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

#' Klustering from dada2 workflow
#'
#' It reads the output matrix from dada2::makeSequenceTable into KTU::klustering or KTU::ktusp.
#'
#' @param seqtab matrix, ouput from dada2::makeSequenceTable.
#' @param subset (optional) 'numeric' for debugging. Limit klustering to a randomly sampled subset of ASVs/OTUs.
#' @param method either 'klustering' or 'ktusp'.
#' @param split_reassemble Reassemble tree tips when subtrees are too small, default is 1000 (see function 'ktusp').
#' @param ...	other parameters passed to 'klustering' or 'ktusp'.
#' @return KTU.table Aggregated KTU table.
#' @return ReqSeq Representative KTU sequences
#' @return kmer.table Tetranucleotide present score table or tetranucleotide frequency table.
#' @return clusters K-clusters of input features.
#' @export
#' @examples
#' data(seqtab)
#' # run klustering on a subset of the ASVs.
#' dada2KTU(seqtab = seqtab, subset = 200, cores = 1)
#' # run ktusp:
#' dada2KTU(seqtab = seqtab, method = "ktusp", subset = 800, cores = 2,
#'          split_tree_init = 5, split_lwrlim = 10000, split_reassemble = 200)
dada2KTU <- function(seqtab = NULL, subset = NULL, method = "klustering", split_reassemble = 1000, ...) {
  require(Biostrings)
  #feature table (ie OTU table) from phyloseq
  stopifnot(any(class(seqtab) %in% "matrix"))
  ft <- data.frame(t(seqtab))
  seqs <-
    {
      hh <- Biostrings::DNAStringSet(dimnames(seqtab)[[2]])
      names(hh) <- as.character(hh)
      hh
    }
  # subset (optional)
  if(is.numeric(subset)){
    seqs <- sample(seqs, subset)
    ft <- ft[match(as.character(seqs), rownames(ft)), ]
  }
  # checks
  if(class(seqs) != "DNAStringSet" | length(seqs) != nrow(ft))
    stop("Input must have sequences as colnames.")
  stopifnot(class(ft) == "data.frame")
  stopifnot(!is.null(names(seqs)))
  # add rownames to first column
  ft <- cbind(asv = rownames(ft), ft)
  # write sequences to temp fasta file
  tempseqs <- tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(seqs, tempseqs)
  # run clustering
  if(method == "klustering") {
    klustering(repseq = tempseqs,
                    feature.table = ft, ...)
  } else if(method == "ktusp") {
    if(nrow(ft) < split_reassemble)
      warning("Number of tips is lower than the threshold to reassamble tips. Consider tuning ktusp or using klustering")
    ktusp(repseq = tempseqs,
          feature.table = ft, ...)
  }
}
#' Klustering using as input a phyloseq object
#'
#' It reads a phyloseq object with the refseq slot.
#'
#' @param ps phyloseq object.
#' @param subset (optional) 'numeric' for debugging. Limit klustering to a randomly sampled subset of ASVs/OTUs.
#' @param dnaseqs otu_table/refseq: slot to search for DNA sequences in phyloseq; either as row/col names in otu_table or as a DNAStringSet in the refseq.
#' @param method either 'klustering' or 'ktusp'.
#' @param split_reassemble Reassemble tree tips when subtrees are too small, default is 1000 (see function 'ktusp').
#' @param ...	other parameters passed to 'klustering' or 'ktusp'.
#' @return KTU.table Aggregated KTU table.
#' @return ReqSeq Representative KTU sequences
#' @return kmer.table Tetranucleotide present score table or tetranucleotide frequency table.
#' @return clusters K-clusters of input features.
#' @export
#' @examples
#' ps <- readRDS(
#' url("https://zenodo.org/record/4155877/files/Schreuder2020_TemporalDynamicsCloacalMicrobiota_phyloseq.rds"))
#' # run klustering on a subset of the ASVs with DNA sequences as ASV name in otu_table:
#' ps2KTU(ps = ps, subset = 200, dnaseqs = "otu_table")
#' # run klustering on a subset of the ASVs with DNA sequences in refseq slot:
#' data(phyloseq_marsh)
#' ps2KTU(ps = phyloseq_marsh, subset = 200, dnaseqs = "refseq")
#' # run ktusp:
#' ps2KTU(ps = ps, method = "ktusp", subset = 800, cores = 2, dnaseqs = "otu_table",
#'          split_tree_init = 5, split_lwrlim = 10000, split_reassemble = 200)
ps2KTU <- function(ps = NULL, subset = NULL, dnaseqs = NULL, method = "klustering", split_reassemble = 1000, ...) {
  require(Biostrings); require(phyloseq)
  #feature table (ie OTU table) from phyloseq
  stopifnot(any(class(ps) %in% "phyloseq"))
  otu1 <- as(otu_table(ps), "matrix")
  # transpose if necessary: rows are taxa
  if(!taxa_are_rows(ps)) otu1 <- t(otu1)
  ft <- as.data.frame(otu1)
  # get sequences
  if(dnaseqs == "otu_table") {
    seqs <- rownames(otu1)
    seqs <- Biostrings::DNAStringSet(seqs)
    names(seqs) <- as.character(seqs)
  } else if(dnaseqs == "refseq") {
    seqs <- phyloseq::refseq(ps)
    rownames(ft) <- seqs[match(rownames(ft), names(seqs))]
    names(seqs) <- as.character(seqs)
  }
  if(class(seqs) != "DNAStringSet" | length(seqs) != ntaxa(ps))
    stop("phyloseq object must have DNA sequences as rownanes in the otu_table slot or as a DNAStringSet in the refseq slot")
  stopifnot(class(ft) == "data.frame")
  stopifnot(!is.null(names(seqs)))
  # add rownames to first column
  ft <- cbind(asv = rownames(ft), ft)
  # subset (optional)
  if(is.numeric(subset)){
    seqs <- sample(seqs, subset)
    ft <- ft[match(as.character(seqs), rownames(ft)), ]
  }
  # write sequences to temp fasta file
  tempseqs <- tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(seqs, tempseqs)
  # run clustering
  if(method == "klustering") {
    KTU::klustering(repseq = tempseqs,
                    feature.table = ft, ...)
  } else if(method == "ktusp") {
    if(nrow(ft) < split_reassemble)
      warning("Number of tips is lower than the threshold to reassamble tips. Consider tuning ktusp or using klustering")
    ktusp(repseq = tempseqs,
          feature.table = ft, ...)
  }
}

#' Subsampled ASVs abundance table from metabarcoding of 16S from marsh soils
#'
#' This dataset contains subsampled most frequent ASVs and their counts as a matrix across 41 samples
#'  from marsh soils in southern Spain.
#'  They come from: Camacho-Sanchez et al. Bacterial assemblage in Mediterranean salt marshes: disentangling
#'  the relative importance of seasonality, zonation and halophytes. In review.
#'
#' @name seqtab
#' @usage data(seqtab)
#' @format A matrix containing 1000 columns, ASVS, and 41 rows, samples.
#' @source \url{https://github.com/csmiguel/marsh_metabarcoding}
"seqtab"

#' Subsampled phyloseq object from metabarcoding of 16S from marsh soils
#'
#' This dataset contains subsampled most frequent ASVs and their counts as a matrix across 41 samples
#'  from marsh soils in southern Spain.
#'  They come from: Camacho-Sanchez et al. Bacterial assemblage in Mediterranean salt marshes: disentangling
#'  the relative importance of seasonality, zonation and halophytes. In review.
#'
#' @name phyloseq_marsh
#' @usage data(phyloseq_marsh)
#' @format A phyloseq object containing 2417 ASVs and 10 samples.
#' @source \url{https://github.com/csmiguel/marsh_metabarcoding}
"phyloseq_marsh"

#' Pre-trained 16S Databases for Kaxonomy annotation
#'
#' The pre-trained databases contain a species+strain levels 16S rRNA gene database retrieved from NCBI (2022/02/10) and SILVA 7-level 16S rRNA gene database (ver. 138).
#'  Two hypervariable-region trimmed formats, V3V4 (341F-805R, primer-trimmed) and V4 (515F-806R, primer-trimmed) regions, are available.
#'  A pair of files (DB and TX, RDS format) are required for annotation by using 'kaxonomy' function.
#'  See how to use a pre-trained database for taxonomy annotation in Examples.
#' Available built-in databases: NCBI_V3V4, NCBI_V4, SILVA138_V3V4, and SILVA138_V4
#' @name kaxaDB
#' @usage NULL
#' @format Pre-trained tetranucleotide database and taxonomy information (path to RDS files) of 16S rRNA genes.
#' @examples
#'
#' kluster <- readRDS("kluster.RDS") # read klustering/ktusp result object from RDS file
#' kaxa <- kaxonomy(dbRDS = "NCBI_V3V4",
#'                  taxaRDS = NULL,
#'                  kmer.table = kluster$kmer.table,
#'                  cores = 4)
#'
#' write.table(kaxa,file = "kaxonomy.tsv",sep="\t") # write out with tab-delimited text file
#' saveRDS(kaxa,"kaxonomy.RDS") # write out with RDS format
kaxaDB <- function(){}
