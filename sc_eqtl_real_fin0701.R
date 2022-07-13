#scRNA-seq analysis
#source("/home/dokada/Dropbox/analysis/2022.5/sc_eqtl_real_fin0701.R") #run 2022.7.5
out_path <- "/home/dokada/work_dir/sc_eqtl_real_fin0701/"
if(file.exists(out_path)){
    unlink(out_path, recursive=TRUE)
    dir.create(out_path)
}else{
    dir.create(out_path)
}

#Processing
data_path <- "/home/dokada/work_dir/bech_data/data_GSE125527/"
files <- sort(list.files(data_path))
n <- length(files)
data_list <- list()
bulk_data <- NULL
for(i in 1:n){
    s <- files[i]
    pbmc.data  <- read.table(paste0(data_path, s),header=T,row.names=1)
    mat <- log1p(as.matrix(pbmc.data))
    mat2 <- apply(mat, 1, function(vec){10^6 * vec/sum(vec)})
    bulk_data <- cbind(bulk_data, rowMeans(mat2))
    data_list[[i]] <-  mat2
    cat(i,"\n")
}
rownames(bulk_data) <- colnames(pbmc.data) #ok
colnames(bulk_data) <- sapply(files, function(x){strsplit(x,"_")[[1]][2]})

#true label
true_label <- sapply(files, function(x){substr(strsplit(x,"_")[[1]][2],1,1)})
gsm_label <-  sapply(files, function(x){strsplit(x,"_")[[1]][1]})


#Bulk data DE and DV analysis
bulk_data <- bulk_data[rowMeans(bulk_data) > 20,]
n_genes <- nrow(bulk_data)
de_p <- de_e <- dv_p <- dv_e <- rep(NA, n_genes)
mean_c <- mean_u <- var_c <- var_u <- rep(NA, n_genes)
for(i in 1:n_genes){
    vec_c <- bulk_data[i, true_label=="C"]
    vec_u <- bulk_data[i, true_label=="U"]
    de_p[i] <- t.test(vec_u, vec_c)$p.value
    de_e[i] <- log2(mean(vec_u)/mean(vec_c))
    dv_p[i] <- var.test(vec_u, vec_c)$p.value
    dv_e[i] <- log2(var(vec_u)/var(vec_c))
    mean_c[i] <- mean(vec_c)
    mean_u[i] <- mean(vec_u)
    var_c[i] <- var(vec_c)
    var_u[i] <-var(vec_u)
}
de_fdr <- p.adjust(de_p, method="BH")
dv_fdr <- p.adjust(dv_p, method="BH")
names(de_p) <- names(dv_p) <- rownames(bulk_data)
names(mean_c) <- names(mean_u) <- names(var_c) <- names(var_u) <- rownames(bulk_data)
names(de_fdr) <- names(dv_fdr) <- names(de_e) <- names(dv_e) <- rownames(bulk_data)

#Cat result
cat("DE gene", sum(de_fdr < 0.05), "\n") #289
cat("DV gene", sum(dv_fdr < 0.05), "\n") #4
intersect(which(dv_fdr < 0.05),  which(de_fdr < 0.05)) #0

#Top DE genes
set.seed(1000)
n_foc <- 1
de_genes <- names(de_p[order(de_p)[1:n_foc]])
dv_genes <- names(dv_p[order(dv_p)[1:n_foc]])
foc_genes <- c(de_genes, dv_genes)
gene_types <- c(rep("DE", n_foc),rep("DV", n_foc))
library(ggplot2)
for(i in 1:length(foc_genes)){

  #Genes
  tmp_gene <- foc_genes[i]
  out_dir <- paste0(out_path, tmp_gene, "/")
  dir.create(out_dir)


  #Bulk gene plot
  vec <- bulk_data[tmp_gene,]
  dat <- data.frame(Group=true_label, Expression=vec)
  g <- ggplot(dat, aes(x=Group, y=Expression))
  g <- g + geom_jitter(height=0, width =0.1) 
  g <- g + stat_summary(fun = mean, geom="point", size =3, color = "red")
  #g <- g + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18))
  #g <- g + theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18))
  g <- g + theme(text = element_text(size = 30))
  plot(g)
  ggsave(file=paste0(out_dir, gene_types[i], ".", tmp_gene, ".png"), plot=g)

  #Clustering
  paras_dat <- NULL
  data_list_gene_C <- list()
  data_list_gene_U <- list()
  for(j in 1:length(data_list)){
    mat2 <- data_list[[j]]
    sc_vec <- mat2[tmp_gene,]
    if(true_label[j] == "C"){
      data_list_gene_C <- c(data_list_gene_C, sc_vec)
    }else{
      data_list_gene_U <- c(data_list_gene_U, sc_vec)
    }
    cls <- kmeans(sc_vec, center=2)$cluster
    mu1 <- mean(sc_vec[cls==1])
    mu2 <- mean(sc_vec[cls==2])
    r1 <- sum(cls==1)/length(cls)
    r2 <- 1 - r1
    mu_plus <- max(mu1, mu2)
    mu_min <- min(mu1, mu2)
    if(mu1 > mu2){
      r <- r1
      cutoff <- max(sc_vec[cls==2])
    }else{
      r <- r2
      cutoff <- max(sc_vec[cls==1])
     }
    paras <- c(mu_plus, mu_min, r)
    names(paras) <- c("mu+", "mu-", "r")
    paras_dat <- rbind(paras_dat, paras)
  }
  write.csv(paras_dat, file=paste0(out_dir, gene_types[i], ".", tmp_gene, ".paras.csv"))

  #Disytibutions of Mu+, Mu-, r
  for(j in 1:ncol(paras_dat)){
    vec <- paras_dat[,j]
    dat <- data.frame(Group=true_label, Expression=vec)
    g <- ggplot(dat, aes(x=Group, y=Expression))
    g <- g + geom_jitter(height=0, width =0.1) 
    g <- g + stat_summary(fun = mean, geom="point", size =3, color = "red")
    #g <- g + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18))
    #g <- g + theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18))
    g <- g + theme(text = element_text(size = 30))
    plot(g)
    ggsave(file=paste0(out_dir, gene_types[i], ".", tmp_gene, ".", colnames(paras_dat)[j],".png"), plot=g)
  }

  #Check
  c_idx <- which(true_label=="C")
  u_idx <- which(true_label=="U")

  #C'sest_e_c
  e_mup_c <- mean(paras_dat[c_idx,1])
  var_mup_c <- var(paras_dat[c_idx,1])
  e_mum_c <- mean(paras_dat[c_idx,2])
  var_mum_c <- var(paras_dat[c_idx,2])
  e_r_c <- mean(paras_dat[c_idx,3])
  var_r_c <- var(paras_dat[c_idx,3])
  est_e_c <- e_mup_c * e_r_c + (1 - e_r_c)*e_mum_c
  est_var_c <- ((e_mup_c - e_mum_c)^2)*var_r_c + (e_r_c^2)*var_mup_c + ((1 - e_r_c)^2)*var_mum_c +var_r_c*(var_mup_c + var_mum_c)

  #U's est_e_c
  e_mup_u <- mean(paras_dat[u_idx,1])
  var_mup_u <- var(paras_dat[u_idx,1])
  e_mum_u <- mean(paras_dat[u_idx,2])
  var_mum_u <- var(paras_dat[u_idx,2])
  e_r_u <- mean(paras_dat[u_idx,3])
  var_r_u <- var(paras_dat[u_idx,3])
  est_e_u <- e_mup_u * e_r_u + (1 - e_r_u)*e_mum_u
  est_var_u <- ((e_mup_u - e_mum_u)^2)*var_r_u + (e_r_u^2)*var_mup_u + ((1 - e_r_u)^2)*var_mum_u + var_r_u*(var_mup_u + var_mum_u)

  #Save, 
  res <- c(e_mup_c, var_mup_c, e_mum_c, var_mum_c, e_r_c, var_r_c, e_mup_u, var_mup_u, e_mum_u, var_mum_u, e_r_u, var_r_u)
  names(res) <- c("e_mup_c", "var_mup_c", "e_mum_c", "var_mum_c", "e_r_c", "var_r_c","e_mup_u", "var_mup_u", "e_mum_u", "var_mum_u", "e_r_u", "var_r_u")
  write.csv(res, file=paste0(out_dir,  gene_types[i], ".", tmp_gene, ".paras.stat.csv"))

  #Save, 
  res <- c(est_e_c, est_var_c, est_e_u, est_var_u)
  names(res) <- c("est_e_c", "est_var_c", "est_e_u", "est_var_u")
  write.csv(res, file=paste0(out_dir,  gene_types[i], ".", tmp_gene, ".est.csv"))


  #Pool
  #tmp <- sample(unlist(data_list_gene_C), 10000)
  #png(paste0(out_dir, gene_types[i], ".", tmp_gene, ".hist.Control.png"), width=960, height=960)
  #hist(unlist(data_list_gene_C), main="Control", cex.axis=4, cex.lab=4, cex.main=4, xlab="", ylab="", freq=F, ylim=c(0, 0.01))
  #dev.off()
  ##tmp <- sample(unlist(data_list_gene_U), 10000)
  #png(paste0(out_dir, gene_types[i], ".", tmp_gene, ".hist.UC.png"), width=960, height=960)
  #hist(unlist(data_list_gene_U), main="UC", cex.axis=4, cex.lab=4, cex.main=4,xlab="", ylab="", freq=F, ylim=c(0, 0.01))
  #dev.off()
  bin_max <- max(c(unlist(data_list_gene_C),unlist(data_list_gene_U))) + 100
  png(paste0(out_dir, gene_types[i], ".", tmp_gene, ".hist.double.png"), width=960, height=960)
  hist(unlist(data_list_gene_C), cex.axis=4, cex.lab=4, cex.main=4, xlab="", ylab="", freq=F, col="#00000020", breaks=seq(from=0,to=bin_max, by=100),ylim=c(0,0.0125), main="")
  hist(unlist(data_list_gene_U), cex.axis=4, cex.lab=4, cex.main=4, xlab="", ylab="", freq=F, add=T, col="#ff000020", breaks=seq(from=0,to=bin_max, by=100),ylim=c(0,0.0125),main="")
  legend("topright", legend = c("C", "U"), col = c("#000000", "#ff0000"), lwd=4, cex=4)
  dev.off()
}

#check
most_de_gene <- names(de_p[order(de_p)[1]]) #POM121
most_de_gene
de_p[most_de_gene]
most_dv_gene <- names(dv_p[order(dv_p)[1]]) #MAP1LC3B2
most_dv_gene
dv_p[most_dv_gene]

#check
parastat <- read.csv(paste0(out_path, "MAP1LC3B2/DV.MAP1LC3B2.paras.stat.csv"), header=T, row.names=1)
e_mup_c   <- parastat["e_mup_c", 1]
var_mup_c <- parastat["var_mup_c", 1]
e_mum_c   <- parastat["e_mum_c", 1]
var_mum_c <- parastat["var_mum_c", 1]
e_r_c     <- parastat["e_r_c", 1]
var_r_c   <- parastat["var_r_c", 1]
e_mup_u   <- parastat["e_mup_u", 1]
var_mup_u <- parastat["var_mup_u", 1]
e_mum_u   <- parastat["e_mum_u", 1]
var_mum_u <- parastat["var_mum_u", 1]
e_r_u     <- parastat["e_r_u", 1]
var_r_u   <- parastat["var_r_u", 1]

#var_r_u's change affect the result
c1 <- ((e_mup_c - e_mum_c)^2)*var_r_c
c2 <- (e_r_c^2)*var_mup_c + ((1 - e_r_c)^2)*var_mum_c
c3 <- var_r_c*(var_mup_c + var_mum_c)
cat("c1:", c1, "\n")
cat("c2:", c2, "\n")
cat("c3:", c3, "\n")
u1 <- ((e_mup_u - e_mum_u)^2)*var_r_u
u2 <- (e_r_u^2)*var_mup_u + ((1 - e_r_u)^2)*var_mum_u
u3 <- var_r_u*(var_mup_u + var_mum_u)
cat("u1:", u1, "\n")
cat("u2:", u2, "\n")
cat("u3:", u3, "\n")
