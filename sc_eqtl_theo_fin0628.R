#theoritical grid (2022.6.28)
#source("/home/dokada/Dropbox/analysis/2022.5/sc_eqtl_theo_fin0628.R") #run 2022.6.30
out_path <- "/home/dokada/work_dir/sc_eqtl_theo_fin0628/"
if(file.exists(out_path)){
    unlink(out_path, recursive=TRUE)
    dir.create(out_path)
}else{
    dir.create(out_path)
}


#3d plot
library(scatterplot3d)
e_mu1 <- 2
e_mu2 <- 1
e_r <- 0.5
v_r <- 0.05
v_mu1 <- 0.3
v_mu2 <- 0.1
d1_grid <- seq(0, 1, by=0.1)
d2_grid <- seq(0, 1, by=0.1)
dr_grid <- c(0.2, 0.5, 0.8)
xyz <- expand.grid(d1_grid, d2_grid, dr_grid)
colnames(xyz) <- c("d1", "d2", "dr")

#mean
scores <- e_r * (xyz[,"dr"] -1) + e_r * xyz[,"dr"] * xyz[,"d1"] + (1 - e_r) * xyz[,"dr"] * xyz[,"d2"]
cols <- ifelse(scores > 0, "red", "blue")
png(paste0(out_path, "mean_theo.png"), width=960, height=960)
par(mar=c(15,15,10,15)) #bot, left, top, roght
par(oma=c(10,8,2,8))
par(mgp=c(15,6,0))
scatterplot3d(xyz[,1], xyz[,2], xyz[,3], color=cols, xlab="d+", ylab="d-", zlab="",cex.axis=3, cex.symbols=4, cex.lab=2)
dev.off()

#variance
v_score <- ((e_mu1 - e_mu2)^2) * v_r + (e_r^2) * v_mu1 + ((1 - e_r)^2) * v_mu2 + v_r * (v_mu1 + v_mu2)
v_score2 <- (((e_mu1+xyz[,"d1"]) - (e_mu2 + xyz[,"d2"]))^2) * v_r + ((xyz[,"dr"]*e_r)^2) * v_mu1 + ((1 - (xyz[,"dr"]*e_r))^2) * v_mu2 + v_r * (v_mu1 + v_mu2)
var_scores <- log2(v_score2/v_score)
cols <- ifelse(var_scores > 0, "red", "blue")
png(paste0(out_path, "var_theo.png"), width=960, height=960)
par(mar=c(15,15,10,15)) #bot, left, top, roght
par(oma=c(10,8,2,8))
par(mgp=c(15,6,0))
scatterplot3d(xyz[,1], xyz[,2], xyz[,3], color=cols, xlab="d+", ylab="d-", zlab="",cex.axis=3, cex.symbols=4, cex.lab=2)
dev.off()
