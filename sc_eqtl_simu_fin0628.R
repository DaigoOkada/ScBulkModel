#Simulation for DVG(2022.6.28)
#source("/home/dokada/Dropbox/analysis/2022.5/sc_eqtl_simu_fin0628.R") #run 2022.7.5
#mean is equal, variance is different
out_path <- "/home/dokada/work_dir/sc_eqtl_simu_fin0628/"
if(file.exists(out_path)){
    unlink(out_path, recursive=TRUE)
    dir.create(out_path)
}else{
    dir.create(out_path)
}

#Simulation1
set.seed(1000)
n_cells <- 10000
sc_v1_fix <- 1
sc_v2_fix <- 1
v_err <- 0.1
e_mu1 <- 2
e_mu2 <- 1
e_r <- 0.9
v_r <- 0.001
v_mu1 <- 1 
v_mu2 <- 0.1 
N <- 1000

#Data generation(Control)
mu1_ind <- rnorm(N, e_mu1, sqrt(v_mu1))
mu2_ind <- rnorm(N, e_mu2, sqrt(v_mu2))
r_ind <- runif(N, e_r - sqrt(12*v_r)/2, e_r + sqrt(12*v_r)/2)
bulk <- rep(NA, N)
for(i in 1:N){
  mu1 <- mu1_ind[i]
  mu2 <- mu2_ind[i]
  r <- r_ind[i]
  n_subset1 <- rbinom(1,prob=r, size=n_cells)
  n_subset2 <- n_cells - n_subset1
  sc_vec <- c(rnorm(n_subset1, mu1, sqrt(sc_v1_fix)), rnorm(n_subset2, mu2, sqrt(sc_v2_fix)))
  bulk[i] <- mean(sc_vec) + rnorm(1, 0, sqrt(v_err))
}

#Check
theo_bulk_mean <- e_mu1 * e_r + (1 - e_r)*e_mu2
theo_bulk_var <- ((e_mu1 - e_mu2)^2) * v_r + (e_r^2) * v_mu1 + 
((1 - e_r)^2) * v_mu2 + v_r * (v_mu1 + v_mu2) + v_err

#Data generation(Case)
d1 <- 0.5
d2 <- 0.5
dr <- (e_r - d2)/(e_r*(1 + d1 - d2))
e_mu1_case <- e_mu1 + d1*(e_mu1 - e_mu2)
e_mu2_case <- e_mu2 + d2*(e_mu1 - e_mu2)
e_r_case <- dr * e_r
mu1_ind <- rnorm(N, e_mu1_case, sqrt(v_mu1))
mu2_ind <- rnorm(N, e_mu2_case, sqrt(v_mu2))
r_ind <- runif(N, e_r_case - sqrt(12*v_r)/2, e_r_case + sqrt(12*v_r)/2)
bulk_case <- rep(NA, N)
for(i in 1:N){
  mu1 <- mu1_ind[i]
  mu2 <- mu2_ind[i]
  r <- r_ind[i]
  n_subset1 <- rbinom(1,prob=r, size=n_cells)
  n_subset2 <- n_cells - n_subset1
  sc_vec <- c(rnorm(n_subset1, mu1, sqrt(sc_v1_fix)), rnorm(n_subset2, mu2, sqrt(sc_v2_fix)))
  bulk_case[i] <- mean(sc_vec) + rnorm(1, 0, sqrt(v_err))
}

theo_bulk_mean_case <- e_mu1_case * e_r_case + (1 - e_r_case)*e_mu2_case
theo_bulk_var_case <- ((e_mu1_case - e_mu2_case)^2) * v_r + (e_r_case^2) * v_mu1 + 
((1 - e_r_case)^2) * v_mu2 + v_r * (v_mu1 + v_mu2) + v_err

#Draw the original single cell distributions
cex.axis <- 4
cex.lab <- 4
set.seed(1000)
sc_v1_fix
sc_v2_fix 
e_mu1
e_mu2 
e_r
e_mu1_case
e_mu2_case
e_r_case
x_grid <- seq(-2, 8, length=1000)
pdf_cont <- e_r * dnorm(x_grid, e_mu1, sqrt(sc_v1_fix)) + (1 - e_r) * dnorm(x_grid, e_mu2, sqrt(sc_v2_fix))
pdf_case <- e_r_case * dnorm(x_grid, e_mu1_case, sqrt(sc_v1_fix)) + (1 - e_r_case) * dnorm(x_grid, e_mu2_case, sqrt(sc_v2_fix))
png(paste0(out_path, "simu1_dist.png"), width=960, height=960)
#par(mar=c(15,15,10,15)) #bot, left, top, roght
#par(oma=c(10,8,2,8))
#par(mgp=c(15,6,0))
matplot(x_grid, cbind(pdf_cont, pdf_case), col=c("black", "red"), type="l", lty=1, xlab="", ylab="", cex.axis=cex.axis, cex.lab=cex.lab, lwd=4)
dev.off()


#Draw the original parameters (r)
r_grid <- seq(0, 1, length=1000)
r_pdf_cont <- dunif(r_grid, e_r - sqrt(12*v_r)/2, e_r + sqrt(12*v_r)/2)
r_pdf_case <- dunif(r_grid, e_r_case - sqrt(12*v_r)/2, e_r_case + sqrt(12*v_r)/2)
png(paste0(out_path, "simu1_r.png"), width=960, height=960)
#par(mar=c(10,10,2,2)) #bot, left, top, roght
#par(oma=c(10,10,2,2))
#par(mgp=c(8,8,0))
matplot(r_grid, cbind(r_pdf_cont, r_pdf_case), col=c("black", "red"), type="l", lty=1, xlab="", ylab="", cex.axis=cex.axis, cex.lab=cex.lab, lwd=4)
dev.off()

#Draw the original parameters (mu)
x_min <- min(e_mu1 - 3*sqrt(v_mu1), e_mu1_case - 3*sqrt(v_mu1), e_mu2 - 3*sqrt(v_mu2), e_mu2_case - 3*sqrt(v_mu2))
x_max <- max(e_mu1 + 3*sqrt(v_mu1), e_mu1_case + 3*sqrt(v_mu1),e_mu2 + 3*sqrt(v_mu2), e_mu2_case + 3*sqrt(v_mu2))
x_grid <- seq(x_min, x_max, length=1000)
mu1_pdf_cont <- dnorm(x_grid, e_mu1, sqrt(v_mu1))
mu1_pdf_case <- dnorm(x_grid, e_mu1_case, sqrt(v_mu1))
mu2_pdf_cont <- dnorm(x_grid, e_mu2, sqrt(v_mu2))
mu2_pdf_case <- dnorm(x_grid, e_mu2_case, sqrt(v_mu2))
png(paste0(out_path, "simu1_mu.png"), width=960, height=960)
#par(mar=c(5,5,10,10)) #bot, left, top, roght
#par(oma=c(4,4,2,2))
#par(mgp=c(10,4,0))
matplot(x_grid, cbind(mu1_pdf_cont, mu1_pdf_case, mu2_pdf_cont, mu2_pdf_case), type="l", col=c("black", "red", "black", "red"), lty=c(1, 1, 2, 2), xlab="", ylab="", cex.axis=cex.axis, cex.lab=cex.lab, lwd=4)
dev.off()

#Check
mean(bulk)
var(bulk)
print(theo_bulk_mean)
print(theo_bulk_var)
mean(bulk_case)
print(theo_bulk_mean_case)
var(bulk_case)
print(theo_bulk_var_case)
var.test(bulk, bulk_case)
t.test(bulk, bulk_case)
png(paste0(out_path, "simu1_boxplot.png"), width=960, height=960)
#par(mar=c(15,15,10,15)) #bot, left, top, roght
#par(oma=c(10,8,2,8))
#par(mgp=c(15,6,0))
boxplot(list(C=bulk, D=bulk_case), xlab="", ylab="", cex.axis=cex.axis, cex.lab=cex.lab)
dev.off()


#Simulation2
set.seed(1000)
n_cells <- 10000
sc_v1_fix <- 1
sc_v2_fix <- 1
v_err <- 0.1
e_mu1 <- 2
e_mu2 <- 1
e_r <- 0.9
v_r <- 0.001
v_mu1 <- 0.01 
v_mu2 <- 0.01 

#Data generation(Control)
mu1_ind <- rnorm(N, e_mu1, sqrt(v_mu1))
mu2_ind <- rnorm(N, e_mu2, sqrt(v_mu2))
r_ind <- runif(N, e_r - sqrt(12*v_r)/2, e_r + sqrt(12*v_r)/2)
bulk <- rep(NA, N)
for(i in 1:N){
  mu1 <- mu1_ind[i]
  mu2 <- mu2_ind[i]
  r <- r_ind[i]
  n_subset1 <- rbinom(1,prob=r, size=n_cells)
  n_subset2 <- n_cells - n_subset1
  sc_vec <- c(rnorm(n_subset1, mu1, sqrt(sc_v1_fix)), rnorm(n_subset2, mu2, sqrt(sc_v2_fix)))
  bulk[i] <- mean(sc_vec) + rnorm(1, 0, sqrt(v_err))
}

#Check
theo_bulk_mean <- e_mu1 * e_r + (1 - e_r)*e_mu2
theo_bulk_var <- ((e_mu1 - e_mu2)^2) * v_r + (e_r^2) * v_mu1 + 
((1 - e_r)^2) * v_mu2 + v_r * (v_mu1 + v_mu2) + v_err


#Data generation(Case)
d1 <- 10
d2 <- 0.1
dr <- (e_r - d2)/(e_r*(1 + d1 - d2))
e_mu1_case <- e_mu1 + d1*(e_mu1 - e_mu2)
e_mu2_case <- e_mu2 + d2*(e_mu1 - e_mu2)
e_r_case <- dr * e_r
mu1_ind <- rnorm(N, e_mu1_case, sqrt(v_mu1))
mu2_ind <- rnorm(N, e_mu2_case, sqrt(v_mu2))
r_ind <- runif(N, e_r_case - sqrt(12*v_r)/2, e_r_case + sqrt(12*v_r)/2)
bulk_case <- rep(NA, N)
for(i in 1:N){
  mu1 <- mu1_ind[i]
  mu2 <- mu2_ind[i]
  r <- r_ind[i]
  n_subset1 <- rbinom(1,prob=r, size=n_cells)
  n_subset2 <- n_cells - n_subset1
  sc_vec <- c(rnorm(n_subset1, mu1, sqrt(sc_v1_fix)), rnorm(n_subset2, mu2, sqrt(sc_v2_fix)))
  bulk_case[i] <- mean(sc_vec) + rnorm(1, 0, sqrt(v_err))
}

theo_bulk_mean_case <- e_mu1_case * e_r_case + (1 - e_r_case)*e_mu2_case
theo_bulk_var_case <- ((e_mu1_case - e_mu2_case)^2) * v_r + (e_r_case^2) * v_mu1 + 
((1 - e_r_case)^2) * v_mu2 + v_r * (v_mu1 + v_mu2) + v_err


#Draw the original single cell distributions
cex.axis <- 4
cex.lab <- 4
set.seed(1000)
sc_v1_fix
sc_v2_fix 
e_mu1
e_mu2 
e_r
e_mu1_case
e_mu2_case
e_r_case
x_grid <- seq(-2, 8, length=1000)
pdf_cont <- e_r * dnorm(x_grid, e_mu1, sqrt(sc_v1_fix)) + (1 - e_r) * dnorm(x_grid, e_mu2, sqrt(sc_v2_fix))
pdf_case <- e_r_case * dnorm(x_grid, e_mu1_case, sqrt(sc_v1_fix)) + (1 - e_r_case) * dnorm(x_grid, e_mu2_case, sqrt(sc_v2_fix))
png(paste0(out_path, "simu2_dist.png"), width=960, height=960)
#par(mar=c(15,15,10,15)) #bot, left, top, roght
#par(oma=c(10,8,2,8))
#par(mgp=c(15,6,0))
matplot(x_grid, cbind(pdf_cont, pdf_case), col=c("black", "red"), type="l", lty=1, xlab="", ylab="", cex.axis=cex.axis, cex.lab=cex.lab, lwd=4)
dev.off()

#Draw the original parameters (r)
r_grid <- seq(0, 1, length=1000)
r_pdf_cont <- dunif(r_grid, e_r - sqrt(12*v_r)/2, e_r + sqrt(12*v_r)/2)
r_pdf_case <- dunif(r_grid, e_r_case - sqrt(12*v_r)/2, e_r_case + sqrt(12*v_r)/2)
png(paste0(out_path, "simu2_r.png"), width=960, height=960)
#par(mar=c(10,10,2,2)) #bot, left, top, roght
#par(oma=c(10,10,2,2))
#par(mgp=c(8,8,0))
matplot(r_grid, cbind(r_pdf_cont, r_pdf_case), col=c("black", "red"), type="l", lty=1, xlab="", ylab="", cex.axis=cex.axis, cex.lab=cex.lab, lwd=4)
dev.off()

#Draw the original parameters (mu)
x_min <- min(e_mu1 - 3*sqrt(v_mu1), e_mu1_case - 3*sqrt(v_mu1), e_mu2 - 3*sqrt(v_mu2), e_mu2_case - 3*sqrt(v_mu2))
x_max <- max(e_mu1 + 3*sqrt(v_mu1), e_mu1_case + 3*sqrt(v_mu1),e_mu2 + 3*sqrt(v_mu2), e_mu2_case + 3*sqrt(v_mu2))
x_grid <- seq(x_min, x_max, length=1000)
mu1_pdf_cont <- dnorm(x_grid, e_mu1, sqrt(v_mu1))
mu1_pdf_case <- dnorm(x_grid, e_mu1_case, sqrt(v_mu1))
mu2_pdf_cont <- dnorm(x_grid, e_mu2, sqrt(v_mu2))
mu2_pdf_case <- dnorm(x_grid, e_mu2_case, sqrt(v_mu2))
png(paste0(out_path, "simu2_mu.png"), width=960, height=960)
#par(mar=c(5,5,10,10)) #bot, left, top, roght
#par(oma=c(4,4,2,2))
#par(mgp=c(10,4,0))
matplot(x_grid, cbind(mu1_pdf_cont, mu1_pdf_case, mu2_pdf_cont, mu2_pdf_case), type="l", col=c("black", "red", "black", "red"), lty=c(1, 1, 2, 2), xlab="", ylab="", cex.axis=cex.axis, cex.lab=cex.lab, lwd=4)
dev.off()

#Check
mean(bulk)
var(bulk)
print(theo_bulk_mean)
print(theo_bulk_var)
mean(bulk_case)
print(theo_bulk_mean_case)
var(bulk_case)
print(theo_bulk_var_case)
var.test(bulk, bulk_case)
t.test(bulk, bulk_case)
png(paste0(out_path, "simu2_boxplot.png"), width=960, height=960)
#par(mar=c(15,15,10,15)) #bot, left, top, roght
#par(oma=c(10,8,2,8))
#par(mgp=c(15,6,0))
boxplot(list(C=bulk, D=bulk_case), xlab="", ylab="", cex.axis=cex.axis, cex.lab=cex.lab)
dev.off()



