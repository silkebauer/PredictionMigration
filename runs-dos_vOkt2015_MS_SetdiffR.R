
# the model for investigating how predictability along a migratory route
# affects migratory behaviour
#
# codes for running the model

# setup {{{
source("c:/A_DropBox/Dropbox/PredictMig/Runs4Paper/pred_mig-dos_vJune2015.R")

library(gplots)
setwd("c:/A_DropBox/Dropbox/PredictMig/Runs4Paper/VaryRho_overSites/")

#}}}


# full set of scenarios
(try_Rfuture <- c(1, 5))    
(try_sigma <- c(seq(0.1, 1.0, by=0.1)))      
(try_rho <- c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
(try_K <- c(5)) # total number of sites; 

rho.k1 <- c(seq(0.0, 0.9, by=0.1),0.99)
rho.kset <- rep(0.99, 11)

rho.k <- matrix(rbind( cbind(rho.kset, rho.k1, rho.k1, rho.k1, rho.k1),
                       cbind(rho.k1, rho.kset, rho.k1, rho.k1, rho.k1),
                       cbind(rho.k1, rho.k1, rho.kset, rho.k1, rho.k1),
                       cbind(rho.k1, rho.k1, rho.k1, rho.kset, rho.k1),
                       cbind(rho.k1, rho.k1, rho.k1, rho.k1, rho.kset) ), ncol=5)

# for test runs
# Rfuture <- 1
# sigma <- 0.5
# rho <- 0.9

count <- 0

for (kkk in 1:length(try_K)) {
for (rrr in 1:length(try_Rfuture)) {
  for(mmm in 1:length(try_sigma)) {
#    for(xxx in 1:length(try_rho)) { 
  for (xxx in 1:nrow(rho.k))  {
    
    count <- count+1
    file.name <- paste("run_", count, ".pdf", sep="")
            

P <- init.param()
P$K <- try_K[kkk]
rho <- try_rho[xxx]
scale.rho <- (rho ^ (1/P$K))
P <- init.param(P, Q = 100,
                mort.max = c(0.0001, rep(0.005, P$K)),          
                mort.min = c(0.0001, rep(0.0001, P$K)),
                alpha   = c(rep(0.05, P$K+1)), 
                beta    = rep(0, P$K+1),           
                tau     = c(seq(0, 70, length.out=P$K+1)), 
                rho     = c(NA, rho.k[xxx,]),            
                sigma   = c(rep(try_sigma[mmm], P$K+1)), 
                R.wide  = c(30),       
                R.future = c(try_Rfuture[rrr])
                )

res <- backward(P)
L <- P$De*apply(res$Action, c(1,2), function(z) min(which( z < 1)) - 1) 
P$tau*P$De
(policy <- L[P$Q+1+seq(-P$Q*P$De, P$Q*P$De, 1)/P$De,])


# calculate forward as an exceptation for all `i`
    res$S <- init.Action(P, 0)
    i <- -P$Q:P$Q
    yy <- dnorm(i, 0, P$sigma[1])
    yy <- yy/sum(yy)
    res$S[,1,1] <- yy
    res <- forward(res, P)
    z0 <- apply(res$S, c(2,3), sum) 
    Surv0 <- apply(res$S, 3, sum)
    S0 <- unlist(as.numeric(Surv0[length(Surv0)]))


# calculate the average reproductive value
    s <- res$S[,1,1]/sum(res$S[,1,1])
    (ReproValue <- sum(s * res$V[,1,1], na.rm=TRUE))



res$S <- init.Action(P, 0)     # run forward for advanced spring
res$S[1,1,1] <- 1
res <- forward(res, P)
z1 <- apply(res$S, c(2,3), sum) 
Surv1 <- apply(res$S, 3, sum)
Surv1[length(Surv1)]
S1 <- unlist(as.numeric(Surv1[length(Surv1)]))

res$S <- init.Action(P, 0)    # run forward for unchanged springs
res$S[101,1,1] <- 1
res <- forward(res, P)
z2 <- apply(res$S, c(2,3), sum) 
Surv2 <- apply(res$S, 3, sum)
Surv2[length(Surv2)]
S2 <- unlist(as.numeric(Surv2[length(Surv2)]))

res$S <- init.Action(P, 0)   # run forward for delayed spring
res$S[201,1,1] <- 1
res <- forward(res, P)
z3 <- apply(res$S, c(2,3), sum) 
Surv3 <- apply(res$S, 3, sum)
Surv3[length(Surv3)]
S3 <- unlist(as.numeric(Surv3[length(Surv3)]))



#------------------- text output to files -----------------------------            

    # Write parameter values in file
    write.table(cbind(P$K, P$R.wide*P$De, P$R.future, P$alpha[1], 
                      P$mort.max[1], P$mort.max[2], P$mort.max[3], 
                      P$sigma[1], P$sigma[2], P$sigma[3],
                      round(P$rho[1],5), round(P$rho[2],5), round(P$rho[3],5), round(P$rho[4],5), 
                      round(P$rho[5],5), round(P$rho[6],5)), 
                append=T, sep=",", file="Parameters.txt", col.names=F)

    # output from backward computation - optimal departures
    #write(paste("Run ", count), append=T, sep=",", file="dep.txt")
    write.table(cbind(count, policy), append=T, sep=",", file="dep.txt", col.names=F)

    
    # output from forward computation - proportion of population at sites
    z0_temp <- cbind(t(z0[,seq(from=1, to=1001, by=10)]), NA)
    z1_temp <- cbind(t(z1[,seq(from=1, to=1001, by=10)]), NA)
    z2_temp <- cbind(t(z2[,seq(from=1, to=1001, by=10)]), NA)
    z3_temp <- cbind(t(z3[,seq(from=1, to=1001, by=10)]), NA)
    

#     write.table(cbind(count, z0_temp[,c(1:11)], z0_temp[,2], z0_temp[,3],  
#                              z1_temp[,1], z1_temp[,2], z1_temp[,3],
#                              z2_temp[,1], z2_temp[,2], z2_temp[,3],
#                              z3_temp[,1], z3_temp[,2], z3_temp[,3]),
#                 append=T, sep=",",file="fw_test_2.txt", col.names=F)

    t0.matrix <- t1.matrix <- t2.matrix <- t3.matrix <- matrix(, nrow = 101, ncol = 13)
    t0.matrix[1:nrow(z0_temp),1:ncol(z0_temp)] <- z0_temp
    t1.matrix[1:nrow(z1_temp),1:ncol(z1_temp)] <- z1_temp
    t2.matrix[1:nrow(z2_temp),1:ncol(z2_temp)] <- z2_temp
    t3.matrix[1:nrow(z3_temp),1:ncol(z3_temp)] <- z3_temp


    write.table(cbind(count, t0.matrix[,1:12], t1.matrix[,1:12], t2.matrix[,c(1:12)], t3.matrix[,c(1:12)]),
                append=T, sep=",",file="fw_test_2.txt", col.names=F)


    # RV and survival values 
    write.table(cbind(count, round(ReproValue,5), round(S0,5), round(S1,5), round(S2,5), round(S3,5)), 
                append=T, sep=",", file="RV_Surv.txt", col.names=F)



#------------------- graphical output to pdf -----------------------------            


pdf(file.name, paper="a4", width=10, height=10)
layout(matrix(c(1,2,3,4,5,6,7,8),ncol=2,byrow=T),height=c(2.5,3,3,3,4.5))


# Plot 1 & 2 - textplot with parameter values 
    par(mar=c(2,2,2,2))
    textplot(data.frame(tau=round(P$tau,1),rho=round(P$rho,3),sigma=round(P$sigma,1),R.wide=P$R.wide,
                        R.future=P$R.future), 
             cex=0.9, mar=c(1,0.2,2,1), show.rownames=F)
    mtext("Parameters: ",3, cex=0.9)


# Plot 2  - Optimal departure from sites
    par(mar=c(4,3,3,3),cex.lab=0.8,cex.axis=0.8)
        #P$tau*P$De
        #L[P$Q+1+seq(-P$Q*P$De, P$Q*P$De, 1)/P$De,]
    matplot(P$De*seq(-P$Q, P$Q, 1), L[,-(P$K+1)], type="l", lty=0:P$K+1, xlim=c(-10, 11.5),
            col=0:P$K+1, xlab="advancement of spring (real time)", ylab="leaving time (real time)")
    legend("topright", legend=0:(P$K-1), lty=1:(P$K), col=1:(P$K), bty="n", title="location", cex=0.55)


# Plot 3 - Textplot with leaving strategy
    textplot(L[P$Q+1+seq(-P$Q*P$De, P$Q*P$De, 1)/P$De,],cex=1,mar=c(1,1,1,1))
    #mtext("Departure strategy:",3,cex=1)


# Plot 4 - result of forward                       
      image(P$De*(0:P$N), 0:P$K, t(1-z1), yaxt="n", xlab="time", ylab="location", main="i = -10")
      axis(side=2, at=0:P$K)
      abline(h=seq(-0.5,P$K+0.5, 1), lwd=2)


# Plot 5 - result of forward
      image(P$De*(0:P$N), 0:P$K, t(1-z2), yaxt="n", xlab="time", ylab="location", main="i = 0")
      axis(side=2, at=0:P$K)
      abline(h=seq(-0.5,P$K+0.5, 1), lwd=2)


# plots 6-8 - reproductive values 
      matplot(P$De*(0:P$N), t(res$V[1,,]), type="l", xlab="time", ylab="reproductive value", main="i = -10")
      legend("bottomleft", legend=0:(P$K), lty=(0:(P$K))+1, col=(0:(P$K))+1, bty="n", title="location")
      
      matplot(P$De*(0:P$N), t(res$V[101,,]), type="l", xlab="time", ylab="reproductive value", main="i = 0")
      legend("bottomleft", legend=0:(P$K), lty=(0:(P$K))+1, col=(0:(P$K))+1, bty="n", title="location")
      
      matplot(P$De*(0:P$N), t(res$V[201,,]), type="l", xlab="time", ylab="reproductive value", main="i = 10")
      legend("bottomleft", legend=0:(P$K), lty=(0:(P$K))+1, col=(0:(P$K))+1, bty="n", title="location")

dev.off()    

rm(z,z0,res,Action,A,L,H.L,H.S,P,V)

}}}}





# screen output - not used for simulations

# layout(matrix(1:3, ncol=1))
# image(P$De*(0:P$N), 0:P$K, t(1-z1), yaxt="n", xlab="time", ylab="location",
#       main="i = -10")
# axis(side=2, at=0:P$K)
# #abline(h=seq(-0.5,P$K+0.5, 1), lwd=2)
# image(P$De*(0:P$N), 0:P$K, t(1-z2), yaxt="n", xlab="time", ylab="location",
#       main="i = 0")
# axis(side=2, at=0:P$K)
# #abline(h=seq(-0.5,P$K+0.5, 1), lwd=2)
# image(P$De*(0:P$N), 0:P$K, t(1-z3), yaxt="n", xlab="time", ylab="location",
#       main="i = 10")
# axis(side=2, at=0:P$K)
# #abline(h=seq(-0.5,P$K+0.5, 1), lwd=2)
# layout(1)
# #}}}
# 
# i <- 0:P$N
# i <- 590:750
# layout(matrix(1:3, ncol=1))
# image(P$De*(i), 0:P$K, t(1-z1[,i+1]), yaxt="n", xlab="time", ylab="location",
#       main="i = -10")
# axis(side=2, at=0:P$K)
# #abline(h=seq(-0.5,P$K+0.5, 1), lwd=2)
# image(P$De*(i), 0:P$K, t(1-z2[,i+1]), yaxt="n", xlab="time", ylab="location",
#       main="i = 0")
# axis(side=2, at=0:P$K)
# #abline(h=seq(-0.5,P$K+0.5, 1), lwd=2)
# image(P$De*(i), 0:P$K, t(1-z3[,i+1]), yaxt="n", xlab="time", ylab="location",
#       main="i = 10")
# axis(side=2, at=0:P$K)
# #abline(h=seq(-0.5,P$K+0.5, 1), lwd=2)
# layout(1)



# vim: foldmethod=marker
