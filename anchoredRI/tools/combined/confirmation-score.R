source("generator.R")
d <- read.csv("confirmation-limits.csv")
for (k in c("lower","upper")) {
 d[[paste0("z_",k)]] <- NA_real_
 for (ii in split(seq_len(nrow(d)),d$family)) {
  i <- ii[1]
  d[[paste0("z_",k)]][ii] <- qnorm(if(k=="lower").025 else .975) - stage2_healthy_z(d[[k]][ii],d$lambda[i],d$mu[i],d$sigma[i])
 }
}
write.csv(d,"confirmation-scored.csv",row.names=FALSE)
