# run pso with shuffled complexes and bad neighborhood
# requires cec2013 package
# run cec2013 benchmark

if (!require(cec2013,character.only = TRUE)){
    stop("cec2013 package not installed")
}else{
  require(cec2013)
}
source("pso_full_vec_sce_bn.R")
rndmint <- read.table("rndmint.dat", quote="\"", comment.char="")
gen=10000*10/100
for(i in 1:30){
  seed=rndmint[100-i,1]
  set.seed(seed)
  PSO=PSO_sce_bn_fv_R(100,2,10,rep(-100,10),rep(100,10),"1",gen=gen,F,1e-6)
  write.csv(PSO,paste('bfnc_1_rep',i,'PSO_sce_bn.csv',sep=''))
}
for(i in 1:30){
  seed=rndmint[100-i,1]
  set.seed(seed)
  PSO=PSO_sce_bn_fv_R(100,2,10,rep(-100,10),rep(100,10),"2",gen=gen,F,0.999)
  write.csv(PSO,paste('bfnc_2_rep',i,'PSO_sce_bn.csv',sep=''))
}
for(i in 1:30){
  seed=rndmint[100-i,1]
  set.seed(seed)
  PSO=PSO_sce_bn_fv_R(100,2,10,rep(-100,10),rep(100,10),"6",gen=gen,F,1e-6)
  write.csv(PSO,paste('bfnc_6_rep',i,'PSO_sce_bn.csv',sep=''))
}
for(i in 1:30){
  seed=rndmint[100-i,1]
  set.seed(seed)
  PSO=PSO_sce_bn_fv_R(100,2,10,rep(-100,10),rep(100,10),"8",gen=gen,F,0.8)
  write.csv(PSO,paste('bfnc_8_rep',i,'PSO_sce_bn.csv',sep=''))
}
for(i in 1:30){
  seed=rndmint[100-i,1]
  set.seed(seed)
  PSO=PSO_sce_bn_fv_R(100,2,10,rep(-100,10),rep(100,10),"9",gen=gen,F,1e-6)#
  write.csv(PSO,paste('bfnc_9_rep',i,'PSO_sce_bn.csv',sep=''))
}
for(i in 1:30){
  seed=rndmint[100-i,1]
  set.seed(seed)
  PSO=PSO_sce_bn_fv_R(100,2,10,rep(-100,10),rep(100,10),"10",gen=gen,F,1e-6)#
  write.csv(PSO,paste('bfnc_10_rep',i,'PSO_sce_bn.csv',sep=''))
}
for(i in 1:30){
  seed=rndmint[100-i,1]
  set.seed(seed)
  PSO=PSO_sce_bn_fv_R(100,2,10,rep(-100,10),rep(100,10),"11",gen=gen,F,1e-6)#
  write.csv(PSO,paste('bfnc_11_rep',i,'PSO_sce_bn.csv',sep=''))
}
for(i in 1:30){
  seed=rndmint[100-i,1]
  set.seed(seed)
  PSO=PSO_sce_bn_fv_R(100,2,10,rep(-100,10),rep(100,10),"12",gen=gen,F,1e-6)#
  write.csv(PSO,paste('bfnc_12_rep',i,'PSO_sce_bn.csv',sep=''))
}
for(i in 1:30){
  seed=rndmint[100-i,1]
  set.seed(seed)
  PSO=PSO_sce_bn_fv_R(100,2,10,rep(-100,10),rep(100,10),"14",gen=gen,F,1e-6)#
  write.csv(PSO,paste('bfnc_14_rep',i,'PSO_sce_bn.csv',sep=''))
}
for(i in 1:30){
  seed=rndmint[100-i,1]
  set.seed(seed)
  PSO=PSO_sce_bn_fv_R(100,2,10,rep(-100,10),rep(100,10),"28",gen=gen,F,1e-6)#
  write.csv(PSO,paste('bfnc_28_rep',i,'PSO_sce_bn.csv',sep=''))
}