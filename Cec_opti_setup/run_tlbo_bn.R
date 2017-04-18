# run tlbo with bad neighbourhood approach
# run cec2013 benchmark
source("tlbo_full_vec_bn.R")
if (!require(cec2013,character.only = TRUE)){
  stop("cec2013 package not installed")
}else{
  require(cec2013)
}
rndmint <- read.table("rndmint.dat", quote="\"", comment.char="")
gen=10000*10/100/2
for(i in 1:30){
  seed=rndmint[100-i,1]
  set.seed(seed)
  TLBO=TLBO_bn_fv_R(100,1,10,rep(-100,10),rep(100,10),"1",gen=gen,F,1e-6)
  write.csv(TLBO,paste('bfnc_1_rep',i,'tlbo_bn.csv',sep=''))
}
for(i in 1:30){
  seed=rndmint[100-i,1]
  set.seed(seed)
  TLBO=TLBO_bn_fv_R(100,1,10,rep(-100,10),rep(100,10),"2",gen=gen,F,0)
  write.csv(TLBO,paste('bfnc_2_rep',i,'tlbo_bn.csv',sep=''))
}
for(i in 1:30){
  seed=rndmint[100-i,1]
  set.seed(seed)
  TLBO=TLBO_bn_fv_R(100,1,10,rep(-100,10),rep(100,10),"6",gen=gen,F,0.5)
  write.csv(TLBO,paste('bfnc_6_rep',i,'tlbo_bn.csv',sep=''))
}
for(i in 1:30){
  seed=rndmint[100-i,1]
  set.seed(seed)
  TLBO=TLBO_bn_fv_R(100,1,10,rep(-100,10),rep(100,10),"8",gen=gen,F,10e13)
  write.csv(TLBO,paste('bfnc_8_rep',i,'tlbo_bn.csv',sep=''))
}
for(i in 1:30){
  seed=rndmint[100-i,1]
  set.seed(seed)
  TLBO=TLBO_bn_fv_R(100,1,10,rep(-100,10),rep(100,10),"9",gen=gen,F,1e-6)#
  write.csv(TLBO,paste('bfnc_9_rep',i,'tlbo_bn.csv',sep=''))
}
for(i in 1:30){
  seed=rndmint[100-i,1]
  set.seed(seed)
  TLBO=TLBO_bn_fv_R(100,1,10,rep(-100,10),rep(100,10),"10",gen=gen,F,1e-6)#
  write.csv(TLBO,paste('bfnc_10_rep',i,'tlbo_bn.csv',sep=''))
}
for(i in 1:30){
  seed=rndmint[100-i,1]
  set.seed(seed)
  TLBO=TLBO_bn_fv_R(100,1,10,rep(-100,10),rep(100,10),"11",gen=gen,F,1e-6)#
  write.csv(TLBO,paste('bfnc_11_rep',i,'tlbo_bn.csv',sep=''))
}
for(i in 1:30){
  seed=rndmint[100-i,1]
  set.seed(seed)
  TLBO=TLBO_bn_fv_R(100,1,10,rep(-100,10),rep(100,10),"12",gen=gen,F,0)#
  write.csv(TLBO,paste('bfnc_12_rep',i,'tlbo_bn.csv',sep=''))
}
for(i in 1:30){
  seed=rndmint[100-i,1]
  set.seed(seed)
  TLBO=TLBO_bn_fv_R(100,1,10,rep(-100,10),rep(100,10),"14",gen=gen,F,0)
  write.csv(TLBO,paste('bfnc_14_rep',i,'tlbo_bn.csv',sep=''))
}
for(i in 1:30){
  seed=rndmint[100-i,1]
  set.seed(seed)
  TLBO=TLBO_bn_fv_R(100,1,10,rep(-100,10),rep(100,10),"28",gen=gen,F,0.2)
  write.csv(TLBO,paste('bfnc_28_rep',i,'tlbo_bn.csv',sep=''))
}
