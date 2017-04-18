# pso no shuffling with bad neighbourhood approach
PSO_fv_bn_R=function(pop,complexes,dim,xmin,xmax,fitness,gen,printall=T,start_shuffle_prob=0.995){
  library(cec2013)
  fit_func<-function(z){
    switch(fitness,
           "1"=return(cec2013(1,z)),
           "2"=return(cec2013(2,z)),
           "3"=return(cec2013(3,z)),
           "4"=return(cec2013(4,z)),
           "5"=return(cec2013(5,z)),
           "6"=return(cec2013(6,z)),
           "7"=return(cec2013(7,z)),
           "8"=return(cec2013(8,z)),
           "9"=return(cec2013(9,z)),
           "10"=return(cec2013(10,z)),
           "11"=return(cec2013(11,z)),
           "12"=return(cec2013(12,z)),
           "13"=return(cec2013(13,z)),
           "14"=return(cec2013(14,z)),
           "15"=return(cec2013(15,z)),
           "16"=return(cec2013(16,z)),
           "17"=return(cec2013(17,z)),
           "18"=return(cec2013(18,z)),
           "19"=return(cec2013(19,z)),
           "20"=return(cec2013(20,z)),
           "21"=return(cec2013(21,z)),
           "22"=return(cec2013(22,z)),
           "23"=return(cec2013(23,z)),
           "24"=return(cec2013(24,z)),
           "25"=return(cec2013(25,z)),
           "26"=return(cec2013(26,z)),
           "27"=return(cec2013(27,z)),
           "28"=return(cec2013(28,z))
    )
  }
  # accepable errors according to cec 2013
  error<- switch(fitness,
                 "1"=-1400+1e-8,
                 "2"=-1300+1e-8, 
                 "3"=-1200+1e-8, 
                 "4"=-1100+1e-8, 
                 "5"=-1000+1e-8, 
                 "6"=-900+1e-8, 
                 "7"=-800+1e-8, 
                 "8"=-700+1e-8, 
                 "9"=-600+1e-8, 
                 "10"=-500+1e-8, 
                 "11"=-400+1e-8, 
                 "12"=-300+1e-8, 
                 "13"=-200+1e-8, 
                 "14"=-100+1e-8, 
                 "15"=100+1e-8, 
                 "16"=200+1e-8, 
                 "17"=300+1e-8, 
                 "18"=400+1e-8, 
                 "19"=500+1e-8, 
                 "20"=600+1e-8, 
                 "21"=700+1e-8, 
                 "22"=800+1e-8, 
                 "23"=900+1e-8, 
                 "24"=1000+1e-8, 
                 "25"=1100+1e-8, 
                 "26"=1200+1e-8, 
                 "27"=1300+1e-8, 
                 "28"=1400+1e-8 
  )
  
  if(length(xmax)!=dim){
    stop('xmax: wrong number of boundaries defined. Dimension and length(xmax) differ')
  }
  if(length(xmin)!=dim){
    stop('xmin: wrong number of boundaries defined. Dimension and length(xmin) differ')
  }
  
  ###################
  #### Initialisation
  ###################
  
  vmax=0.3*(xmax-xmin)
  swarm=matrix(t(runif(pop*dim,xmin,xmax)),nrow=pop,byrow=T)
  vel=matrix(rep(0,dim*pop),nrow=pop) 
  pbest_loc=swarm
  xmax_mat=matrix(rep(xmax,pop),nrow=pop,byrow=T)
  xmin_mat=matrix(rep(xmin,pop),nrow=pop,byrow=T)
  
  # parallel evaluation
  result=apply(swarm,1,function(x) fit_func(x))
  pbest=result
  gbest=min(result)
  pos=which.min(result)
  gbest_loc=swarm[pos,]
  
  ranks=rank(result,ties.method = "random")
  
  bad_hood=matrix(nrow=ceiling(length(ranks)/2),ncol=dim*2)
  centr=swarm[which(ranks>floor(pop/2)),]
  sph=seq(from=1,to=0.3/(ceiling(pop/2)),length=ceiling(pop/2))
  bad_ranks=ranks[which(ranks>floor(pop/2))]
  centr=centr[rev(sort(bad_ranks,index.return=T)$ix),]
  hood=as.matrix(sph)%*%t(as.matrix(vmax))
  bad_hood=cbind(centr-hood,centr+hood)
  ##############
  # Update
  ##############
  k=1 # generation index
  if(printall){
    results=list() # to be returned in the end
  }else{
    results=matrix(ncol=dim,nrow=1)
  }
  reshuffle_prob=1
  fac=(gen-k)/gen
  fac2=k/gen
  samegbest=0
  stuck=F
  stuck_ind=0
  while(k<=gen && gbest>=error && !stuck){
    gbest_loc_mat=matrix(rep(gbest_loc,pop),nrow=pop,byrow=T)
    wmax=(0.9-0.2)*fac+0.2 # based on Suganthan, Roshida and yoshida et al. #0.9
    #wmin=0.2
    vmax=(0.5*(xmax-xmin)-(xmax-xmin)/20)*fac+(xmax-xmin)/20
    #particle specific vmax, low when gbest (improve exploitation) and high when ranked badly (go on exploration)
    vmax_part=matrix(rep(vmax,pop),nrow=pop,byrow=T)/(pop-ranks+1)
    w=wmax#=(wmax-wmin)*ranks[i]/pop+wmin
    c2=0.5+(2.5-0.5)*fac2 ## increasing attraction to global best
    c1=0.1+(1.5-0.1)*fac # decreasing attraction to personal best
    vel=w*vel+c1*runif(pop)*(pbest_loc-swarm)+c2*runif(pop)*(gbest_loc_mat-swarm)
    too_fast=which(vel>vmax_part)
    vel[too_fast]=vmax_part[too_fast]
    too_fast=which(vel<(-vmax_part))
    vel[too_fast]=-vmax_part[too_fast]
    swarm=swarm+vel
    # reflection back into the space when 'hitting' the boundary
    bound_max=which(swarm > xmax_mat)
    bound_min=which(swarm < xmin_mat)
    swarm[bound_max] = xmax_mat[bound_max]-(swarm[bound_max]-xmax_mat[bound_max])
    swarm[bound_min] = xmin_mat[bound_min]-(swarm[bound_min]-xmin_mat[bound_min])
    #
    
    in_bad_hood=rep(FALSE,pop)
    for(n in 1:(round(length(ranks)/2))){
      bool_vec=apply(swarm,1,function(x) sum(x>=bad_hood[n,1:dim]&x<=bad_hood[n,(dim+1):(2*dim)])==dim)
      in_bad_hood[bool_vec]=TRUE
    }
    ln_bn=length(in_bad_hood[in_bad_hood])
    if(ln_bn>0){
      gbest_loc_mat=matrix(rep(gbest_loc,ln_bn),nrow=ln_bn,byrow=T)
     swarm[in_bad_hood,]=swarm[in_bad_hood,]+c2*runif(ln_bn)*(gbest_loc_mat-swarm[in_bad_hood,])
    }
    # parallel evaluation
    result=apply(swarm,1,function(x) fit_func(x))
    #new pbest
    newpbest=which(result<pbest)
    pbest[newpbest]=result[newpbest]
    pbest_loc[newpbest,]=swarm[newpbest,]
    if(gbest>min(result)){
      if((gbest-min(result))<1e-8){
        samegbest=samegbest+1
        stuck_ind=stuck_ind+1
      }else{
        stuck_ind=0
      }
      dif=gbest-min(result)
      gbest=min(result)
      pos=which.min(result)
      gbest_loc=swarm[pos,]
      fac=(gen-k)/gen
      fac2=k/gen
      samegbest=0
    }else{
      # print('here')
      samegbest=samegbest+1
      stuck_ind=stuck_ind+1
    }
    if(samegbest>=10){
      if(fac>=1){
        fac=1
      }else{
        fac=fac*1.1
      }
      fac2=fac2*0.9

    }
    if(stuck_ind>=100){
      stuck=T
    }
    # print('fac')
    # print(fac)
    # print('fac2')
    # print(fac2)
    # print('---------')
    ranks=rank(result,ties.method = "random")
    bad_hood=matrix(nrow=ceiling(length(ranks)/2),ncol=dim*2)
    centr=swarm[which(ranks>floor(pop/2)),]
    sph=seq(from=1,to=0.2/(ceiling(pop/2)),length=ceiling(pop/2))
    bad_ranks=ranks[which(ranks>floor(pop/2))]
    centr=centr[rev(sort(bad_ranks,index.return=T)$ix),]
    hood=as.matrix(sph)%*%t(as.matrix(vmax))
    bad_hood=cbind(centr-hood,centr+hood)
    if(printall){
      results[[k]]=cbind(swarm,result)
    }else{
      if(k==1){
        results=cbind(t(gbest_loc),gbest)
      }
      else{
        resulttemp=cbind(t(gbest_loc),gbest)
        results=rbind(results,resulttemp)
      }
    }
    k=k+1 
    reshuffle_prob=reshuffle_prob*start_shuffle_prob
    a=runif(1)
    # reinitilize worst half if prob reached
    if(a>reshuffle_prob){
      worst=length(ranks[which(ranks>floor(pop/2))])
      swarm[which(ranks>floor(pop/2)),]=matrix(t(runif(worst*dim,xmin,xmax)),nrow=worst,byrow=T) 
      reshuffle_prob=start_shuffle_prob
    }
  }
  if(printall){
    results[[k]]=cbind(t(gbest_loc),gbest)
  }else{
    resulttemp=cbind(t(gbest_loc),gbest)
    results=rbind(results,resulttemp)
  }
  return(results)
}
