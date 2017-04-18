# TLBO with bad neighbourhood approach
# Learning experience with teaching learning based optimization or LETLBO Zou et al. (2015)
# The bad neighbourhood penalizes bad ranks and their region. "Learners" in these bad neighboruhood are attracted towards the global best.

TLBO_bn_fv_R=function(class_size,classes,dim,xmin,xmax,fitness,gen,printall=T,start_shuffle_prob=0.995){
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
    stop('xmax: not enough boundaries defined. Dimension and length(xmax) differ')
  }
  if(length(xmin)!=dim){
    stop('xmin: not enough boundaries defined. Dimension and length(xmin) differ')
  }
  
  if(class_size<20){
    readline("Your class size is low (equivilent to population). There's a high risk to get stuck in while loops with low classsize due to strategy seperation of the population. You've been warned. Press enter if you want to continue but better restart with higher class size.") 
  }
  # Initialise learners 
  school_old=matrix(t(runif(class_size*dim,xmin,xmax)),nrow=class_size,byrow=T) 
  xmax_mat=matrix(rep(xmax,class_size),nrow=class_size,byrow=T)
  xmin_mat=matrix(rep(xmin,class_size),nrow=class_size,byrow=T)
  # evaluate learners and define first teacher
  result_old=apply(school_old,1,function(x) fit_func(x))
  # evaluate teacher
  teacher=min(result_old)
  pos=which.min(result_old)
  teacher_loc=school_old[pos,]
  # calculate mean
  mean_loc=colMeans(school_old)
  # rank population
  ranks=rank(result_old,ties.method = "random")
  
  #create bad neighbourhood according to ranks
  # we find the centers of all neighbourhoods. 
  # They are equal to the position of learners that ranked badly
  centr=school_old[which(ranks>floor(class_size/2)),]
  # here the correct indexing is important to match the correct sphere size to rank
  bad_ranks=ranks[which(ranks>floor(class_size/2))]
  centr=centr[rev(sort(bad_ranks,index.return=T)$ix),]
  # we create sphere sizes. The worst gets the biggest sphere.
  sph=seq(from=1,to=0.3/(ceiling(class_size/2)),length=ceiling(class_size/2))
  hood=as.matrix(sph)%*%t(as.matrix(xmax))
  bad_hood=cbind(centr-hood,centr+hood)
  
  # optimisation algorithm update 
  k=1 # generation count
  if(printall){
    results=list() # to be returned in the end
  }else{
    results=matrix(ncol=dim,nrow=1)
  }
  # helper for indexing 
  class_vec=1:class_size
  
  # define new school matrix as indexes are accessed and different strategies are used
  school_new=matrix(ncol=dim,nrow=class_size)
  reshuffle_prob=1
  stuck=F
  stuck_ind=0
  while((teacher>error) && (k<=gen) && !stuck){
    # two strategies in teachers phase, determined through a and b
    a=runif(class_size)
    b=runif(class_size)
    # index of learners the first strategy will be applied to
    first=which(a<b)
    # index of learners the second strategy will be applied to
    second=which(a>=b)
    # while loop to make sure there are sufficient members in each strategy
    while(length(first)<3|length(second)<3){
      a=runif(class_size)
      b=runif(class_size)
      first=which(a<b)
      second=which(a>=b)
    }
    arun=runif(length(first))
    brun=runif(length(second))
    
    ############
    # strategy a
    ############
    
    # TF: teaching factor. It is either 1 or 2.
    # In Zou et al. 2015 it's described as TF=round(1+rand(0,1)[2-1])
    # Here, 2-1 was willingly ignored as it's equal to 1 and therefore there's no difference to the below
    # also in the same paper the random distribution is not defined further. I assume a uniform distribution.
    TF=round(1+runif(length(first)))
    
    # following two amtrices are important to avoid R repeating the mean_loc and teacher_loc in the wrong order
    # R tends to fill by coloumn not row by default
    
    # teacher loc matrix for vectorized computation
    teacher_loc_mat=matrix(rep(teacher_loc,length(first)),nrow=length(first),byrow=T)
    # mean loc matrix for vectorized computation
    mean_loc_mat=matrix(rep(mean_loc,length(first)),nrow=length(first),byrow=T)
    
    # calculate new school according to first strategy during teachers phase
    school_new[first,]=school_old[first,]+arun*(teacher_loc_mat-TF*mean_loc_mat)
    
    ############
    # strategy b
    ############
    # sample on other student that is part of the second group but not equal index
    # in a for(i in 1:whatever) loop this is is equal to find other != i
    other=sample(class_vec[second])
    while(any(other==class_vec[second])){
      other=sample(class_vec[second])
    }
    # which is better
    other_better=which(result_old[second]>result_old[other])
    # make id based on result
    secnd_id=second # second when other_better == F
    secnd_id[other_better]=other[other_better] # other when other_better == T
    # new teacher loc matrix of nrow of second strategy
    teacher_loc_mat=matrix(rep(teacher_loc,length(second)),nrow=length(second),byrow=T)
    # applying second strategy
    school_new[second,]=school_old[second,]+brun*(teacher_loc_mat-school_old[secnd_id,])
    
    #     ################# End of teacher's strategy
    # 
    # check if in new bad neighbourhood
    # this checks if all design variables (=length dimension) are in any of the previously defined bad neighbourhoods.
    # Here a for loop is used to go through every bad neighbourhood.
    # If a learner is in a bad neighbourhood, no matter which, the bool_vec is set to true for it's ID/index
    # A bad neighbourhood is defined by intervals (min,max) for each variable. For dimension 2 a bad neighbourd has two intervals for each dimension one.
    # The apply function checks if all variables of a learner are in the intervals.
    # R gives each TRUE the value 1 and each FALSE the value 0.
    # A learner is in a bad neighbourhood (or all variable values of a learner are in the intervals of a bad neighboruhood) when the sum is equal to the dimension(:= number of variables).
    # Checking each neighbourhood seperatly is necessary to not mix up intervals of different bad neighbourhoods.
    
    in_bad_hood=rep(FALSE,class_size)
    for(n in 1:(ceiling(length(ranks)/2))){
      bool_vec=apply(school_new,1,function(x) sum(x>=bad_hood[n,1:dim]&x<=bad_hood[n,(dim+1):(2*dim)])==dim)
      in_bad_hood[bool_vec]=TRUE
    }
    # length of learners in bad neighbourhoods or number of learners in bad neighbourhood
    ln_bn=length(in_bad_hood[in_bad_hood])
    if(ln_bn>0){
      rndm_gb=runif(ln_bn)
      teacher_loc_mat=matrix(rep(teacher_loc,ln_bn),nrow=ln_bn,byrow=T)
    # learners in bad neighbourhood go towards the teacher
      school_new[in_bad_hood,]=school_new[in_bad_hood,]-rndm_gb*(teacher_loc_mat-school_new[in_bad_hood,])
    }
    bound_max=which(school_new> xmax_mat)
    bound_min=which(school_new < xmin_mat)
    
    school_new[bound_max] = xmax_mat[bound_max]-(school_new[bound_max]-xmax_mat[bound_max])
    school_new[bound_min] = xmin_mat[bound_min]-(school_new[bound_min]-xmin_mat[bound_min])
    
    # evaluation of objective function
    result_new=apply(school_new,1,function(x) fit_func(x))
    
    # update results, school if better
    newbetter=which(result_new<=result_old)
    result_old[newbetter]=result_new[newbetter]
    school_old[newbetter,]=school_new[newbetter,]
    # rank learners
    ranks=rank(result_old,ties.method = 'random')
    # calculate new bad neighbourhood based on updated result_old
    centr=school_old[which(ranks>floor(class_size/2)),]
    bad_ranks=ranks[which(ranks>floor(class_size/2))]
    centr=centr[rev(sort(bad_ranks,index.return=T)$ix),]
    sph=seq(from=0.1+(1-0.1)*(gen-k)/gen,to=0.3/(ceiling(class_size/2))*(gen-k)/gen,length=ceiling(class_size/2))
    hood=as.matrix(sph)%*%t(as.matrix(xmax))
    bad_hood=cbind(centr-hood,centr+hood)
    
    #calculating new teacher as bad neighbourhood approach is also applied in learners' phase
    teacher=min(result_old)
    pos=which.min(result_old)
    teacher_loc=school_old[pos,]
    
    ############################
    # beginning of learners' phase 
    #############################
    
    # two strategies in learners' phase, determined through a and b
    a=runif(class_size)
    b=runif(class_size)
    first=which(a<b)
    second=which(a>=b)
    while(length(first)<4|length(second)<4){
      a=runif(class_size)
      b=runif(class_size)
      first=which(a<b)
      second=which(a>=b)
    }
    # first strategy, comparison to one other student
    other=sample(class_vec[first])
    while(any(other==class_vec[first])){
      other=sample(class_vec[first])
    }
    # which is better
    other_better=which(result_old[first]>result_old[other])
    # make ids based on result for substraction 
    # when other is better it's school_old[other,]-school_old[first,] and vice versa
    other_id=other
    other_id[other_better]=first[other_better]
    #
    frst_id=first
    frst_id[other_better]=other[other_better]
    school_new[first,]=school_old[first,]+runif(length(first))*(school_old[frst_id,]-school_old[other_id,])
    # 2nd strategy: comparing two other students
    other_one=sample(class_vec[second])
    while(any(other_one==class_vec[second])){
      other_one=sample(class_vec[second])
    }
    other_two=sample(class_vec[second])
    while(any(other_two==class_vec[second]|other_two==other_one)){
      other_two=sample(class_vec[second])
    }
    
    id_other=which(result_old[other_one]<result_old[other_two])
    # make ids based on result for substraction 
    # when other is better it's school_old[other,]-school_old[first,] and vice versa
    one_better=other_two
    one_better[id_other]=other_one[id_other]
    two_better=other_one
    two_better[id_other]=other_two[id_other]
    
    # applying second strategy
    school_new[second,]=school_old[second,]+runif(length(second))*(school_old[one_better,]-school_old[two_better,])
    #bad neighbourhood
    in_bad_hood=rep(FALSE,class_size)
    for(n in 1:(ceiling(length(ranks)/2))){
      bool_vec=apply(school_new,1,function(x) sum(x>=bad_hood[n,1:dim]&x<=bad_hood[n,(dim+1):(2*dim)])==dim)
      in_bad_hood[bool_vec]=TRUE
    }
    ln_bn=length(in_bad_hood[in_bad_hood])
    if(ln_bn>0){
      rndm_gb=runif(ln_bn)
      teacher_loc_mat=matrix(rep(teacher_loc,ln_bn),nrow=ln_bn,byrow=T)
      # bad neighbourhood update
      school_new[in_bad_hood,]=school_new[in_bad_hood,]-rndm_gb*(teacher_loc_mat-school_new[in_bad_hood,])
    }
    bound_max=which(school_new> xmax_mat)
    bound_min=which(school_new < xmin_mat)
    
    school_new[bound_max] = xmax_mat[bound_max]-(school_new[bound_max]-xmax_mat[bound_max])
    school_new[bound_min] = xmin_mat[bound_min]-(school_new[bound_min]-xmin_mat[bound_min])
    
    # parallel evaluation
    result_new=apply(school_new,1,function(x) fit_func(x))
    # replacing school with better solution
    newbetter=which(result_new<=result_old)
    result_old[newbetter]=result_new[newbetter]
    school_old[newbetter,]=school_new[newbetter,]
    
    ranks=rank(result_old,ties.method = 'random')
    
    # new bad neighbourhood
    centr=school_old[which(ranks>floor(class_size/2)),]
    bad_ranks=ranks[which(ranks>floor(class_size/2))]
    centr=centr[rev(sort(bad_ranks,index.return=T)$ix),]
    hood=as.matrix(sph)%*%t(as.matrix(xmax))
    bad_hood=cbind(centr-hood,centr+hood)
    #calculating new teacher
    if((teacher-min(result_old))<1e-8){
      stuck_ind=stuck_ind+1
    }else{
      stuck_ind=0
    }
    if(stuck_ind>=100){
      stuck=T
    }
    teacher=min(result_old)
    pos=which.min(result_old)
    teacher_loc=school_old[pos,]
    
    #new mean
    mean_loc=colMeans(school_old)
    if(printall){
      results[[k]]=cbind(school_old,result_old)
    }else{
      if(k==1){
        results=cbind(t(teacher_loc),teacher)
      }
      else{
        resulttemp=cbind(t(teacher_loc),teacher)
        results=rbind(results,resulttemp)
      }
    }
    reshuffle_prob=reshuffle_prob*start_shuffle_prob
    a=runif(1)
    # reinitilize worst half if prob reached
    if(a>reshuffle_prob){
      worst=length(ranks[which(ranks>floor(class_size/2))])
      school_old[which(ranks>floor(class_size/2)),]=matrix(t(runif(worst*dim,xmin,xmax)),nrow=worst,byrow=T) 
      reshuffle_prob=start_shuffle_prob
    }
    
    k=k+1
  }
  if(printall){
    results[[k]]=as.vector(cbind(t(teacher_loc),teacher))
  }else{
    resulttemp=cbind(t(teacher_loc),teacher)
    results=rbind(results,resulttemp)
  }
  return(results)
}
