# TLBO no shuffling, no bad neighbourhood approach
# Learning experience with teaching learning based optimization or LETLBO Zou et al. (2015)

TLBO_fv_R=function(class_size,classes,dim,xmin,xmax,fitness,gen,printall=T,start_shuffle_prob=0.995){
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
  teacher=c()
  means=c()
  result_old=c()
  teacher_loc=matrix(ncol=dim,nrow=classes)
  mean_loc=matrix(ncol=dim,nrow=classes)
  # evaluate learners and define first teacher
  result_old=apply(school_old,1,function(x) fit_func(x))
  teacher=min(result_old)
  pos=which.min(result_old)
  teacher_loc=school_old[pos,]
  mean_loc=colMeans(school_old)
  ranks=rank(result_old,ties.method = "random")
  result_new=matrix(ncol=classes,nrow=class_size)
  # optimisation algorithm update 
  k=1
  if(printall){
    results=list() # to be returned in the end
  }else{
    results=matrix(ncol=dim,nrow=1)
  }
  class_vec=1:class_size
  school_new=matrix(ncol=dim,nrow=class_size)
  reshuffle_prob=1
  stuck=F
  stuck_ind=0
  while((teacher>error) && (k<=gen) && !stuck){
    # two strategies in teachers phase, determined through a and b
    a=runif(class_size)
    b=runif(class_size)
    first=which(a<b)
    second=which(a>=b)
    while(length(first)<3&length(second)<3){
      a=runif(class_size)
      b=runif(class_size)
      first=which(a<b)
      second=which(a>=b)
    }
    # strategy a
    TF=round(1+runif(length(first)))
    teacher_loc_mat=matrix(rep(teacher_loc,length(first)),nrow=length(first),byrow=T)
    mean_loc_mat=matrix(rep(mean_loc,length(first)),nrow=length(first),byrow=T)           
    school_new[first,]=school_old[first,]+runif(length(first))*(teacher_loc_mat-TF*mean_loc_mat)
    # strategy b
    other=sample(class_vec[second])
    while(any(other==class_vec[second])){
      other=sample(class_vec[second])
    }
    # which is better
    other_better=which(result_old[second]>result_old[other])
    # make id based on result 
    secnd_id=second
    secnd_id[other_better]=other[other_better]
    teacher_loc_mat=matrix(rep(teacher_loc,length(second)),nrow=length(second),byrow=T)
    school_new[second,]=school_old[second,]+runif(length(second))*(teacher_loc_mat-school_old[secnd_id,])
    bound_max=which(school_new> xmax_mat)
    bound_min=which(school_new < xmin_mat)
    
    school_new[bound_max] = xmax_mat[bound_max]-(school_new[bound_max]-xmax_mat[bound_max])
    school_new[bound_min] = xmin_mat[bound_min]-(school_new[bound_min]-xmin_mat[bound_min])
    
    # parallel evaluation
    result_new=apply(school_new,1,function(x) fit_func(x))
    #
    newbetter=which(result_new<result_old)
    result_old[newbetter]=result_new[newbetter]
    school_old[newbetter,]=school_new[newbetter,]
    ranks=rank(result_old,ties.method = 'random')
    # end of teachers phase
    ############################
    #
    #############################
    
    # beginning of learners' phase 
    # two strategies in learners' phase, determined through a and b
    a=runif(class_size)
    b=runif(class_size)
    first=which(a<b)
    second=which(a>=b)
    while(length(first)<4&length(second)<4){
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
    while(any(other_two==class_vec[second]&other_two==other_one)){
      other_two=sample(class_vec[second])
    }
    id_other=which(result_old[other_one]<result_old[other_two])
    # make ids based on result for substraction 
    # when other is better it's school_old[other,]-school_old[first,] and vice versa
    one_better=other_two
    one_better[id_other]=other_one[id_other]
    two_better=other_one
    two_better[id_other]=other_two[id_other]
    #
    school_new[second,]=school_old[second,]+runif(length(second))*(school_old[one_better,]-school_old[two_better,])
    bound_max=which(school_new> xmax_mat)
    bound_min=which(school_new < xmin_mat)
    
    school_new[bound_max] = xmax_mat[bound_max]-(school_new[bound_max]-xmax_mat[bound_max])
    school_new[bound_min] = xmin_mat[bound_min]-(school_new[bound_min]-xmin_mat[bound_min])
    
    # parallel evaluation
    result_new=apply(school_new,1,function(x) fit_func(x))

    # replacing school with better solution
    newbetter=which(result_new<result_old)

    result_old[newbetter]=result_new[newbetter]
    school_old[newbetter,]=school_new[newbetter,]
    ranks=rank(result_old,ties.method = 'random')

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
    ranks=rank(result_old,ties.method = 'random')
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
    k=k+1
    reshuffle_prob=reshuffle_prob*start_shuffle_prob
    a=runif(1)
    # reinitilize worst half if prob reached
    if(a>reshuffle_prob){
      worst=length(ranks[which(ranks>floor(class_size/2))])
      school_old[which(ranks>floor(class_size/2)),]=matrix(t(runif(worst*dim,xmin,xmax)),nrow=worst,byrow=T) 
      reshuffle_prob=start_shuffle_prob
    }
    
  }
  if(printall){
    results[[k]]=cbind(t(teacher_loc),teacher)
  }else{
    resulttemp=cbind(t(teacher_loc),teacher)
    results=rbind(results,resulttemp)
  }

  return(results)
}

