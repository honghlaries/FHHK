


obs.var2.5 = c('Ach1',  'Ach2',  'Amb1',  'Amb2',  'Amb3')
R.prob2.5 = matrix(c(
  1.00 ,  .60  , .30,  .20,   .20,                                               
  .60,  1.00,   .20,   .30,   .10,                                                
  .30,   .20,  1.00,   .70,   .60 ,                                               
  .20,   .30,   .70,  1.00,   .50,                                                
  .20,   .10,   .60,  .50,  1.00), ncol=5,byrow=TRUE)   
model2.5=matrix(c(
  'Ambit ->  Amb1',      'a', NA,
  'Ambit -> Amb2' ,      'b', NA,
  'Ambit -> Amb3' ,      'c', NA,
  'Achieve -> Ach1',     'd', NA,
  'Achieve -> Ach2',     'e', NA,
  'Ambit <-> Achieve',   'f', NA,
  'Amb1 <-> Amb1' ,      'u', NA,
  'Amb2 <-> Amb2' ,      'v', NA,
  'Amb3 <-> Amb3' ,      'w', NA,
  'Ach1 <-> Ach1' ,      'x', NA,
  'Ach2 <-> Ach2' ,      'y', NA,
  'Achieve <-> Achieve',  NA, 1,
  'Ambit <-> Ambit',      NA, 1),
  ncol=3, byrow=TRUE)
sem2.5= sem(model2.5,R.prob2.5,60, obs.var2.5)
summary(sem2.5,digits=3)