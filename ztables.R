pnorm(-3:3,0,1,FALSE)
1-pnorm(0.01,0,1,FALSE)

table1=matrix(0,31,10)
table2=matrix(0,31,10)

for (i in 1:31){
  for(j in 1:10){
    (i-1)/10+(j-1)/100
    table1[i,j]=1-pnorm((i-1)/10+(j-1)/100,0,1,FALSE)
    table2[i,j]=pnorm((i-1)/10+(j-1)/100,0,1,FALSE)
  }
}
