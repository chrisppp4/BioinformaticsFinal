library(class)
setwd("./")
wt <- read.csv("./WT_assay_clinical.csv")
diff <- data.frame(seq(42,21448))
colnames(diff) <- c('X')
for (i in seq(42,21448)) {
  temp = cor.test(wt[[i]],wt$Event.Free.Survival.Time.in.Days)
  diff$p.value[i-41] = temp$p.value
}
expression_set_0.1 <- c(which(abs(diff$p.value) > 0.1) + 41)
expression_set_0.05 <- c(which(abs(diff$p.value) > 0.05) + 41)
es_0.1 <- wt[,expression_set_0.1]
es_0.05 <- wt[,expression_set_0.05]
es_0 <- wt[,seq(42,21448)]
fourcl = unname(quantile(wt$Event.Free.Survival.Time.in.Days,c(0.25,0.5,0.75)))
for (x in seq(1,length(wt$Event.Free.Survival.Time.in.Days))) {
  if( wt$Event.Free.Survival.Time.in.Days[x] <= fourcl[1]) {
    wt$cohort[x] <- "level1"
  } else {
    if( wt$Event.Free.Survival.Time.in.Days[x] <= fourcl[2] ){
      wt$cohort[x] <- "level2"
    } else {
      if( wt$Event.Free.Survival.Time.in.Days[x] <= fourcl[3] ){
        wt$cohort[x] <- "level3"
      } else {
        wt$cohort[x] <- "level4"
      }
    }
  }
}
cl = matrix(nrow = 4, ncol = length(wt$Event.Free.Survival.Time.in.Days))
for (x in seq(1,length(wt$Event.Free.Survival.Time.in.Days))) {
  cl[1,x] = FALSE
  cl[2,x] = FALSE
  cl[3,x] = FALSE
  cl[4,x] = FALSE
  if (wt$cohort[x] == "level1") {
    cl[1,x] = TRUE
  }
  if (wt$cohort[x] == "level2") {
    cl[2,x] = TRUE
  }
  if (wt$cohort[x] == "level3") {
    cl[3,x] = TRUE
  }
  if (wt$cohort[x] == "level4") {
    cl[4,x] = TRUE
  }
}
train = es_0.1[1:66,]
test = es_0.1[67:nrow(es_0.1),]

ovrknn <- function(train,test,cl,kvalue){
  trueresult = wt$cohort[67:132]
  nnlevel1 = as.logical.factor(knn(train,test,cl[1,1:66],k=kvalue))
  nnlevel2 = as.logical.factor(knn(train,test,cl[2,1:66],k=kvalue))
  nnlevel3 = as.logical.factor(knn(train,test,cl[3,1:66],k=kvalue))
  nnlevel4 = as.logical.factor(knn(train,test,cl[4,1:66],k=kvalue))
  testresult = vector(length = 66)
  cllabel = c("level1","level2","level3","level4")
  for (i in seq(1:66)) {
    groupclass = c(nnlevel1[i],nnlevel2[i],nnlevel3[i],nnlevel4[i])
    selection = which(groupclass)
    if (length(selection) == 0) {
      randomselection = ceiling(runif(1,0.0001,4))
      testresult[i] = cllabel[randomselection]
    }
    if (length(selection) == 1) {
      testresult[i] = cllabel[selection[1]]
    }
    if (length(selection) > 1) {
      randomselection = ceiling(runif(1,0.0001,length(selection)))
      testresult[i] = cllabel[selection[randomselection]]
    }
  }
  accuracy = length(which(c(testresult==trueresult))) / 66
  return(accuracy)
}
print("One vs Rest:")
print("k=3,pvalue=0.1")
print(ovrknn(train,test,cl,3))
print("k=3,pvalue=0.05")
print(ovrknn(es_0.05[1:66,],es_0.05[67:132,],cl,3))
print("k=3,pvalue=0")
print(ovrknn(es_0[1:66,],es_0[67:132,],cl,3))
print("k=4,pvalue=0.1")
print(ovrknn(train,test,cl,4))
print("k=4,pvalue=0.05")
print(ovrknn(es_0.05[1:66,],es_0.05[67:132,],cl,4))
print("k=4,pvalue=0")
print(ovrknn(es_0[1:66,],es_0[67:132,],cl,4))
print("k=5,pvalue=0.1")
print(ovrknn(train,test,cl,5))
print("k=5,pvalue=0.05")
print(ovrknn(es_0.05[1:66,],es_0.05[67:132,],cl,5))
print("k=5,pvalue=0")
print(ovrknn(es_0[1:66,],es_0[67:132,],cl,5))

train_group1 = list(matrix(nrow = 17, ncol = ncol(es_0.1)))
cl_group1 = list(matrix(nrow = 17, ncol = 1))
test_group1 = list(matrix(nrow = 16, ncol = ncol(es_0.1)))
trueresult_group1 = list(matrix(nrow = 16, ncol = 1))

train_group1[[1]] = es_0.1[which(wt$cohort=="level1")[1:17],]
cl_group1[[1]] = wt$cohort[which(wt$cohort=="level1")[1:17]]
test_group1[[1]] = es_0.1[which(wt$cohort=="level1")[18:33],]
trueresult_group1[[1]] = wt$cohort[which(wt$cohort=="level1")[18:33]]

train_group1[[2]] = es_0.1[which(wt$cohort=="level2")[1:17],]
cl_group1[[2]] = wt$cohort[which(wt$cohort=="level2")[1:17]]
test_group1[[2]] = es_0.1[which(wt$cohort=="level2")[18:33],]
trueresult_group1[[2]] = wt$cohort[which(wt$cohort=="level2")[18:33]]

train_group1[[3]] = es_0.1[which(wt$cohort=="level3")[1:17],]
cl_group1[[3]] = wt$cohort[which(wt$cohort=="level3")[1:17]]
test_group1[[3]] = es_0.1[which(wt$cohort=="level3")[18:33],]
trueresult_group1[[3]] = wt$cohort[which(wt$cohort=="level3")[18:33]]

train_group1[[4]] = es_0.1[which(wt$cohort=="level4")[1:17],]
cl_group1[[4]] = wt$cohort[which(wt$cohort=="level4")[1:17]]
test_group1[[4]] = es_0.1[which(wt$cohort=="level4")[18:33],]
trueresult_group1[[4]] = wt$cohort[which(wt$cohort=="level4")[18:33]]

train_group2 = list(matrix(nrow = 17, ncol = ncol(es_0.1)))
cl_group2 = list(matrix(nrow = 17, ncol = 1))
test_group2 = list(matrix(nrow = 16, ncol = ncol(es_0.1)))
trueresult_group2 = list(matrix(nrow = 16, ncol = 1))

train_group2[[1]] = es_0.05[which(wt$cohort=="level1")[1:17],]
cl_group2[[1]] = wt$cohort[which(wt$cohort=="level1")[1:17]]
test_group2[[1]] = es_0.05[which(wt$cohort=="level1")[18:33],]
trueresult_group2[[1]] = wt$cohort[which(wt$cohort=="level1")[18:33]]

train_group2[[2]] = es_0.05[which(wt$cohort=="level2")[1:17],]
cl_group2[[2]] = wt$cohort[which(wt$cohort=="level2")[1:17]]
test_group2[[2]] = es_0.05[which(wt$cohort=="level2")[18:33],]
trueresult_group2[[2]] = wt$cohort[which(wt$cohort=="level2")[18:33]]

train_group2[[3]] = es_0.05[which(wt$cohort=="level3")[1:17],]
cl_group2[[3]] = wt$cohort[which(wt$cohort=="level3")[1:17]]
test_group2[[3]] = es_0.05[which(wt$cohort=="level3")[18:33],]
trueresult_group2[[3]] = wt$cohort[which(wt$cohort=="level3")[18:33]]

train_group2[[4]] = es_0.05[which(wt$cohort=="level4")[1:17],]
cl_group2[[4]] = wt$cohort[which(wt$cohort=="level4")[1:17]]
test_group2[[4]] = es_0.05[which(wt$cohort=="level4")[18:33],]
trueresult_group2[[4]] = wt$cohort[which(wt$cohort=="level4")[18:33]]

train_group3 = list(matrix(nrow = 17, ncol = ncol(es_0.1)))
cl_group3 = list(matrix(nrow = 17, ncol = 1))
test_group3 = list(matrix(nrow = 16, ncol = ncol(es_0.1)))
trueresult_group3 = list(matrix(nrow = 16, ncol = 1))

train_group3[[1]] = es_0[which(wt$cohort=="level1")[1:17],]
cl_group3[[1]] = wt$cohort[which(wt$cohort=="level1")[1:17]]
test_group3[[1]] = es_0[which(wt$cohort=="level1")[18:33],]
trueresult_group3[[1]] = wt$cohort[which(wt$cohort=="level1")[18:33]]

train_group3[[2]] = es_0[which(wt$cohort=="level2")[1:17],]
cl_group3[[2]] = wt$cohort[which(wt$cohort=="level2")[1:17]]
test_group3[[2]] = es_0[which(wt$cohort=="level2")[18:33],]
trueresult_group3[[2]] = wt$cohort[which(wt$cohort=="level2")[18:33]]

train_group3[[3]] = es_0[which(wt$cohort=="level3")[1:17],]
cl_group3[[3]] = wt$cohort[which(wt$cohort=="level3")[1:17]]
test_group3[[3]] = es_0[which(wt$cohort=="level3")[18:33],]
trueresult_group3[[3]] = wt$cohort[which(wt$cohort=="level3")[18:33]]

train_group3[[4]] = es_0[which(wt$cohort=="level4")[1:17],]
cl_group3[[4]] = wt$cohort[which(wt$cohort=="level4")[1:17]]
test_group3[[4]] = es_0[which(wt$cohort=="level4")[18:33],]
trueresult_group3[[4]] = wt$cohort[which(wt$cohort=="level4")[18:33]]

ovoknn <- function(train1,test1,cl1,kvalue1,trueresult1){
  test_in_all = rbind(test1[[1]],test1[[2]],test1[[3]],test1[[4]])
  nnsubset1 = as.character.factor(knn(rbind(train1[[1]],train1[[2]]),test_in_all,rbind(cl1[[1]],cl1[[2]]),kvalue1))
  nnsubset2 = as.character.factor(knn(rbind(train1[[1]],train1[[3]]),test_in_all,rbind(cl1[[1]],cl1[[3]]),kvalue1))
  nnsubset3 = as.character.factor(knn(rbind(train1[[1]],train1[[4]]),test_in_all,rbind(cl1[[1]],cl1[[4]]),kvalue1))
  nnsubset4 = as.character.factor(knn(rbind(train1[[2]],train1[[3]]),test_in_all,rbind(cl1[[2]],cl1[[3]]),kvalue1))
  nnsubset5 = as.character.factor(knn(rbind(train1[[2]],train1[[4]]),test_in_all,rbind(cl1[[2]],cl1[[4]]),kvalue1))
  nnsubset6 = as.character.factor(knn(rbind(train1[[3]],train1[[4]]),test_in_all,rbind(cl1[[3]],cl1[[4]]),kvalue1))
  testresult = c()
  for (i in seq(1:64)) {
    temp = c(nnsubset1[i],nnsubset2[i],nnsubset3[i],nnsubset4[i],nnsubset5[i],nnsubset6[i])
    temp_table = as.data.frame(table(temp))
    max_freq = temp_table[which(temp_table$Freq==max(temp_table$Freq)),]
    if (length(max_freq$temp)==1) {
      testresult[i] = as.character(max_freq$temp)
    }
    else {
      testresult[i] = as.character(max_freq$temp[ceiling(runif(1,0.0001,length(max_freq)))])
    }
  }
  trueresult_set = c(trueresult1[[1]],trueresult1[[2]],trueresult1[[3]],trueresult1[[4]])
  accuracy1 = length(which(c(testresult==trueresult_set))) / 64
  return(accuracy1)
}
print("One vs One:")
print("k=3,4,5 pvalue=0.1")
print(ovoknn(train_group1,test_group1,cl_group1,3,trueresult_group1))
print(ovoknn(train_group1,test_group1,cl_group1,4,trueresult_group1))
print(ovoknn(train_group1,test_group1,cl_group1,5,trueresult_group1))
print("k=3,4,5 pvalue=0.05")
print(ovoknn(train_group2,test_group2,cl_group2,3,trueresult_group2))
print(ovoknn(train_group2,test_group2,cl_group2,4,trueresult_group2))
print(ovoknn(train_group2,test_group2,cl_group2,5,trueresult_group2))
print("k=3,4,5 pvalue=0")
print(ovoknn(train_group3,test_group3,cl_group3,3,trueresult_group3))
print(ovoknn(train_group3,test_group3,cl_group3,4,trueresult_group3))
print(ovoknn(train_group3,test_group3,cl_group3,5,trueresult_group3))