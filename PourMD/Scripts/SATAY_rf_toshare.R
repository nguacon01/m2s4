library(caret) #validation croisée
library(randomForest) #forets aleatoires
library(pROC) #courbe ROC
library(stats) #predictions
library(FactoMineR)
library(factoextra)
library(BBmisc)


#generation de la foret
rf=randomForest(as.factor(V8)~V2+V3+V4+V5+V6+V7,data=trainingene)
print(rf) #matrice de confusion

#validation croisée
train.control <- trainControl(method = "cv",number=5)
model=train(as.factor(V8)~V2+V3+V4+V5+V6+V7,data=trainingene,method="rf",trControl=train.control)

#courbe roc
score <-predict(rf,trainingene,type="prob")[,"1"]
roc_obj <-roc(trainingene$V8=="1",score)
plot(1-roc_obj$specificities,roc_obj$sensitivities,type="l",xlab="False Positive",ylab="True Positive")
abline(0,1)
#AUC
roc_obj$auc

#identification du seuil optimal
i=which(max(roc_obj$sensitivities+roc_obj$specificities)==roc_obj$sensitivities+roc_obj$specificities)
roc_obj$thresholds[i] #seuil optimal
points(1-roc_obj$specificities[i],roc_obj$sensitivities[i],col="green") #affichage sur la courbe

1-roc_obj$specificities[i] #FPR
roc_obj$sensitivities[i] #TPR


#predictions sur le set complet
V=predict(rf,genefull,type="response")
genefull[,length(genefull)+1]=V

write.table(CLQCA20184_genes,"/Users/fabien/SATAY/mapping_july2019/CLQCA20184/CLQCA20184_genes_predE.txt", sep="\t",row.names=FALSE,col.names=FALSE)



#subsampling

Genez=BBmisc::normalize(trainingene[,2:7]) #standardisation des variables quantitatives
Genez[,7:8]=c(as.character(trainingene$V1),trainingene$V8) #ajout des variables qualitatives
cha=HCPC(Genez[,1:6],-1,graph=F) #clustering
fviz_cluster(cha,data = Genez[,1:6]) #visualisation

#separation en classes
classe1=Genez[which(cha$data.clust$clust==1),] 
classe2=Genez[which(cha$data.clust$clust==2),]
classe3=Genez[which(cha$data.clust$clust==3),]

#initialisation du trainingset
tsetg=data.frame(Genez[1,],0)
tsetg=tsetg[-1,]

listeclass=list(classe1,classe2,classe3)

completerts=function(tset,listeclass,set){
  p=c()
  it=1
  dist=list(c())
  for(i in 1:length(listeclass)){ #matrice des distances entre chaque individu de chaque classe
    dist[[i]]=as.matrix(get_dist(listeclass[[i]]))
  }
  for(i in 1:length(listeclass)){ #proportion de chaque classe
    p[i]=length(listeclass[[i]][,1])/length(set[,1])
  }
  error=c()
  precision0=-1
  convergence=F
  while(convergence==F){ 
    for(i in 1:length(p)){ #ajout des individus au training set
      for(j in 1:ceiling(30*p[i])){
        ind=which(dist[[i]]==sort(dist[[i]],decreasing=T)[1],arr.ind = T)
        ind1=ind[1,1]
        ind2=ind[1,2]
        gene1=trainingene[which(trainingene$V1==Genez$V7.1[ind1]),]
        gene2=trainingene[which(trainingene$V1==Genez$V7.1[ind2]),]
        tsetg=rbind(tsetg,gene1)
        tsetg=rbind(tsetg,gene2)
        dist[[i]]=dist[[i]][-c(ind1,ind2),-c(ind1,ind2)]
      }
    }
    foret=randomForest(as.factor(V8)~V2+V3+V4+V5+V6+V7,data=tsetg) #prédictions
    precision1=foret$err.rate[length(foret$err.rate[,1]),1] #taux d'erreurs
    print(precision1)
    error[it]=precision1
    if(abs(precision1-precision0)<0.001){convergence=T} #verification de la "stabilité" du taux d'erreur
    precision0=precision1
    it=it+1
    for(k in 1:length(p)){ #pénalisation des individus mal prédits
      pred=predict(foret,listeclass[[k]],type="response")
      ind=which(listeclass[[k]]$V8!=pred)
      listeclass[[k]][ind,2]=listeclass[[k]][ind,2]+1000
      }
    }
  return(list(tsetg,error))
}

rf=randomForest(as.factor(V8)~V2+V3+V4+V5+V6+V7,data=trainingsetgene) #random forest sur le nouveau training set
V=predict(rf,trainingene,type="response")

#courbe ROC etc
score <-predict(rf,trainingene,type="prob")[,"1"]
roc_obj <-roc(trainingene$V8=="1",score)
plot(1-roc_obj$specificities,roc_obj$sensitivities,type="l",xlab="False Positive",ylab="True Positive")
