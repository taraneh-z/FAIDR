#bulk feature calculations for each IDR (in this case evolutionary signatures similar to those calculated in Zarin et al., eLife, 2019)
filename<-"evolsig_unconstrained_10orthol_50perc_non_na.txt"
#target values corresponding to each IDR (based on protein annotations; y=1 if protein is associated with function, y=0 if not)
targetfile<-"targets_unconstrained_evolsig.txt";
#installs the glmnet package if not already installed
if (!require("glmnet")) install.packages("glmnet")
library('glmnet');
data<-read.delim(filename);
data[is.na(data)]<-0; #set missing data to 0. R's regression models don't always work well with missing data.
targets<-read.delim(targetfile);

gene2row<-list();
gene2idr<-list();
for (i in 1:length(data$idr_name)) {
  gene<-strsplit(toString(data$idr_name[[i]]),'_')[[1]][1];
  idr<-as.numeric(strsplit(toString(data$idr_name[[i]]),'_')[[1]][3])
  gene2row[[gene]]<-c(gene2row[[gene]],i);
  gene2idr[[gene]]<-c(gene2idr[[gene]],idr);

}
X<-as.matrix(data[,3:length(data[1,])]);
initw<-NULL; #will store initial weights
multi<-NULL;
for (g in 1:length(gene2row)) {
  initw[gene2row[[g]][1]]<-1;
  if (length(gene2row[[g]])>1) { #if there is more than one idr for this gene ...
    multi<-c(multi,g);
    initw[gene2row[[g]]] <- 1/length(gene2row[[g]]); #downweight them all equally
  }
}
lam<-0.2;
Tstats<-matrix(nr=(length(X[1,])+1), nc=length(targets));
PostIDR<-matrix(nr=length(X[,1]), nc=length(targets));
Pred<-matrix(nr=length(gene2row), nc=length(targets));
geney<-matrix(nr=length(gene2row), nc=length(targets));
for (t in 3:length(targets)) {
  y<-targets[,t];
  w<-initw;
  z<-NULL;
  #infer IDR posterior probabilities in logistic regression model.
  niter<-20;
  pr<-NULL;
  mod<-glmnet(X,y,family="gaussian",weights=w); #intial model
  for (i in 1:niter) {
    for (j in 1:5) { #M-step: fit the weighted logistic regression by iterative least squares
      pr<-1/(1+exp(-predict(mod,X,s=lam)));
      w1<-pr*(1-pr);
      w1[w1==0]<-exp(-20); ### avoid 0s
      z<-predict(mod,X,s=lam) + (y - pr)/w1;
      mod<-glmnet(X,z,family="gaussian",weights=w1*w)
    }
    #E-step: recalculate weights given latest model
    for (g in multi) {
      rows<-gene2row[[g]];
      w[rows] <- y[rows]*pr[rows] + (1-y[rows])*(1-pr[rows])
       w[rows] <- w[rows]/sum(w[rows]);
    }
  }
  #final weights using a formula that doesn't depend on y
  for (g in multi) {
    rows<-gene2row[[g]];
    w[rows] <- (pr[rows]+0.05)/sum(pr[rows]+0.05) #post(x) given that this is a positive example (y=1), smoothed so small pr are ignored
  }
  PostIDR[,t]<-w;
  #store the final predictions
  for (g in 1:length(gene2row)) {
    geney[g,t]<-y[gene2row[[g]][1]];
    #predictions for each gene are the sum of the idrs in the gene
    Pred[g,t]<-sum(pr[gene2row[[g]]]*w[gene2row[[g]]]);
  }

  #find the non-zero coefficients
  cols<-abs(as.numeric(coef(mod,s=lam))[2:(length(X[1,])+1)])>0;

  if (sum(cols)>0) {
    subX<-X[,cols];
    Tstats[c(1,1+which(cols)),t]<-summary(glm(z~subX, weights=w*w1))$coefficients[,3];
  }
}
#Make dataframes from matrices
PostIDR_df <- as.data.frame(PostIDR[,3:ncol(PostIDR)])
colnames(PostIDR_df) <- colnames(targets)[3:ncol(targets)]
rownames(PostIDR_df) <- as.vector(data[,1])
Tstats_df <- as.data.frame(Tstats[2:nrow(Tstats), 3:ncol(Tstats)])
rownames(Tstats_df) <- colnames(data)[3:ncol(data)]
colnames(Tstats_df) <- colnames(targets)[3:ncol(targets)]
Pred_df <- as.data.frame(Pred[,3:ncol(Pred)])
rownames(Pred_df) <- names(gene2row)
colnames(Pred_df) <- colnames(targets)[3:ncol(targets)]
#write dataframes to files
write.csv(PostIDR_df, "PostIDR_df.csv")
write.csv(Tstats_df, "Tstats_df.csv")
write.csv(Pred_df, "Pred_df.csv")
