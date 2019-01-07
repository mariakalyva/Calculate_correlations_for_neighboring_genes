#This is an R function to calculate correlations between genomic neighborhoods
#Data: nucleosome data- 100 bits

nucNei<-function(df1,df2){
  # df1 <- data.frame with names of genes
  # df2 <- data.frame with nucleosome data in 100 bits of info
  # Example run: nucCor(GeneNames,Nucleosome_Data)
  #Nucleosome_data look like
  #GENE CHR START END -100Columns of scores-
  #GeneNames is one column data.frame with gene names
  mm<-c()
  data_correlations <- c()
  m<-nrow(df2)
  n<-nrow(df1)
  
  for(i in 1:n){
    for(j in 6:(m-5)){
      
      #test if name matches
      if(df1[i,1] == df2[j,1]){
        
        #calculate correlations if genes in the same chrom
        if(df2[j,2]==df2[j-1,2]){
          cc<-cor(t(df2[j,5:104]),t(df2[j-1,5:104]))
        }
        else { cc<-0}
        if(df2[j,2]==df2[j-2,2]){
          cc1<-cor(t(df2[j,5:104]),t(df2[j-2,5:104]))
        }
        else { cc1<-0}
        if(df2[j,2]==df2[j-3,2]){
          cc2<-cor(t(df2[j,5:104]),t(df2[j-3,5:104]))
        }
        else { cc2<-0}
        if(df2[j,2]==df2[j-4,2]){
          cc3<-cor(t(df2[j,5:104]),t(df2[j-4,5:104]))
        }
        else { cc3<-0}
        if(df2[j,2]==df2[j-5,2]){
          cc4<-cor(t(df2[j,5:104]),t(df2[j-5,5:104]))
        }
        else { cc4<-0}
        if(df2[j,2]==df2[j+1,2]){
          cc5<-cor(t(df2[j,5:104]),t(df2[j+1,5:104]))
        }
        else { cc5<-0}
        if(df2[j,2]==df2[j+2,2]){
          cc6<-cor(t(df2[j,5:104]),t(df2[j+2,5:104]))
        }
        else { cc6<-0}
        if(df2[j,2]==df2[j+3,2]){
          cc7<-cor(t(df2[j,5:104]),t(df2[j+3,5:104]))
        }
        else { cc7<-0}
        if(df2[j,2]==df2[j+4,2]){
          cc8<-cor(t(df2[j,5:104]),t(df2[j+4,5:104]))
        }
        else { cc8<-0}
        if(df2[j,2]==df2[j+5,2]){
          cc9<-cor(t(df2[j,5:104]),t(df2[j+5,5:104]))
        }
        else { cc9<-0}
        
        d<-c(cc,cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,cc9)
        
        
      }
      else{next()}
    }
    mm[i]<-mean(d)
  }
  data_correlations <<- mm
}

