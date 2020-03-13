library(readr)
library(Matrix)
library(data.tree)
library(viridis)
library(igraph)
library(dplyr)
library(data.table)
options(scipen=999999999)
##############################################################
# funtion to extract the three column (sparse matrix) format table for Hi-C data
df_out<-function(test_aws,res=25000){
  chr_bins<-unique(c(test_aws$X1,test_aws$X2))
  #include gaps in the heatmap
  full_range<-seq(range(chr_bins)[1],range(chr_bins)[2],by=res)
  miss_bins<-full_range[which(!(full_range %in% chr_bins))]
  full_chr_df<-rbind(test_aws,data.frame(X1=miss_bins,X2=miss_bins,X3=NA))
  
  #update the heatmap to include gaps
  full_chr_df<-full_chr_df%>%mutate(X3=ifelse(is.na(X3),0,ifelse(is.nan(X3),0,X3)))
  return(full_chr_df)
}


#############################################################
# Helper function to retrieve the path to root in growing tree
ff = function(x){ 
  if (class(x) == "list" & length(x)>0) 
    lapply(x, ff) 
  else 
    TRUE
}
# Recursive Louvain clustering
cl_bipart<-function(g_chr1){
  options(scipen=999999999)
  #initialisation
  #container to save cluster hierarchy as list of lists
  chr1_tree_list<-list()
  #container for cl members
  chr1_tree_cl<- list()
  
  #Louvain community detection
  modul_g<-cluster_louvain(g_chr1)
  mg_mem<-communities(modul_g)
  names(mg_mem)<-paste(lapply(mg_mem,length),lapply(mg_mem,function(x)min(as.numeric(x))),lapply(mg_mem,function(x)max(as.numeric(x))),sep='_')
  
  #save cluster membership,conductance and expansion
  chr1_tree_cl<-c(chr1_tree_cl,mg_mem)
  
  #temporary list of candidate cluster to further partition
  ok_part<-names(chr1_tree_cl)
  #initiate the tree
  for (i in ok_part){chr1_tree_list[[i]]<-list()}  
  
  #save path to all considered leaves
  lnames <- names(unlist(lapply(chr1_tree_list, ff)))
  names(lnames)<-names(chr1_tree_cl)
  
  #recursive looping
  while(length(ok_part)>0){
    temp_part<-c() 
    for(i in ok_part){
      #create the subnetwork of considered cluster
      sub_g1<- induced.subgraph(g_chr1,V(g_chr1)$name %in% chr1_tree_cl[[i]])
      #compute the louvain community detection
      modul_sg<-cluster_louvain(sub_g1)
      msg_mem<-communities(modul_sg)
      names(msg_mem)<-paste(lapply(msg_mem,length),lapply(msg_mem,function(x)min(as.numeric(x))),lapply(msg_mem,function(x)max(as.numeric(x))),sep='_')
      #skip to next iteration if no sub-clustering generated
      if(length(msg_mem)==1 & all(chr1_tree_cl[[i]] %in% unlist(msg_mem))){next}

      #print('cl append')
      
      #Only consider future cluster partition if they contain more than 1 loci
      ok_part_temp<-names(which(unlist(lapply(msg_mem,length))>1))
      for(k in ok_part_temp){chr1_tree_cl[[k]]<-msg_mem[[k]]}
      
      #print('tree growth')
      for (j in ok_part_temp){chr1_tree_list[[c(unlist(strsplit(lnames[i],split='\\.')),j)]]<-list()} 
      
      temp_part<-c(temp_part,ok_part_temp)
      rm(sub_g1,msg_mem,modul_sg)
    }
    print(length(temp_part))
    ok_part<-temp_part 
    #gather updated paths
    lnames <- names(unlist(lapply(chr1_tree_list, ff)))
    #name each path according to leaf of considered path
    names(lnames)<-lapply(lnames,function(x)unlist(strsplit(x,split='\\.'))[length(unlist(strsplit(x,split='\\.')))])
    
  }
  
  return(list(cl_member=chr1_tree_cl,part_tree=chr1_tree_list))  
  
}
#############################################################
chr_set<-paste0("chr",1:19)

for(i in chr_set){
  print(paste('load',i,'data'))
  # Import the Hi-C data for cell stages considered
  chr_CN<-read_delim(paste0("~/../../media/vipin/DISQUEDUR/research_associate_project/HIC_OE_compare/data/CN/",i,".txt"),"\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)
  chr_ES<-read_delim(paste0("~/../../media/vipin/DISQUEDUR/research_associate_project/HIC_OE_compare/data/ES/",i,".txt"),"\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)
  chr_NPC<-read_delim(paste0("~/../../media/vipin/DISQUEDUR/research_associate_project/HIC_OE_compare/data/NPC/",i,".txt"),"\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)
  
  print(paste('build',i,'dataframe'))
  #compute the matrices for the different datasets
  NPC_chr19<-data.table(as.data.frame(df_out(chr_NPC)))
  names(NPC_chr19)<-c('ego','alter','NPC')
  ES_chr19<-data.table(as.data.frame(df_out(chr_ES)))
  names(ES_chr19)<-c('ego','alter','ES')
  CN_chr19<-data.table(as.data.frame(df_out(chr_CN)))
  names(CN_chr19)<-c('ego','alter','CN')
  rm(chr_CN,chr_NPC,chr_ES)
  
  #########################################
  #build overall interaction dataframe  
  tot_df<-merge(ES_chr19, NPC_chr19, all = TRUE, by = c("ego",'alter'))
  rm(ES_chr19,NPC_chr19)
  tot_df<-merge(tot_df, CN_chr19, all = TRUE, by = c("ego",'alter'))
  rm(CN_chr19)
  #build the correlation wrt the different stereotypical
  #stereotypical profile correlation-> visualisation 
  #purposes only
  #ev_profile<-matrix(c(-1,0,1,1,0,-1,-1,1,-1),ncol=3,byrow = T,dimnames = list(c("up","down","peak"),c("Var1","Var2",'Var3')))

  #interactions exactly fitting the upnward trend
  print(paste('build',i,'deformation matrix'))
  #extract the deformation map containing interactions strictly following the desired profile       
  up_df<-tot_df%>%filter(CN>NPC & NPC>ES & ego!=alter)%>%mutate(weight=CN-ES)%>%dplyr::select(ego,alter,weight)%>%mutate(weight=weight/sum(.$weight))
  print(paste(i,'louvain clustering'))
  g_chr1<-graph_from_data_frame(up_df,directed=F)

  chr_spec_res<- cl_bipart(g_chr1)
  save(chr_spec_res,file = paste0("~/../../media/vipin/DISQUEDUR/research_associate_project/HIC_OE_compare/res/",i,'_clust.rda'))
  rm(chr_f_mat,chr_spec_res,g_chr1,tot_df,up_df,inter_bin)
}
