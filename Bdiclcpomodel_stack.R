#######################################################################################
######################### Best Dic Model = Bdicmod ####################################
# model selection based on WAIc, DIC and LCPO                               ##########
########################### (inla stack) ##############################################
#######################################################################################

#######################################################################################
######## Inputs to be introduced          #############################################
###### * y: Renspose Variable       ###################################################
###### * variables: All the explicative avriables to test                       #####
###### * data: your data      ########################################################
###### * n: number of model that we want to check in the output  ####################
###### * family: distribution family of the model (in case of occurrence data, binomial) #
###### * ...: arguments of       inla ###########################################
#######################################################################################
Bdiclcpomodel_stack<-function(resp, variables, datos, n, family="binomial",...)
{  
  #T?rminos que usaremos 
  sel.terms <- switch('terms',terms=variables)
  
  # todas las combinaciones de los m elementos de v
  comb.terms <- function(m, v=sel.terms) {
    if(m==0) return('resp ~ -1 + beta0')
    else {
      combis <- apply(combn(v, m), 2, paste, collapse=' + ')
      return(paste('resp ~ -1 + beta0', combis, sep=' + '))
    }
  }
  
  #Lista con todos los modelos posibles
  f.list <- unlist(sapply(0:length(sel.terms), comb.terms))
  
  # lanzamos cada uno de los modelos guardados en el objeto 'f.list' y nos quedamos con el DIC
  dic<-numeric()
  LCPO<-numeric()
  waic<-numeric()
  for(i in 1:length(f.list)){
    res =inla(formula = eval(parse(text=f.list[i])), family=family, data=datos, ...)
    dic[i] <- res$dic$dic
    LCPO[i] = -mean(log(res$cpo$cpo))
    waic[i]<-res$waic$waic
    print(c(i, dic[i], waic[i], LCPO[i]))
  }
  
  #mostramos los modelos ordenados SEG??N EL dic
  modelos_dic<-data.frame(f.list[order(dic)[1:n]], dic[order(dic)[1:n]], waic[order(dic)[1:n]], LCPO[order(dic)[1:n]])
  colnames(modelos_dic)<-c("Modelos", "Dic", "Waic", "LCPO")
  
  #mostramos los modelos ordenados SEG??N EL waic
  modelos_waic<-data.frame(f.list[order(waic)[1:n]], dic[order(waic)[1:n]], waic[order(waic)[1:n]], LCPO[order(waic)[1:n]])
  colnames(modelos_waic)<-c("Modelos", "Dic", "Waic", "LCPO")
  
  
  #mostramos los modelos ordenados SEG??N EL LCPO
  modelos_lcpo<-data.frame(f.list[order(LCPO)[1:n]], dic[order(LCPO)[1:n]], waic[order(LCPO)[1:n]], LCPO[order(LCPO)[1:n]])
  colnames(modelos_lcpo)<-c("Modelos", "Dic", "Waic", "LCPO")
  

  modelos<-list(modelos_dic, modelos_waic, modelos_lcpo)
  names(modelos)<-c("Modelos dic", "Modelos waic", "Modelos lcpo")
  modelos
  
}
