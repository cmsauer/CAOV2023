#Code based on Lena Morrill's code modified by Carolin (27/12/2021)

getwd()

library(compositions)
library(TMB)
library(gridExtra)
library(ggplot2)
library(dplyr)

#-------------------------------------------------------------------------------------------#
## the model, in cpp
folder_of_TMB_model <- "tmb_RE/" ## YOU CHANGE THIS! RELATIVE PATH TO TMB MODEL (IN CPP)
TMB::compile(paste0(folder_of_TMB_model, "mvn_beta_no_cor.cpp"), "-std=gnu++17")
dyn.load(dynlib(paste0(folder_of_TMB_model, "../tmb_RE/mvn_beta_no_cor")))
#-------------------------------------------------------------------------------------------#



#-------------------------------------------------------------------------------------------#
## Functions

impute = function(mat, inputation_value){
  mat[mat == 0] = inputation_value
  normalise_rw(mat)
}

normalise_rw <- function(x){
  if(is.null(dim(x))){
    x/sum(x)
  }else{
    ## normalise row-wise
    sweep(x, 1, rowSums(x), '/')
  }
}

python_like_select_name = function(vector, grep_substring){
  vector[grepl(pattern = grep_substring, x = names(vector))]
}

plot_betas <- function(TMB_obj, names_cats=NULL, rotate_axis=T, theme_bw=T, remove_SBS=T, only_slope=F, return_df=F, plot=T,
                       line_zero=T, add_confint=F, return_plot=T, return_ggplot=F, title=NULL, add_median=F){
  if(typeof(TMB_obj) == 'character'){
    .summary_betas <- NA
    if(theme_bw){
      plt <- ggplot()+theme_bw()
      if(plot) print(plt)
    }else{
      plt <- ggplot()
      if(plot) print(plt)
    }
  }else{
    .summary_betas <- summary(TMB_obj)
    .summary_betas <- cbind.data.frame(python_like_select_rownames(.summary_betas, 'beta'),
                                       type_beta=rep(c('Intercept', 'Slope')),
                                       LogR=rep(1:(nrow(python_like_select_rownames(.summary_betas, 'beta'))/2), each=2))
    if(only_slope){
      .summary_betas <- .summary_betas[.summary_betas$type_beta == 'Slope',]
    }
    
    if(!is.null(names_cats)){
      if(remove_SBS){
        names_cats <- gsub("SBS", "", names_cats) 
      }
      if(length(unique(.summary_betas$LogR)) != length(names_cats)){
        stop('Number of beta slope/intercept pairs should be the same as the length of the name of the categories')
      }
      .summary_betas$LogR = names_cats[.summary_betas$LogR]
    }
    plt <- ggplot(.summary_betas, aes(x=LogR, y=`Estimate`))
    
    if(line_zero) plt <- plt + geom_hline(yintercept = 0, lty='dashed', col='blue')
    if(add_median) { plt <- plt +
      geom_hline(yintercept = median(c(0,.summary_betas$Estimate[.summary_betas$type_beta == 'Slope'])),
                 lty='dashed', col='red') }
    
    plt <- plt +
      geom_point()+
      geom_errorbar(aes(ymin=`Estimate`-`Std. Error`, ymax=`Estimate`+`Std. Error`), width=.1)+
      ggtitle('Slopes')+facet_wrap(.~type_beta, scales = "free")
    
    if(theme_bw){
      plt <- plt + theme_bw()
    }
    
    if(add_confint){
      confints <- cbind(.summary_betas, confint=t(give_confidence_interval(.summary_betas[,'Estimate'], .summary_betas[,'Std. Error'])))
      plt <- plt+
        geom_errorbar(data = confints, aes(ymin=confint.1, ymax=confint.2), width=.1 ,col='blue', alpha=0.6)
      
    }
    
    if(rotate_axis){
      plt <- plt + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    }
    
    if(!is.null(title)){
      plt <- plt + ggtitle(title)
    }
    
    if(!TMB_obj$pdHess){
      plt <- plt + annotate("text", x = -Inf, y=Inf, label="not PD", vjust=1, hjust=-.2)+geom_point(col='red')
      if(plot) print(plt)
    }else{
      if(plot) print(plt)
    }
  }
  
  if(return_df){
    .summary_betas
  }else{
    if(return_plot & return_df){stop('<return_plot=T> and <return_df=T> are incompatible')}
    plot_list <- list(plt)
    class(plot_list) <- c("quiet_list", class(plot_list))
    if(return_plot){
      return(cowplot::as_grob(plt))
    }else if(return_ggplot){
      return(plt)
    }
  }
}

softmax = function(x){
  if(is.null(dim(x))){
    ## vector
    sum_x = sum(exp(x))
    exp(x)/sum_x
  }else{
    ## matrix
    sum_x = rowSums(exp(x))
    sweep(exp(x), 1, sum_x, '/')
  }
}

python_like_select_rownames = function(matrix, grep_substring){
  matrix[grepl(pattern = grep_substring, x = rownames(matrix)),]
}

select_slope_2 = function(i, verbatim=TRUE){
  if(is.null(dim(i))){
    i[c(F,T)]
  }else{
    i[,c(F,T)]
  }
}

#-------------------------------------------------------------------------------------------#


#-------------------------------------------------------------------------------------------#
## Load and prepare data to analyse

cns_if <- read.csv("../data/cell_line_genomics_data/cns_IF_data_forALRmodel_CMS20211227.csv") %>% 
  filter(!is.na(X..CA.cells)) %>% 
  arrange(sample_id)

# metadata (groups)
if_data <- cns_if %>% 
  select(X..CA.cells, CA_group, X..MN.cells, MN_group) %>% 
  mutate(CA_group = case_when(CA_group == "high" ~ 1, 
                              CA_group == "low" ~ 0),
         MN_group = case_when(MN_group == "high" ~ 1, 
                              MN_group == "low" ~ 0))
rownames(if_data) <- cns_if$sample_id

#cns exposure data
exposures <- cns_if[,2:8] 
rownames(exposures) <- cns_if$sample_id

n <- nrow(exposures)
d <- ncol(exposures)

image(t(exposures))


## Add imputation, if needed, Alternative ways of imputating data: zCompositions::multLN(), and others
exposures <- impute(exposures, 1e-2)
which_zero = t(apply(exposures, 1, function(i) as.numeric((i==0)) ))

ALR_transformed_data = as(compositions::alr(exposures), 'matrix')
image(t(ALR_transformed_data)) ## absolute exposures

dmin1 = d-1
#-------------------------------------------------------------------------------------------#

TMB_data = list(Y = ALR_transformed_data,
                d = dmin1,
                n = n,
                x = cbind(1,if_data$CA_group))

TMB_params = list(logs_sd=runif(n = dmin1, min = 0, max = 2),
                  beta = (matrix(runif(dmin1*2, min = -4, max = 4),
                                 nrow = 2, byrow=TRUE))
) ## some initial parameters. The parameter that we are most interested in is the second row of beta (the coefficients for the
## change between the two conditions)

#-------------------------------------------------------------------------------------------#
## Run model

obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="mvn_beta_no_cor")
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
rep <- sdreport(obj)
rep

# png("CA_highVSlow_TMB.png")
plot_betas(rep)
# dev.off()

#-------------------------------------------------------------------------------------------#
## If you wanted to run regression with some continuous value, instead of 

TMB_data_continuous = list(Y = ALR_transformed_data,
                           d = dmin1,
                           n = n,
                           x = cbind(1, if_data$X..CA.cells))

# obj <- MakeADFun(data = TMB_data_continuous, parameters = TMB_params, DLL="tmb_MVN_partial_ILR_FEb")
obj <- MakeADFun(data = TMB_data_continuous, parameters = TMB_params, DLL="mvn_beta_no_cor")
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
rep_cont <- sdreport(obj)
rep_cont

png("CA_cont_TMB.png")
plot_betas(rep_cont) 
dev.off()
#-------------------------------------------------------------------------------------------#

