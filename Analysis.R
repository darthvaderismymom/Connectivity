setwd("~/Рабочий стол/R data/resting-state dataset") 

data_O <- read.csv('data_O.csv', sep = ';', header = TRUE, stringsAsFactors = FALSE, dec = ',')
data_C <- read.csv('data_C.csv', sep = ';', header = TRUE, stringsAsFactors = FALSE, dec = ',')

############################## Cw Lw correlations with Raven #######################################
library(Hmisc)
corstarsl <- function(x){ 
  require(Hmisc) 
  x <- as.matrix(x) 
  R <- rcorr(x)$r 
  p <- rcorr(x)$P 
  p<-p.adjust(p, method = 'fdr')
  ## define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", " ")))
  
  ## trunctuate the matrix that holds the correlations to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1] 
  
  ## build a new matrix that includes the correlations with their apropriate stars 
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x)) 
  diag(Rnew) <- paste(diag(R), " ", sep="") 
  rownames(Rnew) <- colnames(x) 
  colnames(Rnew) <- paste(colnames(x), "", sep="") 
  
  ## remove upper triangle
  Rnew <- as.matrix(Rnew)
  Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
  Rnew <- as.data.frame(Rnew) 
  
  ## remove last column and return the matrix (which is now a data frame)
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  return(Rnew) 
}

# cw open
cwo <- data_O[,c('pcrvtot', 'Cw_normalized_alpha', 'Cw_normalized_alpha_1', 'Cw_normalized_alpha_2', 'Cw_normalized_beta_1',
                  'Cw_normalized_beta_2', 'Cw_normalized_theta' )]
CwOpenCorr <- corstarsl(cwo)

# cw closed
cwc <- data_C[,c('pcrvtot', 'Cw_normalized_alpha', 'Cw_normalized_alpha_1', 'Cw_normalized_alpha_2', 'Cw_normalized_beta_1',
             'Cw_normalized_beta_2', 'Cw_normalized_theta' )]
CwClosedCorr <- corstarsl(cwc)

# lw open
lwo <- data_O[,c('pcrvtot', 'Lw_normalized_alpha', 'Lw_normalized_alpha_1', 'Lw_normalized_alpha_2', 'Lw_normalized_beta_1',
                  'Lw_normalized_beta_2', 'Lw_normalized_theta' )]
LwOpenCorr <- corstarsl(lwo)

# lw closed
lwc <- data_C[,c('pcrvtot', 'Lw_normalized_alpha', 'Lw_normalized_alpha_1', 'Lw_normalized_alpha_2', 'Lw_normalized_beta_1',
             'Lw_normalized_beta_2', 'Lw_normalized_theta' )]
LwClosedCorr <- corstarsl(lwc)

################################# Cw Lw Raven correlations adjusted for WM ##############################
# cw open
cwo <- data_O[,c('pure_IQ', 'Cw_normalized_alpha', 'Cw_normalized_alpha_1', 'Cw_normalized_alpha_2', 'Cw_normalized_beta_1',
               'Cw_normalized_beta_2', 'Cw_normalized_theta' )]
CwOpenCorrPure <- corstarsl(cwo)

# cw closed
cwc <- data_C[,c('pure_IQ', 'Cw_normalized_alpha', 'Cw_normalized_alpha_1', 'Cw_normalized_alpha_2', 'Cw_normalized_beta_1',
                 'Cw_normalized_beta_2', 'Cw_normalized_theta' )]
CwClosedCorrPure <- corstarsl(cwc)

# lw open
lwo <- data_O[,c('pure_IQ', 'Lw_normalized_alpha', 'Lw_normalized_alpha_1', 'Lw_normalized_alpha_2', 'Lw_normalized_beta_1',
               'Lw_normalized_beta_2', 'Lw_normalized_theta' )]
LwOpenCorrPure <- corstarsl(lwo)

# lw closed
lwc <- data_C[,c('pure_IQ', 'Lw_normalized_alpha', 'Lw_normalized_alpha_1', 'Lw_normalized_alpha_2', 'Lw_normalized_beta_1',
                 'Lw_normalized_beta_2', 'Lw_normalized_theta' )]
LwClosedCorrPure <- corstarsl(lwc)

############################### здесь будет алгоритм корреляции GFP + Raven ###########################

###################################### Partial Correlations #########################################
library(ggm)
"pcor" <-
  function (u, S) 
  {
    ### Partial correlation between u[1:2], given th rest of u. S: cov matrix.
    k <- solve(S[u,u])
    -k[1,2]/sqrt(k[1,1]*k[2,2])
  }

"pcor.test" <-
  function(r, q, n){
    df = n - 2 - q
    tval <- r * sqrt(df)/sqrt(1-r*r)
    pv <- 2 * pt(-abs(tval), df)
    list(tval = tval, df = df, pvalue = pv)
    
  }

# closed eyes
# функция, output: pvalue, размер корр
partial_test <- function(x, y)
{
  result <- pcor(c(x, "pcrvtot", y), var(data_C))
  pc <- pcor.test(result, 1, 75)
  pc <- c(pc$pvalue, result)
  pc
}

# Lw x Raven controlling for Cw
partial_test("Lw_normalized_alpha", "Cw_normalized_alpha") #alpha
partial_test("Lw_normalized_alpha_1", "Cw_normalized_alpha_1") #alpha_1
partial_test("Lw_normalized_alpha_2", "Cw_normalized_alpha_2") #alpha_2
partial_test("Lw_normalized_beta_1", "Cw_normalized_beta_2") #beta_1
partial_test("Lw_normalized_beta_2", "Cw_normalized_beta_2") #beta_2
partial_test("Lw_normalized_theta", "Cw_normalized_theta") #theta

# Cw x Raven controlling for Lw
partial_test("Cw_normalized_alpha", "Lw_normalized_alpha") #alpha pval=0.056
partial_test("Cw_normalized_alpha_1", "Lw_normalized_alpha_1") #alpha_1 pval=0.030
partial_test("Cw_normalized_alpha_2", "Lw_normalized_alpha_2") #alpha_2 pval=0.034
partial_test("Cw_normalized_beta_1", "Lw_normalized_beta_2") #beta_1
partial_test("Cw_normalized_beta_2", "Lw_normalized_beta_2") #beta_2
partial_test("Cw_normalized_theta", "Lw_normalized_theta") #theta

#open eyes
# функция, output: pvalue, размер корр
partial_test <- function(x, y)
{
  result <- pcor(c(x, "pcrvtot", y), var(data_O))
  pc <- pcor.test(result, 1, 65)
  pc <- c(pc$pvalue, result)
  pc
}

# Lw x Raven controlling for Cw
partial_test("Lw_normalized_alpha", "Cw_normalized_alpha") #alpha pval=0.37
partial_test("Lw_normalized_alpha_1", "Cw_normalized_alpha_1") #alpha_1
partial_test("Lw_normalized_alpha_2", "Cw_normalized_alpha_2") #alpha_2 pval=0.009
partial_test("Lw_normalized_beta_1", "Cw_normalized_beta_2") #beta_1
partial_test("Lw_normalized_beta_2", "Cw_normalized_beta_2") #beta_2
partial_test("Lw_normalized_theta", "Cw_normalized_theta") #theta

# Cw x Raven controlling for Lw
partial_test("Cw_normalized_alpha", "Lw_normalized_alpha") #alpha pval=0.056
partial_test("Cw_normalized_alpha_1", "Lw_normalized_alpha_1") #alpha_1 pval=0.030
partial_test("Cw_normalized_alpha_2", "Lw_normalized_alpha_2") #alpha_2 pval=0.034
partial_test("Cw_normalized_beta_1", "Lw_normalized_beta_2") #beta_1
partial_test("Cw_normalized_beta_2", "Lw_normalized_beta_2") #beta_2
partial_test("Cw_normalized_theta", "Lw_normalized_theta") #theta

