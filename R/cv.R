library(dplyr)

#' Create dataframe 'dfX' with selected features
#'
#' @param ... names of selected features columns
#'
#' @return dataframe dfX
#' @export
#'
#' @examples
new_dataframe_x <- function(df, ...) {
  cols <- c(...) # Convert the input to a character vector
  pattern <- paste(cols, collapse = "|") # Create a regular expression pattern
  new_df <- select(df, matches(pattern))
  return(new_df)
}

#' Create dataframe 'dfY' with Y values
#'
#' @param new_vector vector with y values
#'
#' @return dataframe dfY
#' @export
#'
#' @examples
new_dataframe_y <- function(new_vector) {
  dfY=data.frame()
  length(dfY)=length(new_vector)
  dfY<- cbind( new_vector)
  colnames(dfY)[1] <- "y"
  dfY=data.frame(dfY)
  return(dfY)
}

#' Construct X matrix
#'
#' @param dfX dataframe with selected X features only. Can be created with new_dataframe_x().
#'
#' @return X matrix, first column contains 1s, other rows contain the X features in dfX
#' @export
#'
#' @examples
constructX<-function(dfX){
  X<-dfX
  X<-cbind(1,X)
  colnames(X)[1] <- "x0"
  X<-data.matrix(X)
  X<-unname(X)
  return(X)
}

#' Construct Y matrix
#'
#' @param dfY dataframe with 'y' column. Can be created with new_dataframe_y().
#'
#' @return Y values in vector
#' @export
#'
#' @examples
constructY<-function(dfY){
  Y=c()
  Y=dfY[["y"]]
  return(Y)
}

#' Construct X-transpose X
#'
#' @param X X matrix constructed with constructX()
#'
#' @return transpose of X multiplied by X
#' @export
#'
#' @examples
constructXT_X<-function(X){
  a= t(X) %*% X
  return(a)
}

#' Construct The Inverse of X-transpose X
#'
#' @param X X matrix constructed with constructX()
#' @param XT_X (optional) x-tranposed multiplied by x
#'
#' @return the inverse of the transpose of X multiplied by X
#' @export
#'
#' @examples 
constructinverseXT_X<-function(X, XT_X=NULL){
  
  if(is.null(XT_X)){
    XT_X=constructXT_X(X)
  }
  
  b=solve(XT_X)
  return(b)
}

#' Construct X-transpose Y
#'
#' @param X X matrix constructed with constructX()
#' @param Y Y vector constructed with constructY()
#'
#' @return transpose of X multiplied by Y
#' @export
#'
#' @examples
constructXT_Y<-function(X, Y){
  
  a= t(X) %*% Y
  return(a)
  
}

#' Construct Y-transpose Y
#'
#' @param Y Y vector constructed with constructY()
#'
#' @return transpose of Y multiplied by Y
#' @export
#'
#' @examples
constructYT_Y<-function(Y){
  
  a= t(Y) %*% Y
  return(a)
  
}

#' Estimate Regression Coefficients B
#'
#' @param X X matrix constructed with constructX()
#' @param Y Y vector constructed with constructY()
#' @param Inv_XT_X (optional) inverse of x-tranposed multiplied by x
#' @param XT_Y (optional) x-tranposed multiplied by y
#'
#' @return B vector
#' @export
#'
#' @examples
computeB<-function(X,Y, Inv_XT_X=NULL,XT_Y=NULL){
  
  if(is.null(Inv_XT_X)){
    Inv_XT_X=constructinverseXT_X(X)
  }
  if(is.null(XT_Y)){
    XT_Y=constructXT_Y(X,Y)
  }
  
  b= Inv_XT_X %*% XT_Y
  return(b)
  
}

#' Compute Total Sum of Squares
#'
#' @param X X matrix constructed with constructX()
#' @param Y Y vector constructed with constructY()
#' @param YT_Y transpose of Y multiplied by Y
#'
#' @return total sum of squares
#' @export
#'
#' @examples
computeSST<-function(X,Y,YT_Y=NULL){
  
  if(is.null(YT_Y)){
    YT_Y=constructYT_Y(Y)
  }
  n=length(Y)
  ybar=mean(Y)
  sst=YT_Y-n*ybar*ybar
  return(sst)
  
}

#' Compute Residual Sum of Squares
#'
#' @param X X matrix constructed with constructX()
#' @param Y Y vector constructed with constructY()
#' @param YT_Y transpose of Y multiplied by Y
#' @param XT_Y (optional) x-tranposed multiplied by y
#' @param B (optional) vector of estimated regression coefficients. [The remaining parameters are not needed if B is provided]
#' @param Inv_XT_X (optional) inverse of x-tranposed multiplied by x
#'
#' @return residual sum of squares
#' @export
#'
#' @examples
computeSSE<-function(X,Y,YT_Y=NULL,B=NULL,Inv_XT_X=NULL,XT_Y=NULL){
  
  if(is.null(YT_Y)){
    YT_Y=constructYT_Y(Y)
  }
  if(is.null(XT_Y)){
    XT_Y=constructXT_Y(X,Y)
  }
  if(is.null(B)){
    B=computeB(X,Y,Inv_XT_X=Inv_XT_X,XT_Y=XT_Y)
  }
  sse=YT_Y-t(B)%*%XT_Y
  return(sse)
}

#' Compute Regression Sum of Squares
#'
#' @param X X matrix constructed with constructX()
#' @param Y Y vector constructed with constructY()
#' @param B (optional) vector of estimated regression coefficients. [The remaining parameters are not needed if B is provided]
#' @param Inv_XT_X (optional) inverse of x-tranposed multiplied by x
#' @param XT_Y (optional) x-tranposed multiplied by y
#' 
#' @return regression sum of squares
#' @export
#'
#' @examples
computeSSR<-function(X,Y,B=NULL,Inv_XT_X=NULL,XT_Y=NULL){
  
  if(is.null(B)){
    B=computeB(X,Y,Inv_XT_X=Inv_XT_X,XT_Y=XT_Y)
  }
  if(is.null(XT_Y)){
    XT_Y=constructXT_Y(X,Y)
  }
  n=length(Y)
  ybar=mean(Y)
  ssr=t(B)%*%XT_Y-n*ybar*ybar
  return(ssr)
}

#' Compute Regression degree of freedom
#'
#' @param X X matrix constructed with constructX()
#'
#' @return regression degree of freedom
#' @export
#'
#' @examples
calculateDFR<-function(X){
  P=ncol(X)
  K= P-1
  return(K)
}
#' Compute Residual degree of freedom
#'
#' @param X X matrix constructed with constructX()
#'
#' @return residual degree of freedom
#' @export
#'
#' @examples
calculateDFE<-function(X){
  P=ncol(X)
  K= P-1
  N=nrow(X)
  dfe=N-P
  return(dfe)
}
#' Compute Total degree of freedom
#'
#' @param X X matrix constructed with constructX()
#'
#' @return total degree of freedom
#' @export
#'
#' @examples
calculateDFT<-function(X){
  N=nrow(X)
  dft=N-1
  return(dft)
}

#' Compute the Mean Square Error
#'
#' @param X X matrix constructed with constructX()
#' @param Y Y vector constructed with constructY()
#' @param sse (optional) residual sum of squares
#'
#' @return mean square error
#' @export
#'
#' @examples
calculateMSE<-function(X,Y,sse=NULL){
  if(is.null(sse)){
    sse=computeSSE(X,Y)
  }
  dfe=calculateDFE(X)
  MSE=sse/dfe
  return(MSE)
}

#' Compute the Mean Square Regression
#'
#' @param X X matrix constructed with constructX()
#' @param Y Y vector constructed with constructY()
#' @param ssr (optional) regression sum of squares
#'
#' @return mean square regression
#' @export
#'
#' @examples
calculateMSR<-function(X,Y,ssr=NULL){
  if(is.null(ssr)){
    ssr=computeSSR(X,Y)
  }
  dfr=calculateDFR(X)
  MSR=ssr/dfr
  
  return(MSR)
}

#' Calculate F Value
#'
#' @param X X matrix constructed with constructX()
#' @param Y Y vector constructed with constructY()
#' @param mse (optional) mean square error
#' @param msr (optional) mean square regression
#'
#' @return anova f value
#' @export
#'
#' @examples
calculateFVALUE<-function(X,Y,msr=NULL,mse=NULL){
  if(is.null(msr)){
    msr=calculateMSR(X,Y)
  }
  if(is.null(mse)){
    mse=calculateMSE(X,Y)
  }
  FVALUE=msr/mse
  
  return(FVALUE)
}

#' Calculate Anova Table
#'
#' @param X X matrix constructed with constructX()
#' @param Y Y vector constructed with constructY()
#' @param YT_Y transpose of Y multiplied by Y
#' @param B (optional) vector of estimated regression coefficients. [The remaining parameters are not needed if B is provided]
#' @param Inv_XT_X (optional) inverse of x-tranposed multiplied by x
#' @param XT_Y (optional) x-tranposed multiplied by y
#'
#' @return anova table
#' @export
#'
#' @examples
calculateANOVATABLE<-function(X,Y,YT_Y=NULL,B=NULL,Inv_XT_X=NULL,XT_Y=NULL){
  SST = computeSST(X,Y,YT_Y=YT_Y)
  SSR=computeSSR(X,Y,B,Inv_XT_X=Inv_XT_X,XT_Y=XT_Y)
  SSE=computeSSE(X,Y,B,Inv_XT_X=Inv_XT_X,XT_Y=XT_Y,YT_Y=YT_Y)
  DF_SST=calculateDFT(X)
  DF_SSR=calculateDFR(X)
  DF_SSE=calculateDFE(X)
  MSE=calculateMSE(X,Y,sse=SSE)
  MSR=calculateMSR(X,Y,ssr=SSR)
  FVALUE=calculateFVALUE(X,Y)
  
  tab <- matrix(c(SSR, DF_SSR, MSR, FVALUE, SSE, DF_SSE, MSE,' ', SST,DF_SST,' ',' '), ncol=4, byrow=TRUE)
  colnames(tab) <- c('SS','D.O.F','MS','F')
  rownames(tab) <- c('regression','residual','total')
  tab <- as.table(tab)
  
  return(tab)
}
#' Get a specific B value
#'
#' @param X X matrix constructed with constructX()
#' @param Y Y vector constructed with constructY()
#' @param n index of B to get
#' @param B (optional) vector of estimated regression coefficients. 
#'
#' @return specific B value
#' @export
#'
#' @examples
getBvalues<-function(X,Y,n,B=NULL){
  if(is.null(B)){
    B=computeB(X,Y)
  }
  r = 0
  for (i in 1:length(B) ){
    if(i==n+1){
      r=B[i]
    }
  }
  return(r)
}

#' Get specific Y value
#'
#' @param X X matrix constructed with constructX()
#' @param Y Y vector constructed with constructY()
#' @param n index of Y to get
#' @param y_pred (optional) y-predicted, calculated with compute_y_pred()
#'
#' @return specific Y value
#' @export
#'
#' @examples
get_ypred_values<-function(X,Y,n,y_pred =NULL){
  if(is.null(y_pred)){
    y_pred=compute_y_pred(X,Y)
  }
  for (i in 1:length(y_pred)){
    if(i==n){
      r=y_pred[i]
    }
  }
  return(r)
}

#' Compute Predicted Y values
#'
#' @param X X matrix constructed with constructX()
#' @param Y Y vector constructed with constructY()
#' @param B (optional) vector of estimated regression coefficients. 
#'
#' @return y-predicted
#' @export
#'
#' @examples
compute_y_pred<-function(X,Y,B=NULL){
  if(is.null(B)){
    B=computeB(X,Y)
  }
  y_pred=X %*% B
  return(y_pred)
}

#' calculateCI_B
#'
#'@description This function calculates the confidence interval of any Beta estimator in a multiple regression model
#' @param X X matrix constructed with constructX()
#' @param Y Y vector constructed with constructY()
#' @param index the index of Beta 
#' @param alpha (optional) the value of significance level
#' @param known (optional) a boolean parameter that takes (T/F) as value
#' @param sigma (optional) the standard error of the dataset and it has  default value (NULL)
#' @param mse (optional) mean square error
#' @param inverseXT_X (optional) inverse of x-tranposed multiplied by x
#' @param B (optional) vector of estimated regression coefficients. 
#' 
#' @return the confidence interval of Beta (string)
#' @export
#'
calculateCI_B=function(X,Y,index,alpha = .05,known = F,sigma= NULL,mse=NULL,inverseXT_X=NULL,B=NULL){
  dof=calculateDFE(X)
  if(is.null(mse)){
    mse=calculateMSE(X,Y)
  }
  if(is.null(inverseXT_X)){
    inverseXT_X=constructinverseXT_X(X)
  }
  b=getBvalues(X,Y,index,B=B)
  l=index+1
  if(known){
    z_score <- qnorm(1 - alpha/2)
    Upper=b+(z_score*(((sigma)^2*inverseXT_X[l,l])^.5))
    Lower=b-(z_score*(((sigma)^2*inverseXT_X[l,l])^.5))
    return(paste("The", (1-alpha) * 100, "% confidence interval for beta ", index, " is [", Lower, ", ", Upper, "]."))
  }
  else{
    t_score <- qt(1 - alpha/2, dof)
    Upper=b+(t_score*(mse*inverseXT_X[l,l])^.5)
    Lower=b-(t_score*(mse*inverseXT_X[l,l])^.5)
    return(paste("The", (1-alpha) * 100, "% confidence interval for beta ", index, " is [", Lower, ", ", Upper, "]."))  }
}


#' Get the mean Y value at certain X values
#'
#' @param X X matrix constructed with constructX()
#' @param Y Y vector constructed with constructY()
#' @param X0 X values
#' @param B (optional) vector of estimated regression coefficients. 
#' 
#' @return y0
#' @export
#'
#' @examples
compute_y0<-function(X,Y, X0,B=NULL){
  if(is.null(B)){
    B=computeB(X,Y)
  }
  y_pred= t(X0)%*% B
  return(y_pred)
}

#' Confidence Interval for the Mean Response
#'
#' @param X X matrix constructed with constructX()
#' @param Y Y vector constructed with constructY()
#' @param alpha (optional) the value of significance level
#' @param sigma (optional) the standard error of the dataset and it has  default value (NULL)
#' @param mse (optional) mean square error
#' @param known (optional) a boolean parameter that takes (T/F) as value
#' @param B (optional) vector of estimated regression coefficients. 
#' @param inverse_XTX (optional) inverse of x-tranposed multiplied by x
#'
#' @return mean response confidence interval (string)
#' @export
#'
#' @examples
computeCI_mean_response <- function(X,Y,X0,alpha=0.05, sigma = NULL, mse=NULL, known = F,B=NULL,inverse_XTX=NULL){
  
  y_hat = compute_y0(X,Y,X0=X0,B=B)
  dof = calculateDFE(X)
  if(is.null(mse)){
    mse=calculateMSE(X,Y)
  }
  if(is.null(inverse_XTX)){
    inverse_XTX=constructinverseXT_X(X)
  }
  
  if (known){
    z_score <- qnorm(1 - alpha/2)
    Lower = y_hat - (z_score*(mse*((t(X0) %*% inverse_XTX) %*% X0))^0.5)
    Upper =  y_hat + (z_score*(mse*((t(X0) %*% inverse_XTX) %*% X0))^0.5)
    return(paste("The", (1-alpha) * 100, "% confidence interval for mean response at ",'[',paste(X0,collapse=', '),']', " is [", Lower, ",", Upper, "]."))
    
  }else{
    t_score = qt(1-alpha/2, dof)
    Lower = y_hat - (t_score*(mse*((t(X0) %*% inverse_XTX) %*% X0))^0.5)
    Upper =  y_hat + (t_score*(mse*((t(X0) %*% inverse_XTX) %*% X0))^0.5)
    return(paste("The", (1-alpha) * 100, "% confidence interval for mean response at ",'[',paste(X0,collapse=', '),']', " is [", Lower, ",", Upper, "]."))
  }
}

#' Confidence Interval for a New Observation
#'
#' @param X X matrix constructed with constructX()
#' @param Y Y vector constructed with constructY()
#' @param alpha (optional) the value of significance level
#' @param sigma (optional) the standard error of the dataset and it has  default value (NULL)
#' @param mse (optional) mean square error
#' @param known (optional) a boolean parameter that takes (T/F) as value
#' @param B (optional) vector of estimated regression coefficients. 
#' @param inverse_XTX (optional) inverse of x-tranposed multiplied by x
#' 
#' @return new observation confidence interval (string)
#' @export
#'
#' @examples
computeCI_new_observation <- function(X,Y,X0,alpha=0.05, sigma = NULL, mse=NULL, known = F,B=NULL,inverse_XTX=NULL){
  
  y_hat = compute_y0(X,Y,X0=X0,B=B)
  dof = calculateDFE(X)
  if(is.null(mse)){
    mse=calculateMSE(X,Y)
  }
  if(is.null(inverse_XTX)){
    inverse_XTX=constructinverseXT_X(X)
  }
  
  if (known){
    z_score = qnorm(1 - alpha/2)
    Lower = y_hat - (z_score*(mse*(1+((t(X0)%*%inverse_XTX)%*%X0)))^0.5)
    Upper =  y_hat + (z_score*(mse*(1+((t(X0)%*%inverse_XTX)%*%X0)))^0.5)
    return(paste('The', (1-alpha) * 100, '% confidence interval for mean response at' ,'[',paste(X0,collapse=', '),']', ' is [', Lower, ',', Upper,' ].'))
  }
  else{
    t_score = qt(1-alpha/2, dof)
    Lower = y_hat - (t_score*(mse*(1+((t(X0) %*% inverse_XTX) %*% X0)))^0.5)
    Upper =  y_hat + (t_score*(mse*(1+((t(X0) %*% inverse_XTX)%*%X0)))^0.5)
    return(paste('The', (1-alpha) * 100, '% confidence interval for new observation at' ,'[',paste(X0,collapse=', '),']', ' is [', Lower, ',', Upper,' ].'))
  }
}


#' Compute the Regular Error
#'
#' @param X X matrix constructed with constructX()
#' @param Y Y vector constructed with constructY()
#' @param B (optional) vector of estimated regression coefficients. 
#'
#' @return regular error vector
#' @export
#'
#' @examples
regular_error <- function(X,Y,B=NULL) {
  if(is.null(B)){
    B=computeB(X,Y)
  }
  Y_Hat = compute_y_pred(X=X,Y=Y,B=B)
  return(Y-Y_Hat)
}

#' Compute the Standard Error
#'
#' @param X X matrix constructed with constructX()
#' @param Y Y vector constructed with constructY()
#' @param E (optional) regular error 
#' @param B (optional) vector of estimated regression coefficients
#' @param MSE (optional) mean square error
#'
#' @return standard error matrix
#' @export
#'
#' @examples
standard_error <- function(X,Y,E=NULL,B=NULL,MSE=NULL) {
  if(is.null(E)){
    E = regular_error(X,Y,B)
  }
  if(is.null(MSE)){
    MSE = calculateMSE(X,Y)
  }
  
  D = c()
  for(e_i in E) {
    d_i = e_i/(MSE^0.5)
    D = c(D,d_i)
  }
  return(as.matrix(D))  
  
}

#' Compute the Student Error
#'
#' @param X X matrix constructed with constructX()
#' @param Y Y vector constructed with constructY()
#' @param E (optional) regular error 
#' @param B (optional) vector of estimated regression coefficients
#' @param MSE (optional) mean square error
#' @param inv_XT_X (optional) inverse of x-tranposed multiplied by x
#'
#' @return student error matrix
#' @export
#'
#' @examples
student_error <- function(X,Y,E=NULL,B=NULL,MSE=NULL,inv_XT_X=NULL) {
  if(is.null(E)){
    E = regular_error(X,Y,B)
  }
  if(is.null(MSE)){
    MSE = calculateMSE(X,Y)
  }
  if(is.null(inv_XT_X)){
    inv_XT_X = constructinverseXT_X(X)
  }
  
  H = X %*% inv_XT_X %*% t(X)
  R = c()
  for(i in 1:length(E)) {
    e_i = E[i]
    h_ii = H[i,i]
    r_i = e_i/((MSE*(1-h_ii))^0.5)
    R = c(R,r_i)
  }
  
  return(as.matrix(R))
}


#' Calculate R-Squared value
#'
#' @param X X matrix constructed with constructX()
#' @param Y Y vector constructed with constructY()
#' @param SSR regression sum of squares
#' @param SST total sum of squares
#'
#' @return r-squared
#' @export
#'
#' @examples
calculateRSQUARED<-function(X,Y,SSR=NULL,SST=NULL){
  if(is.null(SSR)){
    SSR=computeSSR(X,Y)
  }
  if(is.null(SST)){
    SST=computeSST(X,Y)
  }
  RSQUARED=SSR/SST
  return(RSQUARED)
}


#' Calculate the Adjusted R-Squared value
#'
#' @param X X matrix constructed with constructX()
#' @param Y Y vector constructed with constructY()
#' @param RSQUARED rsquared, calculated with calculateRSQUARED()
#'
#' @return
#' @export
#'
#' @examples
calculateRSQUAREDADJUSTED<-function(X,Y,RSQUARED=NULL){
  if(is.null(RSQUARED)){
    RSQUARED= calculateRSQUARED(X,Y)
  }
  dft= calculateDFT(X)
  dfe= calculateDFE(X)
  ADJUSTED_RSQUARED <- 1 - ((1 - RSQUARED) * ((dft) / (dfe)))
  
  return(ADJUSTED_RSQUARED)
}


#' Normalize X
#'
#' @param X X matrix constructed with constructX()
#'
#' @return normalization of x
#' @export
#'
#' @examples
normalization <- function(X) {
  
  X <- X[, -1]
  mean_cols <- colMeans(X)
  xstarsxx_vec <- numeric(ncol(X))  # initialize an empty vector for xstarsxx
  
  for (j in 1:ncol(X)) {
    xstarsxx <- sum(X[, j]^2) - (nrow(X) * (mean_cols[j]^2))
    xstarsxx_vec[j] <- xstarsxx  # store xstarsxx in the vector
    for (i in 1:nrow(X)) {
      xstar <- X[i, j]
      xstarmean <- mean_cols[j]
      X[i, j] <- (xstar -mean_cols[j]) / sqrt(xstarsxx)
    }
  }
  
  return(X)
}

#' Calculate the correlation matrix
#'
#' @param X X matrix constructed with constructX()
#'
#' @return correlation matrix
#' @export
#'
#' @examples
coerelation_matrix<-function(X){
  X=normalization(X)
  X_T=t(X)
  return(X_T %*% X)
}
#' Calculate the variance inflation factor
#'
#' @param X X matrix constructed with constructX()
#'
#' @return vif
#' @export
#'
#' @examples
VIF <- function(X) {
  Q <- coerelation_matrix(X)
  Q <- solve(Q)
  diag_X <- diag(Q)
  diag_names <- colnames(Q)[which(diag(Q) != 0)]
  M <- matrix(nrow = 1, ncol = length(diag_X))
  colnames(M) <- diag_names
  M[1, ] <- diag_X
  return(M)
}
#' Calculate the partial F value
#'
#' @param X X matrix constructed with constructX()
#' @param Y Y vector constructed with constructY()
#' @param idx_vec vector with the indices of reduced model (non-tested features)
#' @param SSRT total regression sum of squares
#' @param MSEM model's mean square error
#'
#' @return partial f-value
#' @export
#'
#' @examples
partialF<-function(X,Y,idx_vec,SSRT=NULL,MSEM=NULL){
  if(is.null(SSRT)){
    SSRT=computeSSR(X,Y)
  }
  SSRR=computeSSR(X[,idx_vec+1],Y)
  if(is.null(MSEM)){
    MSEM=calculateMSE(X,Y)
  }
  dfregressor=(ncol(X))-(length(idx_vec))
  F=((SSRT-SSRR)/dfregressor)/MSEM
  return(F)
}
