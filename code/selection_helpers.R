## functions for handling on-the-fly rescaling

imr <- function (x) {dnorm(-1*x)/(1-pnorm(-1*x))}

datamean <- function (x,y,avx,avy) {
  
  (ifelse(is.na(x),0,x)*avx + ifelse(is.na(y),0,y)*avy)/(avx+avy)
  
}

is.wholenumber <- function (x) {x==round(x)}


qformula <- function (a,b,c) {
  c((-1*b + sqrt(b^2 - 4*a*c))/(2*a), (-1*b + sqrt(b^2 - 4*a*c))/(2*a))
}


arrange_uppertri <- function (x) {
  
  if (!is.wholenumber(qformula(1,1,-2*length(x))[1]))
    stop(paste0("Length of matrixvars needs to be coercable into lower triangle (ie 3, 6, 10, etc.),\n  Length is ",length(x)))
  
  if (is.wholenumber(qformula(1,1,-2*length(x))[1])) {
    mat <- matrix(0, nrow = qformula(1,1,-2*length(x))[1], ncol = qformula(1,1,-2*length(x))[1])
    mat[upper.tri(mat, diag = TRUE)] <- x
  }
  
  mat
}

un_cholesky <- function (x) {
  t(arrange_uppertri(x)) %*% arrange_uppertri(x)
}

make_rescaler_str <- function (matrixvars,datavars) {
  
  mat1 <- t(arrange_uppertri(matrixvars))
  mat2 <- arrange_uppertri(matrixvars)
  
  matdims <- qformula(1,1,-2*length(matrixvars))
  
  if (length(datavars)!=matdims[1])
    stop(paste0("Length of datavars is ",length(datavars),"; needs length ",matdims[1]," to be conformable with matrixvars"))
  
  matrixelements <- matrix(data=NA,nrow=nrow(mat1),ncol=ncol(mat2))
  for (i in 1:nrow(mat1)) {
    for (j in 1:ncol(mat2)) {
      matrixelements[i,j] <- paste0(mat1[i,],"*",mat2[,j],collapse=" + ")
    }
  }
  matrixelements
  
  matrixelements <- matrix(paste0("(",matrixelements,")"),nrow=nrow(mat1))
  
  matrixelements <- matrixelements[-1,-1]
  datavars <- datavars[-1]
  
  
  multmat <- c()
  
  for (k in 1:nrow(matrixelements)) {
    multmat <- c(multmat,
                 paste0(datavars[k],"*","(",paste0(datavars,"*",matrixelements[k,],collapse=" + "),")")
    )
  }
  
  multmat <- paste0(multmat,collapse=" + ")
  
  multmat
}

getcovparams <- function (num_params=3) {
  x <- num_params
  one <- expand.grid(letters[1:x],1:x)
  two <- paste0(one[,1],one[,2]) %>% sort()
  three <- matrix(two,nrow=x,ncol=x,byrow = T)
  out <- three[lower.tri(three,diag=T)] %>% sort()
  out
}

build_rescale_function <- function (corvars) {
  
  matrixvars <- getcovparams(length(corvars))
  
  if (sum(duplicated(c(matrixvars,corvars))) > 0)
    stop("Duplicate names")
  
  body1 <- make_rescaler_str(matrixvars,corvars)
  
  args1 <- paste0(matrixvars,collapse=", ")
  args2 <- paste0(corvars,collapse=", ")
  args <- paste0(args1,", ",args2)
  
  doit <- paste0("rescaler <- function (",args,") ","{sqrt(",body1,"+ (pi^2)/3)}")
  
  eval(str2lang(doit),envir = .GlobalEnv)
  
}

startparjitter <- function (x) {
  
  l <-  x[which(names(x)!=c("scale_RP","scale_SP"))] %>% length()
  jitter <- runif(l,min=0.95,max=1.05)
  
  if (sum(c("scale_RP","scale_SP") %in% names(x))) {
    jitter <- c(jitter,1,1)
    
  } 
  
  y <- x*jitter
  y
}

