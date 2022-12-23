convKernel <- function(sigma = 1.4, k = c('gaussian','LoG','sharpen','laplacian','emboss','sobel')) {
  k <- match.arg(k)
  l <- sigma * 7
  # check if odd and if not increas by one
  if (l%%2==0) l <- l + 1
  # dynamic adaptation of kernel size according the value of sigma
  x <- c(-floor(l/2):floor(l/2))
  y <- c(-floor(l/2):floor(l/2))
  if (k=='gaussian')
    M <- outer(X = x, Y = y, FUN = function(X, Y) return(1/(2*pi*sigma^2)*exp(-(X^2+Y^2)/(2*sigma^2))))
  
  if (k=='LoG')
    M <- outer(X = x, Y = y, FUN = function(X, Y) return(-1/(pi*sigma^4)*(1-(X^2+Y^2)/(2*sigma^2))*exp(-(X^2+Y^2)/(2*sigma^2))))
  
  if (k=='sharpen')   M <- matrix(data = c(0,-1,0,-1,5,-1,0,-1,0), nrow = 3)
  if (k=='laplacian') M <- matrix(data = c(.5,1,.5,1,-6,1,.5,1,.5), nrow = 3)
  if (k=='emboss')    M <- matrix(c(2,0,0,0,-1,0,0,0,-1), nrow = 3)
  if (k=='sobel')     M <- matrix(c(1,2,1,0,0,0,-1,-2,-1), nrow = 3)
  # create S3 class
  if ((k == 'LoG') || (k == 'gaussian'))
    output <- list('matrix' = M, 'kernel' = k, 'sigma' = sigma)
  else
    output <- list('matrix' = M, 'kernel' = k)
  class(output) <- 'convKern'
  return(output)
}

#' @S3method  print convKern
print.convKern <- function(x, ...) {
  cat('\nConvolution Kernel:', x$kernel)
  if ((x$kernel == 'LoG') || (x$kernel == 'gaussian')) cat('\n\nSigma =', x$sigma)
  cat('\n\nConvolution Matrix:\n')
  print(x$matrix)
}


#' @S3method plot convKern
#' @import fields grDevices
plot.convKern <- function(x, ...) {
  col1 <- colorRampPalette(c('black', 'midnightblue', 'darkblue','dodgerblue', 'forestgreen', 'darkolivegreen1', 'gold1', 'orange', 'red'))
  ncolors <- 1000  
  if ((x$kernel == 'gaussian') || (x$kernel == 'LoG')) {
    image.plot(x$matrix, col = col1(ncolors), main = paste('Kernel:', x$kernel, ' - Sigma =', x$sigma))  
  }
  else
    image.plot(x$matrix, col = col1(ncolors), main = paste('Kernel:', x$kernel))
}




applyFilter <- function(x, kernel) {
  # check the k entry
  if (class(kernel)!='convKern') stop('kernel MUST be a convKern class object')
  # prepare matrix if x is matrix by adding required extra lines
  # number of extralines to be added to matrix and arrays (extraplanes)
  if (length(dim(x)) > 3) stop('applyFilter function works only on matrices and 3D arrays')
  extralines <- dim(kernel$matrix)[1] %/% 2
  # now resizes array or matrix adding as new extralines (or planes) as the length
  # of half side of kernel - 1 (the center pixel)
  if (class(x)=='matrix') {
    for (n in 1:extralines) x <- cbind(x[,1], x, x[,ncol(x)])
    for (n in 1:extralines) x <- rbind(x[1,], x, x[nrow(x),])
  }
  if (class(x)=='array') {
    for (n in 1:extralines) x <- abind(x[,1,], x, x[,dim(x)[2],], along = 2)
    for (n in 1:extralines) x <- abind(x[1,,], x, x[dim(x)[1],,], along = 1)
  }
  Nrow <- dim(x)[1]
  Ncol <- dim(x)[2]
  if (class(x)=='matrix') Nslices <- 1 else Nslices <- dim(x)[3]
  dataOutput <- x
  # building the index matrix for transporting the kernel on the array or matrix
  kindex <- c()
  for (n in -extralines:extralines)  # index for columns
    for (m in -extralines:extralines)  # index for rows
      kindex <- c(kindex, n*Nrow + m)
  # apply laplacian and emboss convolution kernels
  if ((kernel$kernel == 'laplacian') || (kernel$kernel == 'emboss')) {
    result <- .C('applyKernelWithoutNorm', as.double(x), as.double(kernel$matrix), as.integer(extralines), as.integer(kindex),
                 as.integer(Nrow), as.integer(Ncol), as.integer(Nslices), as.double(dataOutput))
    result <- result[[8]]
  }
  # apply Gaussian, LoG and sharpen convolution kernels
  else if (kernel$kernel != 'sobel')  {# apply Gaussian, LoG and sharpen convolution kernels
    result <- .C('applyKernel', as.double(x), as.double(kernel$matrix), as.integer(extralines), as.integer(kindex),
                 as.integer(Nrow), as.integer(Ncol), as.integer(Nslices), as.double(dataOutput))
    result <- result[[8]]
  }
  # apply Sobel convolution kernels
  if (kernel$kernel == 'sobel') {
    Gx <- .C('applyKernelWithoutNorm', as.double(x), as.double(kernel$matrix), as.integer(extralines), as.integer(kindex),
             as.integer(Nrow), as.integer(Ncol), as.integer(Nslices), as.double(dataOutput))
    Gx <- Gx[[8]]
    dataOutput <- x
    # Gy convolution kernel is given by transposition of sobel default convolution kernel 
    Gy <- .C('applyKernelWithoutNorm', as.double(x), as.double(t(kernel$matrix)), as.integer(extralines), as.integer(kindex),
             as.integer(Nrow), as.integer(Ncol), as.integer(Nslices), as.double(dataOutput))
    Gy <- Gy[[8]]
    result <- sqrt(Gx^2 + Gy^2)
  }
  
  
  if (class(x)=='matrix') {
    output <- matrix(data = result, nrow = Nrow)
    # resize matrix to original size
    output <- output[(extralines+1):(nrow(output) - extralines), (extralines+1):(ncol(output) - extralines)]
  }
  if (class(x)=='array')  {
    output <- array(data = result, dim = c(Nrow, Ncol, Nslices))
    output <- output[(extralines+1):(nrow(output) - extralines), (extralines+1):(ncol(output) - extralines),]
  }
  return(output)
}