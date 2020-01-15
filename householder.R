### series of algorithms to compute qr decomposition and various utilities
### source: https://blogs.mathworks.com/cleve/2016/10/03/householder-reflections-and-the-qr-decomposition/#3c155450-18cf-4e5d-a095-c5c9e12e2e28

A <- matrix(c(12, -51, 4, 6, 167, -68, -4, 24, -41), 3, 3)
x <- c(12, 6, -4)

### to be used with full qr function
A <- matrix(c(35, 1,  6,  26, 19, 24,
              3,  32, 7,  21, 23, 25,
              31, 9,  2,  22, 27, 20,
              8,  28, 33, 17, 10, 15,
              30, 5,  34, 12, 14, 16,
              4,  36, 29, 13, 18, 11), 6, 6, byrow = T)


### function to generate householder reflections
house_ref <- function(x){
  nu <- sqrt(sum(x^2))
  if(nu != 0){
    u <- x/nu
    u[1] <- u[1] + (sign(u[1]) + (u[1] == 0))
    u <- u / sqrt(abs(u[1]))
  } else {
    u <- x
    u[1] <- sqrt(2)
  }
  u
}
house_ref(x)
u <- house_ref(x)

### function to generate householder matrix
### where u is a householder reflection and x is a vector of interest
house_mat <- function(u, x){
  x - u%*%crossprod(u, x)
}
house_mat(u = u, x = x)
house_mat(u = u, x = diag(1, 3))


### computes qr factorization via householder reflections
### stores in compact form a la LINPACK
house_qr <- function(x){
  m <- nrow(x)
  n <- ncol(x)
  k <- min(m, n)
  # U <- matrix(0, m, n)
  qraux <- numeric(length = k)
  R <- x
  for(i in 1:k){
    u <- house_ref(R[i:m, i])
    # U[i:m, i] <- u
    qraux[i] <- u[1]
    R[i:m, i:n] <- house_mat(u, R[i:m, i:n])
    # R[-(1:i), i] <- 0
    R[-(1:i), i] <- u[-1]
  }
  out <- list(qr = R, qraux = qraux)
  out
}
QR <- house_qr(A)
qr(A)


### function which returns R portion of qr decomposition
### accepts as input a qr decomposition
house_qrR <- function(x){
  x <- x$qr
  m <- ncol(x)
  for(j in 1:m){
    x[-(1:j), j] <- 0
  }
  x
}


### function which returns matrix-vector or matrix-matrix product Qy
house_qrqy <- function(x, y){
  r <- x$qr
  qraux <- x$qraux
  m <- nrow(r)
  n <- ncol(r)
  U <- numeric(length = m)
  Z <- y
  for(j in n:1){
    U[j] <- qraux[j]
    U[-(1:j)] <- r[-(1:j), j]
    Z <- house_mat(U, Z)
  }
  Z
}
house_qrqy(QR, diag(1, 6))


### function which returns matrix-vector or matrix-matrix product Qty
house_qrqty <- function(x, y){
  r <- x$qr
  qraux <- x$qraux
  m <- nrow(r)
  n <- ncol(r)
  U <- numeric(length = m)
  Z <- y
  for(j in 1:n){
    U[j] <- qraux[j]
    U[-(1:j)] <- r[-(1:j), j]
    Z <- house_mat(U, Z)
  }
  Z
}
house_qrqy(QR, diag(1, 6))


###function which returns Q portion of qr decomposition
### special case of Qy where y is a diagonal matrix of size == nrow(qr)
house_qrQ <- function(x){
  r <- x$qr
  qraux <- x$qraux
  m <- nrow(r)
  n <- ncol(r)
  U <- numeric(length = m)
  Z <- diag(1, m)
  for(j in n:1){
    U[j] <- qraux[j]
    U[-(1:j)] <- r[-(1:j), j]
    Z <- house_mat(U, Z)
  }
  Z
}



### speed check!
microbenchmark::microbenchmark(
  qr(A), house_qr(A)
)
