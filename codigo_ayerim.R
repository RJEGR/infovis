
options(digits = 9)

# 1. Metodo de la definicion
# funciones
Polbase=function(n,s,t){
   L=0 # hay que ponerlo igual a 0 para hacerlo vector, aunque se use en un productorio
   for(i in 1:n){
     L[i]=1 #ahora si se pone a 1
    
       for(j in 1:n){
        
           if (i!=j){ # !=i significa que j no tomará los valores de i
            
              
               L[i]= L[i]*( t - s [j])/ (s[i]- s[j])
               }
         }
    
       }
   return(L) #Return es el resultado
}

# funcion para determinar polinomio
Polinterp=function(B,L,n){
  p=0
  for (i in 1:n){
    p=p + B[i]*L[i]
  }
  
  # return(p)
  
  res <- list('Approximated value' = p) # , 'Neville iterations table'=q
  
  return(res)
  
  
}


# B=c(10,20,35.33,35.5)
# s=c(1,3.5,5,10)



s <- c(0.0, 0.1, 0.3, 0.6, 1.0, 1.1)
B <- c(-6.00000, -5.89483, -5.65014, -5.17788, -4.28172, -3.99583) # El vector B, que contendrá los valores de la función interpolada


n=length(s)

# Punto 1

T1 = 0.41 #valor a interpolar

L1=Polbase(n,s,T1) # El vector L, que contendrá los polinomios de base evaluados en t

pol1=Polinterp(B,L1,n)

pol1

# 2. Metodo Neville
# https://rstudio-pubs-static.s3.amazonaws.com/286755_5d9464e5f87e4579b5afff23276f0381.html


poly.neville <- function(x, y, x0) {
  
  n <- length(x)
  q <- matrix(data = 0, n, n)
  q[,1] <- y
  
  for (i in 2:n) {
    for (j in i:n) {
      q[j,i] <- ((x0 - x[j-i+1]) * q[j,i-1] - (x0 - x[j]) * q[j-1,i-1]) / (x[j] - x[j-i+1])
    }
  }
  
  res <- list('Approximated value'=q[n,n], 'Neville iterations table' = round(q, 5)) # , 
  return(res)
}


poly.neville(s, B, 0.41)


# de newton

# devtools::install_github("FedericoComoglio/rSymPy")

divided.differences <- function(x, y, x0) {
  require(rSymPy)
  n <- length(x)
  q <- matrix(data = 0, n, n)
  q[,1] <- y
  f <- as.character(round(q[1,1], 5))
  fi <- ''
  
  for (i in 2:n) {
    for (j in i:n) {
      q[j,i] <- (q[j,i-1] - q[j-1,i-1]) / (x[j] - x[j-i+1])
    }
    fi <- paste(fi, '*(x - ', x[i-1], ')', sep = '', collapse = '')
    
    f <- paste(f, ' + ', round(q[i,i], 5), fi, sep = '', collapse = '')
  }
  
  x <- Var('x')
  sympy(paste('e = ', f, collapse = '', sep = ''))
  approx <- sympy(paste('e.subs(x, ', as.character(x0), ')', sep = '', collapse = ''))
  
  #     # 'Approximation from Interpolation'=as.numeric(approx), 
  return(list('Approximation from Interpolation'=as.numeric(approx), 
    'Interpolated Function'=f, 'Divided Differences Table'=q))
}

z1 <- Polinterp(B,L1,n)
z2 <- poly.neville(s, B, 0.41)
z3 <- divided.differences(s, B, 0.41)

y <- c(z1$`Approximated value`, z2$`Approximated value`, z3$`Approximation from Interpolation`)

dff <- data.frame(x  = 0.41, y, color = c('a', 'b', 'c'))

# 

# install.packages("ggplot2")

library(ggplot2)

df <- data.frame(y = B, x = s, color = 'black')

ggplot() + 
  geom_point(data = df, aes(x, y)) +
  geom_point(data = dff, aes(x, y, color = color)) +
  xlim(-1, 2)
