# Intro.... ----
# Código para método de la bisección 
# Estudiante: Ayerim Cinai Montelongo Cano
# 30 de Septiembre de 2021

# Indicaciones: Los programas no van a pedirle al usuario que introduzca ninguna de las funciones, ni tampoco los valores que se necesiten. Todo lo que requieras para que corra o se ejecute el programa ya debe de estar incluido, simplemente vas a comentar y descomentar líneas de codigo.


options(digits=5, nsmall = 5)

raiz <- function(a,b,TOL = 1E-10, N = 100) {
  
  if(f(a)*f(b) > 0)
    stop('\nLa funcion no tiene una raiz en el intervalo [a,b].\n')
  else
    if(abs(f(a)) < TOL)
      cat('\nLa funcion tiene una raiz en a.\n')
  else
    if(abs(f(b)) < TOL)
      cat('\nLa funcion tiene una raiz en b.\n')
  else
    cat('\nLa funcion tiene una raiz en el intervalo (a,b).\n')
  
  n = 1
  
  an = a
  bn = b
  
  # calculando el punto medio
  
  Pn = (bn+an)/2
  
  df[[n]] <- data.frame(n, an, bn, Pn, fp = f(Pn))
  
  # return(df)
  
  while(abs(f(Pn)) > TOL & n < N) {
    
    an_1 = an
    pn_1 = Pn
    bn_1 = bn
    
    n = n + 1
    
    if(f(an_1)*f(pn_1) < 0) {
      an = an_1
      bn = pn_1
    } else {
      an = pn_1
      bn = bn_1
    }
    
    
    Pn = (bn+an)/2
    
    df[[n]] <- data.frame(n, an, bn, Pn, fp = f(Pn))
    
  }
  
  df <- do.call(rbind, df)
  
  Pn <- df[nrow(df), 4]
  
  cat('\nLa funcion tiene una raiz en', Pn)
  
  return(df)
}


# Metodo Biseccion ----
# Mediante el metodo de la biseccion, encuentra la raiz de:

# Ejercicio 1 ----

pi = 3.141592654

TOL = 1E-10

f <- function(x) {tan(pi*x) - 6}

df <- list() # variable de almacenamiento de tabla


tbl <- raiz(a = 0, b = 0.48, TOL = 1E-10)

tbl %>% view()

# Ejercicio 2 ----

f <- function(x) {cos(x) - x}

raiz(a = 0.5, b = pi/4, TOL = 1E-10) -> tbl2

tbl2 %>% view()


# Ejercicio 3 ----

f <- function(x) {sin(sqrt(x)) - x}

raiz(a = 0.4, b = 0.6, TOL = 1E-10, N = 50)

