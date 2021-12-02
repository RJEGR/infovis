# Intro ----
# Código para método de la bisección 
# Estudiante: Ayerim Cinai Montelongo Cano
# 30 de Septiembre de 2021

# Indicaciones: Los programas no van a pedirle al usuario que introduzca ninguna de las funciones, ni tampoco los valores que se necesiten. Todo lo que requieras para que corra o se ejecute el programa ya debe de estar incluido, simplemente vas a comentar y descomentar líneas de codigo.


# options(digits=8)

options(digits=5)

raiz <- function(a,b,TOL = 1E-10, N = 100) {
  
  df <- list() # variable de almacenamiento de tabla
  
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
  
  df[[n]] <- data.frame(n, an, bn, Pn, fp = f(Pn)) #almacenar en tabla
  
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
  
  df <- do.call(rbind, df) # rbind combina filas y columnas de una lista
  
  Pn <- df[nrow(df), 4] # Ir al último valor de la fila y devolver el valor de la columna 4
  
  cat('\nLa funcion tiene una raiz en: ', Pn, '\n')
  
  return(df)
}

# pi = 3.141592654
pi = 3.1415926535897932

TOL = 1E-10


# Método Bisección ----
# Mediante el método de la bisección, encuentra la raíz de:

# Ejercicio 1 ----

f <- function(x) {tan(pi*x) - 6}

raiz(a = 0, b = 0.48, TOL = 1E-10) -> tbl1

# Ejercicio 2 ----

f <- function(x) {cos(x) - x}

raiz(a = 0.5, b = pi/4, TOL = 1E-10) -> tbl2

# Ejercicio 3 ----

f <- function(x) {sin(sqrt(x)) - x}

raiz(a = 0.4, b = 0.6, TOL = 1E-10) -> tbl3



# plots ---

library(tidyverse)

tbl <- tbl3

caption <- paste0('Número de iteraciones: ', 
                  nrow(tbl))

tbl %>%
  mutate(col = ifelse(n == nrow(tbl), 
                      'red','black')) %>%
  ggplot(aes(Pn,fp)) +
  geom_line(linetype = 'dashed', color = 'grey33') +
  geom_point(size = 1, alpha = 0.9,
             aes(color = col)) +
  # scale_y_continuous(limits = c(-6,0)) +
  scale_color_manual(values = c('black','red')) +
  theme_classic(base_size = 16) +
  theme(legend.position = 'none') +
  labs(y = 'f(Pn)', 
       caption = caption)

