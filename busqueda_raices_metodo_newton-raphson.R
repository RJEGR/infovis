# Intro ----
# Código para método de Newton-Raphson
# Estudiante: Ayerim Cinai Montelongo Cano
# 30 de Septiembre de 2021

# Usa TOL = 1 × 10−10, y para cada uno de los ejercicios elabora una tabla donde muestres el número de n, pn, f(pn). Al final, el programa debe indicar claramente el valor de la raíz.

# Funciones

# install.packages('Deriv')

library(Deriv)

options(digits = 5)

df <- list() 

raiz <- function(Po, TOL = 1E-10, N = 100) {
  
  n = 1
  
  # Función Newton - Raphson
  # Para la 1er iteración
  
  fd <- Deriv(f, 'x')
  
  Pn = Po - (f(Po) / fd(Po))
  
  df[[n]] <- data.frame(n, Pn, fp = f(Pn)) #almacenar en tabla
  
  # Iteraciones subsecuentes
  
  while(n < N) {
    
    n = n + 1
    
    if(abs(Po - Pn) >= TOL) {
      
      Po <- Pn
    
      Pn = Po - (f(Po) / fd(Po))
      
      df[[n]] <- data.frame(n, Pn, fp = f(Pn))
      
    } else {
      
      df <- do.call(rbind, df) # rbind combina filas y columnas de una lista
      
      Pn <- df[nrow(df), 2] # Ir al ultimo valor de la fila y devolver el valor de la columna 2
      
      cat('\nLa funcion tiene una raiz en', Pn, '\n')
      
      return(df)
      
      # stop()
    }

  }

}

pi = 3.1415926535897932

# Método de Newton-Raphson
#Mediante el método de Newton-Raphson, encuentra la raíz de:

# Ejercicio 1 ----

f <- function(x) {tan(pi*x) - 6}

Po = 0.23

raiz(Po) -> tbl1

# Ejercicio 2 ----

f <- function(x) {cos(x) - x}

Po = pi / 4

raiz(Po)  -> tbl2

# Ejercicio 3 ----

f <- function(x) {sin(sqrt(x)) - x}

Po = 1

raiz(Po) -> tbl3


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
