library(tidyverse)

# financiamiento de vivienda y perfil de audiencia (2020 - Enero - Julio)
x  <- read.csv("~/Downloads/202007.csv", encoding = "ISO-8859-1")
y  <- read.csv("~/Downloads/201912.csv", encoding = "ISO-8859-1")
z  <- read.csv("~/Downloads/201812.csv", encoding = "ISO-8859-1")

names(x)

# by_state_y <- y %>%
#   filter(CUMULATIVE_MONTH %in% c(1:7)) %>%
#   group_by(STATE) %>%
#   summarise(n = length(YEAR))


# by gender ----

x %>%
  filter(!GENERO %in% "No disponible") %>%
  group_by(GENERO) %>%
  summarise(n = length(YEAR)) %>%
  mutate(pct = n / sum(n) * 100, Year = "2020") -> by_gender_x

y %>%
  filter(CUMULATIVE_MONTH %in% c(1:7)) %>%
  filter(!GENERO %in% "No disponible") %>%
  group_by(GENERO) %>%
  summarise(n = length(YEAR)) %>%
  mutate(pct = n / sum(n) * 100, Year = "2019") -> by_gender_y

z %>%
  filter(CUMULATIVE_MONTH %in% c(1:7)) %>%
  filter(!GENERO %in% "No disponible") %>%
  group_by(GENERO) %>%
  summarise(n = length(YEAR)) %>%
  mutate(pct = n / sum(n) * 100, Year = "2018") -> by_gender_z

by_gender <- rbind(by_gender_x, 
                        by_gender_y, 
                        by_gender_z) %>%
  mutate(pct = round(pct))


# by age ----

x %>%
  filter(CUMULATIVE_MONTH %in% c(1:7)) %>%
  filter(!AGE_RANGE %in% "No disponible") %>%
  group_by(AGE_RANGE) %>%
  summarise(n = length(YEAR)) %>%
  mutate(pct = n / sum(n) * 100, Year = "2020") -> by_age_x

y %>%
  filter(CUMULATIVE_MONTH %in% c(1:7)) %>%
  filter(!AGE_RANGE %in% "No disponible") %>%
  group_by(AGE_RANGE) %>%
  summarise(n = length(YEAR)) %>%
  mutate(pct = n / sum(n) * 100, Year = "2019") -> by_age_y

z %>%
  filter(CUMULATIVE_MONTH %in% c(1:7)) %>%
  filter(!AGE_RANGE %in% "No disponible") %>%
  group_by(AGE_RANGE) %>%
  summarise(n = length(YEAR)) %>%
  mutate(pct = n / sum(n) * 100, Year = "2018") -> by_age_z

by_age <- rbind(by_age_x, 
                by_age_y, 
                by_age_z) %>%
  mutate(pct = round(pct)) %>%
  mutate(AGE_RANGE = as.character(AGE_RANGE)) %>%
  mutate(AGE_RANGE = ifelse(AGE_RANGE %in% "60 o m\xe1s", 
                            "60 o menos", AGE_RANGE)) %>%
  mutate(AGE_RANGE = as.factor(AGE_RANGE))
  

# by uma and gender

x %>%
  as_tibble() %>%
  filter(CUMULATIVE_MONTH %in% c(1:7)) %>%
  filter(!UMA_RANGE %in% "No disponible") %>%
  filter(!GENERO %in% "No disponible") %>%
  group_by(UMA_RANGE, GENERO) %>%
  summarise(n = length(YEAR)) %>%
  group_by(GENERO) %>%
  mutate(pct = n / sum(n) * 100, Year = "2020")
  
# by real_state_value ----

x %>%
  as_tibble() %>%
  filter(CUMULATIVE_MONTH %in% c(1:7)) %>%
  filter(!REAL_STATE_VALUE %in% "No disponible") %>%
  filter(!GENERO %in% "No disponible") %>%
  group_by(REAL_STATE_VALUE, GENERO) %>%
  summarise(n = length(YEAR)) %>%
  group_by(GENERO) %>%
  mutate(pct = n / sum(n) * 100, Year = "2020") -> real_st_x

y %>%
  as_tibble() %>%
  filter(CUMULATIVE_MONTH %in% c(1:7)) %>%
  filter(!REAL_STATE_VALUE %in% "No disponible") %>%
  filter(!GENERO %in% "No disponible") %>%
  group_by(REAL_STATE_VALUE, GENERO) %>%
  summarise(n = length(YEAR)) %>%
  group_by(GENERO) %>%
  mutate(pct = n / sum(n) * 100, Year = "2019") -> real_st_y

z %>%
  as_tibble() %>%
  filter(CUMULATIVE_MONTH %in% c(1:7)) %>%
  filter(!REAL_STATE_VALUE %in% "No disponible") %>%
  filter(!GENERO %in% "No disponible") %>%
  group_by(REAL_STATE_VALUE, GENERO) %>%
  summarise(n = length(YEAR)) %>%
  group_by(GENERO) %>%
  mutate(pct = n / sum(n) * 100, Year = "2018") -> real_st_z

levels <- c("Económica", "Media", 
            "Popular", "Tradicional", 
            "Residencial", "Residencial plus")

by_real_st <- rbind(real_st_x, 
                    real_st_y,
                    real_st_z) %>%
  mutate(pct = round(pct)) %>%
  mutate(REAL_STATE_VALUE = as.character(REAL_STATE_VALUE)) %>%
  mutate(REAL_STATE_VALUE = ifelse(REAL_STATE_VALUE %in% "Econ\xf3mica", 
                            "Económica", REAL_STATE_VALUE)) %>%
  mutate(REAL_STATE_VALUE = factor(REAL_STATE_VALUE, 
                                   levels = levels))

# price = as.numeric(PRICE)
# proyeccion ----

x %>%
  as_tibble() %>%
  filter(!REAL_STATE_VALUE %in% "No disponible") %>%
  filter(!CUMULATIVE_MONTH %in% "No disponible") %>%
  group_by(REAL_STATE_VALUE, CUMULATIVE_MONTH) %>%
  summarise(n = length(YEAR)) %>%
  mutate(Year = "2020") -> real_st_x

y %>%
  as_tibble() %>%
  filter(!REAL_STATE_VALUE %in% "No disponible") %>%
  filter(!CUMULATIVE_MONTH %in% "No disponible") %>%
  group_by(REAL_STATE_VALUE, CUMULATIVE_MONTH) %>%
  summarise(n = length(YEAR)) %>%
  mutate(Year = "2019") -> real_st_y

z %>%
  as_tibble() %>%
  filter(!REAL_STATE_VALUE %in% "No disponible") %>%
  filter(!CUMULATIVE_MONTH %in% "No disponible") %>%
  group_by(REAL_STATE_VALUE, CUMULATIVE_MONTH) %>%
  summarise(n = length(YEAR)) %>%
  mutate(Year = "2018") -> real_st_z

levels <- c("Económica", "Media", 
            "Popular", "Tradicional", 
            "Residencial", "Residencial plus")

by_real_st_proyect <- rbind(real_st_x, 
                    real_st_y,
                    real_st_z) %>%
  mutate(REAL_STATE_VALUE = as.character(REAL_STATE_VALUE)) %>%
  mutate(REAL_STATE_VALUE = ifelse(REAL_STATE_VALUE %in% "Econ\xf3mica", 
                                   "Económica", REAL_STATE_VALUE)) %>%
  mutate(REAL_STATE_VALUE = factor(REAL_STATE_VALUE, 
                                   levels = levels))

# micro-nicho
# table(x$STATE)

x %>%
  as_tibble() %>%
  filter(CUMULATIVE_MONTH %in% c(1:7)) %>%
  filter(STATE %in% "Baja California") %>%
  group_by() %>%
  group_by(MUNICIPIO) %>%
  summarise(n = length(YEAR)) %>%
  mutate(Year = "2020") -> municipio_x

y %>%
  as_tibble() %>%
  filter(CUMULATIVE_MONTH %in% c(1:7)) %>%
  filter(STATE %in% "Baja California") %>%
  group_by() %>%
  group_by(MUNICIPIO) %>%
  summarise(n = length(YEAR)) %>%
  mutate(Year = "2019") -> municipio_y

z %>%
  as_tibble() %>%
  filter(CUMULATIVE_MONTH %in% c(1:7)) %>%
  filter(STATE %in% "Baja California") %>%
  group_by() %>%
  group_by(MUNICIPIO) %>%
  summarise(n = length(YEAR)) %>%
  mutate(Year = "2018") -> municipio_z

by_municipio <- rbind(municipio_x, 
                      municipio_y,
                      municipio_z) %>%
  filter(!MUNICIPIO %in% "No distribuido") %>%
  group_by(Year) %>%
  mutate(pct = n / sum(n) * 100) %>%
  mutate(pct = round(pct))

# Modalidad de credito
x %>%
  as_tibble() %>%
  filter(CUMULATIVE_MONTH %in% c(1:7)) %>%
  group_by(Modalidad) %>%
  summarise(n = length(YEAR)) %>%
  mutate(Year = "2020") -> modalidad_x

y %>%
  as_tibble() %>%
  filter(CUMULATIVE_MONTH %in% c(1:7)) %>%
  group_by(Modalidad) %>%
  summarise(n = length(YEAR)) %>%
  mutate(Year = "2019") -> modalidad_y

z %>%
  as_tibble() %>%
  filter(CUMULATIVE_MONTH %in% c(1:7)) %>%
  group_by(Modalidad) %>%
  summarise(n = length(YEAR)) %>%
  mutate(Year = "2018") -> modalidad_z

by_modalidad <- rbind(modalidad_x, 
                      modalidad_y,
                      modalidad_z) %>%
  group_by(Year) %>%
  mutate(pct = n / sum(n) * 100) %>%
  mutate(pct = round(pct))

# datavis ----

by_gender %>%
  ggplot(aes(x = Year, y = pct)) +
  geom_col(aes(fill = GENERO)) +
  geom_label(aes(label = paste0(pct, "%"))) +
  theme_linedraw(base_family = "GillSans", base_size = 14) +
  ggsci::scale_fill_jama() +
  coord_flip() +
  labs(x = "", y = "", caption = "Financiamiento inmobiliario por género") +
  guides(fill = guide_legend(title = "", label.position = "top"))

setwd("~/Documents/Marketing/")

ggsave("gender_financiamiento.png", width = 6, height = 3)

library(ggrepel)

by_age %>%
  group_by(Year) %>% 
  mutate(lPos = cumsum(pct) - pct/2) %>%
  ggplot(aes(y = Year, x = pct)) +
  geom_col(aes(fill = AGE_RANGE),
           position = position_stack(reverse = T)) +
  geom_label_repel(aes(
    x = lPos,
    label = paste0(pct, "%")),
                   color = "black",
                   fill = "white", 
                   size = 3) + # nudge_y = 0.5
  theme_linedraw(base_family = "GillSans", base_size = 14) +
  ggsci::scale_fill_lancet() +
  labs(x = "", y = "", caption = "Financiamiento inmobiliaria por rango de edad") +
  guides(fill = guide_legend(title = "", label.position = "top"))

ggsave("age_financiamiento.png", width = 6, height = 3)

library(ggrepel)

by_real_st %>%
  ggplot(aes(x = 2, y = pct, fill = REAL_STATE_VALUE)) +
  geom_bar(stat = "identity",color="white") +
  geom_label_repel(aes(label= paste0(pct, "%")),
                   position=position_stack(vjust=0.5),
                   color = "white",size = 3) +
  coord_polar(theta="y") +
  xlim(0.5, 2.5) + # as dona
  labs(x = "", y = "", caption = "Valor inmobiliario a financiamiento") +
  guides(fill = guide_legend(title = "", label.position = "top")) +
  facet_grid(Year~GENERO) +
  theme_void(base_family = "GillSans", base_size = 16) +
  ggsci::scale_fill_uchicago()

ggsave("valor_inmobiliario.png", width = 6.5, height = 10)

by_real_st_proyect %>%
  ungroup() %>%
  mutate(CUMULATIVE_MONTH = factor(CUMULATIVE_MONTH)) %>%
  ggplot(aes(x = CUMULATIVE_MONTH, y = n, color = Year)) +
  geom_point() +
  geom_line(aes(group = Year)) +
  labs(x = "Meses ", y = "Inmobiliario finianciado", caption = "Proyección") +
  
  guides(color = guide_legend(title = "", label.position = "top")) +
  facet_wrap(~ REAL_STATE_VALUE) +
  theme_linedraw(base_family = "GillSans", base_size = 14) +
  ggsci::scale_color_uchicago()

ggsave("proyection.png", width = 10, height = 5)

#

by_real_st_proyect %>%
  group_by(CUMULATIVE_MONTH, Year) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  mutate(CUMULATIVE_MONTH = factor(CUMULATIVE_MONTH)) %>%
  ggplot(aes(x = CUMULATIVE_MONTH, y = n, color = Year)) +
  geom_point() +
  geom_line(aes(group = Year), size = 3) +
  labs(x = "Meses ", y = "Inmobiliario finianciado", caption = "Proyección") +
  
  guides(color = guide_legend(title = "", label.position = "top")) +
  # facet_wrap(~ REAL_STATE_VALUE) +
  theme_linedraw(base_family = "GillSans", base_size = 20) +
  ggsci::scale_color_uchicago()

ggsave("proyection_global.png", width = 10, height = 5)

by_municipio %>%
  mutate(lPos = cumsum(pct) - pct/2) %>%
  ggplot(aes(y = Year, x = pct)) +
  geom_col(aes(fill = MUNICIPIO),
           position = position_stack(reverse = T)) +   
  geom_label_repel(aes(
             x = lPos,
             label = paste0(pct, "%")),
             color = "black",
             fill = "white",
             size = 3) +
  theme_linedraw(base_family = "GillSans", base_size = 20) +
  ggsci::scale_fill_lancet() +
  labs(x = "", y = "", caption = "Financiamiento inmobiliario dentro del Estado de Baja California") +
  guides(fill = guide_legend(title = "", label.position = "left"))


ggsave("municipio_bc.png", width = 8, height = 5)

by_modalidad %>%
  ggplot(aes(x = 2, y = pct, fill = Modalidad)) +
  geom_bar(stat = "identity",color="white") +
  geom_label_repel(aes(label= paste0(pct, "%")),
                   position=position_stack(vjust=0.5),
                   color = "white",size = 3) +
  coord_polar(theta="y") +
  xlim(0.5, 2.5) + # as dona
  labs(x = "", y = "", caption = "Financiamiento dirigido a") +
  guides(fill = guide_legend(title = "", label.position = "top")) +
  facet_grid(~Year) +
  theme_void(base_family = "GillSans", base_size = 16) +
  ggsci::scale_fill_uchicago()

ggsave("modalidad_financiamiento.png", width = 10, height = 6)


# stringi::stri_enc_detect(as.character(x$STATE[2]))


#  financiamiento de vivienda y perfil de audiencia (2019 - Enero - Dic)
