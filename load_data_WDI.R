
# Загрузка пакетов
library('WDI')
library('data.table')
library('dplyr')

# NE.GDI.FTOT.KD         Gross fixed capital formation 
#                          (constant 2000 US$)
# NY.GDP.MKTP.KD         GDP (constant 2000 US$)
# SL.TLF.TOTL.IN         Labor force, total
# SP.POP.TOTL            Population, total
# GC.NFN.TOTL.GD.ZS      Net investment in nonfinancial 
#                          assets (% of GDP)
# NY.GNS.ICTR.ZS         Gross savings (% of GDP)

# Список стран без макрорегионов, групп стран и мира в целом
fileURL <- 'https://raw.githubusercontent.com/aksyuk/R-data/master/KeyBooks/data_countries.csv'
df_states <- read.csv2(fileURL, stringsAsFactors = F)

# стоимость ОПФ
df_K_raw = WDI(country = 'all', indicator = 'NE.GDI.FTOT.KD',
               start = 1960, end = 2019)
dt_K <- data.table(df_K_raw[df_K_raw$iso2c %in% df_states$ISO2, ])

# ВВП
df_Y_raw = WDI(country = 'all', indicator = 'NY.GDP.MKTP.KD', 
               start = 1960, end = 2019)
dt_Y <- data.table(df_Y_raw[df_Y_raw$iso2c %in% df_states$ISO2, ])

# Численность рабочей силы
df_L_raw = WDI(country = 'all', indicator = 'SL.TLF.TOTL.IN', 
               start = 1960, end = 2019)
dt_L <- data.table(df_L_raw[df_L_raw$iso2c %in% df_states$ISO2, ])

# Численность населения
df_pop_raw = WDI(country = 'all', indicator = 'SP.POP.TOTL', 
                 start = 1960, end = 2019)
dt_pop <- data.table(df_pop_raw[df_pop_raw$iso2c %in% df_states$ISO2, ])

# Чистые инвестиции, % от ВВП
df_I_raw = WDI(country = 'all', indicator = 'GC.NFN.TOTL.GD.ZS',
               start = 1960, end = 2019)
dt_I <- data.table(df_I_raw[df_I_raw$iso2c %in% df_states$ISO2, ])

# Норма сбережения, % от ВВП
df_s_raw = WDI(country = 'all', indicator = 'NY.GNS.ICTR.ZS',
               start = 1960, end = 2019)
dt_s <- data.table(df_s_raw[df_s_raw$iso2c %in% df_states$ISO2, ])

# Объединяем показатели в одну таблицу 
dt_all_data <- merge(dt_K, dt_Y, by = c('iso2c', 'country', 'year'))
dt_all_data <- merge(dt_all_data, dt_I, by = c('iso2c', 'country', 'year'))
dt_all_data <- merge(dt_all_data, dt_L, by = c('iso2c', 'country', 'year'))
dt_all_data <- merge(dt_all_data, dt_pop, by = c('iso2c', 'country', 'year'))
dt_all_data <- merge(dt_all_data, dt_s, by = c('iso2c', 'country', 'year'))

# данные по всем странам
str(dt_all_data)
summary(dt_all_data)

# # сохраняем объекты рабочего пространства
# save.image('./data/solow_all_states_data.RData')

# Фильтруем данные: оставить только США
DF <- filter(dt_all_data, iso2c == 'US')

# Приводим имена переменных в соответствие с требованиями MATLAB
colnames(DF) <- gsub('[.]', '_', colnames(DF))

# Считаем пропуски в статистике
sapply(DF, function(x){sum(is.na(x))})

# Записываем данные по США в файл
DF <- na.omit(DF)
write.csv2(DF, file = './data/solow_data_US.csv', row.names = F)
