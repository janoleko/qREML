
energy_price = read.csv("~/Downloads/energy_ger2.csv", sep = ";", header = TRUE)
energy_price = read.csv("~/Downloads/energy_ger3.csv", sep = ";", header = TRUE)
# data from Bundesnetzagentur (https://www.smard.de/home/downloadcenter/download-marktdaten)

energy_gen = read.csv("~/Downloads/power_generation.csv", sep = ";", header = TRUE)
# data from Bundesnetzagentur (https://www.smard.de/home/downloadcenter/download-marktdaten)

# cleaning price data
library(stringr)

energy_price = energy_price[,c(1,2,8)]

colnames(energy_price) = c("date_from", "date_to", "price")

energy_price$price = as.numeric(stringr::str_replace(energy_price$price, ",", "."))

tail(energy_price)

plot(energy_price$price, type = "l")


library(fHMM)

eurdol = download_data("EURUSD=X", from = "2017-01-01", to = "2022-01-01")
eurdol$date = as.Date(eurdol$Date)


library(dplyr)

eurdol = data.frame(date = seq(as.Date(eurdol$date[1]), as.Date(eurdol$date[nrow(eurdol)]), by = "day")) %>% 
  full_join(eurdol, by = "date") %>% 
  mutate(close2 = na.approx(Close))

eurdol$date[1]


energy_price$date = as.Date(energy_price$date_from, format = "%d.%m.%Y")

energy2 = energy_price %>% left_join(eurdol, by = "date")

colnames(energy2)

energy3 = energy2 %>% dplyr::select(date = date, price = price, eurdol = close2)

energy4 = energy3[1:1500,]

plot(energy4$eurdol, energy4$price)


# download oil

oil = download_data("BZ=F", from = "2017-01-01", to = "2022-01-01")
oil$date = as.Date(oil$Date)

oil = data.frame(date = seq(oil$date[1], oil$date[nrow(oil)], by = "day")) %>% 
  full_join(oil, by = "date") %>% 
  mutate(oil = na.approx(Close)) %>% 
  dplyr::select(date, oil)

energy4 = energy3 %>% left_join(oil, by = "date")

energy5 = energy4 %>% filter(date < "2021-10-01")


# download gas
# download oil

gas = download_data("NG=F", from = "2017-01-01", to = "2022-01-01")
gas$date = as.Date(gas$Date)

gas = data.frame(date = seq(gas$date[1], gas$date[nrow(gas)], by = "day")) %>% 
  full_join(gas, by = "date") %>% 
  mutate(gas = na.approx(Close)) %>% 
  dplyr::select(date, gas)

energy6 = energy5 %>% left_join(gas, by = "date")

plot(energy6$date, energy6$price, type = "l")

plot(energy6$gas, energy6$price)


write.csv(energy6, "./data/energy_france.csv", row.names = FALSE)









# cleaning generation data
colnames(energy_gen) = c("date_from", "date_to",
                         "biomass", "hydro", "wind_offshore", "wind_onshore", "solar",
                         "other_renewables", "nuclear", "brown_coal", "hard_coal", "natural_gas", "pumped_storage", "other_conventional")

for(i in 3:ncol(energy_gen)){
  energy_gen[,i] = energy_gen[,i] %>% 
    str_remove("\\.") %>% 
    str_replace(",", ".") %>% 
    as.numeric()
}

energy_gen$wind = energy_gen$wind_offshore + energy_gen$wind_onshore
energy_gen$total = rowSums(energy_gen[,3:ncol(energy_gen)])

plot(energy_gen$total, energy_price$price, ylim = c(0, 200))

