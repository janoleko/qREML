
# Petrel data preprocessing -----------------------------------------------

library(moveHMM)
library(sp)

# Load data from Movebank
# URL <- paste0("https://www.datarepository.movebank.org/bitstream/handle/10255/",
#               "move.568/At-sea%20distribution%20Antarctic%20Petrel%2c%",
#               "20Antarctica%202012%20%28data%20from%20Descamps%20et%20al.%",
#               "202016%29-gps.csv")
# raw <- read.csv(url(URL))
raw = read.csv("~/Downloads/petrels.csv")

# Keep relevant columns: ID, time, lon, lat
data_all <- raw[, c(13, 3, 4, 5)]
colnames(data_all) <- c("ID", "time", "lon", "lat")
data_all$time <- as.POSIXct(data_all$time, tz = "MST")

# Function to compute first-order differences for grouped data
diff_by_ID <- function(x, ID, ...) {
  # Indices of first and last value of each group n <- length(ID)
  i0 <- which(ID[-1] != ID[-n])
  i_first <- c(1, i0 + 1)
  i_last <- c(i0, n)
  # First-order differences
  dx <- rep(NA, n)
  dx[-i_last] <- difftime(time1 = x[-i_first], time2 = x[-i_last], ...)
  return(dx)
}

# Keep only a few tracks for this example (excluding tracks that # have unusually long intervals)
dtimes <- diff_by_ID(data_all$time, data_all$ID)
keep_ids <- setdiff(unique(data_all$ID),
                    unique(data_all$ID[which(dtimes > 30)]))[1:10]
data <- subset(data_all, ID %in% keep_ids)
data <- data[with(data, order(ID, time)),]

# Define centre for each track as first observation
i0 <- c(1, which(data$ID[-1] != data$ID[-nrow(data)]) + 1)
centres <- data[i0, c("ID", "lon", "lat")]
data$centre_lon <- rep(centres$lon, rle(data$ID)$lengths)
data$centre_lat <- rep(centres$lat, rle(data$ID)$lengths)

# Add distance to centre as covariate (based on sp for great circle distance)
data$d2c <- sapply(1:nrow(data), function(i) {
  spDistsN1(pts = matrix(as.numeric(data[i, c("lon", "lat")]), ncol = 2),
            pt = c(data$centre_lon[i], data$centre_lat[i]),
            longlat = TRUE)
})

# Remove unnecessary columns data$centre_lon <- NULL data$centre_lat <- NULL
# Derive step lengths and turning angles using moveHMM
movehmm_data <- prepData(trackData = data,
                         coordNames = c("lon", "lat"),
                         type = "LL")
data$step <- movehmm_data$step
data$angle <- movehmm_data$angle

# Replace zero step length to very small number because it's overkill 
# to use a zero-inflated distribution just for one zero observation 
wh_zero <- which(data$step == 0)
data$step[wh_zero] <- runif(length(wh_zero),
                            min = 0,
                            max = min(data$step[-wh_zero], na.rm = TRUE))

# Shorten track names
data$ID <- factor(data$ID)
levels(data$ID) <- paste0("PET-", LETTERS[1:length(unique(data$ID))])

write.csv(data, file = "./data/petrels.csv", row.names = FALSE)
