zoo <- read.table("zoo.data", sep=",")

zoodata <- zoo[ , -c(1)]
zoodata <- zoodata[-27, ] # delete 2nd frog
rownames(zoodata) <- zoo$V1[-27]
colnames(zoodata) <- c("hair", "feathers", "eggs", "milk", "airborne", "aquatic", "predator",
                         "toothed", "backbone", "breathes", "venomous", "fins", "legs", "tail",
                         "domestic", "catsize", "type")
zoo <- zoodata

usethis::use_data(zoo, overwrite = TRUE)