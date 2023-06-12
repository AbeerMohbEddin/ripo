install.packages("readxl")
library(readxl)
list.files()

data <- read.csv("filtered_file.csv")
data
library(dplyr)

selected_data <- data %>%
  select(contains("NP"))
library(dplyr)


write.csv(selected_data, "CGE_NP.csv",  quote = FALSE, row.names = FALSE)
selected_data1 <- data %>%
       select(contains("PE"))
 write.csv(selected_data1, "CGE_PE1.csv",  quote = FALSE, row.names = FALSE)