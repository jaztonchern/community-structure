library(stringdist)
library(vwr)
library(tidyverse)

#reads txt files
df1 <- read.delim("englishCPdatabase2.txt")
headers <- t(read.delim("clearpondHeaders_EN.txt"))

#extracts first row of data from column names
df1 <- rbind(colnames(df1),df1)

#insert correct headers
colnames(df1) <- headers

#extract only relevant columns, then cleans up entire df
df1 <- df1[,1:5]
df1[1,] <- gsub("X","",df1[1,])
df1[,3:5] <- lapply(df1[,3:5], function(x) as.numeric(as.character(x)))

trials <- 1:nrow(df1)
levenshtein.list = list()

#function to compute edit distance
for(i in 1:nrow(df1)){
  word <- df1$Word[i] #extracts target word
  levenshtein.vec <- levenshtein.distance(word,df1$Word) #computes edit distance between target word and all other words
  levenshtein.list[[i]] <- levenshtein.vec #saves the output into a list
}

