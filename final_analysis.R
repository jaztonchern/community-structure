# set up - load libraries and data 
library(tidyverse)
library(langnetr)
library(mgsub)
library(stringr)
library(igraph)
library(rempsyc)

lexique <- readxl::read_excel("Lexique-xx.xlsx") %>% # where xx is the language
  dplyr::select(c(1,3,5,7,9,11,13))

cp.full <- read_tsv('xxxCPdatabase2.txt', col_names = FALSE) # where xxx is the language
head(cp.full, 2)

cp.col <- read.delim('clearpondHeaders_XX.txt') # where XX is the language
colnames(cp.full) <- cp.col$Headers # merge column headers with the database 

cp <- cp.full %>% dplyr::select(Word, Phono, Length_Phonomes) #select the necessary columns for cleaning
head(cp)

test.df <- cp # clones the data table for manipulation

#separate phono into single syllable
colmn <- paste("phono",1:max(test.df$Length_Phonomes))
test.df <- test.df %>% 
  separate(Phono,into=colmn)

#extract unique phoneme characters
unique.df <- reshape2::melt(test.df[,2:(max(test.df$Length_Phonomes)+1)],measure.vars=colmn,variable.name='position',value.name='phoneme')
unique.vec <- unique.df[,2] %>% unique()

#which phonemes have more than one character
unique2.vec <- unique.vec[which(nchar(unique.vec)>1)]

#which letters have not been used yet
replacement.vec <- c(letters[!(letters%in%unique.vec)],LETTERS[!(LETTERS%in%unique.vec)])

#replace 2-char phonemes with 1-char phonemes
cp$PhonoReplaced <- mgsub(cp$Phono,unique2.vec,replacement.vec[1:length(unique2.vec)])

#remove periods from phonemes
cp$PhonoReplaced <- gsub('.','',cp$PhonoReplaced,fixed=T)

#langnetr with phono
phono.network <- tolangnet(cp$PhonoReplaced)
phono.network.wlabels <- nodeindex(phono.network, cp$PhonoReplaced)
phono.network.edgelist <- toedgelist(phono.network.wlabels)

# Add phonemes back to each node by their matching name
phono.network.wlabels <- set_vertex_attr(phono.network.wlabels, "word", index = cp$PhonoReplaced, value = cp$Word)

#remove hermits and islands
hermits(phono.network.wlabels)
summary(phono.network.wlabels)
comp.g <- components(phono.network.wlabels)
phono.network.lcc <- induced_subgraph(phono.network.wlabels,v=comp.g$membership==1)

#community detection (Louvain)
phono.network.lcc.louvain <- cluster_louvain(phono.network.lcc)
lcc.louvain.map <- data.frame(node = V(phono.network.lcc)$name, community = phono.network.lcc.louvain$membership)
phono.network.lcc <- set_vertex_attr(phono.network.lcc, "membership", index = lcc.louvain.map$node, value = lcc.louvain.map$community)
phono.network.wlabels$louvain <- phono.network.lcc.louvain$membership
unique(phono.network.wlabels$louvain)

df <- data.frame()
for(i in unique(phono.network.wlabels$louvain)) { 
  #create subgraphs for each community 
  subgraph <- induced_subgraph(phono.network.wlabels,v=which(phono.network.wlabels$louvain==i)) 
  #get size of each subgraph 
  size <- gorder(subgraph) 
  #organise data into df
  df <- df %>% 
    bind_rows(data.frame(louvain=i,size=size)) 
}

cp.analysis <- cp.full %>% 
  dplyr::select(Word, Phono, Length_Phonomes, Length_Letters, Frequency, dPTAN, dPTAF) %>% 
  mutate(Names = cp$PhonoReplaced) %>% 
  filter(Names %in% phono.network.lcc.louvain$names) %>% 
  mutate(original_community = phono.network.lcc.louvain$membership) %>% 
  mutate(Frequency = log10(Frequency)) %>% 
  left_join(lexique, by = "Word") #add lexique data

cp.analysis1 <- cp.analysis %>% 
  group_by(original_community) %>% 
  summarise(N=n()) %>% 
  ungroup() %>% 
  arrange(N) %>% 
  rowid_to_column("community")

cp.analysis <- left_join(cp.analysis,cp.analysis1,by="original_community")
cp.analysis$community <- factor(cp.analysis$community)
cp.analysis$original_community <- factor(cp.analysis$original_community)

# comparison with a ER random network 
set.seed(1)
random_lcc <- replicate(10,erdos.renyi.game(gorder(phono.network.lcc), gsize(phono.network.lcc), type = 'gnm'))
random_lcc_louvain <- sapply(random_lcc,cluster_louvain)

# modularity and community structure descriptives -> table 1 with the other languages 
table1 <- data.frame (length(unique(cp.analysis$Word)), # number of words in the lcc
                      length(unique(cp.analysis$Word))/nrow(cp.full),
                      modularity(phono.network.lcc.louvain),
                      degree(phono.network.lcc) %>% mean(),
                      mean(sapply(random_lcc_louvain,modularity)), # Q is much lower
                      mean(sapply(random_lcc,degree)) %>% mean()) # same as phono.network.lcc
colnames(table1) <- c("No. of Words","% of total possible words","Q","degree","Q_rand","degree_rand")
write_csv(table1,"Q_dutch.csv")

# ---------------------------------------------------------
### 5. a. Raw Biphone Counts (Louvain Communities)
# ---------------------------------------------------------

#===============================
### 1. Create Real Communities
#===============================
# Extract only relevant variables from giant component network as data frame
g.lvn.bp <- data.frame(V(phono.network.lcc)$name,V(phono.network.lcc)$word,V(phono.network.lcc)$membership)
names(g.lvn.bp) = c("phonemes","name","membership")

# Split both phonemes and phonemes minus first char into biphones to get all non-overlapping biphones.
g.lvn.bp <- g.lvn.bp %>% mutate(
  phonemes_split = str_extract_all(phonemes, ".."),
  phonemes_minus1 = sub(".", "", phonemes),
  phonemes_minus1_split = str_extract_all(phonemes_minus1, "..")
)
# For str_extract_all: "." matches any character. Repeat twice to split by two chars.
# For sub: "." matches any character. Replace first char of all strings with nothing so we can then get the overlapping biphones of same phonemes.

# Get all phonemes in each community.
# Remember to use unlist after when grouping by membership
g.lvn.phonemes <- g.lvn.bp %>% group_by(membership) %>% 
  summarize(phonemes = c(
    unlist(phonemes_split),
    unlist(phonemes_minus1_split)))
g.lvn.phonemes[g.lvn.phonemes==""]<-NA
g.lvn.phonemes <- g.lvn.phonemes %>% 
  drop_na()
# This splits up the phonemes in each word into a column so they are no longer associated with just the word they came from. This is required because str_extract_all creates a list with a vector for each word.
# Also combines additional phonemes from phonemes_minus1_split!

# Extract data of phonemes by each desired community number (membership)
g.lvn.phonemes.dat <- data.frame(table(g.lvn.phonemes))
# Because raw biphone counts is a categorical variable, what we need to get is the number of times each biphone is seen in any word in each community, so we need to first split the data by membership then use the table function to get frequency for each biphone. Then, convert the table to a data frame so we can plot the graph of frequency of raw biphone counts for each community. Do not extract all communities but only focus on the few largest ones!

#===============================
### 2. Create Random Communities
#===============================
g.rdm <- data.frame(table(phono.network.lcc.louvain$membership))
names(g.rdm) = c("Random","Count")

g.rdm <- as.numeric(uncount(g.rdm, Count)[,1]) # Extract only column from data frame as vector, then convert the factor to numeric so sampling can happen

set.seed(1)
g.rdm <- sample(g.rdm) # when sample is called on a vector without any other parameters, it shuffles the elements inside it
g.rdm <- data.frame(node = V(phono.network.lcc)$name, random = g.rdm)
phono.network.lcc <- set_vertex_attr(phono.network.lcc, "random", index = g.rdm$node, value = g.rdm$random)
cp.analysis <- left_join(cp.analysis,g.rdm,by=c("Names"="node"))
cp.analysis$random_ordered <- cp.analysis1$community[match(cp.analysis$random, cp.analysis1$original_community)]

# Extract only relevant variables from giant component network as data frame.
g.rdm.bp <- data.frame(V(phono.network.lcc)$name,V(phono.network.lcc)$word,V(phono.network.lcc)$random)
names(g.rdm.bp) = c("phonemes","name","random")

# Split both phonemes and phonemes minus first char into biphones to get all non-overlapping biphones.
g.rdm.bp <- g.rdm.bp %>%  mutate(
  phonemes_split = str_extract_all(phonemes, ".."),
  phonemes_minus1 = sub(".", "", phonemes),
  phonemes_minus1_split = str_extract_all(phonemes_minus1, "..")
)

# Get all phonemes in each community.
g.rdm.phonemes <- g.rdm.bp %>% group_by(random) %>% 
  summarize(phonemes = c(
    unlist(phonemes_split),
    unlist(phonemes_minus1_split))) 
g.rdm.phonemes[g.rdm.phonemes==""]<-NA
g.rdm.phonemes <- g.rdm.phonemes %>% 
  drop_na()
# Extract data of phonemes by each desired community number (membership).
g.rdm.phonemes.dat <- data.frame(table(g.rdm.phonemes))

#===============================
### 3. Communities Analysis
#===============================
#descriptive statistics of lexical characteristics
cp.analysis.real <- cp.analysis %>% 
  group_by(community) %>% 
  summarise(N=n(),
            mean_phonemes = mean(Length_Phonomes),
            SD_phonemes = sd(Length_Phonomes),
            mean_letters = mean(Length_Letters),
            SD_letters = sd(Length_Letters),
            mean_freq = mean(Frequency),
            SD_freq = sd(Frequency),
            mean_PN = mean(dPTAN),
            SD_PN = sd(dPTAN),
            mean_PF = mean(dPTAF),
            SD_PF = sd(dPTAF)) %>% 
  ungroup()

cp.analysis.random <- cp.analysis %>% 
  group_by(random_ordered) %>% 
  summarise(N=n(),
            mean_phonemes = mean(Length_Phonomes),
            SD_phonemes = sd(Length_Phonomes),
            mean_letters = mean(Length_Letters),
            SD_letters = sd(Length_Letters),
            mean_freq = mean(Frequency),
            SD_freq = sd(Frequency),
            mean_PN = mean(dPTAN),
            SD_PN = sd(dPTAN),
            mean_PF = mean(dPTAF),
            SD_PF = sd(dPTAF)) %>% 
  ungroup() %>% 
  rename(community=random_ordered)

# multiple regression of lexical characteristics against community size
real.dat <- cp.analysis %>% 
  dplyr::select(3:7,10:15,community) %>% 
  mutate(community = as.numeric(community))
summary(lm(community~.,data=real.dat))
QuantPsyc::lm.beta(lm(community~.,data=real.dat))

rdm.dat <- cp.analysis %>% 
  dplyr::select(3:7,10:15,random)
summary(lm(random~.,data=rdm.dat))
QuantPsyc::lm.beta(lm(random~.,data=rdm.dat))

# k-s test of biphones in communities
ks.df <- data.frame(matrix(ncol=2)) #declare data table to store k-s test results
colnames(ks.df) <- c("D-statistic","p-value") #rename column names of data table
for(i in 1:nrow(cp.analysis1)){
  g.lvn.phonemes.test <- g.lvn.phonemes.dat[g.lvn.phonemes.dat$membership == i,] # extract data for each community
  g.rdm.phonemes.test <- g.rdm.phonemes.dat[g.rdm.phonemes.dat$random == i,]
  set.seed(1) #make results reproducible
  ks_score <- ks.test(g.lvn.phonemes.test$Freq,g.rdm.phonemes.test$Freq) # 2-way k-s test
  ks.df[i,] <- c(ks_score[[c(1,1)]],ks_score[[c(2,1)]]) # store k-s data in data table
}
ks.df <- rowid_to_column(ks.df,"original_community")
ks.df <- left_join(ks.df,cp.analysis1,by="original_community")

# plot of k-s test against community size
attach(ks.df)
plot(N, `p-value`, xaxt = "n",  xlab = "Community Size", ylab = "p-value", main = "Significance of K-S Tests and Community Size")
axis(1, at = seq(0, round(max(ks.df$N),3), by = 200), las=1)
abline(h=0.05, col = "red")
box() # add back box for barplot

# plot of k-s test d-statistic against p-value
plot(`D-statistic`, `p-value`, xaxt = "n",  xlab = "D-Statistic", ylab = "p-value", main = "Significance of K-S Tests and D-Statistic")
axis(1, at = seq(0, round(max(ks.df$`D-statistic`),3), by = 0.05), las=1)
abline(h=0.05, col = "red")
box() # add back box for barplot

#===============================
### 3. Barplots of Phonemes in Real vs Random Community
#===============================
# By first creating a factor for phonemes whereby they are sorted in decreasing order based on their frequencies, we can then plot them in decreasing order later.
reorder_size <- function(x) {
  x.sort = x[order(x$Freq, decreasing = TRUE),]
  x$phonemes = factor(x$phonemes, levels = x.sort$phonemes)
  return(x)
}
# We extract only words in target community.

# Real community 1
g.lvn.phonemes.1 <- g.lvn.phonemes.dat[g.lvn.phonemes.dat$membership == 1,]
g.lvn.phonemes.1$percentage <- g.lvn.phonemes.1$Freq/cp.analysis1$N[which(cp.analysis1$original_community==1)]
g.lvn.phonemes.1 = reorder_size(g.lvn.phonemes.1)

barplot(percentage ~ phonemes, data = g.lvn.phonemes.1, xaxt = "n",  xlab = "", ylab = "Raw biphone counts", main = "Raw biphone counts of real community 43")
title(xlab = "Biphones", line = 1) # move x-axis label closer to graph since x-axis category labels are removed
box() # add back box for barplot
# The plot is now in decreasing order of phonemes based on their frequencies!

# Random community 1
g.rdm.phonemes.1 <- g.rdm.phonemes.dat[g.rdm.phonemes.dat$random == 1,]
g.rdm.phonemes.1$percentage <- g.rdm.phonemes.1$Freq/cp.analysis1$N[which(cp.analysis1$original_community==1)]
# Note that we do not reorder the biphones for the random communities but rather base their order on the real communities because "The x-axis represents the different biphones found within these communities and the biphones (on both x-axes) were arranged based on their frequency of occurrence in the real community in descending order".
g.rdm.phonemes.1$phonemes = factor(g.rdm.phonemes.1$phonemes, levels = g.lvn.phonemes.1$phonemes)

barplot(percentage ~ phonemes, data = g.rdm.phonemes.1, xaxt = "n", xlab = "", ylab = "Raw biphone counts", main = "Raw biphone counts of random community 43")
title(xlab = "Biphones", line = 1)
box() # add back box for barplot


