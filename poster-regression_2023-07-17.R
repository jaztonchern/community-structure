# script to rerun the regression to fix the following bugs:
# 1. the DV is the community SIZE and not community label (ordered by size) 
# 2. the correct phonological neighborhood size and frequency is obtained from CP 
# 3. remove number of letters as it is multicollinear with number of phonemes. 

# 2023-07-06 cogsci 2023 poster 
# note: this script does not resolve the duplicates issue when merging 

# updated 2023-07-17 to report correlations instead of regressions for poster 

# dutch ----

load('dutch_final.RData') # cp.analysis

# get rid of duplicates
cp.analysis.2 <- cp.analysis |> select(Word, Phono, Length_Phonomes, Length_Letters, Frequency, community) |> distinct()

# merge with the correct density and nhf values 
cp.full <- cp.full |> select(Word, dPTAN, dPTAF)
cp.analysis.2 <- cp.analysis.2 |> left_join(cp.full)

# recreate the random communities 
comm_sizes <- table(cp.analysis.2$community) |> as.data.frame()
colnames(comm_sizes) <- c('community', 'comm_size')
comm_sizes$community <- as.numeric(comm_sizes$community)
class(comm_sizes$community)
class(cp.analysis.2$community)
cp.analysis.2$community <- as.numeric(cp.analysis.2$community)

head(comm_sizes)
tail(comm_sizes)

set.seed(1)
test <- sample(cp.analysis.2$community)
test_sizes <- table(test) |> as.data.frame() |> arrange(Freq)

head(test_sizes)
tail(test_sizes)
colnames(test_sizes) <- c('random', 'comm_size')

head(test)
tail(test)

cp.analysis.2$random <- test

summary(cp.analysis.2)

class(cp.analysis.2$random)
test_sizes$random <- as.numeric(test_sizes$random)
class(test_sizes$random)

# analysis 

real.dat <- cp.analysis.2 |> left_join(comm_sizes)

summary(lm(comm_size~Length_Phonomes + Frequency + dPTAN + dPTAF, data=real.dat))
QuantPsyc::lm.beta(lm(comm_size~Length_Phonomes + Frequency + dPTAN + dPTAF,data=real.dat))

rdm.dat <- cp.analysis.2 |> left_join(test_sizes)

summary(lm(comm_size~Length_Phonomes + Frequency + dPTAN + dPTAF,data=rdm.dat))
QuantPsyc::lm.beta(lm(comm_size~Length_Phonomes + Frequency + dPTAN + dPTAF,data=rdm.dat))

### ====

# regression models have an anti- suppression effect - so it may be better to show the straight correlation between community size and the mean of the four lexical characteristics 
# the below cor.tests are at the community level (summarized data)

## real communities 
cor.phonemes <- real.dat |> group_by(community) |> summarize(mean = mean(Length_Phonomes)) |> left_join(comm_sizes)
cor.test(cor.phonemes$mean, cor.phonemes$comm_size)

cor.freq <- real.dat |> group_by(community) |> summarize(mean = mean(Frequency)) |> left_join(comm_sizes)
cor.test(cor.freq$mean, cor.freq$comm_size)

cor.degree <- real.dat |> group_by(community) |> summarize(mean = mean(dPTAN)) |> left_join(comm_sizes)
cor.test(cor.degree$mean, cor.degree$comm_size)

cor.nhf <- real.dat |> group_by(community) |> summarize(mean = mean(dPTAF)) |> left_join(comm_sizes)
cor.test(cor.nhf$mean, cor.nhf$comm_size)

## random communities 
cor.phonemes <- rdm.dat |> group_by(random) |> summarize(mean = mean(Length_Phonomes)) |> left_join(test_sizes)
cor.test(cor.phonemes$mean, cor.phonemes$comm_size)

cor.freq <- rdm.dat |> group_by(random) |> summarize(mean = mean(Frequency)) |> left_join(test_sizes)
cor.test(cor.freq$mean, cor.freq$comm_size)

cor.degree <- rdm.dat |> group_by(random) |> summarize(mean = mean(dPTAN)) |> left_join(test_sizes)
cor.test(cor.degree$mean, cor.degree$comm_size)

cor.nhf <- rdm.dat |> group_by(random) |> summarize(mean = mean(dPTAF)) |> left_join(test_sizes)
cor.test(cor.nhf$mean, cor.nhf$comm_size)

### ====

# french ----

load('french_final.RData')

# get rid of duplicates
cp.analysis.2 <- cp.analysis |> select(Word, Phono, Length_Phonomes, Length_Letters, Frequency, community) |> distinct()

# merge with the correct density and nhf values 
cp.full <- cp.full |> select(Word, fPTAN, fPTAF)
cp.analysis.2 <- cp.analysis.2 |> left_join(cp.full)

# recreate the random communities 
comm_sizes <- table(cp.analysis.2$community) |> as.data.frame()
colnames(comm_sizes) <- c('community', 'comm_size')
comm_sizes$community <- as.numeric(comm_sizes$community)
class(comm_sizes$community)
class(cp.analysis.2$community)
cp.analysis.2$community <- as.numeric(cp.analysis.2$community)

head(comm_sizes)
tail(comm_sizes)

set.seed(1)
test <- sample(cp.analysis.2$community)
test_sizes <- table(test) |> as.data.frame() |> arrange(Freq)

head(test_sizes)
tail(test_sizes)
colnames(test_sizes) <- c('random', 'comm_size')

head(test)
tail(test)

cp.analysis.2$random <- test

summary(cp.analysis.2)

class(cp.analysis.2$random)
test_sizes$random <- as.numeric(test_sizes$random)
class(test_sizes$random)

# analysis 

real.dat <- cp.analysis.2 |> left_join(comm_sizes)

summary(lm(comm_size~Length_Phonomes + Frequency + fPTAN + fPTAF, data=real.dat))
QuantPsyc::lm.beta(lm(comm_size~Length_Phonomes + Frequency + fPTAN + fPTAF,data=real.dat))

model <- lm(comm_size~Length_Phonomes + Frequency + fPTAN + fPTAF, data=real.dat)

rdm.dat <- cp.analysis.2 |> left_join(test_sizes)

summary(lm(comm_size~Length_Phonomes + Frequency + fPTAN + fPTAF,data=rdm.dat))
QuantPsyc::lm.beta(lm(comm_size~Length_Phonomes + Frequency + fPTAN + fPTAF,data=rdm.dat))

### ====

# regression models have an anti- suppression effect - so it may be better to show the straight correlation between community size and the mean of the four lexical characteristics 
# the below cor.tests are at the community level (summarized data)

## real communities 
cor.phonemes <- real.dat |> group_by(community) |> summarize(mean = mean(Length_Phonomes)) |> left_join(comm_sizes)
cor.test(cor.phonemes$mean, cor.phonemes$comm_size)

cor.freq <- real.dat |> group_by(community) |> summarize(mean = mean(Frequency)) |> left_join(comm_sizes)
cor.test(cor.freq$mean, cor.freq$comm_size)

cor.degree <- real.dat |> group_by(community) |> summarize(mean = mean(fPTAN)) |> left_join(comm_sizes)
cor.test(cor.degree$mean, cor.degree$comm_size)

cor.nhf <- real.dat |> group_by(community) |> summarize(mean = mean(fPTAF)) |> left_join(comm_sizes)
cor.test(cor.nhf$mean, cor.nhf$comm_size)

## random communities 
cor.phonemes <- rdm.dat |> group_by(random) |> summarize(mean = mean(Length_Phonomes)) |> left_join(test_sizes)
cor.test(cor.phonemes$mean, cor.phonemes$comm_size)

cor.freq <- rdm.dat |> group_by(random) |> summarize(mean = mean(Frequency)) |> left_join(test_sizes)
cor.test(cor.freq$mean, cor.freq$comm_size)

cor.degree <- rdm.dat |> group_by(random) |> summarize(mean = mean(fPTAN)) |> left_join(test_sizes)
cor.test(cor.degree$mean, cor.degree$comm_size)

cor.nhf <- rdm.dat |> group_by(random) |> summarize(mean = mean(fPTAF)) |> left_join(test_sizes)
cor.test(cor.nhf$mean, cor.nhf$comm_size)

### ====


# german ----

load('german_final.RData') # cp.analysis, dPTAF has a lot of missing data, almost 6k!  

# get rid of duplicates
cp.analysis.2 <- cp.analysis |> select(Word, Phono, Length_Phonomes, Length_Letters, Frequency, community) |> distinct()

# merge with the correct density and nhf values 
cp.full <- cp.full |> select(Word, gPTAN, gPTAF)
cp.analysis.2 <- cp.analysis.2 |> left_join(cp.full)

# recreate the random communities 
comm_sizes <- table(cp.analysis.2$community) |> as.data.frame()
colnames(comm_sizes) <- c('community', 'comm_size')
comm_sizes$community <- as.numeric(comm_sizes$community)
class(comm_sizes$community)
class(cp.analysis.2$community)
cp.analysis.2$community <- as.numeric(cp.analysis.2$community)

head(comm_sizes)
tail(comm_sizes)

set.seed(1)
test <- sample(cp.analysis.2$community)
test_sizes <- table(test) |> as.data.frame() |> arrange(Freq)

head(test_sizes)
tail(test_sizes)
colnames(test_sizes) <- c('random', 'comm_size')

head(test)
tail(test)

cp.analysis.2$random <- test

summary(cp.analysis.2)

class(cp.analysis.2$random)
test_sizes$random <- as.numeric(test_sizes$random)
class(test_sizes$random)

# analysis 

real.dat <- cp.analysis.2 |> left_join(comm_sizes)

summary(lm(comm_size~Length_Phonomes + Frequency + gPTAN + gPTAF, data=real.dat))
QuantPsyc::lm.beta(lm(comm_size~Length_Phonomes + Frequency + gPTAN + gPTAF,data=real.dat))

rdm.dat <- cp.analysis.2 |> left_join(test_sizes)

summary(lm(comm_size~Length_Phonomes + Frequency + gPTAN + gPTAF,data=rdm.dat))
QuantPsyc::lm.beta(lm(comm_size~Length_Phonomes + Frequency + gPTAN + gPTAF,data=rdm.dat))

### ====

# regression models have an anti- suppression effect - so it may be better to show the straight correlation between community size and the mean of the four lexical characteristics 
# the below cor.tests are at the community level (summarized data)

## real communities 
cor.phonemes <- real.dat |> group_by(community) |> summarize(mean = mean(Length_Phonomes)) |> left_join(comm_sizes)
cor.test(cor.phonemes$mean, cor.phonemes$comm_size)

cor.freq <- real.dat |> group_by(community) |> summarize(mean = mean(Frequency)) |> left_join(comm_sizes)
cor.test(cor.freq$mean, cor.freq$comm_size)

cor.degree <- real.dat |> group_by(community) |> summarize(mean = mean(gPTAN)) |> left_join(comm_sizes)
cor.test(cor.degree$mean, cor.degree$comm_size)

cor.nhf <- real.dat |> group_by(community) |> summarize(mean = mean(gPTAF)) |> left_join(comm_sizes)
cor.test(cor.nhf$mean, cor.nhf$comm_size)

## random communities 
cor.phonemes <- rdm.dat |> group_by(random) |> summarize(mean = mean(Length_Phonomes)) |> left_join(test_sizes)
cor.test(cor.phonemes$mean, cor.phonemes$comm_size)

cor.freq <- rdm.dat |> group_by(random) |> summarize(mean = mean(Frequency)) |> left_join(test_sizes)
cor.test(cor.freq$mean, cor.freq$comm_size)

cor.degree <- rdm.dat |> group_by(random) |> summarize(mean = mean(gPTAN)) |> left_join(test_sizes)
cor.test(cor.degree$mean, cor.degree$comm_size)

cor.nhf <- rdm.dat |> group_by(random) |> summarize(mean = mean(gPTAF)) |> left_join(test_sizes)
cor.test(cor.nhf$mean, cor.nhf$comm_size)

### ====

# spanish ----

load('spanish_final.RData') # cp.analysis, dPTAF has a lot of mising data as well. it is also discovered that the dPTAN value (which is presumably phonological density is incorrect... too many zeros)

# get rid of duplicates
cp.analysis.2 <- cp.analysis |> select(Word, Phono, Length_Phonomes, Length_Letters, Frequency, community) |> distinct()

# merge with the correct density and nhf values 
cp.full <- cp.full |> select(Word, sPTAN, sPTAF)
cp.analysis.2 <- cp.analysis.2 |> left_join(cp.full)

# recreate the random communities 
comm_sizes <- table(cp.analysis.2$community) |> as.data.frame()
colnames(comm_sizes) <- c('community', 'comm_size')
comm_sizes$community <- as.numeric(comm_sizes$community)
class(comm_sizes$community)
class(cp.analysis.2$community)
cp.analysis.2$community <- as.numeric(cp.analysis.2$community)

head(comm_sizes)
tail(comm_sizes)

set.seed(1)
test <- sample(cp.analysis.2$community)
test_sizes <- table(test) |> as.data.frame() |> arrange(Freq)

head(test_sizes)
tail(test_sizes)
colnames(test_sizes) <- c('random', 'comm_size')

head(test)
tail(test)

cp.analysis.2$random <- test

summary(cp.analysis.2)

class(cp.analysis.2$random)
test_sizes$random <- as.numeric(test_sizes$random)
class(test_sizes$random)

# analysis 

real.dat <- cp.analysis.2 |> left_join(comm_sizes)

summary(lm(comm_size~Length_Phonomes + Frequency + sPTAN + sPTAF, data=real.dat))
QuantPsyc::lm.beta(lm(comm_size~Length_Phonomes + Frequency + sPTAN + sPTAF,data=real.dat))

rdm.dat <- cp.analysis.2 |> left_join(test_sizes)

summary(lm(comm_size~Length_Phonomes + Frequency + sPTAN + sPTAF,data=rdm.dat))
QuantPsyc::lm.beta(lm(comm_size~Length_Phonomes + Frequency + sPTAN + sPTAF,data=rdm.dat))

### ====

# regression models have an anti- suppression effect - so it may be better to show the straight correlation between community size and the mean of the four lexical characteristics 
# the below cor.tests are at the community level (summarized data)

## real communities 
cor.phonemes <- real.dat |> group_by(community) |> summarize(mean = mean(Length_Phonomes)) |> left_join(comm_sizes)
cor.test(cor.phonemes$mean, cor.phonemes$comm_size)

cor.freq <- real.dat |> group_by(community) |> summarize(mean = mean(Frequency)) |> left_join(comm_sizes)
cor.test(cor.freq$mean, cor.freq$comm_size)

cor.degree <- real.dat |> group_by(community) |> summarize(mean = mean(sPTAN)) |> left_join(comm_sizes)
cor.test(cor.degree$mean, cor.degree$comm_size)

cor.nhf <- real.dat |> group_by(community) |> summarize(mean = mean(sPTAF)) |> left_join(comm_sizes)
cor.test(cor.nhf$mean, cor.nhf$comm_size)

## random communities 
cor.phonemes <- rdm.dat |> group_by(random) |> summarize(mean = mean(Length_Phonomes)) |> left_join(test_sizes)
cor.test(cor.phonemes$mean, cor.phonemes$comm_size)

cor.freq <- rdm.dat |> group_by(random) |> summarize(mean = mean(Frequency)) |> left_join(test_sizes)
cor.test(cor.freq$mean, cor.freq$comm_size)

cor.degree <- rdm.dat |> group_by(random) |> summarize(mean = mean(sPTAN)) |> left_join(test_sizes)
cor.test(cor.degree$mean, cor.degree$comm_size)

cor.nhf <- rdm.dat |> group_by(random) |> summarize(mean = mean(sPTAF)) |> left_join(test_sizes)
cor.test(cor.nhf$mean, cor.nhf$comm_size)

### ====
