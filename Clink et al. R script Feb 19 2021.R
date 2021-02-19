# This is the R script needed to recreate analyses presented in Clink et al 2021
# Limited evidence for individual signatures or site-level patterns of variation in male Northern gray gibbon (Hylobates funereus) duet codas 


# Table of contents -------------------------------------------------------
# Part 1. Table summarizing features
# Part 2. Supervised classification and data visualization
# Part 3. Multivariate, variance components model
# Part 4. Small-scale patterns of geographic variation


# Load necessary packages -------------------------------------------------
library(stringr)
library(dplyr)
library(ggplot2)
library(lme4)
library(bbmle)
library(plyr)
library(flextable)
library(e1071)
library(umap)
library(cowplot)
library(bayesplot)
library(rstan)
rstan_options(auto_write=TRUE)
options(mc.cores = parallel::detectCores())
library(spdep)


# Load in the data --------------------------------------------------------
all.mfcc.combined <- read.csv('all.mfcc.combined.csv')
combined.codas.all.sites <- read.csv('combined.codas.all.sites.csv')
male.coda.gps <- read.csv('male.coda.gps.csv')


# Part 1. Table summarizing features ------------------------------------

min.max <- data.frame(mean=sapply(combined.codas.all.sites[,c(3:22)],mean),
                      sd=sapply(combined.codas.all.sites[,c(3:22)],sd),
                      min=sapply(combined.codas.all.sites[,c(3:22)],min),
                      max=sapply(combined.codas.all.sites[,c(3:22)],max))

min.max$se <- min.max$sd / sqrt(nrow(min.max))

min.max <- round(cbind.data.frame(min.max),2)

min.max$mean.se <- paste(min.max$mean, '±', min.max$se)
min.max$range <- paste(min.max$min, '-', min.max$max,sep='')

NewMinMax <- min.max[,c('mean.se','range')] # paste(min.max$mean.se , '\n', min.max$range )


Feature <- c('Coda duration (s)',
             'Number of notes',
             'Minimum low frequency (Hz)',
             'Minimum high frequency (Hz)',
             'Minimum bandwidth (Hz)',
             'Maximum bandwidth (Hz)',
             'Mean minimum frequency (Hz)',
             'Mean maximum frequency (Hz)',
             'Maximum low frequency (Hz)',
             'Maximum high frequency (Hz)',
             'Mean bandwidth (Hz)',
             'Minimum note duration (s)',
             'Maximum note duration (s)',
             'Note rate (number of notes / duration)',
             'Note 1 duration (s)',
             'Note 1 minimum frequency (Hz)',
             'Note 1 maximum frequency (Hz)',
             'Note 2 duration (s)',
             'Note 2 minimum frequency (Hz)',
             'Note 2 maximum frequency (Hz)')

table.with.features <- cbind.data.frame(Feature,NewMinMax)
colnames(table.with.features) <- c('Feature','Mean ± SEM', 'Range')

myft <- flextable(
  (table.with.features))
myft <- width(myft, width = 1)
myft <- bold(myft, part = "header") 
myft

#save_as_docx(myft,path='Table 2.docx')


# Part 2.  Supervised classification and data visualization ---------------

## Part 2a. Supervised classification
# Convert all grouping variables to factor for SVM
combined.codas.all.sites$Site <- as.factor(combined.codas.all.sites$Site)
combined.codas.all.sites$individual <- as.factor(combined.codas.all.sites$individual)
all.mfcc.combined$site <- as.factor(all.mfcc.combined$site)
all.mfcc.combined$group <- as.factor( all.mfcc.combined$group)

# Supervised clustering using spectrogram features by site
ml.model.specfeatures.svm.site <- e1071::svm(combined.codas.all.sites[,c(3:22)], 
                                combined.codas.all.sites$Site, kernel = "radial", 
                                cross = 25)


ml.model.specfeatures.svm.site$tot.accuracy

# Supervised clustering using spectrogram features by individual
ml.model.specfeatures.svm <- e1071::svm(combined.codas.all.sites[,c(3:22)], 
                           combined.codas.all.sites$individual, kernel = "radial", 
                           cross = 25)


ml.model.specfeatures.svm$tot.accuracy

# Supervised clustering using MFCCs by site
ml.model.mfccs.svm.site <- e1071::svm(all.mfcc.combined[,-c(1,2,180)], 
                                all.mfcc.combined$site, kernel = "radial", 
                                cross = 25)


ml.model.mfccs.svm.site$tot.accuracy


# Supervised clustering using MFCCs by individual
ml.model.mfccs.svm <- e1071::svm(all.mfcc.combined[,-c(1,2,180)], 
                           all.mfcc.combined$group, kernel = "radial", 
                           cross = 25)

ml.model.mfccs.svm$tot.accuracy # LOOCV: 64.03712


## Part 2b. Using UMAP to visualize differences in individual males
male.individual.umap <- 
  umap::umap(combined.codas.all.sites[,c(3:22)],labels=as.factor(combined.codas.all.sites$individual),
             controlscale=TRUE,scale=3)

plot.for.male.individuals <-
  cbind.data.frame(male.individual.umap$layout[,1:2],
                   combined.codas.all.sites$individual)
colnames(plot.for.male.individuals) <-
  c("Dim.1", "Dim.2", "individual")

my_plot_umap_specfeatures_individ <-
  ggplot(data = plot.for.male.individuals, aes(
    x = Dim.1,
    y = Dim.2,
    colour = individual
  )) +
  geom_point(size = 3) +
  scale_color_manual(values = matlab::jet.colors (length(unique(plot.for.male.individuals$individual)))) +
  theme_bw() + ggtitle('Spectrogram features') + xlab('UMAP: Dim 1')+ylab('UMAP: Dim 2')+ theme(legend.position = "none")

my_plot_umap_specfeatures_individ


MaleCoda.umap <- 
  umap::umap(all.mfcc.combined.sub[,-c(1,2,180)],
             labels=as.factor(all.mfcc.combined.sub$group),
             controlscale=TRUE,scale=3)

recorder.id.df.umap <-
  cbind.data.frame(
    MaleCoda.umap$layout[,1:2],
    as.factor(all.mfcc.combined.sub$group)
  )


colnames(recorder.id.df.umap) <-
  c("Comp.1", "Comp.2", "group")


my_plot_umap_mfcc_individ <-
  ggplot(data = recorder.id.df.umap, aes(
    x = Comp.1,
    y = Comp.2,
    colour = group
  )) +
  geom_point(size = 3) +
  scale_color_manual(values = matlab::jet.colors(length(unique( as.factor(recorder.id.df.umap$group)))) ) +
  theme_bw() +xlab('UMAP: Dim 1')+ylab('UMAP: Dim 2')+ggtitle('Mel-frequency cepstral coefficients')+ theme(legend.position = "none")

print(my_plot_umap_mfcc_individ)

# Combine both plots
cowplot::plot_grid(my_plot_umap_specfeatures_individ,my_plot_umap_mfcc_individ)



# Part 3.  Multivariate, variance components model ------------------------

## Part 3a. Prepare data for model and run model
# Isolate relevant features from data set
# maxbw and max95 are highly correlated (0.9) so if we know max95 we have a pretty good idea of maxbw
# rest dur and call.dur also correlated (0.8)
d.manova <- combined.codas.all.sites[,c("call.dur", "note1dur","note1maxfreq",
                                        "note2dur","note2maxfreq",
                                        "noterate")]

# Check the structure of the data
cor(d.manova)

# Log transform data
d.manova <- log(d.manova)

# Create pairs plot for inspection of data
pairs(d.manova)

## Check the structure of the data
str(d.manova)

# Set-up data to pass to Stan. 
# Integer-coded vector of group IDs
group.int <- as.numeric(combined.codas.all.sites$individual)

# Check structure 
table(group.int)

# Integer-coded vector of site IDs
site.int <- as.numeric(as.factor(combined.codas.all.sites$site))

# Check structure 
table(site.int)

# Center data matrix at feature means
col.means <- apply(d.manova, MARGIN=2, FUN="mean")
y.centered <- sweep(d.manova, MARGIN=2, STATS=col.means)

# Create a data list to pass to Stan
data_list <- list(
  K = dim(d.manova)[2],
  J= length(unique(combined.codas.all.sites$individual)),
  M= length(unique(combined.codas.all.sites$site)),
  N= dim(d.manova)[1],
  y= as.matrix(y.centered), ## features centered at zero
  group= group.int,
  site= site.int
)


# Code to run the STAN model
# NOTE: the .stan file must be linked in the code below
# The code is commented out because it takes a long time to run! 
# mfinal.stan = stan(file="MANOVA.male codas.stan", 
#                    model_name = "mfinal", 
#                    data=data_list, iter=3000, warmup=1500, chains=2, 
#                    cores=2, 
#                    control = list(stepsize = 0.5, adapt_delta = 0.99, max_treedepth = 20))


# Optional code to save the output
#save(mfinal.stan, file = "male.coda.2020.rda")
load("male.coda.2020.rda")

# Check model output
# Create traceplots to check for mixing for site level variance
draws <- as.array(mfinal.stan, pars="ICC_site")
# mcmc_trace(draws)
# stan_dens(mfinal.stan, pars=c("ICC_site"))
round(summary(mfinal.stan, pars=c("ICC_site"))$summary, 3)

# Create traceplots to check for mixing for group level variance
draws <- as.array(mfinal.stan, pars="ICC_group")
# mcmc_trace(draws)
# stan_dens(mfinal.stan, pars=c("ICC_group"))
round(summary(mfinal.stan, pars=c("ICC_group"))$summary, 3)

# Check degrees of freedom parameter
# stan_dens(mfinal.stan, pars=c("DF_obs"))
# round(summary(mfinal.stan, pars=c("DF_obs"))$summary, 3)

## Extract the posterior samples
## These are permuted samples; we are talking permuted because output structure is more simple
## Once extracted can't use to check mixing because no longer in order along markov chain
post.samples <- rstan::extract(mfinal.stan,
                               pars= c("ICC_site", "ICC_group", "Maha_sqd", "DF_site","DF_group",
                                       "DF_obs", "Vcov_site", "Vcov_group", "Vcov_obs"),
                               permuted=TRUE)


# Create vector with feature names
features <- c("call.dur","note1dur", "note1maxfreq","note2dur", "note2maxfreq",
              "noterate")

# Add column names to samples
colnames(post.samples$ICC_group) <- features
colnames(post.samples$ICC_site) <- features

# Convert to dataframe
icc.site <- as.data.frame(post.samples$ICC_site)
icc.group <- as.data.frame(post.samples$ICC_group)

# Calculate ICC for observation-level and add column names
icc.obs <- as.data.frame(1 - icc.group - icc.site)
colnames(icc.obs) <- features


## Part 3b. Code to create ICC plots of posterior densities

## Note 1 duration
note1dur.site <- cbind.data.frame(icc.site$note1dur, rep("Site",length(icc.site$note1dur)))
note1dur.group <- cbind.data.frame(icc.group$note1dur,rep("Group",length(icc.group$note1dur)))
note1dur.obs <- cbind.data.frame(icc.obs$note1dur,rep("Obs",length(icc.obs$note1dur)))

colnames(note1dur.site) <- c("note1dur.samples","icc")
colnames(note1dur.group) <- c("note1dur.samples","icc")
colnames(note1dur.obs) <- c("note1dur.samples","icc")

note1dur.dens.df <- rbind.data.frame(note1dur.site,note1dur.group,note1dur.obs)
head(note1dur.dens.df)

note1dur.plot <- ggplot(note1dur.dens.df, aes(x=note1dur.samples, fill=icc)) + geom_density()+xlim(0,1)+#+ylim(0,4)+
  scale_fill_manual(name = "icc",
                    values = c(alpha("#0000FF", .9),
                               alpha("#00FFFF", .9),
                               alpha("#FFFF00", .9)))+  theme_classic() +
  guides(fill=F)+xlab("")+ylab("") +theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                                          axis.text.x  = element_text(size=18))+
  ggtitle("Note 1 Duration") + 
  theme(plot.title = element_text(size = 20, face = "bold"))

note1dur.plot

## Note 1 max freq

note1maxfreq.site <- cbind.data.frame(icc.site$note1maxfreq, rep("Between-site",length(icc.site$note1maxfreq)))
note1maxfreq.group <- cbind.data.frame(icc.group$note1maxfreq,rep("Between-male",length(icc.group$note1maxfreq)))
note1maxfreq.obs <- cbind.data.frame(icc.obs$note1maxfreq,rep("Within-male",length(icc.group$note1maxfreq)))

colnames(note1maxfreq.site) <- c("note1maxfreq.samples","icc")
colnames(note1maxfreq.group) <- c("note1maxfreq.samples","icc")
colnames(note1maxfreq.obs) <- c("note1maxfreq.samples","icc")

note1maxfreq.dens.df <- rbind.data.frame(note1maxfreq.site,note1maxfreq.group,note1maxfreq.obs)
head(note1maxfreq.dens.df)

note1maxfreq.plot <- ggplot(note1maxfreq.dens.df, aes(x=note1maxfreq.samples, fill=icc)) + geom_density()+xlim(0,1)+ #ylim(0,90)+
  scale_fill_manual(name = "ICC",
                    values = c(alpha("#0000FF", .9),
                               alpha("#00FFFF", .9),
                               alpha("#FFFF00", .9)))+  theme_classic() +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        axis.text.x  = element_text(size=18),legend.text=element_text(size=24))+
  xlab("")+ylab("") + 
  theme(legend.position = c(.85,1))+
  #guides(fill=F)+
  guides(fill = guide_legend(size=5, keywidth = 3, title="",keyheight = 1))+
  theme(legend.key.size =  unit(0.5, "in"))+ # Change key size in the legend
  theme(legend.text = element_text(size=16,face = "bold"))+
  ggtitle("Note 1 Max. Frequency") + 
  theme(plot.title = element_text(size = 20, face = "bold"))# Change the size labels in the legend 

note1maxfreq.plot



## Note 2 duration
note2dur.site <- cbind.data.frame(icc.site$note2dur, rep("Site",length(icc.site$note2dur)))
note2dur.group <- cbind.data.frame(icc.group$note2dur,rep("Group",length(icc.group$note2dur)))
note2dur.obs <- cbind.data.frame(icc.obs$note2dur,rep("Obs",length(icc.group$note2dur)))

colnames(note2dur.site) <- c("note2dur.samples","icc")
colnames(note2dur.group) <- c("note2dur.samples","icc")
colnames(note2dur.obs) <- c("note2dur.samples","icc")

note2dur.dens.df <- rbind.data.frame(note2dur.site,note2dur.group,note2dur.obs)
head(note2dur.dens.df)

note2dur.plot <- ggplot(note2dur.dens.df, aes(x=note2dur.samples, fill=icc)) + geom_density()+xlim(0,1)+ #ylim(0,90)+
  scale_fill_manual(name = "ICC",
                    values = c(alpha("#0000FF", .9),
                               alpha("#00FFFF", .9),
                               alpha("#FFFF00", .9)))+  theme_classic() +
  guides(fill=F)+xlab("")+ylab("") +theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                                          axis.text.x  = element_text(size=18))+
  ggtitle("Note 2 Duration") + 
  theme(plot.title = element_text(size = 20, face = "bold"))

note2dur.plot


## Note 2 max freq
note2maxfreq.site <- cbind.data.frame(icc.site$note2maxfreq, rep("Between-site",length(icc.site$note2maxfreq)))
note2maxfreq.group <- cbind.data.frame(icc.group$note2maxfreq,rep("Between-male",length(icc.group$note2maxfreq)))
note2maxfreq.obs <- cbind.data.frame(icc.obs$note2maxfreq,rep("Within-male",length(icc.group$note2maxfreq)))

colnames(note2maxfreq.site) <- c("note2maxfreq.samples","icc")
colnames(note2maxfreq.group) <- c("note2maxfreq.samples","icc")
colnames(note2maxfreq.obs) <- c("note2maxfreq.samples","icc")

note2maxfreq.dens.df <- rbind.data.frame(note2maxfreq.site,note2maxfreq.group,note2maxfreq.obs)
head(note2maxfreq.dens.df)

note2maxfreq.plot <- ggplot(note2maxfreq.dens.df, aes(x=note2maxfreq.samples, fill=icc)) + geom_density()+xlim(0,1)+ #ylim(0,90)+
  scale_fill_manual(name = "ICC",
                    values = c(alpha("#0000FF", .9),
                               alpha("#00FFFF", .9),
                               alpha("#FFFF00", .9)))+    theme_classic() +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        axis.text.x  = element_text(size=18),legend.text=element_text(size=24))+
  xlab("")+ylab("") + 
  guides(fill=F)+
  # theme(legend.position = c(.85,1))+
  # guides(fill = guide_legend(size=5, keywidth = 3, title="",keyheight = 1))+
  # theme(legend.key.size =  unit(0.5, "in"))+ # Change key size in the legend
  # theme(legend.text = element_text(size=16,face = "bold"))+
  ggtitle("Note 2 Max. Frequency") + 
  theme(plot.title = element_text(size = 20, face = "bold"))# Change

note2maxfreq.plot

## Call duration
call.dur.site <- cbind.data.frame(icc.site$call.dur, rep("Site",length(icc.site$call.dur)))
call.dur.group <- cbind.data.frame(icc.group$call.dur,rep("Group",length(icc.group$call.dur)))
call.dur.obs <- cbind.data.frame(icc.obs$call.dur,rep("Obs",length(icc.group$call.dur)))

colnames(call.dur.site) <- c("trill.samples","icc")
colnames(call.dur.group) <- c("trill.samples","icc")
colnames(call.dur.obs) <- c("trill.samples","icc")

call.dur.df <- rbind.data.frame(call.dur.site,call.dur.group,call.dur.obs)
head(call.dur.df)

call.dur.plot <- ggplot(call.dur.df, aes(x=trill.samples, fill=icc)) +xlim(0,1)+
  geom_density()+
  scale_fill_manual(name = "ICC",
                    values = c(alpha("#0000FF", .9),
                               alpha("#00FFFF", .9),
                               alpha("#FFFF00", .9)))+    theme_classic() +
  guides(fill=F)+ylab("")+xlab("Proportion of Variance") + theme(axis.title.x = element_text(size=16, face="bold"))+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        axis.text.x  = element_text(size=18))+ggtitle("Coda Duration") + 
  theme(plot.title = element_text(size = 20, face = "bold"))
call.dur.plot


## Note rate
noterate.site <- cbind.data.frame(icc.site$noterate, rep("Site",length(icc.site$noterate)))
noterate.group <- cbind.data.frame(icc.group$noterate,rep("Group",length(icc.group$noterate)))
noterate.obs <- cbind.data.frame(icc.obs$noterate,rep("Obs",length(icc.group$noterate)))

colnames(noterate.site) <- c("trill.samples","icc")
colnames(noterate.group) <- c("trill.samples","icc")
colnames(noterate.obs) <- c("trill.samples","icc")

trill.dens.df <- rbind.data.frame(noterate.site,noterate.group,noterate.obs)
head(trill.dens.df)

noterate.plot <- ggplot(trill.dens.df, aes(x=trill.samples, fill=icc)) + geom_density()+
  scale_fill_manual(name = "ICC",
                    values = c(alpha("#0000FF", .9),
                               alpha("#00FFFF", .9),
                               alpha("#FFFF00", .9)))+  theme_classic() +
  guides(fill=F)+ylab("")+xlab("Proportion of Variance") + theme(axis.title.x = element_text(size=16, face="bold"))+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        axis.text.x  = element_text(size=18))+xlim(0,1)+ggtitle("Note Rate") + 
  theme(plot.title = element_text(size = 20, face = "bold"))
noterate.plot


## Combine all plots

cowplot::plot_grid(hjust=-1.25, label_size = 15, vjust=1.5,note1dur.plot,note1maxfreq.plot,note2dur.plot,note2maxfreq.plot,call.dur.plot,noterate.plot,
                   #labels=c("Note 1 Dur","Note 1 Max Freq","Note 2 Dur","Note 2 Max Freq", "Call Dur", "Note Rate"),
                   ncol=2)


# Part 4.  Small-scale patterns of geographic variation -------------------
## Calculate geographic distance
gps.data.updated <- merge(combined.codas.all.sites,male.gps.data, by.x = 'individual', by.y = 'pair')
gps.data.updated <- droplevels(subset(gps.data.updated, site=="SA"))
geo_distance <- spDists(cbind(gps.data.updated$lon,gps.data.updated$lat),longlat = T)
geo_matrix <- as.matrix(geo_distance)
(geo_matrix[1,])

geo_matrix2 <- geo_matrix
geo_matrix2[lower.tri(geo_matrix2 , diag=TRUE)] <- NA
geo_vec <- na.exclude(as.vector(geo_matrix2))


##Calculate acoustic distance
# Check which features we want to include
colnames(gps.data.updated[,3:24])

pca.dists <- dist(gps.data.updated[,3:24])
acoustic.dist <- pca.dists
acou_matrix <- as.matrix(acoustic.dist)

aco_matrix2 <- acou_matrix
aco_matrix2[lower.tri(aco_matrix2, diag=TRUE)] <- NA
aco_vec <- na.exclude(as.vector(aco_matrix2))

residual_vectors <- cbind(aco_vec, geo_vec)
resid_vecs <- as.data.frame(residual_vectors)
str(resid_vecs)

resid_vecs <- subset(resid_vecs, subset=geo_vec >0)


### Calculate confidence intervals 
acoustic.dist.dataframe<- cbind.data.frame(resid_vecs$geo_vec,resid_vecs$aco_vec)

# Remove any observations that were closer than 500 m
new.acoustic.data.frame <- subset(acoustic.dist.dataframe, resid_vecs$geo_vec>=0.5)

x <- new.acoustic.data.frame$`resid_vecs$geo_vec`
y <- new.acoustic.data.frame$`resid_vecs$aco_vec`


my.spar <- 1.5 ## this is the smoothing parameter for all spline bootstraps
sp.frame <- data.frame(x=x,y=y)
sp.resampler <- function() {
  n <- nrow(sp.frame)
  resample.rows <- sample(1:n,size=n,replace=TRUE)
  return(sp.frame[resample.rows,])
}

sp.spline.estimator <- function(data,m=300) {
  # Fit spline to data, with cross-validation to pick lambda
  fit <- smooth.spline(x=data[,1],y=data[,2],spar=my.spar)
  # Set up a grid of m evenly-spaced points on which to evaluate the spline
  eval.grid <- seq(from=min(x),to=max(x),length.out=m)
  # Slightly inefficient to re-define the same grid every time we call this,
  # but not a big overhead
  # Do the prediction and return the predicted values
  return(predict(fit,x=eval.grid)$y) # We only want the predicted values
}

sp.spline.cis <- function(B,alpha,m=300) {
  spline.main <- sp.spline.estimator(sp.frame,m=m)
  # Draw B boottrap samples, fit the spline to each
  spline.boots <- replicate(B,sp.spline.estimator(sp.resampler(),m=m))
  # Result has m rows and B columns
  cis.lower <- 2*spline.main - apply(spline.boots,1,quantile,probs=1-alpha/2)
  cis.upper <- 2*spline.main - apply(spline.boots,1,quantile,probs=alpha/2)
  return(list(main.curve=spline.main,lower.ci=cis.lower,upper.ci=cis.upper,
              x=seq(from=min(x),to=max(x),length.out=m)))
}

sp.cis <- sp.spline.cis(B=1000,alpha=0.05)

plot(x,y,xlab="Geographic Distance (km)",
     ylab="Male Coda Dissimilarity",col= "grey80", pch =".", cex=.7)
#smooth.spline(y,x, cv=FALSE)
lines(x=sp.cis$x,y=sp.cis$main.curve, col="red")
lines(x=sp.cis$x,y=sp.cis$lower.ci)
lines(x=sp.cis$x,y=sp.cis$upper.ci)


