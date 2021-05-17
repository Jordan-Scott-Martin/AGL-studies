####################################################################
###Prepare workspace 
####################################################################

#load packages
library(brms);library(rstan)
library(ggplot2);library(tidybayes)

#stan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#set working directory
setwd("")

####################################################################
###Copy model###
####################################################################

#load data
copy<-read.csv("copy data.csv")

#names
grammar<-colnames(copy[,c(17:31)])
# "AFirst" "ALast"  "BFirst" "BLast"  "Copy"   "CpyAS4" "CpyAS6" 
#"CpyAT4" "CpyAT6" "CpyBS4" "CpyBS6" "CpyBT4" "CpyBT6" "SomeA"  "SomeB" 

#build multi-response model for each grammar
AFirst.m <- bf(AFirst ~ Group + (1|cor|ID)) + bernoulli()
ALast.m <- bf(ALast ~ Group + (1|cor|ID)) + bernoulli()
BFirst.m <- bf(BFirst ~ Group + (1|cor|ID)) + bernoulli()
BLast.m <- bf(BLast ~ Group + (1|cor|ID)) + bernoulli()
Copy.m <- bf(Copy ~ Group + (1|cor|ID)) + bernoulli()
CpyAS4.m <- bf(CpyAS4 ~ Group + (1|cor|ID)) + bernoulli()
CpyAS6.m <- bf(CpyAS6 ~ Group + (1|cor|ID)) + bernoulli()
CpyAT4.m <- bf(CpyAT4 ~ Group + (1|cor|ID)) + bernoulli()
CpyAT6.m <- bf(CpyAT6 ~ Group + (1|cor|ID)) + bernoulli()
CpyBS4.m <- bf(CpyBS4 ~ Group + (1|cor|ID)) + bernoulli()
CpyBS6.m <- bf(CpyBS6 ~ Group + (1|cor|ID)) + bernoulli()
CpyBT4.m <- bf(CpyBT4 ~ Group + (1|cor|ID)) + bernoulli()
CpyBT6.m <- bf(CpyBT6 ~ Group + (1|cor|ID)) + bernoulli()
SomeA.m <- bf(SomeA ~ Group + (1|cor|ID)) + bernoulli()
SomeB.m <- bf(SomeB ~ Group + (1|cor|ID)) + bernoulli()

#run model
copymod<- brm(formula=AFirst.m + ALast.m + BFirst.m + BLast.m + Copy.m + CpyAS4.m + CpyAS6.m + CpyAT4.m + 
              CpyAT6.m + CpyBS4.m + CpyBS6.m + CpyBT4.m + CpyBT6.m + SomeA.m + SomeB.m + set_rescor(FALSE),
             data=copy,
             prior=c(prior("normal(0,2)",class="Intercept"),
                     prior("normal(0,2)",class="b"),
                     prior("cauchy(0,2)",class="sd",resp="AFirst"),
                     prior("cauchy(0,2)",class="sd",resp="ALast"),
                     prior("cauchy(0,2)",class="sd",resp="BFirst"),
                     prior("cauchy(0,2)",class="sd",resp="BLast"),
                     prior("cauchy(0,2)",class="sd",resp="Copy"),
                     prior("cauchy(0,2)",class="sd",resp="CpyAS4"),
                     prior("cauchy(0,2)",class="sd",resp="CpyAS6"),
                     prior("cauchy(0,2)",class="sd",resp="CpyAT4"),
                     prior("cauchy(0,2)",class="sd",resp="CpyAT6"),
                     prior("cauchy(0,2)",class="sd",resp="CpyBS4"),
                     prior("cauchy(0,2)",class="sd",resp="CpyBS6"),
                     prior("cauchy(0,2)",class="sd",resp="CpyBT4"),
                     prior("cauchy(0,2)",class="sd",resp="CpyBT6"),
                     prior("cauchy(0,2)",class="sd",resp="SomeA"),
                     prior("cauchy(0,2)",class="sd",resp="SomeB"),
                     prior("lkj(2)", class="cor")),
             
             warmup=1000,iter=3000, chains=4, seed=9,
             control=list(adapt_delta=0.99))

saveRDS(copymod, "copymodels.rds")

copymod<-readRDS("copymodels.rds")

######################################################################
#summary
copy.ef<-fixef(copymod, summary=TRUE, robust=TRUE, probs=c(0.05,0.95))

#post prob for hearing > deaf performance
pmcmc<-c(
  hypothesis(copymod, "AFirst_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(copymod, "ALast_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(copymod, "BFirst_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(copymod, "BLast_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(copymod, "Copy_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(copymod, "CpyAS4_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(copymod, "CpyAS6_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(copymod, "CpyAT4_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(copymod, "CpyAT6_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(copymod, "CpyBS4_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(copymod, "CpyBS6_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(copymod, "CpyBT4_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(copymod, "CpyBT6_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(copymod, "SomeA_GroupHEARING < 0")$hypothesis$Post.Prob,
  hypothesis(copymod, "SomeB_GroupHEARING < 0")$hypothesis$Post.Prob)

#logit sd
logitsd<-pi/sqrt(3)

#hold hearing/deaf comparisons
group.comparison<-data.frame(grammar=grammar, 
                             logO=copy.ef[c(16:30),1],
                             SE=copy.ef[c(16:30),2],
                             d=copy.ef[c(16:30),1]/logitsd,
                             odds=exp(copy.ef[c(16:30),1]),
                             CIl=copy.ef[c(16:30),3],CIu=copy.ef[c(16:30),4],
                             pmcmc=pmcmc)
row.names(group.comparison) <- NULL

group.comparison

#save results
write.csv(group.comparison, "copy_group_results.csv")

######################################################################
#grammar correlations
cor<-VarCorr(copymod, prob = c(0.05,0.95))
#90% CI for target grammar correlations
cor$ID$cor[,3:4,"Copy_Intercept"]

cor<-matrix(cor$ID$cor[1:15,1,1:15],15,15,dimnames=list(grammar,grammar))
cor<-round(cor,2)

write.csv(cor,"copy_cor.csv")

######################################################################
#plot group differences
library(tidyr)

#names
grammar<-colnames(copy[,c(17:31)])

#logistic function (log odds --> probabilities)
logistic <- function(x) 1/(1+exp(-x))

#deaf 
copy.post<-posterior_samples(copymod)
copy.deaf<-data.frame(copy.post[c(1:15)])
colnames(copy.deaf)<-grammar
copy.deafl<-gather(copy.deaf, "grammar", "OR", AFirst:SomeB, factor_key=TRUE)
copy.deafl$Group<-rep("deaf",nrow(copy.deafl))
copy.deafl$prob<-logistic(copy.deafl$OR)

#hearing OR
copy.hear<-data.frame(copy.post[c(16:30)])
colnames(copy.hear)<-grammar
copy.hear<-copy.deaf+copy.hear
copy.hearl<-gather(copy.hear, "grammar", "OR", AFirst:SomeB, factor_key=TRUE)
copy.hearl$Group<-rep("hearing",nrow(copy.hearl))
copy.hearl$prob<-logistic(copy.hearl$OR)

#summarize probabilities
aggregate(prob ~ grammar, FUN = median, data = copy.deafl)
aggregate(prob ~ grammar, FUN = mad, data = copy.deafl)
aggregate(prob ~ grammar, FUN = quantile, c(0.05,0.95), data = copy.deafl)
aggregate(prob ~ grammar, FUN = function(x) sum(x>0.5)/length(x), data = copy.deafl)

aggregate(prob ~ grammar, FUN = median, data = copy.hearl)
aggregate(prob ~ grammar, FUN = mad, data = copy.hearl)
aggregate(prob ~ grammar, FUN = quantile, c(0.05,0.95), data = copy.hearl)
aggregate(prob ~ grammar, FUN = function(x) sum(x>0.5)/length(x), data = copy.hearl)

#get logOR
copy.deaf_target = copy.deaf$Copy
copy.deaf_LOR = copy.deaf_target - copy.deaf
copy.hear_target = copy.hear$Copy
copy.hear_LOR = copy.hear_target - copy.hear

groupLOR =

 round(
  rbind(
  data.frame(
  LOR = apply(copy.deaf_LOR, 2, median),
  mad = apply(copy.deaf_LOR, 2, mad),
  l90 = apply(copy.deaf_LOR, 2, FUN = quantile, c(0.05,0.95))[1,],
  u90 = apply(copy.deaf_LOR, 2, FUN = quantile, c(0.05,0.95))[2,],
  pp = apply(copy.deaf_LOR, 2, FUN = function(x) sum(x>0)/length(x)) ),
  
  data.frame(
  LOR = apply(copy.hear_LOR, 2, median),
  mad = apply(copy.hear_LOR, 2, mad),
  l90 = apply(copy.hear_LOR, 2, FUN = quantile, c(0.05,0.95))[1,],
  u90 = apply(copy.hear_LOR, 2, FUN = quantile, c(0.05,0.95))[2,],
  pp = apply(copy.hear_LOR, 2, FUN = function(x) sum(x>0)/length(x)) )
  ), 2 )

groupLOR$group = rep(c("DEAF","HEARING"), each = length(grammar))

write.csv(groupLOR, "copy_groupLOR.csv")


#organize grammars
copydata<-rbind(copy.deafl,copy.hearl)
copydata$grammar <- factor(copydata$grammar, levels=grammar)

#reverse order of factor
copydata$Group = relevel(as.factor(copydata$Group), ref = "hearing")

#plot
library(tidybayes);library(cowplot)

copyplot<-

ggplot(copydata, aes(x=grammar, y =  prob, group = Group, color = Group))+
stat_pointinterval(.width = c(0.9), position=position_dodge(0.5), size=2, fatten_point=2)+
#coord_flip()+
scale_y_continuous(limits=c(0.2,1), breaks = c(0.25, 0.5, 0.75, 1))+
scale_color_manual(values=c("#000000","#9c9898"), labels = c("Hearing", "Deaf"))+
ggtitle(expression(paste(Copy))) +
ylab("Probability of grammar-consistent response\n")+
geom_hline(yintercept=0.5, linetype="dashed", size=0.5, color="#bab6b6")+
theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),unit="cm"),
        axis.ticks.y =element_blank(),
        axis.text.x =element_text(size=10, angle = 45, hjust = 1),
        axis.text.y =element_text(size=10),
        axis.title.x =element_blank(),
        axis.title.y =element_blank(),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size= 1, 
                                    linetype="solid"))#+
legend = get_legend(copyplot)

copyplot = copyplot + guides(color=FALSE)

save_plot("copy results.png", copyplot,
          base_aspect_ratio=1.4,
          dpi=300)

######################################################################
#individual-level plot
library(tidyr)

id<-unique(copy[,c("ID", "Group") ])
id<-id[order(id$ID),]
df<-data.frame(ID=id[,1],Group=id[,2])

pred<-fitted(copymod, scale = "response",newdata=df, robust=TRUE)
predm<-data.frame(pred[1:30,1,1:15])
predl90<-data.frame(pred[1:30,3,1:15])
predu90<-data.frame(pred[1:30,4,1:15])

#subject info
predm$ind<-id$ID
predm$Group<-id$Group

#wide to long
indml<-gather(predm, "grammar", "re", AFirst:SomeB, factor_key=TRUE)
indl90<-gather(predl90, "grammar", "l90", AFirst:SomeB, factor_key=TRUE)
indu90<-gather(predu90, "grammar", "u90", AFirst:SomeB, factor_key=TRUE)
indml$l90 = indl90$l90
indml$u90 = indu90$u90

#save results
write.csv(indml, "individual predictions copy.csv")

#individual level pp+
pred2<-fitted(copymod, scale = "response", newdata=df, robust=TRUE, summary = FALSE)
pp_target =apply(pred2, 2:3, FUN = function(x) sum(x > 0.5)/length(x) )

#proportion of hearing participants with p > 0.70 at pp+ >= 0.95
sum(pp_target[1:15,"Copy"] >= 0.95)/15

#proportion of deaf participants with p > 0.70 at pp+ >= 0.95
sum(pp_target[16:30,"Copy"] >= 0.95)/15

#reverse order of factor
indml$Group = relevel(as.factor(indml$Group), ref = "HEARING")

#plot
copyidplot<-
  
  ggplot(indml, aes(x=ind, y=re, group = Group, color = Group))+
  geom_point()+
  geom_errorbar(aes(ymin=l90,ymax=u90,width=0))+
  facet_wrap(~grammar, scales="fixed", nrow=3, ncol=5)+
  xlab("\nSubject ID")+
  ylab("Mean probability of consistent response\n")+
  geom_hline(yintercept=0.5, linetype="dashed", size=0.5, color="#bab6b6")+
  theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),unit="cm"),
        axis.ticks.x =element_blank(),
        axis.ticks.y =element_blank(),
        axis.text.x =element_text(size=10),
        axis.text.y =element_text(size=10),
        axis.title.x =element_text(size=12, face="bold"),
        axis.title.y =element_text(size=12, face="bold"),
        legend.title=element_text(size=12,face="bold"),
        panel.border = element_rect(fill=NA,color="black", size= 1, 
                                    linetype="solid"))


save_plot("copy ind results.png", copyidplot,
            base_width=12,base_height=7,
            dpi=300)


#just target grammar
copyid.p<-


  ggplot(indml[indml$grammar=="Copy",], aes(x=ind, y=re, group = Group, color = Group))+
  geom_point()+
  geom_errorbar(aes(ymin=l90,ymax=u90,width=0))+
  scale_color_manual(values=c("#000000","#9c9898"), labels = c("Hearing", "Deaf"))+
  scale_x_continuous(breaks=c(1,10,20,30))+
  scale_y_continuous(limits=c(0.2,1), breaks = c(0.25, 0.5, 0.75, 1))+
  xlab("\nSubject ID")+
  ggtitle(expression(paste(Copy))) +
  ylab("Probability of consistent response\n")+
  geom_hline(yintercept=0.5, linetype="dashed", size=0.5, color="#bab6b6")+
  theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),unit="cm"),
        plot.title=element_text(hjust=0.5, face="bold"),
        axis.ticks.x =element_blank(),
        axis.ticks.y =element_blank(),
        axis.text.x =element_text(size=10),
        axis.text.y =element_text(size=10),
        axis.title.x =element_text(size=12, face="bold"),
        axis.title.y =element_blank(),
        legend.title=element_text(size=12,face="bold"),
        legend.text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.85,0.20),
        panel.border = element_rect(fill=NA,color="black", size= 1, 
                                    linetype="solid"))+
    guides(color=FALSE)

#####################################################################
###Mirror model###
####################################################################

#load data
mirror<-read.csv("mirror data.csv")

#names
grammar<-colnames(mirror[,c(17:31)])
#AEdges4 AEdges6 AFirst ALast BEdges4 BEdges6 BFirst BLast Mid_AA4
#Mid_AA6 Mid_BB4 Mid_BB6 Mirror SomeA SomeB

#build multi-response model
AEdges4.m <- bf(AEdges4 ~ Group + (1|cor|ID)) + bernoulli()
AEdges6.m <- bf(AEdges6 ~ Group + (1|cor|ID)) + bernoulli()
AFirst.m <- bf(AFirst ~ Group + (1|cor|ID)) + bernoulli()
ALast.m <- bf(ALast ~ Group + (1|cor|ID)) + bernoulli()
BEdges4.m <- bf(BEdges4 ~ Group + (1|cor|ID)) + bernoulli()
BEdges6.m <- bf(BEdges6 ~ Group + (1|cor|ID)) + bernoulli()
BFirst.m <- bf(BFirst ~ Group + (1|cor|ID)) + bernoulli()
BLast.m <- bf(BLast ~ Group + (1|cor|ID)) + bernoulli()
Mid_AA4.m <- bf(Mid_AA4 ~ Group + (1|cor|ID)) + bernoulli()
Mid_AA6.m <- bf(Mid_AA6 ~ Group + (1|cor|ID)) + bernoulli()
Mid_BB4.m <- bf(Mid_BB4 ~ Group + (1|cor|ID)) + bernoulli()
Mid_BB6.m <- bf(Mid_BB6 ~ Group + (1|cor|ID)) + bernoulli()
Mirror.m <- bf(Mirror ~ Group + (1|cor|ID)) + bernoulli()
SomeA.m <- bf(SomeA ~ Group + (1|cor|ID)) + bernoulli()
SomeB.m <- bf(SomeB ~ Group + (1|cor|ID)) + bernoulli()

#run model
mirrormod<- brm(formula=AEdges4.m + AEdges6.m + AFirst.m + ALast.m + BEdges4.m + BEdges6.m + BFirst.m + BLast.m + 
                Mid_AA4.m + Mid_AA6.m + Mid_BB4.m + Mid_BB6.m + Mirror.m + SomeA.m + SomeB.m + set_rescor(FALSE),
              data=mirror,
              prior=c(prior("normal(0,2)",class="Intercept"),
                      prior("normal(0,2)",class="b"),
                      prior("cauchy(0,2)",class="sd",resp="AEdges4"),
                      prior("cauchy(0,2)",class="sd",resp="AEdges6"),
                      prior("cauchy(0,2)",class="sd",resp="AFirst"),
                      prior("cauchy(0,2)",class="sd",resp="ALast"),
                      prior("cauchy(0,2)",class="sd",resp="BEdges4"),
                      prior("cauchy(0,2)",class="sd",resp="BEdges6"),
                      prior("cauchy(0,2)",class="sd",resp="BFirst"),
                      prior("cauchy(0,2)",class="sd",resp="BLast"),
                      prior("cauchy(0,2)",class="sd",resp="MidAA4"),
                      prior("cauchy(0,2)",class="sd",resp="MidAA6"),
                      prior("cauchy(0,2)",class="sd",resp="MidBB4"),
                      prior("cauchy(0,2)",class="sd",resp="MidBB6"),
                      prior("cauchy(0,2)",class="sd",resp="Mirror"),
                      prior("cauchy(0,2)",class="sd",resp="SomeA"),
                      prior("cauchy(0,2)",class="sd",resp="SomeB"),
                      prior("lkj(2)", class="cor")),
              
              warmup=1000,iter=3000, chains=4, seed=9,
              control=list(adapt_delta=0.99))

saveRDS(mirrormod, "mirrormodels.rds")

mirrormod<-readRDS("mirrormodels.rds")

######################################################################
#summary
mirror.ef<-fixef(mirrormod, summary=TRUE, robust=TRUE, probs=c(0.05,0.95))

#p-value
pmcmc<-c(
  hypothesis(mirrormod, "AEdges4_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(mirrormod, "AEdges6_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(mirrormod, "AFirst_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(mirrormod, "ALast_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(mirrormod, "BEdges4_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(mirrormod, "BEdges6_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(mirrormod, "BFirst_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(mirrormod, "BLast_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(mirrormod, "MidAA4_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(mirrormod, "MidAA6_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(mirrormod, "MidBB4_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(mirrormod, "MidBB6_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(mirrormod, "Mirror_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(mirrormod, "SomeA_GroupHEARING < 0")$hypothesis$Post.Prob,
  hypothesis(mirrormod, "SomeB_GroupHEARING < 0")$hypothesis$Post.Prob)

#logit sd
logitsd<-pi/sqrt(3)

#hold hearing/deaf comparisons
group.comparison<-data.frame(grammar=grammar, 
                             lOR=mirror.ef[c(16:30),1],
                             SE=mirror.ef[c(16:30),2],
                             d=mirror.ef[c(16:30),1]/logitsd,
                             OR=exp(mirror.ef[c(16:30),1]),
                             CIl=mirror.ef[c(16:30),3],CIu=mirror.ef[c(16:30),4],
                             pmcmc=pmcmc)
row.names(group.comparison) <- NULL

#save results
write.csv(group.comparison, "mirror_group_results.csv")

######################################################################
#grammar correlations
cor<-VarCorr(mirrormod, prob = c(0.05,0.95))
#90% CI for target grammar correlations
cor$ID$cor[,3:4,"Mirror_Intercept"]

cor<-matrix(cor$ID$cor[1:15,1,1:15],15,15,dimnames=list(grammar,grammar))
cor<-round(cor,2)
write.csv(corr,"mirror_cor.csv")


######################################################################
#plot group differences
library(tidyr)

#names
grammar<-colnames(mirror[,c(17:31)])

#logistic function (log odds --> probabilities)
logistic <- function(x) 1/(1+exp(-x))

#deaf 
mir.post<-posterior_samples(mirrormod)
mir.deaf<-data.frame(mir.post[c(1:15)])
colnames(mir.deaf)<-grammar
mir.deafl<-gather(mir.deaf, "grammar", "OR", AEdges4:SomeB, factor_key=TRUE)
mir.deafl$Group<-rep("deaf",nrow(mir.deafl))
mir.deafl$prob<-logistic(mir.deafl$OR)

#hearing OR
mir.hear<-data.frame(mir.post[c(16:30)])
colnames(mir.hear)<-grammar
mir.hear<-mir.deaf+mir.hear
mir.hearl<-gather(mir.hear, "grammar", "OR", AEdges4:SomeB, factor_key=TRUE)
mir.hearl$Group<-rep("hearing",nrow(mir.hearl))
mir.hearl$prob<-logistic(mir.hearl$OR)

#summarize probabilities
aggregate(prob ~ grammar, FUN = median, data = mir.deafl)
aggregate(prob ~ grammar, FUN = mad, data = mir.deafl)
aggregate(prob ~ grammar, FUN = quantile, c(0.05,0.95), data = mir.deafl)
aggregate(prob ~ grammar, FUN = function(x) sum(x>0.5)/length(x), data = mir.deafl)

aggregate(prob ~ grammar, FUN = median, data = mir.hearl)
aggregate(prob ~ grammar, FUN = mad, data = mir.hearl)
aggregate(prob ~ grammar, FUN = quantile, c(0.05,0.95), data = mir.hearl)
aggregate(prob ~ grammar, FUN = function(x) sum(x>0.5)/length(x), data = mir.hearl)

#get logOR
mir.deaf_target = mir.deaf$Mirror
mir.deaf_LOR = mir.deaf_target - mir.deaf
mir.hear_target = mir.hear$Mirror
mir.hear_LOR = mir.hear_target - mir.hear

groupLOR =

 round(
  rbind(
  data.frame(
  LOR = apply(mir.deaf_LOR, 2, median),
  mad = apply(mir.deaf_LOR, 2, mad),
  l90 = apply(mir.deaf_LOR, 2, FUN = quantile, c(0.05,0.95))[1,],
  u90 = apply(mir.deaf_LOR, 2, FUN = quantile, c(0.05,0.95))[2,],
  pp = apply(mir.deaf_LOR, 2, FUN = function(x) sum(x>0)/length(x)) ),
  
  data.frame(
  LOR = apply(mir.hear_LOR, 2, median),
  mad = apply(mir.hear_LOR, 2, mad),
  l90 = apply(mir.hear_LOR, 2, FUN = quantile, c(0.05,0.95))[1,],
  u90 = apply(mir.hear_LOR, 2, FUN = quantile, c(0.05,0.95))[2,],
  pp = apply(mir.hear_LOR, 2, FUN = function(x) sum(x>0)/length(x)) )
  ), 2 )

groupLOR$group = rep(c("DEAF","HEARING"), each = length(grammar))

write.csv(groupLOR, "mirror_groupLOR.csv")

#organize grammars
mirdata<-rbind(mir.deafl,mir.hearl)

mirdata$grammar <- factor(mirdata$grammar, levels=grammar)

#reverse order of factor
mirdata$Group = relevel(as.factor(mirdata$Group), ref = "hearing")

#plot
library(tidybayes);library(cowplot)

mirplot<-

ggplot(mirdata, aes(x=grammar, y =  prob, group = Group, color = Group))+
stat_pointinterval(.width = c(0.9), position=position_dodge(0.5), size=2, fatten_point=2)+
#coord_flip()+
#scale_x_continuous(breaks=c(1,10,20,30))+
scale_y_continuous(limits=c(0.2,1), breaks = c(0.25, 0.5, 0.75, 1))+
scale_color_manual(values=c("#000000","#9c9898"), labels = c("Hearing", "Deaf"))+
ggtitle(expression(paste(Mirror))) +
ylab("Probability of grammar-consistent response\n")+
geom_hline(yintercept=0.5, linetype="dashed", size=0.5, color="#bab6b6")+
theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),unit="cm"),
        axis.ticks.y =element_blank(),
        axis.text.x =element_text(size=10, angle = 45, hjust = 1),
        axis.text.y =element_text(size=10),
        axis.title.x =element_blank(),
        axis.title.y =element_blank(),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size= 1, 
                                    linetype="solid"))+
  guides(color=FALSE)


save_plot("mirror results_2.png", mirplot,
          base_aspect_ratio=1.4,
          dpi=300)

######################################################################
#individual-level plot

id<-unique(mirror[,c("ID", "Group") ])
id<-id[order(id$ID),]
df<-data.frame(ID=id[,1],Group=id[,2])

pred<-fitted(mirrormod, scale = "response",newdata=df, robust=TRUE)
predm<-data.frame(pred[1:30,1,1:15])
predl90<-data.frame(pred[1:30,3,1:15])
predu90<-data.frame(pred[1:30,4,1:15])

#subject info
predm$ind<-id$ID
predm$Group<-id$Group

#wide to long
indml<-gather(predm, "grammar", "re", AEdges4:SomeB, factor_key=TRUE)
indl90<-gather(predl90, "grammar", "l90", AEdges4:SomeB, factor_key=TRUE)
indu90<-gather(predu90, "grammar", "u90", AEdges4:SomeB, factor_key=TRUE)
indml$l90 = indl90$l90
indml$u90 = indu90$u90

#save results
write.csv(indml, "individual predictions mirror.csv")

#individual level pp+
pred2<-fitted(mirrormod, scale = "response", newdata=df, robust=TRUE, summary = FALSE)
pp_target =apply(pred2, 2:3, FUN = function(x) sum(x > 0.5)/length(x) )

#proportion of hearing participants with p > 0.70 at pp+ >= 0.95
sum(pp_target[1:15,"Mirror"] >= 0.95)/15

#proportion of deaf participants with p > 0.70 at pp+ >= 0.95
sum(pp_target[16:30,"Mirror"] >= 0.95)/15

#reverse order of factor
indml$Group = relevel(as.factor(indml$Group), ref = "HEARING")

#plot
miridplot<-
  
  ggplot(indml, aes(x=ind, y=re, group = Group, color = Group))+
  geom_point()+
  geom_errorbar(aes(ymin=l90,ymax=u90,width=0))+
  facet_wrap(~grammar, scales="fixed", nrow=3, ncol=5)+
  xlab("\nSubject ID")+
  ylab("Mean probability of consistent response\n")+
  geom_hline(yintercept=0.5, linetype="dashed", size=0.5, color="#bab6b6")+
  theme(plot.margin=unit(c(1,1,1,1),unit="cm"),
        axis.ticks.x =element_blank(),
        axis.ticks.y =element_blank(),
        axis.text.x =element_text(size=10),
        axis.text.y =element_text(size=10),
        axis.title.x =element_text(size=12, face="bold"),
        axis.title.y =element_text(size=12, face="bold"),
        legend.title=element_text(size=12,face="bold"),
        panel.border = element_rect(fill=NA,color="black", size= 1, 
                                    linetype="solid"))

save_plot("mirror ind results.png", miridplot,
          base_width=12,base_height=7,
          dpi=300)

#just target grammar
mirid.p<-
  
  ggplot(indml[indml$grammar=="Mirror",], aes(x=ind, y=re, group = Group, color = Group))+
  geom_point()+
  geom_errorbar(aes(ymin=l90,ymax=u90,width=0))+
  scale_color_manual(values=c("#000000","#9c9898"), labels = c("Hearing", "Deaf"))+
  scale_x_continuous(breaks=c(1,10,20,30))+
  scale_y_continuous(limits=c(0.2,1), breaks = c(0.25, 0.5, 0.75, 1))+
  xlab("\nSubject ID")+
  ggtitle(expression(paste(Mirror))) +
  ylab("Probability of consistent response\n")+
  geom_hline(yintercept=0.5, linetype="dashed", size=0.5, color="#bab6b6")+
  theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),unit="cm"),
        plot.title=element_text(hjust=0.5, face="bold"),
        axis.ticks.x =element_blank(),
        axis.ticks.y =element_blank(),
        axis.text.x =element_text(size=10),
        axis.text.y =element_text(size=10),
        axis.title.x =element_text(size=12, face="bold"),
        axis.title.y =element_blank(),
        legend.title=element_text(size=12,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size= 1, 
                                    linetype="solid"))+
    guides(color=FALSE)


#####################################################################
###abna model###
####################################################################

#load data
abna<-read.csv("abna data.csv")

#names
grammar<-colnames(abna[,c(17:29)])
#ABFirst ABnA AEdge. AEdge..1

#build multi-response model
ABFirst.m <- bf(ABFirst ~ Group + (1|cor|ID)) + bernoulli()
ABnA.m <- bf(ABnA ~ Group + (1|cor|ID)) + bernoulli()
AEdge.m <- bf(AEdge. ~ Group + (1|cor|ID)) + bernoulli()
AEdge1.m <- bf(AEdge..1 ~ Group + (1|cor|ID)) + bernoulli()
AFirst.m <- bf(AFirst ~ Group + (1|cor|ID)) + bernoulli()
ALast.m <- bf(ALast ~ Group + (1|cor|ID)) + bernoulli()
BALast.m <- bf(BALast ~ Group + (1|cor|ID)) + bernoulli()
BFirst.m <- bf(BFirst ~ Group + (1|cor|ID)) + bernoulli()
BLast.m <- bf(BLast ~ Group + (1|cor|ID)) + bernoulli()
MidB.m <- bf(MidB ~ Group + (1|cor|ID)) + bernoulli()
MidBB.m <- bf(MidBB ~ Group + (1|cor|ID)) + bernoulli()
SomeA.m <- bf(SomeA ~ Group + (1|cor|ID)) + bernoulli()
SomeB.m <- bf(SomeB ~ Group + (1|cor|ID)) + bernoulli()

#run model
abnamod<- brm(formula=ABFirst.m + ABnA.m + AEdge.m + AEdge1.m + AFirst.m + ALast.m + BALast.m + BFirst.m +
                BLast.m + MidB.m + MidBB.m + SomeA.m + SomeB.m + set_rescor(FALSE),
                data=abna,prior=c(prior("normal(0,2)",class="Intercept"),
                          prior("normal(0,2)",class="b"),
                          prior("cauchy(0,2)",class="sd",resp="ABFirst"),
                          prior("cauchy(0,2)",class="sd",resp="ABnA"),
                          prior("cauchy(0,2)",class="sd",resp="AEdge"),
                          prior("cauchy(0,2)",class="sd",resp="AEdge1"),
                          prior("cauchy(0,2)",class="sd",resp="AFirst"),
                          prior("cauchy(0,2)",class="sd",resp="ALast"),
                          prior("cauchy(0,2)",class="sd",resp="BALast"),
                          prior("cauchy(0,2)",class="sd",resp="BFirst"),
                          prior("cauchy(0,2)",class="sd",resp="BLast"),
                          prior("cauchy(0,2)",class="sd",resp="MidB"),
                          prior("cauchy(0,2)",class="sd",resp="MidBB"),
                          prior("cauchy(0,2)",class="sd",resp="SomeA"),
                          prior("cauchy(0,2)",class="sd",resp="SomeB"),
                          prior("lkj(2)", class="cor")),
                
                warmup=1000,iter=3000, chains=4, seed=9,
                control=list(adapt_delta=0.99))

saveRDS(abnamod, "abnamodels.rds")

abnamod<-readRDS("abnamodels.rds")

######################################################################
#summary
abna.ef<-fixef(abnamod, summary=TRUE, robust=TRUE, probs=c(0.05,0.95))

#p-value
pmcmc<-c(
  hypothesis(abnamod, "ABFirst_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(abnamod, "ABnA_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(abnamod, "AEdge_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(abnamod, "AEdge1_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(abnamod, "AFirst_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(abnamod, "ALast_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(abnamod, "BALast_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(abnamod, "BFirst_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(abnamod, "BLast_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(abnamod, "MidB_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(abnamod, "MidBB_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(abnamod, "SomeA_GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(abnamod, "SomeB_GroupHEARING > 0")$hypothesis$Post.Prob)

#logit sd
logitsd<-pi/sqrt(3)

#hold hearing/deaf comparisons
group.comparison<-data.frame(grammar=grammar, 
                             lOR=abna.ef[c(14:26),1],
                             SE=abna.ef[c(14:26),2],
                             d=abna.ef[c(14:26),1]/logitsd,
                             OR=exp(abna.ef[c(14:26),1]),
                             CIl=abna.ef[c(14:26),3],CIu=abna.ef[c(14:26),4],
                             pmcmc=pmcmc)
row.names(group.comparison) <- NULL

#save results
write.csv(group.comparison, "abna_group_results.csv")

######################################################################
#grammar correlations
cor<-VarCorr(abnamod, prob = c(0.05,0.95))
#90% CI for target grammar correlations
cor$ID$cor[,3:4,"ABnA_Intercept"]

cor<-matrix(cor$ID$cor[1:13,1,1:13],13,13,dimnames=list(grammar,grammar))
cor<-round(cor,2)
write.csv(cor,"abna_cor.csv")

######################################################################
#plot group differences
library(tidyr)

#names
grammar<-colnames(abna[,c(17:29)])

#logistic function (log odds --> probabilities)
logistic <- function(x) 1/(1+exp(-x))

#deaf 
abna.post<-posterior_samples(abnamod)
abna.deaf<-data.frame(abna.post[c(1:13)])
colnames(abna.deaf)<-grammar
abna.deafl<-gather(abna.deaf, "grammar", "OR", ABFirst:SomeB, factor_key=TRUE)
abna.deafl$Group<-rep("deaf",nrow(abna.deafl))
abna.deafl$prob<-logistic(abna.deafl$OR)

#hearing OR
abna.hear<-data.frame(abna.post[c(14:26)])
colnames(abna.hear)<-grammar
abna.hear<-abna.deaf+abna.hear
abna.hearl<-gather(abna.hear, "grammar", "OR", ABFirst:SomeB, factor_key=TRUE)
abna.hearl$Group<-rep("hearing",nrow(abna.hearl))
abna.hearl$prob<-logistic(abna.hearl$OR)

#summarize probabilities
aggregate(prob ~ grammar, FUN = median, data = abna.deafl)
aggregate(prob ~ grammar, FUN = mad, data = abna.deafl)
aggregate(prob ~ grammar, FUN = quantile, c(0.05,0.95), data = abna.deafl)
aggregate(prob ~ grammar, FUN = function(x) sum(x>0.5)/length(x), data = abna.deafl)

aggregate(prob ~ grammar, FUN = median, data = abna.hearl)
aggregate(prob ~ grammar, FUN = mad, data = abna.hearl)
aggregate(prob ~ grammar, FUN = quantile, c(0.05,0.95), data = abna.hearl)
aggregate(prob ~ grammar, FUN = function(x) sum(x>0.5)/length(x), data = abna.hearl)

#get logOR
abna.deaf_target = abna.deaf$ABnA
abna.deaf_LOR = abna.deaf_target - abna.deaf
abna.hear_target = abna.hear$ABnA
abna.hear_LOR = abna.hear_target - abna.hear

groupLOR =

 round(
  rbind(
  data.frame(
  LOR = apply(abna.deaf_LOR, 2, median),
  mad = apply(abna.deaf_LOR, 2, mad),
  l90 = apply(abna.deaf_LOR, 2, FUN = quantile, c(0.05,0.95))[1,],
  u90 = apply(abna.deaf_LOR, 2, FUN = quantile, c(0.05,0.95))[2,],
  pp = apply(abna.deaf_LOR, 2, FUN = function(x) sum(x>0)/length(x)) ),
  
  data.frame(
  LOR = apply(abna.hear_LOR, 2, median),
  mad = apply(abna.hear_LOR, 2, mad),
  l90 = apply(abna.hear_LOR, 2, FUN = quantile, c(0.05,0.95))[1,],
  u90 = apply(abna.hear_LOR, 2, FUN = quantile, c(0.05,0.95))[2,],
  pp = apply(abna.hear_LOR, 2, FUN = function(x) sum(x>0)/length(x)) )
  ), 2 )

groupLOR$group = rep(c("DEAF","HEARING"), each = length(grammar))

write.csv(groupLOR, "abna_groupLOR.csv")

#organize grammars for plot
abnadata<-rbind(abna.deafl,abna.hearl)

#capitalize ABNA
abnadata$grammar = as.character(abnadata$grammar)
abnadata$grammar[abnadata$grammar=="ABnA"] = "AB^NA"
grammar<-colnames(abna[,c(17:29)])
grammar[2]="AB^NA"
abnadata$grammar <- factor(abnadata$grammar, levels=grammar)

#reverse order of factor
abnadata$Group = relevel(as.factor(abnadata$Group), ref = "hearing")


#plot
library(tidybayes);library(cowplot)

abnaplot<-

ggplot(abnadata, aes(x=grammar, y =  prob, group = Group, color = Group))+
stat_pointinterval(.width = c(0.9), position=position_dodge(0.5), size=2, fatten_point=2)+
#scale_x_continuous(breaks=c(1,10,20,30))+
scale_y_continuous(limits=c(0.2,1), breaks = c(0.25, 0.5, 0.75, 1))+
scale_color_manual(values=c("#000000","#9c9898"), labels = c("Hearing", "Deaf"))+
ggtitle(expression(paste(AB^N,A))) +
ylab("Probability of consistent response\n")+
geom_hline(yintercept=0.5, linetype="dashed", size=0.5, color="#bab6b6")+
theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),unit="cm"),
        axis.ticks.y =element_blank(),
        axis.text.x =element_text(size=10, angle = 45, hjust = 1),
        axis.text.y =element_text(size=10),
        axis.title.x =element_blank(),
        axis.title.y =element_text(size=12, face="bold"),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size= 1, 
                                    linetype="solid"))+
  guides(color=FALSE)

save_plot("abna results.png", abnaplot,
          base_aspect_ratio=1.4,
          dpi=300)

######################################################################
#individual-level plot

id<-unique(abna[,c("ID", "Group") ])
id<-id[order(id$ID),]
df<-data.frame(ID=id[,1],Group=id[,2])

pred<-fitted(abnamod, scale = "response",newdata=df, robust=TRUE, prob=c(0.05,0.95))
predm<-data.frame(pred[1:30,1,1:13])
predl90<-data.frame(pred[1:30,3,1:13])
predu90<-data.frame(pred[1:30,4,1:13])

#subject info
predm$ind<-id$ID
predm$Group<-id$Group

#wide to long
indml<-gather(predm, "grammar", "re", ABFirst:SomeB, factor_key=TRUE)
indl90<-gather(predl90, "grammar", "l90", ABFirst:SomeB, factor_key=TRUE)
indu90<-gather(predu90, "grammar", "u90", ABFirst:SomeB, factor_key=TRUE)
indml$l90 = indl90$l90
indml$u90 = indu90$u90

#save results
write.csv(indml, "individual predictions abna.csv")

#individual level pp+
pred2<-fitted(abnamod, scale = "response", newdata=df, robust=TRUE, summary = FALSE)
pp_target =apply(pred2, 2:3, FUN = function(x) sum(x > 0.5)/length(x) )

#proportion of hearing participants with p > 0.70 at pp+ >= 0.95
sum(pp_target[1:15,"ABnA"] >= 0.95)/15

#proportion of deaf participants with p > 0.70 at pp+ >= 0.95
sum(pp_target[16:30,"ABnA"] >= 0.95)/15

#reverse order of factor
indml$Group = relevel(as.factor(indml$Group), ref = "HEARING")

#plot
abnaidplot<-
  
  ggplot(indml, aes(x=ind, y=re, group = Group, color = Group))+
  geom_point()+
  geom_errorbar(aes(ymin=l90,ymax=u90,width=0))+
  facet_wrap(~grammar, scales="fixed", nrow=3, ncol=5)+
  xlab("\nSubject ID")+
  ylab("Mean probability of consistent response\n")+
  geom_hline(yintercept=0.5, linetype="dashed", size=0.5, color="#bab6b6")+
  theme(plot.margin=unit(c(1,1,1,1),unit="cm"),
        axis.ticks.x =element_blank(),
        axis.ticks.y =element_blank(),
        axis.text.x =element_text(size=10),
        axis.text.y =element_text(size=10),
        axis.title.x =element_text(size=12, face="bold"),
        axis.title.y =element_text(size=12, face="bold"),
        legend.title=element_text(size=12,face="bold"),
        panel.border = element_rect(fill=NA,color="black", size= 1, 
                                    linetype="solid"))

save_plot("abna ind results.png", abnaidplot,
          base_width=12,base_height=7,
          dpi=300)


#just target grammar

abnaid.p<-
  
  ggplot(indml[indml$grammar=="ABnA",], aes(x=ind, y=re, group = Group, color = Group))+
  geom_point()+
  geom_errorbar(aes(ymin=l90,ymax=u90,width=0))+
  scale_color_manual(values=c("#000000","#9c9898"), labels = c("Hearing", "Deaf"))+
  scale_x_continuous(breaks=c(1,10,20,30))+
  scale_y_continuous(limits=c(0.2,1), breaks = c(0.25, 0.5, 0.75, 1))+
  xlab("\nSubject ID")+
  ggtitle(expression(paste(AB^N,A))) +
  ylab("Probability of consistent response\n")+
  geom_hline(yintercept=0.5, linetype="dashed", size=0.5, color="#bab6b6")+
  theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),unit="cm"),
        plot.title=element_text(hjust=0.5, face="bold"),
        axis.ticks.x =element_blank(),
        axis.ticks.y =element_blank(),
        axis.text.x =element_text(size=10),
        axis.text.y =element_text(size=10),
        axis.title.x =element_text(size=12, face="bold"),
        axis.title.y =element_text(size=12, face="bold"),
        legend.title=element_text(size=12,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size= 1, 
                                    linetype="solid"))+
    guides(color=FALSE)

######################################################################
######################################################################
######################################################################
######################################################################

#####################################################################
###Combined plot
####################################################################

#average probability estimates
fig5A<-
plot_grid(abnaplot, mirplot, copyplot, legend, ncol = 4, rel_widths = c(1,1,1,0.5))

save_plot("Fig 5A2.png", fig5A,
          base_width=14,base_height=4,
          dpi=600)

#individual probability estimates
fig5B<-
plot_grid(abnaid.p, mirid.p, copyid.p, legend,  ncol = 4, rel_widths = c(1,1,1,0.5))

save_plot("Fig 5B.png", fig5B,
          base_width=14,base_height=4,
          dpi=600)

fig5<-plot_grid(fig5A, NULL, fig5B, labels=c("A","","B"), ncol=1, rel_heights=c(1,0.05,1))
fig5

save_plot("Fig 5.png", fig5,
          base_width=14,base_height=8,
          dpi=600)

save_plot("Fig 5.tiff", fig5, compression="lzw",
          base_width=14,base_height=8,
          dpi=600)

######################################################################
######################################################################
######################################################################
######################################################################

#####################################################################
###ABnA stimulus property model
####################################################################

#data
ABNAtrial<-read.csv("ABnA trial type.csv")

#build model
ABnA.m <- bf(ABnA ~ Group*General + Group*Similar + (General + Similar|ID)) + bernoulli()

#run model
abnaint<- brm(formula=ABnA.m,
              data=ABNAtrial,  
              prior=c(prior("normal(0,2)",class="Intercept"),
                prior("normal(0,2)",class="b"),
                prior("cauchy(0,2)",class="sd"),
                prior("lkj(2)", class="cor")),
              warmup=1000,iter=3000, chains=4, seed=9,
              control=list(adapt_delta=0.99))

#save
saveRDS(abnaint,"abna_int.RDS")

#summary
abnafe<-fixef(abnamain,robust=TRUE,probs=c(0.05,0.95))
abnafe

#logit sd
logitsd<-pi/sqrt(3)

#cohens d
abnafe[,1]/logitsd

#pp
c(
  hypothesis(abnaint, "GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(abnaint, "GeneralRecognition > 0")$hypothesis$Post.Prob,
  hypothesis(abnaint, "SimilarSimilar > 0")$hypothesis$Post.Prob,
  hypothesis(abnaint, "GroupHEARING:GeneralRecognition > 0")$hypothesis$Post.Prob,
  hypothesis(abnaint, "GroupHEARING:SimilarSimilar  > 0")$hypothesis$Post.Prob)

#predict values
newdata<-data.frame(Group=c("HEARING","HEARING", "HEARING", "HEARING", 
                            "DEAF", "DEAF", "DEAF", "DEAF"),
                    General=c("Recognition","Generalization","Generalization","Recognition",
                              "Recognition","Generalization","Generalization","Recognition"),
                    Similar=c("Similar","Dissimilar","Similar","Dissimilar",
                              "Similar","Dissimilar","Similar","Dissimilar"))
  
predval = fitted(abnaint,newdata=newdata,scale="linear",probs=c(0.05,0.95),robust=TRUE,re_formula=NA, summary = FALSE)

#columns
newdata

#Hearing - Deaf for recognition , similar
c1 = predval[,1] - predval[,5]
#Hearing - Deaf for recognition , dissimilar
c2 = predval[,4] - predval[,8]
#Hearing - Deaf for generalization , similar
c3 = predval[,3] - predval[,7]
#Hearing - Deaf for generalization , dissimilar
c4 = predval[,2] - predval[,6]

apply(cbind(c1,c2,c3,c4), 2, median)
apply(cbind(c1,c2,c3,c4), 2, mad)
apply(cbind(c1,c2,c3,c4), 2, quantile, c(0.05,0.95))
apply(cbind(c1,c2,c3,c4), 2, FUN = function(x) sum(x>0)/length(x))
apply(cbind(c1,c2,c3,c4), 2, FUN = function(x) median(x/logitsd))


#####################################################################
###Copy stimulus type model
####################################################################

#data
Copytrial<-read.csv("Copy trial type.csv")
  
#build model
Copy.m <- bf(Copy ~ Group*General + Group*Similar + (General + Similar|ID)) + bernoulli()

#run model
copyint<- brm(formula=Copy.m,
              data=Copytrial,  
              prior=c(prior("normal(0,2)",class="Intercept"),
                prior("normal(0,2)",class="b"),
                prior("cauchy(0,2)",class="sd"),
                prior("lkj(2)", class="cor")),
              warmup=1000,iter=3000, chains=4, seed=9,
              control=list(adapt_delta=0.99))

#save
saveRDS(copyint,"copy_int.RDS")

#summary
copyfe<-fixef(copyint,robust=TRUE,probs=c(0.05,0.95))
copyfe

#logit sd
logitsd<-pi/sqrt(3)

#cohens d
copyfe[,1]/logitsd

#pp
c(
  hypothesis(copyint, "GroupHEARING > 0")$hypothesis$Post.Prob,
  hypothesis(copyint, "GeneralRecognition > 0")$hypothesis$Post.Prob,
  hypothesis(copyint, "SimilarSimilar > 0")$hypothesis$Post.Prob,
  hypothesis(copyint, "GroupHEARING:GeneralRecognition > 0")$hypothesis$Post.Prob,
  hypothesis(copyint, "GroupHEARING:SimilarSimilar  > 0")$hypothesis$Post.Prob)

#predict values
newdata<-data.frame(Group=c("HEARING","HEARING", "HEARING", "HEARING", 
                            "DEAF", "DEAF", "DEAF", "DEAF"),
                    General=c("Recognition","Generalization","Generalization","Recognition",
                              "Recognition","Generalization","Generalization","Recognition"),
                    Similar=c("Similar","Dissimilar","Similar","Dissimilar",
                              "Similar","Dissimilar","Similar","Dissimilar"))
  
predval = fitted(copyint,newdata=newdata,scale="linear",probs=c(0.05,0.95),robust=TRUE,re_formula=NA, summary = FALSE)

#columns
newdata

#Hearing - Deaf for recognition , similar
c1 = predval[,1] - predval[,5]
#Hearing - Deaf for recognition , dissimilar
c2 = predval[,4] - predval[,8]
#Hearing - Deaf for generalization , similar
c3 = predval[,3] - predval[,7]
#Hearing - Deaf for generalization , dissimilar
c4 = predval[,2] - predval[,6]

apply(cbind(c1,c2,c3,c4), 2, median)
apply(cbind(c1,c2,c3,c4), 2, mad)
apply(cbind(c1,c2,c3,c4), 2, quantile, c(0.05,0.95))
apply(cbind(c1,c2,c3,c4), 2, FUN = function(x) sum(x>0)/length(x))
apply(cbind(c1,c2,c3,c4), 2, FUN = function(x) median(x/logitsd))
