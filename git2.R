library(adehabitatLT)
library(lubridate)
library(dplyr)
library(lme4)
library(MuMIn)
library(emmeans)
mydata<-readRDS("AreaDF3.rds")
head(mydata)

car<-as.ltraj(xy=mydata[,c("Easting","Northing")], date=mydata$Timestamp, id=mydata$Individual) #make into a burst
refda<-strptime("00:00:00", "%H:%M:%S", tz="EST") #set a reference date for filling in unsampled hours
car2<-setNA(car, refda, 2, units="hour") #Add an NA for the hours not sampled
car3<-ld(car2) #make into a df
car3$hour<-hour(car3$date) # add a column with the hour
car3<-car3[!(car3$hour==10),] #deleting hours that weren't sampled
car3<-car3[!(car3$hour==14),]
car3<-car3[!(car3$hour==16),]
car3<-car3[!(car3$hour==8),]
nrow(car3)
car5<-as.ltraj(xy=car3[,c("x","y")], date=car3$date, id=car3$id)
foo <- function(date) {
  da <- as.POSIXlt(date)
  ho <- da$hour + da$min
  return(ho>11&ho<13) #this is the time when every burst stops
}

deer<-cutltraj(car5, "foo(date)", nextr=TRUE)
dfdeer<-ld(deer)
dfdeer
dd<-data.frame(dfdeer)
dd

mean(dd$dist, na.rm=TRUE)
ds<-dd %>% 
  group_by(burst) %>%
  summarise(diff=last(date)-first(date)) #column with the duration of every burst

d7<-aggregate(dd$dist, by=list(burst=dd$burst), FUN=sum, na.rm=TRUE) #df with total length (distance) of every burst
dg<-merge.data.frame(dd, ds, by.x="burst", by.y="burst" )
dr<-merge.data.frame(dg, d7, by.x="burst", by.y="burst")
dr$diff<-as.numeric(dr$diff)
dw<-dr[(dr$diff>=16),] #only keeping bursts 12 hr or longer
dw<-dw[complete.cases(dw$R2n),] #deletes the NAs
d2<-split(dw, dw$burst)
df1<-lapply(d2, tail, 1) #keeping the last row of every element because the last R2n is needed
bi<-bind_rows(df1)
mdf<-merge.data.frame(mydata, bi, by.x="Individual", by.y="id" )
mdf<-subset(mdf, !duplicated(burst)) #delete duplicates
mdf<-mdf[(mdf$R2n>0),] #delete the zeros
mdf<-mdf[complete.cases(mdf$R2n),] #deletes the NAs
names(mdf)[names(mdf) == "x.y"] <- "distmoved"
mdf<-mdf[(mdf$distmoved>0),] #delete the zeros
mdf$Group <- factor(mdf$Group, levels=c("Control", "Transient", "Resident"))
names(mdf)[names(mdf) == "SEX"] <- "Sex"
mdf$Sex[mdf$Sex == "M"] <- "Male"
mdf$Sex[mdf$Sex == "F"] <- "Female"
mdf2<-mdf

mdf2$sqd<-sqrt(mdf2$R2n) #This value is squared so took the square root
#linear mixed model with animal ID as a random effect
displ<-lmer(log(sqd)~(Group+Sex+Start+End)^2+(1|collarID),data=mdf2, na.action = na.pass, REML=F)
dredge(displ)
dfit<-lmer(log(sqd)~Group+Sex+Start+End+End:Group+Group:Sex+Group:Start+(1|collarID),data=mdf2, na.action = na.pass, REML=F)
r.squaredGLMM(dfit)
emmeans(dfit, pairwise~Start|Group, type="response")
emmeans(dfit, pairwise~SEX|Group, type="response")
emmeans(dfit, pairwise~End|Group, type="response")

#Nightly displacement by group and habitat
emmip(dfit,  ~Start|Group, CIs=TRUE, CIarg=aes(lwd=1, color="black", linetype=1), type="response")+
  facet_wrap(~Group)+scale_linetype_manual(name = "Start", values=c("blank","blank"))+aes(linetype=Start)+
  theme(strip.background=element_rect(color="black", fill="darkslategray3"), panel.border = element_rect(colour='black', fill='NA'), 
        panel.background = element_rect(fill='gray97'), panel.grid = element_line(colour = NA), panel.spacing = unit(1, "lines"),
        legend.key = element_rect(fill='white'),
        axis.text=element_text(size=12), axis.title.x = element_blank(), legend.position="none",
        axis.title.y = element_text(size=12, vjust=3.9), strip.text = element_text(size=12, face="bold"))+
  theme(plot.margin=unit(c(0.2,0.2,0.2,0.5), "cm"))+ylab(expression(atop("Nightly", paste("displacement (km)"))))+scale_y_continuous(labels=function(x) format(x,scientific=FALSE))+scale_x_discrete(name='Start', labels=c('Bottomland
hardwood','Upland
Pine'))

#Nightly displacement by group and sex
emmip(dfit, ~Sex|Group, CIs=TRUE, CIarg=aes(lwd=1, color="black", linetype=1), type="response")+
  facet_wrap(~Group)+scale_linetype_manual(name = "Sex", values=c("blank","blank"))+aes(linetype=Sex)+
  theme(strip.background=element_rect(color="black", fill="darkslategray3"), panel.border = element_rect(colour='black', fill='NA'), 
        panel.background = element_rect(fill='gray97'), panel.grid = element_line(colour = NA), panel.spacing = unit(1, "lines"),
        legend.key = element_rect(fill='white'), legend.position="none", 
        axis.text=element_text(size=12), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=12, vjust=3.9), strip.text = element_text(size=12, face="bold"))+
  theme(plot.margin=unit(c(0.2,0.2,0.2,0.5), "cm"))+
  ylab(expression(atop("Nightly", paste("displacement (m)"))))+
  scale_y_continuous(labels=function(x) format(x,scientific=FALSE))

