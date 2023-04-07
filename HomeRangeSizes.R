rm(list=ls(all=TRUE))
library(ggeffects); library(sjPlot)
library(segclust2d); library(move); library(openair)
library(dplyr); library(lme4); library(MuMIn)
library(ggplot2);library(emmeans)
library(ggpubr); library(ggspatial)
set.seed(34)


##Read in control animals###
controls<-read.csv("precontrolcollar.csv")
unique(controls$collarID)
controls$Timestamp<-as.POSIXct(controls$Timestamp, tz="EST", 
	format="%m/%d/%Y %H:%M:%S") #convert time to POSIX#
###Sampling set period of time (to match transient/resident time periods)###
controls$date<-controls$Timestamp 
	#selectByDate function requires column named 'date' (lowercase!!) #
controlwint<-selectByDate(controls, 
	month = c("December","January","February", 
			"November","October","March"))
newdata <- controlwint[order(controlwint$collarID, controlwint$Timestamp),]
	 #ordering the dataframe#
contlist<-split(newdata, newdata$collarID) 
	#split into list by collarID#
monthlist<-list() #make a list for the loop#
names(contlist)
for (i in 1:length(contlist)) {
  monthlist[[i]]<-cut(contlist[[i]]$Timestamp, breaks="45 days")
  
} #break up each animal into 45 days chunks"


timebreaks<-data.frame(unlist(monthlist)) #unlist and put into dataframe#
timebreaks$rown <- 1:nrow(timebreaks) #add column with row number#
newdata$rown<-1:nrow(newdata) #add column with row number#
dataset<-merge.data.frame(newdata, timebreaks, 
	by.y="rown", by.x="rown") 
	#merge dataset with time breaks based 
	#on row numbers (won't work if both not in same order)#
dataset$ID_Timebreak<-paste(dataset$collarID, 
	dataset$`unlist.monthlist.`) 
	#add a column with every time break for every animal#
duration<-dataset %>% group_by(ID_Timebreak) %>% 
	summarise(firstvis=min(Timestamp), 
	lastvis=max(Timestamp), duration=lastvis-firstvis) 
	#tibble with the duration of every time break (not all div by 45)#
complete<-merge.data.frame(duration, dataset, 
	by.y="ID_Timebreak", by.x="ID_Timebreak") 
	#merge duration with the dataframe#
complete %>% filter(duration > 44) #filter out ones less than 45#
## complete[5000:6000,]
by_cyl <- complete %>% group_by(collarID) 
	#grouping by ID for time selection#
r77<-sample_n(by_cyl, 1, replace = FALSE) 
	#sampling one time break for every animal#
controlfinal<-complete %>%
  filter(ID_Timebreak %in% r77$ID_Timebreak) 
	#filter the dataframe based on the breaks selected#
## controlfinal[5000:5500,] #dataframe merging can get wonky 
	# so spot check to make sure everything's lined up"
coordinates(controlfinal)<-c("Easting","Northing")
proj4string(controlfinal)<-CRS("+proj=longlat +datum=WGS84")
control2<-as.data.frame(spTransform(controlfinal, 
	CRS('+proj=utm +zone=17N +datum=NAD83'))) #convert to UTM#
control2$Timestamp<-as.POSIXct(control2$Timestamp, 
	tz="EST", format="%m/%d/%Y %H:%M:%S")

###Read in translocated and perform segmentation###

moved<-read.csv("collarstranslocation.csv")
moved<- moved %>% distinct() #removing duplicates#
coordinates(moved)<-c("Easting","Northing")
proj4string(moved)<-CRS("+proj=longlat +datum=WGS84")
moveddf<-as.data.frame(spTransform(moved, 
	CRS('+proj=utm +zone=17N +datum=NAD83')))
moveddf$Timestamp<-as.POSIXct(moveddf$Timestamp, 
	tz="EST", format="%m/%d/%Y %H:%M:%S")
movedlist<-split(moveddf, moveddf$ID)
unique(moved$collarID)
finalmoved<-movedlist

finalmoved<-movedlist[names(movedlist) %in% 
                        c("40996 PostMove", "41011 PostMove", "41002 PostMove") == FALSE] 
		#dropping data deficient ones#

shift.list<-list() #create list for the loop#
shift.info<-list() #create list for the loop#
for (i in seq_along(finalmoved)) {
  shift.list[[i]]<-segmentation(finalmoved[[i]],lmin=25,  
		Kmax = 7, seg.var = c("Easting","Northing"), 
		scale.variable=FALSE) 
  shift.info[[i]]<-augment(shift.list[[i]])
  
}

#sidequest- calculating durations of the different behavioral states
segresults <- do.call("rbind", shift.info) 
	#merge list elements back into a single dataframe#
segresults$ID_state<-paste(segresults$collarID, 
	segresults$state) #create column with ID for 
				#every state for every animal#
segduration<-segresults %>% group_by(ID_state) %>% 
  summarise(firstvis=min(Timestamp), 
            lastvis=max(Timestamp), segduration=lastvis-firstvis) 
#This creates a tibble with the duration of every state#
segresultclean<-subset(segresults, !duplicated(ID_state)) 
#keep one row per ID_state to tidy things up#
segresultfinal<-merge.data.frame(segduration, 
                                 segresultclean, by.y="ID_state", by.x="ID_state")

#back to the original task#
trans_res<-segresults %>% filter(state <= 2 ) #only keeping 1st two states#
transkeep=trans_res[c("collarID","SEX","state","Treatment", 
	"PrePost","Timestamp","Easting","Northing","ID_state")]  
	#keep only the needed columns#

control2$state<-0 #Assigning controls state of 0 to 
			#match state of 1 for trans and 2 for res#
control2$ID_state<-paste(control2$collarID, control2$state)

controlkeep=control2[c("collarID","SEX","state","Treatment", 
	"PrePost","Timestamp","Easting","Northing","ID_state")] 
	#make same columns as trans for binding"

finaldata<-rbind(controlkeep, transkeep) #bind the two dataframes"

finaldata$Group[finaldata$state==0]<-"Control" 
	#Making new column and assigning control/trans/res 
	#based on state number#
finaldata$Group[finaldata$state==1]<-"Transient" 
finaldata$Group[finaldata$state==2]<-"Resident"
## head(finaldata)
final2<-finaldata[with(finaldata, order(collarID, Timestamp)),] 
	#need ascending timestamps for every individual"
final2<-final2[!(final2$collarID=="41015"),]
final2<-final2[!(final2$collarID=="41002"),]
final2<-final2[!(final2$collarID=="40996"),]
unique(final2$collarID)
###Home range construction###
raccoons<-move(x=final2$Easting, y=final2$Northing, 
	animal=final2$ID_state, data=final2, time=final2$Timestamp)
datalist<-split(raccoons, final2$ID_state) ## names(datalist)

modlist95<-list()
udlist95<-list()
plotlist95<-list()
arealist95<-list()
modlist60<-list()
udlist60<-list()
plotlist60<-list()
arealist60<-list()

for (i in 1:length(datalist) ){ ## i=names(datalist)[2]
  modlist95[[i]]<-brownian.bridge.dyn(datalist[[i]], location.error=7, raster=50, ext=4, margin=5, window=11)
  flush.console()
  udlist95[[i]]<-getVolumeUD(modlist95[[i]])
  plotlist95[[i]]<-udlist95[[i]]<=0.95
  arealist95[[i]]<-sum(values(plotlist95[[i]]))
  
 
}

#Ouput of arealist is number of raster cells, so
	# multiply area by the size of the raster cells then 
	# divide by 1e06 to get sq km#
Area<-lapply(arealist95, function(x) (x*2500)/1000000) 
Area2<-as.data.frame(unlist(Area)) #get the areas into a df#
Individual<-as.data.frame(names(datalist))
AreaDF<-merge.data.frame(data.frame(cbind(Area2, Individual)), 
	finaldata, by.x="names.datalist.", by.y="ID_state" )

modlist60<-list()
udlist60<-list()
plotlist60<-list()
arealist60<-list()

for (i in names(datalist) ){ ## i=1
  modlist60[[i]]<-brownian.bridge.dyn(datalist[[i]], 
	location.error=7, raster=50, ext=4, margin=5, window=13,verbose=F)
  flush.console()
  udlist60[[i]]<-getVolumeUD(modlist60[[i]])
  plotlist60[[i]]<-udlist60[[i]]<=0.60
  arealist60[[i]]<-sum(values(plotlist60[[i]]))
}

Area60<-lapply(arealist60, function(x) (x*2500)/1000000) 
Area260<-as.data.frame(unlist(Area60))
AreaDF2<-merge.data.frame(data.frame(cbind(Area260, Individual)), 
	AreaDF, by.x="names.datalist.", by.y="names.datalist." )
AreaDF3<-AreaDF2
names(AreaDF3)[names(AreaDF3) == "unlist.Area."] <- "UD95"
names(AreaDF3)[names(AreaDF3) == "unlist.Area60."] <- "UD60"
names(AreaDF3)[names(AreaDF3) == "names.datalist."] <- "Individual"
AreaDF3<-AreaDF3[!(AreaDF3$Individual=="40997 2"),]
AreaDF3<-AreaDF3[!(AreaDF3$Individual=="41013 2"),]
AreaDF3=AreaDF3[order(AreaDF3$Individual,AreaDF3$Timestamp),]

#Assigning everyone high and low density treatments for starting and ending habitat
AreaDF3<- AreaDF3 %>% mutate(Start =
	case_when(
		Treatment=="H2L"  ~ "High",
		Treatment=="H2H"  ~ "High", 
		Treatment=="L2L"  ~ "Low",
		Treatment=="L2H"  ~ "Low",
		collarID=="40993"~ "Low",
		collarID=="41001"~ "Low",
		collarID=="41010"~ "Low",
		collarID=="41012"~ "High",
		collarID=="41015"~ "High",
		collarID=="50000"~ "Low",
		collarID=="70000"~ "Low"))


AreaDF3<- AreaDF3 %>% mutate(End =
	case_when(
		Treatment=="H2L"  ~ "Low",
		Treatment=="H2H"  ~ "Low", 
		Treatment=="L2L"  ~ "Low",
		Treatment=="L2H"  ~ "High",
		collarID=="40993"~ "Low",
		collarID=="41001"~ "Low",
		collarID=="41010"~ "Low",
		collarID=="41012"~ "High",
		collarID=="41015"~ "High",
		collarID=="50000"~ "Low",
		collarID=="70000"~ "Low"))

names(AreaDF3)[names(AreaDF3) == "SEX"] <- "Sex"
names(AreaDF3)[names(AreaDF3) == "Group"] <- "State"
names(AreaDF3)[names(AreaDF3) == "state"] <- "numberstate"
AreaDF3$Sex[AreaDF3$Sex == "M"] <- "Male"
AreaDF3$Sex[AreaDF3$Sex == "F"] <- "Female"
AreaInput2<-subset(AreaDF3, !duplicated(Individual)) 
	#keep one row per ID_state to tidy things up#

AreaInput2 %>%  group_by(State) %>% 
  tally() %>% print(n=27)
nrow(AreaInput2)

#calculating the cumulative duration of the control period in high density habitat
AreaDF3a<-AreaDF3[(AreaDF3$PrePost=="PreMove"),]
AreaDF3a<-AreaDF3[(AreaDF3$Start=="High"),]
duration<-AreaDF3a %>% group_by(Individual) %>% 
  summarise(firstvis=min(Timestamp), 
            lastvis=max(Timestamp), duration=lastvis-firstvis) 
sum(duration$duration)
AreaInput2$Group <- factor(AreaInput2$State, levels=c("Control", "Transient", "Resident"))

#Running linear models with size of UDs as the response
UD95full<-lm(log(UD95)~(Group+Sex+Start+End)^2, data=AreaInput2, na.action = na.pass)
UD60full<-lm(log(UD60)~(Group+Sex+Start+End)^2, data=AreaInput2, na.action = na.pass)
dredge(UD95full)
dredge(UD60full)

#Predictions and interactions from the top models
UD95fit<-lm(log(UD95)~Group+Sex+Start, data=AreaInput2, na.action = na.pass)
ggpredict(UD95fit, ci.lvl=0.95, terms=c("Start"), condition = c(Sex="Male", Group="Resident"))
ggpredict(UD95fit, ci.lvl=0.95, terms=c("Sex"), condition = c(Start="Low", Group="Resident"))
ggpredict(UD95fit, ci.lvl=0.95, terms=c("Group"))
emmeans(UD95fit, pairwise~Group, type="response")

UD60fit<-lm(log(UD60)~Group+Sex+Start+Group:Start, data=AreaInput2, na.action = na.pass)
ggpredict(UD60fit, ci.lvl=0.95, terms=c("Sex"), condition = c(Start="Low", Group="Resident"))
ggpredict(UD60fit, ci.lvl=0.95, terms=c("Group"), condition = c(Sex="Male", Start="Low"))
ggpredict(UD60fit, ci.lvl=0.95, terms=c("Group"), condition = c(Sex="Male", Start="High"))
emmeans(UD60fit, pairwise~Start|Group, type="response")
emmip(UD60fit, ~Start|Group, CIs=TRUE, CIarg=aes(lwd=1, color="black", linetype=1), type="response")

#Table of summary stats
STAT<-AreaInput2 %>% 
  group_by(SEX, Group, Start) %>% 
  summarise(Mean95=mean(UD95),
            sd95=sd(UD95),
            mean60=mean(UD60),
            sd60=sd(UD60),
             )
STAT


#Figure for the 95% UD
plot_model(UD95fit, type="pred", terms = c("Group", "Sex", "Start"), colors = c("chocolate2","darkolivegreen"), dot.size = 3)+
  theme(strip.background=element_rect(color="black", fill="darkslategray3"), plot.title = element_blank(),
        panel.border = element_rect(colour='black', fill='NA'), 
        panel.background = element_rect(fill='gray96'), panel.grid = element_line(colour = NA), panel.spacing = unit(2, "lines"),
        legend.key = element_rect(fill='white'), legend.position="bottom", legend.text = element_text(size=12),
        axis.text=element_text(size=12), axis.title.x = element_text(size=14, vjust=-3),
        axis.title.y = element_text(size=14, vjust=3.9), strip.text = element_text(size=12, face="bold"))+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.9), "cm"))+xlab(NULL)+ylab("95% UD (kilometers squared)")

#Figure for the 60% UD
plot_model(UD60fit, type="pred", terms = c("Group", "Sex", "Start"), colors = c("chocolate2","darkolivegreen"), dot.size = 3)+
  theme(strip.background=element_rect(color="black", fill="darkslategray3"), plot.title = element_blank(),
        panel.border = element_rect(colour='black', fill='NA'), 
        panel.background = element_rect(fill='gray96'), panel.grid = element_line(colour = NA), panel.spacing = unit(2, "lines"),
        legend.key = element_rect(fill='white'), legend.position="bottom", legend.text = element_text(size=12),
        axis.text=element_text(size=12), axis.title.x = element_text(size=14, vjust=-3),
        axis.title.y = element_text(size=14, vjust=3.9), strip.text = element_text(size=12, face="bold"))+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.9), "cm"))+xlab(NULL)+ylab("60% UD (kilometers squared)")

#Figure showing segmentation
r7<-segmap(shift.list[[11]], coord.names=c("Easting", "Northing"))+geom_path(size=0.5)+
  ylab("Northing (m)")+xlab("Easting (m)")+annotation_scale(location="bl", width_hint=0.23, 
                                                            height=unit(0.3, "cm"), pad_x = unit(0.3, "cm"), pad_y=unit(0.5, "cm"),
                                                            text_pad = unit(0.25, "cm"), text_cex=1, tick_height = 1, style="bar", bar_cols = "white")+
  scale_color_manual("State", values = c('coral2', 'turquoise4', "red"),
                     labels=c("Transient", "Resident", "3"))+ theme(plot.margin=unit(c(0.5,0.5,.5,.5),"cm"), panel.border = element_rect(colour='black', fill='NA'), 
                                                                    panel.background = element_rect(fill='gray97'), panel.grid = element_line(colour = NA), 
                                                                    legend.key = element_rect(fill='white'), legend.position="bottom",
                                                                    axis.text=element_text(size=12), axis.title.x = element_text(size=14, vjust=1),
                                                                    axis.title.y = element_text(size=14, vjust=1))+
  geom_point(aes(colour=factor(state)), size=1)+xlim(420000, 450000)

r8<-plot(shift.list[[11]], xcol="Timestamp")+
  theme(plot.margin=unit(c(0.5,0.5,.5,.5),"cm"), axis.text=element_text(size=12), axis.title.x = element_text(size=14, vjust=1),
        axis.title.y = element_text(size=14, vjust=1), legend.position = "bottom", strip.text.x = element_text(size=13))+
  ylab("Location")+xlab("Month")+scale_fill_manual("State", values = c('coral2', 'turquoise4'),
     labels=c("Transient", "Resident"))+scale_color_manual(values=c('coral2', 'turquoise4'),  
name="State", labels=c("Transient", "Resident"))

ggarrange(r7, r8, labels = c("A", "B"), hjust=-4)
