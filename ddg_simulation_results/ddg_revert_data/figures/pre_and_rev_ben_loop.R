library(ggplot2)
library(cowplot)
library(gam)
require(mgcv)
require(base)
require(graphics)
require(plotrix)
library(nls2)
library(segmented)
#draws probability of accepting reversion and preadaption over -15 and 15 markov steps for both ben and delt subs


#color scale to match ggplots default


ggplotColours <- function(n = 4, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
ggcols=ggplotColours(n = 4)
cbbPalette <- c(ggcols[1],ggcols[2], ggcols[3],ggcols[4])
cols <- c("Flexible Backbone & Sidechain"=cbbPalette[1],"Fixed Backbone & Flexible Sidechain"=cbbPalette[2],"Fixed Backbone & Sidechain"=cbbPalette[3],"Flexible Backbone & Fixed Sidechain"=cbbPalette[4])



#gets data and figures out hot and cold spots
getdata<-function(where){

#read in data 
TT.files<-list.files(where,pattern="revert")


TT_by_15_high<-NULL
TT_by_15_low<-NULL
for( i in 1:length(TT.files)){
  name<-paste(where,TT.files[i],sep='')
  TT<-read.csv(name,head=TRUE,sep=",")

  prob<-TT$Prob
  ben<-TT$Prob_target_fixed
  high_prob<-NULL
  low_prob<-NULL


  for( j in 1:(length(prob))){  
     check_ben<-ben[j]
	#if mut wasnt ben when fixed
      if(check_ben < 1){
		low_prob=c(low_prob,prob[j])
		}
	#if mut was ben when it fixed
   	if(check_ben == 1){
		high_prob=c(high_prob,prob[j])
	

    }

   #stick everything in a matrix, dims are based on knowing how many insertion evaluations where done. 
   if(length(high_prob)!=0){
  	mat<-matrix(high_prob,nrow=31)
  	TT_by_15_high<-cbind(TT_by_15_high,mat)
    }

   if(length(low_prob)!=0){
  	mat<-matrix(low_prob,nrow=31)
  	TT_by_15_low<-cbind(TT_by_15_low,mat)
    }

}


#make data frames out of both the ben and the del data 
dat_TT_high<-rowMeans(TT_by_15_high,na.rm = TRUE)
dat_TT_high<-data.frame(x=-15:15,
y=as.vector(as.numeric(dat_TT_high)),
all_mean=as.vector(as.numeric(dat_TT_high)))


dat_TT_low<-rowMeans(TT_by_15_low,na.rm = TRUE)

dat_TT_low<-data.frame(x=-15:15,
y=as.vector(as.numeric(dat_TT_low)),
all_mean=as.vector(as.numeric(dat_TT_low)))


info<-list("high"=dat_TT_high, "low"=dat_TT_low)
return(info)
}


#fits functions and returns the parameters values
getfit<-function(dat,type){

	if(type == "rev"){
		TT_fit_high<-nls2(y~I(c*exp(a*x)+b),data=subset(dat,x>=0),start=list(a=-.5,b=.5,c=.3))
		a<-as.numeric(summary(TT_fit_high)$parameters[1,1])
		b<-as.numeric(summary(TT_fit_high)$parameters[2,1])
		c<-as.numeric(summary(TT_fit_high)$parameters[3,1])
		info<-list("a_r"=a, "b_r"=b, "c_r"=c)
		return(info)
	}

	if(type == "pre"){
		TT_fit_high<-nls2(y~I(c*exp(a*x)+b),data=subset(dat,x<0),start=list(a=.5,b=.5,c=.3))
		a<-as.numeric(summary(TT_fit_high)$parameters[1,1])
		b<-as.numeric(summary(TT_fit_high)$parameters[2,1])
		c<-as.numeric(summary(TT_fit_high)$parameters[3,1])
		info<-list("a_p"=a, "b_p"=b, "c_p"=c)
		return(info)
	}


}


#get TT data
tt_data<-getdata('./bb_T_ch_T/')


#make fits
high_params_rev<-getfit(dat=tt_data$high, type="rev")
high_params_pre<-getfit(dat=tt_data$high, type="pre")


#stick everything in data frame for ggplot
dat_TT_high<-data.frame(tt_data$high,high_params_rev,high_params_pre)


#data from fits, i pasted this in here
print(high_params_rev)
print(high_params_pre)
#$a_r
#[1] -0.2197254
#
#$b_r
#[1] 0.4976108

#$c_r
#[1] 0.1510452

#> print(high_params_pre)
#$a_p
#[1] 0.116294

#$b_p
#[1] 0.798284

#$c_p
#[1] 0.1512913


#set up plotting evn for ben when fixed
lp<-ggplot(data=dat_TT_high, aes(x=x, y=y))
lp<-lp+labs(x="Markov Step",y="Probability of Accepting Mutation")
lp<-lp+guides(color=guide_legend(override.aes=list(fill=NA)))
lp<-lp+scale_x_continuous(breaks=seq(-15,15,5))

lp<-lp+theme(legend.key.size = unit(1, 'lines'), legend.key = element_rect(fill = "transparent", colour = "transparent"))+scale_colour_manual(name="",values=cols) +theme(legend.justification= c(1, 1), legend.position = c(1, .45))

lphot<-lp+ggtitle("Beneficial When Fixed")
lphot<-lphot+coord_cartesian(ylim=c(.45,1))+theme(legend.position="none")
lphot<-lphot+geom_point(aes(x=x,y=all_mean,colour="Flexible Backbone & Sidechain"))
lphot<-lphot+geom_smooth(aes(x=x,y=y,colour="Flexible Backbone & Sidechain"),subset(dat_TT_high,x>=0),method=gam,formula=y~I(0.1510452*exp(-0.2197254
*x)+0.4976108))
lphot<-lphot+geom_smooth(aes(x=x,y=y,colour="Flexible Backbone & Sidechain"),subset(dat_TT_high,x<0),method=gam,formula=y~I(0.1512913*exp( 0.116294*x)+ 0.798284)) 



#fit for TT for delt when fixed
low_params_rev<-getfit(dat=tt_data$low, type="rev")
low_params_pre<-getfit(dat=tt_data$low, type="pre")

print(low_params_rev)
print(low_params_pre)
#$a_r
#[1] -0.1127801

#$b_r
#[1] 0.7839215

#$c_r
#[1] 0.1889705

#> print(low_params_pre)
#$a_p
#[1] 0.467777

#$b_p
#[1] 0.5424723

#$c_p
#[1] 0.1049348


#> print(low_params_rev)

#stick everything in data frame to make ggplot happy
dat_TT_low<-data.frame(tt_data$low,low_params_rev,low_params_pre)

lpcold<-ggplot(data=dat_TT_low, aes(x=x, y=y))
lpcold<-lpcold+labs(x="Markov Step",y="Probability of Accepting Mutation")
lpcold<-lpcold+guides(color=guide_legend(override.aes=list(fill=NA)))
lpcold<-lpcold+scale_x_continuous(breaks=seq(-15,15,5))

lpcold<-lpcold+theme(legend.key.size = unit(1, 'lines'), legend.key = element_rect(fill = "transparent", colour = "transparent"))+scale_colour_manual(name="",values=cols) +theme(legend.justification= c(1, 1), legend.position = c(1, .45))
lpcold<-lpcold+coord_cartesian(ylim=c(.45,1))
lpcold<-lpcold+ggtitle("Deleterious When Fixed")

lpcold<-lpcold+geom_point(aes(x=x,y=all_mean,colour="Flexible Backbone & Sidechain"))
lpcold<-lpcold+geom_smooth(aes(x=x,y=y,colour="Flexible Backbone & Sidechain"),subset(dat_TT_low,x>=0),method=gam,formula=y~I(0.1889705*exp(-0.1127801  *x)+0.7839215)) 
lpcold<-lpcold+geom_smooth(aes(x=x,y=y,colour="Flexible Backbone & Sidechain"),subset(dat_TT_low,x<0),method=gam,formula=y~I(0.1049348*exp( 0.467777*x)+0.5424723 )) 

###EVERYTHING BELOW THIS IS BASICALLY CUT AND PASTE FOR THE OTHER DATA SETS, FT,TF, and FF

####FT CASE

ft_data<-getdata('./bb_F_ch_T/')

high_params_rev<-getfit(dat=ft_data$high, type="rev")
high_params_pre<-getfit(dat=ft_data$high, type="pre")
print(high_params_rev)
print(high_params_pre)
#$a_r
#[1] -0.1778823

#$b_r
#[1] 0.46733

#$c_r
#[1] 0.1813134

#> print(high_params_pre)
#$a_p
#[1] 0.09206659

#$b_p
#[1] 0.8115055

#$c_p
#[1] 0.1521805


dat_FT_high<-data.frame(ft_data$high,high_params_rev,high_params_pre)


lphot<-lphot+geom_point(aes(x=x,y=all_mean,colour="Fixed Backbone & Flexible Sidechain"),dat_FT_high)
lphot<-lphot+geom_smooth(aes(x=x,y=y,colour="Fixed Backbone & Flexible Sidechain"),subset(dat_FT_high,x>=0),method=gam,formula=y~I(0.1813134*exp(-0.1778823*x)+0.46733))
lphot<-lphot+geom_smooth(aes(x=x,y=y,colour="Fixed Backbone & Flexible Sidechain"),subset(dat_FT_high,x<0),method=gam,formula=y~I(0.1521805*exp( 0.09206659*x)+0.8115055))

high_params_rev<-getfit(dat=ft_data$low, type="rev")
high_params_pre<-getfit(dat=ft_data$low, type="pre")
print(high_params_rev)
print(high_params_pre)
#$a_r
#[1] -0.0365823

#$b_r
#[1] 0.6432465

#$c_r
#[1] 0.3274432

#> print(high_params_pre)
#$a_p
#[1] 0.2970286

#$b_p
#[1] 0.2970286

#$c_p
#[1] 0.1152017



dat_FT_low<-data.frame(ft_data$low,high_params_rev,high_params_pre)
lpcold<-lpcold+geom_point(aes(x=x,y=all_mean,colour="Fixed Backbone & Flexible Sidechain"),dat_FT_low)
lpcold<-lpcold+geom_smooth(aes(x=x,y=y,colour="Fixed Backbone & Flexible Sidechain"),subset(dat_FT_low,x>=0),method=gam,formula=y~I(0.3274432*exp(-0.0365823 *x)+0.6432465)) 
lpcold<-lpcold+geom_smooth(aes(x=x,y=y,colour="Fixed Backbone & Flexible Sidechain"),subset(dat_FT_low,x<0),method=gam,formula=y~I(0.1152017*exp(0.2970286 *x)+0.2970286)) 




### F F data

ff_data<-getdata('./bb_F_ch_F/')

high_params_rev<-getfit(dat=ff_data$high, type="rev")
high_params_pre<-getfit(dat=ff_data$high, type="pre")
print(high_params_rev)
print(high_params_pre)

#$a_r
#[1] -0.1800182

#$b_r
#[1] 0.495713

#$c_r
#[1] 0.1531873

#> print(high_params_pre)
#$a_p
#[1] 0.09904182

#$b_p
#[1] 0.8261907

#$c_p
#[1] 0.134038



dat_FF_high<-data.frame(ff_data$high,high_params_rev,high_params_pre)


lphot<-lphot+geom_point(aes(x=x,y=all_mean,colour="Fixed Backbone & Sidechain"),dat_FF_high)
lphot<-lphot+geom_smooth(aes(x=x,y=y,colour="Fixed Backbone & Sidechain"),subset(dat_FF_high,x>=0),method=gam,formula=y~I(0.1531873*exp(-0.1800182*x)+0.495713))
lphot<-lphot+geom_smooth(aes(x=x,y=y,colour="Fixed Backbone & Sidechain"),subset(dat_FF_high,x<0),method=gam,formula=y~I(0.134038*exp(0.09904182*x)+ 0.8261907))

high_params_rev<-getfit(dat=ff_data$low, type="rev")
high_params_pre<-getfit(dat=ff_data$low, type="pre")
print(high_params_rev)
print(high_params_pre)

#> print(high_params_rev)
#$a_r
#[1] -0.09449502

#$b_r
#[1] 0.791287

#$c_r
#[1] 0.1893194

#> print(high_params_pre)
#$a_p
#[1] 0.3741811

#$b_p
#[1] 0.5395779

#$c_p
#[1] 0.1025782


dat_FF_low<-data.frame(ff_data$low,high_params_rev,high_params_pre)
lpcold<-lpcold+geom_point(aes(x=x,y=all_mean,colour="Fixed Backbone & Sidechain"),dat_FF_low)
lpcold<-lpcold+geom_smooth(aes(x=x,y=y,colour="Fixed Backbone & Sidechain"),subset(dat_FF_low,x>=0),method=gam,formula=y~I(0.1893194*exp( -0.09449502*x)+0.791287)) 
lpcold<-lpcold+geom_smooth(aes(x=x,y=y,colour="Fixed Backbone & Sidechain"),subset(dat_FF_low,x<0),method=gam,formula=y~I(0.1025782*exp(0.3741811*x)+0.5395779)) 


### T F data

tf_data<-getdata('./bb_T_ch_F/')

high_params_rev<-getfit(dat=tf_data$high, type="rev")
high_params_pre<-getfit(dat=tf_data$high, type="pre")
print(high_params_rev)
print(high_params_pre)
#$a_r
#[1] -0.1813503

#$b_r
#[1] 0.456604

#$c_r
#[1] 0.1646182

#> print(high_params_pre)
#$a_p
#[1] 0.1342193

#$b_p
#[1] 0.8303031

#$c_p
#[1] 0.1277431



dat_TF_high<-data.frame(tf_data$high,high_params_rev,high_params_pre)


lphot<-lphot+geom_point(aes(x=x,y=all_mean,colour="Flexible Backbone & Fixed Sidechain"),dat_TF_high)
lphot<-lphot+geom_smooth(aes(x=x,y=y,colour="Flexible Backbone & Fixed Sidechain"),subset(dat_TF_high,x>=0),method=gam,formula=y~I( 0.1646182*exp(  -0.1813503*x)+ 0.456604))
lphot<-lphot+geom_smooth(aes(x=x,y=y,colour="Flexible Backbone & Fixed Sidechain"),subset(dat_TF_high,x<0),method=gam,formula=y~I( 0.1277431*exp(0.1342193*x)+0.8303031))


high_params_rev<-getfit(dat=tf_data$low, type="rev")
high_params_pre<-getfit(dat=tf_data$low, type="pre")
print(high_params_rev)
print(high_params_pre)
#$a_r
#[1] -0.08342557

#$b_r
#[1] 0.7494862

#$c_r
#[1] 0.2206749

#> print(high_params_pre)
#$a_p
#[1] 0.3811864

#$b_p
#[1] 0.5284677

#$c_p
#[1] 0.09707901



dat_TF_low<-data.frame(tf_data$low,high_params_rev,high_params_pre)
lpcold<-lpcold+geom_point(aes(x=x,y=all_mean,colour="Flexible Backbone & Fixed Sidechain"),dat_TF_low)
lpcold<-lpcold+geom_smooth(aes(x=x,y=y,colour="Flexible Backbone & Fixed Sidechain"),subset(dat_TF_low,x>=0),method=gam,formula=y~I(0.2206749*exp(-0.08342557 *x)+0.7494862)) 
lpcold<-lpcold+geom_smooth(aes(x=x,y=y,colour="Flexible Backbone & Fixed Sidechain"),subset(dat_TF_low,x<0),method=gam,formula=y~I(0.09707901*exp(0.38118644*x)+ 0.5284677)) 


temp<-plot_grid(lpcold,lphot,nrow=2,labels=c("A","B"))
save_plot("ben_both.pdf",temp,base_aspect_ratio=1,ncol=2,nrow=2)


