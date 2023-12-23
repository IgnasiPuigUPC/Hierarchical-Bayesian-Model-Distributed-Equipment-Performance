# ---
# title: "Plot functions"
# author: "Ignasi Puig"
# date: "17/1/2023"
# ---

## Function plot0

# Plots raw data chart (observations) with weekly printing time in Y-axis

# Input variables
# 1. id: machine id
# 2. obs = XY: data.frame with the details of the posterior predictive
# 3. yLim
# 4. title: logical TRUE print title, FALSE do not print it

plot0<- function(id, obs=XY, y1Lim, y2Lim, title = TRUE ) {
  
  obs = obs[obs$id==id,]
  sn = unique(obs$machineSN[obs$id==id])
  
  b <- y1Lim/y2Lim
  
  plot(ggplot(obs, aes(x = week, y = y))+
         geom_line(lty=1,col='black')+
         geom_point(col='black')+
         geom_col(aes(y = off * b), fill = 'blue', alpha = 0.2, width = 0.5) + 
         scale_y_continuous('no. of errors',limits = c(0,y1Lim), breaks = seq(0,100,2),
                            sec.axis = sec_axis(~ ./b, name = "Printing time (h)")) +
         scale_x_continuous(breaks = seq(0,100,2)) +
         labs(title=ifelse(title, paste('Printer ',id,'.',sep=''), ''), x='week'))
  
}


# Function plotFull

# plots the printer values with the probability assigned to each week of being an outlier.

plotFull <- function(id, obs=XY, outlierThr = 0.5, y1Lim, y2Lim, title = TRUE ) {
  
  obs = XY[XY$id==id,]
  # dotSize <- obs$z
  # dotSize[obs$z < outlierThr] <- 0 #all observations below 0.5 are considered in-control
  # # only 3 levels are considered for OOC observations:
  # # low: [0.5, 0.7]
  # # medium: [0.7,0.9]
  # # high: (0.9, 1.0]
  # dotSize <- cut(dotSize, breaks = c(0, seq(0.5, 1, 0.5/3)), labels = c(0,0.25,1.5,3), include.lowest = TRUE)
  
  b <- y1Lim/y2Lim
  
  plot(ggplot(obs)+
         geom_line(aes(x=week,y=y),lty=1,col='black')+
         geom_point(aes(x=week,y=y),col='black')+
         geom_line(aes(x=week,y=meanThetaPost),lty=2,col="black")+
         geom_point(aes(x=week,y=meanThetaPost),pch=1,col="black")+
         geom_point(data=obs[obs$z>=outlierThr,],
                    # aes(x = week, y = y, stroke = as.numeric(as.character(dotSize[obs$z>=outlierThr]))),
                    aes(x = week, y = y, stroke = 2),
                    pch=21,  colour="black",  size = 5) +
         geom_ribbon(aes(x=week,ymin=QminYPred,ymax=QmaxYPred),fill='grey56',alpha=0.2)+
         geom_col(aes(x = week, y = z * b), fill = 'red', width = 0.5, alpha = 0.3) +
         geom_hline(yintercept = c(0.5,1) * b, size = 0.1, linetype = 5, col = 'red', alpha = 0.4) +
         scale_y_continuous(limits = c(0,y1Lim), breaks = seq(0,100,2), 
                            sec.axis = sec_axis(~ . /b, name = 'P(Zjt=1)', breaks = seq(0,1,0.25)))+
         scale_x_continuous(breaks = seq(0,100,2)) +
         theme(plot.margin=unit(c(1,3,1,1), "cm")) +
         theme(axis.title.y.right = element_text(vjust=3))+
         labs(title= ifelse(title, paste('Printer ',id,'.',sep=''), ''), x='week',y='no. of errors'))
}



