## =============================================================================
## Polar BarPlot
## Original Polar Histogram by Christophe Ladroue
## Source: http://chrisladroue.com/2012/02/polar-histogram-pretty-and-useful/
## Modified from original by Christos Hatzis 3.22.2012 (CH)
## Modified from modified by Christie Haskell 7.25.2014 (CHR)
## =============================================================================
polarBarPlot <-
  function(
    df,
    binSize=1,
    spaceBar=0.05,
    spaceItem=0.2,
    spaceFamily=1.2,
    innerRadius=0.3,
    outerRadius=1,
    nguides=3,
    guides=pretty(range(c(0, df$value)), n=nguides, min.n=2),
    alphaStart=-0,
    circleProportion=1,
    direction="inwards",
    familyLabels=TRUE,
    itemSize=3,
    legLabels=NULL,
    legTitle="Source"){

    require(ggplot2)
    require(plyr)

    # ordering
    df<-arrange(df,family,item,score)
	
    # family and item indices
    df$indexFamily <- as.integer(factor(df$family))
    df$indexItem <- with(df, as.integer(factor(item, levels=item[!duplicated(item)])))        
    df$indexScore <- as.integer(factor(df$score))

    df<-arrange(df,family,item,score)

    # define the bins

    vMax <- max(df$value)

    guides <- guides[guides < vMax]
    df$value <- df$value/vMax

    # linear projection  
    affine<-switch(direction,
                   'inwards'= function(y) (outerRadius-innerRadius)*y+innerRadius,
                   'outwards'=function(y) (outerRadius-innerRadius)*(1-y)+innerRadius,
                   stop(paste("Unknown direction")))

    df<-within(df, {
      xmin <- (binSize + spaceBar) + 
        (indexItem - 1) * (spaceItem + (binSize + spaceBar)) +
        (indexFamily - 1) * (spaceFamily - spaceItem)
      xmax <- xmin + binSize
      ymax <- affine(1 - value)
    }
    )

    df<-df[with(df, order(family,item,value)), ]
    df<-ddply(df,.(item),mutate,ymin=c(1,ymax[1:(length(ymax)-1)]))

    # build the guides
    guidesDF<-data.frame(
      xmin=rep(df$xmin,length(guides)),
      y=rep(guides/vMax,1,each=nrow(df)))

    guidesDF<-within(guidesDF,{
      xend<-xmin+binSize+spaceBar
      y<-affine(1-y)
    })


    # Building the ggplot object

    totalLength<-tail(df$xmin+binSize+spaceBar+spaceFamily,1)/circleProportion-0

	xmins = by(df, df$family, function(x) min(x$xmin)) - 0.4
	xmaxs = by(df, df$family, function(x) max(x$xmax)) + 0.4
	ymins = rep(c(max(df$ymin)), dim(xmins)[1]) + 0.01
	ymaxs = rep(c(min(df$ymax)), dim(xmins)[1]) - 0.01
	items = rep('', dim(xmins)[1])
	
	fullXmin = c(xmins, df$xmin)
	fullXmax = c(xmaxs, df$xmax)
	fullYmin = c(ymins, df$ymin)
	fullYmax = c(ymaxs, df$ymax)
	fullItems = c(items, paste(substr(df$item, 1, 1), '-', substr(df$item, 7, 7)))
	newDf = data.frame(fullXmin, fullXmax, fullYmin, fullYmax, fullItems)
			
    # histograms
    p<-ggplot(newDf)+geom_rect(
      aes(
        xmin=fullXmin,
        xmax=fullXmax,
        ymin=fullYmin,
        ymax=fullYmax,
        fill=fullItems)
    )

    # guides  
    p<-p+geom_segment(
      aes(
        x=xmin,
        xend=xend,
        y=y,
        yend=y),
      colour="white",
      data=guidesDF)

    # label for guides
    guideLabels<-data.frame(
      x=0,
      y=affine(1-guides/vMax),
      label=guides
    )

    p<-p+geom_text(
      aes(x=x,y=y,label=label),
      data=guideLabels,
      angle=-alphaStart*180/pi,
      hjust=1,
      size=4)

    # item labels
    readableAngle<-function(x){
      angle<-x*(-360/totalLength)-alphaStart*180/pi+90
      angle+ifelse(sign(cos(angle*pi/180))+sign(sin(angle*pi/180))==-2,180,0)
    }
    readableJustification<-function(x){
      angle<-x*(-360/totalLength)-alphaStart*180/pi+90
      ifelse(sign(cos(angle*pi/180))+sign(sin(angle*pi/180))==-2,1,0)
    }

    dfItemLabels<-ddply(df,.(item),summarize,xmin=xmin[1])
    dfItemLabels<-within(dfItemLabels,{
      x <- xmin +  (binSize + spaceBar)/2
      angle <- readableAngle(xmin +  (binSize + spaceBar)/2)
      hjust <- readableJustification(xmin +  (binSize + spaceBar)/2)
    })

#    p<-p+geom_text(
#      aes(
#        x=x,
#        label=item,
#        angle=angle,
#        hjust=hjust),
#      y=1.02,
#      size=itemSize,
#      vjust=0.5,
#      data=dfItemLabels)

    # family labels
    if(familyLabels){
      #familyLabelsDF<-ddply(df,.(family),summarise,x=mean(xmin+binSize),angle=mean(xmin+binSize)*(-360/totalLength)-alphaStart*180/pi)
      familyLabelsDF<-aggregate(xmin~family,data=df,FUN=function(s) mean(s+binSize/2))
      familyLabelsDF<-within(familyLabelsDF,{
        x<-xmin
        angle<-xmin*(-360/totalLength)-alphaStart*180/pi
      })
	  
      p<-p+geom_text(
        aes(
          x=x,
          label=family,
          angle=angle),
        data=familyLabelsDF,
        hjust=0.5,
        vjust=0,
        y=1.1)
    }  

    # empty background and remove guide lines, ticks and labels
    p<-p+theme(
      panel.background=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank()
    )

    # x and y limits
    p<-p+xlim(0,tail(df$xmin+binSize+spaceFamily,1)/circleProportion)
    p<-p+ylim(0,outerRadius+0.7)

    # project to polar coordinates
    p<-p+coord_polar(start=alphaStart)

    # nice colour scale
    if(is.null(legLabels)) legLabels <- levels(df$score)
    names(legLabels) <- levels(df$score)
    #p<-p+scale_fill_brewer(name=legTitle, palette='Set1',type='qual', labels=legLabels)
	colorVector = hsv(h=sort(rep(seq(0,1,by=0.3), 4)), s=rep(c(0.25, 0.5, 0.75, 1), 6), v=rep(c(0.8, 0.8, 0.8, 0.8), 6))
	colorVector = c('#EEEEEE', colorVector)	# Add one to facilitate the outline as the first element
	p<-p+scale_fill_manual(name=legTitle, values=colorVector)
    p
  }
