# Run SOM on the monthly change data for a particular SST product

library(kohonen)
library(dplyr)
library(tidyverse)
library(grid)
library(gridExtra)
library(zoo)

# normalize data columns
normalit<-function(m){
  (m - min(m, na.rm=T))/(max(m, na.rm=T)-min(m, na.rm=T))
}


col.par = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#CBD588", "#5F7FC7",
            "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
            "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D",
            "#8A7C64", "#599861")

#load("~/Desktop/SstCP_AVHRR-OI/AVHRR-OI_MonthlyChange.Rdata")
#plot.tit1 <- "AVHRR-OI"

load("~/Desktop/StationCP/Insitu_MonthlyChange.Rdata")
plot.tit1 <- "In Situ"

clust.data <- list()
col.seq <- paste0(as.character(rep(1970:2017, each=2)), c("_median","_sd"))
df = data.frame(matrix(vector(), 12, 98), stringsAsFactors=F)
names(df) <- c("index","month",col.seq)
df$index <- as.factor(df$index)

# combine all the clusters
for (cl in 1:5){
  clust <- stationtrend.list[[cl]]
  clust <- clust[!sapply(clust, is.null)]
  # need to ensure in situ data has equal number of col
  clust <- lapply(clust, function(x) anti_join(df, x, by = "index")%>%
                    bind_rows(x) %>%
                    dplyr::filter(!is.na(index)))
  
  clust.data[[cl]] <- do.call(rbind, clust)
  rm(clust)
}

c1.data <- do.call(rbind, clust.data)

som.stor <- list()

for (m in 1:12){
  data.som <- c1.data %>%
    dplyr::select(index, month, matches("median")) %>%
    dplyr::filter(month==m) %>%
    dplyr::select_if(~ !all(is.na(.)))
  
  
  temp.norm <- t(apply(data.som[,-c(1,2)], 1, function(x) normalit(x)))
  
  # temp.diff <- matrix(NA, nrow=nrow(temp.norm), ncol=ncol(temp.norm))
  # for (j in 1:nrow(temp.norm)){
  #   temp.diff[j,2:ncol(temp.diff)] <- diff(temp.norm[j,], lag=1)
  # }
  # temp.norm <- temp.diff
  
  ind.maxNA <- apply(temp.norm, 1, function(x) sum(is.na(x))/length(x) > .75)
  data.som[,3:ncol(data.som)] <- temp.norm
  rm(temp.norm)
  
  data.som <- data.som[!duplicated(data.som[,3:ncol(data.som)]), ]
  
  names(data.som)[3:ncol(data.som)] <- substr(colnames(data.som[3:ncol(data.som)]),1,4)
  
  plot.tit2 <- paste(plot.tit1, month.name[m], " Cluster ",cl)
  
  X <- as.matrix(data.som[,3:ncol(data.som)])
  res.som <- som(X, grid = somgrid(4, 4, "hexagonal"),
                 rlen=1000, alpha=c(0.05, 0.005), maxNA.fraction=.75, radius=0)
  
  plot(res.som, main = "Cluster 1")
  grp.som <- res.som$unit.classif
  
  plot.som <- cbind(data.frame("grp"=grp.som),data.som[!ind.maxNA,])
  plot.som <- plot.som %>% gather(variable, value, -(grp:month)) %>% arrange(index)
  plot.som$variable <- as.numeric(substr(plot.som$variable,start = 1,stop = 4))
  
  som.stor[[m]] <- plot.som
  
  name.grp <- sort(unique(grp.som))
  n.grp <- length(name.grp)
  plot.list <- list()
  cnt = 0
  for (i in 1:n.grp){
    cnt = cnt+1
    df.cent <- data.frame(x=plot.som$variable,y=as.vector(res.som$codes[[1]][name.grp[i],]))
    
    plot.list[[cnt]] <- ggplot(plot.som %>% dplyr::filter(grp==name.grp[i]), aes(colour=index)) +
      geom_crossbar(aes(x=variable, ymin=value, ymax=value, y=value, group=index),
                    width = 1, fatten = 6) +
      geom_line(aes(x=variable, y=value, group=index), lwd=1.2) +
      geom_crossbar(data=df.cent, aes(x=x, ymin=y, ymax=y, y=y), width=1, col='black') +
      geom_line(data=df.cent, aes(x=x, y=y), color="black", lwd=1.5) +
      ylim(c(0,1)) + ylab("Temperature (°C)") + xlab("Year") +
      scale_color_manual(name="SOM Cluster Stations",
                         values=col.par[1:24],
                         guide=guide_legend(override.aes=list(alpha=.5, size=2))) +
      scale_x_continuous(breaks=min(plot.som$variable):max(plot.som$variable)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size=9),
            panel.grid.major.x=element_line(colour="gray90",size=0.75))
  }
  
  plot(res.som, type="changes")
  plot(res.som, type="dist.neighbours")
  
  ggsave(file = paste0("~/Desktop/StationCP/",month.abb[m],".png"),
         arrangeGrob(grobs = plot.list, ncol = 4),  width = 40, height=16)
  
}

save(som.stor, file="~/Desktop/StationCP/SOM_results.Rdata")

# Plot of SOM series from monthly data ------------------------------------
load("~/Desktop/StationCP/Insitu_DailyChange.Rdata")
load("~/R/Data/site_list_v4.2.RData")

daily.list <- unlist(statdaily.list, recursive = F)

for (k in 1:12){
  plot.som <- som.stor[[k]]
  name.grp <- sort(unique(plot.som$grp))
  for (i in 1:length(name.grp)){
    subset.names <- unique(plot.som$index[which(plot.som$grp==name.grp[i])])
    stat.loc <- list()
    for (n in 1:length(subset.names)){
      stat.name <- subset.names[n]
      stat.ind <- which(site_list$index==stat.name)
      stat.loc[[n]] <- data.frame("station"=stat.name, 
                                  "lon"=site_list$lon[stat.ind],
                                  "lat"=site_list$lat[stat.ind])
      
      stat.ind <- which(names(daily.list)==stat.name)
      t <- daily.list[[stat.ind]][k,]
      t <- t[,-((ncol(t)-1):ncol(t))]             # drop the last two columns median and sd
      
      freq <- max(as.numeric(substr(colnames(t),start=6, stop=7)),na.rm=T)
      st.dom.yr <- c(as.numeric(substr(colnames(t)[1],start = 1, stop = 4)),
                     as.numeric(substr(colnames(t)[1],start = 6, stop = 7)))
      
      yrs.all <- rep(1970:2017, each=freq)
      dom.all <- rep(1:freq, (2017-1970+1))
      
      if (!exists("date.all.lu")){
        date.all.lu <- data.frame("years"=yrs.all, "dom"=dom.all)
        GrA <- date.all.lu
      }
      
      data.vec <- rep(NA, nrow(date.all.lu))
      ind.st <- which(date.all.lu$years==st.dom.yr[1] & date.all.lu$dom==st.dom.yr[2])
      data.vec[ind.st:(ind.st+length(t)-1)] <- as.numeric(t)
      data.vec <- data.frame("replace"=data.vec, stringsAsFactors = F)
      colnames(data.vec) <- stat.name
      GrA <- cbind(GrA, data.vec)
    }
    rm(date.all.lu)
    
    GrA_plot <- GrA %>%
      gather(station, temp, 3:ncol(.)) %>% 
      mutate(date=paste(years,dom,sep="."))
    
    tick.st <- which(GrA_plot$years==min(GrA_plot$years) & GrA_plot$dom==min(GrA_plot$dom))[1]
    tick.en <- which(GrA_plot$years==max(GrA_plot$years) & GrA_plot$dom==max(GrA_plot$dom))[1]

    ggplot(data=GrA_plot, aes(x=date, y=temp, colour=station)) +
      geom_line(aes(group=station)) +
      geom_smooth(aes(x=date, y=temp, group=station), na.rm=T, method="loess") +
      scale_x_discrete('Date', breaks=GrA_plot$date[seq(tick.st,tick.en,31)],
                       labels=substr(GrA_plot$date[seq(tick.st,tick.en,31)],start=1,stop=4)) +
      facet_grid(station ~ ., margins = F, scales = "fixed") +
      theme(panel.background = element_rect(fill = 'gray90', colour = 'gray30'),
            axis.text.x = element_text(angle = 45, hjust = 1, size=9))
  }
}