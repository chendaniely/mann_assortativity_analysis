#V2 for Final
search()
library(igraph)
library(rgexf)

writing.dir <- c("/home/sdal/mann2/mann/analysis/")
wd <- c("/home/sdal/mann2/mann/multidisciplinary-diffusion-model-experiments/results/simulations")
setwd(wd)

#ITERATION LISTS
batch.dir.list <- c("02-lens_batch_2015-11-30_16-09-28_457_750",
                    "02-lens_batch_2015-12-01_19-09-26_495_750",
                    "02-lens_batch_2015-12-02_08-30-37_740_750",
                    "02-lens_batch_2015-12-02_17-25-33_776_750",
                    "02-lens_batch_2015-12-02_22-45-10_646_750",
                    "02-lens_batch_2015-12-04_19-59-57_776_750",
                    "02-lens_batch_2015-12-05_15-26-02_776_750",
                    "02-lens_batch_2015-12-06_02-25-29_756_750",
                    "02-lens_batch_2015-12-06_13-14-04_770_750",
                    "02-lens_batch_2015-12-07_13-38-26",
                    "02-lens_batch_2015-12-07_23-06-53",
                    "02-lens_batch_2015-12-08_21-42-32")

btwn.list <- c(21:29,3,31:39,4,41:45)

wn.list <- c(5,51:59,6,61:69,7,71:79,8)

master.matrix <- matrix(rep(NA,(length(btwn.list)*length(wn.list))),nrow=length(btwn.list),ncol=length(wn.list))
 
for(i in 1:length(btwn.list)){ #BETWEEN

    for(j in 1:length(wn.list)){ #WITHIN

        rm(catch.across.runs)
        catch.across.runs <- rep(NA,length(batch.dir.list))
        
        for(k in 1:length(batch.dir.list)){
            setwd(wd)
            if (file.exists(paste(batch.dir.list[k],"/a250_bm-0.",btwn.list[i],"_bs0.1_wm0.",wn.list[j],"_ws0.2_c0.25_r000/output",sep=""))) { #IF N0 1

                setwd(paste(batch.dir.list[k],"/a250_bm-0.",btwn.list[i],"_bs0.1_wm0.",wn.list[j],"_ws0.2_c0.25_r000/output",sep=""))

                if (file.exists("network_of_agents.pout")) { #IF NO 2
                    d <- read.csv("./network_of_agents.pout",header=F)

                    if(tail(d$V1,n=1)==99 & tail(d$V2,n=1)==249) { #IF NO 3
                         #GET MEAN, SD AND RANGE
                        mean.for.d <- apply(d[,4:13],1,mean)
                        sd.for.d <- apply(d[,4:13],1,sd)
                        min.for.d <- apply(d[,4:13],1,min)
                        max.for.d <- apply(d[,4:13],1,max)
                        mean.for.d.neg <- apply(d[,4:8],1,mean)
                        mean.for.d.pos <- apply(d[,9:13],1,mean)
    
                        d <- cbind(d,mean.for.d,mean.for.d.neg,mean.for.d.pos,sd.for.d,min.for.d,max.for.d)
    
                        d$diff.means.pos.neg <- d$mean.for.d.pos-d$mean.for.d.neg
                        
                        d$node.name <- NA
                        d$node.name <- paste("A",d$V2,sep="")
                        
                        for(v in 99:99){ #v is tick number
                            assign(paste("d.t.",v,sep=""), d[d$V1==v,])
                        
                            #LOAD NETWORK
                            e.list.in <- read.table(gzfile("edge_list.gz"))
                            e.list.use <- e.list.in[,1:2]
                            gg <- graph.data.frame(e.list.use,directed=FALSE)
                            
                            V(gg)$state <- get(paste("d.t.",v,sep=""))$diff.means.pos.neg[match(V(gg)$name,get(paste("d.t.",v,sep=""))$node.name)]
                            
                            #DATA WRITING HERE #WITNIN RUNS
                            write(c(date(),"",c(i,"",j,"",k),"ERR: NONE"),file=paste(writing.dir,"AnalysisTracker.txt",sep=""),ncolumns=8,append=TRUE)
                            catch.across.runs[k] <- assortativity(gg,V(gg)$state,directed=FALSE)
                        } #v
                        
                    } #IF NO 3
                    else {
                        write(c(date(),"",c(i,"",j,"",k),"ERR: NO 3"),file=paste(writing.dir,"AnalysisTracker.txt",sep=""),ncolumns=8,append=TRUE) 
                    } #ELSE NO 3
                } #IF NO 2
                 
                else {
                    write(c(date(),"",c(i,"",j,"",k),"ERR: NO 2"),file=paste(writing.dir,"AnalysisTracker.txt",sep=""),ncolumns=8,append=TRUE)
                } #ELSE NO 2
                 
            } #IF NO 1
             
            else {
                write(c(date(),"",c(i,"",j,"",k),"ERR: NO 1"),file=paste(writing.dir,"AnalysisTracker.txt",sep=""),ncolumns=8,append=TRUE)
            } #ELSE NO 1

        } #k

        master.matrix[i,j] <- mean(catch.across.runs,na.rm=TRUE)
        
    } #j
} #i


setwd(writing.dir)
save.image()

#EOF
