library(getopt)

spec = matrix(c('input','i',1,"character",'work_dir','w',1,"character",'function_file','f',1,"character",'sample','s',1,"character"),byrow = TRUE,ncol = 4)
opt = getopt(spec)
setwd(opt$work_dir)
source(opt$function_file)

map_df <- as.data.frame(t(na.omit(t(read.table(opt$input,header = FALSE,sep=" ")))))
colnames(map_df) <- c("count","perc","type")
map_df$count <- as.numeric(gsub(" ","",map_df$count))

map_df$perc <- paste(map_df$type,"times (",map_df$count,",",map_df$perc,"%",")",sep = " ")


map <- data.frame(count=map_df$count,type=map_df$type)
map$sample <- rep(opt$sample,nrow(map))

pdf("tmp_mappie.pdf",width=6,height=5)
pieplot(map_df)
dev.off()

write.table(map,"tmp_map.txt",quote = FALSE,sep=" ",row.names = FALSE,col.names = FALSE)
