output <- read.table("digitalorganisms/allgenomesaminoacids/miscellaneous_genes/v2/blast_all_raw.txt", sep="\t")
names(output)<- c("Gene", "Hit Name", "E-Value", "Score", "Coverage")
summary(output$`E-Value`)
all(2e-15 >= output$`E-Value`) ## Check your data is filtered by e-value
output <- output %>% filter(`E-Value`<= 1e-10)
output <- output %>% filter(Coverage >= 75)
all(2e-15 >= output$`E-Value`)
summary(output$`E-Value`) ## evalues between 0 and 1.2e-15)

output <- rbind(output[grep(pattern = "K03781", x = output$Gene),],output[grep(pattern = "K03782", x = output$Gene),])

output$Genome <- output$`Hit Name`
output$KO <- output$Gene

#need to remove the "-contigs" in the dataframe
output <- data.frame(lapply(output, function(x) {gsub("-contigs", "", x)}))
#can now parse to get our unique IDs (last part of the hit name)
library("tidyverse")
for (i in 1:nrow(output)){
  string  <- strsplit(output[i,6], "_")
  end <- length(string[[1]])
  output$Genome[i]<-string[[1]][end]
}

#here we are going to work with some of the catalase genes and plot them as histograms in relation to the phylogeny

file <- output
file <- file %>% filter(KO %in% c("K03781", "K03782"))
for (a in 1:nrow(file)){
  name <- file[a,6]
  index <- which(genomesnames$contigs_db_path == name)
  if (length(index)){file[a,6]<-(genomesnames$name[index])}
  else {file[a,6] <- file[a,6]}
}
file<-distinct(file)
#now we want to calculate how many copies per genome
catalase <- as.data.frame(genomesnames$name)
catalase$K03781 <- 0
catalase$K03782 <- 0
for (genome in 1:length(genomesnames$name)){
  catalase[genome,2] <- nrow(filter(filter(file, Genome == genomesnames$name[genome]),KO=="K03781"))
  catalase[genome,3] <- nrow(filter(filter(file, Genome == genomesnames$name[genome]),KO=="K03782"))
}
catalase$total <- catalase$K03781 + catalase$K03782
catalaselong <- gather(catalase[,1:3], key = "Catalase",value="count", K03781:K03782)
names(catalaselong)[1]<-"Genomes"

write.csv(x = file, file = "newblast/Poster/Catalase_Raw.csv", quote = F)
write.csv(x = catalase, file = "newblast/Poster/Catalase_Totals.csv", quote = F)

#we want to then import the tree we want
tree <- read.tree("phylogeny/midpoint_trees/RAxML_midpoint_mags_2_1")
# get the labels on the tree properly
clade_position <- data.frame(tree$tip.label,
                             c(1:339)
)
colnames(clade_position)<-c("label", "order")
genomes <- read.table("digitalorganisms/allgenomesdatabases/view.txt", sep ="\t", header=T
)
new_order<-match(clade_position$label,genomes$genome_id)
label2_names <- genomes[new_order,2]
d <- data.frame(label = tree$tip.label, label2=label2_names)
tree2 <- full_join(tree, d, by = "label")
test <- ggtree(tree2)
#subs <- length(unique(row.names(clean)))
#xaxiswidth <- 2*max(test$data$x) + (0.05*ncol(combined) + 0.05*subs + 0.05*3)*max(test$data$x)
tree3 <- groupClade(tree2, c(527, 497, 414, 402, 394))
#add , color=group to aes for colored tips and tip labels
p<-ggtree(tree3) + theme_tree2() + geom_tiplab(aes(label=label2),size=2,align=TRUE) + ggtree::vexpand(.1, -1) +theme(legend.position = "none",
                                                                                                                                             axis.title.y = element_blank(),
                                                                                                                                             plot.title = element_text(size = 12, 
                                                                                                                                                                       face = "bold",
                                                                                                                                                                       hjust = 0.5,
                                                                                                                                                                       vjust = -5)) + hexpand(.15)  + ylim(-1,340) + 
  geom_hilight(node=527, fill=(brewer.pal(10, "Paired"))[1], alpha=.6, extendto=2.6) + 
  geom_hilight(node=497, fill=(brewer.pal(10, "Paired"))[3], alpha=.6, extendto=2.575) + 
  geom_hilight(node=414, fill=(brewer.pal(10, "Paired"))[5], alpha=.6, extendto=2.5)+ 
  geom_hilight(node=402, fill=(brewer.pal(10, "Paired"))[7], alpha=.6, extendto=2.65) + 
  geom_hilight(node=394, fill=(brewer.pal(10, "Paired"))[9], alpha=.6, extendto=2.355) + 
  geom_cladelabel(node=527, label="A. macleodii", color=(brewer.pal(10, "Paired"))[2], offset=1.22, angle=90, hjust='center', offset.text=.05, barsize=1.5, fontsize=4) + 
  geom_cladelabel(node=497, label="A. mediterranea", color=(brewer.pal(10, "Paired"))[4], offset=1.15, angle=90, hjust='center', offset.text=.05, barsize=1.5, fontsize=4) + 
  geom_cladelabel(node=414, label="A. australica", color=(brewer.pal(10, "Paired"))[6], offset=1.17, angle=90, hjust='center', offset.text=.05, barsize=1.5, fontsize=4)  + 
  geom_cladelabel(node=402, label="A. naphthalenivorans", color=(brewer.pal(10, "Paired"))[8], offset=1.3, angle=90, hjust='center', offset.text=.05, barsize=1.5, fontsize=4)  + 
  geom_cladelabel(node=394, label="A. stellipolaris", color=(brewer.pal(10, "Paired"))[10], offset=1.02, hjust='center', angle=90, offset.text=.05, barsize=1.5, fontsize=4) +
  geom_tiplab(aes(label=label2),size=2,align=TRUE) 

h2 <- gheatmap(p, completeness, offset = 1.09,                               # offset shifts the heatmap to the right,
               width = 0.22,                              # width defines the width of the heatmap column,
               # color defines the boarder of the heatmap columns
               colnames = T, colnames_angle=0, font.size=2, color=F, hjust=.5) + theme(legend.position = "bottom", legend.title = element_text(size = 10), legend.text = element_text(size = 5), legend.box = "horizontal", legend.margin = margin()) + scale_fill_viridis_c( direction=-1, name = "Completeness")
h2 <- h2 + new_scale_fill()
h2 <- gheatmap(h2, numcontigs, offset = 1.52,                               # offset shifts the heatmap to the right,
               width = 0.22,                              # width defines the width of the heatmap column,
               # color defines the boarder of the heatmap columns
               colnames = T, colnames_angle=0, font.size=2, color=F, hjust=.5) + theme(legend.position = "bottom", legend.title = element_text(size = 10), legend.text = element_text(size = 5), legend.box = "horizontal", legend.margin = margin()) + scale_fill_viridis_c(option="rocket", direction=-1, name = "Contigs")
h2 <- h2 + new_scale_fill()
h2 <- gheatmap(h2, mags, offset = 1.94,                               # offset shifts the heatmap to the right,
               width = 0.22,                              # width defines the width of the heatmap column,
               # color defines the boarder of the heatmap columns
               colnames = T, colnames_angle=0, font.size=2, color=F, hjust=.5) + theme(legend.position = "bottom", legend.title = element_text(size = 10), legend.text = element_text(size = 5), legend.box = "horizontal", legend.margin = margin()) + scale_fill_manual(values = c(brewer.pal(7, name = "Greys" )[2], brewer.pal(7, name = "Greys" )[7]), name = "Type of Genome")
h2 <- h2 + new_scale_fill()

# catalaseplot <- filter(catalaselong, Catalase == "K03781")
# catalaseplot <- catalaseplot[,c(1,3)]
# catalaseplot2 <- as.data.frame(catalaseplot[,c(2)])
# rownames(catalaseplot2) <- catalaseplot$Genomes
# names(catalaseplot2)<- "K03781"
# catalaseplot2$K03781 <- as.factor(catalaseplot2$K03781)
# h2 <- gheatmap(h2, catalaseplot2, offset = 7,                               # offset shifts the heatmap to the right,
#                width = 2,                              # width defines the width of the heatmap column,
#                # color defines the boarder of the heatmap columns
#                colnames = T, colnames_angle=0, font.size=2, color=F, hjust=.5) + theme(legend.position = "bottom", legend.title = element_text(size = 10), legend.text = element_text(size = 5), legend.box = "horizontal", legend.margin = margin()) + scale_fill_manual(values = c("white", "dodgerblue"), name = "K03781")
# h2 <- h2 + new_scale_fill()
# 
# catalaseplot <- filter(catalaselong, Catalase == "K03782")
# catalaseplot <- catalaseplot[,c(1,3)]
# catalaseplot2 <- as.data.frame(catalaseplot[,c(2)])
# rownames(catalaseplot2) <- catalaseplot$Genomes
# names(catalaseplot2)<- "K03782"
# catalaseplot2$K03782 <- as.factor(catalaseplot2$K03782)
# h2 <- gheatmap(h2, catalaseplot2, offset =10.5,                               # offset shifts the heatmap to the right,
#                width = 2,                              # width defines the width of the heatmap column,
#                # color defines the boarder of the heatmap columns
#                colnames = T, colnames_angle=0, font.size=2, color=F, hjust=.5) + theme(legend.position = "bottom", legend.title = element_text(size = 10), legend.text = element_text(size = 5), legend.box = "horizontal", legend.margin = margin()) + scale_fill_manual(values = c("white", "yellowgreen"), name = "K03782")
# h2 <- h2 + new_scale_fill()

p4 <- facet_plot(h2, panel = 'Catalase', data = catalaselong,
                 geom = geom_barh,
                 mapping = aes(x = count, fill = as.factor(Catalase)),
                 stat='identity', color = "#3F3939") + theme(legend.position = "bottom", legend.text = element_text(size = 8)) + scale_fill_manual(values = c("yellowgreen", "dodgerblue"))
p4 <- facet_widths(p4, widths = c(2, 1))
ggsave(paste0("newblast/Poster/Catalase_genes.pdf"), p4, width=12, height = 35, units = "in")

#here we just want to make a total catalase figure
catalaselong <- gather(catalase[,c(1,4)], key = "Catalase",value="count", total)
names(catalaselong)[1]<-"Genomes"

h2 <- gheatmap(p, completeness, offset = 1.09,                               # offset shifts the heatmap to the right,
               width = 0.22,                              # width defines the width of the heatmap column,
               # color defines the boarder of the heatmap columns
               colnames = T, colnames_angle=0, font.size=2, color=F, hjust=.5) + theme(legend.position = "bottom", legend.title = element_text(size = 10), legend.text = element_text(size = 5), legend.box = "horizontal", legend.margin = margin()) + scale_fill_viridis_c( direction=-1, name = "Completeness")
h2 <- h2 + new_scale_fill()
h2 <- gheatmap(h2, numcontigs, offset = 1.52,                               # offset shifts the heatmap to the right,
               width = 0.22,                              # width defines the width of the heatmap column,
               # color defines the boarder of the heatmap columns
               colnames = T, colnames_angle=0, font.size=2, color=F, hjust=.5) + theme(legend.position = "bottom", legend.title = element_text(size = 10), legend.text = element_text(size = 5), legend.box = "horizontal", legend.margin = margin()) + scale_fill_viridis_c(option="rocket", direction=-1, name = "Contigs")
h2 <- h2 + new_scale_fill()
h2 <- gheatmap(h2, mags, offset = 1.94,                               # offset shifts the heatmap to the right,
               width = 0.22,                              # width defines the width of the heatmap column,
               # color defines the boarder of the heatmap columns
               colnames = T, colnames_angle=0, font.size=2, color=F, hjust=.5) + theme(legend.position = "bottom", legend.title = element_text(size = 10), legend.text = element_text(size = 5), legend.box = "horizontal", legend.margin = margin()) + scale_fill_manual(values = c(brewer.pal(7, name = "Greys" )[2], brewer.pal(7, name = "Greys" )[7]), name = "Type of Genome")
h2 <- h2 + new_scale_fill()

p4 <- facet_plot(h2, panel = 'Catalase', data = catalaselong,
                 geom = geom_barh,
                 mapping = aes(x = count, fill = as.factor(Catalase)),
                 stat='identity', color = "#3F3939") + theme(legend.position = "bottom", legend.text = element_text(size = 8)) + scale_fill_manual(values = "#777777")
p4 <- facet_widths(p4, widths = c(2, 1))
ggsave(paste0("newblast/Poster/Catalase_genes_2.pdf"), p4, width=12, height = 35, units = "in")

catalaseK03781 <- gather(catalase[,c(1,2)], key = "Catalase",value="count", K03781)
names(catalaseK03781)[1]<-"Genomes"
catalaseK03782 <- gather(catalase[,c(1,3)], key = "Catalase",value="count", K03782)
names(catalaseK03782)[1]<-"Genomes"

h2 <- gheatmap(p, completeness, offset = 1.09,                               # offset shifts the heatmap to the right,
               width = 0.22,                              # width defines the width of the heatmap column,
               # color defines the boarder of the heatmap columns
               colnames = T, colnames_angle=0, font.size=2, color=F, hjust=.5) + theme(legend.position = "bottom", legend.title = element_text(size = 10), legend.text = element_text(size = 5), legend.box = "horizontal", legend.margin = margin()) + scale_fill_viridis_c( direction=-1, name = "Completeness")
h2 <- h2 + new_scale_fill()
h2 <- gheatmap(h2, numcontigs, offset = 1.52,                               # offset shifts the heatmap to the right,
               width = 0.22,                              # width defines the width of the heatmap column,
               # color defines the boarder of the heatmap columns
               colnames = T, colnames_angle=0, font.size=2, color=F, hjust=.5) + theme(legend.position = "bottom", legend.title = element_text(size = 10), legend.text = element_text(size = 5), legend.box = "horizontal", legend.margin = margin()) + scale_fill_viridis_c(option="rocket", direction=-1, name = "Contigs")
h2 <- h2 + new_scale_fill()
h2 <- gheatmap(h2, mags, offset = 1.94,                               # offset shifts the heatmap to the right,
               width = 0.22,                              # width defines the width of the heatmap column,
               # color defines the boarder of the heatmap columns
               colnames = T, colnames_angle=0, font.size=2, color=F, hjust=.5) + theme(legend.position = "bottom", legend.title = element_text(size = 10), legend.text = element_text(size = 5), legend.box = "horizontal", legend.margin = margin()) + scale_fill_manual(values = c(brewer.pal(7, name = "Greys" )[2], brewer.pal(7, name = "Greys" )[7]), name = "Type of Genome")
h2 <- h2 + new_scale_fill()

p4 <- facet_plot(h2, panel = 'Catalase K03781', data = catalaseK03781,
                 geom = geom_barh,
                 mapping = aes(x = count, fill = as.factor(Catalase)),
                 stat='identity', color = "#3F3939") + theme(legend.position = "bottom", legend.text = element_text(size = 8)) + scale_fill_manual(values = "yellowgreen")

p5 <- facet_plot(p4, panel = 'Catalase K03782', data = catalaseK03782,
                 geom = geom_barh,
                 mapping = aes(x = count, fill = as.factor(Catalase)),
                 stat='identity', color = "#3F3939") + theme(legend.position = "bottom", legend.text = element_text(size = 8))  + scale_fill_manual(values = c("yellowgreen", "dodgerblue"))
p5 <- facet_widths(p5, widths = c(4, 1,1))
ggsave(paste0("newblast/Poster/Catalase_genes_3.pdf"), p5, width=12, height = 35, units = "in")

########################
## Here we want to recreate the figures above, but this time only with specific genomes. In this case, we are only using the genomes isolated from depths 40m and above. 
# We are going to include metadata columns, with the completion, contig, and type, then we want a set of environmental data, which will be depth and latitude, and lastly we want a column for catalase like we did above (the two different KOs)

catalase <- read.csv("/pool001/demers/newblast/Poster/Catalase_Totals.csv", row.names = 1)
catalaselong <- gather(catalase[,1:3], key = "Catalase",value="count", K03781:K03782)
names(catalaselong)[1]<-"Genomes"
catalaselong$Catalase <- as.factor(catalaselong$Catalase)

completenesscat <- data.frame("Completeness"=completeness[,3])
row.names(completenesscat) <- row.names(completeness)
new_order<-match(clade_position$label,row.names(completenesscat))
new_names <- row.names(completenesscat)[new_order]
completenesscat <- data.frame("Completeness"=completenesscat[new_order,])
row.names(completenesscat) <- new_names

magscat <- mags
new_order<-match(clade_position$label,row.names(magscat))
new_names <- row.names(magscat)[new_order]
magscat <- data.frame("MAG or Isolate"=magscat[new_order,])
row.names(magscat) <- new_names

numcontigscat <- numcontigs
new_order<-match(clade_position$label,row.names(numcontigscat))
new_names <- row.names(numcontigscat)[new_order]
numcontigscat <- data.frame("Number of Contigs"=numcontigscat[new_order,])
row.names(numcontigscat) <- new_names

depthcat <- depth
new_order<-match(clade_position$label,row.names(depthcat))
new_names <- row.names(depthcat)[new_order]
depthcat <- data.frame("Depth of Isolation"=depthcat[new_order,])
row.names(depthcat) <- new_names

latitude <- read.csv("/pool001/demers/digitalorganisms/alteromonas_metadata.csv")[,c(1,5,33)]
latitude$latitude <- abs(latitude$latitude)
test <- latitude$taxon_oid %in% genomesnames$contigs_db_path
latitude <- latitude[test,]
row.names(latitude)<-latitude$taxon_oid
new1<-match(latitude$taxon_oid,genomesnames$contigs_db_path)
names1 <- genomesnames[new1,1]
row.names(latitude) <- names1
latitudecat <- data.frame("latitude"=latitude[,3])
row.names(latitudecat) <- row.names(latitude)
new_order<-match(clade_position$label,row.names(latitudecat))
new_names <- row.names(latitudecat)[new_order]
latitudecat <- data.frame("ABS(Latitude)"=latitudecat[new_order,])
row.names(latitudecat) <- new_names

sizefraction <- read.csv("/pool001/demers/digitalorganisms/alteromonas_metadata.csv")[,c(1,5,44)]
sizefraction$sizefraction <- abs(sizefraction$sizefraction)
test <- sizefraction$taxon_oid %in% genomesnames$contigs_db_path
sizefraction <- sizefraction[test,]
row.names(sizefraction)<-sizefraction$taxon_oid
new1<-match(sizefraction$taxon_oid,genomesnames$contigs_db_path)
names1 <- genomesnames[new1,1]
row.names(sizefraction) <- names1
sizefractioncat <- data.frame("sizefraction"=sizefraction[,3])
row.names(sizefractioncat) <- row.names(sizefraction)
new_order<-match(clade_position$label,row.names(sizefractioncat))
new_names <- row.names(sizefractioncat)[new_order]
sizefractioncat <- data.frame("ABS(sizefraction)"=sizefractioncat[new_order,])
row.names(sizefractioncat) <- new_names

#we want to then import the tree we want
tree <- read.tree("/pool001/demers/newblast/Poster/Catalase_SurfaceGenomes_RAxML_midpointRootedTree")
# get the labels on the tree properly
clade_position <- data.frame(tree$tip.label,
                             c(1:122)
)
colnames(clade_position)<-c("label", "order")
genomes <- read.table("digitalorganisms/allgenomesdatabases/view.txt", sep ="\t", header=T
)
new_order<-match(clade_position$label,genomes$genome_id)
label2_names <- genomes[new_order,2]
d <- data.frame(label = tree$tip.label, label2=label2_names)
tree2 <- full_join(tree, d, by = 'label')
test <- ggtree(tree2)
#subs <- length(unique(row.names(clean)))
#xaxiswidth <- 2*max(test$data$x) + (0.05*ncol(combined) + 0.05*subs + 0.05*3)*max(test$data$x)
tree3 <- groupClade(tree2, c(184,172,130,161,156))
#add , color=group to aes for colored tips and tip labels
p<-ggtree(tree3) + theme_tree2() + geom_tiplab(aes(label=label2),size=4,align=TRUE) + ggtree::vexpand(.1, -1) +theme(legend.position = "none",
                                                                                                                     axis.title.y = element_blank(),
                                                                                                                     plot.title = element_text(size = 12, 
                                                                                                                                               face = "bold",
                                                                                                                                               hjust = 0.5,
                                                                                                                                               vjust = -5)) + hexpand(.25)  + ylim(-1,123) + 
  geom_hilight(node=184, fill=(brewer.pal(10, "Paired"))[1], alpha=.6, extendto=1.3125) + 
  geom_hilight(node= 172, fill=(brewer.pal(10, "Paired"))[3], alpha=.6, extendto=1.305) + 
  geom_hilight(node=130, fill=(brewer.pal(10, "Paired"))[5], alpha=.6, extendto=1.2175)+ 
  geom_hilight(node=156, fill=(brewer.pal(10, "Paired"))[7], alpha=.6, extendto=1.345) + 
  geom_hilight(node=161, fill=(brewer.pal(10, "Paired"))[9], alpha=.6, extendto=1.16) + 
  geom_cladelabel(node=184, label="A. macleodii", color=(brewer.pal(10, "Paired"))[2], offset=.81, angle=90, hjust='center', offset.text=.02, barsize=1.5, fontsize=4) + 
  geom_cladelabel(node=172, label="A. mediterranea", color=(brewer.pal(10, "Paired"))[4], offset=.825, angle=90, hjust='center', offset.text=.02, barsize=1.5, fontsize=4) + 
  geom_cladelabel(node=130, label="A. australica", color=(brewer.pal(10, "Paired"))[6], offset=.725, angle=90, hjust='center', offset.text=.02, barsize=1.5, fontsize=4)  + 
  geom_cladelabel(node=156, label="A. naphthalenivorans", color=(brewer.pal(10, "Paired"))[8], offset=.8255, angle=90, hjust='center', offset.text=.02, barsize=1.5, fontsize=4)  + 
  geom_cladelabel(node=161, label="A. stellipolaris", color=(brewer.pal(10, "Paired"))[10], offset=.645, hjust='center', angle=90, offset.text=.02, barsize=1.5, fontsize=4) +
  geom_tiplab(aes(label=label2),size=4,align=TRUE) #+ geom_text2(aes(subset=!isTip, label=node), hjust=-.3)

h2 <- gheatmap(p, completenesscat, offset = .79,                               # offset shifts the heatmap to the right,
               width = 0.22,                              # width defines the width of the heatmap column,
               # color defines the boarder of the heatmap columns
               colnames = T, colnames_angle=0, font.size=2, color=F, hjust=.5) + theme(legend.position = "bottom", legend.title = element_text(size = 10), legend.text = element_text(size = 5), legend.box = "horizontal", legend.margin = margin()) + scale_fill_viridis_c( direction=-1, name = "Completeness")
h2 <- h2 + new_scale_fill()
h2 <- gheatmap(h2, numcontigscat, offset = .97,                               # offset shifts the heatmap to the right,
               width = 0.22,                              # width defines the width of the heatmap column,
               # color defines the boarder of the heatmap columns
               colnames = T, colnames_angle=0, font.size=2, color=F, hjust=.5) + theme(legend.position = "bottom", legend.title = element_text(size = 10), legend.text = element_text(size = 5), legend.box = "horizontal", legend.margin = margin()) + scale_fill_viridis_c(option="rocket", direction=-1, name = "Contigs")
h2 <- h2 + new_scale_fill()
h2 <- gheatmap(h2, magscat, offset = 1.14,                               # offset shifts the heatmap to the right,
               width = 0.22,                              # width defines the width of the heatmap column,
               # color defines the boarder of the heatmap columns
               colnames = T, colnames_angle=0, font.size=2, color=F, hjust=.5) + theme(legend.position = "bottom", legend.title = element_text(size = 10), legend.text = element_text(size = 5), legend.box = "horizontal", legend.margin = margin()) + scale_fill_manual(values = c(brewer.pal(7, name = "Greys" )[2], brewer.pal(7, name = "Greys" )[7]), name = "Type of Genome")
h2 <- h2 + new_scale_fill()
h2 <- gheatmap(h2, sizefractioncat, offset = 1.31,                               # offset shifts the heatmap to the right,
               width = 0.22,                              # width defines the width of the heatmap column,
               # color defines the boarder of the heatmap columns
               colnames = T, colnames_angle=0, font.size=2, color=F, hjust=.5) + theme(legend.position = "bottom", legend.title = element_text(size = 10), legend.text = element_text(size = 5), legend.box = "horizontal", legend.margin = margin()) 
h2 <- h2 + new_scale_fill()
h2 <- gheatmap(h2, depthcat, offset = 1.62,                               # offset shifts the heatmap to the right,
               width = 0.22,                              # width defines the width of the heatmap column,
               # color defines the boarder of the heatmap columns
               colnames = T, colnames_angle=0, font.size=2, color=F, hjust=.5) + theme(legend.position = "bottom", legend.title = element_text(size = 10), legend.text = element_text(size = 5), legend.box = "horizontal", legend.margin = margin()) + scale_fill_viridis_c(option="magma", direction=-1, name = "Depth")
h2 <- h2 + new_scale_fill()
h2 <- gheatmap(h2, latitudecat, offset = 1.79,                               # offset shifts the heatmap to the right,
               width = 0.22,                              # width defines the width of the heatmap column,
               # color defines the boarder of the heatmap columns
               colnames = T, colnames_angle=0, font.size=2, color=F, hjust=.5) + theme(legend.position = "bottom", legend.title = element_text(size = 10), legend.text = element_text(size = 5), legend.box = "horizontal", legend.margin = margin()) + scale_fill_viridis_c(option="cividis", direction=-1, name = "Latitude")
h2 <- h2 + new_scale_fill()



# catalaseplot <- filter(catalaselong, Catalase == "K03781")
# catalaseplot <- catalaseplot[,c(1,3)]
# catalaseplot2 <- as.data.frame(catalaseplot[,c(2)])
# rownames(catalaseplot2) <- catalaseplot$Genomes
# names(catalaseplot2)<- "K03781"
# catalaseplot2$K03781 <- as.factor(catalaseplot2$K03781)
# h2 <- gheatmap(h2, catalaseplot2, offset = 7,                               # offset shifts the heatmap to the right,
#                width = 2,                              # width defines the width of the heatmap column,
#                # color defines the boarder of the heatmap columns
#                colnames = T, colnames_angle=0, font.size=2, color=F, hjust=.5) + theme(legend.position = "bottom", legend.title = element_text(size = 10), legend.text = element_text(size = 5), legend.box = "horizontal", legend.margin = margin()) + scale_fill_manual(values = c("white", "dodgerblue"), name = "K03781")
# h2 <- h2 + new_scale_fill()
# 
# catalaseplot <- filter(catalaselong, Catalase == "K03782")
# catalaseplot <- catalaseplot[,c(1,3)]
# catalaseplot2 <- as.data.frame(catalaseplot[,c(2)])
# rownames(catalaseplot2) <- catalaseplot$Genomes
# names(catalaseplot2)<- "K03782"
# catalaseplot2$K03782 <- as.factor(catalaseplot2$K03782)
# h2 <- gheatmap(h2, catalaseplot2, offset =10.5,                               # offset shifts the heatmap to the right,
#                width = 2,                              # width defines the width of the heatmap column,
#                # color defines the boarder of the heatmap columns
#                colnames = T, colnames_angle=0, font.size=2, color=F, hjust=.5) + theme(legend.position = "bottom", legend.title = element_text(size = 10), legend.text = element_text(size = 5), legend.box = "horizontal", legend.margin = margin()) + scale_fill_manual(values = c("white", "yellowgreen"), name = "K03782")
# h2 <- h2 + new_scale_fill()

p4 <- facet_plot(h2, panel = 'Catalase', data = catalaselong,
                 geom = geom_barh,
                 mapping = aes(x = count, fill = Catalase),
                 stat='identity', color = "#3F3939") + theme(legend.position = "bottom", legend.text = element_text(size = 8)) + scale_fill_manual(values = c("yellowgreen", "dodgerblue"))
p4 <- facet_widths(p4, widths = c(2, 1))
ggsave(paste0("/pool001/demers/newblast/Poster/Catalase_genes_4.pdf"), p4, width=22, height = 35, units = "in")


# do all depths, just with those genomes whose completeness - contamination is greater athan 70%
catalase <- read.csv("/pool001/demers/newblast/Poster/Catalase_Totals.csv", row.names = 1)
catalaselong <- gather(catalase[,1:3], key = "Catalase",value="count", K03781:K03782)
names(catalaselong)[1]<-"Genomes"
catalaselong$Catalase <- as.factor(catalaselong$Catalase)

compcon <- read.csv("/pool001/demers/digitalorganisms/alteromonas_metadata.csv")[,c(1,5,26,38,39,44)]
test <- compcon$taxon_oid %in% genomesnames$contigs_db_path
compcon <- compcon[test,]
new1<-match(compcon$taxon_oid,genomesnames$contigs_db_path)
names1 <- genomesnames[new1,1]
compcon$GenomeName_raw <- names1
compcon$Contamination <- as.numeric(compcon$Contamination)
compcon <- compcon %>%
  filter((Completeness-Contamination)>=70)
compcon <- merge(compcon, catalase, by.x="GenomeName_raw", by.y = "genomesnames.name")
genomes <- read.table("/pool001/demers/digitalorganisms/allgenomesdatabases/view.txt", sep ="\t", header=T
)
compcon$Genome.Name...Sample.Name<-genomes$genome[match(compcon$GenomeName_raw,genomes$genome_id)]
write.csv(compcon, "/pool001/demers/digitalorganisms/catalase/allgenomes_comconabove70.csv")

#we want to then import the tree we want
tree <- read.tree("/pool001/demers/digitalorganisms/catalase/catalase_tree_2")
# get the labels on the tree properly
clade_position <- data.frame(tree$tip.label,
                             c(1:214)
)
colnames(clade_position)<-c("label", "order")
genomes <- read.table("/pool001/demers/digitalorganisms/allgenomesdatabases/view.txt", sep ="\t", header=T
)
new_order<-match(clade_position$label,genomes$genome_id)
label2_names <- genomes[new_order,2]
d <- data.frame(label = tree$tip.label, label2=label2_names)
tree2 <- full_join(tree, d, by = 'label')
test <- ggtree(tree2)
#subs <- length(unique(row.names(clean)))
#xaxiswidth <- 2*max(test$data$x) + (0.05*ncol(combined) + 0.05*subs + 0.05*3)*max(test$data$x)
tree3 <- groupClade(tree2, c(334,311,250,287,281))
#add , color=group to aes for colored tips and tip labels
p<-ggtree(tree3) + theme_tree2() + geom_tiplab(aes(label=label2),size=4,align=TRUE) + ggtree::vexpand(.1, -1) +theme(legend.position = "none",
                                                                                                                     axis.title.y = element_blank(),
                                                                                                                     plot.title = element_text(size = 12, 
                                                                                                                                               face = "bold",
                                                                                                                                               hjust = 0.5,
                                                                                                                                               vjust = -5)) + hexpand(.25)  + ylim(-1,215) + 
  geom_hilight(node=334, fill=(brewer.pal(10, "Paired"))[1], alpha=.6, extendto=2.15625) + 
  geom_hilight(node= 311, fill=(brewer.pal(10, "Paired"))[3], alpha=.6, extendto=2.1565) + 
  geom_hilight(node=250, fill=(brewer.pal(10, "Paired"))[5], alpha=.6, extendto=2.075)+ 
  geom_hilight(node=287, fill=(brewer.pal(10, "Paired"))[7], alpha=.6, extendto=2.175) + 
  geom_hilight(node=281, fill=(brewer.pal(10, "Paired"))[9], alpha=.6, extendto=1.975) + 
  geom_cladelabel(node=334, label="A. macleodii", color=(brewer.pal(10, "Paired"))[2], offset=.81, angle=90, hjust='center', offset.text=.02, barsize=1.5, fontsize=4) + 
  geom_cladelabel(node=311, label="A. mediterranea", color=(brewer.pal(10, "Paired"))[4], offset=.765, angle=90, hjust='center', offset.text=.02, barsize=1.5, fontsize=4) + 
  geom_cladelabel(node=250, label="A. australica", color=(brewer.pal(10, "Paired"))[6], offset=.725, angle=90, hjust='center', offset.text=.02, barsize=1.5, fontsize=4)  + 
  geom_cladelabel(node=287, label="A. naphthalenivorans", color=(brewer.pal(10, "Paired"))[8], offset=.8355, angle=90, hjust='center', offset.text=.02, barsize=1.5, fontsize=4)  + 
  geom_cladelabel(node=281, label="A. stellipolaris", color=(brewer.pal(10, "Paired"))[10], offset=.645, hjust='center', angle=90, offset.text=.02, barsize=1.5, fontsize=4) +
  geom_tiplab(aes(label=label2),size=4,align=TRUE) #+ geom_text2(aes(subset=!isTip, label=node), hjust=-.3)

p4 <- facet_plot(p, panel = 'Catalase', data = catalaselong,
                 geom = geom_barh,
                 mapping = aes(x = count, fill = Catalase),
                 stat='identity', color = "#3F3939") + theme(legend.position = "bottom", legend.text = element_text(size = 8)) + scale_fill_manual(values = c("yellowgreen", "dodgerblue"))
p4 <- facet_widths(p4, widths = c(2, 1))
ggsave(paste0("/pool001/demers/newblast/Poster/Catalase_genes_4.pdf"), p4, width=18, height = 35, units = "in")

p<-ggtree(tree3) + theme_tree2() + geom_tiplab(aes(label=label2),size=4,align=TRUE) + ggtree::vexpand(.1, -1) +theme(legend.position = "none",
                                                                                                                     axis.title.y = element_blank(),
                                                                                                                     plot.title = element_text(size = 12, 
                                                                                                                                               face = "bold",
                                                                                                                                               hjust = 0.5,
                                                                                                                                               vjust = -5)) + hexpand(.25)  + ylim(-1,215) + 
  geom_hilight(node=334, fill=(brewer.pal(10, "Paired"))[1], alpha=.6, extendto=2.15625) + 
  geom_hilight(node= 311, fill=(brewer.pal(10, "Paired"))[3], alpha=.6, extendto=2.1565) + 
  geom_hilight(node=250, fill=(brewer.pal(10, "Paired"))[5], alpha=.6, extendto=2.075)+ 
  geom_hilight(node=287, fill=(brewer.pal(10, "Paired"))[7], alpha=.6, extendto=2.175) + 
  geom_hilight(node=281, fill=(brewer.pal(10, "Paired"))[9], alpha=.6, extendto=1.975) + 
  geom_cladelabel(node=334, label="A. macleodii", color=(brewer.pal(10, "Paired"))[2], offset=.81, angle=90, hjust='center', offset.text=.02, barsize=1.5, fontsize=4) + 
  geom_cladelabel(node=311, label="A. mediterranea", color=(brewer.pal(10, "Paired"))[4], offset=.765, angle=90, hjust='center', offset.text=.02, barsize=1.5, fontsize=4) + 
  geom_cladelabel(node=250, label="A. australica", color=(brewer.pal(10, "Paired"))[6], offset=.725, angle=90, hjust='center', offset.text=.02, barsize=1.5, fontsize=4)  + 
  geom_cladelabel(node=287, label="A. naphthalenivorans", color=(brewer.pal(10, "Paired"))[8], offset=.8355, angle=90, hjust='center', offset.text=.02, barsize=1.5, fontsize=4)  + 
  geom_cladelabel(node=281, label="A. stellipolaris", color=(brewer.pal(10, "Paired"))[10], offset=.645, hjust='center', angle=90, offset.text=.02, barsize=1.5, fontsize=4) +
  geom_tiplab(aes(label=label2),size=4,align=TRUE) #+ geom_text2(aes(subset=!isTip, label=node), hjust=-.3)

p4 <- facet_plot(p, panel = 'Catalase', data = catalaselong,
                 geom = geom_barh,
                 mapping = aes(x = count, fill = Catalase),
                 stat='identity', color = "#3F3939") + theme(legend.position = "bottom", legend.text = element_text(size = 8)) + scale_fill_manual(values = c("yellowgreen", "dodgerblue"))
p4 <- facet_widths(p4, widths = c(2, 1))
ggsave(paste0("/pool001/demers/newblast/Poster/Catalase_genes_6.pdf"), p4, width=18, height = 35, units = "in")

p<-ggtree(tree3, layout = "circular") + theme_tree2() + geom_tiplab(aes(label=label2),size=3,align=TRUE) +theme(legend.position = "none",
                                                                                                                axis.title.y = element_blank(),
                                                                                                                plot.title = element_text(size = 12, 
                                                                                                                                          face = "bold",
                                                                                                                                          hjust = 0.5,
                                                                                                                                          vjust = -5)) + 
  geom_hilight(node=334, fill=(brewer.pal(10, "Paired"))[1], alpha=.6, extendto=2.15625) + 
  geom_hilight(node= 311, fill=(brewer.pal(10, "Paired"))[3], alpha=.6, extendto=2.1565) + 
  geom_hilight(node=250, fill=(brewer.pal(10, "Paired"))[5], alpha=.6, extendto=2.075)+ 
  geom_hilight(node=287, fill=(brewer.pal(10, "Paired"))[7], alpha=.6, extendto=2.175) + 
  geom_hilight(node=281, fill=(brewer.pal(10, "Paired"))[9], alpha=.6, extendto=1.975) + 
  geom_cladelabel(node=334, label="A. macleodii", color=(brewer.pal(10, "Paired"))[2], offset=.81, angle=11,hjust='center', offset.text=.02, barsize=1.5, fontsize=6) + 
  geom_cladelabel(node=311, label="A. mediterranea", color=(brewer.pal(10, "Paired"))[4], offset=.765, angle=85, hjust='center', offset.text=.02, barsize=1.5, fontsize=6) + 
  geom_cladelabel(node=250, label="A. australica", color=(brewer.pal(10, "Paired"))[6], offset=.725, angle=21, hjust='center', offset.text=.02, barsize=1.5, fontsize=6)  + 
  geom_cladelabel(node=287, label="A. naphthalenivorans", color=(brewer.pal(10, "Paired"))[8], offset=.8355, angle=-13, hjust='center', offset.text=.02, barsize=1.5, fontsize=6)  + 
  geom_cladelabel(node=281, label="A. stellipolaris", color=(brewer.pal(10, "Paired"))[10], offset=.645, angle=-27, hjust='center',  offset.text=.02, barsize=1.5, fontsize=6) +
  geom_tiplab(aes(label=label2),size=3,align=TRUE) #+ geom_text2(aes(subset=!isTip, label=node), hjust=-.3)
p <- p + new_scale_fill()
h3 <- p + geom_fruit(data = catalaselong, geom = geom_bar, mapping=aes(y=Genomes, x=count, fill = Catalase), offset=.6, pwidth=0.5, orientation="y", stat="identity") + scale_fill_discrete(name = "Catalase Type", type = c("forestgreen", "dodgerblue2"))
ggsave(paste0("/pool001/demers/newblast/Poster/Catalase_genes_7.pdf"), h3, width=31, height = 31, units = "in")

p<-ggtree(tree3, layout = "circular") + theme_tree2() + geom_tippoint() +theme(legend.position = "none",
                                                                               axis.title.y = element_blank(),
                                                                               plot.title = element_text(size = 12, 
                                                                                                         face = "bold",
                                                                                                         hjust = 0.5,
                                                                                                         vjust = -5)) + 
  geom_hilight(node=334, fill=(brewer.pal(10, "Paired"))[1], alpha=.6, extendto=1.4) + 
  geom_hilight(node= 311, fill=(brewer.pal(10, "Paired"))[3], alpha=.6, extendto=1.4) + 
  geom_hilight(node=250, fill=(brewer.pal(10, "Paired"))[5], alpha=.6, extendto=1.4)+ 
  geom_hilight(node=287, fill=(brewer.pal(10, "Paired"))[7], alpha=.6, extendto=1.4) + 
  geom_hilight(node=281, fill=(brewer.pal(10, "Paired"))[9], alpha=.6, extendto=1.4) + 
  geom_cladelabel(node=334, label="A. macleodii", color=(brewer.pal(10, "Paired"))[2], offset=.052, angle=11,hjust='center', offset.text=.02, barsize=1.5, fontsize=5) + 
  geom_cladelabel(node=311, label="A. mediterranea", color=(brewer.pal(10, "Paired"))[4], offset=.005, angle=85, hjust='center', offset.text=.02, barsize=1.5, fontsize=5) + 
  geom_cladelabel(node=250, label="A. australica", color=(brewer.pal(10, "Paired"))[6], offset=.055, angle=21, hjust='center', offset.text=.02, barsize=1.5, fontsize=5)  + 
  geom_cladelabel(node=287, label="A. naphthalenivorans", color=(brewer.pal(10, "Paired"))[8], offset=.063, angle=-12, hjust='center', offset.text=.02, barsize=1.5, fontsize=5)  + 
  geom_cladelabel(node=281, label="A. stellipolaris", color=(brewer.pal(10, "Paired"))[10], offset=.075, angle=-27, hjust='center',  offset.text=.02, barsize=1.5, fontsize=5)
p <- p + new_scale_fill()
h3 <- p + geom_fruit(data = catalaselong, geom = geom_bar, mapping=aes(y=Genomes, x=count, fill = Catalase), offset=.05, pwidth=0.4, orientation="y", stat="identity") + scale_fill_discrete(name = "Catalase Type", type = c("forestgreen", "dodgerblue2"))
ggsave(paste0("/pool001/demers/newblast/Poster/Catalase_genes_8.pdf"), h3, width=31, height = 31, units = "in")
p<-ggtree(tree3, layout = "circular") + theme_tree2() + geom_tippoint() +theme(legend.position = "none",
                                                                               axis.title.y = element_blank(),
                                                                               plot.title = element_text(size = 12, 
                                                                                                         face = "bold",
                                                                                                         hjust = 0.5,
                                                                                                         vjust = -5)) + 
  geom_hilight(node=334, fill=(brewer.pal(10, "Paired"))[1], alpha=.6, extendto=1.4) + 
  geom_hilight(node= 311, fill=(brewer.pal(10, "Paired"))[3], alpha=.6, extendto=1.4) + 
  geom_hilight(node=250, fill=(brewer.pal(10, "Paired"))[5], alpha=.6, extendto=1.4)+ 
  geom_hilight(node=287, fill=(brewer.pal(10, "Paired"))[7], alpha=.6, extendto=1.4) + 
  geom_hilight(node=281, fill=(brewer.pal(10, "Paired"))[9], alpha=.6, extendto=1.4) + 
  geom_cladelabel(node=334, label="A. macleodii", color=(brewer.pal(10, "Paired"))[2], offset=.052, angle=11,hjust='center', offset.text=.02, barsize=1.5, fontsize=5) + 
  geom_cladelabel(node=311, label="A. mediterranea", color=(brewer.pal(10, "Paired"))[4], offset=.005, angle=85, hjust='center', offset.text=.02, barsize=1.5, fontsize=5) + 
  geom_cladelabel(node=250, label="A. australica", color=(brewer.pal(10, "Paired"))[6], offset=.055, angle=21, hjust='center', offset.text=.02, barsize=1.5, fontsize=5)  + 
  geom_cladelabel(node=287, label="A. naphthalenivorans", color=(brewer.pal(10, "Paired"))[8], offset=.063, angle=-12, hjust='center', offset.text=.02, barsize=1.5, fontsize=5)  + 
  geom_cladelabel(node=281, label="A. stellipolaris", color=(brewer.pal(10, "Paired"))[10], offset=.075, angle=-27, hjust='center',  offset.text=.02, barsize=1.5, fontsize=5)
p <- p + new_scale_fill()
df_katE <- catalase[,1:2]
df_KatG <- catalase[,c(1,3)]
names(df_katE)[1]<-"Genomes"
names(df_KatG)[1]<-"Genomes"
h3 <- p + geom_fruit_list(geom_fruit(data = df_katE, geom = geom_bar, mapping=aes(y=Genomes, x=K03781), offset=.05, pwidth=0.15, orientation="y", stat="identity", fill="forestgreen"), 
                          new_scale_fill(),geom_fruit(data = df_KatG, geom = geom_bar, mapping=aes(y=Genomes, x=K03782), offset=.22, pwidth=0.15, orientation="y", stat="identity", fill="dodgerblue2"))
ggsave(paste0("/pool001/demers/newblast/Poster/Catalase_genes_9.pdf"), h3, width=31, height = 31, units = "in")


##### open the ocean microbiomics dataset and filter just to the TARA genomes and then to the genomes that have completeness-contamincation >= 70% and to the small size fraction, MAGs from the 0.2-3µm size fraction, and SRF samples
genomessums <- read.csv("/pool001/demers/DOSTraits/genomes-summary.csv")
genomessums <- genomessums %>% filter((Mean.Completeness-Mean.Contamination)>= 70)
genomessums <- genomessums %>% filter(size.fraction == "0.22-3")
genomessums <- genomessums[grep(x = genomessums$Genome, pattern = "TARA_"),]
genomessums <- genomessums %>% filter(depth <= 10)
write.csv(x=genomessums, file = "/pool001/demers/digitalorganisms/catalase/catalase_genomes.csv")

# here going to open the txt result files and limit them to the genomes in genomessums
katg <- read.table("/pool001/demers/digitalorganisms/catalase/KatG.txt", sep = "\t", header = 1)
katg <- katg %>% filter(sequence.E.value <= 1e-10)
for (i in 1:nrow(katg)){
  katg$scaffold[i] <- paste0(strsplit(x = katg$target.name[i], split = "-")[[1]][1],"-",strsplit(x = katg$target.name[i], split = "-")[[1]][2], "-")
}

katg$genome <- genomesscaffolds$V2[match(katg$scaffold, genomesscaffolds$V1)]
katg <- katg %>% filter(genome %in% genomessums$Genome)
write.table(katg, "/pool001/demers/digitalorganisms/catalase/KatG_filtered.txt", sep = "\t", row.names = F, quote = F)

kate <- read.table("/pool001/demers/digitalorganisms/catalase/KatE.txt", sep = "\t", header = 1)
kate <- kate %>% filter(sequence.E.value <= 1e-10)
for (i in 1:nrow(kate)){
  kate$scaffold[i] <- paste0(strsplit(x = kate$target.name[i], split = "-")[[1]][1],"-",strsplit(x = kate$target.name[i], split = "-")[[1]][2], "-")
}
kate$genome <- genomesscaffolds$V2[match(kate$scaffold, genomesscaffolds$V1)]
kate <- kate %>% filter(genome %in% genomessums$Genome)
write.table(kate, "/pool001/demers/digitalorganisms/catalase/KatE_filtered.txt", sep = "\t", row.names = F, quote = F)


catmn <- read.table("/pool001/demers/digitalorganisms/catalase/CatMn.txt", sep = "\t", header = 1)
catmn <- catmn %>% filter(sequence.E.value <= 1e-10)
for (i in 1:nrow(catmn)){
  catmn$scaffold[i] <- paste0(strsplit(x = catmn$target.name[i], split = "-")[[1]][1],"-",strsplit(x = catmn$target.name[i], split = "-")[[1]][2], "-")
}
catmn$genome <- genomesscaffolds$V2[match(catmn$scaffold, genomesscaffolds$V1)]
catmn <- catmn %>% filter(genome %in% genomessums$Genome)
write.table(catmn, "/pool001/demers/digitalorganisms/catalase/CatMn_filtered.txt", sep = "\t", row.names = F, quote = F)

#############
# Use the hmm results and plot the same thing that we did above, just to see - 
output <- read.table("digitalorganisms/genehitshmms/moreStringent/allhitsforhmms.out")[,c(1,3,5,6)]
names(output)<- c("Hit Name","Gene",  "E-Value", "Score")
output <- rbind(output[grep(pattern = "K03781", x = output$Gene),],output[grep(pattern = "K03782", x = output$Gene),])
summary(output$`E-Value`)
all(1e-10 <= output$`E-Value`) ## Check your data is filtered by e-value
output <- output %>% filter(`E-Value`<= 1e-10)
all(1e-10 <= output$`E-Value`)
summary(output$`E-Value`) ## evalues between 0 and 1.2e-15)


output$Genome <- output$`Hit Name`
output$KO <- output$Gene

#need to remove the "-contigs" in the dataframe
output <- data.frame(lapply(output, function(x) {gsub("-contigs", "", x)}))
#can now parse to get our unique IDs (last part of the hit name)
library("tidyverse")
for (i in 1:nrow(output)){
  string  <- strsplit(output[i,5], "_")
  end <- length(string[[1]])
  output$Genome[i]<-string[[1]][end]
}

output$KO <- gsub(pattern = "hits_moreStringent.align", replacement = "", output$KO)

#here we are going to work with some of the catalase genes and plot them as histograms in relation to the phylogeny

file <- output
for (a in 1:nrow(file)){
  name <- file[a,5]
  index <- which(genomesnames$contigs_db_path == name)
  if (length(index)){file[a,5]<-(genomesnames$name[index])}
  else {file[a,5] <- file[a,5]}
}
file<-distinct(file)
#now we want to calculate how many copies per genome
catalase <- as.data.frame(genomesnames$name)
catalase$K03781 <- 0
catalase$K03782 <- 0
for (genome in 1:length(genomesnames$name)){
  catalase[genome,2] <- nrow(filter(filter(file, Genome == genomesnames$name[genome]),KO=="K03781"))
  catalase[genome,3] <- nrow(filter(filter(file, Genome == genomesnames$name[genome]),KO=="K03782"))
}
catalase$total <- catalase$K03781 + catalase$K03782


write.csv(x = file, file = "newblast/Poster/Catalase_Raw_fromHMMs.csv", quote = F)
write.csv(x = catalase, file = "newblast/Poster/Catalase_Totals_fromHMMs.csv", quote = F)

compcon <- read.csv("/pool001/demers/digitalorganisms/alteromonas_metadata.csv")[,c(1,5,26,38,39,44)]
test <- compcon$taxon_oid %in% genomesnames$contigs_db_path
compcon <- compcon[test,]
new1<-match(compcon$taxon_oid,genomesnames$contigs_db_path)
names1 <- genomesnames[new1,1]
compcon$GenomeName_raw <- names1
compcon$Contamination <- as.numeric(compcon$Contamination)
compcon <- compcon %>%
  filter((Completeness-Contamination)>=70)
compcon <- merge(compcon, catalase, by.x="GenomeName_raw", by.y = "genomesnames$name")
catalase <- catalase %>%
  filter(`genomesnames$name` %in% compcon$GenomeName_raw)
catalaselong <- gather(catalase[,1:3], key = "Catalase",value="count", K03781:K03782)
names(catalaselong)[1]<-"Genomes"
catalaselong$Catalase <- as.factor(catalaselong$Catalase)

genomes <- read.table("/pool001/demers/digitalorganisms/allgenomesdatabases/view.txt", sep ="\t", header=T
) 
compcon$Genome.Name...Sample.Name<-genomes$genome[match(compcon$GenomeName_raw,genomes$genome_id)]
write.csv(compcon, "/pool001/demers/digitalorganisms/catalase/allgenomes_comconabove70_fromHMMs.csv")

#we want to then import the tree we want
tree <- read.tree("/pool001/demers/digitalorganisms/catalase/catalase_tree_2")
# get the labels on the tree properly
clade_position <- data.frame(tree$tip.label,
                             c(1:214)
)
colnames(clade_position)<-c("label", "order")
genomes <- read.table("/pool001/demers/digitalorganisms/allgenomesdatabases/view.txt", sep ="\t", header=T
)
new_order<-match(clade_position$label,genomes$genome_id)
label2_names <- genomes[new_order,2]
d <- data.frame(label = tree$tip.label, label2=label2_names)
tree2 <- full_join(tree, d, by = 'label')
test <- ggtree(tree2)
#subs <- length(unique(row.names(clean)))
#xaxiswidth <- 2*max(test$data$x) + (0.05*ncol(combined) + 0.05*subs + 0.05*3)*max(test$data$x)
tree3 <- groupClade(tree2, c(334,311,250,287,281))
#add , color=group to aes for colored tips and tip labels
p<-ggtree(tree3) + theme_tree2() + geom_tiplab(aes(label=label2),size=4,align=TRUE) + ggtree::vexpand(.1, -1) +theme(legend.position = "none",
                                                                                                                     axis.title.y = element_blank(),
                                                                                                                     plot.title = element_text(size = 12, 
                                                                                                                                               face = "bold",
                                                                                                                                               hjust = 0.5,
                                                                                                                                               vjust = -5)) + hexpand(.25)  + ylim(-1,215) + 
  geom_hilight(node=334, fill=(brewer.pal(10, "Paired"))[1], alpha=.6, extendto=2.15625) + 
  geom_hilight(node= 311, fill=(brewer.pal(10, "Paired"))[3], alpha=.6, extendto=2.1565) + 
  geom_hilight(node=250, fill=(brewer.pal(10, "Paired"))[5], alpha=.6, extendto=2.075)+ 
  geom_hilight(node=287, fill=(brewer.pal(10, "Paired"))[7], alpha=.6, extendto=2.175) + 
  geom_hilight(node=281, fill=(brewer.pal(10, "Paired"))[9], alpha=.6, extendto=1.975) + 
  geom_cladelabel(node=334, label="A. macleodii", color=(brewer.pal(10, "Paired"))[2], offset=.81, angle=90, hjust='center', offset.text=.02, barsize=1.5, fontsize=4) + 
  geom_cladelabel(node=311, label="A. mediterranea", color=(brewer.pal(10, "Paired"))[4], offset=.765, angle=90, hjust='center', offset.text=.02, barsize=1.5, fontsize=4) + 
  geom_cladelabel(node=250, label="A. australica", color=(brewer.pal(10, "Paired"))[6], offset=.725, angle=90, hjust='center', offset.text=.02, barsize=1.5, fontsize=4)  + 
  geom_cladelabel(node=287, label="A. naphthalenivorans", color=(brewer.pal(10, "Paired"))[8], offset=.8355, angle=90, hjust='center', offset.text=.02, barsize=1.5, fontsize=4)  + 
  geom_cladelabel(node=281, label="A. stellipolaris", color=(brewer.pal(10, "Paired"))[10], offset=.645, hjust='center', angle=90, offset.text=.02, barsize=1.5, fontsize=4) +
  geom_tiplab(aes(label=label2),size=4,align=TRUE) #+ geom_text2(aes(subset=!isTip, label=node), hjust=-.3)

p4 <- facet_plot(p, panel = 'Catalase', data = catalaselong,
                 geom = geom_barh,
                 mapping = aes(x = count, fill = Catalase),
                 stat='identity', color = "#3F3939") + theme(legend.position = "bottom", legend.text = element_text(size = 8)) + scale_fill_manual(values = c("yellowgreen", "dodgerblue"))
p4 <- facet_widths(p4, widths = c(2, 1))
ggsave(paste0("/pool001/demers/newblast/Poster/Catalase_genes_6_hmm.pdf"), p4, width=18, height = 35, units = "in")

p<-ggtree(tree3, layout = "circular") + theme_tree2() + geom_tiplab(aes(label=label2),size=3,align=TRUE) +theme(legend.position = "none",
                                                                                                                     axis.title.y = element_blank(),
                                                                                                                     plot.title = element_text(size = 12, 
                                                                                                                                               face = "bold",
                                                                                                                                               hjust = 0.5,
                                                                                                                                               vjust = -5)) + 
  geom_hilight(node=334, fill=(brewer.pal(10, "Paired"))[1], alpha=.6, extendto=2.15625) + 
  geom_hilight(node= 311, fill=(brewer.pal(10, "Paired"))[3], alpha=.6, extendto=2.1565) + 
  geom_hilight(node=250, fill=(brewer.pal(10, "Paired"))[5], alpha=.6, extendto=2.075)+ 
  geom_hilight(node=287, fill=(brewer.pal(10, "Paired"))[7], alpha=.6, extendto=2.175) + 
  geom_hilight(node=281, fill=(brewer.pal(10, "Paired"))[9], alpha=.6, extendto=1.975) + 
  geom_cladelabel(node=334, label="A. macleodii", color=(brewer.pal(10, "Paired"))[2], offset=.81, angle=11,hjust='center', offset.text=.02, barsize=1.5, fontsize=6) + 
  geom_cladelabel(node=311, label="A. mediterranea", color=(brewer.pal(10, "Paired"))[4], offset=.765, angle=85, hjust='center', offset.text=.02, barsize=1.5, fontsize=6) + 
  geom_cladelabel(node=250, label="A. australica", color=(brewer.pal(10, "Paired"))[6], offset=.725, angle=21, hjust='center', offset.text=.02, barsize=1.5, fontsize=6)  + 
  geom_cladelabel(node=287, label="A. naphthalenivorans", color=(brewer.pal(10, "Paired"))[8], offset=.8355, angle=-13, hjust='center', offset.text=.02, barsize=1.5, fontsize=6)  + 
  geom_cladelabel(node=281, label="A. stellipolaris", color=(brewer.pal(10, "Paired"))[10], offset=.645, angle=-27, hjust='center',  offset.text=.02, barsize=1.5, fontsize=6) +
  geom_tiplab(aes(label=label2),size=3,align=TRUE) #+ geom_text2(aes(subset=!isTip, label=node), hjust=-.3)
p <- p + new_scale_fill()
h3 <- p + geom_fruit(data = catalaselong, geom = geom_bar, mapping=aes(y=Genomes, x=count, fill = Catalase), offset=.6, pwidth=0.5, orientation="y", stat="identity") + scale_fill_discrete(name = "Catalase Type", type = c("forestgreen", "dodgerblue2"))
ggsave(paste0("/pool001/demers/newblast/Poster/Catalase_genes_7_hmm.pdf"), h3, width=31, height = 31, units = "in")

p<-ggtree(tree3, layout = "circular") + theme_tree2() + geom_tippoint() +theme(legend.position = "none",
                                                                                                                axis.title.y = element_blank(),
                                                                                                                plot.title = element_text(size = 12, 
                                                                                                                                          face = "bold",
                                                                                                                                          hjust = 0.5,
                                                                                                                                          vjust = -5)) + 
  geom_hilight(node=334, fill=(brewer.pal(10, "Paired"))[1], alpha=.6, extendto=1.4) + 
  geom_hilight(node= 311, fill=(brewer.pal(10, "Paired"))[3], alpha=.6, extendto=1.4) + 
  geom_hilight(node=250, fill=(brewer.pal(10, "Paired"))[5], alpha=.6, extendto=1.4)+ 
  geom_hilight(node=287, fill=(brewer.pal(10, "Paired"))[7], alpha=.6, extendto=1.4) + 
  geom_hilight(node=281, fill=(brewer.pal(10, "Paired"))[9], alpha=.6, extendto=1.4) + 
  geom_cladelabel(node=334, label="A. macleodii", color=(brewer.pal(10, "Paired"))[2], offset=.052, angle=11,hjust='center', offset.text=.02, barsize=1.5, fontsize=5) + 
  geom_cladelabel(node=311, label="A. mediterranea", color=(brewer.pal(10, "Paired"))[4], offset=.005, angle=85, hjust='center', offset.text=.02, barsize=1.5, fontsize=5) + 
  geom_cladelabel(node=250, label="A. australica", color=(brewer.pal(10, "Paired"))[6], offset=.055, angle=21, hjust='center', offset.text=.02, barsize=1.5, fontsize=5)  + 
  geom_cladelabel(node=287, label="A. naphthalenivorans", color=(brewer.pal(10, "Paired"))[8], offset=.063, angle=-12, hjust='center', offset.text=.02, barsize=1.5, fontsize=5)  + 
  geom_cladelabel(node=281, label="A. stellipolaris", color=(brewer.pal(10, "Paired"))[10], offset=.075, angle=-27, hjust='center',  offset.text=.02, barsize=1.5, fontsize=5)
p <- p + new_scale_fill()
h3 <- p + geom_fruit(data = catalaselong, geom = geom_bar, mapping=aes(y=Genomes, x=count, fill = Catalase), offset=.05, pwidth=0.4, orientation="y", stat="identity") + scale_fill_discrete(name = "Catalase Type", type = c("forestgreen", "dodgerblue2"))
ggsave(paste0("/pool001/demers/newblast/Poster/Catalase_genes_8_hmm.pdf"), h3, width=31, height = 31, units = "in")
p<-ggtree(tree3, layout = "circular") + theme_tree2() + geom_tippoint() +theme(legend.position = "none",
                                                                               axis.title.y = element_blank(),
                                                                               plot.title = element_text(size = 12, 
                                                                                                         face = "bold",
                                                                                                         hjust = 0.5,
                                                                                                         vjust = -5)) + 
  geom_hilight(node=334, fill=(brewer.pal(10, "Paired"))[1], alpha=.6, extendto=1.4) + 
  geom_hilight(node= 311, fill=(brewer.pal(10, "Paired"))[3], alpha=.6, extendto=1.4) + 
  geom_hilight(node=250, fill=(brewer.pal(10, "Paired"))[5], alpha=.6, extendto=1.4)+ 
  geom_hilight(node=287, fill=(brewer.pal(10, "Paired"))[7], alpha=.6, extendto=1.4) + 
  geom_hilight(node=281, fill=(brewer.pal(10, "Paired"))[9], alpha=.6, extendto=1.4) + 
  geom_cladelabel(node=334, label="A. macleodii", color=(brewer.pal(10, "Paired"))[2], offset=.052, angle=11,hjust='center', offset.text=.02, barsize=1.5, fontsize=5) + 
  geom_cladelabel(node=311, label="A. mediterranea", color=(brewer.pal(10, "Paired"))[4], offset=.005, angle=85, hjust='center', offset.text=.02, barsize=1.5, fontsize=5) + 
  geom_cladelabel(node=250, label="A. australica", color=(brewer.pal(10, "Paired"))[6], offset=.055, angle=21, hjust='center', offset.text=.02, barsize=1.5, fontsize=5)  + 
  geom_cladelabel(node=287, label="A. naphthalenivorans", color=(brewer.pal(10, "Paired"))[8], offset=.063, angle=-12, hjust='center', offset.text=.02, barsize=1.5, fontsize=5)  + 
  geom_cladelabel(node=281, label="A. stellipolaris", color=(brewer.pal(10, "Paired"))[10], offset=.075, angle=-27, hjust='center',  offset.text=.02, barsize=1.5, fontsize=5)
p <- p + new_scale_fill()
df_katE <- catalase[,1:2]
df_KatG <- catalase[,c(1,3)]
names(df_katE)[1]<-"Genomes"
names(df_KatG)[1]<-"Genomes"
h3 <- p + geom_fruit_list(geom_fruit(data = df_katE, geom = geom_bar, mapping=aes(y=Genomes, x=K03781), offset=.05, pwidth=0.15, orientation="y", stat="identity", fill="forestgreen"), 
                           new_scale_fill(),geom_fruit(data = df_KatG, geom = geom_bar, mapping=aes(y=Genomes, x=K03782), offset=.22, pwidth=0.15, orientation="y", stat="identity", fill="dodgerblue2"))
ggsave(paste0("/pool001/demers/newblast/Poster/Catalase_genes_9_hmm.pdf"), h3, width=31, height = 31, units = "in")
