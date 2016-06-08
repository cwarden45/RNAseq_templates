library("Vennerable")
#install manually

set1.up = as.character(table1$symbol[table1$status == "Trt1 Up"])
set1.down = as.character(table1$symbol[table1$status == "Trt1 Down"])

set2.up = as.character(table2$symbol[table2$status == "Trt2 Up"])
set2.down = as.character(table2$symbol[table2$status == "Trt2 Down"])

gene.list = list(UP1=set1.up,
				DOWN1=set1.down,
				UP2=set2.up,
				DOWN2=set2.down)
vennObj = Venn(gene.list)

png("venn.png")
plot(vennObj, type="ellipses")
dev.off()