# Feb 28 2024
# Xochitl Diaz
# Convert between gene types

#From the website
#https://www.genenames.org/tools/multi-symbol-checker/
# Paste all the genes you want to look for and downlad a csv

gene_names = read.csv("hgnc-symbol-check.csv")

table(gene_names$Match.type)
#Alias symbol Approved symbol Entry withdrawn Previous symbol       Unmatched 
#678           23380              18             901           13851 

check=which(gene_names$Match.type=="Entry withdrawn")
gene_names[check,]

check=which(gene_names$Match.type=="Previous symbol")
gene_names[check,]

check=which(gene_names$Match.type=="Unmatched")
gene_names[check,]
transform = gene_names$Input[check]

modified_list <-gsub("\\.", "-", gene_names$Input[check])
gene_names$Input[check] = modified_list

write.table(gene_names$Input, "new_gene_list.csv", quote = FALSE, row.names = FALSE)

gene_names = read.csv("hgnc-symbol-check2.csv")

table(gene_names$Match.type)
#Alias symbol Approved symbol Entry withdrawn Previous symbol       Unmatched 
#678           23437              18             910           13019 

check=which(gene_names$Match.type=="Unmatched")

sum(grepl("AL", gene_names$Input[check])) # example AL133516-1
# 3381
sum(grepl("AC", gene_names$Input[check])) # example AC138393-3
# 8474
sum(grepl("FO", gene_names$Input[check])) # example FO704657-1
# 11
sum(grepl("BX", gene_names$Input[check])) # example BX284668-6
# 69
sum(grepl("AP", gene_names$Input[check])) # example AP001816-1
#745
sum(grepl("AF", gene_names$Input[check])) # example AF213884-3
#80
sum(grepl("AJ", gene_names$Input[check])) # example AJ239322-2
#12
sum(grepl("LINC", gene_names$Input[check])) # example LINC01238-1
#2