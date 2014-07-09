library(data.table)
library("data.table")
install.packages("data.table")
install.packages("ggplot2")
library(ggplot2)
library(plyr)
library(data.table)
f <- read.table("matches.processed.tab")
f
plot(f$V3,
f <- read.table("matches.processed.tab", header=T)
plot(f$id, f$interface_id)
?setkey
setkey(f, id)
setkey(f, 'id')
setkey(f, c('id'))
?setkey
setkey(f, id)
id
?id
f <- read.table("matches.processed.tab", header=T)
plot(f$protein_id, f$interface_id)
plot(f$protein_id, f$interface_id)
setkey(f, protein_id)
f <- fread("matches.processed.tab", header=T)
setkey(f, protein_id)
f[58[
f[58]
f[59]
f[60]
f
f[60,]
?setkey
f["60"]
f[1]
f[2]
?setkey
f[5]
setkey(f, query)
f[protein_id == 59,]
f[protein_id == 60,]
f[protein_id == 58,]
b <- fread("deltablast_report.tab")
summary(b$identical/b$align_len)
summary(round(b$identical/b$align_len, 2))
g <- round(b$identical/b$align_len, 2)
g[g == 58,]
g == 58
g == 59
?any
?count
?where
table(g == 58)
table(g == 59)
table(g == 0.59)
table(g == 0.58)
table(g == 0.56)
table(g == 0.55)
table(g == 0.54)
table(g == 0.58 and )
table(round(b$identical/b$align_len, 2))
table(round(b$identical/b$align_len, 2)*100)
f <- fread("matches.processed.tab", header=T)
table(f$id)
plot(f$id, f$interface_id)
?heatmap
?table
?heatmap
cor(f$id, f$interface_id)
?table
table(f$id, f$interface_id)
?table
as.matrix(table(f$id, f$interface_id))
image(as.matrix(table(f$id, f$interface_id)))
?image
m <- as.matrix(table(f$id, f$interface_id))
m[m == 0]
m
m[0
m[0]
m[1]
m[2]
m[100]
m[400]
m[90]
m[40]
m[40,40]
m[40][40]
m[40]
n
?matrix
m[40,40]
m[40,44]
m[40,39]
m[39,39]
m[4`,39]
m[44,39]
s <- table(f$id, f$interface_id)
s[1,1]
s[1,0]
s[1,4]
s[1,10]
s[s==0]
s[3025]
as.data.frame(s)
d <- as.data.frame(s)
d[d==0]
summary(d)
d$25
d[Freq == 0]
d[Freq == 0,]
d[d$Freq == 0]
d[d$Freq == 0,]
t
typeof(s)
s
id(s)
?table
?table
s <- table(f$id, f$interface_id, useNA="always")
s
image(as.matrix(log(table(f$id, f$interface_id))))
f <- fread("matches.processed.tab", header=T)
image(as.matrix(log(table(f$id, f$interface_id))))
image(as.matrix(log(table(f$id, f$interface_id))))
f <- fread("matches.processed.tab", header=T)
image(as.matrix(log(table(f$id, f$interface_id))))
f <- fread("matches.processed.tab", header=T)
image(as.matrix(log(table(f$id, f$interface_id))))
image(as.matrix(log(table(f$interface_size, f$interface_id))))
log(table(f$interface_size, f$interface_id))
history(1000)
savehistory(file="history.R")
