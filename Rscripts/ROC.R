# R ROC

# 1 - Black - score
# 2 - Red - score without HPA 
# 3 - Blue - score without Ig
# 4 - Orange - identity without HPA and Ig
# 5 - Cyan - weighted score without HPA and Ig


f1 <- read.table("b1.out")
f2 <- read.table("b2.out")
f3 <- read.table("b3.out")
f4 <- read.table("b4.out")
f5 <- read.table("b5.out")
f6 <- read.table("b6.out")
f7 <- read.table("b7.out")

png("ROC1.png")
plot(f3, ylab="TP rate", xlab="FP rate", col="blue")
# points(f1, col="black") 
# points(f2, col="gray") 
points(f4, col="orange") 
# points(f5, col="cyan") 
points(f6, col="red") 
points(f7, col="green") 
dev.off()

plot(f7, ylab="TP rate", xlab="FP rate", col="blue")
points(f6, col="red") 
