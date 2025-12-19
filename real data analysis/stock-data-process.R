library(huge)
#--------------------------------------------------
data("stockdata")
data = data.frame(stockdata$info[,2],t(stockdata$data))
colnames(data)[1] = "sector"
sector_name = unique(data$sector)

sub_1 = subset(data,sector == "Consumer Staples")
sub_2 = subset(data,sector == "Utilities")
sub_3 = subset(data,sector == "Energy")

W = t(rbind(sub_1,sub_2,sub_3)[,-1])
n = nrow(W);p = ncol(W)
X = matrix(0,n - 1,p)

for (i in 2:n) {
  
  X[i - 1,] = log(W[i,]/W[i - 1,])
}

n1 = nrow(X)
Xt = X[(n1 - 251 + 1):n1,]
XS = X[-((n1 - 251 + 1):n1),]

# remove the outliers
abs_vals <- abs(Xt)
threshold <- quantile(abs_vals, probs = 0.999)
rows_to_remove <- apply(abs_vals, 1, function(row) any(row > threshold))
Xt <- Xt[!rows_to_remove, , drop = FALSE]
dim(Xt)

abs_vals <- abs(XS)
threshold <- quantile(abs_vals, probs = 0.999)
rows_to_remove <- apply(abs_vals, 1, function(row) any(row > threshold))
XS <- XS[!rows_to_remove, , drop = FALSE]
dim(XS)


