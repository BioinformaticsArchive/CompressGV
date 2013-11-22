library(ROCR)

pos <- scan("data/analysis/P")
tp = length(which(pos>=0))
fn = length(which(pos<0))
neg <- scan("data/analysis/N")
fp = length(which(neg>=0))
tn = length(which(neg<0))

sen = tp / (tp + fn)
spec = tn / (tn + fp)

all <- c(pos, neg)
lbl <- c(rep(1,length(pos)), rep(0,length(neg)))

pred <- prediction(all, lbl)

performance(pred, measure="auc")@y.values[[1]]
sen
spec
