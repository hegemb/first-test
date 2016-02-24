## This is a function made to try to apply combat to our data.
## Try combat
batch = labInfo[,"New_plate"]
edata = log2(exprs(data))
mod = model.matrix(~as.factor(Case_ctrl),data=background)

## Correct for batch using ComBat
combat_edata = ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
## End(Not run)

log2.exprs <- combat_edata

pp <- prcomp(log2.exprs)
pcData <- data.frame(pp$x)
