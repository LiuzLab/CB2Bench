library(ScreenBEAM)

args <- commandArgs(trailingOnly = TRUE)
load(args[1])
nx <- ncol(dat)-4
control.samples <- colnames(dat)[5:(5+nx/2-1)]
case.samples <- colnames(dat)[(5+nx/2):(5+nx-1)]
tmp.name <- tempfile()
write.table(dat, file=tmp.name, sep="\t")

df.ret <- ScreenBEAM(
		     input.file = tmp.name,
		     control.samples = control.samples, case.samples = case.samples,
		     control.groupname = 'LOW', case.groupname = 'HIGH',
		     data.type = 'NGS',
		     do.normalization = TRUE,
		     filterLowCount = FALSE,
		     filterBy = 'control',
		     gene.columnId = 2,
		     nitt = 1500,
		     burnin = 500
		     )
#df.gene <- data.frame(gene=df.ret$gene, score=-df.ret[,4])
save(df.ret, file=args[2])
