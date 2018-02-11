
library(fda)

#use Growth data that are saved in library "fda"

dim(growth$hgtf)

data_fem  = growth$hgtf
data_man  = growth$hgtm
data_both =cbind(data_fem,data_man)

#plot all data
matplot(growth $age, data_both, col=c(rep("red",54), rep("blue",39)), type="l", 
			lty=c(rep(1,54), rep(3,39)), xlab="age", ylab= "height")

###  GIRLS ###
################
data_mtx = data_fem
tau	     = list(0.05,0.5,0.95)
alpha	 = list(-0.45,0,0.45)
#PEC
set.seed(1234)
output.pec0     = pec.k(data_mtx,nk=2,alpha=-0.43)
output.pec      = lapply(alpha,pec.k,Y= data_mtx, nk=2,reset.tol=50,lab.ini=-output.pec0[[4]])
basis.pec_girls = lapply(output.pec,listBasis,basis0=output.pec[[1]][[2]],num=2)


###  BOYS ###
################
data_mtx = data_man
tau	     = list(0.05,0.5,0.95)
alpha	 = list(-0.45,0,0.45)
#PEC
set.seed(1234)
output.pec0    = pec.k(data_mtx,nk=2,alpha=-0.43)
output.pec     = lapply(alpha,pec.k,Y= data_mtx, nk=2,reset.tol=50,lab.ini=-output.pec0[[4]])
basis.pec_boys = lapply(output.pec,listBasis,basis0=output.pec[[1]][[2]],num=2)

###  BOTH ###
################
data_mtx = data_both
tau	     = list(0.05,0.5,0.95)
alpha	 = list(-0.45,0,0.45)
#PEC
set.seed(1234)
output.pec0     = pec.k(data_mtx,nk=2,alpha=-0.43)
output.pec     = lapply(alpha,pec.k,Y= data_mtx, nk=2,reset.tol=50,lab.ini=-output.pec0[[4]])
basis.pec_both = lapply(output.pec,listBasis,basis0=output.pec[[1]][[2]],num=2)


# plots
########
plot_path = paste0("/Users/PetraB/Desktop/PEC_growth_data",".pdf")
pdf(plot_path, width=12, height=8)
# pdf("/Users/Petra/Desktop/growth_data_PECs.pdf")
par(mfrow=c(2,3)) #par(mfrow=c(1,1))

plot(basis.pec_girls[[2]][,1], type="l", ylim=c(-0.35,0), xlab="age", ylab= "height", main= "1st PC - girls")
lines(basis.pec_girls[[1]][,1], col="blue")
lines(basis.pec_girls[[3]][,1], col="red")

plot(basis.pec_boys[[2]][,1], type="l", ylim=c(-0.35,0), xlab="age", ylab= "height", main= "1st PC - boys")
lines(basis.pec_boys[[1]][,1], col="blue")
lines(basis.pec_boys[[3]][,1], col="red")

plot(basis.pec_both[[2]][,1], type="l", ylim=c(-0.35,0), xlab="age", ylab= "height", main= "1st PC - all")
lines(basis.pec_both[[1]][,1], col="blue")
lines(basis.pec_both[[3]][,1], col="red")

plot(basis.pec_girls[[2]][,2], type="l", ylim=c(-0.5,0.5), xlab="age", ylab= "height", main= "2nd PC - girls")
lines(basis.pec_girls[[1]][,2], col="blue")
lines(basis.pec_girls[[3]][,2], col="red")

plot(basis.pec_boys[[2]][,2], type="l", ylim=c(-1,1), xlab="age", ylab= "height", main= "2nd PC - boys")
lines(-basis.pec_boys[[1]][,2], col="blue")
lines(basis.pec_boys[[3]][,2], col="red")

plot(-basis.pec_both[[2]][,2], type="l", ylim=c(-0.5,0.5), xlab="age", ylab= "height", main= "2nd PC - all")
lines(-basis.pec_both[[1]][,2], col="blue")
lines(-basis.pec_both[[3]][,2], col="red")
#box()
dev.off()
