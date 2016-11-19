#-----------------------------------     
#    Makefile per SPLASH 
#-----------------------------------     
F77=gfortran 
############## Linux
# SWI= -c -O3 -Wall -fcheck=bounds -llapack -lblas
# COM= -O3 -Wall -fcheck=bounds -llapack -lblas
SWI= -c -O3 -fcheck=all -ffpe-trap=invalid,zero -fbacktrace -C -g3 -fcheck=bounds -llapack -lblas
COM= -O3 -fcheck=all -ffpe-trap=invalid,zero -fbacktrace -C -g3 -fcheck=bounds -llapack -lblas

moduli=	main.o\
	input.o\
	initial.o\
	filiniz.o\
	spline.o\
	splint.o
       

spl: $(moduli)
	$(F77) $(COM) $(moduli) -o spl -llapack -lblas 

main.o     : main.f
	$(F77) $(SWI) main.f
input.o     : input.f
	$(F77) $(SWI) input.f
initial.o     : initial.f
	$(F77) $(SWI) initial.f
filiniz.o     : filiniz.f
	$(F77) $(SWI) filiniz.f
spline.o     : spline.f
	$(F77) $(SWI) spline.f
splint.o     : splint.f
	$(F77) $(SWI) splint.f
