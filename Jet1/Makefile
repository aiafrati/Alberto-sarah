#-----------------------------------     
#    Makefile per SPLASH 
#-----------------------------------     
#.SUFFIXES: .o

F77=pgf77
#F77=f77

#####   Linux
 SWI= -c #-Wall
# COM= #-Wall
 COM= -llapack -lblas

#####   Digital
#SWI= -c -g -O5 -C -arch host -tune host -inline speed
#SWI= -c -g  -C -arch host -tune host -inline speed
#COM= -O5 -g -C -arch host -tune host -inline speed -non_shared
#COM=  -g -C -arch host -tune host -inline speed -non_shared

########## IBM
#SWI= -c -O3 -qstrict -qipa -qarch=pwr2 -Pk
#SWI= -c -O3 
#SWI= -c  
#COM= -lessl

moduli=	main.o\
	caveta.o\
	cel.o\
	check.o\
	cinfax.o\
	curv.o\
	delta.o\
	der.o\
	der2.o\
	dinam.o\
	distan.o\
	disuni.o\
	doldfil1.o\
	filiniz.o\
	finte.o\
	gauleg.o\
	gauss.o\
	geoini.o\
	geopla.o\
        initia.o\
	initial.o\
	input.o\
	intret.o\
	izero.o\
	izero2.o\
        milan.o\
        nortan.o\
        poli3.o\
        prebem.o\
        press.o\
        ridis.o\
        ridis2.o\
        ropn.o\
        sava.o\
        solpre.o\
	solver.o\
	splone.o\
	splont.o\
	splver.o\
	stampa.o\
	taglio.o\
	tagli2.o\
	wigl.o
#
#spl: $(moduli)
#	$(F77) $(COM) $(moduli) -ldxml -o spl 
spl: $(moduli) 
	$(F77) $(moduli) $(COM)   -o spl 

main.o     : main.f
	$(F77) $(SWI) main.f
caveta.o     : caveta.f
	$(F77) $(SWI) caveta.f
cel.o     : cel.f
	$(F77) $(SWI) cel.f
check.o     : check.f
	$(F77) $(SWI) check.f
cinfax.o     : cinfax.f
	$(F77) $(SWI) cinfax.f
curv.o     : curv.f
	$(F77) $(SWI) curv.f
delta.o     : delta.f
	$(F77) $(SWI) delta.f
der.o     : der.f
	$(F77) $(SWI) der.f
der2.o     : der2.f
	$(F77) $(SWI) der2.f
dinam.o     : dinam.f
	$(F77) $(SWI) dinam.f
distan.o     : distan.f
	$(F77) $(SWI) distan.f
disuni.o     : disuni.f
	$(F77) $(SWI) disuni.f
doldfil1.o     : doldfil1.f
	$(F77) $(SWI) doldfil1.f
filiniz.o     : filiniz.f
	$(F77) $(SWI) filiniz.f
finte.o     : finte.f
	$(F77) $(SWI) finte.f
gauleg.o     : gauleg.f
	$(F77) $(SWI) gauleg.f
gauss.o     : gauss.f
	$(F77) $(SWI) gauss.f
geoini.o     : geoini.f
	$(F77) $(SWI) geoini.f
geopla.o     : geopla.f
	$(F77) $(SWI) geopla.f
initia.o     : initia.f
	$(F77) $(SWI) initia.f
initial.o     : initial.f
	$(F77) $(SWI) initial.f
input.o     : input.f
	$(F77) $(SWI) input.f
intret.o     : intret.f
	$(F77) $(SWI) intret.f
izero.o     : izero.f
	$(F77) $(SWI) izero.f
izero2.o     : izero2.f
	$(F77) $(SWI) izero2.f
milan.o     : milan.f
	$(F77) $(SWI) milan.f
nortan.o     : nortan.f
	$(F77) $(SWI) nortan.f
poli3.o     : poli3.f
	$(F77) $(SWI) poli3.f
prebem.o     : prebem.f
	$(F77) $(SWI) prebem.f
press.o     : press.f
	$(F77) $(SWI) press.f
ridis.o     : ridis.f
	$(F77) $(SWI) ridis.f
ridis2.o     : ridis2.f
	$(F77) $(SWI) ridis2.f
ropn.o     : ropn.f
	$(F77) $(SWI) ropn.f
sava.o     : sava.f
	$(F77) $(SWI) sava.f
solpre.o     : solpre.f
	$(F77) $(SWI) solpre.f
solver.o     : solver.f
	$(F77) $(SWI) solver.f
splone.o     : splone.f
	$(F77) $(SWI) splone.f
splont.o     : splont.f
	$(F77) $(SWI) splont.f
splver.o     : splver.f
	$(F77) $(SWI) splver.f
stampa.o     : stampa.f
	$(F77) $(SWI) stampa.f
taglio.o     : taglio.f
	$(F77) $(SWI) taglio.f
tagli2.o     : tagli2.f
	$(F77) $(SWI) tagli2.f
wigl.o       : wigl.f
	$(F77) $(SWI) wigl.f

