#
# Unix Makefile for couple. The list of 57 subroutines is
# correct, but it has not been tested. The may be some missing
# tab characters.
#
couple: couple.o \
	abortc.o amat.o backsub.o bathy.o bdf.o \
        bessj0.o bessy0.o bmat.o botdeps.o bounds.o \
	cdfun.o cdj0.o cdsrt.o cdtan.o cexib.o cexibd.o \
        cexiw.o cexiwb.o cexiwd.o chdh.o \
        choldc.o cjdj.o \
        clineqs.o cpres.o \
	crscup.o csds.o csymeig.o cupb.o \
	cupbw.o cupw.o cupwb.o cvmult.o cwns.o dreal.o \
	file.o foresub.o egnval.o galrkn.o gama0.o \
	gama00.o genabc.o header.o hk01.o hk02.o \
	indep.o indexx.o input.o mgs.o mmult.o \
	propm.o prtmat.o prttl.o rcf.o rwprof.o \
	tmmult.o zcexib.o zcexiw.o
#
	gfortran -o couple couple.o \
	abortc.o amat.o backsub.o bathy.o bdf.o \
        bessj0.o bessy0.o bmat.o botdeps.o bounds.o \
	cdfun.o cdj0.o cdsrt.o cdtan.o cexib.o cexibd.o \
        cexiw.o cexiwb.o cexiwd.o chdh.o \
        choldc.o cjdj.o \
        clineqs.o cpres.o \
	crscup.o csds.o csymeig.o cupb.o \
	cupbw.o cupw.o cupwb.o cvmult.o cwns.o dreal.o \
	file.o foresub.o egnval.o galrkn.o gama0.o \
	gama00.o genabc.o header.o hk01.o hk02.o \
	indep.o indexx.o input.o mgs.o mmult.o \
	propm.o prtmat.o prttl.o rcf.o rwprof.o \
	tmmult.o zcexib.o zcexiw.o
#

clean:
	rm *.o

couple.o: couple.for
	gfortran -c couple.for
abortc.o: abortc.for
	gfortran -c abortc.for
amat.o: amat.for 
	gfortran -c amat.for
backsub.o: backsub.for 
	gfortran -c backsub.for
bathy.o: bathy.for 
	gfortran -c bathy.for
bdf.o: bdf.for 
	gfortran -c bdf.for
bessj0.o: bessj0.for 
	gfortran -c bessj0.for
bessy0.o: bessy0.for 
	gfortran -c bessy0.for
bmat.o: bmat.for 
	gfortran -c bmat.for
botdeps.o: botdeps.for
	gfortran -c botdeps.for
bounds.o: bounds.for 
	gfortran -c bounds.for
cdfun.o: cdfun.for
	gfortran -c cdfun.for
cdj0.o: cdj0.for
	gfortran -c cdj0.for
cdsrt.o: cdsrt.for
	gfortran -c cdsrt.for
cdtan.o: cdtan.for
	gfortran -c cdtan.for
cexib.o: cexib.for
	gfortran -c cexib.for
cexibd.o: cexibd.for
	gfortran -c cexibd.for
cexiw.o: cexiw.for
	gfortran -c cexiw.for
cexiwb.o: cexiwb.for
	gfortran -c cexiwb.for
cexiwd.o: cexiwd.for
	gfortran -c cexiwd.for
chdh.o: chdh.for
	gfortran -c chdh.for
choldc.o: choldc.for
	gfortran -c choldc.for
cjdj.o: cjdj.for
	gfortran -c cjdj.for
clineqs.o: clineqs.for
	gfortran -c clineqs.for
cpres.o: cpres.for
	gfortran -c cpres.for
crscup.o: crscup.for
	gfortran -c crscup.for
csds.o: csds.for
	gfortran -c csds.for
csymeig.o: csymeig.for
	gfortran -c csymeig.for
cupb.o: cupb.for 
	gfortran -c cupb.for
cupbw.o: cupbw.for
	gfortran -c cupbw.for
cupw.o: cupw.for
	gfortran -c cupw.for
cupwb.o: cupwb.for
	gfortran -c cupwb.for
cvmult.o: cvmult.for
	gfortran -c cvmult.for
cwns.o: cwns.for
	gfortran -c cwns.for
dreal.o: dreal.for
	gfortran -c dreal.for
file.o: file.for
	gfortran -c file.for
foresub.o: foresub.for
	gfortran -c foresub.for
egnval.o: egnval.for
	gfortran -c egnval.for
galrkn.o: galrkn.for
	gfortran -c galrkn.for
gama0.o: gama0.for
	gfortran -c gama0.for
gama00.o: gama00.for
	gfortran -c gama00.for
genabc.o: genabc.for
	gfortran -c genabc.for
header.o: header.for
	gfortran -c header.for
hk01.o: hk01.for
	gfortran -c hk01.for
hk02.o: hk02.for
	gfortran -c hk02.for
indep.o: indep.for
	gfortran -c indep.for
indexx.o: indexx.for
	gfortran -c indexx.for
input.o: input.for
	gfortran -c input.for
mgs.o: mgs.for
	gfortran -c mgs.for
mmult.o: mmult.for
	gfortran -c mmult.for
propm.o: propm.for
	gfortran -c propm.for
prtmat.o: prtmat.for
	gfortran -c prtmat.for
prttl.o: prttl.for
	gfortran -c prttl.for
rcf.o: rcf.for
	gfortran -c rcf.for
rwprof.o: rwprof.for
	gfortran -c rwprof.for
tmmult.o: tmmult.for
	gfortran -c tmmult.for
zcexib.o: zcexib.for
	gfortran -c zcexib.for
zcexiw.o: zcexiw.for
	gfortran -c zcexiw.for
