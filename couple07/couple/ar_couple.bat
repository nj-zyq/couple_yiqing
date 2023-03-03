rem
rem Use ar command that came with g95 to create libcouple.a
rem containing 57 subroutines
rem
ar -rcvf libcouple.a abortc.o amat.o backsub.o bathy.o bdf.o
ar -rcvf libcouple.a bessj0.o bessy0.o bmat.o botdeps.o bounds.o
ar -rcvf libcouple.a cdfun.o cdj0.o cdsrt.o cdtan.o cexib.o cexibd.o
ar -rcvf libcouple.a cexiw.o cexiwb.o cexiwd.o chdh.o 
ar -rcvf libcouple.a choldc.o cjdj.o
ar -rcvf libcouple.a clineqs.o cpres.o
ar -rcvf libcouple.a crscup.o csds.o csymeig.o cupb.o
ar -rcvf libcouple.a cupbw.o cupw.o cupwb.o cvmult.o cwns.o dreal.o 
ar -rcvf libcouple.a file.o foresub.o egnval.o galrkn.o gama0.o
ar -rcvf libcouple.a gama00.o genabc.o header.o hk01.o hk02.o  
ar -rcvf libcouple.a indep.o indexx.o input.o mgs.o mmult.o 
ar -rcvf libcouple.a propm.o prtmat.o prttl.o rcf.o rwprof.o
ar -rcvf libcouple.a tmmult.o zcexib.o zcexiw.o
