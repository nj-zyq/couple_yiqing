rem
rem Compile and link couple07
rem
del COUPLE.EXE
rem
lf90 -c COUPLE.FOR
rem lf90 -c ABORTC.FOR
rem lf90 -c AMAT.FOR
rem lf90 -c BACKSUB.FOR 
rem lf90 -c BATHY.FOR 
rem lf90 -c BDF.FOR 
rem lf90 -c BESSJ0.FOR 
rem lf90 -c BESSY0.FOR 
rem lf90 -c BMAT.FOR
rem lf90 -c BOTDEPS.FOR 
rem lf90 -c BOUNDS.FOR 
rem lf90 -c CDFUN.FOR
rem lf90 -c CDJ0.FOR 
rem lf90 -c CDSRT.FOR 
rem lf90 -c CDTAN.FOR 
rem lf90 -c CEXIB.FOR
rem lf90 -c CEXIBD.FOR
rem lf90 -c CEXIW.FOR
rem lf90 -c CEXIWB.FOR
rem lf90 -c CEXIWD.FOR
rem lf90 -c CHDH.FOR 
rem lf90 -c CHOLDC.FOR
rem lf90 -c CJDJ.FOR 
rem lf90 -c CLINEQS.FOR 
rem lf90 -c CPRES.FOR 
rem lf90 -c CRSCUP.FOR
rem lf90 -c CSDS.FOR 
rem lf90 -c CSYMEIG.FOR 
rem lf90 -c CUPB.FOR 
rem lf90 -c CUPBW.FOR 
rem lf90 -c CUPW.FOR 
rem lf90 -c CUPWB.FOR 
rem lf90 -c CVMULT.FOR 
rem lf90 -c CWNS.FOR 
rem lf90 -c DREAL.FOR 
rem lf90 -c FILE.FOR 
rem lf90 -c FORESUB.FOR
rem lf90 -c EGNVAL.FOR
rem lf90 -c GALRKN.FOR 
rem lf90 -c GAMA0.FOR 
rem lf90 -c GAMA00.FOR 
rem lf90 -c GENABC.FOR 
lf90 -c HEADER.FOR 
rem lf90 -c HK01.FOR 
rem lf90 -c HK02.FOR 
rem lf90 -c INDEP.FOR 
rem lf90 -c INDEXX.FOR 
rem lf90 -c INPUT.FOR 
rem lf90 -c MGS.FOR 
rem lf90 -c MMULT.FOR 
rem lf90 -c PROPM.FOR 
rem lf90 -c PRTMAT.FOR 
rem lf90 -c PRTTL.FOR
rem lf90 -c RCF.FOR 
rem lf90 -c RWPROF.FOR 
rem lf90 -c TMMULT.FOR 
rem lf90 -c ZCEXIB.FOR
rem lf90 -c ZCEXIW.FOR
lf90 @COUPLE.LNK
