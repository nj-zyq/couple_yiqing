#
#Batch file to run couple
#use"%s/oldname/newname/g" to change input data name
#
cd /home/oalib/evans/couple97
cp test1.dat COUPLE.DAT
couple
cp COUPLE.LOG test1.log
cp COUPLE.PRT test1.prt
rm COUPLE.LOG
rm COUPLE.PRT
cp COUPLE.TL test1.tl
cp COUPLES.TL test1s.tl
cp COUPLEO.TL test1o.tl
cp COUPLEI.TL test1i.tl
cp COUPLE.CPR test1.cpr
cp COUPLES.CPR test1s.cpr
cp COUPLEO.CPR test1o.cpr
cp COUPLEI.CPR test1i.cpr
rm COUPLE.TL
rm COUPLES.TL
rm COUPLEO.TL
rm COUPLEI.TL
rm COUPLE.CPR
rm COUPLES.CPR
rm COUPLEO.CPR
rm COUPLEI.CPR
rm COUPLE.DAT
#
#rm SCRACH8.DAT
#rm SCRACH23.DAT
#rm SCRACH24.DAT
#rm SCRACH31.DAT
#rm SCRACH33.DAT
#rm SCRACH34.DAT
#rm SCRACH41.DAT
#rm SCRACH42.DAT
#rm SCRACH43.DAT
#rm SCRACH44.DAT
