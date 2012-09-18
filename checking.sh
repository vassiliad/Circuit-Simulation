./zice $1 CORRECT 0
cc=$?
./zice $1 CHECKING 1
rc=$?
./check CORRECT CHECKING $cc $rc
#rm CORRECT
#rm CHECKING