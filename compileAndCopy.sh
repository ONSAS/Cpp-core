clear
cd src
make
echo " compilation done."
cp timeStepIteration.lnx $ONSAS_PATH/src
echo " binary copied."
make clean
cd ..
