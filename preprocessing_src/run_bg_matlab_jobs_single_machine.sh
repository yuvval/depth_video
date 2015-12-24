#!/bin/tcsh
echo ==========================================
set MatlabScript=$argv[1]

source ../theano-env/bin/activate
cat $MatlabScript
matlab -nosplash -nodisplay -logfile log_$MatlabScript.txt -r "run $MatlabScript" < /dev/null > &  /dev/null 
