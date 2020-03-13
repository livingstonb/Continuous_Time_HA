path0="/home/livingstonb/GitHub/Continuous_Time_HA/output/"

path1="${path0}*.csv"
path2="${path0}*.xlsx"
path3="${path0}*.out"
path4="${path0}*.err"

dlpath="${path0}download/"
cp $path1 $dlpath
cp $path2 $dlpath
cp $path3 $dlpath
cp $path4 $dlpath 
