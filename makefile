
.PHONY : all batch combine check_combine copy_for_download

all :

batch :
	sbatch "code/batch/server.sbatch"

combine :
	sbatch "code/batch/combine_runs.sbatch"

check_combine :
	cat "output/combine.err"

dldir = output/download/
copy_for_download :
	-rm -rf "$(dldir)"
	-mkdir -p "$(dldir)"
	-cp output/*.csv "$(dldir)"
	-cp output/*.xlsx "$(dldir)"
	-cp output/*.out "$(dldir)"
	-cp output/*.err "$(dldir)"

spath := "$$MW:/home/livingstonb/GitHub/Continuous_Time_HA/output/download/*"
download :
	-mkdir -p output/server
	-scp $(spath) output/server/

readme :
	-pandoc readme.md -o readme.pdf
	-xdg-open readme.pdf

clean :
	-rm -rf output/*
	-rm -rf temp/*