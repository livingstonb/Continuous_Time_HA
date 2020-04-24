
.PHONY : all batch combine check_combine copy_for_download

all :

batch :
	sbatch "code/batch/server.sbatch"

combine :
	sbatch "code/batch/combine_runs.sbatch"

check_combine :
	cat "output/combine.err"

dldir = output/download
copy_for_download :
	-cp "output/*.csv" $(dldir)
	-cp "output/*.xlsx" $(dldir)
	-cp "output/*.out" $(dldir)
	-cp "output/*.err" $(dldir)

readme :
	-pandoc readme.md -o readme.pdf
	-xdg-open readme.pdf

clean :
	-rm -rf output/*
	-rm -rf temp/*