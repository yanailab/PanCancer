#!/usr/bin/bash
#SBATCH --job-name=Gyn_part1 # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=db3562@nyulangone.org
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=100gb # Job memory request
#SBATCH --time=10:00:00 # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log # Standard output and error log
#SBATCH --partition cpu_short

source ~/.bashrc

if [ "$#" == 3 ]; then
	cd ~/Documents/STARCH
	conda activate STARCH
	python run_STARCH.py -i ~/Documents/Gyn/${1}${2}/${1}${2}${3}/Clones --output Clones --n_clusters 3 --outdir ~/Documents/Gyn/${1}${2}/${1}${2}${3}/Clones
	conda deactivate
	cd ~/Documents/Gyn/${1}${2}/${1}${2}${3}
	conda activate seurat_functions
	module load r/4.0.3
	Rscript ~/Documents/Gyn/Runs/part1.R ${1} ${2} ${3}

elif [ "$#" == 4 ]; then
	cd ~/Documents/STARCH
	conda activate STARCH
	python run_STARCH.py -i ~/Documents/Gyn/${1}${2}${3}/${1}${2}${3}${4}/Clones --output Clones --n_clusters 3 --outdir ~/Documents/Gyn/${1}${2}${3}/${1}${2}${3}${4}/Clones
	conda deactivate
	cd ~/Documents/Gyn/${1}${2}${3}/${1}${2}${3}${4}
	conda activate seurat_functions
	module load r/4.0.3
	Rscript ~/Documents/Gyn/Runs/part1.R ${1} ${2} ${3} ${4}

fi

