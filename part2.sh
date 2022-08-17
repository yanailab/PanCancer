#!/usr/bin/bash
#SBATCH --job-name=Tumors_part2 # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=db3562@nyulangone.org
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=200gb # Job memory request
#SBATCH --time=20:00:00 # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log # Standard output and error log
#SBATCH --partition cpu_short

source ~/.bashrc

if [ "$#" == 5 ]; then
	cd ~/Documents/Tumors/${1}${2}/${1}${2}${3}
	conda activate seurat_functions
	module load r/4.0.3
	Rscript ~/Documents/Tumors/Runs/part2.R ${4} ${5}
	cp Subclones/infercnv.png ./infercnv_subclones.png
	cp Subclones/infercnv.median_filtered.png ./infercnv_subclones_filtered.png
	conda deactivate
	
elif [ "$#" == 6 ]; then
	cd ~/Documents/Tumors/${1}${2}${3}/${1}${2}${3}${4}
	conda activate seurat_functions
	module load r/4.0.3
	Rscript ~/Documents/Tumors/Runs/part2.R ${5} ${6}
	cp Subclones/infercnv.png ./infercnv_subclones.png
	cp Subclones/infercnv.median_filtered.png ./infercnv_subclones_filtered.png
	conda deactivate
fi




conda activate seurat_functions
module load r/4.0.3
Rscript ~/Documents/Tumors/Runs/scenic.malignant1.R
conda deactivate

cd Scenic.Malignant
module load python/cpu/3.6.5
cp 1.1_exprMatrix_filtered_t.txt 1.1_exprMatrix_filtered_t.tsv
pyscenic grn -o grnboost_output.csv 1.1_exprMatrix_filtered_t.tsv 1.1_inputTFs.txt
module unload python/cpu/3.6.5
cd ../

conda activate seurat_functions
module load r/4.0.3
Rscript ~/Documents/Tumors/Runs/scenic.malignant2.R
conda deactivate


