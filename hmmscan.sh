#!/bin/bash
#SBATCH -p long # Partition or queue. In this case, short!
#SBATCH --job-name=hmmscan_ncbi_chaperones_sans_NASP2 # Job name
#SBATCH --mail-type=END # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=shla9937@colorado.edu
#SBATCH --nodes=1 # Only use a single node
#SBATCH --ntasks=1 # Run on a single node
#SBATCH --cpus-per-task=60 # cpus
#SBATCH --mem=400gb # Memory limit
#SBATCH --time=100:00:00 # Time limit hrs:min:sec
#SBATCH --output=slurmfiles/slurm_%j.out # Standard output and error log
#SBATCH --error=slurmfiles/slurm_%j.err # %j inserts job number
pwd; hostname; date
module purge
module load hmmer
# commands go here
echo "Hello"
name="ncbi_chaperones_sans_NASP2"
mkdir ${name}
hmmscan -o ${name}/${name}.out --tblout ${name}/${name}_tbl --domtblout ${name}/${name}_domtbl --pfamtblout ${name}/${name}_pfamtbl hmmdb/chaperones_sans_NASP queries/ncbi_archaeal_proteins.fasta
done 
wait
