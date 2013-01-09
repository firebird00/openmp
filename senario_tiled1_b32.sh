## Give the Job a descriptive name
#PBS -N par-lab-ask2_ni 

## Output and error files
#PBS -o run_omp_tiled_b32.out
#PBS -e run_omp_tiled_b32.err

## Limit memory, runtime etc.
##PBS -l walltime=01:00:00
##PBS -l pmem=1gb

## How many nodes:processors_per_node should we get?
#PBS -l nodes=1:clover:highmem:ppn=8

## Start
## Run the job (use full paths to make sure we execute the correct things
## Just replace the path with your local path to openmp file
openmp_exe=/home/parallel/parlab11/ask2_nik/tiled

for N in 512 1024 2048 4096
do
	#Execute OpenMP executable
	for t in 1 2 3 4 5 6 7 8
	do
		#export OMP_NUM_THREADS=$t
		$openmp_exe $N 32   $t   0
	done
done 
