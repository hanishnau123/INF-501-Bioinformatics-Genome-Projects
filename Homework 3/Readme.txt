Copy the following files into a dir in monsoon:
	homework.cpp
	homework.h
	homework_cpp.h
	Makefile

Step 01: Download the files into the dir

	Monsoon: cp /common/contrib/classroom/inf503/hw_dataset.fa .
	Monsoon: cp /common/contrib/classroom/inf503/test_genome.fasta .
	
Step 02: Monsoon: make

Step 03: To execute part A
 
	Monsoon: srun --mem=10GB -t 00:60:00 ./homework A /hw3_dataset.fa ./test_genome.fasta
 
Step 05: To execute part B
 
	Monsoon: srun --mem=10GB -t 00:60:00 ./homework B /hw3_dataset.fa ./test_genome.fasta
	
Step 06: To execute part C
 
	Monsoon: srun --mem=10GB -t 00:60:00 ./homework C /hw3_dataset.fa ./test_genome.fasta

Step 07: To execute part D
 
	Monsoon: srun --mem=10GB -t 00:60:00 ./homework D /hw3_dataset.fa ./test_genome.fasta

