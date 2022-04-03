The folder contains the following files:
homework.cpp, homework.h, homework_BLAST.h, hw4_dataset.fa, test_genome.fasta, readme
	
Following commands are needed to execute the program:

Step 01: move into the directory

Step 02: make (to make the makefile)

Step 03: To execute part 1 A in monsoon:
 
	 srun --mem=10GB -t 00:60:00 ./homework A /hw4_dataset.fa ./test_genome.fasta
	
Step 04: To check the job status for time:

	 jobstats -j <job_id>

	 if job number is not generated then run

	 jobstats -r
 
Step 05: To execute part 1 B of program:
 
	 srun --mem=10GB -t 00:60:00 ./homework B /hw4_dataset.fa ./test_genome.fasta
	
Step 06: To execute part 2 A of program: 
 
	 srun --mem=10GB -t 00:60:00 ./homework C /hw4_dataset.fa ./test_genome.fasta

Step 07: To execute part 2 B of program:
 
	 srun --mem=10GB -t 00:60:00 ./homework D /hw4_dataset.fa ./test_genome.fasta

