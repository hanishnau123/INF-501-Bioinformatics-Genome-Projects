The folder contains the following files:
homework.cpp, homework.h, hw6_dataset.fa
	
Following commands are needed to execute the program:

Step 01: move into the directory

Step 02: make (to make the makefile)
	 monsoon: make

Step 03: To execute part A in monsoon for searching sequences in part A:
 
	 monsoon: srun --mem=10GB -t 00:60:00 ./homework A ./hw6_dataset.fa
	 this can be executed without srun too

	 monsoon: ./homework A ./hw6_dataset.fa	

Step 04: To check the job status for time:

	 jobstats -j <job_id>
 
Step 05: To execute part A II in monsoon for suffix trie size:
 
	 srun --mem=10GB -t 00:60:00 ./homework B /hw6_dataset.fa

	 this can be executed without srun too as:

	 monsoon: ./homework B ./hw6_dataset.fa


Step 06: To check the job status for time:

	 jobstats -j <job_id>