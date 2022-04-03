The folder contains the following files:
homework.cpp, homework.h, hw5_dataset.fa
	
Following commands are needed to execute the program:

Step 01: move into the directory

Step 02: make (to make the makefile)
	 monsoon: make

Step 03: To execute part A in monsoon:
 
	 monsoon: srun --mem=10GB -t 00:60:00 ./homework A ./hw5_dataset.fa
	 this can be executed without srun too

	 monsoon: ./homework A ./hw5_dataset.fa	

Step 04: To check the job status for time:

	 jobstats -j <job_id>
 
Step 05: To execute part B in monsoon:
 
	 srun --mem=10GB -t 00:60:00 ./homework B /hw5_dataset.fa

	 this can be executed without srun too as:

	 monsoon: ./homework B ./hw5_dataset.fa


Step 06: To check the job status for time:

	 jobstats -j <job_id>