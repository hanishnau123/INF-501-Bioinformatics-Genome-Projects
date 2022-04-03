#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string.h>
#include "homework.h"
#include "homework_BLAST.h"
#include <time.h>

using namespace std;
	
int main(int argc, char ** argv){ 

	if (argc != 4){
		cout << endl << endl << "=========================="<< endl;
		cout << "Error: Total 4 input parameters expected" <<endl;
		cout << "Proper usage is:" <<endl;
		cout << "./homework <problem-flag> <filepath01> <filepath02>" << endl;
		cout << "Example:" << endl;
		cout << "./homework A file1.fa file2.fa"  << endl;
		cout << "=========================="<< endl << endl;
		cout << "exiting..." << endl;
		exit(-1);
	}else{

		cout << "\n\nThe number of arguments passed: " << argc << endl;
		cout << "The first argument is: " << argv[0] << endl; 	
		cout << "The second argument is: " << argv[1] << endl; 
		cout << "The third argument is: " << argv[2] << endl;
		cout << "The fourth argument is: " << argv[3] << endl;
	}

	int inputValue = 1000000;
	char* secondArgument = argv[1];

	FASTAreadset_DA fasta(argv[2]);
	FASTAreadset_BLAST fasta_BLAST(argv[2]);

	switch(*secondArgument)
	{
		case 'A': 
			cout << "\nPart 1 A: Smith-Waterman algorithm implementation between two genomic sequences: \n" << endl;
			fasta.readQuery();
			fasta.readGenomicData(argv[3]);
			fasta.swFunction();
			break;

		case 'B': 
			cout << "\nPart 1 B: Performing alignment of queries to the subject sequence: \n" << endl;
			fasta.randomSequenceGenerator(inputValue);
			fasta.readGenomicData(argv[3]);
			fasta.swFunction();
			break;

		case 'C': 
			cout << "\nPart 2 A: Seed based Smith-Waterman: \n" << endl;
			fasta_BLAST.read11mers();
			fasta_BLAST.read11merGenome(argv[3]);
			fasta_BLAST.hashtableCreation(argv[3]);
			fasta_BLAST.searchQuerySeedInHashtable(argv[3]);
			break;
			
		case 'D': 
			cout << "\nPart 2 B: Testing the code, generating completely random 50mers and aligning them to SARS-COV2 data: \n" << endl;
			fasta_BLAST.randomSequenceGenerator(inputValue);
			fasta_BLAST.read11merGenome(argv[3]);
			fasta_BLAST.hashtableCreation(argv[3]);
			fasta_BLAST.searchRandomQuerySeedInHashtable(argv[3]);
			break;

		default:
			cout << "Enter second argument from A to D !!!" << endl;
			break;
	}

	return (0);
}
