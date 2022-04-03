#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string.h>
#include "homework.h"
#include "homework_chain.h"

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

	int userHashTableSize = 10000000;
	char* argument = argv[1];

	FASTAreadset_DA fasta(argv[2]);
	FASTAreadset_Chain fastaTwo(argv[2], userHashTableSize);

	switch(*argument)
	{
		case 'A':
			fasta.readData();
			fasta.radixFun();
			break;

		case 'B':
			fasta.genomicData(argv[3]);
			fasta.genomicDataSearchFun();
			break;

		case 'C': 
			fastaTwo.readData();
			fastaTwo.hashTableCreation();
			break;

		case 'D':
			fastaTwo.readGenomeData(argv[3]);
			fastaTwo.searchGenomeData();
			break;

		default:
			cout << "\nPlease select another option" << endl;
			break;
	}
	return (0);
}
