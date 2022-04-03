#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string.h>
#include "homework.h"
#include <time.h>

using namespace std;
	
int main(int argc, char ** argv){ 

	if (argc != 3){
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
	}

	char* secondArgument = argv[1];

	switch(*secondArgument)
	{
		case 'A':
		{
			cout << "\nPart A: \n" << endl;
            char* copyGenomeChar = readGenome(argv[2]);
            questionA(argv[2], copyGenomeChar);
			break;
		}	
		case 'B':
		    cout << "\nPart B: \n" << endl;
		    char* copyGenomeChar = readGenome(argv[2]);
		    questionB(argv[2], copyGenomeChar);
		    break;
	}
	
	return (0);
}
