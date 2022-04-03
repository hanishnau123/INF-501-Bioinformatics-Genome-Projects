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
			cout << "\nPart A: Reading the SARS-COV2 genome sequence and storing it in the suffix trie\n" << endl;
            string copyGenomeString = readGenome(argv[2]);
			copyGenomeString += "$";
			suffix_trie suffixTrie(copyGenomeString);
			
			cout << "Searching sequence 1: TTTTAAGTGTTATGGAGTGTCTCCTACTAAATTAAA \n" << endl;
			suffixTrie.Search("TTTTAAGTGTTATGGAGTGTCTCCTACTAAATTAAA");
			cout << "Searching sequence 2: TCTACCAGTGTCTATGACCAAGACATCAGTAGATTG \n" << endl;
			suffixTrie.Search("TCTACCAGTGTCTATGACCAAGACATCAGTAGATTG");
			cout << "Searching sequence 3: TTTACAAGACTTCAGAGTTTAGATAATGTGGCTTTT \n" << endl;
			suffixTrie.Search("TTTACAAGACTTCAGAGTTTAGATAATGTGGCTTTT");
			cout << "Searching sequence 4: TACCAATTTACCTTTACAGCTAGTTTTTTCTACAGG \n" << endl;
			suffixTrie.Search("TACCAATTTACCTTTACAGCTAGTTTTTTCTACAGG");
			cout << "Searching sequence 5: CCTTACCGCAGAGACAGAAGAAACAGCAAACTGTGA \n" << endl;
            suffixTrie.Search("CCTTACCGCAGAGACAGAAGAAACAGCAAACTGTGA");
			break;
		}	
		case 'B':
		    cout << "Part B: \n" << endl;
			string copyGenomeString = readGenome(argv[2]);
			copyGenomeString += "$";
			suffix_trie suffixTrie(copyGenomeString);

			cout << "Size of the suffix trie is: \n" << endl;
			cout << suffixTrie.printSize() << endl;
			
		    break;
	}
	
	return (0);
}
