#ifndef FASTA_H
#define FASTA_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctime>
#include <cstdlib>

using namespace std;
     
struct Node 
{
	struct Node *parent;
	struct Node *children[4];
	int index;
};

struct stackValues{
    int mismatch = 0;
    Node* StackNode = nullptr;
};

class Stack{
	public:
	int top;
    stackValues Stackvalue[100000];

	Stack()
	{
		top = -1;
	}

	bool push(stackValues Val)
	{
		if (top >= 100000 - 1)
		{
			cout << "Stack Overflow";
			return false;
		}
		else
		{
			Stackvalue[++top] = Val;
			return true;
		}
	}

	stackValues pop()
	{
		if (top < 0)
		{
			cout << "Stack Underflow";
			stackValues nullStruct;
			nullStruct.mismatch = 0;
			nullStruct.StackNode = nullptr;
			return nullStruct;
		}
		else
		{
			stackValues val = Stackvalue[top--];
			return val;
		}
	}

	bool isEmpty()
	{
		return (top < 0);
	}

};

class prefix_trie 
{
    private:
        string filePath;
        ifstream input;
		
    public:
		Node* root;
		string path;
		char* genomeChar;
        
        prefix_trie()
        {
			root = (Node*)calloc(1, sizeof(Node));
        }    
        
        prefix_trie(string filepath, char* Genome)
        {
			root = (Node*)calloc(1, sizeof(Node));
			this->path = filepath;
			this->genomeChar = Genome;
        }

		prefix_trie(prefix_trie &obj)
		{
			this->root = obj.root;
    		this->path = obj.path;
		}

		int indexReturn(char charVal)
		{
			int returnValue;
			switch (charVal)
			{
				case 'A':
					returnValue = 0;
					break;
				case 'C':
					returnValue = 1;
					break;
				case 'G':
					returnValue = 2;
					break;
				case 'T':
					returnValue = 3;
					break;
				default:
					returnValue = -1;
					break;
			}
			return returnValue;
		}

		char returnCharFunction(int index)
		{
			char returnChar;
			switch (index) {
				case 0:
					returnChar = 'A';
					break;
				case 1:
					returnChar = 'C';
					break;
				case 2:
					returnChar = 'G';
					break;
				case 3:
					returnChar = 'T';
					break;
				default:
					returnChar = '\0';
					break;
			}
			return returnChar;
		}

		void insertFunc(char *seq)
		{
			this->root->index = -1;
			Node* temp = this -> root;
			int length = strlen(seq);

			for (int i = 0; i < length; i++)
			{
				if (temp->children[indexReturn(seq[i])] == nullptr)
				{
					temp->children[indexReturn(seq[i])] = (Node*)calloc(1, sizeof(Node));
					temp->children[indexReturn(seq[i])]->index = i;
					temp->children[indexReturn(seq[i])]->parent = temp;
				}
				temp = temp->children[indexReturn(seq[i])];
			}
		}

		int prefixTrieSize(Node* Node)
		{
			if (Node == nullptr)
				return 0;
			else
				return prefixTrieSize(Node->children[0]) + prefixTrieSize(Node->children[1]) + prefixTrieSize(Node->children[2]) + prefixTrieSize(Node->children[3]) + 1;
		}
		
		int fuzzyMatching(Node* node, char* query, int mismatch)
		{
			int queryLength = (int)strlen(query);
			int mismatchValue = 0;
			Stack stack;

			stackValues stackstruct;
			stackstruct.mismatch = mismatchValue;
			stackstruct.StackNode = node;
			stack.push(stackstruct);

			while (stack.top >= 0)
			{
				stackValues tempStack = stack.pop();
				Node* TempNode = tempStack.StackNode;
				mismatchValue = tempStack.mismatch;

				int i = TempNode->index + (int)1;
				while(i < queryLength)
				{
					if (TempNode->children[indexReturn(query[i])] == nullptr)
					{
						mismatchValue++;
						if (mismatchValue > mismatch)
							break;
						for (int j = 0; j < 4; j++)
						{
							if (TempNode->children[j] != nullptr)
							{
								stackValues nextLeaf;
								nextLeaf.mismatch  = mismatchValue;
								nextLeaf.StackNode = TempNode->children[j];
								stack.push(nextLeaf);
							}
						}
						tempStack = stack.pop();
						TempNode = tempStack.StackNode;
						mismatchValue = tempStack.mismatch;
					}
					else
						TempNode = TempNode->children[indexReturn(query[i])];
					i++;
				}
			}

			if (mismatchValue == 0)
				return 0;
			else if (mismatchValue == 1)
				return 1;
			else
				return 2;
		}
   
        ~prefix_trie()
        {

        }
};

char* readGenome(string Path)
{
	ifstream myFile(Path);

	if (!myFile)
	{
		cout << "Could not open the file" << endl;
		exit(1);
	}
	else
	{
		string skipString;
		char charValue;
		int index  = 0;
		getline(myFile, skipString);
		int charSize = 0;
		while (myFile.get(charValue))
		{
			if ((charValue == '\n') || (charValue == '\r'))
			{
				continue;
			}

			charSize++;
		}
		myFile.clear();
		myFile.seekg(0, ios::beg);
		char* copyChar;
		copyChar = new char[charSize];
		getline(myFile,skipString);
		while (myFile.get(charValue))
		{
			if ((charValue == '\n') || (charValue == '\r'))
				continue;
			copyChar[index] = charValue;
			index++;
		}
		copyChar[index] = '\0';
		if (myFile.is_open())
		{
			myFile.close();
		}
		return copyChar;
	}
}

char** randomGenerator(int numb, char* genome, float errorValue)
{
	char charA[3] = {'C', 'G', 'T'};
	char charC[3] = {'A', 'G', 'T'};
	char charG[3] = {'A', 'C', 'T'};
	char charT[3] = {'A', 'C', 'G'};
	char** returnVal;
	returnVal = new char*[numb];
	for (int i = 0; i < numb; i++)
	{
		returnVal[i] = new char[36 + 1];
		int randomIndex = rand() % ((int)strlen(genome) - 1 - 36 + 1);
		for(int j = 0; j < 36; j++)
		{
			returnVal[i][j] = genome[j + randomIndex];
			if (errorValue > 0)
			{
				float prob = ((float)rand() / (RAND_MAX));
				if (prob < errorValue)
				{
					int Rand;
					switch(returnVal[i][j])
					{
						case 'A':
							Rand = rand() % (3);
							returnVal[i][j] = charA[Rand];
							break;
						case 'C':
							Rand = rand() % (3);
							returnVal[i][j] = charC[Rand];
							break;
						case 'G':
							Rand = rand() % (3);
							returnVal[i][j] = charG[Rand];
							break;
						case 'T':
							Rand = rand() % (3);
							returnVal[i][j] = charT[Rand];
							break;
						default:
						break;
					}
				}
			}
		}
		returnVal[i][36] = '\0';
	}
	return returnVal;
}

void questionA(string path, char* genome)
{
	int randomTargets[3] = {5000, 50000, 100000};
	int mismatchNumber;
	int mismatchFound[3] = {0, 0, 0};
	int totalMatches = 0;

	for (int j = 0; j < 3; j++)
	{
		prefix_trie Trie(path, genome);
		cout << "Added " << randomTargets[j] << " Random 36mers to the Prefix-trie\n" << endl;
		char** randomTarg = randomGenerator(randomTargets[j], Trie.genomeChar, 0);
		for (int i = 0; i < (int)randomTargets[j]; i++)
		{
			Trie.insertFunc(randomTarg[i]);
		}
		cout << "Size of the Trie = " << Trie.prefixTrieSize(Trie.root) << endl;
		for(int k = 0; k < (int)strlen(genome) - 36 + 1; k++)
		{
			char Query[36 + 1];
			for (int m = 0; m < 36; m++)
			{
				Query[m] = genome[m + k];
			}
			Query[36] = '\0';
			mismatchNumber = Trie.fuzzyMatching(Trie.root, Query, 1);

			if (mismatchNumber == 0)
				mismatchFound[0]++;
			else if (mismatchNumber == 1)
				mismatchFound[1]++;
			else
				mismatchFound[2]++;
		}

		cout << "\nNumber of 0 Mismatches = " << mismatchFound[0] << "\n" <<
				"Number of 1 Mismatches = " << mismatchFound[1] << "\n" <<
				"Number of 2 or more Mismatches = " << mismatchFound[2] << "\n" << endl;
        totalMatches = mismatchFound[0] + mismatchFound[1] + mismatchFound[2];
		cout << "Number of Matches = " << totalMatches << endl;
		cout << "****************************************" << endl;
		cout << "\n";
	}
}

void questionB(string path, char* genome)
{
	int randomTargets[3] = {5000, 50000, 100000};
	int mismatchNumber;
	int mismatchFound[3] = {0, 0, 0};
	int totalMatches = 0;

	for (int j = 0; j < 3; j++)
	{
		prefix_trie Trie(path, genome);
		cout << "Added " << randomTargets[j] << " Random 36mers to the Prefix-trie\n" << endl;
		char** randomTarg = randomGenerator(randomTargets[j], Trie.genomeChar, 0.01);
		for (int i = 0; i < (int)randomTargets[j]; i++)
		{
			Trie.insertFunc(randomTarg[i]);
		}
		cout << "Size of the Trie = " << Trie.prefixTrieSize(Trie.root) << endl;
		for(int k = 0; k < (int)strlen(genome) - 36 + 1; k++)
		{
			char queryArr[36 + 1];
			for (int m = 0; m < 36; m++)
			{
				queryArr[m] = genome[m + k];
			}
			queryArr[36] = '\0';
			mismatchNumber = Trie.fuzzyMatching(Trie.root, queryArr, 1);

			if (mismatchNumber == 0)
				mismatchFound[0]++;
			else if (mismatchNumber == 1)
				mismatchFound[1]++;
			else
				mismatchFound[2]++;
		}

		cout << "\nNumber of 0 Mismatches = " << mismatchFound[0] << "\n" << "Number of 1 Mismatches = " << mismatchFound[1] << "\n" << "Number of 2 or more Mismatches = " << mismatchFound[2] << "\n" << endl;
        totalMatches = mismatchFound[0] + mismatchFound[1] + mismatchFound[2];
		cout << "Number of Matches = " << totalMatches << endl;
		cout << "****************************************" << endl;
		cout << "\n";
	}
}
#endif
