#ifndef FASTA_CHAIN_H
#define FASTA_CHAIN_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "homework.h"

using namespace std;
     
struct chainNode 
{
    char *head;
    unsigned int chain;
	struct chainNode *next;
};

//main class 
class FASTAreadset_Chain
{
    private:
		char ** header;
		char ** dataChain;
		string filePath;
        ifstream input;
		unsigned int totalChain, valueChain, fragments, collisions = 0;
		int hashtableSize;
		chainNode* headChain;
		chainNode* dataheadChain;
		bool boolArray[4294967295];

    public:

		chainNode ** hashTable;
        FASTAreadset_Chain()
        {
			headChain = NULL;
			dataheadChain = NULL;
        }
        
        FASTAreadset_Chain(string filepath, int userSize)
        {
            filePath = filepath;
			hashtableSize = userSize;

            if(input.is_open())
			{
				input.close();
			}

			input.open(filePath.c_str());
			headChain = NULL;
			dataheadChain = NULL;
        }
		
		int sizeOfTable()
		{
			string x;
			int counter = 0;
			ifstream input;
			input.open(filePath.c_str()); 
	
			while(getline(input,x))
			{
				counter++;
			} 

			return counter;
		}

        void readData()
        {
			int counter = sizeOfTable();
			header = new char*[counter];
			dataChain = new char*[counter];			
					
			for(int i=0;i<counter;i++)
			{
				header[i] = new char[16];
				dataChain[i] = new char[16];
			}

			chainNode * current = new chainNode;

			for(int i=0;i<counter;i++)
			{
				chainNode * newNode = new chainNode;

				input >> header[i];	
				input >> dataChain[i];

				newNode->next = NULL;

				if(current != NULL)
				{
					current->next = newNode;
				}

				current = newNode;

				if(i == 0)
				{
					headChain = newNode;
				}

			}
			input.close();
        }

		void hashTableCreation()
		{
			
			hashTable =  new chainNode*[hashtableSize];
 			
			for(int i = 0 ; i < hashtableSize ;i++)
			{
				hashTable[i] = NULL;
 			}

			unsigned int size = sizeOfTable();
			
			int power = 0;

 			for(int i = 0; i < hashtableSize ; i++)
			{
				power = 1;
				for (int j = 0; j < 16; j++)
				{
					if(dataChain[i][j]=='A')
					{
						valueChain = 0;
					}
					else if(dataChain[i][j]=='C')
					{
						valueChain = 1;
					}
					else if(dataChain[i][j]=='G')
					{
						valueChain = 2;
					}
					else if(dataChain[i][j]=='T')
					{
						valueChain = 3;
					}

					totalChain = totalChain + (valueChain * power);
					power = power * 4;
				}

				totalChain = totalChain % hashtableSize;
				
				chainNode *current = hashTable[totalChain];
				
             	chainNode *temp = NULL;
             	
				while (current != NULL)
				{
					temp = current;
					current = current->next;
				}

				if (current == NULL)
				{
					current = new chainNode();
					current->chain = totalChain;
					current->head = new char[16];

					for(int k = 0 ; k < 16 ; k++)
					{
						current->head[k] = dataChain[i][k];
					}

					if (temp == NULL)
					{
						hashTable[totalChain] = current;
					}

					else
					{
						temp->next = current;
					}

				}
				else
				{
					for(int m = 0 ; m < 16 ; m++)
					{
	                	(current->head)[m] = dataChain[i][m];
					}
				}
			}
		}

	   	void readGenomeData(string filePath)
		{
			ifstream input;
			char temp[16];
			char charChain = '\0';
            if(input.is_open())
			{
				input.close();
			}
			input.open(filePath.c_str());

			int count = 0;
			int rec = 0;

			chainNode* current;
			
			while(charChain != '\n')
			{
				input.get(charChain);
			}

			while (input.get(charChain))
			{
  				if(charChain == 'A' || charChain == 'C' || charChain == 'G' || charChain == 'T' || charChain == 'N')
				{
					temp[count] = charChain;
  					count++;
				
  					if(count == 16)
					{
  						rec ++;
						count = 0;
  						chainNode* newNode = new chainNode;
						newNode ->head = new char[16];

  						for(int k = 0; k < 16; k++)
						{
			    			(newNode -> head)[k] = temp[k];
						}

			    		(newNode -> head)[16] = '\0';

			    		newNode -> next = NULL;

			    		if(rec == 1)
						{
			    			dataheadChain = newNode;
						}
						else
						{
							current -> next = newNode;
						}
						current = newNode;
					}

					if(dataheadChain != NULL)
					{
						for(int i = 1; i < 16; i++)
						{
							temp[i - 1] = temp[i];
						}
					}
				}
			}         
			
			cout << "\nNumber of 16-mer fragments found in the readset: " << rec << endl;
			input.close();
		}
		
		void searchGenomeData()
		{
			int count;
			chainNode* current = dataheadChain;
			chainNode* temp = dataheadChain;
			
			while(current->next != NULL)
			{
				count++;
				current=current->next;
			}

			unsigned int size = count;
			int power = 0;

			while (temp!=NULL)
			{
				totalChain = 0;
				power = 1;

				for (int i = 0; i < 16; i++)
				{
					
					if(hashTable[i] != NULL)
					{
						temp = hashTable[i];

						while(temp->next != NULL)
						{ 
							if(hashTable[i]->chain == i)
							{
								fragments++;
 							 	temp = temp->next;
							}
						}

					}
				}
				
			}
			cout << "\nGenome 16-mer fragments found in read set: " << fragments << endl;
		}

		void deleteChain(chainNode** headRef)
        {
            chainNode* current = *headRef;
            chainNode* next = NULL;
 
            while (current != NULL) 
            {
                next = current->next;
                free(current);
                current = next;
            }
            *headRef = NULL;
        }
		
		void deleteHash()
		{
			for(int i=0;i<hashtableSize;i++)
			{
				chainNode* current = hashTable[i];
				chainNode* temp = new chainNode();
				
				while (current != NULL)
				{
					temp = current;
					free(temp);
					current = current->next;
				}

				hashTable[i]= NULL;
			}
		}
   
        ~FASTAreadset_Chain()
        {
            deleteChain(&headChain);
			deleteChain(&dataheadChain);
			deleteHash();
		}       
};

#endif