#ifndef FASTA_H
#define FASTA_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <math.h>

using namespace std;
     
struct Node 
{
    char dataLL[16];
    char headLL[16];
	struct Node *next;
};

//main class 
class FASTAreadset_DA
{
    private:
		bool boolDataArr[4294967295];
		ifstream input;
		Node* headLL;
		Node* dataHead;
		char ** header;
		char ** data;
		unsigned int totalCount, value, hashSize, numbOfCollisions = 0;
        string filePath;

    public:
        
        FASTAreadset_DA()
        {
			headLL = NULL;
			dataHead = NULL;
        }    
        
        FASTAreadset_DA(string filepath)
        {
            filePath = filepath;

            if(input.is_open())
			{
				input.close();
			}
			input.open(filePath.c_str());
        }
		
		int tableSize()
		{
			string x;
			int count = 0;
			ifstream input;
			input.open(filePath.c_str()); 

			while(getline(input, x))
			{
				count++;
			}
			return count;
		}

        void readData()
        {
			int count = tableSize();
			data = new char*[count];
			header = new char*[count];
			
			for(int i=0; i < count; i++)
			{
				data[i] = new char[16];
				header[i] = new char[16];
			}

			for(int i=0; i < count; i++)
			{
				Node * newNode = new Node;

				input >> newNode -> headLL;		
				input >> newNode -> dataLL;

				data[i] = newNode->dataLL;

				newNode->next = NULL;

				if(i == 0)
				{
					headLL = newNode;
				}
			}
			input.close();
        }

		void radixFun()
		{
			unsigned int size = tableSize();
			int power = 0;
			int unique = 0;
			
			for (int i = 0; i < size; i++)
			{
				totalCount = 0;
				power = 1;

				for (int j = 0; j < 16; j++)
				{
					if(data[i][j]=='A')
					{
						value = 0;
					}
					else if(data[i][j]=='C')
					{
						value = 1;
					}
					else if(data[i][j]=='G')
					{
						value = 2;
					}
					else if(data[i][j]=='T')
					{
						value = 3;
					}
					totalCount = totalCount + (value * power);
					power = power * 4;
				}

				if(hashSize < totalCount)
				{
					hashSize = totalCount;
				}

				if(boolDataArr[totalCount]==true)
				{
					numbOfCollisions++;
				}
				else
				{
					boolDataArr[totalCount]=true;
				}
			}

			for(int i = 0; i < hashSize; i++)
			{
				
				if(boolDataArr[i] == true)
				{
					unique++;
				}
			}
			
			cout << "Size of mgchar hash table is: " << hashSize << endl;
			cout << "Number of numbOfCollisions observed: " << numbOfCollisions << endl;
			cout << "Number of unique sequences observed: " << unique<< endl;	
		}

	   	void genomicData(string filePath)
		{
			ifstream input;
			char temp[16];
			char gchar = '\0';
			Node* current;
			int count = 0;
			int records = 0;

            if(input.is_open())
			{
				input.close();
			}

			input.open(filePath.c_str());

			while(gchar != '\n')
			{
				input.get(gchar);
			}

			while (input.get(gchar))
			{
  				if(gchar == 'A' || gchar == 'C' || gchar == 'G' || gchar == 'T' || gchar == 'N')
				{
					temp[count] = gchar;
  					count++;

  					if(count == 16)
					{
						count = 0;
  						records++;
  						
  						Node* newNode = new Node;

  						for(int x = 0; x < 16; x++)
						{
			    			(newNode -> dataLL)[x] = temp[x];
						}

			    		(newNode -> dataLL)[16] = '\0';

			    		newNode -> next = NULL;

			    		if(records == 1)
						{
			    			dataHead = newNode;
						}
						else.
						{
							current -> next = newNode;
						}

						current = newNode;
					}

					if(dataHead != NULL)
					{
						for(int i = 1; i < 16; i++)
						{
							temp[i - 1] = temp[i];
						}
					}
				}  
			}         

			cout << "Total number of 16-character fragments observed are: " << records << endl;
			input.close();
		}
	   
		void genomicDataSearchFun()
		{
			unsigned int size = tableSize();
			unsigned int totalCount, value, hashSize, x = 0;
			int power = 0;

			Node* temp = dataHead;
			
			while (temp!=NULL)
			{
				totalCount = 0;
				power = 1;
				for (int j = 0; j < 16; j++)
				{
					if(temp->dataLL[j]=='A')
					{
						value = 0;
					}
					else if(temp->dataLL[j]=='C')
					{
						value = 1;
					}
					else if(temp->dataLL[j]=='G')
					{
						value = 2;
					}
					else if(temp->dataLL[j]=='T')
					{
						value = 3;
					}

					totalCount = totalCount + (value * power);
					power = power * 4;
				}
				
				temp = temp -> next;

				if(hashSize < totalCount)
				{
					hashSize = totalCount;
				}

				if(boolDataArr[totalCount]==true)
				{
					x++;
				}
			}
		}

		void deleteList(Node** head_ref)
        {
 
            Node* current = *head_ref;
            Node* next = NULL;
 
            while (current != NULL) 
            {
                next = current->next;
                free(current);
                current = next;
            }
 
            *head_ref = NULL;
        }
   
        ~FASTAreadset_DA()
        {
            deleteList(&headLL);
			deleteList(&dataHead);
			delete[] header;
			delete[] data;

        }
};

#endif