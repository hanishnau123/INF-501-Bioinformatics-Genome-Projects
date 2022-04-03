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
	char headLL[50];
    char datatLL[50];
	char genomeDataLL[70];
	struct Node *next;
};

class FASTAreadset_DA 
{
    private:
		int index;
		int gapPenalty = -3;
		int match = 2;
		int mismatch = -1;
        string filePath;
        ifstream input;
		Node* head;
		Node* data;
		
    public:
        
        FASTAreadset_DA()
        {
			head = NULL;
			data = NULL;
        }    
        
        FASTAreadset_DA(string filepath)
        {
            filePath = filepath;

            if(input.is_open())
			{
				input.close();
			}

			input.open(filePath.c_str()); 
			head = NULL;
			data = NULL;
        }
		
		int sizeOfTable(string filePath)
		{
			string inputLine;
			int counter = 0;
			ifstream input;

			input.open(filePath.c_str()); 

			while(getline(input, inputLine))
			{
				counter++;
			} 

			return counter;
		}

        void readQuery()
        {
			
			int counter = sizeOfTable(filePath);
			Node * currentNode = new Node;

			cout << "Number of reads/query lines: " << counter << endl;

			for(int i=0; i<counter/2; i++)
			{
				Node * newNode = new Node;

				input >> newNode->headLL;	
				input >> newNode->datatLL;

				newNode->next = NULL;

				if(currentNode != NULL)
				{
					currentNode->next = newNode;
				}

				currentNode = newNode;

				if(i == 0)
				{
					head = newNode;
				}

			}
			input.close();
        }

	   	void readGenomicData(string filepath)	
		{
			int count = 0;
			int recordNumb = 0;
			char temp[70];
			char genomeChar = '\0';
			Node* current;
			ifstream input2;

            if(input2.is_open()) 
			{
				input2.close();
			}

			input2.open(filepath.c_str());
			
			while(genomeChar != '\n')
			{
				input2.get(genomeChar);
			}

			while(input2.get(genomeChar))
			{
  				if(genomeChar == 'A' || genomeChar == 'C' || genomeChar == 'G' || genomeChar == 'T' || genomeChar == 'N')
				{
					temp[count] = genomeChar;
  					count++;
				
  					if(count == 70)
					{
						count--;
  						recordNumb++;
  						Node* newNode = new Node;

  						for(int k = 0; k < 70; k++)
						{
			    			(newNode -> genomeDataLL)[k] = temp[k];
						}

			    		(newNode -> genomeDataLL)[70] = '\0';

			    		newNode -> next = NULL;
						count = 0;
			    		
			    		if(recordNumb == 1)
						{
			    			data = newNode;
						}
						else
						{
							current -> next = newNode;
						}

						current = newNode;
					}

					if(data != NULL)
					{
						for(int i = 1; i < 70; i++)
						{
							temp[i - 1] = temp[i];
						}
					}
				}  
			}         

			cout << "\nTotal number of Genomers: " << recordNumb << endl;
			input2.close();
		}

		// random sequence generator
		void randomSequenceGenerator(int userVal)
		{
			Node * currentNode = new Node;
			srand(time(0));

			for (int i = 0; i < userVal; i++)
			{
				Node * newNode = new Node;
				
				for (int j = 0; j < 50; j++)
				{
					int result = ( rand() % 4 );

					if(result==0)
					{
						newNode->datatLL[j] = 'A';
					}
					else if(result==1)
					{
						newNode->datatLL[j] = 'C';
					}
					else if(result==2)
					{
						newNode->datatLL[j] = 'G';
					}
					else if(result==3)
					{
						newNode->datatLL[j] = 'T';
					}

				}

				newNode->next = NULL;

				if(currentNode != NULL)
				{
					currentNode->next = newNode;
				}

				currentNode = newNode;

				if(i == 0)
				{
					head = newNode;
				}
			}
		}

		double similarityScore(char a, char b) 
		{
			double result;
			if(a==b)
			{
				result=match;
			}
			else
			{
				result=mismatch;
			}
			return result;
		}

		double findMax(double array[], int length)
		{
			double max = array[0];
			index = 0;

			for(int i=1; i<length; i++)
			{
				if(array[i] > max)
				{
					max = array[i];
					index=i;
				}
			}
			return max;
		}

		// Smith Waterman Algorithm
		void swFunction()
		{
			string sequenceA;
			string sequenceB;
			
			Node* tempRead = head;
			Node* tempGenome = data;

			while(tempRead != NULL)
			{
				sequenceA = tempRead -> datatLL;

				while(tempGenome != NULL)
				{					
					sequenceB = tempGenome -> genomeDataLL;
					tempGenome = tempGenome -> next;

					cout<<"Sequence 01 and Sequence 02: \n"<<endl;
					cout << sequenceA << endl;
					cout << sequenceB << endl;
			
					int lengthOfSeqA = sequenceA.length();
					int lengthOfSeqB = sequenceB.length();
					double matrix[lengthOfSeqA+1][lengthOfSeqB+1];

					for(int i=0;i<=lengthOfSeqA;i++)
					{
						for(int j=0;j<=lengthOfSeqB;j++)
						{
							matrix[i][j]=0;
						}
					}
					
					double tracebackArr[4];
					int I_i[lengthOfSeqA+1][lengthOfSeqB+1];
					int I_j[lengthOfSeqA+1][lengthOfSeqB+1];

					for (int i=1;i<=lengthOfSeqA;i++)
					{
						for(int j=0;j<=lengthOfSeqB;j++)
						{
							tracebackArr[0] = matrix[i-1][j-1] + similarityScore(sequenceA[i-1],sequenceB[j]);
							tracebackArr[1] = matrix[i-1][j]+gapPenalty;
							tracebackArr[2] = matrix[i][j-1]+gapPenalty;
							tracebackArr[3] = 0;
							matrix[i][j] = findMax(tracebackArr,4);
							
							if(index==0)
							{							
								I_i[i][j] = i-1;
								I_j[i][j] = j-1;
							}
							else if(index==1)
							{	
								I_i[i][j] = i-1;
								I_j[i][j] = j;
							}
							else if(index==2)
							{
								I_i[i][j] = i;
								I_j[i][j] = j-1;
							}
							else if(index==3)
							{
								I_i[i][j] = i;
								I_j[i][j] = j;
							}
							
						}
					}
					cout<<endl;

					double matrixMaxValue = 0;
					int iMax=0, jMax=0;

					for(int i=1;i<lengthOfSeqA;i++)
					{
						for(int j=1;j<lengthOfSeqB;j++)
						{
							if(matrix[i][j]>matrixMaxValue)
							{
								matrixMaxValue = matrix[i][j];
								iMax=i;
								jMax=j;
							}
						}
					}

					cout << "Max score in the matrix is: " << matrixMaxValue << endl;

					int current_i_value=iMax, current_j_value=jMax; 

					int next_i_value = I_i[current_i_value][current_j_value];
					int next_j_value = I_j[current_i_value][current_j_value];

					int tick=0;

					char consensus_a_array[lengthOfSeqA + lengthOfSeqB + 2];
					char consensus_b_array[lengthOfSeqA + lengthOfSeqB + 2];

					while(((current_i_value!=next_i_value) || (current_j_value!=next_j_value)) && (next_j_value!=0) && (next_i_value!=0))
					{
						if(next_i_value == current_i_value)
						{
							consensus_a_array[tick] = '_';                  	
						}  
						else
						{
							consensus_a_array[tick] = sequenceA[current_i_value - 1]; 

						}

						if(next_j_value==current_j_value)  
						{
							consensus_b_array[tick] = '_';                 
						}
						else
						{
		                   consensus_b_array[tick] = sequenceB[current_j_value - 1]; 
						}   

						current_i_value = next_i_value;
						current_j_value = next_j_value;
						next_i_value = I_i[current_i_value][current_j_value];
						next_j_value = I_j[current_i_value][current_j_value];
						tick++;
					}

					cout<<"\nAlignment of sequence 01 and sequence 02: \n"<<endl;
					
					string markings[lengthOfSeqA+lengthOfSeqB+2];

					for(int i=tick-1; i>=0; i--){
						cout<<consensus_a_array[i]; 
					}
					cout<<endl;

					for(int i=tick-1; i>=0; i--)
					{
						if(consensus_a_array[i]==consensus_b_array[i])
						{
							markings[i]="|";
						}
						else if (consensus_a_array[i]== '_' || consensus_b_array[i] == '_')
						{
							markings[i]=" ";
						}
						else if (consensus_a_array[i] != consensus_b_array[i])
						{
							markings[i]="x";
						}
					}

					for(int i=tick-1; i>=0; i--){
						cout<<markings[i]; 
					}
					cout<<endl; 

					for(int j=tick-1; j>=0; j--){
						cout<<consensus_b_array[j]; 
					}
					cout<<endl;
				}
				tempRead = tempRead -> next;
			}
		}

		void deleteList(Node** headRef)
        {
            Node* current = *headRef;
            Node* next = NULL;
 
            while (current != NULL) 
            {
                next = current->next;
                free(current);
                current = next;
            }
            *headRef = NULL;
        }
   
        ~FASTAreadset_DA()
        {
            deleteList(&head);
        }
        
};

#endif
