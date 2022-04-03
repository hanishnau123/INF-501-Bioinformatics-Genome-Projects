#ifndef FASTA_BLAST_H
#define FASTA_BLAST_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <cstring>
#include <math.h>
#include "homework.h"

using namespace std;

struct BLASTNode 
{
    char dataLL[11];
    char randomDataCharArray[50];
	char genomeData[11];
	char *hashData;
    unsigned int hashHeader;
	struct BLASTNode *BLASTNext;
};
 
class FASTAreadset_BLAST
{
    private:
		
		char ** headerArray;
		char ** dataArray;
		char ** dataArrayGenome;
		char genomeHeadArray[38];

		unsigned int totalVal = 0;
		unsigned int powerVal = 0;
		unsigned int value = 0;
		unsigned int hashtableSize = 0;
		unsigned int fragmentPresent = 0;
		unsigned int collisions = 0;
		int usersizeofhashtable;

		int index;
		double gapPenalty = -3;
		double match = 2;
		double mismatch = -1;

		bool boolArray[4200000];

        string filePath;
        string filePath02;
        ifstream inputFile;
        BLASTNode* headBLAST;
		BLASTNode* dataBLAST;

    public:

		BLASTNode **hashTable;

        FASTAreadset_BLAST()
        {
			headBLAST = NULL;
			dataBLAST = NULL;
        }    
        
        FASTAreadset_BLAST(string filepath)
        {
            filePath = filepath;

            if(inputFile.is_open())
			{
				inputFile.close();
			}

			inputFile.open(filePath.c_str()); 

			headBLAST = NULL;
			dataBLAST = NULL;
        }
        
        double similarityScore(char a, char b)
		{
			double result;

			if(a==b){
				result = match;
			}else{
				result = mismatch;
			}

			return result;
		}

		double findMax(double array[], int length)
		{
			index = 0;
			double max = array[0];

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

		void swFunction(string query, string genome)
		{
			string sequenceA;
			string sequenceB;
			sequenceA = query;
			sequenceB = genome;
					
			cout<<"Sequence 01 and Sequence 02: \n"<<endl;
			cout << sequenceA << endl;
			cout << sequenceB << endl;
			
			int sequenceALength = sequenceA.length();
			int sequenceBLength = sequenceB.length();

			double matrix[sequenceALength+1][sequenceBLength+1];

			for(int i=0;i<=sequenceALength;i++)
			{
				for(int j=0;j<=sequenceBLength;j++)
				{
					matrix[i][j]=0;
				}
			}

			double tracebackArray[4];
			int I_i[sequenceALength+1][sequenceBLength+1];
			int I_j[sequenceALength+1][sequenceBLength+1];

			for(int i=1;i<=sequenceALength;i++)
			{
				for(int j=0;j<=sequenceBLength;j++)
				{
					tracebackArray[0] = matrix[i-1][j-1] + similarityScore(sequenceA[i-1],sequenceB[j]);
					tracebackArray[1] = matrix[i-1][j]+gapPenalty;
					tracebackArray[2] = matrix[i][j-1]+gapPenalty;
					tracebackArray[3] = 0;
					matrix[i][j] = findMax(tracebackArray,4);
					
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

			double maxMatrix = 0;
			int i_max_value = 0, j_max_value = 0;

			for(int i=1; i<sequenceALength; i++)
			{
				for(int j=1; j<sequenceBLength; j++)
				{
					if(matrix[i][j]>maxMatrix)
					{
						maxMatrix = matrix[i][j];
						i_max_value=i;
						j_max_value=j;
					}
				}
			}

			cout << "\nMax matrix score is: " << maxMatrix << endl;
					
			int current_i_value = i_max_value, current_j_value = j_max_value; 
			int next_i = I_i[current_i_value][current_j_value];
			int next_j = I_j[current_i_value][current_j_value];
			int tick=0;
			char consensusA[sequenceALength + sequenceBLength + 2];
			char consensusB[sequenceALength + sequenceBLength + 2];

			while(((current_i_value != next_i) || (current_j_value != next_j)) && (next_j != 0) && (next_i != 0))
			{

				if(next_i==current_i_value)
				{
					consensusA[tick] = '_';                  			
				}  
				else
				{
					consensusA[tick] = sequenceA[current_i_value - 1]; 
				}

				if(next_j==current_j_value)  
				{
					consensusB[tick] = '_';                 
				}
				else
				{
					consensusB[tick] = sequenceB[current_j_value - 1]; 
				}   

				current_i_value = next_i;
				current_j_value = next_j;
				next_i = I_i[current_i_value][current_j_value];
				next_j = I_j[current_i_value][current_j_value];
				tick++;
			}

			cout<<"\nAlignment of sequence 01 and sequence 02: \n"<<endl;

			for(int i=tick-1; i>=0; i--){
				cout<<consensusA[i]; 
				cout<<endl;
			} 

			string markingsArray[sequenceALength+sequenceBLength+2];

			for(int i=tick-1; i>=0; i--)
			{
				if(consensusA[i] == consensusB[i])
				{
					markingsArray[i]="|";
				}
				else if(consensusA[i]== '_' || consensusB[i] == '_')
				{
					markingsArray[i]=" ";
				}
				else if (consensusA[i] != consensusB[i])
				{
					markingsArray[i]="x";
				}		
			}
					
			for(int i=tick-1; i>=0; i--){
				cout<<markingsArray[i]; 
				cout<<endl;

			} 

			for(int j=tick-1;j>=0;j--){
				cout<<consensusB[j]; cout<<endl;
			}
		}
		
		int sizeOfTable()
		{
			string inputLine;
			int counter = 0;
			ifstream inputFile;

			inputFile.open(filePath.c_str()); 
					
			while(getline(inputFile,inputLine))
			{
				counter++;
			} 

			return counter;
		}

        void read11mers() 
        {
			int BLASTCounter = sizeOfTable();
				
			cout << "Number of lines in read dataset: " << BLASTCounter << endl;

			headerArray = new char*[BLASTCounter/2];
			dataArray = new char*[BLASTCounter/2];	

			for(int i=0; i<BLASTCounter/2; i++)
			{
				headerArray[i] = new char[50];
				dataArray[i] = new char[50];
			}

			for(int i=0; i<BLASTCounter/2; i++)
			{
				inputFile >> headerArray[i];		
				inputFile >> dataArray[i];
			}

			char querySeed = '\0';
			char temp[50];

			int headerCounter = 0;

			BLASTNode * currentBLAST = new BLASTNode;

			for (int i = 0; i < BLASTCounter/2; i++)
			{
				int count = 0;

				for (int m = 0; m < 41; m++)
				{
				
					for (int j = m; j <= (m + 11); j++)
					{
						BLASTNode * newNodeBLAST = new BLASTNode;
						querySeed = dataArray[i][j];	
						
			  				if(querySeed == 'A' || querySeed == 'C' || querySeed == 'G' || querySeed == 'T')
							{
								temp[count] = querySeed;
								count++;
			  					
			  					if(count == 12)
								{

									for (int k = 0; k < 11; k++)
									{
										(newNodeBLAST->dataLL)[k] = temp[k];
									}
									
									newNodeBLAST->BLASTNext = NULL;

									if(currentBLAST != NULL)
									{
										currentBLAST->BLASTNext = newNodeBLAST;
									}

									currentBLAST = newNodeBLAST;

									if(headerCounter == 0)
									{
										headBLAST = newNodeBLAST;
									}
									headerCounter++;
									count = 0;
									
								}
								
							}
					}
					
				}
				
        	}
        }

        void displayReadQueries()
		{
			BLASTNode* temp;
			temp = headBLAST;
			while(temp != NULL)
			{
				cout << temp -> dataLL << "\n";
				temp = temp -> BLASTNext;
			}
		}

		int radixFunction(string filepath02)
		{
			unsigned int size = sizeOfGenomeTable(filepath02);
			int powerVal = 0;
			
			BLASTNode* current = dataBLAST;

			while(current != NULL)
			{
				totalVal = 0;
				powerVal = 1;

					for (int j = 0; j < 11; j++)
					{

						if((current -> genomeData)[j]=='A')
						{
							value = 0;
						}

						if((current -> genomeData)[j]=='C')
						{
							value = 1;
						}
						
						if((current -> genomeData)[j]=='G')
						{
							value = 2;
						}

						if((current -> genomeData)[j]=='T')
						{
							value = 3;
						}

						totalVal = totalVal + (value * powerVal);
						powerVal = powerVal * 4;
					}
					
					current = current -> BLASTNext;

					if(hashtableSize > totalVal)
					{
						//do nothing
					}
					else
					{
						hashtableSize = totalVal;
					}

					if(boolArray[totalVal]==true)
					{
						collisions++;
					}
					else
					{
						boolArray[totalVal]=true;

					}
			}
			
			int uniqueSeq = 0;

			for(int i = 0; i < hashtableSize; i++)
			{
				
				if(boolArray[i] == true)
				{
					uniqueSeq++;
				}
			}
			
			cout << "\nSize of the HASH table is:" << hashtableSize << endl;
			cout << "\nTotal number of collisions is:" << collisions << endl;
			cout << "\nTotal number of UNIQUE sequences is: " << uniqueSeq<< endl;
			return hashtableSize;
		}
	
		int sizeOfGenomeTable(string filepath_genome)
		{
			string inputLine;
			int counter = 0;
			ifstream inputFile;
			filePath02 = filepath_genome;
			inputFile.open(filePath02.c_str()); 

			while(getline(inputFile,inputLine))
			{
				counter++;
			} 
			
			return counter;
		}

	   	void read11merGenome(string filepath2)
		{
			ifstream input;
			int BLASTCounter = sizeOfGenomeTable(filepath2);
			dataArrayGenome = new char * [BLASTCounter];
			
			cout << "Number of lines in read genome dataset: " << BLASTCounter << endl;


            if(input.is_open()) 
			{
				input.close();
			}

			input.open(filepath2.c_str());

			for (int i = 0; i < BLASTCounter-1; ++i)
			{
				dataArrayGenome[i] = new char[70];
			}

			input >> genomeHeadArray;
			
			for(int i=0;i<BLASTCounter-1;i++)
			{
				input >> dataArrayGenome[i];
			}

			char querySeed = '\0';
			
			char temp[70];

			int headerCounter = 0;

			BLASTNode * currentBLAST = new BLASTNode;

			for (int i = 0; i < BLASTCounter-1; i++)
			{
				int count = 0;

				for (int m = 0; m < 61; m++)
				{
				
					for (int j = m; j <= (m + 11); j++)
					{
						BLASTNode * newNodeBLAST = new BLASTNode;
						querySeed = dataArrayGenome[i][j];	
						
			  				if(querySeed == 'A' || querySeed == 'C' || querySeed == 'G' || querySeed == 'T')
							{
								temp[count] = querySeed;
								count++;
			  					
			  					if(count == 12)
								{
									for (int k = 0; k < 11; k++)
									{
										(newNodeBLAST->genomeData)[k] = temp[k];
									}

									newNodeBLAST->BLASTNext = NULL;

									if(currentBLAST != NULL)
									{
										currentBLAST->BLASTNext = newNodeBLAST;
									}

									currentBLAST = newNodeBLAST;

									if(headerCounter == 0)
									{
										dataBLAST = newNodeBLAST;
									}
									headerCounter++;
									count = 0;
								}
							}
					}
				}
        	}
			input.close();
		}
		
		void hashtableCreation(string filepath2)
		{
		    int hashtableReturnValue = radixFunction(filepath2);
			hashTable =  new BLASTNode*[hashtableReturnValue]; //size is calculated from highest radix value
 			
			for(int i = 0 ; i < hashtableReturnValue ;i++)
			{
				hashTable[i] = NULL;
 			}

			BLASTNode* current = dataBLAST;

     		while(current != NULL)
			{
				totalVal = 0;
				powerVal = 1;

					for (int j = 0; j < 11; j++)
					{
						if((current -> genomeData)[j]=='A')
						{
							value = 0;
						}

						if((current -> genomeData)[j]=='C')
						{
							value = 1;
						}
						
						if((current -> genomeData)[j]=='G')
						{
							value = 2;
						}

						if((current -> genomeData)[j]=='T')
						{
							value = 3;
						}

						totalVal = totalVal + (value * powerVal);
						powerVal = powerVal * 4;
					}
									
				BLASTNode *current02 = hashTable[totalVal];
             	BLASTNode *tempChain = NULL;
             	
				while (current02 != NULL)
				{
					tempChain = current02;
					current02 = current02->BLASTNext;
				}

				if (current02 == NULL)
				{
					current02 = new BLASTNode();
					current02->hashHeader = totalVal;
					current02->hashData = new char[11];
					
					for(int k = 0 ; k < 11 ; k++)
					{
						current02->hashData[k] = (current -> genomeData)[k];
						
					}

					if (tempChain == NULL)
					{
						hashTable[totalVal] = current02;
					}else{
						tempChain->BLASTNext = current02;
					}
				}
				else
				{
					for(int m = 0 ; m < 11 ; m++)
					{
	                	(current02->hashData)[m] = current -> genomeData[m];
					}
				}
				current = current -> BLASTNext;
			}
		}
		
		void searchQuerySeedInHashtable(string filepath2)
		{
			unsigned int size = sizeOfTable();
			int BLASTCounter = sizeOfGenomeTable(filepath2);
			unsigned int totalVal, value, hashsize, fragmentPresent = 0;
			int powerVal = 0;

			BLASTNode* temp = headBLAST;
			
			while (temp!=NULL)
			{
				totalVal = 0;
				powerVal = 1;
				for (int j = 0; j < 11; j++)
				{
					if(temp->dataLL[j]=='A')
					{
						value = 0;
					}

					if(temp->dataLL[j]=='C')
					{
						value = 1;
					}
					
					if(temp->dataLL[j]=='G')
					{
						value = 2;
					}

					if(temp->dataLL[j]=='T')
					{
						value = 3;
					}

					totalVal = totalVal + (value * powerVal);
					powerVal = powerVal * 4;
				}

				if(hashsize > totalVal)
				{
					//do nothing
				}
				else
				{
					hashsize = totalVal;
				}
				
				if(boolArray[totalVal]==true)
				{
					fragmentPresent++;
					
					BLASTNode* current02 = hashTable[totalVal];
					char dataseq[50];
					int counter = 0;
					int val = 1;
					while(current02 != NULL)
					{

    					if(strcmp(temp -> dataLL, current02 -> hashData)==0)
    					{
    					    BLASTNode* current03 = dataBLAST;

    					    while(current03 != NULL)
    					    {
    					        if(strcmp(temp -> dataLL, current03 -> genomeData)==0)
    					        {

    					                for(int l = 0; l< 11 ; l++)
                					    {
                					        dataseq[l] = (temp -> dataLL)[l];
                					    }
                					    
    					                for(int z = 11; z < 50; z++)
    					                {
    					                    if(current03 -> BLASTNext != NULL)
    					                    {
    					                        current03 = current03 -> BLASTNext;
    					                        dataseq[z] = (current03 -> genomeData)[0];
    					                    } 

    					                }
    					                
    					            }
    					          
    					        current03 = current03 -> BLASTNext;
    					    }
    					}
						   
    				swFunction(temp -> dataLL, dataseq);	    
					current02 = current02->BLASTNext;
					}
				}
				else
				{
					//do nothing
				}
				temp = temp -> BLASTNext;
			}
		}
		
		void randomSequenceGenerator(int userValue)
		{
			BLASTNode * current_ll = new BLASTNode;
			srand(time(0));

			for (int i = 0; i < userValue; i++)
			{
				BLASTNode * newNode = new BLASTNode;
				
				for (int j = 0; j < 50; j++)
				{
					int result = ( rand() % 4 );

					if(result==0)
					{
						newNode->randomDataCharArray[j] = 'A';
					}
					else if(result==1)
					{
						newNode->randomDataCharArray[j] = 'C';
					}
					else if(result==2)
					{
						newNode->randomDataCharArray[j] = 'G';
					}
					else if(result==3)
					{
						newNode->randomDataCharArray[j] = 'T';
					}

				}
				
				newNode->BLASTNext = NULL;

				if(current_ll != NULL)
				{
					current_ll->BLASTNext = newNode;
				}

				current_ll = newNode;

				if(i == 0)
				{
					headBLAST = newNode;
				}
			}
		}
		
		void searchRandomQuerySeedInHashtable(string filePath)
		{
			unsigned int size = sizeOfTable();
			int BLASTCounter = sizeOfGenomeTable(filePath);
			unsigned int totalVal,value,hashsize,fragmentPresent = 0;
			int powerVal = 0;

			BLASTNode* temp = dataBLAST;
			BLASTNode* randomTemp = headBLAST;
			int counter = 0;
			char dataSeq[50];
			
			
			while (randomTemp != NULL)
			{
			    	while (temp!=NULL)
			        {
			            for(int l = 0; l< 11 ; l++)
					    {
					        dataSeq[l] = (temp -> genomeData)[l];
					    }
    					 for(int z = 11; z < 50; z++)
		                {
		                    temp = temp -> BLASTNext;
		                    dataSeq[z] = (temp -> genomeData)[0];
		                    
		                   
		                }
				
			            swFunction(randomTemp -> randomDataCharArray, dataSeq);
    			
				        temp = temp -> BLASTNext;
				        counter++;
			        }
			        randomTemp = randomTemp -> BLASTNext;
			}
			cout << "\nNumber of 11mer fragments found: " << fragmentPresent << endl;
		}
		
		void printGenomicDataBLAST()
		{
			BLASTNode* current = dataBLAST;
			while(current != NULL)
			{
				cout << current -> genomeData << "\n";
				current = current -> BLASTNext;
			}
		}
	
		void deleteHash(BLASTNode** headRef)
        {
 
            BLASTNode* current = *headRef;
            BLASTNode* BLASTNext = NULL;
 
            while (current != NULL) 
            {
                BLASTNext = current->BLASTNext;
                free(current);
                current = BLASTNext;
            }
 
            *headRef = NULL;
        }
        
        void deleteListBLAST(BLASTNode** headRef)
        {
            BLASTNode* current = *headRef;
            BLASTNode* BLASTNext = NULL;
 
            while (current != NULL) 
            {
                BLASTNext = current->BLASTNext;
                free(current);
                current = BLASTNext;
            }
 
            *headRef = NULL;
        }
		
		void deleteHashtable()
		{
			for(int i=0;i<usersizeofhashtable;i++)
			{
				BLASTNode* current = hashTable[i];
				BLASTNode* temp = new BLASTNode();
				
				while (current != NULL)
				{
					temp = current;
					free(temp);
					current = current->BLASTNext;
				}

				hashTable[i]= NULL;
			}

		}
   
        ~FASTAreadset_BLAST() //destructor 
            {
                deleteListBLAST(&headBLAST);
				deleteListBLAST(&dataBLAST);
				deleteHashtable();
				delete[] dataArray;
				delete[] dataArrayGenome;
			}       
};

#endif