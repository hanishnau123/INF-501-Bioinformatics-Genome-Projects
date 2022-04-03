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
#include<vector>

using namespace std;

class Node{
	char charValue;
	int suffixRow;
	public:
	vector<Node*> childNode;
	
	Node(){

	}

	Node(char v , int row){
		charValue=v;
		suffixRow = row;
	}

	void addChildNode(Node* x){
		childNode.push_back(x);
	}

	Node* searchChildNode(char start){
		for(int i = 0; i < childNode.size(); i++)
		{
			if(childNode.at(i) -> charValue==start)
			{
				return childNode.at(i);
			}
		}
		return NULL;
	}

	int getRow() {
		return suffixRow;
	}

	char getValue() {
		return charValue;
	}
};

class suffix_trie {

	int suffixSize;	
	char* charArray;
	char** Suffixes;
	Node* root;

	void insertIntoSuffix(){
		for (int i = 0 ; i < suffixSize; i++) {
			Node* temp = root;
			for (int j = 0; j < suffixSize; j++) {
				if (Suffixes[i][j - 1] == '$')
				{
					break;
				} 
				if (temp->searchChildNode(Suffixes[i][j]) != NULL) 
				{
					temp = temp->searchChildNode(Suffixes[i][j]);
					continue;
				}
				else 
				{
					Node* newnode=new Node(Suffixes[i][j],i);
					temp->addChildNode(newnode);
					temp = newnode;
				}
			}
		}
	}

	void buildSuffix() {
		string suffix;
		int rowSuffix = 0;
		for(int j=0; j < suffixSize; j++ , rowSuffix++ )
		{
			suffix = "";
			for (int i = j; i < suffixSize; i++) 
			{
				suffix += charArray[i];
			}
			for(int k=0, colSuffix=0;k<suffix.size();k++, colSuffix++)
			{
				Suffixes[rowSuffix][colSuffix]=suffix[k];
			}
		}
	}

	void Traverse(Node* temp, int length,int StringSize) {
		if (temp->childNode.size()==0) 
		{
			cout << "The sequence is present in the suffix trie" << endl;
			return;
		}
		for (int i = 0; i < temp->childNode.size(); i++) 
		{
			Traverse(temp->childNode.at(i), length + 1, StringSize);
		}
	}
	
public:
	suffix_trie(string a) 
	{
		suffixSize = a.size();
		charArray =new char[suffixSize];
		for (int i = 0; i <suffixSize; i++) 
		{
			charArray[i] = a[i];
		}
		Suffixes=new char*[suffixSize];
		for(int i=0;i<suffixSize;i++){
			Suffixes[i]=new char[suffixSize];
		}
			
		root = new Node();
		buildSuffix();
		insertIntoSuffix();	
	}
	
	void Search(string s)
	{
		Node* temp = root;
		for (int i = 0; i < s.size(); i++) 
		{
			temp = temp->searchChildNode(s[i]);
			if(temp== NULL)
			{
				cout << "This suffix isn't in the tree" << endl;
				return;
			}
		}

		Node* traverse = temp;
		Traverse(traverse, s.size(), suffixSize);
		cout << endl;
	}

	int printSize(){
		return suffixSize;
	}
};

string convertToString(char* a, int size)
{
    int i;
    string convertedString = "";
    for (i = 0; i < size; i++) 
	{
        convertedString = convertedString + a[i];
    }

    return convertedString;
}

string* randomerGenerator(int n)
{
    char alphabet[4] = { 'A', 'C', 'T', 'G'};
  
    string res = "";
	string *randomerArray = new string[n];

	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < 36; i++)
		{
			res = res + alphabet[rand() % 4];
		}

		randomerArray[j] = res;
		res = "";
	}

    return randomerArray;
}

string readGenome(string Path)
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
		
		string convertedString;
		convertedString = convertToString(copyChar, index);
		return convertedString;
	}
}

#endif
