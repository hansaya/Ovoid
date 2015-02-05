//============================================================================
// Name        : ovoid3.cpp
// Author      : Hans
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <stdio.h>
#include <omp.h>
#include <sys/types.h>
#include <fstream>
#include <string>

#define SIZE 8    /* size of the field */
//#define FIELD 8
#define P 2       /* prime characteristic */
#define DIM 6
#define VCOUNT 3530
#define N P-1
//#define MAX 4097


//int Q(), f(), mult(), add();

#include "qbf.h"

using namespace std;

void findPoints(bool, int, int &, int, int, int[], int, vector<vector<int> > &);
vector<int> modCheck(vector<int> , int);//passing in vector of ints, number we will be modding with
void calculation(vector< vector<int> >,vector< vector< vector <int> > >,int,int &,vector< vector< vector<int> > > &,int &,int,int,int, int *);
int bilinearFormOdd(vector<int> ,vector<int>  ,vector<vector<int> > multi, vector<vector<int> > add,int , int);
int bilinearFormAny(vector<int> ,vector<int> ,vector<vector<int> > multi, vector<vector<int> > add,int size, int );
int quadraticOdd(vector<int> array,vector<vector<int> > multi,vector<vector<int> > add);
int quadraticAny(vector<int> array,vector<vector<int> > multi,vector<vector<int> > add);
vector<int> generateGrammian(vector<vector<int> > input,vector<vector<int> > multi, vector<vector<int> > add,int mod);
vector<vector<int> > grammianMatrix(vector<vector<int> > ,vector<vector<int> > , vector<vector<int> >, int );
int findMultiple(int number,vector<vector<int> > multi);
vector<int> multiplyVector(int numMult,vector<int> input,vector<vector<int> > multi ,int mod);
void permute(vector<vector<int> > input, int i, int n,vector<vector<int> > cap1,vector<vector<int> > multi, vector<vector<int> > add,int mod, bool &failed,bool lastElement);
bool capAreIsomorphic(vector<vector<int> > cap1,vector<vector<int> > cap2, vector<vector<int> > multi, vector<vector<int> > add,int mod);
bool capAreEqual(vector<vector<int> > cap1,vector<vector<int> > cap2, vector<vector<int> > multi, vector<vector<int> > add,int mod);
vector<vector<vector<int> > > findCaps(vector<vector<int> > vectors,vector<vector<int> > multi, vector<vector<int> > add,int mod,vector<vector<int> > input,vector<vector<vector<int> > > oldResults);
int s_to_i(string);
vector<vector<vector<int> > > findAlltheCaps(vector<vector<int> > vectors,vector<vector<int> > multi, vector<vector<int> > add,int mod,vector<vector<int> > input);
void fileWriteCaps(vector< vector< vector <int> > > arrayVector,string fileName);
int f(vector<int>  x, vector<int>  y, vector<vector<int> > mult, vector<vector<int> > add);

int main() {
	//******************Reading the file
	vector<vector<int> > multiA(SIZE, vector<int>(SIZE));
	vector<vector<int> > addA(SIZE, vector<int>(SIZE));
	vector<vector<int> > vectors(VCOUNT, vector<int>(DIM));

	readfield();

	for(int i=0;i<SIZE;i++){
		for(int j=0;j<SIZE;j++){
			multiA[i][j]=mult(i,j);
			addA[i][j]=add(i,j);
			cout << multiA[i][j] << " ";
		}
		cout << endl;
	}

	readVectors();
	for(int i=0;i<VCOUNT;i++){
		for(int j=0;j<DIM;j++){
			vectors[i][j]=vect(i,j);
			//cout << vectors[i][j];
		}
		//cout << endl;
	}

	//Picking caps
	vector<vector<vector<int> > > capResults;
	vector<vector<int> > cap3Pool(3, vector<int>(DIM));
	cap3Pool[0] = vectors[0];
	cap3Pool[1] = vectors[1];
	cap3Pool[2] = vectors[2];
	vector<vector<vector<int> > > capPool = findAlltheCaps(vectors,multiA,addA,SIZE,cap3Pool);

	cout << "Cap4 pool " << capPool.size() << "*************" << endl;
	for(int i=0;i<capPool.size();i++){
		for(int j=0;j<capPool[i].size();j++){
			for(int k=0;k<capPool[i][j].size();k++)
				cout << capPool[i][j][k];
			cout << endl;
		}
		cout << "cap"<< i <<" ends " << endl;
	}

	return 0;
}
vector<vector<vector<int> > > findAlltheCaps(vector<vector<int> > vectors,vector<vector<int> > multi, vector<vector<int> > add,int mod,vector<vector<int> > input){

	//Cap size 4 work *******************************************************
	cout << "vector size" << vectors.size() << endl;
	//Creating singular points for cap size 4
	vector<vector<int> > list;
	for(int j=0;j<vectors.size();j++){
		bool notSingular = false;
		for(int k=0;k<input.size();k++){
			if(f(input[k],vectors[j],multi,add)==0){
				notSingular =true;
				break;
			}
		}
		if(!notSingular)
			list.push_back(vectors[j]);
	}
	cout << "cap3 list size " << list.size() << endl;
	vector<vector<vector<int> > > cap5pool;

	//finding cap size 4 non isomorphic sets of vectors
	vector<vector<vector<int> > > cap4Pool = findCaps(list,multi,add,SIZE,input,cap5pool);
	cout << "cap4 pool size " << cap4Pool.size() << endl;

	//Cap size 5 work *******************************************************
	//Creating singular points for cap size 4
	vector<vector<vector <int> > > cap4List;
	for(int i=0;i<cap4Pool.size();i++){
		vector<vector<int> > list;
		for(int j=0;j<vectors.size();j++){
			bool notSingular = false;
			for(int k=0;k<cap4Pool[i].size();k++){
				if(f(cap4Pool[i][k],vectors[j],multi,add)==0){
					notSingular =true;
					break;
				}
			}
			if(!notSingular)
				list.push_back(vectors[j]);
		}
		cap4List.push_back(list);
	}
	cout << "cap4 list size " << cap4List[0].size() << " " << cap4List[1].size() << endl;

	//finding cap size 5 non isomorphic sets of vectors
	for(int i=0;i<cap4Pool.size();i++){
		vector<vector<vector<int> > > temp = findCaps(cap4List[i],multi,add,SIZE,cap4Pool[i],cap5pool);
		cap5pool.insert(cap5pool.end(), temp.begin(), temp.end());
		cout << "Ran through cap4 " << i << " size of " << temp.size() << " and pool size is "<< cap5pool.size() << endl;
		fileWriteCaps(cap5pool,"cap5");
	}
	cout << "cap5 size " << cap5pool.size() << endl;

	//Cap size 6 work *******************************************************
	vector<vector<vector<int> > > cap6pool;
	//Creating singular points for cap size 5
	vector<vector<vector <int> > > cap5List;
	for(int i=0;i<cap5pool.size();i++){
		vector<vector<int> > list;
		for(int j=0;j<vectors.size();j++){
			bool notSingular = false;
			for(int k=0;k<cap5pool[i].size();k++){
				if(f(cap5pool[i][k],vectors[j],multi,add)==0){
					notSingular =true;
					break;
				}
			}
			if(!notSingular)
				list.push_back(vectors[j]);
		}
		cap5List.push_back(list);
	}
	cout << "cap4 list size " << cap5List[0].size() << " " << cap5List[1].size() << endl;

	//finding cap size 6 non isomorphic sets of vectors
	for(int i=0;i<cap5pool.size();i++){
			vector<vector<vector<int> > > temp = findCaps(cap5List[i],multi,add,SIZE,cap5pool[i],cap6pool);
			cap6pool.insert(cap6pool.end(), temp.begin(), temp.end());

			cout << "Ran through cap5 " << i << " size of " << temp.size() << " and pool size is "<< cap6pool.size() << endl;
			cout << " Example run " << endl;
			//printing
			for(int p=0;p<cap6pool[cap6pool.size()-1].size();p++){
				for(int k=0;k<cap6pool[cap6pool.size()-1][p].size();k++){
					cout << cap6pool[cap6pool.size()-1][p][k] << " ";
				}
				cout << endl;
			}
			fileWriteCaps(cap6pool,"cap6");
		}
		cout << "cap6 size " << cap6pool.size() << endl;

	return cap6pool;
}

//recursion part for the finding ovoid
vector<vector<vector<int> > > findCaps(vector<vector<int> > vectors,vector<vector<int> > multi, vector<vector<int> > add,int mod,vector<vector<int> > input,vector<vector<vector<int> > > oldResults){

	vector< vector<int> > capResults;
	vector< vector<vector<int> > > allCaps;

	for(int i=0;i<vectors.size();i++){
//		bool notBil = false;
//		for(int j=0;j<input.size();j++){
//			if(bilinearFormAny(input[j],vectors[i],multi,add,vectors[i].size(), mod)==0){
//				notBil = true;
//				break;
//			}
//		}
//		if(!notBil){
			capResults = input;
			capResults.push_back(vectors[i]);
			bool isoFailed = false;
			#pragma omp parallel for
			for(int j=0;j<allCaps.size();j++){
				#pragma omp flush (isoFailed)
				if(!isoFailed && capAreIsomorphic(allCaps[j],capResults,multi,add,SIZE)){
					isoFailed = true;
					#pragma omp flush (isoFailed)
				}
			}
			if(!isoFailed){
				#pragma omp parallel for
				for(int p=0;p<oldResults.size();p++){
					#pragma omp flush (isoFailed)
					if(!isoFailed && capAreIsomorphic(oldResults[p],capResults,multi,add,SIZE)){
						isoFailed = true;
						#pragma omp flush (isoFailed)
					}
				}
			}
			if(!isoFailed){
//				cout << "found one " << endl;
				allCaps.push_back(capResults);
			}
			capResults.clear();
		}
//	}
	return allCaps;
}
vector<int> modCheck(vector<int> abc, int modulo) {
	for (int i = 0; i < abc.size(); i++) {
		abc[i] = abc[i] % modulo;
	}
	return abc;

}//end modCheck function
bool capAreIsomorphic(vector<vector<int> > cap1,vector<vector<int> > cap2, vector<vector<int> > multi, vector<vector<int> > add,int mod){
	//for()
	bool failed = false;
	permute(cap2,0,cap2.size()-1,cap1,multi,add,mod,failed,false); //capSize-1 = 3 and starts at 0
//	cout << "condition is " << failed << endl;
	return failed;
}


void permute(vector<vector<int> > input, int i, int n,vector<vector<int> > cap1,vector<vector<int> > multi, vector<vector<int> > add,int mod, bool &failed,bool lastElement)
{
   int j;
   if (i == n){
	   failed = capAreEqual(cap1,input,multi,add,mod);
//	   for(int k=0;k<input.size();k++){
//		   for(int p=0;p<input[k].size();p++){
//			   cout << input[k][p] << " ";
//		   }
//		   cout << endl;
//	   }
//	   cout <<"******************* "<< failed << endl;
   }
   else if(!failed)
   {
        for (j = i; j <= n; j++)
       {
        if(failed)
        	break;
        else if(lastElement || j==(n-1)){
        	if(j==(n-1)) 	lastElement = true;
			input[i].swap(input[j]);
			permute(input, i+1, n,cap1,multi,add,mod,failed,lastElement);
			input[i].swap(input[j]); //backtrack
        }
       }
   }
}

bool capAreEqual(vector<vector<int> > cap1,vector<vector<int> > cap2, vector<vector<int> > multi, vector<vector<int> > add,int mod){
	vector<int> cap1Results=generateGrammian(cap1,multi,add,mod);
	vector<int> cap2Results=generateGrammian(cap2,multi,add,mod);
//	cout << "cap1 abc : ";
//	for(int i=0;i<cap1Results.size();i++)
//		 cout << cap1Results[i] << " ";
//	cout << endl << "cap2 abc : ";
//	for(int i=0;i<cap1Results.size();i++)
//		cout << cap2Results[i] << " ";
//	cout << endl;
//	cout << "###################################################################### "<<endl;
	if(cap1Results == cap2Results){
		return true;
	}
	return false;

}
vector<int> generateGrammian(vector<vector<int> > input,vector<vector<int> > multi, vector<vector<int> > add,int mod){

	vector< vector<int> > array = grammianMatrix(input,multi,add,mod);
	for(int i=1; i<input.size();i++){
			int multiple = findMultiple(array[0][i],multi);
			input[i]=multiplyVector(multiple,input[i],multi,mod);
			array = grammianMatrix(input,multi,add,mod);
	}

//	cout << "checkCap " << endl;
	vector<int> returnArray;
	#pragma omp parallel for
	for(int i=1;i<array.size()-1;i++){
		for(int j=i+1;j<array[i].size();j++){
			returnArray.push_back(array[i][j]);
//			cout << " &&&&&&&&&&&& " << i << " " << j << " ";
		}
	}
//	cout << endl;
//			cout << array[i][j] << " ";
//			if(i!=0 && j!=0 && i!=j)
//				returnArray.push_back(array[i][j]);
//		}
//		cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<< endl;
//	}

	if(returnArray[0]!=1){
		int multiplicationNumber=findMultiple(returnArray[0],multi);
		#pragma omp parallel for
		for(int i=0;i<returnArray.size();i++){
			returnArray[i]=multi[returnArray[i]][multiplicationNumber];
		}
	}

	return returnArray;
}
int bilinearFormAny(vector<int> arrayA,vector<int> arrayB ,vector<vector<int> > multi,vector<vector<int> > add,int size, int mod){
	int sum = 0;
	for (int j = 0; j < size; j++) {
		sum = add[multi[arrayA[j]][arrayB[(size-1) - j]]][sum];
	}
	return sum;
}
int bilinearFormOdd(vector<int> arrayA,vector<int> arrayB ,vector<vector<int> > multi,vector<vector<int> > add,int size, int mod){
	int sum = 0;
	for (int j = 0; j < size; j++) {
		sum = add[multi[arrayA[j]][arrayB[j]]][sum];
	}
	return sum;
}
int quadraticOdd(vector<int> array,vector<vector<int> > multi,vector<vector<int> > add){
	int sum = 0;
	for (int j = 0; j < array.size(); j++) {
		sum = add[multi[array[j]][array[j]]][sum];
	}
	return sum;
}
int quadraticAny(vector<int> array,vector<vector<int> > multi,vector<vector<int> > add){
	int sum = 0;
	for (int j = 0; j < array.size()/2; j++) {
		sum = add[multi[array[j]][array[(array.size()-1)-j]]][sum];
	}
	return sum;
}
vector<vector<int> > grammianMatrix(vector<vector<int> > matrix,vector<vector<int> > multi, vector<vector<int> > add, int mod){

	vector<vector<int> > results(matrix.size(),vector<int>(matrix.size(),0));
	#pragma omp parallel for
	for(int i=0; i<matrix.size();i++){
		for(int j=i+1; j<matrix.size();j++){
			results[i][j]=results[j][i]=bilinearFormAny(matrix[i],matrix[j],multi,add,matrix[i].size(), mod);
		}
	}
	return results;
}
int findMultiple(int number,vector<vector<int> > multi){
	for(int i=0;i<multi.size();i++){
		if(multi[number][i]==1){
			return i;
		}
	}
	return -1;
}
vector<int> multiplyVector(int numMult,vector<int> input, vector<vector<int> > multi, int mod){

		#pragma omp parallel for
		for(int i=0;i<input.size();i++){
			input[i]=multi[input[i]][numMult];
		}
		return input;
}
//converting any string to int
int s_to_i(string s )
{
  istringstream i(s);
  int x;
  if (!(i >> x))
    return 0;
  return x;
}


int f(vector<int>  x, vector<int>  y, vector<vector<int> > mult, vector<vector<int> > add)
  {
    vector<int> z(x.size());
    int sum=0;

   for(int i=0; i < x.size(); i++)
       z[i] = add[x[i]][y[i]];

    sum = quadraticAny(z,mult,add);
    sum = add[sum][mult[N][quadraticAny(x,mult,add)]];
    sum = add[sum][mult[N][quadraticAny(y,mult,add)]];

    return(sum);
  }
void fileWriteCaps(vector< vector< vector <int> > > arrayVector,string fileName){

	ofstream outFile(fileName.c_str(),ios::out|ios::binary);
	if(outFile){
		for(int z=0;z<arrayVector.size();z++){
			for(int i=0;i<arrayVector[z].size();i++){
				for(int k=0;k<arrayVector[i].size();k++){
					outFile << arrayVector[z][i][k]<< " " ;
				}
				outFile << endl;
			}
			outFile << "|" << endl;
		}
		outFile.close();
	}
	else cerr << "Something really went wrong!!" << endl;

}
