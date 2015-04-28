//============================================================================
// Name        : ovoid3.cpp
// Author      : Hans
// Version     :
// Copyright   : Your copyright notice
// Description : This program written to find number of Ovoids(not the completed
//				 Ovoid just the amount) at field 8 and Dim 6
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
#include <iomanip>

//these defines are for qbf file
#define SIZE 8    /* size of the field */
#define P 2       /* prime characteristic */
#define DIM 6
#define VCOUNT 3530
#define N P-1



#include "qbf.h" //This is written by Athula Gunawardena. I'm using some function out of it to start the program. This is just a helper and I'm lazy to rewrite it :-)

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
void permute(vector<vector<int> > input, int i, int n,vector<vector<int> > cap1,vector<vector<int> > multi, vector<vector<int> > add,int mod, bool &failed,bool lastElement,bool normal);
bool capAreIsomorphic(vector<vector<int> > cap1,vector<vector<int> > cap2, vector<vector<int> > multi, vector<vector<int> > add,int mod,bool normal);
bool capAreEqual(vector<vector<int> > cap1,vector<vector<int> > cap2, vector<vector<int> > multi, vector<vector<int> > add,int mod);
vector<vector<vector<int> > > findCaps(vector<vector<int> > vectors,vector<vector<int> > multi, vector<vector<int> > add,int mod,vector<vector<int> > input,vector<vector<vector<int> > > oldResults);
int s_to_i(string);
vector<vector<vector<int> > > findAlltheCaps(vector<vector<int> > vectors,vector<vector<int> > multi, vector<vector<int> > add,int mod,vector<vector<int> > input);
void fileWriteCaps(vector< vector< vector <int> > > arrayVector,string fileName);
int f(vector<int>  x, vector<int>  y, vector<vector<int> > mult, vector<vector<int> > add);
static inline void loadbar(unsigned int x, unsigned int n, unsigned int w);

int main() {
	//******************Reading the file
	vector<vector<int> > multiA(SIZE, vector<int>(SIZE));
	vector<vector<int> > addA(SIZE, vector<int>(SIZE));
	vector<vector<int> > vectors(VCOUNT, vector<int>(DIM));

	readfield();
	//converting everything to vectors (me being lazy)
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
	cap3Pool[0] = vectors[0]; //picking up arbitrary first 3 vectors (some math behind it)
	cap3Pool[1] = vectors[1];
	cap3Pool[2] = vectors[2];
	vector<vector<vector<int> > > capPool = findAlltheCaps(vectors,multiA,addA,SIZE,cap3Pool);

	//if this algorithm ever finishes then its going to print out the results
	cout << "Cap pool " << capPool.size() << "*************" << endl;
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

	/*
	 * Cap size 4 work *******************************************************
	*/
	cout << "vector list size" << vectors.size() << endl;
	//Creating singular points for cap size 4 by going through cap size 3
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
	vector<vector<vector<int> > > cap5pool; //this is empty at this point coz there is no calculations have been done before. Using this to just send an empty vector list for the function

	//finding cap size 4 non isomorphic sets of vectors
	vector<vector<vector<int> > > cap4Pool = findCaps(list,multi,add,SIZE,input,cap5pool);
	cout << "\n cap4 pool size " << cap4Pool.size() << endl;

	/*
	 * Cap size 5 work *******************************************************
	*/
	//Creating singular points for cap size 5 using caps found previously
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
	cout << "cap4 list size " << cap4List[0].size() << " " << cap4List[1].size() << endl; //Debugging line

	//finding cap size 5 non isomorphic sets of vectors for each cap4 set of vector
	for(int i=0;i<cap4Pool.size();i++){
		vector<vector<vector<int> > > temp = findCaps(cap4List[i],multi,add,SIZE,cap4Pool[i],cap5pool); //doing all the calculations
		cap5pool.insert(cap5pool.end(), temp.begin(), temp.end()); //add it to the pool
		cout << "Ran through cap4 " << i << " size of " << temp.size() << " and pool size is "<< cap5pool.size() << endl;
		fileWriteCaps(cap5pool,"cap5"); //Save the pool into a file
	}
	cout << "\n cap5 size " << cap5pool.size() << endl;

	/*
	 * Cap size 6 work *******************************************************
	*/
	vector<vector<vector<int> > > cap6pool;
	//Creating singular points for cap size 6 using the previous results
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
	cout << "cap5 list size " << cap5List[0].size() << " " << cap5List[1].size() << endl; //debugging, to check how many singular point are there

	//finding cap size 6 non isomorphic sets of vectors for each cap5 vector set
	for(int i=0;i<cap5pool.size();i++){
			vector<vector<vector<int> > > temp = findCaps(cap5List[i],multi,add,SIZE,cap5pool[i],cap6pool); //calculation
			cap6pool.insert(cap6pool.end(), temp.begin(), temp.end());//adding it into the pool

			cout << "Ran through cap5 " << i << " size of " << temp.size() << " and pool size is "<< cap6pool.size() << endl;
			cout << " sample out " << endl;
			//printing (debugging step to check things)
			for(int p=0;p<cap6pool[cap6pool.size()-1].size();p++){
				for(int k=0;k<cap6pool[cap6pool.size()-1][p].size();k++){
					cout << cap6pool[cap6pool.size()-1][p][k] << " ";
				}
				cout << endl;
			}
			fileWriteCaps(cap6pool,"cap6"); //Writing to the file
		}
		cout << "\n cap6 size " << cap6pool.size() << endl; //printing the size of pool

	return cap6pool; //finally done
}

//This where all the magic happens
vector<vector<vector<int> > > findCaps(vector<vector<int> > vectors,vector<vector<int> > multi, vector<vector<int> > add,int mod,vector<vector<int> > input,vector<vector<vector<int> > > oldResults){

	vector< vector<int> > capResults;
	vector< vector<vector<int> > > allCaps;

	for(int i=0;i<vectors.size();i++){ //going through individual singular points
		loadbar(i,vectors.size(),50); //printing the status to the terminal
		//looks like i dont need this part
//		bool notBil = false;
//		for(int j=0;j<input.size();j++){
//			if(bilinearFormAny(input[j],vectors[i],multi,add,vectors[i].size(), mod)==0){
//				notBil = true;
//				break;
//			}
//		}
//		if(!notBil){
			capResults = input; //grabbing the vector set for particular cap. This is from the old results from previous cap
			capResults.push_back(vectors[i]);
			bool isoFailed = false;
			#pragma omp parallel for //for compiler directive for for OpenMp
			for(int j=0;j<allCaps.size();j++){
				#pragma omp flush (isoFailed)
				if(!isoFailed && capAreIsomorphic(allCaps[j],capResults,multi,add,SIZE,false)){ //check for isomorphism. AKA, check whether vector sets are the same
					isoFailed = true;
					#pragma omp flush (isoFailed)
				}
			}
			if(!isoFailed){ //this where if the vector pass the isomorphism with its own results, check agents the pool of results that been calculated previously.
				#pragma omp parallel for //for compiler directive for for OpenMp
				for(int p=0;p<oldResults.size();p++){
					#pragma omp flush (isoFailed)
					if(!isoFailed && capAreIsomorphic(oldResults[p],capResults,multi,add,SIZE,false)){ //old results here coming from the pool (this takes long time) I think i can do performance increase here by keep the old results in memory
						isoFailed = true;
						#pragma omp flush (isoFailed)
					}
				}
			}
			if(!isoFailed){ //if all the conditions pass up there, then this is worthy of going into the selected results for particular set
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
bool capAreIsomorphic(vector<vector<int> > cap1,vector<vector<int> > cap2, vector<vector<int> > multi, vector<vector<int> > add,int mod,bool normal){
	//This function permute all the vectors and check whether they are isomorphic
	bool failed = false;
	permute(cap2,0,cap2.size()-1,cap1,multi,add,mod,failed,false,normal); //capSize-1 = 3 and starts at 0
//	cout << "condition is " << failed << endl;
	return failed;
}


void permute(vector<vector<int> > input, int i, int n,vector<vector<int> > cap1,vector<vector<int> > multi, vector<vector<int> > add,int mod, bool &failed,bool lastElement,bool normal)
{
   int j;
   if (i == n){
	   failed = capAreEqual(cap1,input,multi,add,mod); //checking whether they are equal, AKA isomorphic or not
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
        	//there are lots of math optimizations going on here but they failed and gave bad results.
//        if(failed)
//        	break;
//        else if(normal){
//        	input[i].swap(input[j]);
//        	permute(input, i+1, n,cap1,multi,add,mod,failed,lastElement,normal);
//        	input[i].swap(input[j]); //backtrack
//        }
//        else if((lastElement || j==(n-1))){
//        	if(j==(n-1)) 	lastElement = true;
//			input[i].swap(input[j]);
//			permute(input, i+1, n,cap1,multi,add,mod,failed,lastElement,normal);
//			input[i].swap(input[j]); //backtrack
//        	}
        	//brute force method. math optimizations didnt add up
        	if(failed)
        	     break;
        	input[i].swap(input[j]);
			permute(input, i+1, n,cap1,multi,add,mod,failed,lastElement,normal);
			input[i].swap(input[j]); //backtrack

       }
   }
}

bool capAreEqual(vector<vector<int> > cap1,vector<vector<int> > cap2, vector<vector<int> > multi, vector<vector<int> > add,int mod){
	//get the grammian out of the vectors and check agents each other. AKA isomorphic or not
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

	//Grammian matrix (lots of math)
	vector< vector<int> > array = grammianMatrix(input,multi,add,mod);
	for(int i=1; i<input.size();i++){
			int multiple = findMultiple(array[0][i],multi);
			input[i]=multiplyVector(multiple,input[i],multi,mod);
			array = grammianMatrix(input,multi,add,mod);
	}

//	cout << "checkCap " << endl;
	//I only need the half of the matrix
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

	if(returnArray[0]!=1){ //if it is 1 then there is nothing to do
		int multiplicationNumber=findMultiple(returnArray[0],multi); //dont have to multiply manually, just have to find multiple from multiplication array
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

	//generate grammianMatrix, this is just raw form of it
	vector<vector<int> > results(matrix.size(),vector<int>(matrix.size(),0));
	#pragma omp parallel for
	for(int i=0; i<matrix.size();i++){
		for(int j=i+1; j<matrix.size();j++){
			results[i][j]=results[j][i]=bilinearFormAny(matrix[i],matrix[j],multi,add,matrix[i].size(), mod); // take the bilinear to generate the grammain
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
static inline void loadbar(unsigned int x, unsigned int n, unsigned int w = 50)
{
    if ( (x != n) && (x % (n/100+1) != 0) ) return;

    float ratio  =  x/(float)n;
    int   c      =  ratio * w;

    cout << setw(3) << (int)(ratio*100) << "% [";
    for (int x=0; x<c; x++) cout << "=";
    for (int x=c; x<w; x++) cout << " ";
    cout << "]\r" << flush;
}
