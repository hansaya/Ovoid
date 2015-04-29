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

//void findPoints(bool, int, int &, int, int, int[], int, vector<vector<int> > &);
vector<int> modCheck(vector<int> , int);//passing in vector of ints, number we will be modding with
//void calculation(vector< vector<int> >,vector< vector< vector <int> > >,int,int &,vector< vector< vector<int> > > &,int &,int,int,int, int *);
int bilinearFormOdd(vector<int> ,vector<int>  ,vector<vector<int> > multi, vector<vector<int> > add,int , int);
void bilinearFormAny(const vector<int> &arrayA,const vector<int> &arrayB ,const vector<vector<int> > &multi,const vector<vector<int> > &add,const int &size,const int &mod,int &results);
int quadraticOdd(vector<int> array,vector<vector<int> > multi,vector<vector<int> > add);
int quadraticAny(const vector<int> &array,const vector<vector<int> > &multi,const vector<vector<int> > &add);
void generateGrammian(vector<vector<int> > input,const vector<vector<int> > &multi,const vector<vector<int> > &add,const int &mod,vector<int> &returnArray);
void grammianMatrix(const vector<vector<int> > &matrix,const vector<vector<int> > &multi,const vector<vector<int> > &add,const int &mod,vector<vector<int> > &results);
void findMultiple(const int &number,const vector<vector<int> > &multi,int &results);
void multiplyVector(const int &numMult,const vector<int> &input,const vector<vector<int> > &multi,const int &mod,vector<int> &results);
void permute(vector<vector<int> > input, int i,const int &n,const vector<vector<int> > &cap1,const vector<vector<int> > &multi, const vector<vector<int> > &add,const int &mod, bool &failed,bool lastElement,bool normal);
bool capAreIsomorphic(const vector<vector<int> > &cap1, const vector<vector<int> > &cap2, const vector<vector<int> > &multi, const vector<vector<int> > &add, const int &mod, const bool &normal);
bool capAreEqual(const vector<vector<int> > &cap1,const vector<vector<int> > &cap2,const vector<vector<int> > &multi,const vector<vector<int> > &add,const int &mod);
void findCaps(const vector<vector<int> > &vectors,const vector<vector<int> > &multi, const vector<vector<int> > &add,const int &mod,const vector<vector<int> > &input,const vector<vector<vector<int> > > &oldResults, vector< vector<vector<int> > > &allCapsResults);
int s_to_i(string);
vector<vector<vector<int> > > findAlltheCaps(const vector<vector<int> > &vectors,const vector<vector<int> > &multi,const vector<vector<int> > &add,const int &mod,const vector<vector<int> > &input);
void fileWriteCaps(const vector< vector< vector <int> > > &arrayVector,const string &fileName,const int &i);
int fileReadCaps(vector< vector< vector <int> > > &result,const string &fileName);
bool file_empty(const string &fileName);
int f(const vector<int>  &x,const vector<int>  &y,const vector<vector<int> > &mult,const vector<vector<int> > &add);
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
	cout << "Cap pool " << capPool.size() << "************* " << capPool[1].size() << " --- " << capPool[1][0].size() <<endl;
	for(int i=0;i<capPool.size();i++){
		for(int j=0;j<capPool[i].size();j++){
			for(int k=0;k<capPool[i][j].size();k++)
				cout << capPool[i][j][k];
			cout << endl;
		}
		cout << "cap "<< i <<" ends " << endl;
	}

	return 0;
}
vector<vector<vector<int> > > findAlltheCaps(const vector<vector<int> > &vectors,const vector<vector<int> > &multi,const vector<vector<int> > &add,const int &mod,const vector<vector<int> > &input){

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
	vector<vector<vector<int> > > cap4Pool;
	findCaps(list,multi,add,SIZE,input,cap5pool,cap4Pool);
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

	//Resume the program if a file exists
	int startIndex = 0;
	if(!file_empty("cap5")){
		startIndex = fileReadCaps(cap5pool,"cap5")+1;
		cout << "resuming at pool size: " << cap5pool.size() << " and starting at : " << startIndex << endl;
	}

	//finding cap size 5 non isomorphic sets of vectors for each cap4 set of vector
	for(int i=startIndex;i<cap4Pool.size();i++){
		vector<vector<vector<int> > > temp;
		findCaps(cap4List[i],multi,add,SIZE,cap4Pool[i],cap5pool,temp); //doing all the calculations
		cap5pool.insert(cap5pool.end(), temp.begin(), temp.end()); //add it to the pool
		cout << endl << "Ran through cap4 " << i << " size of " << temp.size() << " and pool size is "<< cap5pool.size() << endl;
		fileWriteCaps(cap5pool,"cap5",i); //Save the pool into a file
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

	//Resume the program if a file exists
	startIndex = 0;
	if(!file_empty("cap6")){
		startIndex = fileReadCaps(cap6pool,"cap6")+1;
		cout << "resuming at pool size: " << cap6pool.size() << " and starting at : " << startIndex << endl;
	}

	//finding cap size 6 non isomorphic sets of vectors for each cap5 vector set
	for(int i=startIndex;i<cap5pool.size();i++){
			vector<vector<vector<int> > > temp;
			findCaps(cap5List[i],multi,add,SIZE,cap5pool[i],cap6pool,temp); //calculation
			cap6pool.insert(cap6pool.end(), temp.begin(), temp.end());//adding it into the pool

			cout << endl << "Ran through cap5 " << i << " size of " << temp.size() << " and pool size is "<< cap6pool.size() << endl;
			cout << " sample out " << endl;
			//printing (debugging step to check things)
			for(int p=0;p<cap6pool[cap6pool.size()-1].size();p++){
				for(int k=0;k<cap6pool[cap6pool.size()-1][p].size();k++){
					cout << cap6pool[cap6pool.size()-1][p][k] << " ";
				}
				cout << endl;
			}
			fileWriteCaps(cap6pool,"cap6",i); //Writing to the file
		}
		cout << "\n cap6 size " << cap6pool.size() << endl; //printing the size of pool

	return cap6pool; //finally done
}

//This where all the magic happens
void findCaps(const vector<vector<int> > &vectors,const vector<vector<int> > &multi, const vector<vector<int> > &add,const int &mod,const vector<vector<int> > &input,const vector<vector<vector<int> > > &oldResults, vector< vector<vector<int> > > &allCapsResults){

	vector< vector<int> > capResults;
//	vector< vector<vector<int> > > allCaps;

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
			for(int j=0;j<allCapsResults.size();j++){
				#pragma omp flush (isoFailed)
				if(!isoFailed && capAreIsomorphic(allCapsResults[j],capResults,multi,add,SIZE,false)){ //check for isomorphism. AKA, check whether vector sets are the same
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
				allCapsResults.push_back(capResults);
			}
			capResults.clear();
		}
//	}
//	return allCapsResults;
}
vector<int> modCheck(vector<int> abc, int modulo) {
	for (int i = 0; i < abc.size(); i++) {
		abc[i] = abc[i] % modulo;
	}
	return abc;

}//end modCheck function
bool capAreIsomorphic(const vector<vector<int> > &cap1, const vector<vector<int> > &cap2, const vector<vector<int> > &multi, const vector<vector<int> > &add, const int &mod, const bool &normal){
	//This function permute all the vectors and check whether they are isomorphic
	bool failed = false;
	permute(cap2,0,cap2.size()-1,cap1,multi,add,mod,failed,false,normal); //capSize-1 = 3 and starts at 0
	return failed;
}


//This is a recursive method
void permute(vector<vector<int> > input, int i,const int &n,const vector<vector<int> > &cap1,const vector<vector<int> > &multi, const vector<vector<int> > &add,const int &mod, bool &failed,bool lastElement,bool normal)
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

bool capAreEqual(const vector<vector<int> > &cap1,const vector<vector<int> > &cap2,const vector<vector<int> > &multi,const vector<vector<int> > &add,const int &mod){
	//get the grammian out of the vectors and check agents each other. AKA isomorphic or not
	vector<int> cap1Results;
	generateGrammian(cap1,multi,add,mod,cap1Results);
	vector<int> cap2Results;
	generateGrammian(cap2,multi,add,mod,cap2Results);

	if(cap1Results == cap2Results){
		return true;
	}
	return false;

}
void generateGrammian(vector<vector<int> > input,const vector<vector<int> > &multi,const vector<vector<int> > &add,const int &mod,vector<int> &returnArray){

	//Grammian matrix (lots of math)
	vector< vector<int> > array(input.size(),vector<int>(input.size(),0));
	grammianMatrix(input,multi,add,mod,array);
	for(int i=1; i<input.size();i++){
		if(array[0][i]!=1){
			int multiple;
			findMultiple(array[0][i],multi,multiple);
			multiplyVector(multiple,input[i],multi,mod,input[i]);
//			array = grammianMatrix(input,multi,add,mod);
		}

	}
	grammianMatrix(input,multi,add,mod,array);

//	#pragma omp parallel for
	for(int i=1;i<array.size()-1;i++){
		for(int j=i+1;j<array[i].size();j++){
			returnArray.push_back(array[i][j]);
		}
	}



	if(returnArray[0]!=1){ //if it is 1 then there is nothing to do
		int multiplicationNumber;
		findMultiple(returnArray[0],multi,multiplicationNumber); //dont have to multiply manually, just have to find multiple from multiplication array
//		#pragma omp parallel for
		for(int i=0;i<returnArray.size();i++){
			returnArray[i]=multi[returnArray[i]][multiplicationNumber];
		}
	}

}
void bilinearFormAny(const vector<int> &arrayA,const vector<int> &arrayB ,const vector<vector<int> > &multi,const vector<vector<int> > &add,const int &size,const int &mod,int &results){
	results = 0;
	for (int j = 0; j < size; j++) {
		results = add[multi[arrayA[j]][arrayB[(size-1) - j]]][results];
	}
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
int quadraticAny(const vector<int> &array,const vector<vector<int> > &multi,const vector<vector<int> > &add){
	int sum = 0;
	for (int j = 0; j < array.size()/2; j++) {
		sum = add[multi[array[j]][array[(array.size()-1)-j]]][sum];
	}
	return sum;
}
void grammianMatrix(const vector<vector<int> > &matrix,const vector<vector<int> > &multi,const vector<vector<int> > &add,const int &mod,vector<vector<int> > &results){

	//generate grammianMatrix, this is just raw form of it
//	#pragma omp parallel for
	for(int i=0; i<matrix.size();i++){
		for(int j=i+1; j<matrix.size();j++){
			bilinearFormAny(matrix[i],matrix[j],multi,add,matrix[i].size(), mod,results[i][j]); // take the bilinear to generate the grammian =results[j][i]
		}
	}
}
void findMultiple(const int &number,const vector<vector<int> > &multi,int &results){
	for(int i=0;i<multi.size();i++){
		if(multi[number][i]==1){
			results= i;
			break;
		}
	}
}
void multiplyVector(const int &numMult,const vector<int> &input,const vector<vector<int> > &multi,const int &mod,vector<int> &results){

//		#pragma omp parallel for
		for(int i=0;i<input.size();i++){
			results[i]=multi[input[i]][numMult];
		}
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


int f(const vector<int>  &x,const vector<int>  &y,const vector<vector<int> > &mult,const vector<vector<int> > &add)
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
void fileWriteCaps(const vector< vector< vector <int> > > &arrayVector,const string &fileName,const int &o){

	ofstream outFile(fileName.c_str(),ios::out|ios::binary);
	if(outFile){
		for(int z=0;z<arrayVector.size();z++){
			for(int i=0;i<arrayVector[z].size();i++){
				for(int k=0;k<arrayVector[z][i].size();k++){
					outFile << arrayVector[z][i][k]<< " " ;
				}
				outFile << endl;
			}
			outFile << "|" << endl;
		}
		outFile << "***" << endl << o << endl;
		outFile.close();
	}
	else cerr << "Something really went wrong!!" << endl;

}
int fileReadCaps(vector< vector< vector <int> > > &result,const string &fileName){

	ifstream redFile(fileName.c_str(),ios::in|ios::binary);
		if(redFile){
			int count =0;
			//Going though each line
			for(string ovoid,value;getline(redFile, ovoid,'|');){
				//Dividing the ovoid
				istringstream ovd(ovoid);
				vector< vector<int> > data;\
				bool hitTheEnd = false;
				for(string lines,value;getline(ovd, lines);){
					if(lines=="***"|| hitTheEnd){
						if(hitTheEnd)
							return s_to_i(lines);
						hitTheEnd = true;
					}
					else if(lines!=""){
						//Dividing the line
						istringstream line(lines);
						vector<int> array;
						while(getline(line,value,' ')){
							array.push_back(s_to_i(value));
						}
						data.push_back(array);
					}
				}
				result.push_back(data);
				count++;
			}
		}
		else cerr << "Something really went wrong!!" << endl;
		redFile.close();
		return -1;

}
//Checking whether file is there
bool file_empty(const string &fileName)
{
	ifstream pFile(fileName.c_str(),ios::in|ios::binary);
    return pFile.peek() == ifstream::traits_type::eof();
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
