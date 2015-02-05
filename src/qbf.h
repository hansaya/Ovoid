/* This file contains the subroutines for reading the field gf(9) 
   quadratic form and bilinear form  */


int a[SIZE][SIZE], m[SIZE][SIZE],v[VCOUNT][DIM];

int readfield(){

	int i, k;


	FILE  *ifp;

	ifp = fopen("field8", "r");
	for(k=0; k < SIZE; k++)
		for(i=0; i < SIZE; i++)
		{
			fscanf(ifp, "%d", &m[k][i]);
		}

	for(k=0; k < SIZE; k++)
		for(i=0; i < SIZE; i++)
		{
			fscanf(ifp, "%d", &a[k][i]);
		}
	fclose(ifp);

	return(0);
}
 int readVectors() {

	FILE  *ifpV;

	ifpV = fopen("list2", "r");
	for(int k=0; k < VCOUNT ; k++)
		for(int i=0; i < DIM ; i++)
		{
		fscanf(ifpV, "%d", &v[k][i]);
		}

	fclose(ifpV);

return(0);
 }



int add(int x, int y)

  {
    return(a[x][y]);
  }

int mult(int x, int y)

  {
    return(m[x][y]);
  }
int vect(int x, int y)

  {
    return(v[x][y]);
  }
