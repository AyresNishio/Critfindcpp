#include <iostream>
using namespace std;



double* multiply(double mat1[], double mat2[], int m, int n, int p)
{
	int i, j, k;
	double *mul;
	mul = (double*)malloc(m*p*sizeof(double));
	for (i = 0; i<m; i++)
	{
		for (j = 0; j<p; j++)
		{
			mul[i*p + j] = 0;
			for (k = 0; k<n; k++)
			{
				//printf("%3.4lf * %3.4lf ; ", mat1[i*n + k], mat2[k*p + j]);
				mul[i*p + j] = mul[i*p + j] + mat1[i*n + k] * mat2[k*p + j];
			}
			//printf("\n");
		}

	}

	return mul;
}

double* invert(double mat[], int m,int *inv) {
	//monta inv eye(m)
	double *Inv;
	Inv = (double*)malloc(m*m * sizeof(double));

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {

			if (i == j) Inv[i*m + j] = 1;
			else Inv[i*m + j] = 0;

		}

	}
	double pivo = 0.;
	double li = 0.;
	for (int i = 0; i < m; i++) {
		pivo = mat[i*m + i];
		li = 0.;
		for (int k = i; k < m; k++) {
			li = li + mat[i*m + k];
		}

		if (abs(pivo)<0.0001||abs(li)<0.0001) {
			*inv = 0;
			printf("matriz não possui inversa");
			return Inv;
		}
		
		for (int j = 0; j < m; j++) {
			mat[i*m + j] = mat[i*m + j] / pivo;
			Inv[i*m + j] = Inv[i*m + j] / pivo;
			
		}


		for (int j = 0; j < m; j++) {
			if (j!=i){
				pivo = mat[j*m + i];
				for (int k = 0; k < m; k++) {
					mat[j*m + k] = mat[j*m + k] - pivo*mat[i*m + k];
					Inv[j*m + k] = Inv[j*m + k] - pivo*Inv[i*m + k];
					
				}
			}
		}



	}

	*inv = 1;
	return Inv;
}

void printmat(double mat[], int m, int n) {
	for (int i = 0; i < m; i++) {
		printf("[");
		for (int j = 0; j < n; j++) {
			printf("%3.4lf ", mat[i*n + j]);
		}
		printf("]\n");
	}
}

int main()
{
	const int nmed = 9;
	const int nbar = 6;
	const int ref = 1;

	//const int nmed = 6, nbar=3;
	int A[nbar*nbar] = {
		0, 1, 1, 0, 1, 0,
		1, 0, 1, 0, 0, 0,
		1, 1, 0, 1, 0, 0,
		0, 0, 1, 0, 1, 1,
		1, 0, 0, 1, 0, 0,
		0, 0, 0, 1, 0, 0 };

	int medida[nmed * 3] = { 1, 2, 1,
		                     2, 3, 1,
							 4, 5, 1,
		                     4, 6, 1,
		                     5, 4, 1,
		                     1, 1, 1,
		                     3, 3, 1,
		                     5, 5, 1,
		                     6, 6, 1};
	
	const int nTmed = nmed + 1;
	double H[nTmed* nbar];
	double Ht[nTmed*nbar];

	//CALCULO JACOBIANA
	for (int i = 0; i < nmed; i++) {
		for (int j = 0; j < nbar; j++) {
			int bar = j + 1;
			int de = medida[i * 3];
			int para = medida[i * 3 + 1];
			int ind = i*nbar + j;
			if (de != para) {
				if (de == bar) {
					H[ind] = (double)1;

				}
				else if (para == bar) {
					H[ind] = (double)-1;

				}
				else {
					H[ind] = 0;

				}
			}
			else if (de==para) {
				
				if (de == bar) {
					H[ind] = 0;
					for (int k = 0; k < nbar; k++) {
						H[ind] += A[(bar-1)*nbar + k];
					}

				}
				else if (A[(bar-1)*nbar + de-1] == 1) {
					H[ind] = (double)-1;
				}
				else {
					H[ind] = 0;
				}
			
			}

		}

	}
	for (int i = 0; i < nbar; i++) {
		if (i+1 == ref) { H[nmed*nbar + i] = 1; }
		else { H[nmed*nbar + i] = 0; }

	}

	printf("-------------------------\nJacobiana(H)\n");
	printmat(H, nTmed, nbar);


	//-----------

	double Z[nTmed];//Matriz de medidas
	for (int i = 0; i < nTmed; i++) {
		Z[i] = 1;
	}

	double R[nTmed*nTmed];//Matriz de desvios
	for (int i = 0; i < nTmed; i++) {
		for (int j = 0; j < nTmed; j++) {

			if (i == j) R[i*(nTmed) + j] = 1;
			else R[i*(nTmed) + j] = 0;

		}

	}


	//------TRANSPOSE H------
	for (int i = 0; i < nbar; i++) {
		for (int j = 0; j < nTmed; j++) {
			Ht[i*nTmed + j] = H[j*nbar + i];
		}
	}
	printf("-------------------------\nJacobiana transposta(Ht)\n");
	printmat(Ht, nbar, nTmed);
	//-----------------------
	
	double *E,*G;
	int inv=1;
	int *pinv;
	pinv = &inv;

	E = (double*)malloc(nTmed*nTmed * sizeof(double));
	G = (double*)malloc(nbar*nbar * sizeof(double));
	G = multiply(Ht, H, nbar, nTmed, nbar);
	printf("-------------------------\nMatriz de Ganho G\n");
	printmat(G, nbar, nbar);
	G = invert(G, nbar, pinv);
	if (inv) {
		E = multiply(H, G, nTmed, nbar, nbar);
		
		E = multiply(E, Ht, nTmed, nbar, nTmed);
		for (int i = 0; i < nTmed*nTmed; i++) {
			E[i] = R[i] - E[i];
		}
		
		printf("-------------------------\nE=R-H*((transpose(H)*H) ^ -1)*transpose(H)\n");
		printmat(E, nTmed, nTmed);
	}
	
	return 0;
}

