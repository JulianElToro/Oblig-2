#include<iostream > 
#include<vector > 
#include<cstdlib > 
#include<math.h> 
#include<complex> 
#include <fstream>

using namespace std;

#define PI 3.141592653589793238462643

const int q = 3;	//	q spin states
const int L = 16;	//	linear system size
const double  T = 0.25;	//	temperature in units of J 


const int  N = L*L;	// total number of spins
const double   pconnect = 1 - exp(-1 / T);  // connection probability 


const int  NCLUSTERS = 1;  // #  of   cluster builds in one MC step .
const int NESTEPS = 10000; // # of equilibrium MC step . 
const int NMSTEPS = 10000; // # of measurement MC step . 
const int NBINS = 10; // #of measurement bins


int S[L][L];	//	the spin  array matrix
vector <int> M(q);	//	number of spins in each state (for q=3, the states 0,1,2) 

vector <complex<double>  > W(q);   //  order parameter weights

//	lattice handling :
enum dirs { RIGHT, LEFT, UP, DOWN };
int  indx(int  x) { return  x; }   //  make an indx on every site
int  xpos(int	i) { return i % L; }
int  ypos(int	i) { return i % L; }

int Nbr(int  i, int  dir)
{
	int  x = xpos(i);
	int y = ypos(i);
	switch (dir)
	{
	case  RIGHT:  return  indx((x + 1) % L);
	case LEFT:	return indx((x - 1 + L) % L);
	case  UP:  return  indx((y + 1) % L);
	case DOWN:	return indx((y - 1 + L) % L);
	}
}


void  FlipandBuildFrom(int s, int e)
{
	int  oldstate(S[s][e]), newstate((S[s][e] + 1) % q);

	S[s][e] = newstate;   //	flip spin

	M[oldstate]--;

	M[newstate]++;	//  update spin counts

	for (int dir = 0; dir < 4; dir++)  //  go thru neighbors
	{
		int j1 = Nbr(s, dir);
		int j2 = Nbr(e, dir);

		if (S[j1][j2] == oldstate)
			if (rand() / (RAND_MAX + 1.) < pconnect) { FlipandBuildFrom(j1,j2); }

		// rand() is a pseudorandom number between 0 and RAND_MAX
	}
}

int  main()
{
	//	initialize	order parameter weights
	for (int s = 0; s < q; s++) {
		W[s] = complex<double>(cos(2 * PI * s / q), sin(2 * PI * s / q));
	}
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			S[i][j] = 0; //	initialize to the spin=0  state
		}
	}
	for (int s = 1; s < q; s++) {
		M[s] = 0; //	initialize	counters .
	}

	M[0] = N;
	// All spins start with state 0
	srand((unsigned)time(0));   //	initialize random number gen .


	//  equilibriate
	for (int t = 0; t < NESTEPS; t++)
	{
		for (int c = 0; c < NCLUSTERS; c++)
		{
			FlipandBuildFrom(rand() % L , rand() % L);
		}
	}

	ofstream ofile;
	ofile.open("mvsTJ.txt");
	ofile << scientific;

	// measure NBINS times
	for (int n = 0; n < NBINS; n++)
	{
		complex<double> m(0., 0.);

		double  m1 = 0, C = 0, m2 = 0, m4 = 0;  //  measurement results

		for (int t = 0; t < NMSTEPS; t++)
		{
			for (int c = 0; c < NCLUSTERS; c++) { FlipandBuildFrom(rand() % L, rand() % L); }
				
			complex<double> tm(0., 0.);

			for (int s = 0; s < q; s++) {
				
				tm += W[s] * double(M[s]);

				tm /= N;
				
				double tm1 = abs(tm); // modulus
				double tm2 = tm1 * tm1;
				m += tm; // average magnetization (complex)
				m1 += tm1; // average magnetization (real)	
				m2 += tm2; // magnetization squared
				m4 += tm2 * tm2; // average magnetization to the forth
			}
		}

		m /= NMSTEPS; m1 /= NMSTEPS;  m2 /= NMSTEPS; m4 /= NMSTEPS;

		//cout << m << " " << m1 << " " << m2 << " " << m4 << endl;

		ofile << m << " " << m1 << " " << m2 << " " << m4 << endl;

	}

	ofile.close();
}