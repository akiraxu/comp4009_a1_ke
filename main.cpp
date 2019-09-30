
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cilk/cilk.h>
#include "hwtimer.h"
using namespace std;

bool DEBUG = false;

int DIM_X;
int DIM_Y;
int DIM_Z;

int fib(int n) 
	{
		if (n < 2)
			return n;
		int a = cilk_spawn fib(n-1);
		int b = cilk_spawn fib(n-2);
		cilk_sync;
		return a + b;
	};

int modulo(int a, int q){
	return a < 0 ? q + a % q : a % q;
}

int arrpos(int x, int y, int z){
	return modulo(x,DIM_X) + modulo(y,DIM_Y) * DIM_X + modulo(z,DIM_Z) * DIM_X * DIM_Y;
}

int stencil(double * grid_a, double * grid_b, double alpha, double beta){
	for(int z = 0; z < DIM_Z; ++z){
		for(int y = 0; y < DIM_Y; ++y){
			for(int x = 0; x < DIM_X; ++x){
				grid_b[arrpos(x,y,z)] = alpha * grid_a[arrpos(x,y,z)] 
				+ beta * (grid_a[arrpos(x + 1,y,z)] 
				+ grid_a[arrpos(x - 1,y,z)] 
				+ grid_a[arrpos(x,y + 1,z)] 
				+ grid_a[arrpos(x,y - 1,z)] 
				+ grid_a[arrpos(x,y,z + 1)] 
				+ grid_a[arrpos(x,y,z - 1)]);
			}
		}
	}
}

string printArray(double * arr, int len){
	ostringstream oss;
	if(DEBUG)cout << "Enter printArray" << endl;
	for(int i = 0; i < len - 1; i++){
		if(DEBUG)cout << "Element: " << arr[i] << " at: " << i << endl;
		oss << fixed << setw(6) << setprecision(3) << arr[i] << " ";
	}
	if(len > 1){
		if(DEBUG)cout << "Last element: " << arr[len - 1] << endl;
		oss << fixed << setw(6) << setprecision(3) << arr[len - 1];
	}
	if(DEBUG)cout << "Exit printArray" << endl;
	return oss.str();
}


int main(int argc, char* argv[])
	{
		if (argc != 2) {
			cout << "Usage: fib <num>" << endl;
			return 1;
		}

		hwtimer_t timer;
		initTimer(&timer);

		int param = atoi(argv[1]);

		startTimer(&timer);
		int answer = fib(param);
		stopTimer(&timer);
		int fibTime = getTimerNs(&timer);

		cout << "fib(" << param << ") = " << answer << endl;
		cout << "Total time: " << fibTime << "ns" << endl;

		double *a = new double[27]{0,0,0,0,0,0,0,0,0,0,0,0,0,27,0,0,0,0,0,0,0,0,0,0,0,0,0};
		double *b = new double[27]{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

		cout << printArray(a, 27) << endl;
		cout << printArray(b, 27) << endl;

		DIM_X = DIM_Y = DIM_Z = 3;

		for(int i = 0; i < 100; i++){
			stencil(a, b, 0.5, 1.0/12);
			a = b;
		}

		double total = 0;
		for(int i = 0; i < DIM_X * DIM_Y * DIM_Z + 1; i++){
			total += b[i];
		}

		cout << printArray(b, 27) << endl;
		cout << "total: " << total << endl;

		return 0;
	};
