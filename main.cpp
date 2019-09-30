
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cilk/cilk.h>
#include "hwtimer.h"
using namespace std;

bool DEBUG = false;

int DIM_X;
int DIM_Y;
int DIM_Z;

long timecost = 0;

struct a1_data {
	double * a;
	double * b;
	double alpha;
	double beta;
	double ax;
	double ay;
	double bx;
	double by;
	int dim_x;
	int dim_y;
	int dim_z;
	int n;
	int round;
	int sum_x;
	int sum_y;
	int sum_z;
};

struct pos {
	int x;
	int y;
	int z;
};

void updateDIM(a1_data * data){
	DIM_X = data->dim_x;
	DIM_Y = data->dim_y;
	DIM_Z = data->dim_z;
}

void freeMem(a1_data * data){
	delete data->a;
	delete data->b;
}

void sawpData(a1_data * data){
	double * temp = NULL;
	temp = data->a;
	data->a = data->b;
	data->b = temp;
}

double calcHeat(double * arr, int len){
	double total = 0;
	for(int i = 0; i < len; i++){
		total += arr[i];
	}
	return total;
}

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

bool inRange(int n, int a, int b){
	if(a == b) return false;
	if(a > b) return n < a && n >= b;
	return n >= a && n < b;
}

int arrayPos(int x, int y, int z){
	if(DEBUG)cout << "AP: x:" << x << " y:" << y << " z:" << z << endl;
	return modulo(x,DIM_X) + modulo(y,DIM_Y) * DIM_X + modulo(z,DIM_Z) * DIM_X * DIM_Y;
}

void cartesianPos(int i, int * x, int * y, int * z){
	*z = i / (DIM_X * DIM_Y);
	*y = (i - *z * DIM_X * DIM_Y) / DIM_X;
	*x = i - *z * DIM_X * DIM_Y - *y * DIM_X;
	int pos = arrayPos(*x,*y,*z);
	if(DEBUG)cout << "CP: x:" << *x << " y:" << *y << " z:" << *z << " (" << i << " --> " << pos << ")" << endl;
}

int arrayTrans(int i, int dx, int dy, int dz){
	int x, y, z;
	cartesianPos(i, &x, &y, &z);
	return arrayPos(modulo(x + dx,DIM_X), modulo(y + dy,DIM_Y), modulo(z + dz,DIM_Z));
}
int arrayTrans_v2(int i, int k, int a, int b){
	return (modulo(modulo(i, b) + k * a, b)) + (i / b);
}

void fillData(a1_data * data, pos a, pos b, double num){
	for(int z = 0; z < data->dim_z; ++z){
		if(inRange(z, a.z, b.z)){
			for(int y = 0; y < data->dim_y; ++y){
				if(inRange(y, a.y, b.y)){
					for(int x = 0; x < data->dim_x; ++x){
						if(inRange(x, a.x, b.x)){
							data->b[arrayPos(x,y,z)] = num;
						}
					}
				}
			}
		}
	}
}

double findMin(double * arr, int len){
	double min = arr[0];
	for(int i = 0; i < len; i++){
		if(arr[i] < min){
			min = arr[i];
		}
	}
	return min;
}

double findMax(double * arr, int len){
	double max = arr[0];
	for(int i = 0; i < len; i++){
		if(arr[i] > max){
			max = arr[i];
		}
	}
	return max;
}

string printLite(a1_data * data, int q, double min, double max){
	char scale[15];
	strcpy(scale, " .:-=+*#%@");
	ostringstream oss;
	oss << "----------ROUND: " << data->round << "----------" << endl;
	oss << "----------SIZE: " << data->n << " TOTAL HEAT: " << calcHeat(data->a, data->n) << "----------\n";
	oss << "----------MIN " << min << "  .:-=+*#%@ " << max << " MAX----------\n";
	for(int z = 0; z < data->dim_z; ++z){
		if(z % q == 0){
			for(int y = 0; y < data->dim_y; ++y){
				if(y % q == 0){
					for(int x = 0; x < data->dim_x; ++x){
						if(x % q == 0){
							int a = int((((data->b[arrayPos(x,y,z)] - min) / (max - min)) - 0.01) * 10);
							oss << scale[a];
						}
					}
					oss << "\n";
				}
			}
			oss << "-------------Z-Axis=" << z << "-------------\n";
		}
	}
	return oss.str();
}

void stencil(double * grid_a, double * grid_b, double alpha, double beta){
	for(int z = 0; z < DIM_Z; ++z){
		for(int y = 0; y < DIM_Y; ++y){
			for(int x = 0; x < DIM_X; ++x){
				grid_b[arrayPos(x,y,z)] = alpha * grid_a[arrayPos(x,y,z)] 
				+ beta * (grid_a[arrayPos(x + 1,y,z)] 
				+ grid_a[arrayPos(x - 1,y,z)] 
				+ grid_a[arrayPos(x,y + 1,z)] 
				+ grid_a[arrayPos(x,y - 1,z)] 
				+ grid_a[arrayPos(x,y,z + 1)] 
				+ grid_a[arrayPos(x,y,z - 1)]);
			}
		}
	}
}

void stencil_v2(a1_data * data){
	for(int i = 0; i < DIM_X * DIM_Y * DIM_Z; ++i){
		data->b[i] = data->alpha * data->a[i] 
		+ data->beta * (data->a[arrayTrans(i,+1,0,0)] 
		+ data->a[arrayTrans(i,-1,0,0)] 
		+ data->a[arrayTrans(i,0,+1,0)] 
		+ data->a[arrayTrans(i,0,-1,0)] 
		+ data->a[arrayTrans(i,0,0,+1)] 
		+ data->a[arrayTrans(i,0,0,-1)]);
	}
}

void stencil_naive(double * grid_a, double * grid_b, double alpha, double beta){
	for(int i = 0; i < DIM_X * DIM_Y * DIM_Z; ++i){
		grid_b[i] = alpha * grid_a[i] 
		+ beta * (grid_a[arrayTrans(i,+1,0,0)] 
		+ grid_a[arrayTrans(i,-1,0,0)] 
		+ grid_a[arrayTrans(i,0,+1,0)] 
		+ grid_a[arrayTrans(i,0,-1,0)] 
		+ grid_a[arrayTrans(i,0,0,+1)] 
		+ grid_a[arrayTrans(i,0,0,-1)]);
	}
}

void stencil_kernel(a1_data * data, int i){
	data->b[i] = data->alpha * data->a[i] 
		+ data->beta * (data->a[arrayTrans(i,+1,0,0)] 
		+ data->a[arrayTrans(i,-1,0,0)] 
		+ data->a[arrayTrans(i,0,+1,0)] 
		+ data->a[arrayTrans(i,0,-1,0)] 
		+ data->a[arrayTrans(i,0,0,+1)] 
		+ data->a[arrayTrans(i,0,0,-1)]);
}

void stencil_kernel_v2(a1_data * data, int i){
	data->b[i] = data->alpha * data->a[i] 
		+ data->beta * (data->a[arrayTrans_v2(i,+1,1,data->sum_x)] 
		+ data->a[arrayTrans_v2(i,-1,1,data->sum_x)] 
		+ data->a[arrayTrans_v2(i,+1,data->sum_x,data->sum_y)] 
		+ data->a[arrayTrans_v2(i,-1,data->sum_x,data->sum_y)] 
		+ data->a[arrayTrans_v2(i,+1,data->sum_y,data->sum_z)] 
		+ data->a[arrayTrans_v2(i,-1,data->sum_y,data->sum_z)]);
}

void stencil_parallel_naive(a1_data * data, int op, int ed, int sub){
	if(ed - op <= sub){
		for(int i = op; i < ed; i++){
			stencil_kernel_v2(data,i);
		}
	}else{
		int mid = op + (ed - op) / 2;
		cilk_spawn stencil_parallel_naive(data, op, mid, sub);
		cilk_spawn stencil_parallel_naive(data, mid, ed, sub);
		cilk_sync;
	}
}

void stencil_parallel_better(a1_data * data, int op, int ed, int sub){
	if(ed == op){
		return;
	}else if(ed - op == 1){
		for(int i = op; i < data->n; i += sub){
			stencil_kernel_v2(data,i);
		}
	}else{
		int mid = op + (ed - op) / 2;
		cilk_spawn stencil_parallel_better(data, op, mid, sub);
		cilk_spawn stencil_parallel_better(data, mid, ed, sub);
		cilk_sync;
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

void init(a1_data * data){
	string line;
	ifstream inputFile ("a1_input.txt");
	if(inputFile.is_open()){
		int n = -1;
		int k = 0;
		while(getline(inputFile,line)){
			istringstream aLine(line);
			if(n == -1){
				double ax, ay, bx, by;
				aLine >> data->dim_x >> data->dim_y >> data->dim_z;
				data -> n = data->dim_x * data->dim_y * data->dim_z;
				updateDIM(data);
				aLine >> data->ax >> data->ay >> data->bx >> data->by;
				data -> a = new double[data->n]();
				data -> b = new double[data->n]();
				data -> alpha = data->ax / data->ay;
				data -> beta = data->bx / data->by;
				data -> round = 0;
				data -> sum_x = data->dim_x;
				data -> sum_x = data->dim_x * data->dim_y;
				data -> sum_x = data->dim_x * data->dim_y * data->dim_z;
			}else{
				for(int i = 0; i < data->dim_x; i++){
					aLine >> data->a[k];
					k++;
				}
			}
			n++;
		}
		inputFile.close();
	}
}

void result(a1_data * data){
	string line;
	ofstream outputFile ("a1_output.txt");
	ofstream visualFile ("a1_visual.txt");
	if(outputFile.is_open()){
		outputFile << data->dim_x << " " << data->dim_y << " " << data->dim_z << " ";
		outputFile << data->ax << " " << data->ay << " " << data->bx << " " << data->by << "\n";
		for(int i = 0; i < DIM_X * DIM_Y * DIM_Z; ++i){
			outputFile << data->b[i];
			if(i % DIM_X == DIM_X - 1){
				outputFile << "\n";
			}else{
				outputFile << " ";
			}
		}
		outputFile.close();
	}
	if(visualFile.is_open()){
		double min = findMin(data->b, data->n);
		double max = findMax(data->b, data->n);
		cout << "Min:" << min << " Max:" << max << endl;
		visualFile << printLite(data, 3, min, max) << endl;
		visualFile.close();
	}
}

void fillQ2(a1_data * data){
	freeMem(data);
	data -> alpha = 0.5;
	data -> beta = 1.0/12;
	data -> ax = 1;
	data -> ay = 2;
	data -> bx = 1;
	data -> by = 12;
	data -> dim_x = 100;
	data -> dim_y = 100;
	data -> dim_z = 100;
	data -> n = 1000000;
	data -> round = 0;
	data -> sum_x = data->dim_x;
	data -> sum_x = data->dim_x * data->dim_y;
	data -> sum_x = data->dim_x * data->dim_y * data->dim_z;
	data -> a = new double[data->n]();
	data -> b = new double[data->n]();

	updateDIM(data);

	pos allZero = {.x = 0, .y = 0, .z = 0};
	pos allMax = {.x = DIM_X, .y = DIM_Y, .z = DIM_Z};
	fillData(data, allZero, allMax, 20.0);

	pos hotA = {.x = 75, .y = 75, .z = 75};
	pos hotB = {.x = 85, .y = 85, .z = 85};
	fillData(data, hotA, hotB, 100.0);

	sawpData(data);
}

string printData(a1_data * data){
	ostringstream oss;
	oss << "alpha:" << data->alpha << endl;
	oss << "beta:" << data->beta << endl;
	oss << "ax:" << data->ax << endl;
	oss << "ay:" << data->ay << endl;
	oss << "bx:" << data->bx << endl;
	oss << "by:" << data->by << endl;
	oss << "dim_x:" << data->dim_x << endl;
	oss << "dim_y:" << data->dim_y << endl;
	oss << "dim_z:" << data->dim_z << endl;
	oss << "n:" << data->n << endl;
	oss << "round:" << data->round << endl;
	//cout << printArray(data->a, data->n);
	return oss.str();
}

int main(int argc, char* argv[])
	{
		if (argc != 4) {
			cout << "Usage: main <1 for q1, 2 for q3> <iteration> <parallelism>" << endl;
			return 1;
		}

		hwtimer_t timer;
		initTimer(&timer);

		int method = atoi(argv[1]);

		int iteration = atoi(argv[2]);

		int parallelism = atoi(argv[3]);
/*
		startTimer(&timer);
		int answer = fib(param);
		stopTimer(&timer);
		int fibTime = getTimerNs(&timer);

		cout << "fib(" << param << ") = " << answer << endl;
		cout << "Total time: " << fibTime << "ns" << endl;
*/
		a1_data * data = new a1_data();
/*
		data -> a = new double[27]{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,100};
		data -> b = new double[27]{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		data -> alpha = 0.5;
		data -> beta = 1.0/12;
*/
		init(data);

		cout << printData(data);

		cout << "Total heat before: " << calcHeat(data->a, data->n) << endl;
		
		//cout << printArray(data->a, data->n) << endl;
		//cout << printArray(data->b, data->n) << endl;

//Start calculation
		startTimer(&timer);
		for(int i = 0; i < iteration; i++){
			if(method == 1){
				stencil_parallel_naive(data, 0, data->n, parallelism);
			}else{
				stencil_parallel_better(data, 0, parallelism, parallelism);
			}
			sawpData(data);
			(data->round)++;
		}
		stopTimer(&timer);
//End calculation

		long runTime = getTimerNs(&timer);
		timecost = runTime;
		cout << "Total heat after: " << calcHeat(data->b, data->n) << endl;
		cout << "Total time: " << setprecision(3) << runTime/1000000 << "ms" << endl;

		//cout << printArray(data->b, 27) << endl;

		//fillQ2(data);
		result(data);
		freeMem(data);
		delete data;
		return 0;
	};
