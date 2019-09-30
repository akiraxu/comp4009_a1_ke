comp4009_a1_ke
100955433

Usage:

main <1 for q1, 2 for q3> <iteration> <parallelism>

e.g.
./main 1 1440 100
run question one way for question 2 input for 1440 iteration, with sub case of size 100

./main 2 100 500
run question one way for question 2 input for 100 iteration, with 500 sub cases


Files:

_compile: compile code
_run: not used
main.cpp: all code here

For question 2
a1_input.txt
a1_sample_output.txt
a1_sample_visual.txt

Methods:
struct a1_data: Store all need information
struct pos: Position
void updateDIM(a1_data * data): Update global var
void freeMem(a1_data * data): Delete memory
void sawpData(a1_data * data): Swap two grid
double calcHeat(double * arr, int len): Calcuate total heat for varification
int modulo(int a, int q): Better modulo than c %
bool inRange(int n, int a, int b): Does n in [a, b) where a < b
int arrayPos(int x, int y, int z): Transform xyz to i
void cartesianPos(int i, int * x, int * y, int * z): Transform i to xyz
int arrayTrans(int i, int dx, int dy, int dz): Get new i with delta
int arrayTrans_v2(int i, int k, int a, int b): Better than two-step transform
void fillData(a1_data * data, pos a, pos b, double num): Fill data with given cornor point
double findMin(double * arr, int len): Find minimum value
double findMax(double * arr, int len): Find maximum value
string printLite(a1_data * data, int q, double min, double max): Generate readable string with modulo index
void stencil(double * grid_a, double * grid_b, double alpha, double beta): The standard way, using xyz
void stencil_v2(a1_data * data): Using i
void stencil_naive(double * grid_a, double * grid_b, double alpha, double beta): Another
void stencil_kernel(a1_data * data, int i): The kernel;
void stencil_kernel_v2(a1_data * data, int i): The kernel with better calcuation
void stencil_parallel_naive(a1_data * data, int op, int ed, int sub): Parallel into chunks
void stencil_parallel_better(a1_data * data, int op, int ed, int sub): Parallel with cache-friendly
string printArray(double * arr, int len): Output array into string
void init(a1_data * data): Read and init form file
void result(a1_data * data): Write result into file
void fillQ2(a1_data * data): Fill question 2 data
string printData(a1_data * data): Print status
