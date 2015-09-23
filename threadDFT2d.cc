// Threaded two-dimensional Discrete FFT transform
// Guru Das Srinagesh
// ECE6122 Project 2


#include <iostream>
#include <string>
#include <math.h>

#include "Complex.h"
#include "InputImage.h"

#include <stdio.h>
#include <pthread.h>

// You will likely need global variables indicating how
// many threads there are, and a Complex* that points to the
// 2d image being transformed.

Complex* ImageData;
int ImageWidth;
int ImageHeight;

#define N_THREADS 16

int test_variable = 0;              // To test how pThreads works

int N=8; // Number of points in the 1-D transform

// pThreads variables
pthread_mutex_t startCountMutex;    // To test how pThreads works
pthread_mutex_t exitMutex;          // To test how pThreads works
pthread_mutex_t elementCountMutex;  // To test how pThreads works
pthread_mutex_t printfMutex;        // To test how pThreads works
pthread_cond_t  exitCond;           // To test how pThreads works
int             startCount;         // To test how pThreads works

Complex* W;     // Twiddle bits

using namespace std;

// Function to reverse bits in an unsigned integer
// This assumes there is a global variable N that is the
// number of points in the 1D transform.
unsigned ReverseBits(unsigned v)
{ //  Provided to students
  unsigned n = N; // Size of array (which is even 2 power k value)
  unsigned r = 0; // Return value
   
  for (--n; n > 0; n >>= 1)
    {
      r <<= 1;        // Shift return value
      r |= (v & 0x1); // Merge in next bit
      v >>= 1;        // Shift reversal value
    }
  return r;
}

// GRAD Students implement the following 2 functions.
// Undergrads can use the built-in barriers in pthreads.

// Call MyBarrier_Init once in main
void MyBarrier_Init()// you will likely need some parameters)
{
}

// Each thread calls MyBarrier after completing the row-wise DFT
void MyBarrier() // Again likely need parameters
{
}
                    
void precomputeW()
{
  W = new Complex[ImageWidth]; 

  /* Compute W only for first half */
  for(int n=0; n<(ImageWidth/2); n++){
    W[n].real = cos(2*M_PI*n/ImageWidth);
    W[n].imag = -sin(2*M_PI*n/ImageWidth);
  }

  ///* W[n + N/2] = -W[n] */
  //for(int n=(ImageWidth/2); n < N; n++){
  //  W[n].real = - W[n - (N/2)].real;
  //  W[n].imag = - W[n - (N/2)].imag;
  //}

}

void Transform1D(Complex* h, int N)
{
  // Implement the efficient Danielson-Lanczos DFT here.
  // "h" is an input/output parameter
  // "N" is the size of the array (assume even power of 2)
  
  /* Reorder array based on bit reversing */
  for(int i=0; i<(8/2); i++){
    int rev_i = ReverseBits(i);
    printf("%d swapped with %d\n", i, rev_i);
    Complex temp = h[i];
    h[i] = h[rev_i];
    h[rev_i] = temp;
  }

  //for(int i=0; i<ImageWidth; i++)
  //  cout << h[i].Mag() << " ";
  //printf("\n");
 
  for(int pt=2; pt <= N; pt *= 2){    // Number of points in transform
    for(int j=0; j < (N/*/pt*/); j/*++*/+=pt){    // Number of times to execute the mults
      for(int k=0; k < (pt/2); k++){  // Number of times to execute the mults
        int offset = pt/2;
        //printf("H[%d] W[%d] H[%d]\t H[%d] W[%d] H[%d]  |  ", (j*pt+k), (k*N/pt)%(N/pt), (j*pt+k+offset), (j*pt+k), (((k+offset)*N/pt)%(N/pt)), (j*pt+k+offset));
        printf("H[%d] W[%d] H[%d]\t H[%d] W[%d] H[%d]  |  ", (j+k), (k*N/pt), (j+k+offset), (j+k), (k*N/pt), (j+k+offset));
        Complex oldfirst = h[j+k];
        Complex oldsecond = h[j+k+offset];
        h[j+k] = oldfirst + W[k*N/pt]*oldsecond;
        h[j+k+offset] = oldfirst - W[k*N/pt]*oldsecond;
        
        //printf("H[%d] W[%d] H[%d]\t H[%d] W[%d] H[%d]  |  ", (j*pt+k), ((j*pt+k)*N/pt)%(N/pt), (j*pt+k+offset), (j*pt+k), (int)((j*pt+k+offset)*N/pt)%(N/pt), (j*pt+k+offset));
      }
      printf("\n");
    }
    printf("\n");
  }

  /* Print out the 1-D transformed array */
  for(int i=0; i<8; i++)
    cout<<"("<<h[i].real<<", "<<h[i].imag<<")"<<endl;
}

void* Transform2DTHread(void* v)
{ // This is the thread starting point.  "v" is the thread number
  // Calculate 1d DFT for assigned rows
  // wait for all to complete
  // Calculate 1d DFT for assigned columns
  // Decrement active count and signal main if all complete

  /* Determine thread ID */
  unsigned long thread_id = (unsigned long)v;

  /* Determine starting row and number of rows per thread */
  int rowsPerThread = ImageHeight / N_THREADS;
  int startingRow = thread_id * rowsPerThread;

  //pthread_mutex_lock(&elementCountMutex);
  //test_variable++;
  //pthread_mutex_unlock(&elementCountMutex);

  pthread_mutex_lock(&printfMutex);
  printf("  Thread %2ld: Just chillin' \n", thread_id);
  pthread_mutex_unlock(&printfMutex);

  pthread_mutex_lock(&startCountMutex);
  startCount--;

  if (startCount == 0)
  { // Last to exit, notify calling function
    pthread_mutex_unlock(&startCountMutex);
    pthread_mutex_lock(&exitMutex);
    pthread_cond_signal(&exitCond);
    pthread_mutex_unlock(&exitMutex);
  }
  else
  {
    pthread_mutex_unlock(&startCountMutex);
  }

  return 0; 
}

void Transform2D(const char* inputFN) 
{ 
  /* Do the 2D transform here. */

  InputImage image(inputFN);  // Create the helper object for reading the image
  ImageWidth = image.GetWidth();
  ImageHeight = image.GetHeight();

  // All mutex and condition variables must be "initialized"
  pthread_mutex_init(&exitMutex,0);
  pthread_mutex_init(&startCountMutex,0);
  pthread_mutex_init(&elementCountMutex,0);
  pthread_mutex_init(&printfMutex,0);
  pthread_cond_init(&exitCond, 0);

  // Create the global pointer to the image array data
  ImageData = image.GetImageData();

  // Precompute W values
  precomputeW();

  /* Perform the 1-D transform on the first row */
  Transform1D(&ImageData[0], 8);

  //// Hold the exit mutex until waiting for exitCond condition
  //pthread_mutex_lock(&exitMutex);

  //startCount = N_THREADS;
  //
  //pthread_t threads[N_THREADS];

  //int i = 0;

  //// Create 16 threads
  //for(i=0; i < N_THREADS; ++i){
  //  pthread_create(&threads[i], 0, Transform2DTHread, (void *)i);
  //}

  //// Wait for all threads complete
  //pthread_cond_wait(&exitCond, &exitMutex);

  //// Write the transformed data
  //printf("\n  Final result: test_variable = %d\n", test_variable);

}

int main(int argc, char** argv)
{
  string fn("testimage.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line

  Transform2D(fn.c_str()); // Perform the transform.
}  
