// from https://www.olcf.ornl.gov/tutorials/cuda-monte-carlo-pi/
// 03/06/18

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand.h>
#include <time.h>
 
//Declare the CUDA kernel
__global__ void kernel(int* count_d, float* randomnums)
{
        int i;
        double x,y,z;
        //Find the overall ID of the thread
        int tid = blockDim.x * blockIdx.x + threadIdx.x;
        i = tid;
        int xidx = 0, yidx = 0;
 
        //Do the MonteCarlo!
        xidx = (i+i);
        yidx = (xidx+1);
 
        //Get the random x,y points
        x = randomnums[xidx];
        y = randomnums[yidx];
        z = ((x*x)+(y*y));
 
        if (z<=1)
                count_d[tid] = 1;
        else
                count_d[tid] = 0;
}
 
//Used to check if there are any errors launching the kernel
void CUDAErrorCheck()
{
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess)
        {
                printf("CUDA error : %s (%d)\n", cudaGetErrorString(error), error);
                exit(0);
        }
}
 
//int main(int argc,char* argv[])
//calling C++ from fortran solution from http://theochem.mercer.edu/interlanguage/testdot.cc
extern "C" {
        int picuda_(void);
        }

int picuda_(void)
{
        //NOTE: if threads and/or blocks is changed, niter needs to be changed to reflect
        //that change (niter=threads*blocks)
        int niter = 100000;
        float *randomnums;
        double pi;
        //Allocate the array for the random numbers
        cudaMalloc((void**)&randomnums, (2*niter)*sizeof(float));
        //Use CuRand to generate an array of random numbers on the device
        int status;
        curandGenerator_t gen;
        status = curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MRG32K3A);
        status |= curandSetPseudoRandomGeneratorSeed(gen, 4294967296ULL^time(NULL));
        status |= curandGenerateUniform(gen, randomnums, (2*niter));
        status |= curandDestroyGenerator(gen);
        //Check to see if there was any problem launching the CURAND kernels and generating
        //the random numbers on the device
        if (status != CURAND_STATUS_SUCCESS)
        {
                printf("CuRand Failure\n");
                exit(EXIT_FAILURE);
        }
 
        //Threads per thread block to be launched
        int threads = 1000;
        //Number of thread blocks launched
        int blocks = 100;
        int* count_d;
        int *count = (int*)malloc(blocks*threads*sizeof(int));
        unsigned int reducedcount = 0;
        //Allocate the array to hold a value (1,0) whether the point in is the circle (1) or not (0)
        cudaMalloc((void**)&count_d, (blocks*threads)*sizeof(int));
        CUDAErrorCheck();
        //Launch the kernel
        kernel <<<blocks, threads>>> (count_d, randomnums);
        //Acts as a kind of code Barrier until the kernel is finished. Kernel calls are nonblocking so
        //the code would continue regardless of whether the kernel succeeded or not without the Sync
        cudaDeviceSynchronize();
        CUDAErrorCheck();
        //Copy the resulting array back
        cudaMemcpy(count, count_d, blocks*threads*sizeof(int), cudaMemcpyDeviceToHost);
        int i = 0;
 
        //Reduce array into int
        for(i = 0; i<niter; i++)
                reducedcount += count[i];
 
        //Free the cudaMalloc()'d arrays
        cudaFree(randomnums);
        cudaFree(count_d);
        free(count);
 
        //Find the ratio
        pi = ((double)reducedcount/niter)*4.0;
        printf("Pi cuda: %f\n", pi);
 
        return 0;
}