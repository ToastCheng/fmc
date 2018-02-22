#include "cuda_runtime.h"  //#include <cutil.h>
#include "device_launch_parameters.h"

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <curand.h>
#include <curand_kernel.h>
extern char* NUM_SPECTRUM;

using namespace std;


// DEFINES

#define NUM_BLOCKS 56 //Keep numblocks a multiple of the #MP's of the GPU (8800GT=14MP)
//The register usage varies with platform. 64-bit Linux and 32.bit Windows XP have been tested.

#ifdef __linux__ //uses 25 registers per thread (64-bit)
	#define NUM_THREADS_PER_BLOCK 320 //Keep above 192 to eliminate global memory access overhead However, keep low to allow enough registers per thread
	#define NUM_THREADS 17920
#endif

#ifdef _WIN32 //uses 26 registers per thread
	#define NUM_THREADS_PER_BLOCK 288 //Keep above 192 to eliminate global memory access overhead However, keep low to allow enough registers per thread
	#define NUM_THREADS 16128
#endif




#define NUMSTEPS_GPU       10000
#define PI                 3.141592654f
#define RPI                0.318309886f
#define MAX_LAYERS         10
// no used
#define STR_LEN            200
// 1: normal, 0: oblique
#define NORMAL             1
#define NUM_OF_DETECTOR    (NORMAL ? 6:9)			// normal: 6 fibers, oblique: 9 fibers
#define ANGLE              (NORMAL ? 0:45)			// normal: 0 degree, oblique: 45 degree
#define NAOfSource         (NORMAL ? 0.26:0.22)     // normal: 0.26, oblique: 0.22
#define NAOfDetector       (NORMAL ? 0.26:0.22)     // normal: 0.26, oblique: 0.22
#define n_detector         1.457f
#define n_source           1.457f
#define illumination_r     0.01f
#define collect_r          0.01f
#define NUMBER_PHOTONS     10000000				//1000000 (20151029)
#define NUM_OF_SIMULATION  26						//21 (20151029)
#define WEIGHTI			   0.0001f					// weight defined too small
#define CHANCE			   0.1f


// TYPEDEFS
typedef struct __align__(16)
{
	float z_min;		// Layer z_min [cm]
	float z_max;		// Layer z_max [cm]
	float mutr;			// Reciprocal mu_total [cm]
	float mua;			// Absorption coefficient [1/cm]
	float mutrE;        // Reciprocal mu_total (excitation) [cm]
	float muaE;         // Absorption coefficient (excitation)[1/cm]
	float g[181];		// Anisotropy factor [-]  //excitation
	float gE[181];		// Anisotropy factor [-]  //emission
	float n;			// Refractive index [-]
	float f_mua;        // fluoro's mua at excitation wavelength
	float quantum;      // quantum yield of the fluorophore in this layer
	float emission_p;   // emission probability at particular wavelength
}LayerStruct;

typedef struct __align__(16)
{
	float x;				// Global x coordinate [cm]
	float y;				// Global y coordinate [cm]
	float z;				// Global z coordinate [cm]
	float dx;				// (Global, normalized) x-direction
	float dy;				// (Global, normalized) y-direction
	float dz;				// (Global, normalized) z-direction
	float weight;			// Photon weight
	int	  layer;			// Current layer
	bool  fluoro;			// true: if it's a fluoro
	int   fluoro_layer;		// fluoro is generated in which layer
}PhotonStruct;

typedef struct
{
	unsigned long number_of_photons;
	unsigned int n_layers;
	float start_weight;
	LayerStruct* layers;
}SimulationStruct;

typedef struct
{
	// at most, for oblique configuration, there are 13 fibers
	float radius[13];
	float NA[13];
	float position[13];
	float angle[13];
	float up_data[13];    // fluoro comes from upper layer
	float down_data[13];  // fluoro comes from bottom layer
	float Exphoton[1];
	float Fluoro[1];

}Fibers;

typedef struct
{
	Fibers* f;
	PhotonStruct* p;						// Pointer to structure array containing all the photon data
	unsigned int* thread_active;			// Pointer to the array containing the thread active status
	unsigned int* num_terminated_photons;	// Pointer to a scalar keeping track of the number of terminated photons
	curandState*  state;
}MemStruct;

__device__ __constant__ unsigned int num_photons_dc[1];
__device__ __constant__ unsigned int n_layers_dc[1];
__device__ __constant__ float start_weight_dc[1];
__device__ __constant__ LayerStruct layers_dc[MAX_LAYERS];
