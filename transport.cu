#include "header.h"
#include <float.h> //for FLT_MAX


int InitMemStructs(MemStruct* HostMem, MemStruct* DeviceMem, SimulationStruct* sim);
void FreeMemStructs(MemStruct* HostMem, MemStruct* DeviceMem);
void FreeSimulationStruct(SimulationStruct* sim, int n_simulations);
__global__ void MCd(MemStruct DeviceMem, unsigned long long seed, int *n);
__global__ void LaunchPhoton_Global(MemStruct DeviceMem);
int InitDCMem(SimulationStruct* sim);
int Write_Simulation_Results(MemStruct* HostMem, SimulationStruct* sim, clock_t simulation_time);
int read_simulation_data(char* filename, SimulationStruct** simulations, int ignoreAdetection);
int interpret_arg(int argc, char* argv[], unsigned long long* seed, int* ignoreAdetection);

__device__ void LaunchPhoton(PhotonStruct* p, curandState *state);
__global__ void LaunchPhoton_Global(MemStruct DeviceMem, unsigned long long seed);
__device__ void fluoro_MC(PhotonStruct* p, curandState *state,Fibers *f);
__device__ void Spin(PhotonStruct*, float*,curandState* state);
__device__ unsigned int Reflect(PhotonStruct*, int, curandState* state);
__device__ unsigned int PhotonSurvive(PhotonStruct*, curandState* state);
__device__ void AtomicAddULL(unsigned long long* address, unsigned int add);
__device__ void detect(PhotonStruct* p, Fibers* f);
__device__ int binarySearch(float *data, float value);
void fiber_initialization(Fibers* f);
void output_fiber(SimulationStruct* sim, float *up, float* down, float *Exp, float *F ,int n);
void calculate_reflectance(Fibers* f, float *up, float *down, float *Exp, float *F);

__device__ float rn_gen(curandState *s)
{
	float x = curand_uniform(s);
    return x;
}

void DoOneSimulation(SimulationStruct* simulation, int index)
{
	unsigned long long seed = time(NULL);
	float up_fluorescence[NUM_OF_DETECTOR] = {0};
	float down_fluorescence[NUM_OF_DETECTOR] = {0};
	float ExPhoton[1] = {0};
	float Fluoro[1] = {0};

	MemStruct DeviceMem;
	MemStruct HostMem;
	unsigned int threads_active_total=1;
	unsigned int i,ii;
	int  H_num_of_fluoro[1] = {0};
	int  *D_num_of_fluoro;
	cudaMalloc(&D_num_of_fluoro, sizeof(int));
	cudaMemcpy(D_num_of_fluoro, H_num_of_fluoro, sizeof(int),cudaMemcpyHostToDevice);

    cudaError_t cudastat;

	InitMemStructs(&HostMem,&DeviceMem,simulation);
	InitDCMem(simulation);

    dim3 dimBlock(NUM_THREADS_PER_BLOCK);
    dim3 dimGrid(NUM_BLOCKS);

	LaunchPhoton_Global<<<dimGrid,dimBlock>>>(DeviceMem, seed);
	cudaThreadSynchronize(); //CUDA_SAFE_CALL( cudaThreadSynchronize() ); // Wait for all threads to finish
	cudastat=cudaGetLastError(); // Check if there was an error
	if(cudastat)printf("Error code=%i, %s.\n",cudastat,cudaGetErrorString(cudastat));

	i=0;

	while(threads_active_total>0)
	{
		i++;
		fiber_initialization(HostMem.f);
	    cudaMemcpy(DeviceMem.f,HostMem.f,NUM_THREADS*sizeof(Fibers),cudaMemcpyHostToDevice);

		//run the kernel
		seed = time(NULL);
		MCd<<<dimGrid,dimBlock>>>(DeviceMem, seed, D_num_of_fluoro);
		cudaThreadSynchronize(); //CUDA_SAFE_CALL( cudaThreadSynchronize() ); // Wait for all threads to finish
		cudastat=cudaGetLastError(); // Check if there was an error
		if(cudastat)printf("Error code=%i, %s.\n",cudastat,cudaGetErrorString(cudastat));

		// Copy thread_active from device to host, later deleted
		cudaMemcpy(HostMem.thread_active, DeviceMem.thread_active, NUM_THREADS*sizeof(unsigned int), cudaMemcpyDeviceToHost); //CUDA_SAFE_CALL(cudaMemcpy(HostMem.thread_active,DeviceMem.thread_active,NUM_THREADS*sizeof(unsigned int),cudaMemcpyDeviceToHost) );
		threads_active_total = 0;
		for(ii=0;ii<NUM_THREADS;ii++) threads_active_total+=HostMem.thread_active[ii];

		cudaMemcpy(HostMem.num_terminated_photons, DeviceMem.num_terminated_photons, sizeof(unsigned int), cudaMemcpyDeviceToHost); //CUDA_SAFE_CALL(cudaMemcpy(HostMem.num_terminated_photons,DeviceMem.num_terminated_photons,sizeof(unsigned int),cudaMemcpyDeviceToHost) );

		//printf("Run %u, Number of photons terminated %u, Threads active %u\n",i,*HostMem.num_terminated_photons,threads_active_total);

		cudaMemcpy(HostMem.f, DeviceMem.f, NUM_THREADS*sizeof(Fibers), cudaMemcpyDeviceToHost); //CUDA_SAFE_CALL(cudaMemcpy(HostMem.f,DeviceMem.f,NUM_THREADS*sizeof(Fibers),cudaMemcpyDeviceToHost));
		calculate_reflectance(HostMem.f,up_fluorescence, down_fluorescence, ExPhoton, Fluoro);
		cudaMemcpy(&H_num_of_fluoro, D_num_of_fluoro, sizeof(int), cudaMemcpyDeviceToHost);
	}
	//cout << "#" << index << " Simulation done!\n";
	//cout << *H_num_of_fluoro << "@@" << endl;

	output_fiber(simulation,up_fluorescence, down_fluorescence, ExPhoton, Fluoro, *H_num_of_fluoro);
	cudaFree(D_num_of_fluoro);
	FreeMemStructs(&HostMem,&DeviceMem);
}

void calculate_reflectance(Fibers* f, float *up, float *down, float *ExP, float * F)
{
	for(int i = 0; i < NUM_THREADS; i++)
	{
		// normal configuration
		if(NORMAL)
		{
			ExP[0] += f[i].Exphoton[0];
			F[0] += f[i].Fluoro[0];

			for(int k = 1; k <= NUM_OF_DETECTOR; k++)
			{
				up[k-1] += f[i].up_data[k];
				down[k-1] += f[i].down_data[k];
			}

		}
		// oblique configuration
		else
		{
			ExP[0] += f[i].Exphoton[0];
			F[0] += f[i].Fluoro[0];

			up[0] += f[i].up_data[1];
			up[1] += f[i].up_data[2];
			up[2] += f[i].up_data[3];
			up[3] += f[i].up_data[4] + f[i].up_data[7];
			up[4] += f[i].up_data[5] + f[i].up_data[8];
			up[5] += f[i].up_data[6] + f[i].up_data[9];
			down[0] += f[i].down_data[1];
			down[1] += f[i].down_data[2];
			down[2] += f[i].down_data[3];
			down[3] += f[i].down_data[4] + f[i].down_data[7];
			down[4] += f[i].down_data[5] + f[i].down_data[8];
			down[5] += f[i].down_data[6] + f[i].down_data[9];
		}
	}
}

//Device function to add an unsigned integer to an unsigned long long using CUDA Compute Capability 1.1
__device__ void AtomicAddULL(unsigned long long* address, unsigned int add)
{
	if(atomicAdd((unsigned int*)address,add)+add<add)
		atomicAdd(((unsigned int*)address)+1,1u);
}

__global__ void MCd(MemStruct DeviceMem, unsigned long long seed, int* n)
{
    //Block index
    int bx = blockIdx.x;

    //Thread index
    int tx = threadIdx.x;

    //First element processed by the block
    int begin = NUM_THREADS_PER_BLOCK * bx;

	float s;	//step length

	float w_temp;

	PhotonStruct p = DeviceMem.p[begin+tx];
	//PhotonStruct fluoro = DeviceMem.fluoro[begin+tx];
	Fibers f = DeviceMem.f[begin+tx];

	int new_layer;

	curandState state = DeviceMem.state[begin+tx];
    curand_init(seed, begin+tx, 0, &state);

	//First, make sure the thread (photon) is active
	unsigned int ii = 0;
	if(!DeviceMem.thread_active[begin+tx]) ii = NUMSTEPS_GPU;

	for(;ii<NUMSTEPS_GPU;ii++) //this is the main while loop
	{
		if(layers_dc[p.layer].mutr!=FLT_MAX)
		{
			if (p.fluoro == false)
			s = -__logf(rn_gen(&state))*layers_dc[p.layer].mutr;//sample step length [cm] //HERE AN OPEN_OPEN FUNCTION WOULD BE APPRECIATED
			else
			s = -__logf(rn_gen(&state))*layers_dc[p.layer].mutrE;//sample step length [cm] //HERE AN OPEN_OPEN FUNCTION WOULD BE APPRECIATED
		}
		else
			s = 100.0f;//temporary, say the step in glass is 100 cm.

		//Check for layer transitions and in case, calculate s
		new_layer = p.layer;
		if(p.z+s*p.dz<layers_dc[p.layer].z_min){new_layer--; s = __fdividef(layers_dc[p.layer].z_min-p.z,p.dz);} //Check for upwards reflection/transmission & calculate new s
		if(p.z+s*p.dz>layers_dc[p.layer].z_max){new_layer++; s = __fdividef(layers_dc[p.layer].z_max-p.z,p.dz);} //Check for downward reflection/transmission

		p.x += p.dx*s;
		p.y += p.dy*s;
		p.z += p.dz*s;

		// 20150313
		if(p.z>layers_dc[p.layer].z_max)p.z=layers_dc[p.layer].z_max;
		if(p.z<layers_dc[p.layer].z_min)p.z=layers_dc[p.layer].z_min;

		if(new_layer!=p.layer)
		{
			// set the remaining step length to 0
			s = 0.0f;

			if(Reflect(&p,new_layer,&state)==0u)//Check for reflection
			{
				if (p.fluoro == false)
				{
					if(new_layer == 0)
					{ //Diffuse reflectance
						detect(&p,&f);
						p.weight = 0; // Set the remaining weight to 0, effectively killing the photon
					}
					if(new_layer > *n_layers_dc)
					{	//Transmitted
						p.weight = 0; // Set the remaining weight to 0, effectively killing the photon
					}
				}

				else
				{
					if(new_layer == 0)
					{ //Diffuse reflectance
						detect(&p,&f);
						p.weight = 0; // Set the remaining weight to 0, effectively killing the photon
					}
					if(new_layer > *n_layers_dc)
					{	//Transmitted
						p.weight = 0; // Set the remaining weight to 0, effectively killing the photon
					}
				}
			}
		}

		if(s > 0.0f)
		{
			if (p.fluoro == false)
			{
				if(rn_gen(&state) < layers_dc[p.layer].mua*layers_dc[p.layer].mutr)  // absorption
				{
					// if absorbed by fluorophore
					// if(rn_gen(&state) < layers_dc[p.layer].f_mua/(layers_dc[p.layer].f_mua+layers_dc[p.layer].mua)) //20150303
					if(rn_gen(&state) < layers_dc[p.layer].f_mua/(layers_dc[p.layer].mua))
					{
						// if emitting fluorescence at particular wavelength
						if(rn_gen(&state) < layers_dc[p.layer].quantum * layers_dc[p.layer].emission_p)
						{
							*n = *n + 1;                             // count the number of fluorophore
							fluoro_MC(&p, &state, &f);
						}

						else
							p.weight = 0;
					}

					else
						p.weight = 0;
				}

				Spin(&p,layers_dc[p.layer].g,&state);
			}

			else
			{
				// weighted

				w_temp = layers_dc[p.layer].muaE*layers_dc[p.layer].mutrE*p.weight;
				p.weight -= w_temp;

				// fix weight
				/*
				if(rn_gen(&state) < layers_dc[p.layer].muaE*layers_dc[p.layer].mutrE)  // absorption
				{
					p.weight = 0;
				}
				*/
				Spin(&p,layers_dc[p.layer].gE,&state);
			}
		}

		if(!PhotonSurvive(&p,&state)) //if the photon doesn't survive
		{
			f.Exphoton[0] ++;

			if(atomicAdd(DeviceMem.num_terminated_photons,1u) < (*num_photons_dc-NUM_THREADS))
			{	// Ok to launch another photon
				LaunchPhoton(&p,&state);//Launch a new photon
			}
			else
			{	// No more photons should be launched.
				DeviceMem.thread_active[begin+tx] = 0u; // Set thread to inactive
				ii = NUMSTEPS_GPU;				// Exit main loop
			}
		}

	}//end main for loop!

	__syncthreads();	//necessary?

	//save the state of the MC simulation in global memory before exiting
	DeviceMem.p[begin+tx] = p;	//This one is incoherent!!!
	DeviceMem.f[begin+tx] = f;

}//end MCd

__device__ void fluoro_MC(PhotonStruct* p, curandState *state, Fibers *f)
{
	f->Fluoro[0] ++;

	float theta = 2 * PI * rn_gen(state);
    p->dz = -1 + 2 * rn_gen(state);
    p->dx = sqrt(1-p->dz*p->dz)*cos(theta);
    p->dy = sqrt(1-p->dz*p->dz)*sin(theta);

	p->fluoro = true;
	p->fluoro_layer = p->layer;
}

__device__ void LaunchPhoton(PhotonStruct* p, curandState *state)
{
	float rnd_Azimuth, rnd_direction, rnd_rotated;
	float AzimuthAngle;
	float launchPosition;
	float theta_direction;
	float rotated_angle;
	float uxprime, uyprime, uzprime;
	float angle = -ANGLE * PI / 180;

	rnd_Azimuth    = rn_gen(state);
	rnd_direction  = rn_gen(state);
	rnd_rotated    = rn_gen(state);
	AzimuthAngle   = 2 * PI * rnd_Azimuth;
	rotated_angle  = 2 * PI * rnd_rotated;

	float beam_width = illumination_r;  // 200 um, Gaussian beam profile

	launchPosition = beam_width*sqrt(-log(rn_gen(state))/2.0);

	p->x = launchPosition*cos(AzimuthAngle)/cos(angle);
	p->y = launchPosition*sin(AzimuthAngle);
	p->z = 0.0;

	theta_direction = asin(NAOfSource/n_source)*rnd_direction;
	p->dz = cos(theta_direction);
	p->dx = sin(theta_direction) * cos(rotated_angle);
	p->dy = sin(theta_direction) * sin(rotated_angle);

	uxprime = cos(angle)*p->dx - sin(angle)*p->dz;
	uyprime = sin(theta_direction)*sin(rotated_angle);
	uzprime = sin(angle)*p->dx + cos(angle)*p->dz;

	p->dx = uxprime, p->dy = uyprime, p->dz = uzprime;

	p->layer = 1;
	p->weight = *start_weight_dc; //specular reflection!
	p->fluoro = false;
	p->fluoro_layer = 0;
}

__global__ void LaunchPhoton_Global(MemStruct DeviceMem, unsigned long long seed)
{
	int bx=blockIdx.x;
    int tx=threadIdx.x;

    //First element processed by the block
    int begin=NUM_THREADS_PER_BLOCK*bx;

	PhotonStruct p;

	curandState state = DeviceMem.state[begin+tx];
    curand_init(seed, 0, 0, &state);

	LaunchPhoton(&p,&state);

	//__syncthreads();
	DeviceMem.p[begin+tx]=p;		//incoherent!?
}

/*
__device__ void Spin(PhotonStruct* p, float g, curandState *state)
{
	float theta, cost, sint;	// cosine and sine of the
						// polar deflection angle theta.
	float cosp, sinp;	// cosine and sine of the
						// azimuthal angle psi.
	float temp;
	float tempdir=p->dx;

	//This is more efficient for g!=0 but of course less efficient for g==0
	temp = __fdividef((1.0f-(g)*(g)),(1.0f-(g)+2.0f*(g)*rn_gen(state)));//Should be close close????!!!!!
	cost = __fdividef((1.0f+(g)*(g) - temp*temp),(2.0f*(g)));
	if(g==0.0f)
		cost = 2.0f*rn_gen(state)-1.0f;//Should be close close??!!!!!

	sint = sqrtf(1.0f - cost*cost);

	__sincosf(2.0f*PI*rn_gen(state),&sinp,&cosp);// spin psi [0-2*PI)

	temp = sqrtf(1.0f - p->dz*p->dz);

	if(temp==0.0f) //normal incident.
	{
		p->dx = sint*cosp;
		p->dy = sint*sinp;
		p->dz = copysignf(cost,p->dz*cost);
	}
	else // regular incident.
	{
		p->dx = __fdividef(sint*(p->dx*p->dz*cosp - p->dy*sinp),temp) + p->dx*cost;
		p->dy = __fdividef(sint*(p->dy*p->dz*cosp + tempdir*sinp),temp) + p->dy*cost;
		p->dz = -sint*cosp*temp + p->dz*cost;
	}

	//normalisation seems to be required as we are using floats! Otherwise the small numerical error will accumulate
	temp=rsqrtf(p->dx*p->dx+p->dy*p->dy+p->dz*p->dz);
	p->dx = p->dx*temp;
	p->dy = p->dy*temp;
	p->dz = p->dz*temp;
}// end Spin
*/
__device__ int binarySearch(float *data, float value)
{
    int middle;
	int left = 0, right = 180;
    while (left <= right)
    {
        middle = (right + left) / 2;

        if (data[middle] == value)
            return middle;

        if (data[middle] > value)
            right = middle - 1;
        else
            left = middle + 1;
    }
	if (data[middle] > value)
	    return middle;
	else
		return middle + 1;
}

__device__ void Spin(PhotonStruct* p, float *g, curandState *state)
{
	float theta, cost, sint;	// cosine and sine of the
						// polar deflection angle theta.
	float cosp, sinp;	// cosine and sine of the
						// azimuthal angle psi.
	float temp;
	float tempdir=p->dx;

	float rn = rn_gen(state);
    int sample;

	sample = binarySearch(g,rn);

	theta = sample-1+__fdividef((rn-g[sample-1]),(g[sample]-g[sample-1]));
    theta = __fdividef(theta*PI,180);
	cost = cos(theta);

	sint = sqrtf(1.0f - cost*cost);

	__sincosf(2.0f*PI*rn_gen(state),&sinp,&cosp);// spin psi [0-2*PI)

	temp = sqrtf(1.0f - p->dz*p->dz);

	if(temp==0.0f) //normal incident.
	{
		p->dx = sint*cosp;
		p->dy = sint*sinp;
		p->dz = copysignf(cost,p->dz*cost);
	}
	else // regular incident.
	{
		p->dx = __fdividef(sint*(p->dx*p->dz*cosp - p->dy*sinp),temp) + p->dx*cost;
		p->dy = __fdividef(sint*(p->dy*p->dz*cosp + tempdir*sinp),temp) + p->dy*cost;
		p->dz = -sint*cosp*temp + p->dz*cost;
	}

	//normalisation seems to be required as we are using floats! Otherwise the small numerical error will accumulate
	temp=rsqrtf(p->dx*p->dx+p->dy*p->dy+p->dz*p->dz);
	p->dx = p->dx*temp;
	p->dy = p->dy*temp;
	p->dz = p->dz*temp;
}// end Spin



__device__ unsigned int Reflect(PhotonStruct* p, int new_layer, curandState *state)
{
	//Calculates whether the photon is reflected (returns 1) or not (returns 0)
	// Reflect() will also update the current photon layer (after transmission) and photon direction (both transmission and reflection)

	float n1 = layers_dc[p->layer].n;
	float n2 = layers_dc[new_layer].n;
	float r;
	float cos_angle_i = fabsf(p->dz);

	if(n1==n2)//refraction index matching automatic transmission and no direction change
	{
		p->layer = new_layer;
		return 0u;
	}

	if(n1>n2 && n2*n2<n1*n1*(1-cos_angle_i*cos_angle_i))//total internal reflection, no layer change but z-direction mirroring
	{
		p->dz *= -1.0f;
		return 1u;
	}

	if(cos_angle_i==1.0f)//normal incident
	{
		r = __fdividef((n1-n2),(n1+n2));
		if(rn_gen(state)<=r*r)
		{
			//reflection, no layer change but z-direction mirroring
			p->dz *= -1.0f;
			return 1u;
		}
		else
		{	//transmission, no direction change but layer change
			p->layer = new_layer;
			return 0u;
		}
	}

	//gives almost exactly the same results as the old MCML way of doing the calculation but does it slightly faster
	// save a few multiplications, calculate cos_angle_i^2;
	float e = __fdividef(n1*n1,n2*n2)*(1.0f-cos_angle_i*cos_angle_i); //e is the sin square of the transmission angle
	r=2*sqrtf((1.0f-cos_angle_i*cos_angle_i)*(1.0f-e)*e*cos_angle_i*cos_angle_i);//use r as a temporary variable
	e=e+(cos_angle_i*cos_angle_i)*(1.0f-2.0f*e);//Update the value of e
	r = e*__fdividef((1.0f-e-r),((1.0f-e+r)*(e+r)));//Calculate r

	if(rn_gen(state)<=r)
	{
		// Reflection, mirror z-direction!
		p->dz *= -1.0f;
		return 1u;
	}
	else
	{
		// Transmission, update layer and direction
		r = __fdividef(n1,n2);
		e = r*r*(1.0f-cos_angle_i*cos_angle_i); //e is the sin square of the transmission angle
		p->dx *= r;
		p->dy *= r;
		p->dz = copysignf(sqrtf(1-e) ,p->dz);
		p->layer = new_layer;
		return 0u;
	}

}

__device__ unsigned int PhotonSurvive(PhotonStruct* p, curandState *state)
{
	//Calculate wether the photon survives (returns 1) or dies (returns 0)

	if(p->weight>WEIGHTI) return 1u; // No roulette needed
	if(p->weight==0) return 0u;	// Photon has exited slab, i.e. kill the photon

	if(rn_gen(state) < CHANCE)
	{
		//p->weight = __float2uint_rn(__fdividef((float)p->weight,CHANCE));
		p->weight = __fdividef((float)p->weight,CHANCE);
		return 1u;
	}
	return 0u;
}

__device__ void detect(PhotonStruct* p, Fibers* f)
{
	float angle = ANGLE*PI/180;
	float critical = asin(f->NA[1]/ n_detector);
    float uz_rotated=(p->dx*sin(angle))+(p->dz*cos(angle));
	float uz_angle = acos(fabs(uz_rotated));
	float distance;

	if(uz_angle <= critical)  // successfully detected
	{
		if(NORMAL)   // normal configuration
		{
			if(p->fluoro == true)
			{
				// ISS circle
				/*
				for(int i = 1; i <= 3 ; i++)
				{
					if(pow((p->x-f->position[i])*cos(angle),2) + pow(p->y,2) <= f->radius[i]*f->radius[i])
					{
						if(p->fluoro_layer==1)
							f->up_data[i] += p->weight;
						else
							f->down_data[i] += p->weight;
					}
				}

				for(int i = 4; i <= 6 ; i++)
				{
					if(pow((p->y-f->position[i]),2) + pow(p->x*cos(angle),2) <= f->radius[i]*f->radius[i])
					{
						if(p->fluoro_layer==1)
							f->up_data[i] += p->weight;
						else
							f->down_data[i] += p->weight;
					}

				}
				*/

				// ISS annular
				distance = sqrt(p->x * p->x + p->y * p->y);

				for(int i = 1; i <= 6 ; i++)
				{
					if((distance>=(f->position[i]-f->radius[i])) && (distance<=(f->position[i]+f->radius[i])))
					{
						float temp;
						temp = (distance*distance + f->position[i]*f->position[i] - f->radius[i]*f->radius[i])/(2*distance*f->position[i]);
						// check for rounding error!
						if(temp > 1.0f)
							temp = 1.0f;

						if(p->fluoro_layer==1)
							f->up_data[i] += p->weight * acos(temp) * RPI;
						else
							f->down_data[i] += p->weight * acos(temp) * RPI;

					}
				}

			}

		}
		// oblique configuration
		else
		{

		}
	}
    return;
}

int InitDCMem(SimulationStruct* sim)
{
	// Copy num_photons_dc to constant device memory
	cudaMemcpyToSymbol(n_layers_dc,&(sim->n_layers),sizeof(unsigned int));

	// Copy start_weight_dc to constant device memory
	cudaMemcpyToSymbol(start_weight_dc,&(sim->start_weight),sizeof(float));

	// Copy layer data to constant device memory
	cudaMemcpyToSymbol(layers_dc,sim->layers,(sim->n_layers+2)*sizeof(LayerStruct));

	// Copy num_photons_dc to constant device memory
	cudaMemcpyToSymbol(num_photons_dc,&(sim->number_of_photons),sizeof(unsigned int));

	return 0;
}

int InitMemStructs(MemStruct* HostMem, MemStruct* DeviceMem, SimulationStruct* sim)
{
	// Allocate p on the device!!
	cudaMalloc((void**)&DeviceMem->p,NUM_THREADS*sizeof(PhotonStruct));

	// Allocate fluoro on the device!!
	// cudaMalloc((void**)&DeviceMem->fluoro,NUM_THREADS*sizeof(PhotonStruct));

	// Allocate thread_active on the device and host
	HostMem->thread_active = (unsigned int*) malloc(NUM_THREADS*sizeof(unsigned int));
	if(HostMem->thread_active==NULL){printf("Error allocating HostMem->thread_active"); exit (1);}
	for(int i=0;i<NUM_THREADS;i++)HostMem->thread_active[i]=1u;

	cudaMalloc((void**)&DeviceMem->thread_active,NUM_THREADS*sizeof(unsigned int));
	cudaMemcpy(DeviceMem->thread_active,HostMem->thread_active,NUM_THREADS*sizeof(unsigned int),cudaMemcpyHostToDevice);

	//Allocate num_launched_photons on the device and host
	HostMem->num_terminated_photons = (unsigned int*) malloc(sizeof(unsigned int));
	if(HostMem->num_terminated_photons==NULL){printf("Error allocating HostMem->num_terminated_photons"); exit (1);}
	*HostMem->num_terminated_photons=0;

	cudaMalloc((void**)&DeviceMem->num_terminated_photons,sizeof(unsigned int));
	cudaMemcpy(DeviceMem->num_terminated_photons,HostMem->num_terminated_photons,sizeof(unsigned int),cudaMemcpyHostToDevice);

	//Allocate and initialize fiber f on the device and host
	HostMem->f = (Fibers*) malloc(NUM_THREADS*sizeof(Fibers));
	cudaMalloc((void**)&DeviceMem->f,NUM_THREADS*sizeof(Fibers));
	fiber_initialization(HostMem->f);
	cudaMemcpy(DeviceMem->f,HostMem->f,NUM_THREADS*sizeof(Fibers),cudaMemcpyHostToDevice);

	//Allocate states on the device and host
	cudaMalloc((void**)&DeviceMem->state,NUM_THREADS*sizeof(curandState));


	return 1;
}

void FreeMemStructs(MemStruct* HostMem, MemStruct* DeviceMem)
{
	free(HostMem->thread_active);
	free(HostMem->num_terminated_photons);
	free(HostMem->f);

	cudaFree(DeviceMem->thread_active);
	cudaFree(DeviceMem->num_terminated_photons);
	cudaFree(DeviceMem->f);
	cudaFree(DeviceMem->state);
}
