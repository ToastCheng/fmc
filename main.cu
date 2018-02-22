#include "header.h"

void FreeSimulationStruct(SimulationStruct* sim, int n_simulations);
int read_data(SimulationStruct** simulations);
void DoOneSimulation(SimulationStruct* simulation, int index);
char* NUM_SPECTRUM;
int main(int argc,char* argv[])
{

	NUM_SPECTRUM = (char*)malloc(50);
	strcat(NUM_SPECTRUM,argv[1]);
	printf("#%s spectrum simulating\n",NUM_SPECTRUM);
	SimulationStruct* simulations;
	int n_simulations;
	unsigned long long seed = (unsigned long long) time(NULL);// Default, use time(NULL) as seed
	n_simulations = read_data(&simulations); // read the input file

	if(n_simulations == 0)
	{
		printf("Something wrong with read_simulation_data!\n");
		return 1;
	}

	clock_t time1,time2;

	// Start the clock
    time1=clock();

	//perform all the simulations
	for(int i = 0; i < n_simulations; i++)
	{
		// Run a simulation
		//printf("simulation:%d/%d\n",i+1,n_simulations);
		DoOneSimulation(&simulations[i],i);
	}

	time2=clock();
	printf("Simulation time: %.2f sec\n",(double)(time2-time1)/CLOCKS_PER_SEC);

	FreeSimulationStruct(simulations, n_simulations);

    //system("PAUSE");
	return 0;
}

void FreeSimulationStruct(SimulationStruct* sim, int n_simulations)
{
	for(int i = 0;i < n_simulations; i++) free(sim[i].layers);
	free(sim);
}
