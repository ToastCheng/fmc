#include "header.h"
#include <cstdio>
#include <float.h>


void output_fiber(SimulationStruct* sim, float *up, float* down, float *Exp, float *F, int n)
{

	ofstream myfile, myfile_up, myfile_bt;

	char *sim_, *sim_up, *sim_bt;
	sim_ = (char*)malloc(30);
	sim_up = (char*)malloc(30);
	sim_bt = (char*)malloc(30);
	sim_[0] ='\0';
	sim_up[0] ='\0';
	sim_bt[0] ='\0';
	strcat(sim_,NUM_SPECTRUM);
	strcat(sim_,"/simulation");
	strcat(sim_,NUM_SPECTRUM);
	strcat(sim_,".txt");
	strcat(sim_up,NUM_SPECTRUM);
	strcat(sim_up,"/simulation_up");
	strcat(sim_up,NUM_SPECTRUM);
	strcat(sim_up,".txt");
	strcat(sim_bt,NUM_SPECTRUM);
	strcat(sim_bt,"/simulation_bt");
	strcat(sim_bt,NUM_SPECTRUM);
	strcat(sim_bt,".txt");
	myfile.open (sim_,ios::app);
	myfile_up.open(sim_up, ios::app);
	myfile_bt.open(sim_bt, ios::app);
	if(!myfile){
		 printf("failed to open %s\n",sim_);
		 exit(1);
	 }
	if(!myfile_up){
		 printf("failed to open %s\n",sim_);
		 exit(1);
	 }
	if(!myfile_bt){
		 printf("failed to open %s\n",sim_);
		 exit(1);
	 }

	double tmp;

	float scale1 = *Exp;    // #exitation photo

	// normal
	if(NORMAL == 1)
	{

		//fout2 << *Exp << " " << *F << endl;

		for(int i = 0; i < NUM_OF_DETECTOR; i++)
		{
			tmp = (double)(up[i]+down[i])/scale1;

			//fout << (double)(up[i]/scale1)/tmp << "\t" << (double)(down[i]/scale1)/tmp << "\t";
			//fout2 << (double)up[i]/scale1 << "\t" << (double)down[i]/scale1 << "\t";

			if (i==NUM_OF_DETECTOR-1){
				myfile << (double)(up[i] + down[i]) / scale1;
				myfile_up << (double)(up[i] / scale1);
				myfile_bt << (double)(down[i] / scale1);
			}
			else{
				myfile << (double)(up[i] + down[i]) / scale1 << "\t";
				myfile_up << (double)(up[i] / scale1) << "\t";
				myfile_bt << (double)(down[i] / scale1) << "\t";
			}
		}
	}
	// oblique
	else
	{
		for(int i = 0; i < 6; i++)
		{
			tmp = (double)(up[i]+down[i])/scale1;
			//fout << (double)(up[i]/scale1)/tmp << " " << (double)(down[i]/scale1)/tmp << " ";
			//fout2 << (double)up[i]/scale1 << " " << (double)down[i]/scale1 << " ";
			if (i == NUM_OF_DETECTOR - 1){ myfile << (double)(up[i] + down[i]) / scale1; }
			else{ myfile << (double)(up[i] + down[i]) / scale1 << "\t"; }
		}
	}

	//fout << endl;
	//fout2 << endl;
	myfile << endl;
	myfile_up << endl;
	myfile_bt << endl;

	myfile.close();
	myfile_up.close();
	myfile_bt.close();
	//fout.close();
	//fout2.close();
}

void read_anglepattern(float* gcumf, int nn)
{
	float gall[181] = {};
	for (int k = 0; k<181; k++)
		gcumf[k] = 0;
	float gprob[181] = {};

	// can be modified
	//char *address = "C:\\Users\\user\\Desktop\\GPU\\Anglepattern\\0.5\\";
	char filename[100];

	ifstream infile;
	switch (nn)
	{
		case 1:
			sprintf(filename, "input/tissue_anglepattern.txt"); //sprintf(filename,"%s0.5_%dnm.txt",address,index);  // can be modified
			//printf("Read tissue_anglepattern.txt\n");
			break;
		case 2:
			sprintf(filename, "input/tissue_anglepattern2.txt"); //sprintf(filename,"%s0.5_%dnm.txt",address,index);  // can be modified
			//printf("Read tissue_anglepattern2.txt\n");
			break;
	}
	infile.open(filename);
	if(!infile) printf("can't open anglepattern.txt!\n");
	float sum = 0;

	for (int i = 0; i < 181; i++)
	{
		infile >> gall[i];
		sum += gall[i] * sin(i*PI / 180);
	}

	for (int i = 0; i < 181; i++)
		gprob[i] = gall[i] * sin(i*PI / 180) / sum;

	for (int i = 0; i < 181; i++)
	{
		if (i == 0)
			gcumf[i] = gprob[i];
		else
			gcumf[i] += gcumf[i - 1] + gprob[i];
	}

	//printf("%e\n",gcumf[50]);

	infile.close();
}

//void read_anglepattern(int index, float* gcumf)
//{
//	float gall[181]={};
//	for(int k=0; k<181;k++)
//		gcumf[k] = 0;
//	float gprob[181]={};
//
//	// can be modified
//	//char *address = "C:\\Users\\user\\Desktop\\GPU\\Anglepattern\\0.5\\";
//	// complete address for the file
//	char filename[100];
//
//	ifstream infile;
//	sprintf(filename, "0.5_%dnm.txt", index); //sprintf(filename,"%s0.5_%dnm.txt",address,index);  // can be modified
//	infile.open(filename);
//
//	float sum = 0;
//
//	for(int i = 0; i < 181; i++)
//	{
//		infile >> gall[i];
//		sum += gall[i]*sin(i*PI/180);
//	}
//
//	for(int i = 0; i < 181; i++)
//		gprob[i] = gall[i]*sin(i*PI/180)/sum;
//
//	for(int i = 0; i < 181; i++)
//	{
//		if(i == 0)
//			gcumf[i] = gprob[i];
//		else
//			gcumf[i] += gcumf[i-1] + gprob[i];
//	}
//
//	//printf("%e\n",gcumf[50]);
//
//	infile.close();
//}

int read_data(SimulationStruct** simulations)
{
	// parameters to be modified
	unsigned long number_of_photons = NUMBER_PHOTONS;
	const int n_simulations = NUM_OF_SIMULATION;
	// double layer, default value = 2
	int n_layers = 2;
	// refractive index of outer medium	// water:1.33 tissue:1.58
	float medium_n = 1.60;
	// refractive index of tissue
	float tissue_n = 1.33;

	float start_weight;
	float upper_thickness;

	// read the file
	// anisotropy
	// directory address, maybe needs to be modified
	// can be modified
	//int excitation[1] = {365};
	// can be modified
	//int emmition[n_simulations] = {400,410,420,430,440,450,460,470,480,490,500,510,520,530,540,550,560,570,580,590,600,610,620,630,640,650};

	fstream myfile;
	char *para;
	para = (char*)malloc(30);
	para[0] = '\0';
	strcat(para,"DRSresult/");
	strcat(para,NUM_SPECTRUM);
	strcat(para,".txt");
	myfile.open (para,ios::in);
	if(!myfile) printf("failed to open %s!\n",para);
	float up_mua, up_mus, down_mua, down_mus;                      // excitation
	float up_mua_E[n_simulations],up_mus_E[n_simulations],         // emission
		  down_mua_E[n_simulations],down_mus_E[n_simulations];
	float up_quantum, down_quantum;                                // quantum yield
	float up_f_mua,   down_f_mua;                                  // mua of fluoro at excitation
	float up_emission_p[n_simulations], down_emission_p[n_simulations];
	float wavelength[n_simulations+1];



	myfile >> upper_thickness;     // input the thickness of upper layer
	myfile >> wavelength[0] >> up_mua >> up_mus >> down_mua >> down_mus;   // excitation parameters
	up_f_mua = up_mua;
	down_f_mua = down_mua;
	for(int i = 0; i < n_simulations; i++)   // emission parameters
	{
	    myfile >> wavelength[i+1] >> up_mua_E[i] >> up_mus_E[i]
		       >> down_mua_E[i] >> down_mus_E[i];
	}

	up_quantum = 1;
	down_quantum = 1;
	for(int i = 0; i < n_simulations; i++){   // emission probability at particular wavelength
	    up_emission_p[i] = 1;
		down_emission_p[i] = 1;
	}
	myfile.close();

	// Allocate memory for the SimulationStruct array
	*simulations = (SimulationStruct*) malloc(sizeof(SimulationStruct)*n_simulations);
	if(*simulations == NULL){perror("Failed to malloc simulations.\n");return 0;}

	for(int i = 0; i < n_simulations; i++)
	{
		(*simulations)[i].number_of_photons=number_of_photons;
		(*simulations)[i].n_layers = n_layers;

		// Allocate memory for the layers (including one for the upper and one for the lower)
		(*simulations)[i].layers = (LayerStruct*) malloc(sizeof(LayerStruct)*(n_layers+2));
		if((*simulations)[i].layers == NULL){perror("Failed to malloc layers.\n");return 0;}//{printf("Failed to malloc simulations.\n");return 0;}

		// Set upper refractive index (medium)
		(*simulations)[i].layers[0].n = medium_n;

		// Set the parameters of tissue (upper layer)
		(*simulations)[i].layers[1].n     = tissue_n;
		(*simulations)[i].layers[1].mua   = up_mua;
		(*simulations)[i].layers[1].muaE  = up_mua_E[i];
		// Set angle pattern
		read_anglepattern((*simulations)[i].layers[1].g, 1);   //read_anglepattern(excitation[0], (*simulations)[i].layers[1].g); //(*simulations)[i].layers[1].g     = read_anglepattern(excitation[0]);
		read_anglepattern((*simulations)[i].layers[1].gE, 1);  //read_anglepattern(emmition[i], (*simulations)[i].layers[1].gE); //(*simulations)[i].layers[1].gE    = read_anglepattern(emmition[i]);
		// Set other parameters
		(*simulations)[i].layers[1].z_min = 0;
		(*simulations)[i].layers[1].z_max = upper_thickness;
		(*simulations)[i].layers[1].mutr  = 1.0f/(up_mua+up_mus);
		(*simulations)[i].layers[1].mutrE = 1.0f/(up_mua_E[i]+up_mus_E[i]);
		(*simulations)[i].layers[1].f_mua = up_f_mua;
		(*simulations)[i].layers[1].quantum = up_quantum;
		(*simulations)[i].layers[1].emission_p = up_emission_p[i];

		// Set the parameters of tissue (lower layer)
		(*simulations)[i].layers[2].n     = tissue_n;
		(*simulations)[i].layers[2].mua   = down_mua;
		(*simulations)[i].layers[2].muaE  = down_mua_E[i];
		// Set angle pattern
		read_anglepattern((*simulations)[i].layers[2].g, 2);   //read_anglepattern(excitation[0], (*simulations)[i].layers[2].g);//(*simulations)[i].layers[2].g     = read_anglepattern(excitation[0]);
		read_anglepattern((*simulations)[i].layers[2].gE, 2);  //read_anglepattern(emmition[i], (*simulations)[i].layers[2].gE);//(*simulations)[i].layers[2].gE    = read_anglepattern(emmition[i]);
		// Set other parameters
		(*simulations)[i].layers[2].z_min = upper_thickness;
		//////////////////////////////////////////////////////////////////////////////
		//TOAST MODIFY 1 -> FLT_MAX
		(*simulations)[i].layers[2].z_max = FLT_MAX;				// set as infinity
		(*simulations)[i].layers[2].mutr  = 1.0f/(down_mua+down_mus);
		(*simulations)[i].layers[2].mutrE = 1.0f/(down_mua_E[i]+down_mus_E[i]);
		(*simulations)[i].layers[2].f_mua = down_f_mua;
		(*simulations)[i].layers[2].quantum = down_quantum;
		(*simulations)[i].layers[2].emission_p = down_emission_p[i];

		// Set lower refractive index (medium)
		(*simulations)[i].layers[n_layers+1].n = tissue_n;		// use "tissue_n" for no reflectance(assume that semi-infinity); use "medium_n" for layers

		//calculate start_weight
		float n1=n_source;
		float n2=(*simulations)[i].layers[1].n;
		float r = (n1-n2)/(n1+n2);
		r = r*r;
		start_weight = 1.0 * (1.0-r);
		//printf("Start weight=%e\n",start_weight);
		(*simulations)[i].start_weight=start_weight;

	}
	//system("pause");

	return n_simulations;
}
