#include "header.h"

void fiber_initialization(Fibers* f)
{
	for(int i = 0; i < NUM_THREADS; i++)
	{
	    // initialization to zeros for data deposition 
		f[i].Exphoton[0] = 0;
		f[i].Fluoro[0] = 0;
		
		for(int j = 0; j <= NUM_OF_DETECTOR; j++)
		{
			f[i].up_data[j] = 0;
			f[i].down_data[j] = 0;

		}
		// normal configuration
		if(NORMAL)   
		{	
			// source fiber	
			f[i].radius[0]   = illumination_r;          		
			f[i].NA[0]       = NAOfSource;				
			f[i].angle[0]    = ANGLE*PI/180;
			f[i].position[0] = 0.0;			
			//first fiber,  SDS = 0.026 cm
			f[i].radius[1]   = collect_r;             
			f[i].NA[1]       = NAOfDetector;				
			f[i].position[1] = 0.022;                    
			f[i].angle[1]    = ANGLE*PI/180;
			//second fiber, SDS = 0.054 cm(before 20151029)
			f[i].radius[2]   = collect_r;          
			f[i].NA[1]       = NAOfDetector;				
			f[i].position[2] = 0.041;                    
			f[i].angle[2]    = ANGLE*PI/180;		
			//third fiber, SDS = 0.078 cm(before 20151103)
			f[i].radius[3]   = collect_r;             
			f[i].NA[3]       = NAOfDetector;				
			f[i].position[3] = 0.061;                    
			f[i].angle[3]    = ANGLE*PI/180;
			//fourth fiber, SDS = 0.049 cm
			f[i].radius[4]   = collect_r;             
			f[i].NA[4]       = NAOfDetector;				
			f[i].position[4] = 0.0215;                    
			f[i].angle[4]    = ANGLE*PI/180;
			//fourth fiber, SDS = 0.076 cm
			f[i].radius[5]   = collect_r;             
			f[i].NA[5]       = NAOfDetector;				
			f[i].position[5] = 0.045;                    
			f[i].angle[5]    = ANGLE*PI/180;
			//fourth fiber, SDS = 0.024 cm
			f[i].radius[6]   = collect_r;             
			f[i].NA[6]       = NAOfDetector;				
			f[i].position[6] = 0.073;                    
			f[i].angle[6]    = ANGLE*PI/180;
		}
		// oblique configuration
		else
		{
			for(int j = 1; j <= 3; j++)
			{
				f[i].radius[j]   = collect_r;	
				f[i].NA[j]       = NAOfDetector;		
				f[i].position[j] = 0.032*j;	
				f[i].angle[j]    = ANGLE*PI/180;
			}
		

			for(int j = 4; j <= 6; j++)
			{
				f[i].radius[j]   = collect_r;	
				f[i].NA[j]       = NAOfDetector;		
				f[i].position[j] = 0.022*(j-3);	
				f[i].angle[j]    = ANGLE*PI/180;
			}

			for(int j = 7; j <= 9; j++)
			{
				f[i].radius[j]   = collect_r;	
				f[i].NA[j]       = NAOfDetector;		
				f[i].position[j] = -0.022*(j-6);	
				f[i].angle[j]    = ANGLE*PI/180;
			}
		}
	}
}