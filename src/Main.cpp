// This file is part of Swept2D
// Copyright (C) 2015 Qiqi Wang, qiqi@mit.edu AND Maitham Alhubail, hubailmm@mit.edu
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) an later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT An WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "Main.h"

int main(int argc,char *argv[])
{
	const int n       = atoi(argv[1]);
	const int xNodes  = atoi(argv[2]);
	const int yNodes  = atoi(argv[3]);
	const int remoteC = atoi(argv[4]);
	cycles            = 10000000;//atoi(argv[4]);
	fileOutputEvery   = -1;
	conservationCheckEvery = -1;
	totalGlobals = 0;

	initMPI(argc,argv,xNodes,yNodes,true);

	Substeps2D::setInitParameters(xNodes*n,yNodes*n);
	if(totalGlobals > 0)
	{
		globals = new double[totalGlobals];
		Substeps2D::setGlobalVariables(xNodes*n,yNodes*n);
	}
	
	for(int i=0;i<substeps;i++)
	{
		localCount[i] = 0;
		stepCount[i]  = 0;
	}
	double ellapsed;
	SweptDiscretization2D swept2d(n,substeps,dataPointSize,constants,totalConservedQuantities,remoteConstants);
	swept2d.setInitCondition(Substeps2D::init);
	swept2d.setOutputDirectory(outputDirectory);
	swept2d.setFileOutput(fileOutputEvery);
	swept2d.setConservationCheck(conservationCheckEvery);
	swept2d.setCheckForRemoteCommands(remoteC);
	swept2d.setOutputLength(4);

	ClassicDiscretization2D classic2d(n,substeps,dataPointSize,constants,totalConservedQuantities);
	classic2d.setInitCondition(Substeps2D::init);
	//ellapsed = classic2d.calculate(cycles*n);	
	ellapsed = classic2d.calculate(1);	
	
	//printf("\nEnergy:\n");
	//classic2d.printFoundation();

	MPI_Allreduce(&localCount,&stepCount,substeps,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
	
	//if(pg.rank == 0)
	{
			//printf("Classic: Total Substeps: %d - us per substep: %f\n",cycles*n,ellapsed*1000000/(cycles*n));
			//for(int i=0;i<substeps;i++)
			//	printf("Substep %d Function Calls: %ld\n",i+1,stepCount[i]);
			
	}
	
	for(int i=0;i<substeps;i++)
	{
		localCount[i] = 0;
		stepCount[i]  = 0;
	}
	
	ellapsed = swept2d.calculate(cycles);	
	
	//printf("\nEnergy:\n");
	//swept2d.printFoundation();

	MPI_Allreduce(&localCount,&stepCount,substeps,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);	
	//if(pg.rank == 0)
	{
			printf("Swept: Total Substeps: %d - us per substep: %f\n",cycles*n,ellapsed*1000000/(cycles*n));
			for(int i=0;i<substeps;i++)
				printf("Substep %d Function Calls: %ld\n",i+1,stepCount[i]);			
	}
	
	finMPI();
	return 0;
}
