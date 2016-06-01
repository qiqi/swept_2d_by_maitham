// This file is part of Swept2D
// Copyright (C) 2015 Qiqi Wang, qiqi@mit.edu AND Maitham Alhubail, hubailmm@mit.edu
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) an later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT Any WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef H_SUBSTEPS2D
#define H_SUBSTEPS2D

#include <math.h>
#include "Swept2DGlobals.h"

extern void printHex(double a);

class Substeps2D
{

public:

	static inline void executeStepFnc(int executeFnc,PointStruct2D *point)
	{
		const double DISS_COEFF = 0.07;
		switch(executeFnc)
		{
			case(0):
			{
				int inputShift   = (executeFnc + 0) * 4;
				int outputShift  = (executeFnc + 1) * 4;
				int forwardIndex = 4;
				for (int i = 0; i < forwardIndex; i++) 
				{
					point->output[i] = point->C_input[i];
				}
				
				// compute the dc * rho * Laplace(u)
				double rc = point->C_input[0 + inputShift];
				double rn = point->N_input[0 + inputShift];
				double rs = point->S_input[0 + inputShift];
				double re = point->E_input[0 + inputShift];
				double rw = point->W_input[0 + inputShift];
				double laplaceRho = (re*re + rw*rw + rn*rn + rs*rs) * .25 - (rc*rc);
				point->ConstantsOutput[2] = laplaceRho;

				for (int i = 0; i < 2; ++i) 
				{
					double uc = point->C_input[i + 1 + inputShift] / rc;
					double un = point->N_input[i + 1 + inputShift] / rn;
					double us = point->S_input[i + 1 + inputShift] / rs;
					double ue = point->E_input[i + 1 + inputShift] / re;
					double uw = point->W_input[i + 1 + inputShift] / rw;
					double laplaceU = (ue + uw + un + us) * .25 - uc;					
					point->ConstantsOutput[i + 1 + 2] = laplaceU;
				}

				double pc = point->C_input[3 + inputShift];
				double pn = point->N_input[3 + inputShift];
				double ps = point->S_input[3 + inputShift];
				double pe = point->E_input[3 + inputShift];
				double pw = point->W_input[3 + inputShift];
				double laplaceP = (pe + pw + pn + ps) * .25 - pc;
				point->ConstantsOutput[5] = laplaceP;

				localCount[executeFnc]++;
			}
			break;
			case(1):
			{
				int inputShift  = (executeFnc - 1) * 4;
				int outputShift = (executeFnc + 0) * 4;
				int forwardIndex = 4;
				for (int i = 0; i < forwardIndex; i++) 
				{
					point->output[i] = point->C_input[i];
				}

				double rhoV_n = point->N_input[0 + inputShift] * point->N_input[2 + inputShift];
				double rhoV_s = point->S_input[0 + inputShift] * point->S_input[2 + inputShift];
				double rhoU_e = point->E_input[0 + inputShift] * point->E_input[1 + inputShift];
				double rhoU_w = point->W_input[0 + inputShift] * point->W_input[1 + inputShift];
				
				double mass = (rhoU_e - rhoU_w) / (2 * dx) + (rhoV_s - rhoV_n) / (2 * dy);
				
				point->C_conserved[0] = mass;
				
				//////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////////////////////////////

				double rhoU_c = point->C_input[0 + inputShift] * point->C_input[1 + inputShift];
				double rhoV_c = point->C_input[0 + inputShift] * point->C_input[2 + inputShift];

				double un = point->N_input[1 + inputShift] / point->N_input[0 + inputShift];
				double us = point->S_input[1 + inputShift] / point->S_input[0 + inputShift];
				double ue = point->E_input[1 + inputShift] / point->E_input[0 + inputShift];
				double uw = point->W_input[1 + inputShift] / point->W_input[0 + inputShift];

				double rhoUV_n = point->N_input[1 + inputShift] * point->N_input[2 + inputShift];
				double rhoUV_s = point->S_input[1 + inputShift] * point->S_input[2 + inputShift];
				double rhoUU_e = point->E_input[1 + inputShift] * point->E_input[1 + inputShift];
				double rhoUU_w = point->W_input[1 + inputShift] * point->W_input[1 + inputShift];

				double pe = point->E_input[3 + inputShift];
				double pw = point->W_input[3 + inputShift];

				
				double momentum_x = ((rhoUU_e - rhoUU_w) / (2 * dx) + rhoU_c * (ue - uw) / (2 * dx)) / 2.0
					              + ((rhoUV_s - rhoUV_n) / (2 * dy) + rhoV_c * (us - un) / (2 * dy)) / 2.0
								  + (pe - pw) / (2 * dx);
							
				point->C_conserved[1] = momentum_x;
				
				double vn = point->N_input[2 + inputShift] / point->N_input[0 + inputShift];
				double vs = point->S_input[2 + inputShift] / point->S_input[0 + inputShift];
				double ve = point->E_input[2 + inputShift] / point->E_input[0 + inputShift];
				double vw = point->W_input[2 + inputShift] / point->W_input[0 + inputShift];

				double rhoVV_n = point->N_input[2 + inputShift] * point->N_input[2 + inputShift];
				double rhoVV_s = point->S_input[2 + inputShift] * point->S_input[2 + inputShift];
				double rhoUV_e = point->E_input[1 + inputShift] * point->E_input[2 + inputShift];
				double rhoUV_w = point->W_input[1 + inputShift] * point->W_input[2 + inputShift];

				double pn = point->N_input[3 + inputShift];
				double ps = point->S_input[3 + inputShift];

				
				double momentum_y = ((rhoUV_e - rhoUV_w) / (2 * dx) + rhoU_c * (ve - vw) / (2 * dx)) / 2.0
					              + ((rhoVV_s - rhoVV_n) / (2 * dy) + rhoV_c * (vs - vn) / (2 * dy)) / 2.0
								  + (ps - pn) / (2 * dy);
				
				point->C_conserved[2] = momentum_y;
				
				//////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////////////////////////////
				
				
				double gamma = 1.4;
				double uc = point->C_input[1 + inputShift] / point->C_input[0 + inputShift];
				double vc = point->C_input[2 + inputShift] / point->C_input[0 + inputShift];
				
				double diffpu = (pe*ue - pw*uw) / (2 * dx);
				double diffpv = (ps*vs - pn*vn) / (2 * dy);
				double udiffp = uc * ((pe - pw) / (2 * dx));
				double vdiffp = vc * ((ps - pn) / (2 * dy));

				double energy = gamma     * (diffpu + diffpv)
					          - (gamma-1) * (udiffp + vdiffp);

				point->C_conserved[3] = energy;
				
				//////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////////////////////////////
				// Compute dissipation
				
				// compute the Laplacian of dc * rho * laplace(u)
				double laplaceDiss[4];

				double diss_c = DISS_COEFF * (1 * 1) * point->C_constants[2];
				double diss_n = DISS_COEFF * (1 * 1) * point->N_constants[2];
				double diss_s = DISS_COEFF * (1 * 1) * point->S_constants[2];
				double diss_e = DISS_COEFF * (1 * 1) * point->E_constants[2];
				double diss_w = DISS_COEFF * (1 * 1) * point->W_constants[2];
				laplaceDiss[0] = (diss_w + diss_e + diss_n + diss_s) * .25 - diss_c;

				for (int i = 0; i < 2; ++i) 
				{
					diss_c = DISS_COEFF * (point->C_input[0 + inputShift] * point->C_input[0 + inputShift]) * point->C_constants[i + 1 + 2];
					diss_n = DISS_COEFF * (point->N_input[0 + inputShift] * point->N_input[0 + inputShift]) * point->N_constants[i + 1 + 2];
					diss_s = DISS_COEFF * (point->S_input[0 + inputShift] * point->S_input[0 + inputShift]) * point->S_constants[i + 1 + 2];
					diss_e = DISS_COEFF * (point->E_input[0 + inputShift] * point->E_input[0 + inputShift]) * point->E_constants[i + 1 + 2];
					diss_w = DISS_COEFF * (point->W_input[0 + inputShift] * point->W_input[0 + inputShift]) * point->W_constants[i + 1 + 2];

					laplaceDiss[i+1] = (diss_w + diss_e + diss_n + diss_s) * .25 - diss_c;
				}

				diss_c = DISS_COEFF * (1 * 1) * point->C_constants[5];
				diss_n = DISS_COEFF * (1 * 1) * point->N_constants[5];
				diss_s = DISS_COEFF * (1 * 1) * point->S_constants[5];
				diss_e = DISS_COEFF * (1 * 1) * point->E_constants[5];
				diss_w = DISS_COEFF * (1 * 1) * point->W_constants[5];
				laplaceDiss[3] = (diss_w + diss_e + diss_n + diss_s) * .25 - diss_c;
				

				const double R_air = 287., T0 = 300.;
				const double c0 = sqrt(gamma * R_air * T0);
				double dissipation_r = laplaceDiss[0] * c0 / dx;
				double dissipation_x = laplaceDiss[1] * c0 / dx;
				double dissipation_y = laplaceDiss[2] * c0 / dy;
				double dissipation_p = laplaceDiss[3] * c0 / dx;
				
				//dissipation_r = 0;
				//dissipation_p = 0;
				mass       += dissipation_r;
				momentum_x += dissipation_x;
				momentum_y += dissipation_y;
				energy     += dissipation_p - (gamma - 1) * (uc * dissipation_x + vc * dissipation_y);
				
				double rhs_w[4];
				double rc = point->C_input[0 + inputShift];
				rhs_w[0] = 0.5 * mass / rc;
				rhs_w[1] = momentum_x / rc;
				rhs_w[2] = momentum_y / rc;
				rhs_w[3] = energy;

				//rhs_w[i_obstacle,j_obstacle,1:3] += 0.1*c0 * w[i_obstacle, j_obstacle, 1:3]
				// rhs_w[1] += 0.5 * c0 * point->C_input[1 + inputShift] * point->C_constants[0];	// changed			
				// rhs_w[2] += 0.5 * c0 * point->C_input[2 + inputShift] * point->C_constants[0];  //changed
				
				//rhs_w += 0.1*c0 * (w - w0) * dc[:,:,newaxis]
				//for (int i = 0; i < 4; i++) 
				//{
				//	rhs_w[i] += 0.1 * c0 * (point->C_input[i + inputShift] - globals[i]) * point->C_constants[1];//changed
				//}

				double dw0[4];
				for (int i = 0; i < 4; i++)
				{
					dw0[i] = -1. * dt * rhs_w[i];
					point->ConstantsOutput[i + 6] = dw0[i];
					point->output[i + outputShift] = point->C_input[i] + .5 * dw0[i];
				}

				localCount[executeFnc]++;
			}
			break;

			// 2 & 3

			case(2):
			{
				int inputShift  = (executeFnc - 1) * 4;
				int outputShift = (executeFnc + 0) * 4;
				int forwardIndex = 8;
				for (int i = 0; i < forwardIndex; i++) 
				{
					point->output[i] = point->C_input[i];
				}
				
				// compute the dc * rho * Laplace(u)
				double rc = point->C_input[0 + inputShift];
				double rn = point->N_input[0 + inputShift];
				double rs = point->S_input[0 + inputShift];
				double re = point->E_input[0 + inputShift];
				double rw = point->W_input[0 + inputShift];
				double laplaceRho = (re*re + rw*rw + rn*rn + rs*rs) * .25 - (rc*rc);
				point->ConstantsOutput[2] = laplaceRho;

				for (int i = 0; i < 2; ++i) 
				{
					double uc = point->C_input[i + 1 + inputShift] / rc;
					double un = point->N_input[i + 1 + inputShift] / rn;
					double us = point->S_input[i + 1 + inputShift] / rs;
					double ue = point->E_input[i + 1 + inputShift] / re;
					double uw = point->W_input[i + 1 + inputShift] / rw;
					double laplaceU = (ue + uw + un + us) * .25 - uc;					
					point->ConstantsOutput[i + 1 + 2] = laplaceU;
				}

				double pc = point->C_input[3 + inputShift];
				double pn = point->N_input[3 + inputShift];
				double ps = point->S_input[3 + inputShift];
				double pe = point->E_input[3 + inputShift];
				double pw = point->W_input[3 + inputShift];
				double laplaceP = (pe + pw + pn + ps) * .25 - pc;
				point->ConstantsOutput[5] = laplaceP;
				
				localCount[executeFnc]++;
			}
			break;
			case(3):
			{
				int inputShift  = (executeFnc - 2) * 4;
				int outputShift = (executeFnc - 1) * 4;
				int forwardIndex = 8;
				for (int i = 0; i < forwardIndex; i++) 
				{
					point->output[i] = point->C_input[i];
				}

				double rhoV_n = point->N_input[0 + inputShift] * point->N_input[2 + inputShift];
				double rhoV_s = point->S_input[0 + inputShift] * point->S_input[2 + inputShift];
				double rhoU_e = point->E_input[0 + inputShift] * point->E_input[1 + inputShift];
				double rhoU_w = point->W_input[0 + inputShift] * point->W_input[1 + inputShift];
				
				double mass = (rhoU_e - rhoU_w) / (2 * dx) + (rhoV_s - rhoV_n) / (2 * dy);
				
				point->C_conserved[0] = mass;
				
				//////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////////////////////////////

				double rhoU_c = point->C_input[0 + inputShift] * point->C_input[1 + inputShift];
				double rhoV_c = point->C_input[0 + inputShift] * point->C_input[2 + inputShift];

				double un = point->N_input[1 + inputShift] / point->N_input[0 + inputShift];
				double us = point->S_input[1 + inputShift] / point->S_input[0 + inputShift];
				double ue = point->E_input[1 + inputShift] / point->E_input[0 + inputShift];
				double uw = point->W_input[1 + inputShift] / point->W_input[0 + inputShift];

				double rhoUV_n = point->N_input[1 + inputShift] * point->N_input[2 + inputShift];
				double rhoUV_s = point->S_input[1 + inputShift] * point->S_input[2 + inputShift];
				double rhoUU_e = point->E_input[1 + inputShift] * point->E_input[1 + inputShift];
				double rhoUU_w = point->W_input[1 + inputShift] * point->W_input[1 + inputShift];

				double pe = point->E_input[3 + inputShift];
				double pw = point->W_input[3 + inputShift];

				
				double momentum_x = ((rhoUU_e - rhoUU_w) / (2 * dx) + rhoU_c * (ue - uw) / (2 * dx)) / 2.0
					              + ((rhoUV_s - rhoUV_n) / (2 * dy) + rhoV_c * (us - un) / (2 * dy)) / 2.0
								  + (pe - pw) / (2 * dx);
							
				point->C_conserved[1] = momentum_x;
				
				double vn = point->N_input[2 + inputShift] / point->N_input[0 + inputShift];
				double vs = point->S_input[2 + inputShift] / point->S_input[0 + inputShift];
				double ve = point->E_input[2 + inputShift] / point->E_input[0 + inputShift];
				double vw = point->W_input[2 + inputShift] / point->W_input[0 + inputShift];

				double rhoVV_n = point->N_input[2 + inputShift] * point->N_input[2 + inputShift];
				double rhoVV_s = point->S_input[2 + inputShift] * point->S_input[2 + inputShift];
				double rhoUV_e = point->E_input[1 + inputShift] * point->E_input[2 + inputShift];
				double rhoUV_w = point->W_input[1 + inputShift] * point->W_input[2 + inputShift];

				double pn = point->N_input[3 + inputShift];
				double ps = point->S_input[3 + inputShift];

				
				double momentum_y = ((rhoUV_e - rhoUV_w) / (2 * dx) + rhoU_c * (ve - vw) / (2 * dx)) / 2.0
					              + ((rhoVV_s - rhoVV_n) / (2 * dy) + rhoV_c * (vs - vn) / (2 * dy)) / 2.0
								  + (ps - pn) / (2 * dy);
				
				point->C_conserved[2] = momentum_y;
				
				//////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////////////////////////////
				
				
				double gamma = 1.4;
				double uc = point->C_input[1 + inputShift] / point->C_input[0 + inputShift];
				double vc = point->C_input[2 + inputShift] / point->C_input[0 + inputShift];
				
				double diffpu = (pe*ue - pw*uw) / (2 * dx);
				double diffpv = (ps*vs - pn*vn) / (2 * dy);
				double udiffp = uc * ((pe - pw) / (2 * dx));
				double vdiffp = vc * ((ps - pn) / (2 * dy));

				double energy = gamma     * (diffpu + diffpv)
					          - (gamma-1) * (udiffp + vdiffp);

				point->C_conserved[3] = energy;
				
				//////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////////////////////////////
				// Compute dissipation
				
				// compute the Laplacian of dc * rho * laplace(u)
				double laplaceDiss[4];

				double diss_c = DISS_COEFF * (1 * 1) * point->C_constants[2];
				double diss_n = DISS_COEFF * (1 * 1) * point->N_constants[2];
				double diss_s = DISS_COEFF * (1 * 1) * point->S_constants[2];
				double diss_e = DISS_COEFF * (1 * 1) * point->E_constants[2];
				double diss_w = DISS_COEFF * (1 * 1) * point->W_constants[2];
				laplaceDiss[0] = (diss_w + diss_e + diss_n + diss_s) * .25 - diss_c;

				for (int i = 0; i < 2; ++i) 
				{
					diss_c = DISS_COEFF * (point->C_input[0 + inputShift] * point->C_input[0 + inputShift]) * point->C_constants[i + 1 + 2];
					diss_n = DISS_COEFF * (point->N_input[0 + inputShift] * point->N_input[0 + inputShift]) * point->N_constants[i + 1 + 2];
					diss_s = DISS_COEFF * (point->S_input[0 + inputShift] * point->S_input[0 + inputShift]) * point->S_constants[i + 1 + 2];
					diss_e = DISS_COEFF * (point->E_input[0 + inputShift] * point->E_input[0 + inputShift]) * point->E_constants[i + 1 + 2];
					diss_w = DISS_COEFF * (point->W_input[0 + inputShift] * point->W_input[0 + inputShift]) * point->W_constants[i + 1 + 2];

					laplaceDiss[i+1] = (diss_w + diss_e + diss_n + diss_s) * .25 - diss_c;
				}

				diss_c = DISS_COEFF * (1 * 1) * point->C_constants[5];
				diss_n = DISS_COEFF * (1 * 1) * point->N_constants[5];
				diss_s = DISS_COEFF * (1 * 1) * point->S_constants[5];
				diss_e = DISS_COEFF * (1 * 1) * point->E_constants[5];
				diss_w = DISS_COEFF * (1 * 1) * point->W_constants[5];
				laplaceDiss[3] = (diss_w + diss_e + diss_n + diss_s) * .25 - diss_c;
				

				const double R_air = 287., T0 = 300.;
				const double c0 = sqrt(gamma * R_air * T0);
				double dissipation_r = laplaceDiss[0] * c0 / dx;
				double dissipation_x = laplaceDiss[1] * c0 / dx;
				double dissipation_y = laplaceDiss[2] * c0 / dy;
				double dissipation_p = laplaceDiss[3] * c0 / dx;
				
				//dissipation_r = 0;
				//dissipation_p = 0;
				mass       += dissipation_r;
				momentum_x += dissipation_x;
				momentum_y += dissipation_y;
				energy     += dissipation_p - (gamma - 1) * (uc * dissipation_x + vc * dissipation_y);
				
				double rhs_w[4];
				double rc = point->C_input[0 + inputShift];
				rhs_w[0] = 0.5 * mass / rc;
				rhs_w[1] = momentum_x / rc;
				rhs_w[2] = momentum_y / rc;
				rhs_w[3] = energy;

				//rhs_w[i_obstacle,j_obstacle,1:3] += 0.1*c0 * w[i_obstacle, j_obstacle, 1:3]
				// rhs_w[1] += 0.5 * c0 * point->C_input[1 + inputShift] * point->C_constants[0];;	//changed
				// rhs_w[2] += 0.5 * c0 * point->C_input[2 + inputShift] * point->C_constants[0];; // changed
				
				//rhs_w += 0.1*c0 * (w - w0) * dc[:,:,newaxis]
				// for (int i = 0; i < 4; i++) 
				// {
				// 	rhs_w[i] += 0.1 * c0 * (point->C_input[i + inputShift] - globals[i]) * point->C_constants[1];//changed
				// }
				
				double dw1[4];
				for (int i = 0; i < 4; i++)
				{
					dw1[i] = -1. * dt * rhs_w[i];
					point->output[i + inputShift] = point->C_constants[i + 6];
					point->ConstantsOutput[i + 6] = dw1[i];
					point->output[i + outputShift] = point->C_input[i] + .5 * dw1[i];
				}
				localCount[executeFnc]++;
			}
			break;

			//

			//4 & 5

			case(4):
			{
				int inputShift  = (executeFnc - 2) * 4;
				int outputShift = (executeFnc - 1) * 4;
				int forwardIndex = 12;
				for (int i = 0; i < forwardIndex; i++) 
				{
					point->output[i] = point->C_input[i];
				}
				
				// compute the dc * rho * Laplace(u)
				double rc = point->C_input[0 + inputShift];
				double rn = point->N_input[0 + inputShift];
				double rs = point->S_input[0 + inputShift];
				double re = point->E_input[0 + inputShift];
				double rw = point->W_input[0 + inputShift];
				double laplaceRho = (re*re + rw*rw + rn*rn + rs*rs) * .25 - (rc*rc);
				point->ConstantsOutput[2] = laplaceRho;

				for (int i = 0; i < 2; ++i) 
				{
					double uc = point->C_input[i + 1 + inputShift] / rc;
					double un = point->N_input[i + 1 + inputShift] / rn;
					double us = point->S_input[i + 1 + inputShift] / rs;
					double ue = point->E_input[i + 1 + inputShift] / re;
					double uw = point->W_input[i + 1 + inputShift] / rw;
					double laplaceU = (ue + uw + un + us) * .25 - uc;					
					point->ConstantsOutput[i + 1 + 2] = laplaceU;
				}

				double pc = point->C_input[3 + inputShift];
				double pn = point->N_input[3 + inputShift];
				double ps = point->S_input[3 + inputShift];
				double pe = point->E_input[3 + inputShift];
				double pw = point->W_input[3 + inputShift];
				double laplaceP = (pe + pw + pn + ps) * .25 - pc;
				point->ConstantsOutput[5] = laplaceP;
				
				localCount[executeFnc]++;
			}
			break;
			case(5):
			{
				int inputShift  = (executeFnc - 3) * 4;
				int outputShift = (executeFnc - 2) * 4;
				int forwardIndex = 12;
				for (int i = 0; i < forwardIndex; i++) 
				{
					point->output[i] = point->C_input[i];
				}

				double rhoV_n = point->N_input[0 + inputShift] * point->N_input[2 + inputShift];
				double rhoV_s = point->S_input[0 + inputShift] * point->S_input[2 + inputShift];
				double rhoU_e = point->E_input[0 + inputShift] * point->E_input[1 + inputShift];
				double rhoU_w = point->W_input[0 + inputShift] * point->W_input[1 + inputShift];
				
				double mass = (rhoU_e - rhoU_w) / (2 * dx) + (rhoV_s - rhoV_n) / (2 * dy);
				
				point->C_conserved[0] = mass;
				
				//////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////////////////////////////

				double rhoU_c = point->C_input[0 + inputShift] * point->C_input[1 + inputShift];
				double rhoV_c = point->C_input[0 + inputShift] * point->C_input[2 + inputShift];

				double un = point->N_input[1 + inputShift] / point->N_input[0 + inputShift];
				double us = point->S_input[1 + inputShift] / point->S_input[0 + inputShift];
				double ue = point->E_input[1 + inputShift] / point->E_input[0 + inputShift];
				double uw = point->W_input[1 + inputShift] / point->W_input[0 + inputShift];

				double rhoUV_n = point->N_input[1 + inputShift] * point->N_input[2 + inputShift];
				double rhoUV_s = point->S_input[1 + inputShift] * point->S_input[2 + inputShift];
				double rhoUU_e = point->E_input[1 + inputShift] * point->E_input[1 + inputShift];
				double rhoUU_w = point->W_input[1 + inputShift] * point->W_input[1 + inputShift];

				double pe = point->E_input[3 + inputShift];
				double pw = point->W_input[3 + inputShift];

				
				double momentum_x = ((rhoUU_e - rhoUU_w) / (2 * dx) + rhoU_c * (ue - uw) / (2 * dx)) / 2.0
					              + ((rhoUV_s - rhoUV_n) / (2 * dy) + rhoV_c * (us - un) / (2 * dy)) / 2.0
								  + (pe - pw) / (2 * dx);
							
				point->C_conserved[1] = momentum_x;
				
				double vn = point->N_input[2 + inputShift] / point->N_input[0 + inputShift];
				double vs = point->S_input[2 + inputShift] / point->S_input[0 + inputShift];
				double ve = point->E_input[2 + inputShift] / point->E_input[0 + inputShift];
				double vw = point->W_input[2 + inputShift] / point->W_input[0 + inputShift];

				double rhoVV_n = point->N_input[2 + inputShift] * point->N_input[2 + inputShift];
				double rhoVV_s = point->S_input[2 + inputShift] * point->S_input[2 + inputShift];
				double rhoUV_e = point->E_input[1 + inputShift] * point->E_input[2 + inputShift];
				double rhoUV_w = point->W_input[1 + inputShift] * point->W_input[2 + inputShift];

				double pn = point->N_input[3 + inputShift];
				double ps = point->S_input[3 + inputShift];

				
				double momentum_y = ((rhoUV_e - rhoUV_w) / (2 * dx) + rhoU_c * (ve - vw) / (2 * dx)) / 2.0
					              + ((rhoVV_s - rhoVV_n) / (2 * dy) + rhoV_c * (vs - vn) / (2 * dy)) / 2.0
								  + (ps - pn) / (2 * dy);
				
				point->C_conserved[2] = momentum_y;
				
				//////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////////////////////////////
				
				
				double gamma = 1.4;
				double uc = point->C_input[1 + inputShift] / point->C_input[0 + inputShift];
				double vc = point->C_input[2 + inputShift] / point->C_input[0 + inputShift];
				
				double diffpu = (pe*ue - pw*uw) / (2 * dx);
				double diffpv = (ps*vs - pn*vn) / (2 * dy);
				double udiffp = uc * ((pe - pw) / (2 * dx));
				double vdiffp = vc * ((ps - pn) / (2 * dy));

				double energy = gamma     * (diffpu + diffpv)
					          - (gamma-1) * (udiffp + vdiffp);

				point->C_conserved[3] = energy;
				
				//////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////////////////////////////
				// Compute dissipation
				
				// compute the Laplacian of dc * rho * laplace(u)
				double laplaceDiss[4];

				double diss_c = DISS_COEFF * (1 * 1) * point->C_constants[2];
				double diss_n = DISS_COEFF * (1 * 1) * point->N_constants[2];
				double diss_s = DISS_COEFF * (1 * 1) * point->S_constants[2];
				double diss_e = DISS_COEFF * (1 * 1) * point->E_constants[2];
				double diss_w = DISS_COEFF * (1 * 1) * point->W_constants[2];
				laplaceDiss[0] = (diss_w + diss_e + diss_n + diss_s) * .25 - diss_c;

				for (int i = 0; i < 2; ++i) 
				{
					diss_c = DISS_COEFF * (point->C_input[0 + inputShift] * point->C_input[0 + inputShift]) * point->C_constants[i + 1 + 2];
					diss_n = DISS_COEFF * (point->N_input[0 + inputShift] * point->N_input[0 + inputShift]) * point->N_constants[i + 1 + 2];
					diss_s = DISS_COEFF * (point->S_input[0 + inputShift] * point->S_input[0 + inputShift]) * point->S_constants[i + 1 + 2];
					diss_e = DISS_COEFF * (point->E_input[0 + inputShift] * point->E_input[0 + inputShift]) * point->E_constants[i + 1 + 2];
					diss_w = DISS_COEFF * (point->W_input[0 + inputShift] * point->W_input[0 + inputShift]) * point->W_constants[i + 1 + 2];

					laplaceDiss[i+1] = (diss_w + diss_e + diss_n + diss_s) * .25 - diss_c;
				}

				diss_c = DISS_COEFF * (1 * 1) * point->C_constants[5];
				diss_n = DISS_COEFF * (1 * 1) * point->N_constants[5];
				diss_s = DISS_COEFF * (1 * 1) * point->S_constants[5];
				diss_e = DISS_COEFF * (1 * 1) * point->E_constants[5];
				diss_w = DISS_COEFF * (1 * 1) * point->W_constants[5];
				laplaceDiss[3] = (diss_w + diss_e + diss_n + diss_s) * .25 - diss_c;
				

				const double R_air = 287., T0 = 300.;
				const double c0 = sqrt(gamma * R_air * T0);
				double dissipation_r = laplaceDiss[0] * c0 / dx;
				double dissipation_x = laplaceDiss[1] * c0 / dx;
				double dissipation_y = laplaceDiss[2] * c0 / dy;
				double dissipation_p = laplaceDiss[3] * c0 / dx;
				
				//dissipation_r = 0;
				//dissipation_p = 0;
				mass       += dissipation_r;
				momentum_x += dissipation_x;
				momentum_y += dissipation_y;
				energy     += dissipation_p - (gamma - 1) * (uc * dissipation_x + vc * dissipation_y);
				
				double rhs_w[4];
				double rc = point->C_input[0 + inputShift];
				rhs_w[0] = 0.5 * mass / rc;
				rhs_w[1] = momentum_x / rc;
				rhs_w[2] = momentum_y / rc;
				rhs_w[3] = energy;

				//rhs_w[i_obstacle,j_obstacle,1:3] += 0.1*c0 * w[i_obstacle, j_obstacle, 1:3]
				// rhs_w[1] += 0.5 * c0 * point->C_input[1 + inputShift] * point->C_constants[0];; //changed				
				// rhs_w[2] += 0.5 * c0 * point->C_input[2 + inputShift] * point->C_constants[0];; //changed
				
				//rhs_w += 0.1*c0 * (w - w0) * dc[:,:,newaxis]
				//for (int i = 0; i < 4; i++) 
				//{
				//	rhs_w[i] += 0.1 * c0 * (point->C_input[i + inputShift] - globals[i]) * point->C_constants[1];//changed
				//}
				for (int i = 0; i < 4; i++)
					point->output[i + outputShift] = rhs_w[i];

				double dw2[4];
				for (int i = 0; i < 4; i++)
				{
					dw2[i] = -1. * dt * rhs_w[i];
					point->output[i + inputShift] = point->C_constants[i + 6];
					point->ConstantsOutput[i + 6] = dw2[i];
					point->output[i + outputShift] = point->C_input[i] + dw2[i];
				}
				localCount[executeFnc]++;
			}
			break;

			//

			// 6 & 7

			case(6):
			{
				int inputShift  = (executeFnc - 3) * 4;
				int outputShift = (executeFnc - 6) * 4;
				int forwardIndex = 16;
				for (int i = 0; i < forwardIndex; i++) 
				{
					point->output[i] = point->C_input[i];
				}
				
				// compute the dc * rho * Laplace(u)
				double rc = point->C_input[0 + inputShift];
				double rn = point->N_input[0 + inputShift];
				double rs = point->S_input[0 + inputShift];
				double re = point->E_input[0 + inputShift];
				double rw = point->W_input[0 + inputShift];
				double laplaceRho = (re*re + rw*rw + rn*rn + rs*rs) * .25 - (rc*rc);
				point->ConstantsOutput[2] = laplaceRho;

				for (int i = 0; i < 2; ++i) 
				{
					double uc = point->C_input[i + 1 + inputShift] / rc;
					double un = point->N_input[i + 1 + inputShift] / rn;
					double us = point->S_input[i + 1 + inputShift] / rs;
					double ue = point->E_input[i + 1 + inputShift] / re;
					double uw = point->W_input[i + 1 + inputShift] / rw;
					double laplaceU = (ue + uw + un + us) * .25 - uc;					
					point->ConstantsOutput[i + 1 + 2] = laplaceU;
				}

				double pc = point->C_input[3 + inputShift];
				double pn = point->N_input[3 + inputShift];
				double ps = point->S_input[3 + inputShift];
				double pe = point->E_input[3 + inputShift];
				double pw = point->W_input[3 + inputShift];
				double laplaceP = (pe + pw + pn + ps) * .25 - pc;
				point->ConstantsOutput[5] = laplaceP;
				
				localCount[executeFnc]++;
			}
			break;
			case(7):
			{
				int inputShift  = (executeFnc - 4) * 4;
				int outputShift = (executeFnc - 7) * 4;
				int forwardIndex = 0;
				for (int i = 0; i < forwardIndex; i++) 
				{
					point->output[i] = point->C_input[i];
				}

				double rhoV_n = point->N_input[0 + inputShift] * point->N_input[2 + inputShift];
				double rhoV_s = point->S_input[0 + inputShift] * point->S_input[2 + inputShift];
				double rhoU_e = point->E_input[0 + inputShift] * point->E_input[1 + inputShift];
				double rhoU_w = point->W_input[0 + inputShift] * point->W_input[1 + inputShift];
				
				double mass = (rhoU_e - rhoU_w) / (2 * dx) + (rhoV_s - rhoV_n) / (2 * dy);
				
				point->C_conserved[0] = mass;
				
				//////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////////////////////////////

				double rhoU_c = point->C_input[0 + inputShift] * point->C_input[1 + inputShift];
				double rhoV_c = point->C_input[0 + inputShift] * point->C_input[2 + inputShift];

				double un = point->N_input[1 + inputShift] / point->N_input[0 + inputShift];
				double us = point->S_input[1 + inputShift] / point->S_input[0 + inputShift];
				double ue = point->E_input[1 + inputShift] / point->E_input[0 + inputShift];
				double uw = point->W_input[1 + inputShift] / point->W_input[0 + inputShift];

				double rhoUV_n = point->N_input[1 + inputShift] * point->N_input[2 + inputShift];
				double rhoUV_s = point->S_input[1 + inputShift] * point->S_input[2 + inputShift];
				double rhoUU_e = point->E_input[1 + inputShift] * point->E_input[1 + inputShift];
				double rhoUU_w = point->W_input[1 + inputShift] * point->W_input[1 + inputShift];

				double pe = point->E_input[3 + inputShift];
				double pw = point->W_input[3 + inputShift];

				
				double momentum_x = ((rhoUU_e - rhoUU_w) / (2 * dx) + rhoU_c * (ue - uw) / (2 * dx)) / 2.0
					              + ((rhoUV_s - rhoUV_n) / (2 * dy) + rhoV_c * (us - un) / (2 * dy)) / 2.0
								  + (pe - pw) / (2 * dx);
							
				point->C_conserved[1] = momentum_x;
				
				double vn = point->N_input[2 + inputShift] / point->N_input[0 + inputShift];
				double vs = point->S_input[2 + inputShift] / point->S_input[0 + inputShift];
				double ve = point->E_input[2 + inputShift] / point->E_input[0 + inputShift];
				double vw = point->W_input[2 + inputShift] / point->W_input[0 + inputShift];

				double rhoVV_n = point->N_input[2 + inputShift] * point->N_input[2 + inputShift];
				double rhoVV_s = point->S_input[2 + inputShift] * point->S_input[2 + inputShift];
				double rhoUV_e = point->E_input[1 + inputShift] * point->E_input[2 + inputShift];
				double rhoUV_w = point->W_input[1 + inputShift] * point->W_input[2 + inputShift];

				double pn = point->N_input[3 + inputShift];
				double ps = point->S_input[3 + inputShift];

				
				double momentum_y = ((rhoUV_e - rhoUV_w) / (2 * dx) + rhoU_c * (ve - vw) / (2 * dx)) / 2.0
					              + ((rhoVV_s - rhoVV_n) / (2 * dy) + rhoV_c * (vs - vn) / (2 * dy)) / 2.0
								  + (ps - pn) / (2 * dy);
				
				point->C_conserved[2] = momentum_y;
				
				//////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////////////////////////////
				
				
				double gamma = 1.4;
				double uc = point->C_input[1 + inputShift] / point->C_input[0 + inputShift];
				double vc = point->C_input[2 + inputShift] / point->C_input[0 + inputShift];
				
				double diffpu = (pe*ue - pw*uw) / (2 * dx);
				double diffpv = (ps*vs - pn*vn) / (2 * dy);
				double udiffp = uc * ((pe - pw) / (2 * dx));
				double vdiffp = vc * ((ps - pn) / (2 * dy));

				double energy = gamma     * (diffpu + diffpv)
					          - (gamma-1) * (udiffp + vdiffp);

				point->C_conserved[3] = energy;
				
				//////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////////////////////////////
				// Compute dissipation
				
				// compute the Laplacian of dc * rho * laplace(u)
				double laplaceDiss[4];

				double diss_c = DISS_COEFF * (1 * 1) * point->C_constants[2];
				double diss_n = DISS_COEFF * (1 * 1) * point->N_constants[2];
				double diss_s = DISS_COEFF * (1 * 1) * point->S_constants[2];
				double diss_e = DISS_COEFF * (1 * 1) * point->E_constants[2];
				double diss_w = DISS_COEFF * (1 * 1) * point->W_constants[2];
				laplaceDiss[0] = (diss_w + diss_e + diss_n + diss_s) * .25 - diss_c;

				for (int i = 0; i < 2; ++i) 
				{
					diss_c = DISS_COEFF * (point->C_input[0 + inputShift] * point->C_input[0 + inputShift]) * point->C_constants[i + 1 + 2];
					diss_n = DISS_COEFF * (point->N_input[0 + inputShift] * point->N_input[0 + inputShift]) * point->N_constants[i + 1 + 2];
					diss_s = DISS_COEFF * (point->S_input[0 + inputShift] * point->S_input[0 + inputShift]) * point->S_constants[i + 1 + 2];
					diss_e = DISS_COEFF * (point->E_input[0 + inputShift] * point->E_input[0 + inputShift]) * point->E_constants[i + 1 + 2];
					diss_w = DISS_COEFF * (point->W_input[0 + inputShift] * point->W_input[0 + inputShift]) * point->W_constants[i + 1 + 2];

					laplaceDiss[i+1] = (diss_w + diss_e + diss_n + diss_s) * .25 - diss_c;
				}

				diss_c = DISS_COEFF * (1 * 1) * point->C_constants[5];
				diss_n = DISS_COEFF * (1 * 1) * point->N_constants[5];
				diss_s = DISS_COEFF * (1 * 1) * point->S_constants[5];
				diss_e = DISS_COEFF * (1 * 1) * point->E_constants[5];
				diss_w = DISS_COEFF * (1 * 1) * point->W_constants[5];
				laplaceDiss[3] = (diss_w + diss_e + diss_n + diss_s) * .25 - diss_c;
				

				const double R_air = 287., T0 = 300.;
				const double c0 = sqrt(gamma * R_air * T0);
				double dissipation_r = laplaceDiss[0] * c0 / dx;
				double dissipation_x = laplaceDiss[1] * c0 / dx;
				double dissipation_y = laplaceDiss[2] * c0 / dy;
				double dissipation_p = laplaceDiss[3] * c0 / dx;
				
				//dissipation_r = 0;
				//dissipation_p = 0;
				mass       += dissipation_r;
				momentum_x += dissipation_x;
				momentum_y += dissipation_y;
				energy     += dissipation_p - (gamma - 1) * (uc * dissipation_x + vc * dissipation_y);
				
				double rhs_w[4];
				double rc = point->C_input[0 + inputShift];
				rhs_w[0] = 0.5 * mass / rc;
				rhs_w[1] = momentum_x / rc;
				rhs_w[2] = momentum_y / rc;
				rhs_w[3] = energy;

				//rhs_w[i_obstacle,j_obstacle,1:3] += 0.1*c0 * w[i_obstacle, j_obstacle, 1:3]
				// rhs_w[1] += 0.5 * c0 * point->C_input[1 + inputShift] * point->C_constants[0]; //changed
				// rhs_w[2] += 0.5 * c0 * point->C_input[2 + inputShift] * point->C_constants[0]; //changed
				
				//rhs_w += 0.1*c0 * (w - w0) * dc[:,:,newaxis]
				// for (int i = 0; i < 4; i++) 
				// {
				// 	rhs_w[i] += 0.1 * c0 * (point->C_input[i + inputShift] - globals[i]) * point->C_constants[1];//changed
				// }
				for (int i = 0; i < 4; i++)
					point->output[i + outputShift] = rhs_w[i];

				double dw3[4];
				for (int i = 0; i < 4; i++)
				{
					dw3[i] = -1. * dt * rhs_w[i];
					point->output[i] = point->C_input[i] + (((point->C_input[i + 4] + dw3[i]) / 6.) + ((point->C_input[i + 8] + point->C_constants[i + 6]) / 3.));
				}

				point->output[1] *= (1 - point->C_constants[0]);			
				point->output[2] *= (1 - point->C_constants[0]);

				for (int i = 0; i < 4; i++)
				{
					point->output[i] += (globals[i] - point->output[i]) * point->C_constants[1];
				}

				localCount[executeFnc]++;
			}
			break;

			//

		}
	}

	static void conservationCheck(double *previous,double *current)
	{
		
	}
	static void linearSystemUpdate(int i,int j,InitPoint2D *point)
	{
		printf("Updating(%d,%d)...\n",i,j);
	}


	static void init(int i,int j,InitPoint2D *point)
	{
		const double PI = atan(1.0) * 4;
		double x = (i + 0.5) * dx - 0.2 * lx;
		double y = (j + 0.5) * dy - 0.5 * ly;		
		double dc = pow(cos((x / lx + 0.2) * PI),256);

		point->U_constants[1] = dc;
		point->U_constants[0] = 0.0; //exp(-pow((x*x + y*y)/1,8));
		point->U_constants[0] = exp(-pow((x*x + y*y)/1,8));
		//point->U_input[16]    = exp(-pow((x*x + y*y)/1,8));

		for (int ii = 0; ii < 4; ++ii) 
		{
			point->U_input[ii] = globals[ii];			
		}
				
		point->U_conserved[0] = 0;
		point->U_conserved[1] = 0;
		point->U_conserved[2] = 0;
		point->U_conserved[3] = 0;
		
	}

	static void setGlobalVariables(int nx,int ny)
	{
		double gamma = 1.4;
		double R     = 287.;
		double T0    = 300.;
		double p0    = 101325.;
		double M0    = 0.2;
		double c0    = sqrt(gamma * R * T0);
		const double rho0 = p0 / (R * T0);  
		const double u0   = c0 * M0;  
		
		globals[0] = sqrt(rho0);
		globals[1] = sqrt(rho0) * u0;
		globals[2] = sqrt(rho0) * 0.0;
		globals[3] = p0;

		int xObstacleStart = 0;
		int xObstacleEnd   = 0;
		int yObstacleStart = 0;
		int yObstacleEnd   = 0;

		for(int i=0;i<nx;i++)
		{
			double x = (i + 0.5) * dx - 0.2 * lx;
			if(x < -1.) xObstacleStart++;
			if(x <  1.) xObstacleEnd++;			
		}

		for(int j=0;j<ny;j++)
		{
			double y = (j + 0.5) * dy - 0.5 * ly;
			if(y < -.5) yObstacleStart++;
			if(y <  .5) yObstacleEnd++;
		}
		
		globals[4] = xObstacleStart;
		globals[5] = xObstacleEnd;
		globals[6] = yObstacleStart;
		globals[7] = yObstacleEnd;
		if(pg.rank == 0)
		printf("Global Variables Setup DONE! - X_Obstacle(%d,%d) , Y_Obstacle(%d,%d)\n",xObstacleStart,xObstacleEnd,yObstacleStart,yObstacleEnd);
	}

	static void setInitParameters(int nx,int ny)
	{
		//PARAMETERS SETUP
		totalGlobals = 8;
		constants = 10;
		remoteConstants = 1;
		substeps = 8;
		dataPointSize = 2;
		
		//OUTPUT SETUP AND CONTROL
		outputDirectory = "c://output";
		//fileOutputEvery = 10;
		
		//CONSERVATION SETUP AND CONTROL
		totalConservedQuantities = 4;
		//conservationCheckEvery = 10;

		//SIMULATION SETUP
		lx = 15.;
		ly = 10.;
		//cycles = 100;
		dx = lx/nx;
		dy = ly/ny;
		
		double c0 = sqrt(1.4 * 287. * 300.);
		dt = .000001 ;
		dt = dx / c0 * 0.5;
	}	
};

#endif
