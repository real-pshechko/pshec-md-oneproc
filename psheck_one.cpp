#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <cassert>
#include <vector>
#include <random>

#include "mpi.h" 
#include "defines_vars_and_structs.h"
#include "defines_calc.h"
#include "my_func.h"
#include "windows.h"

using namespace std;

//ÑÒÐÓÊÒÓÐÛ
CalcSpace Space;

//ÊËÀÑÑÛ
CalcConstsLD Consts;//(6.2E-10, 3.16E-21, 4.055E-10, 0.24E-12, 21.75E-26);


//ÏÅÐÅÌÅÍÍÛÅ Ñ ÏËÀÂÀÞÙÅÉ ÒÎ×ÊÎÉ
real dt, lx, ly, lz, //äëèíà ÿ÷åéêè â ñïèñêå Âåðëå
	rCut, rView; //ðàäèóñ îáðåçêè è ðàäèóñ íàáëþäåíèÿ
real BuferM;

real R2cut, X2min, rVr, rVr1, rVr2, aVr, aVr2, aLJ3, bLJ2, aLJ, bLJ, XX, FV, Ur;

real m1, m2, m3;

//ÖÅËÎ×ÈÑËÅÍÍÛÅ È ÏÐÎ×ÈÅ ÏÅÐÅÌÅÍÍÛÅ
int my_num, num_proc;

int NGen, MyN, NpartW, NAsnd, NArcv, NpartA, NAinf0;
int NeibsM,  Neibs0, Neib, Neibs;
int Ncall;
int Xi, Yi, Zi;
int Np1, Np2;
int k;


int main(int argc, char **argv)
{	
	
	ofstream file("output(eobK).txt"), impfile("impulse(eobK).txt"), forces("force.txt");
	ofstream fff("sosed.txt");
	double rty;
	real lx, ly, lz;
	real Fx, Fy, Fz;
	int timestep;
	X2min = 1.2599210498948732;
	//int i, j, k;
	//file.write('VARIABLES="x, nm","y, nm","p-max, GPa"\n')
	file << "VARIABLES= time,sec , T,K , Eob" << endl;

	Ncall = 0; NAinf0 = 0;

	Space.Lx = 5E-9; Space.Ly = 5E-9; Space.Lz = 5E-9;
	Space.n1 = 10; Space.n2 = 10; Space.n3 = 10;

	rCut = 4 * Consts.sig;
	rView = 1.5 * Consts.sig;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_num);
	MPI_Status status;
	MPI_Request request;
	
	NullPtrsInt TestSend(num_proc), TestRecv(num_proc), MCHKL(num_proc); // ÂÛÄÅËßÒÜ ÈÇ ÑÒÅÊÀ ÈËÈ ÄÈÍÀÌÈ×ÅÑÊÈ ÈÇ ÊÓ×È(÷åðåç pointer) ???

	NullPtrsReal Bufer(num_proc), rMax(num_proc);

	NullPtrsVecR Cent(num_proc);

	Period prd;

	MPI_Request *ReqR = new MPI_Request [num_proc]; MPI_Request *ReqS = new MPI_Request [num_proc];
	MPI_Request *ReqR1 = new MPI_Request [num_proc]; MPI_Request *ReqS1 = new MPI_Request [num_proc];
	

	Neibs0 = -1;
	m1 = Consts.m;
	m2 = 0; m3 = 0;
	dt = 4E-15;//0.24E-12;


	if (my_num == 0)
	{

	}

	else 
	{
		
		NGen = 0;
		CalculateNGen_skal(Space, Consts);
		//NGen = 2;
		
		NullPtrsVecR dr(NGen), Vel(NGen), Accel(NGen), halfVel(NGen);
		NullPtrsInt EnergyPot(NGen);
		
		NullPtrsVecR LanForce(NGen);

		Period prd;

		InitValuesCoord_skal(Space, Consts, dr); 
		
		/*dr.setValueX(0, 0);dr.setValueY(0, 0);  dr.setValueZ(0, 0); 
		dr.setValueX(1, 3.4E-10); dr.setValueY(1, 0); dr.setValueZ(1, 0);*/
		cout << NGen << endl;
		InitValuesOthers_skal(NGen, Vel, Accel, EnergyPot, LanForce);
		
		
		
		int NxBox = (int)ceil(Space.Lx/rCut); //ðàçáèâàåì îáëàñòü íà ÿ÷åéêè, â êàæäîì íàïðàâëåíèè ÿ÷åéêè
		int NyBox = (int)ceil(Space.Ly/rCut);
		int NzBox = (int)ceil(Space.Lz/rCut);

		int NAllBox = NxBox * NyBox * NzBox;


		lx = Space.Lx/(NxBox);
		ly = Space.Ly/(NyBox);
		lz = Space.Lz/(NzBox);

		int ***Addr = new int ** [NzBox];
		for (int i = 0; i < NzBox; i++) {
			Addr[i] = new int * [NyBox];
			for (int j = 0; j < NyBox; j++) {
				Addr[i][j] = new int [NxBox];
			}
		}

		int ***Addr1 = new int ** [NzBox];
		for (int i = 0; i < NzBox; i++) {
			Addr1[i] = new int * [NyBox];
			for (int j = 0; j < NyBox; j++) {
				Addr1[i][j] = new int [NxBox];
			}
		}

		vector<int> neighb(NGen, -1);
		vector<int> verle, notverle;
		vector<int> cutNeighb, cutVerle;

		vector<int> neighb1;
		vector<int> verle1;
		vector<int> cutNeighb1, cutVerle1;

		vector<double> GaussDist(NGen, 0); 
		

		int *cell = new int [NGen];
		int *cell1 = new int [NGen];

		int timestep;
		int III, ind, n;

		rVr = Consts.sig; aVr = Consts.eps; // ìîæíî âûíåñòè èç öèêëà timestep
		R2cut = rCut * rCut;

		rVr2 = rVr * rVr;//Consts.sig * Consts.sig;
		rVr1 = 1.0 / rVr2;
		aVr2 = R2cut / rVr2 - X2min;
		aLJ3 = (rVr2 / R2cut) * (rVr2 / R2cut) * (rVr2 / R2cut);
		bLJ2 = (aLJ3 - 0.5) * aVr2* rVr2 / R2cut;
		aLJ = (aLJ3 - 1.0 + 2.0 * bLJ2) * aLJ3 / (aVr2 * aVr2);
		bLJ =- (aLJ3 - 1.0 + 3.0 * bLJ2) * aLJ3 / (aVr2 * aVr2 * aVr2);

		aLJ3 = 3.0 * aLJ;
		bLJ2 = 2.0 * bLJ;

		aVr2 = -48.0 * aVr / rVr2;
		aVr =   2.0 * aVr; 

		double impX0 = 0, impY0 = 0, impZ0 = 0;
		

//   ÑÍÀ×ÀËÀ Ñ×ÈÒÀÞ ÑÈËÛ ÂÇÀÈÌÎÄÅÉÑÒÂÈß ÌÅÆÄÓ ÀÒÎÌÀÌÈ È ÂÛÑ×ÈÒÛÂÀÞ ÓÑÊÎÐÅÍÈß (ÓÑÊÎÐÅÍÈß ÍÀ×ÀËÜÍÛÅ ÍÅ ÍÓËÅÂÛÅ!!!)
		for (int i = 0; i < NGen; i++)
		{
			Fx=0.; Fy=0.; Fz=0.;
			for (int j = 0; j < NGen; j++)
			{
				if (i != j)
				{
					real r2 = 
								(dr.getValueX(i) - dr.getValueX(j)) * //verle[cutVerle[i]]
								(dr.getValueX(i) - dr.getValueX(j)) + 
								(dr.getValueY(i) - dr.getValueY(j)) * 
								(dr.getValueY(i) - dr.getValueY(j)) +
								(dr.getValueZ(i) - dr.getValueZ(j)) * 
								(dr.getValueZ(i) - dr.getValueZ(j));
					if (r2 < rCut * rCut) //âûâîäèòü òåõ, ñ êåì âçàèìîäåéñòâóåò
					{							
						XX = rVr1 * r2 - X2min;
						FV = rVr2 / r2;
						Ur = FV * FV * FV;

						real RR = aVr * ( Ur * (Ur - 1.0) - (aLJ3 + bLJ2 * XX) * XX * XX); //V(r)/2     êàêîé ðàçìåðíîñòè ìàññèâ? RR(k)
						real FR = aVr2 * ( Ur * (Ur - 0.5) * FV + (aLJ + bLJ * XX) * XX ); //V'(r)/r ÑÈËÀ ÍÀ ÝÐ(ìåæàòîìíîå ðàññòîÿíèå)

						Fx += FR * (dr.getValueX(j) - dr.getValueX(i));						
						Fy += FR * (dr.getValueY(j) - dr.getValueY(i));
						Fz += FR * (dr.getValueZ(j) - dr.getValueZ(i));									
					}						
				}				
			}			
			Accel.setValueX(i, (Fx/m1));
			Accel.setValueY(i, (Fy/m1));
			Accel.setValueZ(i, (Fz/m1));
		}

		int ActivateThermLan = 1;

		for (timestep = 0; timestep < 2000; timestep++)
		{
			
			cout << timestep << "   step" << endl;
			if (timestep > 250) ActivateThermLan = 0;
			if (ActivateThermLan == 1 )//&& timestep < 9)
			{
				ThermLan(LanForce, Consts, Vel);
			}

			


			for (int i = 0; i < NGen; i++)
			{
				double val = Vel.getValueX(i) + 0.5 * Accel.getValueX(i) * dt;
				halfVel.setValueX(i, val);
					val = Vel.getValueY(i) + 0.5 * Accel.getValueY(i) * dt;
				halfVel.setValueY(i, val);
					val = Vel.getValueZ(i) + 0.5 * Accel.getValueZ(i) * dt;
				halfVel.setValueZ(i, val);
			}

			for (int i = 0; i < NGen; i++) 
			{
				real xi = dt * halfVel.getValueX(i);//Vel.getValueX(i) * dt + 0.5 * Accel.getValueX(i) * dt * dt;
				xi += dr.getValueX(i);
				prd.X(xi, Space.Lx, 0, Space.Lx);
				dr.setValueX(i, xi);
			

				real yi = dt * halfVel.getValueY(i);//Vel.getValueY(i) * dt + 0.5 * Accel.getValueY(i) * dt * dt;
				yi += dr.getValueY(i);
				prd.Y(yi, Space.Ly, 0, Space.Ly);
				dr.setValueY(i, yi);
				
								
				real zi = dt * halfVel.getValueZ(i);//Vel.getValueZ(i) * dt + 0.5 * Accel.getValueZ(i) * dt * dt;
				zi += dr.getValueZ(i);
				prd.Z(zi, Space.Lz, 0, Space.Lz);
				dr.setValueZ(i, zi);
			}

			//int aa = 0;
			
			//vector<real> RR, FR;
			vector<real> U(NGen, 0);
			
			//ïåðåñòðîåíèå ñïèñêà Âåðëå
			


			//ðàñ÷åò ñèë ËÄ
			
			
			
			

			for (int i = 0; i < NGen; i++)
			{
			
			
				Fx=0.; Fy=0.; Fz=0.; double u = 0.;
				for (int j = 0; j < NGen; j++)
				{
					if (i != j)
					{
						real r2 = 
									(dr.getValueX(i) - dr.getValueX(j)) * 
									(dr.getValueX(i) - dr.getValueX(j)) + 
									(dr.getValueY(i) - dr.getValueY(j)) * 
									(dr.getValueY(i) - dr.getValueY(j)) +
									(dr.getValueZ(i) - dr.getValueZ(j)) * 
									(dr.getValueZ(i) - dr.getValueZ(j));
						if (r2 < rCut * rCut)
						{
							
							XX = rVr1 * r2 - X2min;

							FV = rVr2 / r2;

							Ur = FV * FV * FV;

							real RR = aVr * ( Ur * (Ur - 1.0) - (aLJ3 + bLJ2 * XX) * XX * XX); //V(r)/2 

							real FR = aVr2 * ( Ur * (Ur - 0.5) * FV + (aLJ+ bLJ * XX) * XX ); //V'(r)/r ÑÈËÀ ÍÀ ÝÐ(ìåæàòîìíîå ðàññòîÿíèå)

							Fx += FR * (dr.getValueX(j) - dr.getValueX(i));
								
							Fy += FR * (dr.getValueY(j) - dr.getValueY(i));
							
							Fz += FR * (dr.getValueZ(j) - dr.getValueZ(i));

							u += RR;
						}
					}
						
				}

				if (ActivateThermLan == 1) 
				{
					Fx += LanForce.getValueX(i); 
					Fy += LanForce.getValueY(i); 
					Fz += LanForce.getValueZ(i);
				}
				
				
				
			

				Accel.setValueX(i, (Fx/m1));
				Accel.setValueY(i, (Fy/m1));
				Accel.setValueZ(i, (Fz/m1));
				U[i] = u;
			}

			
			//ñêîðîñòè
			for (int i = 0; i < NGen; i++)
			{
				double val = halfVel.getValueX(i) + 0.5 * dt * Accel.getValueX(i);
				Vel.setValueX(i, val);
					val = halfVel.getValueY(i) + 0.5 * dt * Accel.getValueY(i);
				Vel.setValueY(i, val);
					val = halfVel.getValueZ(i) + 0.5 * dt * Accel.getValueZ(i);
				Vel.setValueZ(i, val);
			}

			double impX = 0, impY = 0, impZ = 0;
		
			for (int j = 0; j < NGen; j++)
			{
				impX += Vel.getValueX(j);
				impY += Vel.getValueY(j);
				impZ += Vel.getValueZ(j);
			}

			if (timestep == 0) 
			{
				impX0 = impX / NGen;
				impY0 = impY / NGen;
				impZ0 = impZ / NGen;
			}
			impX = impX / NGen - impX0;
			impY = impY / NGen - impY0;
			impZ = impZ / NGen - impZ0;
			double imp = (impX + impY + impZ) / 3;
			
			

			//impfile << timestep << " step|   " << impX << "   impX|   " << impY << "   impY|   " << impZ << "   impZ|   " << imp << "   imp|" << endl;

			double Eob = 0; double kin;
			double kin1 = 0;
			for (int j = 0; j < NGen; j++)
			{
				kin = 0.5 * m1 * (Vel.getValueX(j) * Vel.getValueX(j) + Vel.getValueY(j) * Vel.getValueY(j) + Vel.getValueZ(j) * Vel.getValueZ(j));
				kin1 += Vel.getValueX(j) * Vel.getValueX(j) + Vel.getValueY(j) * Vel.getValueY(j) + Vel.getValueZ(j) * Vel.getValueZ(j);
				Eob += kin + U[j];
			}
			//cout << Eob/NGen << "   Eob" << endl;
			double Targ = (m1 * (kin1/NGen)) / (3 * 1.38E-23);
			//file << timestep << " step| " << kn/NGen/(1.38E-23) << "   kin|  " << u/NGen/(1.38E-23) << "   pot|  " << Eob/NGen/(1.38E-23) << "  Eob|  " << Targ << "   temperature|" << endl;
			U.clear();

			
			file << timestep * dt << " " << Targ << " " << Eob/NGen/(1.38E-23) << endl;

			//file << timestep << " step|   " << verle.size() << "   verle.size|   " << Eob/NGen/(1.38E-23) << "   Eob|   " << Targ << "   temperature|" << endl;

		}
		


	}




	int workTime = clock();


	cout << ((float)workTime) / CLOCKS_PER_SEC << endl;

	MPI_Finalize();
	return 0;
	getchar();
}

//taskkill /f /im md_pshe_pshe_one_oneproc.exe
