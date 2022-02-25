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

//СТРУКТУРЫ
CalcSpace Space;

//КЛАССЫ
CalcConstsLD Consts;//(6.2E-10, 3.16E-21, 4.055E-10, 0.24E-12, 21.75E-26);


//ПЕРЕМЕННЫЕ С ПЛАВАЮЩЕЙ ТОЧКОЙ
real dt, lx, ly, lz, //длина ячейки в списке Верле
	rCut, rView; //радиус обрезки и радиус наблюдения
real BuferM;

real R2cut, X2min, rVr, rVr1, rVr2, aVr, aVr2, aLJ3, bLJ2, aLJ, bLJ, XX, FV, Ur;

real m1, m2, m3;

//ЦЕЛОЧИСЛЕННЫЕ И ПРОЧИЕ ПЕРЕМЕННЫЕ
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
	
	NullPtrsInt TestSend(num_proc), TestRecv(num_proc), MCHKL(num_proc); // ВЫДЕЛЯТЬ ИЗ СТЕКА ИЛИ ДИНАМИЧЕСКИ ИЗ КУЧИ(через pointer) ???

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
		
		
		
		int NxBox = (int)ceil(Space.Lx/rCut); //разбиваем область на ячейки, в каждом направлении ячейки
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

		rVr = Consts.sig; aVr = Consts.eps; // можно вынести из цикла timestep
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
		

		
/*
		bool flag;
		ind = 0;

		for (int k = 0; k < NzBox; k++)
		{
			for (int j = 0; j < NyBox; j++)
			{
				for (int i = 0; i < NxBox; i++)
				{
					flag = false;
					int i1 = 0, i2 = 0;

					for (int n = 0; n < NGen; n++)
					{
						if (dr.getValueX(n) >= (lx * i) && dr.getValueX(n) <= (lx * (i+1)) &&
							dr.getValueY(n) >= (ly * j) && dr.getValueY(n) <= (ly * (j+1)) &&
							dr.getValueZ(n) >= (lz * k) && dr.getValueZ(n) <= (lz * (k+1)))
						{
							cell1[i1] = n;
							i1 += 1;									
						}
					}

					if (i1 > 0)
					{
						flag = true;
						i2 = i1 - 1;
					}

					if (flag == true)
					{
						qsort(cell1, i2, sizeof(int), Comp1);
						Addr1[k][j][i] = cell1[0];

						for (i1 = 0; i1 < i2; i1++)
						{	
							neighb1[cell1[i1]] = cell1[i1];
						}

						neighb1[cell1[i2]] = -1;
					}
						
					else 
					{
						Addr1[k][j][i] = -1;
					}
				}
			}
		}//адресный список 
		//		
		ind = 0;
		int ind_old = 0; int a = 0;
		cutNeighb1.push_back(0); cutVerle1.push_back(0);

		for (int i = 0; i < NGen; i++)//начало списка Верле (общий)
		{
			if (neighb1[i] == -1) cutNeighb1.push_back(i);
		}
	
		for (int i = 0; i < (cutNeighb1.size() - 1); i++)
		{
			for (int j = cutNeighb1[i]; j < cutNeighb1[i + 1]; j++)
			{
				for (int k = cutNeighb1[i]; k < cutNeighb1[i + 1]; k++)
				{
					if (j != k)
					{
						if (neighb1[j] > -1 && neighb1[k] > -1)
						{
							real r = 
								(dr.getValueX(neighb1[j]) - dr.getValueX(neighb1[k])) * (dr.getValueX(neighb1[j]) - dr.getValueX(neighb1[k])) + 
								(dr.getValueY(neighb1[j]) - dr.getValueY(neighb1[k])) * (dr.getValueY(neighb1[j]) - dr.getValueY(neighb1[k])) +
								(dr.getValueZ(neighb1[j]) - dr.getValueZ(neighb1[k])) * (dr.getValueZ(neighb1[j]) - dr.getValueZ(neighb1[k]));
							if (r <= (rCut + 2 * rView))
							{
								verle1.push_back(neighb1[k]);
								a += 1;	
								ind = 1;
							}
						}
					}
				}
				if (ind == 1) cutVerle1.push_back(a);
				ind = 0;
			}
		}
		cutVerle1.pop_back(); 
		
		for (int i = 0; i < (cutVerle1.size() - 1); i++)
		{
			int j;
			Fx=0.; Fy=0.; Fz=0.;
			int a = 0;
			for (int k = cutVerle1[i]; k < cutVerle1[i + 1]; k++)
			{
				j = cutVerle1[i];
				if (verle1[j] != verle1[k])
				{
					if (verle1[j] != -1 && verle1[k] != -1)
					{
						real r2 = 
								(dr.getValueX(verle1[j]) - dr.getValueX(verle1[k])) * //verle[cutVerle[i]]
								(dr.getValueX(verle1[j]) - dr.getValueX(verle1[k])) + 
								(dr.getValueY(verle1[j]) - dr.getValueY(verle1[k])) * 
								(dr.getValueY(verle1[j]) - dr.getValueY(verle1[k])) +
								(dr.getValueZ(verle1[j]) - dr.getValueZ(verle1[k])) * 
								(dr.getValueZ(verle1[j]) - dr.getValueZ(verle1[k]));


						XX = rVr1 * r2 - X2min;

						FV = rVr2 / r2;

						Ur = FV * FV * FV;

						real RR = aVr * ( Ur * (Ur - 1.0) - (aLJ3 + bLJ2 * XX) * XX * XX); //V(r)/2     какой размерности массив? RR(k)

						real FR = aVr2 * ( Ur * (Ur - 0.5) * FV + (aLJ+ bLJ * XX) * XX ); //V'(r)/r СИЛА НА ЭР(межатомное расстояние)

						Fx += FR * (dr.getValueX(verle1[j]) - dr.getValueX(verle1[k])) * 
							(dr.getValueX(verle1[j]) - dr.getValueX(verle1[k]));
						Fy += FR * (dr.getValueY(verle1[j]) - dr.getValueY(verle1[k])) * 
								(dr.getValueY(verle1[j]) - dr.getValueY(verle1[k]));
						Fz += FR * (dr.getValueZ(verle1[j]) - dr.getValueZ(verle1[k])) * 
								(dr.getValueZ(verle1[j]) - dr.getValueZ(verle1[k])); 

					}
				}		
			}
			
			Accel.setValueX(verle1[j], (Fx/m1));
			Accel.setValueY(verle1[j], (Fy/m1));
			Accel.setValueZ(verle1[j], (Fz/m1));
		}
		cutNeighb1.clear();
		cutVerle1.clear();
		verle1.clear();
		neighb1.clear();
		delete[] Addr1;
		delete[] cell1;*/

		
		

//   СНАЧАЛА СЧИТАЮ СИЛЫ ВЗАИМОДЕЙСТВИЯ МЕЖДУ АТОМАМИ И ВЫСЧИТЫВАЮ УСКОРЕНИЯ (УСКОРЕНИЯ НАЧАЛЬНЫЕ НЕ НУЛЕВЫЕ!!!)
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
					if (r2 < rCut * rCut) //выводить тех, с кем взаимодействует
					{							
						XX = rVr1 * r2 - X2min;
						FV = rVr2 / r2;
						Ur = FV * FV * FV;

						real RR = aVr * ( Ur * (Ur - 1.0) - (aLJ3 + bLJ2 * XX) * XX * XX); //V(r)/2     какой размерности массив? RR(k)
						real FR = aVr2 * ( Ur * (Ur - 0.5) * FV + (aLJ + bLJ * XX) * XX ); //V'(r)/r СИЛА НА ЭР(межатомное расстояние)

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

			//новые координаты, используя старые скорости и ускорения
			/*for (int i = 0; i < NGen; i++) 
			{
				real xi = Vel.getValueX(i) * dt + 0.5 * Accel.getValueX(i) * dt * dt;
				xi += dr.getValueX(i);
				prd.X(xi, Space.Lx, 0, Space.Lx);
				dr.setValueX(i, xi);
			

				real yi = Vel.getValueY(i) * dt + 0.5 * Accel.getValueY(i) * dt * dt;
				yi += dr.getValueY(i);
				prd.Y(yi, Space.Ly, 0, Space.Ly);
				dr.setValueY(i, yi);
				
								
				real zi = Vel.getValueZ(i) * dt + 0.5 * Accel.getValueZ(i) * dt * dt;
				zi += dr.getValueZ(i);
				prd.Z(zi, Space.Lz, 0, Space.Lz);
				dr.setValueZ(i, zi);
			}
			
			//скорости на половине шага с использованием старых скоростей и ускорений
			for (int i = 0; i < NGen; i++)
			{
				double val = Vel.getValueX(i) + 0.5 * Accel.getValueX(i) * dt;
				halfVel.setValueX(i, val);
					val = Vel.getValueY(i) + 0.5 * Accel.getValueY(i) * dt;
				halfVel.setValueY(i, val);
					val = Vel.getValueZ(i) + 0.5 * Accel.getValueZ(i) * dt;
				halfVel.setValueZ(i, val);
			}*/


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
			
			//перестроение списка Верле
			/*if (timestep % 10 == 0)
			{
				if (timestep != 0)
				{
					cutNeighb.clear();
					cutVerle.clear();
					verle.clear();
					notverle.clear();
				}

				bool flag;
				III = 0;
				ind = 0;

				for (int k = 0; k < NzBox; k++)
				{
					for (int j = 0; j < NyBox; j++)
					{
						for (int i = 0; i < NxBox; i++)
						{
							flag = false; 
							int i1 = 0, i2 = 0;

							aa += 1;

							for (int n = 0; n < NGen; n++) //сделать всё одним циклом по NGen, внутри проверка по индексам, координату делить на длину(lx), получается число =>
								//это есть координата i j k, всё, что далее будет с подобными появляться - отдавать в связанный список
							{
								if (dr.getValueX(n) >= (lx * i) && dr.getValueX(n) <= (lx * (i + 1)) &&
									dr.getValueY(n) >= (ly * j) && dr.getValueY(n) <= (ly * (j + 1)) &&
									dr.getValueZ(n) >= (lz * k) && dr.getValueZ(n) <= (lz * (k + 1)))
								{
									cell[i1] = n;
									i1 += 1;									
								}
							}
							
							if (i1 > 0)
							{
								flag = true;
								i2 = i1 - 1;
							}

							if (flag == true)
							{
								qsort(cell, i2, sizeof(int), Comp1);
								
								Addr[k][j][i] = cell[0];
								//neighb1.push_back(cell[0]);
								//neighb[cell[0]] = Addr[k][j][i];
								for (i1 = 1; i1 < i2; i1++)
								{	
									neighb[cell[i1-1]] = cell[i1];
									//forces << neighb[cell[i1-1]] << endl;
									//neighb1.push_back(cell[i1]);
									ind += 1;
								}

								neighb[cell[i2]] = -1;
								//neighb1.push_back(-1);
								cutNeighb.push_back(ind);
								//cout << ind << endl;
								//forces << neighb[cell[i2]] << endl;
								ind += 1;
							}
						
							else 
							{
								Addr[k][j][i] = -1;
							}
							//cout << Addr[k][j][i] << "   " << aa << endl; //ТУТ МОЖНО ПОСМОТРЕТЬ ВЕСЬ АДРЕСНЫЙ СПИСОК
						}
					}
				}//адресный список 
				//cout << neighb1.size() << endl;
				//cout << cutNeighb1.size() << endl;
				//cout << NzBox << "  " << NyBox << "   " << NxBox << endl;
				for (int i = 0; i < neighb.size(); i++)
				{
					
					//if (neighb[i] == 0) neighb[i] = -1;
					//forces << neighb[i] << endl;
				}

				int atm = 0;
				int ind1 = 0;
				for (int k = 0; k < NzBox; k++)
				{
					for (int j = 0; j < NyBox; j++)
					{
						for (int i = 0; i < NxBox; i++)
						{	
							//ind1 = 0;
							int b, d;
							n = Addr[k][j][i];
							if (n != -1) b = neighb[n];
							
							while (n != -1) //бесконечный цикл иногда
							{
								while (b != -1)
								{
									if (n != b)
									{
										real r2 = 
										(dr.getValueX(n) - dr.getValueX(b)) * 
										(dr.getValueX(n) - dr.getValueX(b)) + 
										(dr.getValueY(n) - dr.getValueY(b)) * 
										(dr.getValueY(n) - dr.getValueY(b)) +
										(dr.getValueZ(n) - dr.getValueZ(b)) * 
										(dr.getValueZ(n) - dr.getValueZ(b));
									
										if (r2 < rCut * rCut)
										{
											forces << b << endl;
											verle.push_back(b);
										}

										b = neighb[b];
									}
									
								}
								
								for (int k1 = (k - 1); k1 <= (k + 1); k1++)
								{
									int k2 = k1;
									if (k1 == -1) k2 += NzBox;
									if (k1 == (NzBox)) k2 -= NzBox;
				
									for (int g1 = (j - 1); g1 <= (j + 1); g1++)
									{
										int g2 = g1;
										if (g1 == -1) g2 += NyBox;
										if (g1 == (NyBox)) g2 -= NyBox;
						
										for (int i1 = (i - 1); i1 <= (i + 1); i1++)
										{
											int i2 = i1;
											if (i1 == -1) i2 += NxBox;
											if (i1 == (NxBox)) i2 -= NxBox;
											
											int c = Addr[k2][g2][i2];
											
											if (c != -1) d = neighb[c];
							
											while (d != -1)// && (Addr[k][j][i] != Addr[k2][g2][i2]))
											{
												if (n != d)
												{
													real r2 = 
													(dr.getValueX(n) - dr.getValueX(d)) * 
													(dr.getValueX(n) - dr.getValueX(d)) + 
													(dr.getValueY(n) - dr.getValueY(d)) * 
													(dr.getValueY(n) - dr.getValueY(d)) +
													(dr.getValueZ(n) - dr.getValueZ(d)) * 
													(dr.getValueZ(n) - dr.getValueZ(d));

													if (r2 < rCut * rCut)
													{
														forces << d << endl;
														verle.push_back(d);
													}

													d = neighb[d];
													
												}
											}
											
											
										}
									}
								}
							n = neighb[n];
							}
							verle.push_back(-1);
						}
					}
				}





				/*ind = 0;
				int a = 0;
				for (int i = 0; i < (cutNeighb.size() - 1); i++)
				{
					for (int j = cutNeighb[i]; j < cutNeighb[i+1]; j++)
					{
						for (int k = cutNeighb[i]; k < cutNeighb[i+1]; k++)
						{
							if (j != k) 
							{
								if (neighb[j] != -1 && neighb[k] != -1)
								{
									real r = 
									(dr.getValueX(j) - dr.getValueX(k)) * (dr.getValueX(j) - dr.getValueX(k)) + 
									(dr.getValueY(j) - dr.getValueY(k)) * (dr.getValueY(j) - dr.getValueY(k)) +
									(dr.getValueZ(j) - dr.getValueZ(k)) * (dr.getValueZ(j) - dr.getValueZ(k));
									//cout << r << endl;
									if (r < (rCut + 2 * rView))
									{
										verle.push_back(k);
										a += 1;
										ind = 1;
										//forces << k << "  " << i << "  " << neighb[k] << endl;
									}
								}
							}
						}
						if (ind == 1) { cutVerle.push_back(a); }//forces << a << "  " << i << endl; }
						//cout << j << endl;
						else if (ind != 1 && neighb1[j] > -1) { notverle.push_back(neighb1[j]); }
						ind = 0;

						
					}
				}
				cout << verle.size() << endl;
				//cout << cutVerle1.size() << endl;
			}*/

			//упрощенный список Верле
			/*if (timestep % 1 == 0)
			{
				if (timestep != 0)
				{
					verle.clear();
					cutVerle.clear();
				}
				double rview = 2 * rCut * rCut;
				int a = 0;

				for (int i = 0; i < NGen; i++)
				{
					verle.push_back(i);
					a += 1;
					for (int j = 0; j < NGen; j++)
					{
						if (i != j)
						{verle.push_back(j); a += 1;
							//real r2 = 
							//		(dr.getValueX(i) - dr.getValueX(j)) * 
							//		(dr.getValueX(i) - dr.getValueX(j)) + 
							//		(dr.getValueY(i) - dr.getValueY(j)) * 
							//		(dr.getValueY(i) - dr.getValueY(j)) +
							//		(dr.getValueZ(i) - dr.getValueZ(j)) * 
							//		(dr.getValueZ(i) - dr.getValueZ(j));
							//if (r2 < rview)
							//{
							//	verle.push_back(j);
							//	a += 1;
							//	//forces << j << endl;
							//}
						}

						
					}
					//verle.push_back(-1);
					//a += 1;
					cutVerle.push_back(a);
				}
				//cutVerle.push_back(verle.size());
			}*/

			//cout << verle.size() << endl;
			
			

			/*for (int i = 0; i < (NGen - 1); i++)
			{
				int k = cutVerle[i];
				int a = verle[k];
				Fx=0.; Fy=0.; Fz=0.; double u = 0.;
				for (int j = cutVerle[i]; j < cutVerle[i + 1]; j++)
				{
					int b = verle[j];
					if (b != -1)
					{
						if (a != b)
						{//forces << verle[j] << endl;
							real r2 = 
										(dr.getValueX(a) - dr.getValueX(b)) * 
										(dr.getValueX(a) - dr.getValueX(b)) + 
										(dr.getValueY(a) - dr.getValueY(b)) * 
										(dr.getValueY(a) - dr.getValueY(b)) +
										(dr.getValueZ(a) - dr.getValueZ(b)) * 
										(dr.getValueZ(a) - dr.getValueZ(b));

							if (r2 < rCut * rCut)
							{
								XX = rVr1 * r2 - X2min;
								FV = rVr2 / r2;
								Ur = FV * FV * FV;

								real RR = aVr * ( Ur * (Ur - 1.0) - (aLJ3 + bLJ2 * XX) * XX * XX); //V(r)/2 
								real FR = aVr2 * ( Ur * (Ur - 0.5) * FV + (aLJ+ bLJ * XX) * XX ); //V'(r)/r СИЛА НА ЭР(межатомное расстояние)

								Fx += FR * (dr.getValueX(b) - dr.getValueX(a));								
								Fy += FR * (dr.getValueY(b) - dr.getValueY(a));							
								Fz += FR * (dr.getValueZ(b) - dr.getValueZ(a));

								u += RR;
							}
						}
					}
				}

				Accel.setValueX(a, (Fx/m1));
				Accel.setValueY(a, (Fy/m1));
				Accel.setValueZ(a, (Fz/m1));
				U[a] = u;

			}*/

			//расчет сил ЛД
			
			
			
			/*for (int i = 0; i < (cutVerle.size() - 1); i++)
			{
				for (int j = cutVerle[i]; j < cutVerle[i+1]; j++)
				{
					Fx=0.; Fy=0.; Fz=0.; int a = neighb1[verle[j]]; //cout << a << endl;
					for (int k = cutVerle[i]; k < cutVerle[i+1]; k++)
					{
						int b = neighb1[verle[k]]; cout << b << endl;
						if (j != k) 
						{
							
							if (a != -1 && b != -1)
							{
							real r2 = 
									(dr.getValueX(a) - dr.getValueX(b)) * //(dr.getValueX(verle[j]) - dr.getValueX(verle[k]))
									(dr.getValueX(a) - dr.getValueX(b)) + 
									(dr.getValueY(a) - dr.getValueY(b)) * 
									(dr.getValueY(a) - dr.getValueY(b)) +
									(dr.getValueZ(a) - dr.getValueZ(b)) * 
									(dr.getValueZ(a) - dr.getValueZ(b));

								if (r2 < rCut * rCut)
								{
									XX = rVr1 * r2 - X2min;

									FV = rVr2 / r2;

									Ur = FV * FV * FV;

									real RR = aVr * ( Ur * (Ur - 1.0) - (aLJ3 + bLJ2 * XX) * XX * XX); //V(r)/2     какой размерности массив? RR(k)

									real FR = aVr2 * ( Ur * (Ur - 0.5) * FV + (aLJ+ bLJ * XX) * XX ); //V'(r)/r СИЛА НА ЭР(межатомное расстояние)

									Fx += FR * (dr.getValueX(a) - dr.getValueX(b)) * 
										(dr.getValueX(a) - dr.getValueX(b));
									Fy += FR * (dr.getValueY(a) - dr.getValueY(b)) * 
											(dr.getValueY(a) - dr.getValueY(b));
									Fz += FR * (dr.getValueZ(a) - dr.getValueZ(b)) * 
											(dr.getValueZ(a) - dr.getValueZ(b)); 

									U[a] += RR;

								}	
							}

						}

					}*/

			for (int i = 0; i < NGen; i++)
			{
			
			//for (int i = 0; i < (cutVerle.size() - 1); i++)
			//{
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

							real FR = aVr2 * ( Ur * (Ur - 0.5) * FV + (aLJ+ bLJ * XX) * XX ); //V'(r)/r СИЛА НА ЭР(межатомное расстояние)

							Fx += FR * (dr.getValueX(j) - dr.getValueX(i));
								
							Fy += FR * (dr.getValueY(j) - dr.getValueY(i));
							
							Fz += FR * (dr.getValueZ(j) - dr.getValueZ(i));

							u += RR;
						}
					}
				/*int j;
				Fx=0.; Fy=0.; Fz=0.;
				for (int k = cutVerle[i]; k < cutVerle[i + 1]; k++)
				{
					j = cutVerle[i];
					if (verle[j] != verle[k])
					{
						if (verle[j] != -1 && verle[k] != -1)
						{
							real r2 = 
									(dr.getValueX(verle[j]) - dr.getValueX(verle[k])) * //verle[cutVerle[i]]
									(dr.getValueX(verle[j]) - dr.getValueX(verle[k])) + 
									(dr.getValueY(verle[j]) - dr.getValueY(verle[k])) * 
									(dr.getValueY(verle[j]) - dr.getValueY(verle[k])) +
									(dr.getValueZ(verle[j]) - dr.getValueZ(verle[k])) * 
									(dr.getValueZ(verle[j]) - dr.getValueZ(verle[k]));

							if (r2 < rCut * rCut)
							{
								XX = rVr1 * r2 - X2min;

								FV = rVr2 / r2;

								Ur = FV * FV * FV;

								real RR = aVr * ( Ur * (Ur - 1.0) - (aLJ3 + bLJ2 * XX) * XX * XX); //V(r)/2     какой размерности массив? RR(k)

								real FR = aVr2 * ( Ur * (Ur - 0.5) * FV + (aLJ+ bLJ * XX) * XX ); //V'(r)/r СИЛА НА ЭР(межатомное расстояние)

								Fx += FR * (dr.getValueX(verle[k]) - dr.getValueX(verle[j]));
									
								Fy += FR * (dr.getValueY(verle[k]) - dr.getValueY(verle[j]));

								Fz += FR * (dr.getValueZ(verle[k]) - dr.getValueZ(verle[j]));

								u += RR;
								//U[verle[j]] += RR;
							}
							

						}
					}*/		
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

			/*if (ActivateThermLan == 1) 
			{
				for (int j = 0; j < notverle.size(); j++)
				{
					double Fx = LanForce.getValueX(notverle[j]); 
					double Fy = LanForce.getValueY(notverle[j]); 
					double Fz = LanForce.getValueZ(notverle[j]);
					Accel.setValueX(notverle[j], (Fx/m1));
					Accel.setValueY(notverle[j], (Fy/m1));
					Accel.setValueZ(notverle[j], (Fz/m1));
					//ForX += Fx; ForY += Fy; ForZ += Fz;
				}
				
			}*/

			//
			//
			//скорости
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
			
			//постановка задачи, графики температуры, расчет потенциала, методика ускорения расчета, количество времени на расчет с и без списка Верле


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




	

	/*else
	{

		Space.delta_x = Space.Lx / num_proc;
		Space.delta_y = Space.Ly / num_proc;
		Space.delta_z = Space.Lz / num_proc;
		

		if (timestep == 0) //либо через Ncall = 0
		{
			Cent.setValueX(my_num, (Space.delta_x * my_num)); //центры Вороного, пока через центр масс
			Cent.setValueY(my_num, (Space.delta_y * my_num)); 
			Cent.setValueZ(my_num, (Space.delta_z * my_num));

			NGen = 0;

			CalculateNGen(Space, Consts, Cent.getValueX(my_num));

			MyN = NGen / 1.75;

			

		}
		

		else
		{
			//дополнить пересчет центра масс при обновлении
		}

		NpartW = MyN;
		
		
		NullPtrsVecR dr(NGen), Vel(NGen), Accel(MyN);

		NullPtrsInt EnergyPot(MyN), KindA(NGen);

		OperatingWithNeibs OperNeibs(my_num, num_proc, rView, Space.Lx, Space.Ly, Space.Lz);

		NullPtrsInt NeibsNum(num_proc), ContactM(num_proc), NCHKL(num_proc), MCHKL(num_proc);

		//NullPtrsInt NCHKL(num_proc), MCHKL(num_proc);
		vector<NatTo> NatToAll, NatToMe;

		InitValuesCoord(Space, Consts, Cent.getValueX(my_num), dr);
		InitValuesOthers(NGen, MyN, Vel, Accel, EnergyPot, KindA);

		rMax.setValue(my_num, RMax(MyN, dr));
		
		//SendRecvRMax(rMax, ReqR, ReqS);
		

		//передача соседним процессам максимальных радиусов сфер
		real r_max_recv, r_max_send; 
		r_max_send = rMax.getValue(my_num);
		for (i = 1; i < num_proc; i++) 
		{			
			real r_max_recv, r_max_send;
			r_max_send = rMax.getValue(my_num);
			if (i != my_num)
			{
				MPI_Irecv(&r_max_recv, 1, MPI_DOUBLE, i, 4000, MPI_COMM_WORLD, &ReqR[i]);// 4000 - радиусы
				MPI_Isend(&r_max_send, 1, MPI_DOUBLE, i, 4000, MPI_COMM_WORLD, &ReqS[i]);
				MPI_Wait(&ReqR[i], &status);
				MPI_Wait(&ReqS[i], &status);
				rMax.setValue(i, r_max_recv);
			}			
		}


		//передача центров другим процессам
		real Xsys_send, Ysys_send, Zsys_send;
		real Xsys_recv, Ysys_recv, Zsys_recv;

		Xsys_send = Cent.getValueX(my_num); Ysys_send = Cent.getValueY(my_num); Zsys_send = Cent.getValueZ(my_num);

		for (i = 1; i < num_proc; i++) 
		{	
			if (i != my_num)
			{
				MPI_Irecv(&Xsys_recv, 1, MPI_DOUBLE, i, 1000, MPI_COMM_WORLD, &ReqR[i]);// 1000 - центры по x
				MPI_Isend(&Xsys_send, 1, MPI_DOUBLE, i, 1000, MPI_COMM_WORLD, &ReqS[i]);
				MPI_Wait(&ReqR[i], &status);
				MPI_Wait(&ReqS[i], &status);
				Cent.setValueX(i, Xsys_recv);

				MPI_Irecv(&Ysys_recv, 1, MPI_DOUBLE, i, 2000, MPI_COMM_WORLD, &ReqR[i]);// 2000 - центры по y
				MPI_Isend(&Ysys_send, 1, MPI_DOUBLE, i, 2000, MPI_COMM_WORLD, &ReqS[i]);
				MPI_Wait(&ReqR[i], &status);
				MPI_Wait(&ReqS[i], &status);
				Cent.setValueY(i, Ysys_recv);

				MPI_Irecv(&Zsys_recv, 1, MPI_DOUBLE, i, 3000, MPI_COMM_WORLD, &ReqR[i]);// 3000 - центры по z
				MPI_Isend(&Zsys_send, 1, MPI_DOUBLE, i, 3000, MPI_COMM_WORLD, &ReqS[i]);
				MPI_Wait(&ReqR[i], &status);
				MPI_Wait(&ReqS[i], &status);
				Cent.setValueZ(i, Zsys_recv);
			}
		}

		//считаем число и список соседей Вороного
		

		NeibsM = OperNeibs.SearchNeibAmount(rMax, Cent, prd);

		*OperNeibs.SearchNeibNumb(rMax, Cent, prd, NeibsNum);
		
		

		//NatTo *NatToAll = nullptr; NatTo *NatToMe = nullptr;
		
		CheckList(NeibsM, MCHKL, NeibsNum, ReqS, ReqR, status);


		if (Neibs0 < NeibsM || Neibs0 > (NeibsM*3+16))
		{
			Neibs0 = (NeibsM*7)/3;

			//if(NatToAll)
			//{
				
				//delete[] NatToAll;
				//delete[] NatToMe;
			//}
			NatToAll.resize(Neibs0); //= new NatTo [Neibs0];
			NatToMe.resize(Neibs0); //= new NatTo [Neibs0];
		}
		
		

		for (i = 0; i < NeibsM; i++)
		{
			int ii = NeibsNum.getValue(i);
			//if (ii != my_num)
			//{
				for (j = 0; j < NeibsM; j++)
				{
					int jj = NeibsNum.getValue(j);
					NatToAll[i].MpPair(NatToAll, ii, jj, Cent, Space);
				}
			//}
		}	


		LookingNeibs1(NeibsNum, ContactM, NCHKL, KindA, rMax, dr, Space, Cent, prd, NatToAll);


		LookingNeibs2(NeibsNum, ContactM, MCHKL, NCHKL, ReqS, ReqR, status, Bufer);


		CheckList(Neibs, MCHKL, NeibsNum, ReqS, ReqR, status);


		if (Neibs > 0)
		{
			AtomExchange(NatToAll, NatToMe, Space, prd, NeibsNum, KindA, ContactM, Cent, dr, Vel, status, ReqS, ReqR);
		}

		
		//место для функции балансировки


		//начинать со строки ...





		//cout << rMax.getValue(my_num) << endl;
		
		//if (my_num == 2) cout << rMax.getValue(1) << endl;




	}*/
	


	int workTime = clock();


	cout << ((float)workTime) / CLOCKS_PER_SEC << endl;



	//cout << ((float)workTime) / CLOCKS_PER_SEC << endl;
	MPI_Finalize();
	return 0;
	getchar();
}




//taskkill /f /im md_pshe_pshe_one_oneproc.exe










	//	ind = 0; a = 0;
			//	int ind_old = 0; 
			//	cutNeighb.push_back(0); cutVerle.push_back(0);

			//	for (int i = 0; i < NGen; i++)//начало списка Верле (общий)
			//	{
			//		if (neighb[i] == -1) cutNeighb.push_back(i); //cout << i << endl;}
			//		//forces << neighb[i] << "  " << i << endl;
			//	}
			//	
			//	
			//	
			//	for (int i = 0; i < (cutNeighb.size() - 1); i++)
			//	{
			//		for (int j = cutNeighb[i]; j < cutNeighb[i + 1]; j++)
			//		{
			//			for (int k = cutNeighb[i]; k < cutNeighb[i + 1]; k++)
			//			{
			//				if (j != k)
			//				{
			//					if (neighb[j] > -1 && neighb[k] > -1)
			//					{
			//						real r = 
			//							(dr.getValueX(neighb[j]) - dr.getValueX(neighb[k])) * (dr.getValueX(neighb[j]) - dr.getValueX(neighb[k])) + 
			//							(dr.getValueY(neighb[j]) - dr.getValueY(neighb[k])) * (dr.getValueY(neighb[j]) - dr.getValueY(neighb[k])) +
			//							(dr.getValueZ(neighb[j]) - dr.getValueZ(neighb[k])) * (dr.getValueZ(neighb[j]) - dr.getValueZ(neighb[k]));
			//						if (r < (rCut + 2 * rView))
			//						{
			//							verle.push_back(neighb[k]);
			//							//forces << neighb[j] << "  " << i << endl;
			//							//verle[a] = neighb[k];
			//							//cout << verle[a] << endl;
			//							a += 1;	
			//							ind = 1;
			//						}
			//					}
			//				}
			//			}
			//			//verle.push_back(-1);
			//			if (ind == 1) { cutVerle.push_back(a); }//a += 1; }							
			//			else if (ind != 1 && neighb[j] > -1) {notverle.push_back(neighb[j]);}//a += 1; cutVerle.push_back(a);} a+1
			//			//a += 1;
			//			
			//			ind = 0;
			//		}
			//	}
			//	cutVerle.pop_back();
			//}//конец timestep для перестроения списка
			////cout << verle.size() << "   verle.size" << endl; //общий список из всех частиц (список Верле)
			////cout << cutVerle.size() << endl; //массив из номеров, с которых начинается новая частица в списке Верле (размер меньше NGen)
			////cout << cutNeighb.size() << endl; //массив по ячейкам
			
			
			/*//расчет сил ЛД
			double ForX = 0, ForY = 0, ForZ = 0;
			//cout << timestep << endl;
			for (int i = 0; i < (cutVerle.size() - 1); i++)
			{
				int j;
				Fx=0.; Fy=0.; Fz=0.;
				int a = 0;
				for (int k = cutVerle[i]; k < cutVerle[i + 1]; k++)
				{
					j = cutVerle[i];
					if (verle[j] != verle[k])
					{
						if (verle[j] != -1 && verle[k] != -1)
						{
							real r2 = 
									(dr.getValueX(verle[j]) - dr.getValueX(verle[k])) * //verle[cutVerle[i]]
									(dr.getValueX(verle[j]) - dr.getValueX(verle[k])) + 
									(dr.getValueY(verle[j]) - dr.getValueY(verle[k])) * 
									(dr.getValueY(verle[j]) - dr.getValueY(verle[k])) +
									(dr.getValueZ(verle[j]) - dr.getValueZ(verle[k])) * 
									(dr.getValueZ(verle[j]) - dr.getValueZ(verle[k]));

							if (r2 < rCut * rCut)
							{
								XX = rVr1 * r2 - X2min;

								FV = rVr2 / r2;

								Ur = FV * FV * FV;

								real RR = aVr * ( Ur * (Ur - 1.0) - (aLJ3 + bLJ2 * XX) * XX * XX); //V(r)/2     какой размерности массив? RR(k)

								real FR = aVr2 * ( Ur * (Ur - 0.5) * FV + (aLJ+ bLJ * XX) * XX ); //V'(r)/r СИЛА НА ЭР(межатомное расстояние)

								Fx += FR * (dr.getValueX(verle[j]) - dr.getValueX(verle[k])) * 
									(dr.getValueX(verle[j]) - dr.getValueX(verle[k]));
								Fy += FR * (dr.getValueY(verle[j]) - dr.getValueY(verle[k])) * 
										(dr.getValueY(verle[j]) - dr.getValueY(verle[k]));
								Fz += FR * (dr.getValueZ(verle[j]) - dr.getValueZ(verle[k])) * 
										(dr.getValueZ(verle[j]) - dr.getValueZ(verle[k])); 

								U[verle[j]] += RR;
							}
							

						}
					}		
				}

				if (ActivateThermLan == 1) 
				{
					Fx += LanForce.getValueX(verle[j]); 
					Fy += LanForce.getValueY(verle[j]); 
					Fz += LanForce.getValueZ(verle[j]);
				}
				ForX += Fx; ForY += Fy; ForZ += Fz;
				Accel.setValueX(verle[j], (Fx/m1));
				Accel.setValueY(verle[j], (Fy/m1));
				Accel.setValueZ(verle[j], (Fz/m1));
			}*/