#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <vector>
#include <random>


#include "mpi.h"
#include "defines_calc.h"
#include "defines_vars_and_structs.h"
#include "my_func.h"



int Comp1(const void * a, const void * b)
{
	return (*(int*)a - *(int*)b);
}
int Round (real value)
{
   return (int) floor(value + 0.5);
}


void CalculateNGen(CalcSpace Space, CalcConstsLD Consts, real Cent)
{
	for (Space.i1 = -Space.n1; Space.i1 < Space.n1; Space.i1++)
	{
		for (Space.i2 = -Space.n2; Space.i2 < Space.n2; Space.i2++)
		{
			for (Space.i3 = -Space.n3; Space.i3 < Space.n3; Space.i3++)
			{				
				Space.dx = 0.5 * Consts.a * Space.i1 + 0.5 * Consts.a * Space.i3;
				Space.dy = 0.5 * Consts.a * Space.i1 + 0.5 * Consts.a * Space.i2;
				Space.dz = 0.5 * Consts.a * Space.i2 + 0.5 * Consts.a * Space.i3;
				
				if (Space.dx >= (Cent - Space.delta_x) && Space.dx <= (Cent + Space.delta_x) &&
					Space.dy >= 0 && Space.dy <= Space.Ly &&
					Space.dz >= 0 && Space.dz <= Space.Lz) NGen +=1;

			}
		}
	}
	NGen = NGen*1.75;
}


void InitValuesCoord(CalcSpace Space, CalcConstsLD Consts, real Cent, NullPtrsVecR &dr)
{
	int i = 0;
	for (Space.i1 = -Space.n1; Space.i1 < Space.n1; Space.i1++)
	{
		for (Space.i2 = -Space.n2; Space.i2 < Space.n2; Space.i2++)
		{
			for (Space.i3 = -Space.n3; Space.i3 < Space.n3; Space.i3++)
			{				
				Space.dx = 0.5 * Consts.a * Space.i1 + 0.5 * Consts.a * Space.i3;
				Space.dy = 0.5 * Consts.a * Space.i1 + 0.5 * Consts.a * Space.i2;
				Space.dz = 0.5 * Consts.a * Space.i2 + 0.5 * Consts.a * Space.i3;
				
				if (Space.dx >= (Cent - Space.delta_x) && Space.dx <= (Cent + Space.delta_x) &&
					Space.dy >= 0 && Space.dy <= Space.Ly &&
					Space.dz >= 0 && Space.dz <= Space.Lz)
				{
					dr.setValueX(i, Space.dx);
					dr.setValueY(i, Space.dy);
					dr.setValueZ(i, Space.dz);
					i +=1;
				}
			}
		}
	}
}


void InitValuesOthers(int NGen, int MyN, NullPtrsVecR &Vel, NullPtrsVecR &Accel, NullPtrsInt &EnergyPot, NullPtrsInt &KindA, NullPtrsVecR &LanForce)
{
	int i;
	for (i = 0; i < NGen; i++)
	{
		Vel.setValueX(i, 0);
		Vel.setValueY(i, 0);
		Vel.setValueZ(i, 0);

		LanForce.setValueX(i, 0);
		LanForce.setValueY(i, 0);
		LanForce.setValueZ(i, 0);
	}
	
	for (i = 0; i < MyN; i++)
	{
		Accel.setValueX(i, 0);
		Accel.setValueY(i, 0);
		Accel.setValueZ(i, 0);
		EnergyPot.setValue(i, 0);
		KindA.setValue(i, 1);
	}
}


real RMax(int MyN, NullPtrsVecR &dr) // поиск максимального радиуса для построения сферы для грубого списка соседей вороных
{
	int i, k;
	real r;
	real r_max = 0;
	for (i = 0; i < MyN; i++)
	{
		for (k = 0; k < MyN; k++)
		{
			if (k != i)
			{	//r=exp(0.5*log(((dr[tid].x-dr[k].x)*(dr[tid].x-dr[k].x)+(dr[tid].y-dr[k].y)*(dr[tid].y-dr[k].y)+(dr[tid].z-dr[k].z)*(dr[tid].z-dr[k].z))));
				//r = exp(0.5 * log(( Sqr(dr[tid].x - dr[k].x) + Sqr(dr[tid].y - dr[k].y) + Sqr(dr[tid].z - dr[k].z))));
				r = exp(0.5 * log((dr.getValueX(i) - dr.getValueX(k)) * (dr.getValueX(i) - dr.getValueX(k)) + 
									(dr.getValueY(i) - dr.getValueY(k)) * (dr.getValueY(i) - dr.getValueY(k)) + 
										(dr.getValueZ(i) - dr.getValueZ(k)) * (dr.getValueZ(i) - dr.getValueZ(k))));

				if (abs(r) > abs(r_max)) r_max = r;
			}
		}
	}
	return r_max;
}


void LookingNeibs1(NullPtrsInt &NeibsNum, NullPtrsInt &ContactM, NullPtrsInt &NCHKL, NullPtrsInt &KindA, NullPtrsReal &rMax, NullPtrsVecR &dr, CalcSpace Space, NullPtrsVecR &Cent, Period &prd,
				   std::vector<NatTo>&NatToAll
				  /*NatTo NatToAll[]*/)
{
	int i;
	for (i = 0; i < NeibsM; i++) //рассматриваем атомы в приграничной области
	{
		ContactM.setValue(i, 0);
	}

	if (Ncall != 0)
	{
		for (i = 0; i < NeibsM; i++)
		{
			NCHKL.setValue(i, k);
		}
	}
	else
	{
		Ncall = 1;
		for (i = 0; i < NeibsM; i++)
		{
			int ii = NeibsNum.getValue(i);
			real bb = rMax.getValue(my_num) + rMax.getValue(ii) + rView*6; //поменять

			k = 1 + std::min(1, Round(bb/Space.Lx));
			k = (1 + std::min(1, Round(bb/Space.Ly))) * k;
			k = (1 + std::min(1, Round(bb/Space.Lz))) * k;

			NCHKL.setValue(i, k);
		}
	}

	
	if (NeibsM > 0)//&& NpartW > 0
	{
		real VDmin = Space.Lx*Space.Lx + Space.Ly*Space.Ly + Space.Lz*Space.Lz;
		real DL = 0;
		for (i = 0; i < NeibsM; i++)
		{
			int ii = NeibsNum.getValue(i);
			if (ii != my_num)
			{
				DL = prd.Xa(Cent.getValueX(ii), Cent.getValueX(my_num), Space.Lx, 0, Space.Lx) 
					* prd.Xa(Cent.getValueX(ii), Cent.getValueX(my_num), Space.Lx, 0, Space.Lx);
				DL = prd.Ya(Cent.getValueY(ii), Cent.getValueY(my_num), Space.Ly, 0, Space.Ly) 
					* prd.Ya(Cent.getValueY(ii), Cent.getValueY(my_num), Space.Ly, 0, Space.Ly) + DL;
				DL = prd.Za(Cent.getValueZ(ii), Cent.getValueZ(my_num), Space.Lz, 0, Space.Lz) 
					* prd.Za(Cent.getValueZ(ii), Cent.getValueZ(my_num), Space.Lz, 0, Space.Lz) + DL;
			}
			if (VDmin > DL) VDmin = DL;
		}
		
		VDmin = 0.5 * sqrt(VDmin) - 6*rView;

		if (VDmin > 0) VDmin = VDmin * VDmin;
		else VDmin = 0;
		//NpartW = NGen; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		for (i = 0; i < NpartW; i++)
		{
			if (KindA.getValue(i) > 0)
			{
				Xi = prd.Xa(prd.X(dr.getValueX(i), Space.Lx, 0, Space.Lx), Cent.getValueX(my_num), Space.Lx, 0, Space.Lx);
				DL = Xi * Xi;
				Yi = prd.Ya(prd.Y(dr.getValueY(i), Space.Ly, 0, Space.Ly), Cent.getValueY(my_num), Space.Ly, 0, Space.Ly);
				DL = Yi * Yi + DL;
				Zi = prd.Za(prd.Z(dr.getValueZ(i), Space.Lz, 0, Space.Lz), Cent.getValueZ(my_num), Space.Lz, 0, Space.Lz);
				DL = Zi * Zi + DL;
			}
		}
		//строка 464
		if (DL > VDmin)
		{
			for (i = 0; i < NeibsM; i++)
			{
				if (ContactM.getValue(i) < NCHKL.getValue(i))
				{
					for (k = 0; k < NCHKL.getValue(i); k++)
					{
						int ii = NeibsNum.getValue(i);
						if ((NatToAll[i].getValueCx(k)*Xi + NatToAll[i].getValueCy(k)*Yi + NatToAll[i].getValueCz(k)*Zi + NatToAll[i].getValueRR(k)) < 6*rView)
						{
							if (ContactM.getValue(i) < k)
							{
								ContactM.setValue(i, k);

							}
						}
					}
				}
			}
		}
	}
	
	Np1 = my_num - 1;
	Np2 = my_num + 1;

	if (Np1 == 0) Np1 = num_proc - 1;
	if (Np2 > num_proc - 1) Np2 = 1;


}


void LookingNeibs2(NullPtrsInt &NeibsNum, NullPtrsInt &ContactM, NullPtrsInt &MCHKL, NullPtrsInt &NCHKL, MPI_Request ReqS[], MPI_Request ReqR[], MPI_Status &status, NullPtrsReal &Bufer)
{
	int i;
	for (i = 0; i < NeibsM; i++)
	{
		Neib = NeibsNum.getValue(i);

		if (Neib == Np1 || Neib == Np2)
		{
			if (ContactM.getValue(i) < 1) ContactM.setValue(i, 1);
		}
		int mchkl;
		//int mchkl = MCHKL.getValue(i);
		int contactm = ContactM.getValue(i);

		MPI_Irecv(&mchkl, 4, MPI_BYTE, Neib, 5000, MPI_COMM_WORLD, &ReqR[i]);// 5000 - контакты
		MPI_Isend(&contactm, 4, MPI_BYTE, Neib, 5000, MPI_COMM_WORLD, &ReqS[i]);
		MPI_Wait(&ReqR[i], &status);
		MPI_Wait(&ReqS[i], &status);
		MCHKL.setValue(i, mchkl);

	}


	Bufer.setValue(my_num, (rView/rCut)*10E-9);
	BuferM = Bufer.getValue(my_num);
	//std::cout << BuferM << std::endl;
	real bufer_recv;

	for (i = 1; i < num_proc; i++) //передача буферов другим процессам
	{		
		if (i != my_num)
		{
			MPI_Irecv(&bufer_recv, 1, MPI_DOUBLE, i, 6000, MPI_COMM_WORLD, &ReqR[i]);// 6000 - буферы
			MPI_Isend(&BuferM, 1, MPI_DOUBLE, i, 6000, MPI_COMM_WORLD, &ReqS[i]);
			MPI_Wait(&ReqR[i], &status);
			MPI_Wait(&ReqS[i], &status);
			Bufer.setValue(i, bufer_recv);
		}			
	}

	Neibs = 0;

	for (i = 0; i < NeibsM; i++)
	{
		if (ContactM.getValue(i) > 0 || MCHKL.getValue(i) > 0) 
		{
			ContactM.setValue(i, std::max(ContactM.getValue(i), MCHKL.getValue(i)));

			int ii = NeibsNum.getValue(i);
			if (BuferM < Bufer.getValue(ii)) BuferM = Bufer.getValue(ii);
			Neibs+=1;

			if (Neibs < i)
			{
				Neib = NeibsNum.getValue(Neib);//[Neib];
				NeibsNum.setValue(Neib, NeibsNum.getValue(ii));
				NeibsNum.setValue(i, Neib);
				//вызов сет мб строка 535
			}
			else
			{
				ContactM.setValue(i, NCHKL.getValue(i));
				NeibsNum.setValue(i, -1*NeibsNum.getValue(i));
			}

		}
	}
	BuferM = 0.25*BuferM + 0.75*Bufer.getValue(my_num);
	//std::cout << BuferM << std::endl;
	
}


void CheckList(int NeibsM, NullPtrsInt &MCHKL, NullPtrsInt &NeibsNum, MPI_Request ReqS[], MPI_Request ReqR[], MPI_Status &status)
{
	int i, j;
	int neibsnum, mchkl;
	for (i = 1; i < num_proc; i++)
	{
		if (i != my_num)
		{
			MCHKL.setValue(i, 0);
			if (NeibsM < (num_proc - 1))
			{
				for (j = 0; j < NeibsM; j++)
				{
					if (NeibsNum.getValue(j) == i)
					{
						MCHKL.setValue(i, i);
						break;
					}
				}
			}
			else MCHKL.setValue(my_num, i);
			//std::cout << MCHKL.getValue(i) << "    " << NeibsNum.getValue(i);

			neibsnum = NeibsNum.getValue(i);
			mchkl = MCHKL.getValue(i);

			MPI_Irecv(&neibsnum, 1, MPI_INT, i, 100, MPI_COMM_WORLD, &ReqR[i]);
			MPI_Isend(&mchkl, 1, MPI_INT, i, 100, MPI_COMM_WORLD, &ReqS[i]);
		}
	}

	for (i = 1; i < num_proc; i++)
	{
		if (i != my_num)
		{
			MPI_Wait(&ReqR[i], &status);
			MPI_Wait(&ReqS[i], &status);
		}
	}
		//return &MCHKL;
}


void AtomExchange(std::vector<NatTo>&NatToAll, std::vector<NatTo>&NatToMe, CalcSpace &Space, Period &prd, 
				  NullPtrsInt &NeibsNum, NullPtrsInt &KindA, NullPtrsInt &ContactM, NullPtrsVecR &Cent, NullPtrsVecR &dr, NullPtrsVecR &Vel,
				  MPI_Status &status, MPI_Request reqS[], MPI_Request reqR[])
{
	int i, m, mm, N2, N1, n, Np1, Np2;
	double qq;

	double DL = 0;
	int NeibsX = 0;
	int tmp;

	std::vector<int> NeiN(NGen), NeiA(NGen);


	for (i = 0; i < NeibsM; i++)
	{
		NatToAll[i].setValueN1(0); NatToMe[i].setValueN1(0);
		NatToAll[i].setValueN2(0); NatToMe[i].setValueN2(0);
	}

	NAsnd = -1;
	real bb = Space.Lx * Space.Lx + Space.Ly * Space.Ly + Space.Lz * Space.Lz;

	for (i = 0; i < NeibsM; i++)
	{
		int ii = NeibsNum.getValue(i);
		//std::cout << NeibsNum.getValue(i) << std::endl;
		//int ii = i;
		if (ii > 0)
		{
			if (ii != my_num)
			{
				DL = prd.Xa(Cent.getValueX(ii), Cent.getValueX(my_num), Space.Lx, 0, Space.Lx) *
					prd.Xa(Cent.getValueX(ii), Cent.getValueX(my_num), Space.Lx, 0, Space.Lx);
				DL = prd.Ya(Cent.getValueY(ii), Cent.getValueY(my_num), Space.Ly, 0, Space.Ly) *
					prd.Ya(Cent.getValueY(ii), Cent.getValueY(my_num), Space.Ly, 0, Space.Ly) + DL;
				DL = prd.Za(Cent.getValueZ(ii), Cent.getValueZ(my_num), Space.Lz, 0, Space.Lz) *
					prd.Za(Cent.getValueZ(ii), Cent.getValueZ(my_num), Space.Lz, 0, Space.Lz) + DL;
			}

			if (bb > DL) bb = DL;
		}
		
		
		
	}
	
	bb = 0.25 * bb;
	//std::cout << bb << "  " << my_num << std::endl;
	//NpartW = MyN;
	for (i = 0; i < NpartW; i++)
	{
		if (KindA.getValue(i) > 0)
		DL = prd.Xa(dr.getValueX(i), Cent.getValueX(my_num), Space.Lx, 0, Space.Lx) * prd.Xa(dr.getValueX(i), Cent.getValueX(my_num), Space.Lx, 0, Space.Lx) + 
			prd.Ya(dr.getValueY(i), Cent.getValueY(my_num), Space.Ly, 0, Space.Ly) * prd.Ya(dr.getValueY(i), Cent.getValueY(my_num), Space.Ly, 0, Space.Ly) + 
			prd.Za(dr.getValueZ(i), Cent.getValueZ(my_num), Space.Lz, 0, Space.Lz) * prd.Za(dr.getValueZ(i), Cent.getValueZ(my_num), Space.Lz, 0, Space.Lz);
		int ii;
		//std::cout << DL << std::endl;
		if (DL > bb)
		{
			mm = 0;
			for (m = 0; m < NeibsM; m++)
			{
				ii = NeibsNum.getValue(m);
				//std::cout << my_num << "  " << ii << std::endl;
				if (ii > 0)
				{
					qq = prd.Xa(dr.getValueX(i), Cent.getValueX(ii), Space.Lx, 0, Space.Lx) * prd.Xa(dr.getValueX(i), Cent.getValueX(ii), Space.Lx, 0, Space.Lx) + 
						prd.Ya(dr.getValueY(i), Cent.getValueY(ii), Space.Ly, 0, Space.Ly) * prd.Ya(dr.getValueY(i), Cent.getValueY(ii), Space.Ly, 0, Space.Ly) +
						prd.Za(dr.getValueZ(i), Cent.getValueZ(ii), Space.Lz, 0, Space.Lz) * prd.Za(dr.getValueZ(i), Cent.getValueZ(ii), Space.Lz, 0, Space.Lz);

					if (DL > qq)
					{
						DL = qq;
						mm = m;
					}
				}
				
			}
			if (mm > 0)
			{
				NAsnd += 1;
				NeiN[NAsnd] = mm;
				//std::cout << NeiN[NAsnd] << std::endl;
				NeiA[NAsnd] = i;
				NatToAll[mm].setValueN2(NatToAll[mm].getValueN2() + 1);
				

			}
		}
	}


	int nat2al1S, nat2al2S;
	int nat2men1R, nat2men2R;

	NArcv = 0;
	for (i = 0; i < NeibsM; i++)
	{
		int ii = NeibsNum.getValue(i);
		
		if (ii > 0)
		{
			nat2al2S = NatToAll[i].getValueN2();
			MPI_Irecv(&nat2men2R, 1, MPI_INT, ii, 102, MPI_COMM_WORLD, &reqR[i]);
			MPI_Isend(&nat2al2S, 1, MPI_INT, ii, 102, MPI_COMM_WORLD, &reqS[i]);
			MPI_Wait(&reqR[i], &status);
			MPI_Wait(&reqS[i], &status);
			
			NatToMe[i].setValueN2(nat2men2R);
			NArcv += nat2men2R;
		}
		//std::cout << my_num << "  " << nat2al2S << std::endl;
		//std::cout << my_num << "  " << ii << std::endl;
		//}
		
	}
	//std::cout << my_num << "  " << NArcv << std::endl;

	struct MPexc
	{
		real X, Y, Z, Vx, Vy, Vz;
		int Kind;
		MPexc(): X(0.), Y(0.), Z(0.), Vx(0.), Vy(0.), Vz(0.), Kind(0) {}
	};
	

	//std::vector<MPexc> MPrecvA(NArcv), MPsendA(NAsnd);
	std::vector<MPexc> MPrecvA, MPsendA;
	for (i = 0; i < 1.25*NArcv; i++) MPrecvA.push_back(MPexc());
	for (i = 0; i < NAsnd; i++) MPsendA.push_back(MPexc());

	//MPrecvA.resize(NArcv, 0);
	//MPrecvA.resize(NArcv, 0); MPsendA.resize(NAsnd, 0);


	int d = 6*8+4;
	//int d = 8;
	//MPI_Request *tmpR = new MPI_Request [NArcv];
	//MPI_Request *tmpS = new MPI_Request [NAsnd];

	if (NArcv > 0) //строка 991 - 1010
	{
		N2 = 0;
		for (i = 0; i < NeibsM; i++)
		{
			//std::cout << NatToMe[i].getValueN2() << "   " << my_num << "   " << i << std::endl;
			int ii = NeibsNum.getValue(i);
			if (ii > 0)
			{
				//if (ii != my_num)
				//{
					if (NatToMe[i].getValueN2() > 0)//убрать комменты
					{
						
						MPI_Irecv(&MPrecvA[N2], d*NatToMe[i].getValueN2(), MPI_BYTE, ii, 401, MPI_COMM_WORLD, &reqR[i]);
						N2 += NatToMe[i].getValueN2();

						//std::cout << N2 << "   " << my_num << "   " << i << std::endl;
					}
				//}
			}			
		}
	}

	N1 = 0;
	//std::cout << NArcv << "  111  " << my_num << std::endl;
	//std::cout << N2 << "   " << my_num << std::endl;
	if (NAsnd > 0)
	{
		for (m = 0; m < NeibsM; m++)
		{
			int ii = NeibsNum.getValue(m); 
			if (ii > 0)
			{
				if (ii != my_num)
				{
					if (NatToAll[m].getValueN2() > 0)
					{
						N2 = 0; int N2_1 = 0;
						for (n = 0; n < NAsnd; n++)
						{
							if (NeiN[n] == m) //обратить внимание сюда
							{
								i = NeiA[n];
								//std::cout << ii << std::endl;
								
								//MPsendA[1].X = dr.getValueX(i);
								MPsendA[N1].X = dr.getValueX(i);
								//std::cout << MPsendA[N1].X << std::endl;
								MPsendA[N1].Y = dr.getValueY(i);
								MPsendA[N1].Z = dr.getValueZ(i);
								MPsendA[N1].Vx = Vel.getValueX(i);
								MPsendA[N1].Vy = Vel.getValueY(i);
								MPsendA[N1].Vz = Vel.getValueZ(i);
								MPsendA[N1].Kind = KindA.getValue(i);
								KindA.setValue(i, -5);
								N1+=1;
								N2+=1;
								//N2_1+=1;
								NpartA = NpartA - 1;
							}
					
						}
						
						N2 = N1 - N2;
						
						MPI_Isend(&MPsendA[N2], d*NatToAll[m].getValueN2(), MPI_BYTE, ii, 401, MPI_COMM_WORLD, &reqS[m]);

					}
				}
			}
		}
	}

	int NpartA1 = 0;  int NpartA2 = 0;  int NpartA3 = 0; int NpartW = 0;
	//int Npart = MyN; 
	int NpartP;
	int Npend = NGen;


	for (i = 0; i < MyN; i++)
	{
		if (KindA.getValue(i) > 0)
		{
			NpartW = i;
			if (KindA.getValue(i) == m1) NpartA1+=1;
			if (KindA.getValue(i) == m2) NpartA2+=1;
			if (KindA.getValue(i) == m3) NpartA3+=1;
		}
	}
	int NpartA = NpartA1 + NpartA2 + NpartA3;
	int Nalien = 1.25*N1;


	if (NArcv > 0)
	{
		if ((NArcv + NpartA) > MyN)
		{
			Np1 = 1.25*MyN + 3*NArcv;
			if (Nalien < 2*(Npend - MyN)) Np2 = (5*Nalien)/4;
			else Np2 = Nalien;
			
			NeiA.resize(Np1 + Np2, 0); NeiN.resize(Np1 + Np2, 0); 
			//тут функция allcar
		}
	
		NpartP = std::min(NpartW + 1, MyN);

		for (i = 0; i < NpartW; i++)
		{
			if (KindA.getValue(i) <= 0)
			{
				NpartP = i;
				break;
			}
		}

		N2 = 0;

		for (m = 0; m < NeibsM; m++)
		{
			int ii = NeibsNum.getValue(m);
			if (ii > 0)
			{
				if (ii != my_num)
				{
					if (NatToAll[m].getValueN2() > 0)
					{
						MPI_Wait(&reqR[m], &status);
						MPI_Wait(&reqS[m], &status);

						real Rx = prd.Xa(Cent.getValueX(my_num), Cent.getValueX(ii), Space.Lx, 0, Space.Lx);//по коду ВВ локальные центры масс
						real Ry = prd.Ya(Cent.getValueY(my_num), Cent.getValueY(ii), Space.Ly, 0, Space.Ly);
						real Rz = prd.Za(Cent.getValueZ(my_num), Cent.getValueZ(ii), Space.Lz, 0, Space.Lz);
						N1 = 0;

						for (i = NpartP; i < MyN; i++)
						{
							if (KindA.getValue(i) <= 0)
							{
								

								dr.setValueX(i, prd.Xa(MPrecvA[N2].X, Rx, Space.Lx, 0, Space.Lx));
								//std::cout<<MPrecvA[N2].X<<std::endl;
								dr.setValueY(i, prd.Ya(MPrecvA[N2].Y, Ry, Space.Ly, 0, Space.Ly));
								dr.setValueZ(i, prd.Za(MPrecvA[N2].Z, Rz, Space.Lz, 0, Space.Lz));

								Vel.setValueX(i, MPrecvA[N2].Vx);
								Vel.setValueY(i, MPrecvA[N2].Vy);
								Vel.setValueZ(i, MPrecvA[N2].Vz);

								NeiN[i] = 0;  NeiA[i] = 0; KindA.setValue(i, MPrecvA[N2].Kind);

								N1+=1;
								N2+=1;

								if (KindA.getValue(i) == m1) NpartA1+=1;
								if (KindA.getValue(i) == m2) NpartA2+=1;
								if (KindA.getValue(i) == m3) NpartA3+=1;

								if( N1 == NatToMe[m].getValueN2())
								{
									NpartA = NpartA + N1;
									NpartP = i + 1;
									if ( NpartW < i ) NpartW = i;
								}
							}
						}
					}
				}
			}
		}
	}

	//кусок с переаллокациями

}



void MPInf_0()
{

}


void ThermLan(/*std::vector<double>&GaussDist,*/ NullPtrsVecR &LanForce, CalcConstsLD Consts, NullPtrsVecR &Vel)
{
	std::random_device rdx, rdy, rdz;
	std::mt19937 genx(rdx()), geny(rdy()), genz(rdz());
	std::normal_distribution<> distrx(0, 1), distry(0, 1), distrz(0, 1);
	
	double Tmpr2 = 60 * 1.38E-23;
	double Beta2 = 1E+12;//1E+12;
	double Sigm2 = Beta2 * Tmpr2 * 2 / dt;
	double Rx = 1; //тут страшная формула в будущем будет
	double Vx = 0; //тут страшная формула в будущем будет
	std::vector<double> GaussDistX(NGen, 0);
	std::vector<double> GaussDistY(NGen, 0);
	std::vector<double> GaussDistZ(NGen, 0);

	//double x = 0, x2 = 0;

	for (int i = 0; i < NGen; i++) //у каждой оси своё распредление
	{
		GaussDistX[i] = distrx(genx);
		//x += GaussDistX[i];
		//x2 += GaussDistX[i] * GaussDistX[i];
		GaussDistY[i] = distry(geny);
		GaussDistZ[i] = distrz(genz); //определять для тех, кто в списке Верле
	}
	//int i = 0;
	//std::cout << x/NGen << std::endl;
	//std::cout << x2/NGen << std::endl;
	/*while (i < NGen)
	{
		if (distr(gen) > -1 && distr(gen) < 1)
		{
			GaussDist[i] = distr(gen); i += 1;
		}
	}*/


	for (int i = 0; i < NGen; i++)
	{
		double Beta = m1;//Consts.m; //в данном случае масса атома
        double Sigm = sqrt( Sigm2 * Beta );// * Rx; //rx единица //sigm2 - как сигма0 считать
        Beta = Beta2 * Beta;// * Rx;

		//fx += (Sigm * GaussDist[i] - Beta * (Vel.getValueX(i) - Vx));
		//fy += (Sigm * GaussDist[i] - Beta * (Vel.getValueY(i) - Vx));
		//fz += (Sigm * GaussDist[i] - Beta * (Vel.getValueZ(i) - Vx));
		LanForce.setValueX(i, (Sigm * GaussDistX[i] - Beta * (Vel.getValueX(i))) );
		//std::cout << LanForce.getValueX(i) << std::endl;
		LanForce.setValueY(i, (Sigm * GaussDistY[i] - Beta * (Vel.getValueY(i))) );
		LanForce.setValueZ(i, (Sigm * GaussDistZ[i] - Beta * (Vel.getValueZ(i))) );
	}


}


//***********************************СКАЛЯРНОЕ******************************************

void CalculateNGen_skal(CalcSpace Space, CalcConstsLD Consts)
{
	for (Space.i1 = -Space.n1; Space.i1 < Space.n1; Space.i1++)
	{
		for (Space.i2 = -Space.n2; Space.i2 < Space.n2; Space.i2++)
		{
			for (Space.i3 = -Space.n3; Space.i3 < Space.n3; Space.i3++)
			{				
				Space.dx = 0.5 * Consts.a * Space.i1 + 0.5 * Consts.a * Space.i3;
				Space.dy = 0.5 * Consts.a * Space.i1 + 0.5 * Consts.a * Space.i2;
				Space.dz = 0.5 * Consts.a * Space.i2 + 0.5 * Consts.a * Space.i3;
				
				if (Space.dx >= 0 && Space.dx <= Space.Lx &&
					Space.dy >= 0 && Space.dy <= Space.Ly &&
					Space.dz >= 0 && Space.dz <= Space.Lz) NGen +=1;

				/*if (Space.dx >= Space.Lx * 0.25 && Space.dx <= Space.Lx * 0.75 &&
					Space.dy >= Space.Lx * 0.25 && Space.dy <= Space.Ly * 0.75 &&
					Space.dz >= Space.Lx * 0.25 && Space.dz <= Space.Lz * 0.75) NGen +=1;*/

			}
		}
	}
	//NGen = NGen*1.75;
}

void InitValuesCoord_skal(CalcSpace Space, CalcConstsLD Consts, NullPtrsVecR &dr)
{
	int i = 0;
	for (Space.i1 = -Space.n1; Space.i1 < Space.n1; Space.i1++)
	{
		for (Space.i2 = -Space.n2; Space.i2 < Space.n2; Space.i2++)
		{
			for (Space.i3 = -Space.n3; Space.i3 < Space.n3; Space.i3++)
			{				
				Space.dx = 0.5 * Consts.a * Space.i1 + 0.5 * Consts.a * Space.i3;
				Space.dy = 0.5 * Consts.a * Space.i1 + 0.5 * Consts.a * Space.i2;
				Space.dz = 0.5 * Consts.a * Space.i2 + 0.5 * Consts.a * Space.i3;
				
				/*if (Space.dx >= Space.Lx * 0.25 && Space.dx <= Space.Lx * 0.75 &&
					Space.dy >= Space.Lx * 0.25 && Space.dy <= Space.Ly * 0.75 &&
					Space.dz >= Space.Lx * 0.25 && Space.dz <= Space.Lz * 0.75)*/
				if (Space.dx >= 0 && Space.dx <= Space.Lx &&
					Space.dy >= 0 && Space.dy <= Space.Ly &&
					Space.dz >= 0 && Space.dz <= Space.Lz)
				{
					dr.setValueX(i, Space.dx);
					dr.setValueY(i, Space.dy);
					dr.setValueZ(i, Space.dz);
					i +=1;
				}
			}
		}
	}
}

void InitValuesOthers_skal(int NGen, NullPtrsVecR &Vel, NullPtrsVecR &Accel, NullPtrsInt &EnergyPot, NullPtrsVecR &LanForce)
{
	int i;
	std::random_device rdx, rdy, rdz;
	//std::mt19937 genx(rdx()), geny(rdy()), genz(rdz()); //поставить начальное число
	std::mt19937 genx(5), geny(5), genz(5);
	std::normal_distribution<> distrx(0, 1), distry(0, 1), distrz(0, 1);

	/*std::vector<double> GaussDistX(NGen, 0);
	std::vector<double> GaussDistY(NGen, 0);
	std::vector<double> GaussDistZ(NGen, 0);*/
	for (i = 0; i < NGen; i++)
	{
		//double GaussDistX = distrx(genx) * 125;
		//double GaussDistY = distry(geny) * 125;
		//double GaussDistZ = distrz(genz) * 125;
		//
		//Vel.setValueX(i, GaussDistX); //задать скорость через (Гауссово распределение * 100) //термостат выключен сразу
		//Vel.setValueY(i, GaussDistY); 
		//Vel.setValueZ(i, GaussDistZ);

		Vel.setValueX(i, 0); //задать скорость через (Гауссово распределение * 100) //термостат выключен сразу
		Vel.setValueY(i, 0); //посмотреть координаты, выгрузить
		Vel.setValueZ(i, 0);
		Accel.setValueX(i, 0);
		Accel.setValueY(i, 0);
		Accel.setValueZ(i, 0);
		EnergyPot.setValue(i, 0);

		LanForce.setValueX(i, 0);
		LanForce.setValueY(i, 0);
		LanForce.setValueZ(i, 0);
	}
}

