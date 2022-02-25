#ifndef __my_func_h__
#define __my_func_h__

int Comp1(const void * a, const void * b);


int Round (real value);

//real prd_Xa(real Ri, real Rj, real Right, real Left, real L);
//real prd_Ya(real Ri, real Rj, real Right, real Left, real L);
//real prd_Za(real Ri, real Rj, real Right, real Left, real L);


//Nat2Al DirCos(int k, int n, real Rx, real Ry, real Rz, Nat2Al Nat_2_Al[]);
//Nat2Al MpPair(Nat2Al Nat_2_Al[], int my_num, int n, real Lx, real Ly, real Lz, VecR Cent[]);

int SearchNeibAmount(int my_num, int num_proc, 
				real r_view, real Lx, real Ly, real Lz, 
				real r_max[], VecR Cent[]);
int SearchNeibNumb(int my_num, int num_proc, real r_view, real Lx, real Ly, real Lz, 
					 int Neibs_num[], real r_max[], VecR Cent[]);


void CalculateNGen(CalcSpace Space, CalcConstsLD Consts, real Cent);


void InitValuesCoord(CalcSpace Space, CalcConstsLD Consts, real Cent, NullPtrsVecR &dr);


void InitValuesOthers(int NGen, int MyN, NullPtrsVecR &Vel, NullPtrsVecR &Accel, NullPtrsInt &EnergyPot, NullPtrsInt &KindA, NullPtrsVecR &LanForce);


real RMax(int NGen, NullPtrsVecR &dr);


void LookingNeibs1(NullPtrsInt &NeibsNum, NullPtrsInt &ContactM, NullPtrsInt &NCHKL, NullPtrsInt &KindA, NullPtrsReal &rMax, NullPtrsVecR &dr, 
				   CalcSpace Space, NullPtrsVecR &Cent, Period &prd, std::vector<NatTo>&NatToAll/*NatTo NatToAll[]*/);
//void SendRecvRMax(NullPtrsReal rMax, MPI_Request ReqR[], MPI_Request ReqS[]);


void LookingNeibs2(NullPtrsInt &NeibsNum, NullPtrsInt &ContactM, NullPtrsInt &MCHKL, NullPtrsInt &NCHKL, MPI_Request ReqS[], MPI_Request ReqR[], MPI_Status &status, NullPtrsReal &Bufer);


void CheckList(int NeibsM, NullPtrsInt &MCHKL, NullPtrsInt &NeibsNum, MPI_Request ReqS[], MPI_Request ReqR[], MPI_Status &status);


void AtomExchange(std::vector<NatTo>&NatToAll, std::vector<NatTo>&NatToMe, CalcSpace &Space, Period &prd, 
				  NullPtrsInt &NeibsNum, NullPtrsInt &KindA, NullPtrsInt &ContactM, NullPtrsVecR &Cent, NullPtrsVecR &dr, NullPtrsVecR &Vel,
				  MPI_Status &status, MPI_Request reqS[], MPI_Request reqR[]);



void ThermLan(/*std::vector<double>&GaussDist,*/ NullPtrsVecR &LanForce, CalcConstsLD Consts, NullPtrsVecR &Vel);


//***********************************— ¿Àﬂ–ÕŒ≈******************************************

void CalculateNGen_skal(CalcSpace Space, CalcConstsLD Consts);


void InitValuesCoord_skal(CalcSpace Space, CalcConstsLD Consts, NullPtrsVecR &dr);


void InitValuesOthers_skal(int NGen, NullPtrsVecR &Vel, NullPtrsVecR &Accel, NullPtrsInt &EnergyPot, NullPtrsVecR &LanForce);















#endif