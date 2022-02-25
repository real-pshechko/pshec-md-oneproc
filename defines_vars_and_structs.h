#ifndef __defines_vars_and_structs__
#define __defines_vars_and_structs__

typedef double real;

//    СТРУКТУРЫ
struct VecR
{
	real x, y, z;
};


struct Mol
{
	VecR r, rv, ra;
};


struct CalcSpace
{
	int n1, n2, n3,
	i1, i2, i3;
	real Lx, Ly, Lz;//длина расчетной области
	real dx, dy, dz,
	delta_x, delta_y, delta_z;
};


//    КЛАССЫ
class CalcConstsLD
{
public:
	real a, eps, sig, dt, m;
	CalcConstsLD():a(6.2E-10), eps(1.65E-21), sig(3.465E-10), dt(4E-15), m(6.63E-26)
	{

	}

	void de_dimension()
	{
		//ОБЕЗРАЗМЕРИВАНИЕ СДЕЛАТЬ
	}

};




// на векторах
class NullPtrsInt
{
private:
	//int *ptr;
	std::vector<int> ptr;
	
public:
	
	NullPtrsInt(int leight) //конструктор
	{
		ptr.resize(leight, 0);//(leight);//добавлять через append или resize
	}

	virtual ~NullPtrsInt() //деструктор
	{
		ptr.clear();//удалять через pop или .... //погуглить очистку вектора
	}
	
	void setValue(int index, int value) { ptr.at(index) = value; }//ptr[index] = value; } //установить значение

	int getValue(int index) { return ptr.at(index); }//return ptr[index]; } //вернуть значение
};


class NullPtrsReal
{
private:
	//real *ptr;
	std::vector<real> ptr;
	
public:

	NullPtrsReal(int leight) //конструктор
	{
		ptr.resize(leight, 0);
	}

	virtual ~NullPtrsReal() //деструктор
	{
		ptr.clear();
	}
	
	void setValue(int index, real value) { ptr.at(index) = value; } //установить значение

	real getValue(int index) { return ptr.at(index); } //вернуть значение
};


class NullPtrsVecR
{
private:
	//VecR *ptr;
	std::vector<VecR> ptr;
	
public:

	NullPtrsVecR(int leight) //конструктор
	{
		ptr.resize(leight);
	}

	virtual ~NullPtrsVecR() //деструктор
	{
		ptr.clear();
	}
	
	void setValueX(int index, real value) { ptr.at(index).x = value; } //установить значение x
	void setValueY(int index, real value) { ptr.at(index).y = value; } //установить значение y
	void setValueZ(int index, real value) { ptr.at(index).z = value; } //установить значение z

	real getValueX(int index) { return ptr.at(index).x; } //вернуть значение x
	real getValueY(int index) { return ptr.at(index).y; } //вернуть значение y
	real getValueZ(int index) { return ptr.at(index).z; } //вернуть значение z
};


class Period
{
private:
	real prd_X, prd_Y, prd_Z;

public:
	Period()
	{
	}

	virtual ~Period()
	{
	}

	real Xa(real Ri, real Rj, real Right, real Left, real L)
	{
		real prd_X = Ri - Rj;

		if (prd_X > Right)
		{
			prd_X = prd_X - L;
		}
		else if (prd_X < Left)
		{
			prd_X = prd_X + L;
		}

		return prd_X;
	}
	real Ya(real Ri, real Rj, real Right, real Left, real L)
	{
		real prd_Y = Ri - Rj;

		if (prd_Y > Right)
		{
			prd_Y = prd_Y - L;
		}
		else if (prd_Y < Left)
		{
			prd_Y = prd_Y + L;
		}

		return prd_Y;
	}
	real Za(real Ri, real Rj, real Right, real Left, real L)
	{
		real prd_Z = Ri - Rj;

		if (prd_Z > Right)
		{
			prd_Z = prd_Z - L;
		}
		else if (prd_Z < Left)
		{
			prd_Z = prd_Z + L;
		}

		return prd_Z;
	}

	real X(real x, real right, real left, real L)
	{
		real prd_x = x;
		if (prd_x > right) prd_x = prd_x - L;
		if (prd_x < left) prd_x = prd_x + L;
		return prd_x;
	}
	real Y(real y, real right, real left, real L)
	{
		real prd_y = y;
		if (prd_y > right) prd_y = prd_y - L;
		if (prd_y < left) prd_y = prd_y + L;
		return prd_y;
	}
	real Z(real z, real right, real left, real L)
	{
		real prd_z = z;
		if (prd_z > right) prd_z = prd_z - L;
		if (prd_z < left) prd_z = prd_z + L;
		return prd_z;
	}


};


class OperatingWithNeibs
{
private:
	int m_my_num, m_num_proc, m_NeibsM;
	real m_r_view, m_Lx, m_Ly, m_Lz;

public:
	OperatingWithNeibs(int my_num, int num_proc, real rView, real Lx, real Ly, real Lz)//(): // m_my_num(my_num), m_num_proc(num_proc), m_r_view(rView), m_Lx(Lx), m_Ly(Ly), m_Lz(Lz)//переделать как вверху
	{
		m_my_num = my_num;
		m_num_proc = num_proc;
		m_r_view = rView;
		m_Lx = Lx;
		m_Ly = Ly;
		m_Lz = Lz;
	}

	virtual ~OperatingWithNeibs()
	{
	}

	int SearchNeibAmount(NullPtrsReal &rMax, NullPtrsVecR &Cent, Period &prd) //this использовать
	{
		int Np1, Np2, i;
		real RM, Rx, Ry, Rz;
		real bb;

		m_NeibsM = 0;

		Np1 = m_my_num - 1;
		Np2 = m_my_num + 1;
	
		if (Np1 == 0) Np1 = m_num_proc - 1;
		if (Np2 > (m_num_proc-1)) Np2 = 1;

		for (i=1; i < m_num_proc; i++)
		{
			if (i != m_my_num)
			{
				RM = rMax.getValue(m_my_num) + rMax.getValue(i) + m_r_view*6;
				RM = RM*RM;

				Rx = prd.Xa(Cent.getValueX(m_my_num), Cent.getValueX(i), m_Lx, 0, m_Lx);
				Ry = prd.Ya(Cent.getValueY(m_my_num), Cent.getValueY(i), m_Ly, 0, m_Ly);
				Rz = prd.Za(Cent.getValueZ(m_my_num), Cent.getValueZ(i), m_Lz, 0, m_Lz);
				
				bb = Rx * Rx + Ry * Ry + Rz * Rz;

				if (bb < RM || i == Np1 || i == Np2) m_NeibsM+=1;
			}
		}
		return m_NeibsM;
	}

	NullPtrsInt *SearchNeibNumb(NullPtrsReal &rMax, NullPtrsVecR &Cent, Period &prd, NullPtrsInt &Neibs_num)
	{
		int Np1, Np2, i, j;
		real RM, Rx, Ry, Rz;
		real bb;
	
		Np1 = m_my_num - 1;
		Np2 = m_my_num + 1;
	
		if (Np1 == 0) Np1 = m_num_proc - 1;
		if (Np2 > (m_num_proc-1)) Np2 = 1;
	
		j = 0;
	
		for (i = 1; i < m_num_proc; i++)
		{
			if (i != m_my_num)
			{
				RM = rMax.getValue(m_my_num) + rMax.getValue(i) + m_r_view*6;
				RM = RM*RM;
					
				Rx = prd.Xa(Cent.getValueX(m_my_num), Cent.getValueX(i), m_Lx, 0, m_Lx);
				Ry = prd.Ya(Cent.getValueY(m_my_num), Cent.getValueY(i), m_Ly, 0, m_Ly);
				Rz = prd.Za(Cent.getValueZ(m_my_num), Cent.getValueZ(i), m_Lz, 0, m_Lz);
	
				bb = Rx*Rx + Ry*Ry + Rz*Rz;
	
				if (bb < RM || i == Np1 || i == Np2) {Neibs_num.setValue(j, i); j+=1;}
			}
		}
		return &Neibs_num;
	}


};


class NatTo
{
private:
	real RR[7], Cx[7], Cy[7], Cz[7];
	int N1, N2;

public:

	NatTo()
	{
	}

	virtual ~NatTo()
	{
	}

	void setValueRR(int index, real value) { RR[index] = value; }
	void setValueCx(int index, real value) { Cx[index] = value; }
	void setValueCy(int index, real value) { Cy[index] = value; }
	void setValueCz(int index, real value) { Cz[index] = value; }

	void setValueN1(int value) { N1 = value; }
	void setValueN2(int value) { N2 = value; }


	real getValueRR(int index) { return RR[index]; }
	real getValueCx(int index) { return Cx[index]; }
	real getValueCy(int index) { return Cy[index]; }
	real getValueCz(int index) { return Cz[index]; }

	int getValueN1() { return N1; }
	int getValueN2() { return N2; }

	void DirCos(int k, int n, real Rx, real Ry, real Rz, std::vector<NatTo>&Nat)//(int k, int n, real Rx, real Ry, real Rz, NatTo Nat[]) //сделать void
	{
		real qq;
		qq = sqrt(Rx*Rx + Ry*Ry + Rz*Rz);
		Nat[n].setValueRR(k, qq);

		if (qq > 0)
		{
			Nat[n].setValueCx(k, -Rx/qq);
			Nat[n].setValueCy(k, -Ry/qq);
			Nat[n].setValueCz(k, -Rz/qq);
		}
		else
		{
			Nat[n].setValueCx(k, -1);
			Nat[n].setValueCy(k,  0);
			Nat[n].setValueCz(k,  0);
		}
	}

	NatTo MpPair(std::vector<NatTo>&Nat, int my_num, int n, NullPtrsVecR &Cent, CalcSpace Space)//(NatTo Nat[], int my_num, int n, NullPtrsVecR &Cent, CalcSpace Space)
	{
		int k; 
		real Wx, Wy, Wz; 
		real Rx, Ry, Rz, qq;
	
		for (k=0; k<8; k++)
		{
			Nat[n].setValueRR(k, (Space.Lx*Space.Lx+Space.Ly*Space.Ly+Space.Lz*Space.Lz));
			Nat[n].setValueCx(k, -1);
			Nat[n].setValueCy(k,  0);
			Nat[n].setValueCz(k,  0);
		}

		Wx=-(Cent.getValueX(n) - Cent.getValueX(my_num))/2;
		Wy=-(Cent.getValueY(n) - Cent.getValueY(my_num))/2;
		Wz=-(Cent.getValueZ(n) - Cent.getValueZ(my_num))/2;

		if (Wx > 0)  Rx = Wx;				else Rx = Wx + Space.Lx;	//1 quadrant
		if (Wy >= 0) Ry = Wy;				else Ry = Wy + Space.Ly;
		if (Wz >= 0) Rz = Wz;				else Rz = Wz + Space.Lz;
		//Nat_2_Al[n] = DirCos(0, n, Rx, Ry, Rz, Nat_2_Al);
		DirCos(0, n, Rx, Ry, Rz, Nat);

		if (Wx > 0)  Rx = Wx - Space.Lx;	else Rx = Wx;		//2 quadrant
		if (Wy >= 0) Ry = Wy;				else Ry = Wy + Space.Ly;
		if (Wz >= 0) Rz = Wz;				else Rz = Wz + Space.Lz;
		//Nat_2_Al[n] = DirCos(1, n, Rx, Ry, Rz, Nat_2_Al);
		DirCos(1, n, Rx, Ry, Rz, Nat);

		if (Wx > 0)  Rx = Wx;				else Rx = Wx + Space.Lx;	//3 quadrant
		if (Wy >= 0) Ry = Wy - Space.Ly;	else Ry = Wy;
		if (Wz >= 0) Rz = Wz;				else Rz = Wz + Space.Lz;
		//Nat_2_Al[n] = DirCos(2, n, Rx, Ry, Rz, Nat_2_Al);
		DirCos(2, n, Rx, Ry, Rz, Nat);

		if (Wx > 0)  Rx = Wx - Space.Lx;	else Rx = Wx;		//4 quadrant
		if (Wy >= 0) Ry = Wy - Space.Ly;	else Ry = Wy;
		if (Wz >= 0) Rz = Wz;				else Rz = Wz + Space.Lz;
		//Nat_2_Al[n] = DirCos(3, n, Rx, Ry, Rz, Nat_2_Al);
		DirCos(3, n, Rx, Ry, Rz, Nat);

		if (Wx > 0)  Rx = Wx;				else Rx = Wx + Space.Lx;	//5 quadrant
		if (Wy >= 0) Ry = Wy;				else Ry = Wy + Space.Ly;
		if (Wz >= 0) Rz = Wz - Space.Lz;	else Rz = Wz;
		//Nat_2_Al[n] = DirCos(4, n, Rx, Ry, Rz, Nat_2_Al);
		DirCos(4, n, Rx, Ry, Rz, Nat);

		if (Wx > 0)  Rx = Wx - Space.Lx;	else Rx = Wx;	//6 quadrant
		if (Wy >= 0) Ry = Wy;				else Ry = Wy + Space.Ly;
		if (Wz >= 0) Rz = Wz - Space.Lz;	else Rz = Wz;
		//Nat_2_Al[n] = DirCos(5, n, Rx, Ry, Rz, Nat_2_Al);
		DirCos(5, n, Rx, Ry, Rz, Nat);

		if (Wx > 0)  Rx = Wx;				else Rx = Wx + Space.Lx;	//7 quadrant
		if (Wy >= 0) Ry = Wy - Space.Ly;	else Ry = Wy;
		if (Wz >= 0) Rz = Wz - Space.Lz;	else Rz = Wz;
		//Nat_2_Al[n] = DirCos(6, n, Rx, Ry, Rz, Nat_2_Al);
		DirCos(6, n, Rx, Ry, Rz, Nat);

		if (Wx > 0)  Rx = Wx - Space.Lx;	else Rx = Wx;	//8 quadrant
		if (Wy >= 0) Ry = Wy - Space.Ly;	else Ry = Wy;
		if (Wz >= 0) Rz = Wz - Space.Lz;	else Rz = Wz;
		//Nat_2_Al[n] = DirCos(7, n, Rx, Ry, Rz, Nat_2_Al);
		DirCos(7, n, Rx, Ry, Rz, Nat);

		for (k = 0; k < 7; k++)
		{
			if (Nat[n].getValueRR(k) > Nat[n].getValueRR(k+1))
			{
				qq = Nat[n].getValueRR(k);
				Rx = Nat[n].getValueCx(k);
				Ry = Nat[n].getValueCy(k);
				Rz = Nat[n].getValueCz(k);

				Nat[n].setValueRR(k, Nat[n].getValueRR(k+1));
				Nat[n].setValueCx(k, Nat[n].getValueCx(k+1));
				Nat[n].setValueCy(k, Nat[n].getValueCy(k+1)); 
				Nat[n].setValueCz(k, Nat[n].getValueCy(k+1));

				Nat[n].setValueRR(k+1, qq);
				Nat[n].setValueCx(k+1, Rx);
				Nat[n].setValueCy(k+1, Ry);
				Nat[n].setValueCy(k+1, Rz);
			}
		}
		return Nat[n];
	}
};
































//ПЕРЕМЕННЫЕ С ПЛАВАЮЩЕЙ ТОЧКОЙ
extern real dt, //длина ячейки в списке Верле
	rCut, rView, rMax; //радиус обрезки и радиус наблюдения

extern real R2cut, X2min, rVr, rVr1, rVr2, aVr, aVr2, aLJ3, bLJ2, aLJ, bLJ, aLJ3, bLJ2, XX, FV, Ur;

extern real BuferM;

extern real m1, m2, m3;



//ЦЕЛОЧИСЛЕННЫЕ И ПРОЧИЕ ПЕРЕМЕННЫЕ
extern int my_num, num_proc, k;
extern int NeibsM, Neib, Neibs;
extern int Ncall;
extern int NGen, MyN, NpartW, NAsnd, NArcv, NpartA, NAinf0;
extern int Xi, Yi, Zi;
extern int Np1, Np2;
//extern int NxBox, NyBox, NzBox;
//extern MPI_Status status;





#endif