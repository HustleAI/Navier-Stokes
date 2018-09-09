#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

vector<double> progonka(vector<double> b, vector<double> c, vector<double> d, vector<double> f){

	int n = f.size();
	vector<double> x(n), delta(n), l(n);

	delta[0] = -d[0] / c[0];
	l[0] = f[0] / c[0];
	for (int i = 1; i < n; i++){
		delta[i] = -d[i] / (c[i] + b[i] * delta[i - 1]);
		l[i] = (f[i] - b[i] * l[i - 1]) / (c[i] + b[i] * delta[i - 1]);
	}
	x[n - 1] = l[n - 1];
	for (int i = n - 2; i >= 0; i--)
		x[i] = delta[i] * x[i + 1] + l[i];

	return x;
}

int main(){
	//количество узлов вдоль пластины, перпендикулярно пластине
	int XM = 401, YM = 101, XYM, NTIME = 5000;
	XM += 1; YM += 1;
	//длина области,высота расчетной области, входная скорость,кинематическая вязкость
	double L = 0.2, Y = 0.02, Uo = 0.4, Nu = 2e-5, L1 = 0.05, L2 = 0.15;
	double dx, dy, tau, eps = 1e-3, Ro = 1.225, A = 0.25/(Ro*Uo*Uo),mu=Nu*Ro;
	if (XM >= YM) XYM = XM;
	else XYM = YM;
	vector <double>   fui(XM - 1), fuj(YM), fvi(XM), fvj(YM - 1), tempui(XM - 1), tempuj(YM), tempvi(XM), tempvj(YM - 1), b(XYM), c(XYM), d(XYM);
	vector<vector<double> >  UN(XM - 1, vector<double>(YM)), U_N(XM - 1, vector<double>(YM)), UN1(XM - 1, vector<double>(YM)),
		VN(XM, vector<double>(YM - 1)), V_N(XM, vector<double>(YM - 1)), VN1(XM, vector<double>(YM - 1)),
		PN(XM, vector<double>(YM)), P_N(XM, vector<double>(YM)), PN1(XM, vector<double>(YM));

	dx = L / (XM - 2);
	dy = Y / (YM - 2);
	tau = dx / (Uo * 2);
	ofstream fout, f1out, f2out;
	fout.open("Tecplot.plt");
	f1out.open("test.txt");
	f2out.open("Pressure.plt");
	//Нужно задать н.у., чтобы быстрее разрешать погран слой
	for (int q = 0; q < NTIME; q++){											// цикл по временным шагам

		//U
		//Прогонка вдоль
		for (int j = 1; j < YM-1; j++){
			b[0] = 0; c[0] = 1; d[0] = 0; fui[0] = Uo;
			b[XM - 2] = -1; c[XM - 2] = 1; d[XM - 2] = 0; fui[XM - 2] = 0;	 

			for (int i = 1; i < XM - 2; i++){
				b[i] = tau * (-Nu / (dx * dx) - A*tau / (Ro * dx * dx) - UN[i][j] / dx);
				c[i] = 1 + 2 * tau * Nu / (dx * dx) + 2 * A * tau * tau / (Ro * dx * dx);
				d[i] = tau * (-Nu / (dx * dx) - A*tau / (Ro * dx * dx) + UN[i][j] / dx);
				fui[i] = UN[i][j] - tau / (Ro*dx) * (PN[i + 1][j] - PN[i][j]);
			}

			tempui = progonka(b, c, d, fui);
			for (int i = 0; i < XM - 1; i++)
				U_N[i][j] = tempui[i];
		}
		// Расчет давления 
		for (int i = 1; i < XM - 1; i++)
			for (int j = 1; j < YM - 1; j++)
				P_N[i][j] = PN[i][j] - tau*A*(U_N[i][j] - U_N[i - 1][j]) / dx;
		
		for (int i = 0; i < XM - 1; i++){
			U_N[i][YM - 1] = U_N[i][YM - 2];
			if (i * dx < L1 || i * dx > L2)
				U_N[i][0] = U_N[i][1];
		}

		UN1 = U_N;

		// Прогонка поперек
		for (int i = 1; i < XM - 2; i++){
			if (i * dx >= L1 && i * dx <= L2){
				b[0] = 0; c[0] = 1; d[0] = 1;  fuj[0] = 0;
			}
			else {
				b[0] = 0; c[0] = -1; d[0] = 1;  fuj[0] = 0;
			}
			b[YM-1] = -1; c[YM-1] = 1; d[YM-1] = 0; fuj[YM - 1] = 0;

			for (int j = 1; j < YM - 1; j++){
				b[j] = tau *(-Nu / (dy * dy) - 1 / (4 * dy) * (VN[i][j - 1] + VN[i + 1][j - 1]));
				c[j] = 1 + 2 * tau * Nu / (dy * dy) + tau / (4 * dy) * (VN[i][j] + VN[i + 1][j] - VN[i][j - 1] - VN[i + 1][j - 1]);
				d[j] = tau *(-Nu / (dy * dy) + 1 / (4 * dy) * (VN[i][j] + VN[i + 1][j]));
				fuj[j] = U_N[i][j];
			}

			tempuj = progonka(b, c, d, fuj);
			for (int j = 0; j < YM; j++)
				UN1[i][j] = tempuj[j];
		}

		for (int j = 0; j < YM; j++)
			UN1[XM - 2][j] = UN1[XM - 3][j];
			

		//V
		// Прогонка вдоль
		for (int j = 1; j < YM - 2; j++){
			b[0] = 0; c[0] = 1; d[0] = 1; fvi[0] = 0;
			b[XM - 1] = -1; c[XM - 1] = 1; d[XM - 1] = 0; fvi[XM - 1] = 0;
			
			for (int i = 1; i < XM - 1; i++){
				b[i] = tau *(-Nu / (dx * dx) - 1 / (4 * dx) * (UN1[i - 1][j + 1] + UN1[i - 1][j]));
				c[i] = 1 + 2 * tau * Nu / (dx * dx) + tau * 1 / (4 * dx) * (UN1[i][j + 1] + UN1[i][j] - UN1[i - 1][j + 1] - UN1[i - 1][j]);
				d[i] = tau *(-Nu / (dx * dx) + 1 / (4 * dx) * (UN1[i][j + 1] + UN1[i][j]));
				fvi[i] = VN[i][j];
			}
	
			tempvi = progonka(b, c, d, fvi);
			for (int i = 0; i < XM; i++)
				V_N[i][j] = tempvi[i];
		}

		for (int i = 0; i < XM; i++)
			V_N[i][YM - 2] = V_N[i][YM - 3];

		VN1 = V_N;

		// Прогонка поперек
		for (int i = 1; i < XM-1; i++){
			b[0] = 0; c[0] = 1; d[0] = 0; fvj[0] = 0;
			b[YM - 2] = -1; c[YM - 2] = 1; d[YM - 2] = 0; fvj[YM - 2] = 0;
			
			for (int j = 1; j < YM - 2; j++){
				b[j] = tau *(-Nu / (dy * dy) - A*tau / (Ro * dy * dy) - V_N[i][j] / dy);
				c[j] = 1 + 2 * tau * Nu / (dy * dy) + 2 * A * tau * tau / (Ro * dy * dy);
				d[j] = tau *(-Nu / (dy * dy) - A*tau / (Ro * dy * dy) + V_N[i][j] / dy);
				fvj[j] = V_N[i][j] - tau / (Ro * dy) * (P_N[i][j + 1] - P_N[i][j]);
			}

			tempvj = progonka(b, c, d, fvj);
			for (int j = 0; j < YM - 1; j++)
				VN1[i][j] = tempvj[j];
		}
		// Расчет давления
		for (int i = 1; i < XM - 1; i++)
			for (int j = 1; j < YM - 1; j++)
				PN1[i][j] = P_N[i][j] - tau*A*(VN1[i][j] - VN1[i][j - 1]) / dy;

		for (int j = 0; j < YM-1; j++)
			VN1[XM - 1][j] = VN1[XM - 2][j];

		UN = UN1;
		VN = VN1;
		PN = PN1;
	}

	for (int i = 1; i < XM - 2; i++)
		if (i * dx >= L1 && i * dx <= L2)
			f1out << mu*(-3 * UN1[i][0] + 4 * UN1[i][1] - UN1[i][2]) / 2 / dy << endl;	//mu*(UN1[i][1] - UN1[i][0]) / dy << endl;

	fout << "Variables = \"x\",\"y\", \"U\", \"V\" " << endl;
	fout << "Zone i=" << XM - 1 << " j=" << YM - 1 << endl;
	fout.setf(ios::scientific);
	for (int j = 0; j < YM - 1; j++){
		for (int i = 0; i < XM - 1; i++){
			fout << dx * i << "  ";
			fout << dy * j << "  ";
			fout << (UN1[i][j + 1] + UN1[i][j]) / 2 << "  ";
			fout << (VN1[i + 1][j] + VN1[i][j]) / 2 << endl;
		}
	}

	f2out << "Variables = \"x\",\"y\", \"p\" " << endl;
	f2out << "Zone i=" << XM - 2 << " j=" << YM - 2 << endl;
	f2out.setf(ios::scientific);
	for (int j = 1; j < YM - 1; j++){
		for (int i = 1; i < XM - 1; i++){
			f2out << dx * i << "  ";
			f2out << dy * j << "  ";
			f2out << PN1[i][j] << endl;
		}
	}

	system("pause");

	return 0;
}
