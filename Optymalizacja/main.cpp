/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
Data ostatniej modyfikacji: 30.09.2025
*********************************************/

#include"opt_alg.h"
#include"user_funs.h"

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

matrix df1(double t, matrix Y, matrix ud1, matrix ud2) {
	matrix dY(Y);

	double a = 0.98;
	double b = 0.63;
	double g = 9.81;
	double PA = 0.5;
	double PB = 1;
	double DB = 0.00365665;
	double Fin = 0.01;
	double Tin = 20;
	double TA = 90.0;
	double DA = ud1(0);

	double VA = Y(0);
	double VB = Y(1);
	double TB = Y(2);

	double FAout = VA > 0 ? a * b * DA * sqrt(2 * g * VA / PA) : 0;
	double FBout = VB > 0 ? a * b * DB * sqrt(2 * g * VB / PB) : 0;

	dY(0) = -FAout;
	dY(1) = FAout + Fin - FBout;
	dY(2) = FAout / VB * (TA - TB) + Fin / VB * (Tin - TB);

	return dY;
}

int main()
{
	try
	{
		lab1();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	return 0;
}

void lab0()
{
	/*
	//Funkcja testowa
	double epsilon = 1e-2;									// dok�adno��
	int Nmax = 10000;										// maksymalna liczba wywo�a� funkcji celu
	matrix lb(2, 1, -5), ub(2, 1, 5),						// dolne oraz g�rne ograniczenie
		a(2, 1);											// dok�adne rozwi�zanie optymalne
	solution opt;											// rozwi�zanie optymalne znalezione przez algorytm
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);			// wywo�anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie licznik�w

	//Wahadlo
	Nmax = 1000;											// dok�adno��
	epsilon = 1e-2;											// maksymalna liczba wywo�a� funkcji celu
	lb = 0, ub = 5;											// dolne oraz g�rne ograniczenie
	double teta_opt = 1;									// maksymalne wychylenie wahad�a
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);		// wywo�anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie licznik�w

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz�tkowe
		MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);	// rozwi�zujemy r�wnanie r�niczkowe
	ofstream Sout("symulacja_lab0.csv");					// definiujemy strumie� do pliku .csv
	Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
	Sout.close();											// zamykamy strumie�
	Y[0].~matrix();											// usuwamy z pami�ci rozwi�zanie RR
	Y[1].~matrix();
	*/
}

void lab1()
{

	std::ofstream Sout("symulacja_lab1.csv");

	//zadanie teoretyczne

	double* res = new double[2] { 0, 0 };
	double x0 = 50, d = 5, alpha = 1.5;
	int Nmax = 10000;

	//double a = 50, b = 70;
	double epsilon = 0.0001;
	double gamma = 0.000001;
	solution wynik;

	for (int i = 0; i < 1; i++) {
		res = expansion(ff1T, x0, d, alpha, Nmax);
		cout << res[0] << endl << res[1] << endl << solution::f_calls << endl << endl;
		//Sout << "x" << res[0] << ";" << "x" << res[1] << ";" << "x" << solution::f_calls << "\n";

		wynik = fib(ff1T, res[0], res[1], epsilon);
		//Sout << "x" << wynik.x << "x" << wynik.y << "x" << wynik.f_calls << "\n";
		cout << wynik << endl;

		wynik = lag(ff1T, res[0], res[1], epsilon, gamma, Nmax);
		//Sout << "x" << wynik.x << "x" << wynik.y << "x" << wynik.f_calls << "\n";
		cout << wynik << endl;

		x0 = x0 + 2;
	}


	//zadanie praktyczne

	double* res1 = new double[2] { 0, 0 };
	double da = 0.005, delta_da = 0.002;
	double alpha1 = 1.5, epsilon1 = 0.0001, gamma1 = 0.000001;
	int nmax = 1000;

	res1 = expansion(ff1T, da, delta_da, alpha1, nmax);
	solution wynik2;
	//wynik = fib(ff2T, res[0], res[1], epsilon);
	cout << "Metoda Fibonacciego: " << endl;
	cout << "Optymalna wielosc otwou D_A: " << wynik.x << "\nMaksymalna temperatura wody w zbiorniku: " << wynik.y + 50 << "|" << wynik.y << "\nLiczna wywolan fukcji: " << wynik.f_calls << "\nExit flag: " << wynik.flag << endl;;
	//cout << wynik;
	wynik = lag(ff1T, res[0], res[1], epsilon1, gamma1, nmax);
	cout << "Metoda Lagrangea: " << endl;
	cout << "Optymalna wielosc otworu D_A: " << wynik.x << "\nMaksymalna temperatura wody w zbiorniku: " << wynik.y + 50 << "|" << wynik.y << "\nLiczna wywolan fukcji: " << wynik.f_calls << "\nExit flag: " << wynik.flag << endl;
	//cout << wynik;

	//symlacja

	// Warunki początkowe
	matrix Y0(3, 1);
	Y0(0) = 5.0;   // Początkowa objętość w zbiorniku A (VA)
	Y0(1) = 1.0;   // Początkowa objętość w zbiorniku B (VB)
	Y0(2) = 20.0;  // Początkowa temperatura w zbiorniku B (TB)

	// Czas symulacji
	double t0 = 0.0;            // Początkowy czas
	double tend = 2000.0;       // Końcowy czas symulacji
	double dt = 1.0;            // Krok czasowy

	matrix ud1(1, 1), ud2;
	ud1(0) = m2d(wynik.x);

	matrix* S = solve_ode(df1, t0, dt, tend, Y0, ud1, ud2);
	int N = static_cast<int>(floor((tend - t0) / dt) + 1);
	std::vector<double> t_values(N);
	std::vector<double> VA_values(N);
	std::vector<double> VB_values(N);
	std::vector<double> TB_values(N);

	for (int i = 0; i < N; ++i) {
		t_values[i] = S[0](i);
		VA_values[i] = S[1](i, 0); // Pierwsza kolumna to VA
		VB_values[i] = S[1](i, 1); // Druga kolumna to VB
		TB_values[i] = S[1](i, 2); // Trzecia kolumna to TB
	}

	std::ofstream file("symulacja_lab1.csv");
	file << "Czas;VA;VB;TB\n";
	for (int i = 0; i < N; ++i) {
		file << t_values[i] << "x;" << VA_values[i] << "x;" << VB_values[i] << "x;" << TB_values[i] << "x\n";
	}
	file.close();

}

void lab2()
{

}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}