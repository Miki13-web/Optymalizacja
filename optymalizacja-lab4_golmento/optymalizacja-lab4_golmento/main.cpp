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
#define PRAKTYCZNE4

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();
bool check_constraints(matrix x, double a);
matrix read_matrix_from_file(int rows, int cols, std::string filename);

int main()
{
	try
	{
		lab4();
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
	cout << "Optymalizacja - Projekt nr 1" << endl;
	cout << "Jakub Bednarczyk" << endl;
	cout << "Mikolaj Bednarczyk" << endl;
	cout << "Piotr Drozynski" << endl;
	cout << "Kacper Golmento" << endl;

	// Zadanie teoretyczne

	// Wartości początkowe
	double* res = new double[2] { 0, 0 };
	double x0 = -100.0, d = 5, alpha = 1.3, epsilon = 0.0001, gamma = 0.000001;
	int Nmax = 10000;
	//double a = 50, b = 70;

	solution wynik;

	std::ofstream filetest("wyniki_lab1_zadtest.csv");
	filetest << "x0;a;b;Expansion_f_calls;Fibonacci_x;Fibonacci_y;Fibonacci_f_calls;Lagrange_x;Lagrange_y;Lagrange_f_calls\n";

	for (int i = 0; i < 3; i++) {
		x0 = -100.0;
		filetest << "Alpha=" << alpha << "\n";
		for (int j = 0; j < 100; j++) {
			cout << "Zadanie teoretyczne - iteracja " << i + 1 << ":" << endl;
			cout << "x0: " << x0 << endl;
			res = expansion(ff1T, x0, d, alpha, Nmax);
			cout << "Metoda ekspansji:" << endl;
			cout << "a: " << res[0] << endl << "b: " << res[1] << endl << "Liczba wywolan funkcji: " << solution::f_calls << endl;
			cout << "b - a = " << res[1] - res[0] << endl << endl;
			filetest << x0 << ";" << res[0] << ";" << res[1] << ";" << solution::f_calls << ";";

			cout << "Metoda Fibonacciego:" << endl;
			wynik = fib(ff1T, res[0], res[1], epsilon);
			cout << wynik << endl;
			filetest << wynik.x << ";" << wynik.y << ";" << wynik.f_calls << ";";

			cout << "Metoda Lagrange'a:" << endl;
			wynik = lag(ff1T, res[0], res[1], epsilon, gamma, Nmax);
			cout << wynik << endl;
			filetest << wynik.x << ";" << wynik.y << ";" << wynik.f_calls << "\n";

			x0 += 2.0;
		}
		alpha += 0.5;
	}

	/*
	cout << "Zadanie teoretyczne - bez metody ekspansji:" << endl;

	wynik = fib_cout(ff1T, -100.0, 100.0, epsilon);
	cout << "Metoda Fibonacciego:" << endl;
	cout << wynik << endl;
	//filetest << "Bez ekspansji;-100;100;;" << wynik.x << ";" << wynik.y << ";" << wynik.f_calls << ";";

	cout << "Metoda Lagrange'a:" << endl;
	wynik = lag_cout(ff1T, -100.0, 100.0, epsilon, gamma, Nmax);
	cout << "Metoda Lagrange'a:" << endl;
	cout << wynik << endl;
	//filetest << wynik.x << ";" << wynik.y << ";" << wynik.f_calls << "\n";
	*/
	/*
	// Zadanie praktyczne
	double* res1 = new double[2] { 0, 0 };
	double da = 0.005, delta_da = 0.002, alpha1 = 1.5, epsilon1 = 0.0001, gamma1 = 0.000001;
	int nmax = 1000;
	cout << "----------------------------------------" << endl;
	cout << "Zadanie praktyczne:" << endl;
	res1 = expansion(ff1R, da, delta_da, alpha1, nmax);
	cout << "Metoda ekspansji:" << endl;
	cout << "a: " << res1[0] << endl << "b: " << res1[1] << endl << "Liczba wywolan funkcji: " << solution::f_calls << endl << endl;

	wynik = fib(ff1R, res1[0], res1[1], epsilon1);
	cout << "Metoda Fibonacciego: " << endl;
	cout << "Optymalna wielkosc otworu DA: " << wynik.x << endl;
	cout << "Maksymalna temperatura wody w zbiorniku: " << wynik.y + 50 << endl; // "| " << wynik.y << endl;
	cout << "Liczba wywolan funkcji: " << wynik.f_calls << endl;
	cout << "Exit flag: " << wynik.flag << endl;

	wynik = lag(ff1R, res1[0], res1[1], epsilon1, gamma1, nmax);
	cout << "Metoda Lagrange'a: " << endl;
	cout << "Optymalna wielosc otworu DA: " << wynik.x << endl;
	cout << "Maksymalna temperatura wody w zbiorniku: " << wynik.y + 50 << endl; // << " | " << wynik.y << endl;
	cout << "Liczna wywolan funkcji: " << wynik.f_calls << endl;
	cout << "Exit flag: " << wynik.flag << endl;

	// Symulacja

	// Warunki początkowe
	matrix Y0(3, 1);
	Y0(0) = 5.0;   // VA początkowe
	Y0(1) = 1.0;   // VB początkowe
	Y0(2) = 20.0;  // TB początkowe

	// Czas symulacji
	double t0 = 0.0;            // t początkowe
	double tend = 2000.0;       // t końcowe
	double dt = 1.0;            // krok

	matrix ud1(1, 1), ud2;
	ud1(0) = m2d(wynik.x) * 0.0001;

	matrix* S = solve_ode(df1, t0, dt, tend, Y0, ud1, ud2);
	int N = static_cast<int>(floor((tend - t0) / dt) + 1);
	std::vector<double> t_values(N);
	std::vector<double> VA_values(N);
	std::vector<double> VB_values(N);
	std::vector<double> TB_values(N);

	for (int i = 0; i < N; ++i) {
		t_values[i] = S[0](i);
		VA_values[i] = S[1](i, 0); // VA
		VB_values[i] = S[1](i, 1); // VB
		TB_values[i] = S[1](i, 2); // TB
	}

	std::ofstream file("symulacja_lab1.csv");
	file << "Czas;VA;VB;TB\n";
	for (int i = 0; i < N; ++i) {
		file << t_values[i] << ";" << VA_values[i] << ";" << VB_values[i] << ";" << TB_values[i] << "\n";
	}
	file.close();
	*/
}

void lab2()
{
	cout << "Optymalizacja - Projekt nr 2" << endl;
	cout << "Jakub Bednarczyk" << endl;
	cout << "Mikolaj Bednarczyk" << endl;
	cout << "Piotr Drozynski" << endl;
	cout << "Kacper Golmento" << endl;

	srand(time(NULL));

	// Pliki wynikowe
	std::ofstream SoutTeoria("wyniki_teoretyczne.csv");
	std::ofstream SoutSymulacja("symulacja_lab2.csv");

	// Nagłówki kolumn
	SoutTeoria << "Typ;Start_x1;Start_x2;Opt_x1;Opt_x2;Min_f;Liczba_wywolan\n";
	SoutSymulacja << "Czas;Kat;Predkosc_katowa\n";

	matrix X;
	double step = 0.041, alpha = 0.8, beta = 0.1, epsilon = 0.0001;
	double a, b;
	int Nmax = 1000;

	// zadanie teoretyczne
	cout << "Problem teoretyczny:" << endl;

	// inicjalizacja wartosci dla obu metod teoretycznych
	double a_vals[100], b_vals[100];
	for (int i = 0; i < 100; i++) {
		a_vals[i] = ((rand() % 200) / 100.0) - 1;
		b_vals[i] = ((rand() % 200) / 100.0) - 1;
	}

	// --- Najpierw HOOKE–JEEVES ---
	cout << "\n[Hooke-Jeeves]" << endl;
	for (int i = 0; i < 100; i++)
	{
		double a = a_vals[i];
		double b = b_vals[i];
		alpha = 0.8;
		X = matrix(2, new double[2]{ a, b });

		solution hooke = HJ(ff2T, X, step, alpha, epsilon, Nmax);
		cout << hooke;

		SoutTeoria << "Hooke-Jeeves;"
			<< a << ";" << b << ";"
			<< hooke.x(0) << ";" << hooke.x(1) << ";"
			<< hooke.y << ";" << hooke.f_calls << "\n";
	}

	// --- HOOKE–JEEVES z zapisem trajektorii dla jednego wybranego przypadku ---
	cout << "\n[Hooke-Jeeves - zapis trajektorii]" << endl;

	// wybierz dowolny punkt startowy z zakresu [-1,1], np. z pierwszej pętli:
	double a_traj = a_vals[0];
	double b_traj = b_vals[0];

	// utworzenie macierzy startowej
	matrix Xtraj(2, new double[2]{ a_traj, b_traj });

	// wywołanie funkcji zapisującej trajektorię do pliku CSV
	solution traj = HJ_with_trace(ff2T, Xtraj, step, alpha, epsilon, Nmax);

	cout << "Zapisano trajektorię do pliku: trajektoria_HJ.csv" << endl;

	// --- ROSENBROCK z zapisem trajektorii dla tego samego punktu startowego ---
	cout << "\n[Rosenbrock - zapis trajektorii]" << endl;

	// punkt startowy taki sam jak w HJ
	matrix XtrajRosen(2, new double[2]{ a_traj, b_traj });
	matrix StepRosen(2, new double[2]{ step, step });

	solution trajR = Rosen_with_trace(ff2T, XtrajRosen, StepRosen, alpha, beta, epsilon, Nmax);

	cout << "Zapisano trajektorię do pliku: trajektoria_Rosen.csv" << endl;

	//cout << traj;

	// --- Potem ROSENBROCK ---
	cout << "\n[Rosenbrock]" << endl;
	for (int i = 0; i < 100; i++)
	{
		double a = a_vals[i];
		double b = b_vals[i];
		alpha = 1.8;
		X = matrix(2, new double[2]{ a, b });
		matrix Step = matrix(2, new double[2]{ step, step });

		solution rosen = Rosen(ff2T, X, Step, alpha, beta, epsilon, Nmax);
		cout << rosen;

		SoutTeoria << "Rosenbrock;"
			<< a << ";" << b << ";"
			<< rosen.x(0) << ";" << rosen.x(1) << ";"
			<< rosen.y << ";" << rosen.f_calls << "\n";
	}

	SoutTeoria.close();

	//problem rzeczywisty
	cout << "\nProblem rzeczywisty:";

	//sprawdzenie poprawnosci
	matrix x(2, 1, 5); // macierz 2x1 wypelniona wartoscia 5
	cout << ff2R(x) << endl;

	X = matrix(2, new double[2]{ 5, 5 }); // tutaj chyba sa wartosci k
	solution wynikHJ = HJ(ff2R, x, step, alpha, epsilon, Nmax);
	cout << "\n[HJ] Wynik optymalizacji:" << endl;
	cout << "k1 = " << wynikHJ.x(0) << endl;
	cout << "k2 = " << wynikHJ.x(1) << endl;
	cout << "Minimalna wartosc funkcjonalu Q = " << wynikHJ.y << endl;
	cout << "Liczba wywolan funkcji: " << wynikHJ.f_calls << endl;
	cout << "Exit flag: " << wynikHJ.flag << endl;

	matrix Step = matrix(2, new double[2]{ step, step });
	solution wynikR = Rosen(ff2R, X, Step, alpha, beta, epsilon, Nmax);
	cout << "\n[Rosenbrock] Wynik optymalizacji:" << endl;
	cout << "k1 = " << wynikR.x(0) << endl;
	cout << "k2 = " << wynikR.x(1) << endl;
	cout << "Minimalna wartosc funkcjonalu Q = " << wynikR.y << endl;
	cout << "Liczba wywolan funkcji: " << wynikR.f_calls << endl;
	cout << "Exit flag: " << wynikR.flag << endl;

	//symulacja

	double t0 = 0.0;
	double tend = 100.0;
	double dt = 0.1;

	// Initial conditions for the state vector Y
	matrix Y0(2, 1);
	Y0(0) = 0.0;
	Y0(1) = 0.0;

	// Reference values for the desired angle and angular velocity
	matrix ud1(2, 1);
	ud1(0) = 3.14;
	ud1(1) = 0.0;

	// Gain parameters, k1 and k2, within the range [0, 10]
	matrix ud2(2, 1);
	//ud2(0) = wynikHJ.x(0);
	ud2(0) = wynikR.x(0);
	//ud2(1) = wynikHJ.x(1);
	ud2(1) = wynikR.x(1);

	matrix* result = solve_ode(df2, t0, dt, tend, Y0, ud1, ud2);

	int n = get_len(result[0]);
	std::cout << "\nTime\tAngle\tAngular Velocity" << std::endl;
	for (int i = 0; i < n; i++) {
		std::cout << result[0](i) << "\t"
			<< result[1](i, 0) << "\t"
			<< result[1](i, 1) << std::endl;
	}

	for (int i = 0; i < n; ++i) {
		SoutSymulacja << result[0](i) << ";" << result[1](i, 0) << ";" << result[1](i, 1) << "\n";
	}

	SoutSymulacja.close();

	cout << "\nWyniki zapisane do plikow:" << endl;
	cout << " - wyniki_teoretyczne.csv" << endl;
	cout << " - symulacja_lab2.csv" << endl;
}

void lab3()
{
    double epsilon = 1E-3;
    int Nmax = 10000;

    double c_inside = 1.0; 
    double dc_inside = 0.5; 

    double c_outside = 1.0; 
    double dc_outside = 2.0;

    std::ofstream Sout("symulacja_lab3.csv");

#ifdef TEORETYCZNE3
    Sout << "x0_1;x0_2;x1_out;x2_out;norm_out;y_out;f_calls_out;x1_in;x2_in;norm_in;y_in;f_calls_in\n";

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> x0_dist(1.0, 5.5); // Rozszerzony zakres losowania

    std::stringstream test_ss;
    solution test_sol;
    matrix a_val(1.0);
    matrix test_x0{};

    for (int i = 0; i < 3; ++i) {
        double current_a = 0.0;
        if (i == 0) current_a = 4.0;
        else if (i == 1) current_a = 4.4934;
        else current_a = 5.0;
        
        a_val = matrix(current_a);

        for (int j = 0; j < 100; ++j) {
            
            // Punkt musi być wewnątrz obszaru dopuszczalnego dla metody wewnętrznej
            do {
                test_x0 = matrix(2, new double[2] {x0_dist(gen), x0_dist(gen)});
            } while (!check_constraints(test_x0, current_a));

            test_ss << test_x0(0) << ";" << test_x0(1) << ";";

            // Zewnętrzne rozwiązanie
            solution::clear_calls();
            // true = 'is_external' (w pen)
            // lub sterujemy parametrami (dla zewnętrznej dc > 1)
            test_sol = pen(ff3T_outside, test_x0, c_outside, dc_outside, epsilon, Nmax, a_val);
            
            test_ss << test_sol.x(0) << ";" << test_sol.x(1) << ";"
                << sqrt(pow(test_sol.x(0), 2) + pow(test_sol.x(1), 2)) << ";"
                << test_sol.y(0) << ";" << test_sol.f_calls << ";";
            
            // Wewnętrzne rozwiązanie
            solution::clear_calls();
            // Przekazujemy poprawione c_inside
            test_sol = pen(ff3T_inside, test_x0, c_inside, dc_inside, epsilon, Nmax, a_val);
            
            test_ss << test_sol.x(0) << ";" << test_sol.x(1) << ";"
                << sqrt(pow(test_sol.x(0), 2) + pow(test_sol.x(1), 2)) << ";"
                << test_sol.y(0) << ";" << test_sol.f_calls << "\n";
        }
    }
    Sout << test_ss.str();
    Sout.close();
#endif

#ifdef PRAKTYCZNE3
    Sout << "time;x_position;y_position\n";
    
    // Dane zadania - bez zmian
    matrix ud1 = matrix(5, new double[5] {0.47, 1.2, 0.12, 0.6, 9.81});
    matrix x0 = matrix(2, new double[2] {-5.0, 5.0}); // Przykładowy punkt startowy

    // Optymalizacja zewnętrzna
    solution::clear_calls();
    solution opt = pen(ff3R, x0, c_outside, dc_outside, epsilon, Nmax, ud1);
    std::cout << "Optimum: " << opt.x(0) << ", " << opt.x(1) << " | y: " << opt.y(0) << "\n";

    // Symulacja zoptymalizowana
    matrix Y0(4, new double[4] {0.0, opt.x(0), 100.0, 0.0}); 
    matrix* Y = solve_ode(df3, 0.0, 0.01, 7.0, Y0, ud1, opt.x(1));

    int num_rows = get_len(Y[0]);
    for (int i = 0; i < num_rows; ++i) {
        Sout << Y[0](i) << ";" << Y[1](i, 0) << ";" << Y[1](i, 2) << "\n";
    }
    Sout.close();
    delete[] Y;
#endif
}

void lab4()
{
	try
	{
		// Parametry ogólne
		double epsilon = 1e-6;
		int Nmax = 10000;

		// ========================================================================
		// ZADANIE 5a: Testowa funkcja celu (Teoretyczne)
		// ========================================================================
#ifdef TEORETYCZNE4
		std::cout << "Rozpoczynam zadanie 5a (funkcja testowa)..." << std::endl;
		std::cout << "Generowanie 100 punktow startowych..." << std::endl;

		// 1. Pre-generowanie 100 punktów startowych
		std::vector<matrix> start_points;
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> dist(-2.0, 2.0);

		for (int i = 0; i < 100; ++i) {
			double x1 = dist(gen);
			double x2 = dist(gen);
			start_points.push_back(matrix(2, new double[2] {x1, x2}));
		}

		std::ofstream file_5a("lab4_zad5a.csv");
		if (!file_5a.is_open()) throw std::string("Nie udalo sie otworzyc pliku lab4_zad5a.csv");

		file_5a << "Metoda;Krok_Typ;h0;Iteracja;x0_1;x0_2;x_opt_1;x_opt_2;y_opt;f_calls;g_calls;H_calls;Status\n";

		double step_sizes[] = { 0.05, 0.25, 0.0 };
		std::string step_names[] = { "Staly 0.05", "Staly 0.25", "Zmienny" };

		// --- Metoda Najszybszego Spadku (SD) ---
		std::cout << "Obliczenia dla metody SD..." << std::endl;
		for (int k = 0; k < 3; ++k)
		{
			for (int i = 0; i < 100; ++i)
			{
				solution::clear_calls();
				solution res = SD(ff4T, gf4T, start_points[i], step_sizes[k], epsilon, Nmax);

				file_5a << "SD;" << step_names[k] << ";" << step_sizes[k] << ";" << (i + 1) << ";"
					<< start_points[i](0) << ";" << start_points[i](1) << ";"
					<< res.x(0) << ";" << res.x(1) << ";"
					<< res.y(0, 0) << ";"
					<< res.f_calls << ";" << res.g_calls << ";" << 0 << ";"
					<< res.flag << "\n";
			}
		}

		// --- Metoda Gradientów Sprzężonych (CG) ---
		std::cout << "Obliczenia dla metody CG..." << std::endl;
		for (int k = 0; k < 3; ++k)
		{
			for (int i = 0; i < 100; ++i)
			{
				solution::clear_calls();
				solution res = CG(ff4T, gf4T, start_points[i], step_sizes[k], epsilon, Nmax);

				file_5a << "CG;" << step_names[k] << ";" << step_sizes[k] << ";" << (i + 1) << ";"
					<< start_points[i](0) << ";" << start_points[i](1) << ";"
					<< res.x(0) << ";" << res.x(1) << ";"
					<< res.y(0, 0) << ";"
					<< res.f_calls << ";" << res.g_calls << ";" << 0 << ";"
					<< res.flag << "\n";
			}
		}

		// --- Metoda Newtona ---
		std::cout << "Obliczenia dla metody Newtona..." << std::endl;
		for (int k = 0; k < 3; ++k)
		{
			for (int i = 0; i < 100; ++i)
			{
				solution::clear_calls();
				solution res = Newton(ff4T, gf4T, hf4T, start_points[i], step_sizes[k], epsilon, Nmax);

				file_5a << "Newton;" << step_names[k] << ";" << step_sizes[k] << ";" << (i + 1) << ";"
					<< start_points[i](0) << ";" << start_points[i](1) << ";"
					<< res.x(0) << ";" << res.x(1) << ";"
					<< res.y(0, 0) << ";"
					<< res.f_calls << ";" << res.g_calls << ";" << res.H_calls << ";"
					<< res.flag << "\n";
			}
		}

		file_5a.close();
		std::cout << "Zakonczono zadanie 5a. Wyniki w pliku lab4_zad5a.csv" << std::endl;
#endif

		// ========================================================================
		// ZADANIE 5b: Problem rzeczywisty (Praktyczne)
		// ========================================================================
#ifdef PRAKTYCZNE4
		std::cout << "Rozpoczynam zadanie 5b (klasyfikacja)..." << std::endl;

		std::ofstream file_5b("lab4_zad5b.csv");
		if (!file_5b.is_open()) throw std::string("Nie udalo sie otworzyc pliku lab4_zad5b.csv");
		file_5b << "Krok(h0);theta_0;theta_1;theta_2;f_kosztu;f_calls;g_calls;Skutecznosc(%);Status\n";

		matrix X = read_matrix_from_file(3, 100, "XData.txt");
		matrix Y = read_matrix_from_file(1, 100, "YData.txt");

		std::cout << "\n--- WERYFIKACJA POPRAWNOSCI FUNKCJI ---" << std::endl;
		matrix theta_test(3, new double[3] {0.1, 0.1, 0.1});
		matrix J_test = ff4R(theta_test, X, Y);
		matrix G_test = gf4R(theta_test, X, Y);

		std::cout << "Dla theta = [0.1, 0.1, 0.1]:" << std::endl;
		std::cout << "J(theta) = " << J_test(0) << " (Oczekiwane: ~2.72715)" << std::endl;
		std::cout << "Gradient = \n" << G_test << "(Oczekiwane: ~0.299, ~13.60, ~13.35)" << std::endl;
		std::cout << "---------------------------------------\n" << std::endl;
		// ------------------------------------------------------------------------

		// Punkt startowy do WŁAŚCIWEJ optymalizacji musi być [0, 0, 0]
		matrix theta_start(3, new double[3] {0.0, 0.0, 0.0});

		double steps_real[] = { 0.01, 0.001, 0.0001 };

		for (double h : steps_real)
		{
			solution::clear_calls();

			// Tu używamy theta_start = [0,0,0]
			solution theta_opt = CG(ff4R, gf4R, theta_start, h, epsilon, Nmax, X, Y);

			// Obliczanie skuteczności klasyfikacji
			int correct_predictions = 0;
			int m = 100;

			for (int i = 0; i < m; ++i)
			{
				matrix xi = X[i];
				double yi = Y(0, i);
				double z = (trans(theta_opt.x) * xi)(0);
				double h_val = 1.0 / (1.0 + exp(-z));

				int prediction = (h_val >= 0.5) ? 1 : 0;
				if (prediction == (int)yi) correct_predictions++;
			}

			double accuracy = (double)correct_predictions / m * 100.0;

			std::cout << "Krok h=" << h << " | Koszt=" << theta_opt.y(0) << " | Acc=" << accuracy << "%" << std::endl;

			file_5b << h << ";"
				<< theta_opt.x(0) << ";" << theta_opt.x(1) << ";" << theta_opt.x(2) << ";"
				<< theta_opt.y(0) << ";"
				<< theta_opt.f_calls << ";" << theta_opt.g_calls << ";"
				<< accuracy << ";"
				<< theta_opt.flag << "\n";
		}

		file_5b.close();
		std::cout << "Zakonczono zadanie 5b. Wyniki w pliku lab4_zad5b.csv" << std::endl;
#endif

	}
	catch (std::string ex_info)
	{
		std::cerr << "WYSTAPIL BLAD: " << ex_info << std::endl;
	}
}

void lab5()
{

}

void lab6()
{

}

bool check_constraints(matrix x, double a) {
    // g1: -x1 + 1 <= 0  => x1 >= 1
    if (-x(0) + 1.0 > 0) return false;
    // g2: -x2 + 1 <= 0  => x2 >= 1
    if (-x(1) + 1.0 > 0) return false;
    // g3: sqrt(x1^2 + x2^2) - a <= 0
    double r = sqrt(pow(x(0), 2) + pow(x(1), 2));
    if (r - a > 0) return false;
    return true;
}

matrix read_matrix_from_file(int rows, int cols, std::string filename) {
	std::ifstream file(filename);
	if (!file.is_open()) {
		throw std::string("Nie mozna otworzyc pliku: " + filename);
	}

	matrix M(rows, cols);
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			double val;
			file >> val; // Wczytaj liczbę

			if (file.fail()) {
				// Jeśli wczytywanie się nie powiodło (np. koniec pliku), przerwij lub rzuć błąd
				break;
			}

			M(i, j) = val;

			// Pętla usuwająca separatory (średniki, przecinki) i białe znaki po liczbie
			while (file.peek() == ';' || file.peek() == ',' || isspace(file.peek())) {
				file.get(); // "Zjedz" znak separatora/spacji
			}
		}
	}
	file.close();
	return M;
}