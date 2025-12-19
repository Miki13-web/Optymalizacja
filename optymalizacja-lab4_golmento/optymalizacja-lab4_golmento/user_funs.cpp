#include"user_funs.h"

#define M_PI 3.14159265358979323846

matrix ff0T(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla przypadku testowego
{
	matrix y;												// y zawiera wartość funkcji celu
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);		// ud1 zawiera współrzędne szukanego optimum
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla problemu rzeczywistego
{
	matrix y;												// y zawiera wartość funkcji celu
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki początkowe
		MT = matrix(2, new double[2] { m2d(x), 0.5 });		// MT zawiera moment siły działający na wahadło oraz czas działania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);	// rozwiązujemy równanie różniczkowe
	int n = get_len(Y[0]);									// długość rozwiązania
	double teta_max = Y[1](0, 0);							// szukamy maksymalnego wychylenia wahadła
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));							// wartość funkcji celu (ud1 to założone maksymalne wychylenie)
	Y[0].~matrix();											// usuwamy z pamięci rozwiązanie RR
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);										// definiujemy wektor pochodnych szukanych funkcji
	double m = 1, l = 0.5, b = 0.5, g = 9.81;				// definiujemy parametry modelu
	double I = m * pow(l, 2);
	dY(0) = Y(1);																// pochodna z położenia to prędkość
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;	// pochodna z prędkości to przyspieszenie
	return dY;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2) {
	return -cos(0.1 * x(0)) * exp(-pow((0.1 * x(0) - 2 * 3.14), 2)) + 0.002 * pow((0.1 * x(0)), 2);
}

matrix ff1R(matrix x, matrix ud1, matrix ud2) {
	double DA = x(0) * 1e-4; // cm^2 na m^2

	// Warunki początkowe
	matrix Y0(3, 1);
	Y0(0) = 5.0; // VA początkowe
	Y0(1) = 1.0; // VB początkowe
	Y0(2) = 20.0; // TB początkowe

	matrix params(1, 1);
	params(0) = DA; // aktualne DA

	// Rozwiąż równanie różniczkowe
	matrix* Y = solve_ode(dff1R, 0, 1, 2000, Y0, NAN, params);

	// Znajdź temperaturę maksymalną w zbiorniku B
	int n = get_len(Y[0]);
	double T_max = 0;
	for (int i = 0; i < n; ++i) {
		if (Y[1](i, 2) > T_max) {
			T_max = Y[1](i, 2);
		}
	}

	// Zwolnienie pamięci
	delete[] Y;

	return abs(T_max - 50.0);
}

matrix dff1R(double t, matrix Y, matrix ud1, matrix ud2) {
	// Wartości początkowe
	double PA = 2.0; 
	double PB = 1.0; 
	double g = 9.81;
	double a = 0.98;
	double b = 0.63;
	double TA_in = 95.0;
	double TB_in_temp = 20.0;
	double Fin_B = 10.0 / 1000.0;
	double DB = 36.5665 * 1e-4;

	double DA = ud2(0); 

	// Aktualne stany
	double VA = Y(0);
	double VB = Y(1);
	double TB = Y(2);

	double hA = VA / PA;
	double hB = VB / PB;

	double Fout_A = (hA > 0) ? a * b * DA * sqrt(2 * g * hA) : 0;
	double Fout_B = (hB > 0) ? a * b * DB * sqrt(2 * g * hB) : 0;

	// Równania różniczkowe
	matrix dY(3, 1);

	dY(0) = -Fout_A; // dVA
	dY(1) = Fout_A + Fin_B - Fout_B; // dVB

	if (abs(VB) < 1e-9) { // Precyzja
		dY(2) = 0;
	}
	else {
		dY(2) = (Fout_A * (TA_in - TB) + Fin_B * (TB_in_temp - TB)) / VB; // dt
	}

	return dY;
}

matrix df1(double t, matrix Y, matrix ud1, matrix ud2) {
	matrix dY(Y);

	double a = 0.98;
	double b = 0.63;
	double g = 9.81;
	double PA = 2;
	double PB = 1;
	double DB = 0.00365665;
	double Fin = 0.01;
	double Tin = 20;
	double TA = 95.0;
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

matrix ff2T(matrix x, matrix ud1, matrix ud2) {
	double x1 = x(0);
	double x2 = x(1);
	double PI = acos(-1);

	double result = pow(x1, 2) + pow(x2, 2) - cos(2.5 * PI * x1) - cos(2.5 * PI * x2) + 2;

	return matrix(1, 1, result);
}

matrix ff2R(matrix x, matrix ud1, matrix ud2)
{
	// Zmienne optymalizacji
	// x(0) = k1
	// x(1) = k2

	// --- Ustawienia symulacji ---

	// Parametry czasu symulacji
	double t0 = 0.0;
	double tend = 100.0;
	double dt = 0.1;

	// Warunki początkowe (ramię w spoczynku w pozycji 0)
	matrix Y0(2, 1);
	Y0(0) = 0.0; // Kąt początkowy alpha(0)
	Y0(1) = 0.0; // Prędkość kątowa początkowa omega(0)

	// --- Przygotowanie parametrów dla df2 ---

	// ud1_sim: Wartości referencyjne (zadane)
	matrix ud1_sim(2, 1);
	ud1_sim(0) = 3.14; //1592653589793; // alpha_ref = pi
	ud1_sim(1) = 0.0;                  // omega_ref = 0

	// ud2_sim: Współczynniki wzmocnienia (zmienne optymalizacji)
	matrix ud2_sim = x;

	// --- Uruchomienie symulacji ---
	// Rozwiązanie równania różniczkowego
	matrix* Y_solution = solve_ode(df2, t0, dt, tend, Y0, ud1_sim, ud2_sim);

	matrix T = Y_solution[0]; // Wektor czasu
	matrix S = Y_solution[1]; // Wektor stanów [alpha, omega]
	int n = get_len(T);       // Liczba kroków czasowych

	// --- Obliczenie funkcjonału jakości Q --- 
	double Q = 0.0;

	// Pobranie wartości dla czytelności
	double k1 = x(0);
	double k2 = x(1);
	double alpha_ref = ud1_sim(0);
	double omega_ref = ud1_sim(1);

	// Całkowanie metodą prostokątów
	for (int i = 0; i < n - 1; ++i)
	{
		double alpha_t = S(i, 0); // alpha(t) w kroku i
		double omega_t = S(i, 1); // omega(t) w kroku i

		// Obliczenie kwadratów błędów
		double e_alpha_sq = pow(alpha_ref - alpha_t, 2);
		double e_omega_sq = pow(omega_ref - omega_t, 2);

		// Obliczenie M(t) w kroku i
		double M_t = k1 * (alpha_ref - alpha_t) + k2 * (omega_ref - omega_t);
		double M_sq = pow(M_t, 2);

		// Wartość funkcji podcałkowej 
		double integrand = 10.0 * e_alpha_sq + e_omega_sq + M_sq;

		// Sumowanie do całki: Q = integral(integrand) dt ≈ sum(integrand * dt)
		Q += integrand * dt;
	}

	// Zwolnienie pamięci
	Y_solution[0].~matrix();
	Y_solution[1].~matrix();

	// Zwrócenie wartości Q jako macierzy 1x1
	return matrix(1, 1, Q);
}

matrix df2(double t, matrix Y, matrix ud1, matrix ud2) {
	// Y(0) = alpha (kąt)
	// Y(1) = omega (prędkość kątowa)
	matrix dY(2, 1);

	// Stałe modelu
	double l = 2.0;     // długość ramienia
	double m_r = 1.0;   // masa ramienia
	double m_c = 5.0;   // masa ciężarka
	double b = 0.25;    // współczynnik tarcia

	// Obliczenie momentu bezwładności I
	double I = (1.0 / 3.0) * m_r * pow(l, 2) + m_c * pow(l, 2);

	// Wartości referencyjne (zadane) z ud1
	double alpha_ref = ud1(0);
	double omega_ref = ud1(1);

	// Współczynniki wzmocnienia (zmienne optymalizacji) z ud2
	double k1 = ud2(0);
	double k2 = ud2(1);

	// Aktualny stan
	double alpha_t = Y(0);
	double omega_t = Y(1);

	// Obliczenie momentu siły M(t)
	double M = k1 * (alpha_ref - alpha_t) + k2 * (omega_ref - omega_t);

	// Równania dynamiki systemu
	dY(0) = omega_t; // d(alpha)/dt = omega
	dY(1) = (M - b * omega_t) / I; // d(omega)/dt = (M - b*omega) / I

	return dY;
}

matrix ff3T_outside(matrix x, matrix a, matrix ud2)
{
	double x1 = x(0);
	double x2 = x(1);
	double aa = a(0);
	double c = ud2(0); // współczynnik kary

	// Funkcja celu
	double r = sqrt((x1*x1 + x2*x2) / (M_PI*M_PI));
	double f = sin(M_PI * r) / (M_PI * r);

	// Ograniczenia
	double g1 = -x1 + 1.0;
	double g2 = -x2 + 1.0;
	double g3 = sqrt(x1*x1 + x2*x2) - aa;

	double S =
		pow(std::max(0.0, g1), 2) +
		pow(std::max(0.0, g2), 2) +
		pow(std::max(0.0, g3), 2);

	return matrix(f + c * S); // kara
}

matrix ff3T_inside(matrix x, matrix a, matrix ud2)
{
	double x1 = x(0);
	double x2 = x(1);
	double aa = a(0);
	double c = ud2(0); // współczynnik kary

	double r = sqrt((x1*x1 + x2*x2) / (M_PI*M_PI));
	double f = sin(M_PI * r) / (M_PI * r);

	double g1 = -x1 + 1.0;
	double g2 = -x2 + 1.0;
	double g3 = sqrt(x1*x1 + x2*x2) - aa;

	// Jeżeli punkt narusza ograniczenia -> bardzo duża kara
	if (g1 >= 0 || g2 >= 0 || g3 >= 0)
		return matrix(1e12);

	// Wewnętrzna funkcja kary
	double eps = 1e-12; // Tolerancja
	double S = -(1.0/(g1 - eps) + 1.0/(g2 - eps) + 1.0/(g3 - eps));

	return matrix(f + c * S);
}

matrix df3(double t, matrix Y, matrix ud1, matrix omega_m)
{
	double C   = ud1(0);
	double rho = ud1(1);
	double r   = ud1(2);
	double m   = ud1(3);
	double g   = ud1(4);

	double x  = Y(0);
	double vx = Y(1);
	double y  = Y(2);
	double vy = Y(3);

	double omega = omega_m(0);

	double S = M_PI * r * r;

	// Siły oporu
	double Dx = 0.5 * C * rho * S * vx * fabs(vx);
	double Dy = 0.5 * C * rho * S * vy * fabs(vy);

	// Siły Magnusa
	double FMx = rho * vy * omega * M_PI * pow(r,3);
	double FMy = rho * vx * omega * M_PI * pow(r,3);

	matrix dY(4,1);
	dY(0) = vx;
	dY(1) = -(Dx + FMx)/m;
	dY(2) = vy;
	dY(3) = -(Dy + FMy)/m - g;

	return dY;
}

matrix ff3R(matrix x, matrix ud1, matrix ud2)
{
	double v0x = x(0);
	double omega = x(1);
	double c = ud2(0); // Współczynnik kary

	double S = 0.0;

	// Kary za przekroczenie zakresu zmiennych decyzyjnych [-10, 10]
	if (v0x < -10) S += pow(-10 - v0x, 2);
	if (v0x > 10) S += pow(v0x - 10, 2);
	if (omega < -10) S += pow(-10 - omega, 2);
	if (omega > 10) S += pow(omega - 10, 2);

	matrix Y0(4, new double[4] {0.0, v0x, 100.0, 0.0});
	matrix omega_m(1, new double[1] {omega});

	matrix* Y = solve_ode(df3, 0.0, 0.01, 7.0, Y0, ud1, omega_m);

	int N = get_len(Y[0]);

	// Znalezienie x_end (moment uderzenia w ziemię y=0)
	double x_end = 0.0;
	for (int i = 0; i < N - 1; i++) {
		double y_curr = Y[1](i, 2);
		double y_next = Y[1](i + 1, 2);

		// Wykrycie przejścia przez 0
		if (y_curr >= 0 && y_next < 0) {
			// Interpolacja liniowa dla dokładniejszego x_end
			double frac = (y_curr - 0.0) / (y_curr - y_next);
			double x_curr = Y[1](i, 0);
			double x_next = Y[1](i + 1, 0);
			x_end = x_curr + frac * (x_next - x_curr);
			break;
		}
	}

	bool found_50 = false;
	for (int i = 0; i < N - 1; i++) {
		double y_curr = Y[1](i, 2);
		double y_next = Y[1](i + 1, 2);

		// Sprawdzamy, czy w tym kroku przekroczyliśmy 50m
		if (y_curr >= 50.0 && y_next < 50.0) {
			// Interpolacja liniowa, aby znaleźć x dokładnie na wysokości 50m
			double frac = (y_curr - 50.0) / (y_curr - y_next);
			double x_curr = Y[1](i, 0);
			double x_next = Y[1](i + 1, 0);
			double x_at_50 = x_curr + frac * (x_next - x_curr);

			// Naliczanie kary, jeśli x nie jest w przedziale [3, 7]
			if (x_at_50 < 3) S += pow(3 - x_at_50, 2);
			if (x_at_50 > 7) S += pow(x_at_50 - 7, 2);

			found_50 = true;
			break;
		}
	}

	if (!found_50) S += 1000.0;

	delete[] Y;

	// Wynik
	return matrix(-x_end + c * S);
}

// Funkcja celu testowa: f(x1, x2) = 1/6*x1^6 - 1.05*x1^4 + 2*x1^2 + x2^2 + x1*x2
matrix ff4T(matrix x, matrix ud1, matrix ud2)
{
	double x1 = x(0);
	double x2 = x(1);
	double y = (1.0 / 6.0) * pow(x1, 6) - 1.05 * pow(x1, 4) + 2.0 * pow(x1, 2) + pow(x2, 2) + x1 * x2;
	return matrix(1, 1, y);
}

// Gradient funkcji testowej
matrix gf4T(matrix x, matrix ud1, matrix ud2)
{
	matrix g(2, 1);
	double x1 = x(0);
	double x2 = x(1);

	// Pochodna po x1: x1^5 - 4.2*x1^3 + 4*x1 + x2
	g(0) = pow(x1, 5) - 4.2 * pow(x1, 3) + 4.0 * x1 + x2;

	// Pochodna po x2: 2*x2 + x1
	g(1) = 2.0 * x2 + x1;

	return g;
}

// Hesjan funkcji testowej
matrix hf4T(matrix x, matrix ud1, matrix ud2)
{
	matrix H(2, 2);
	double x1 = x(0);
	double x2 = x(1);

	// d^2 f / dx1^2 = 5*x1^4 - 12.6*x1^2 + 4
	H(0, 0) = 5.0 * pow(x1, 4) - 12.6 * pow(x1, 2) + 4.0;

	// d^2 f / dx1 dx2 = 1
	H(0, 1) = 1.0;

	// d^2 f / dx2 dx1 = 1
	H(1, 0) = 1.0;

	// d^2 f / dx2^2 = 2
	H(1, 1) = 2.0;

	return H;
}

// --- Funkcje problemu rzeczywistego (Lab 4 - Klasyfikacja) ---

// Funkcja pomocnicza sigmoid
double sigmoid(double z) {
	return 1.0 / (1.0 + exp(-z));
}

// Funkcja kosztu logistycznego
// x - wektor theta (parametry)
// ud1 - macierz X (dane wejściowe, 3 wiersze, m kolumn)
// ud2 - wektor Y (etykiety, 1 wiersz, m kolumn)
matrix ff4R(matrix x, matrix ud1, matrix ud2)
{
	matrix X = ud1;
	matrix Y = ud2; // Przyjmujemy Y jako wektor poziomy (zgodnie z formatem w pliku)
	int m = get_len(trans(Y)); // Liczba próbek (liczba kolumn)
	double J = 0.0;

	for (int i = 0; i < m; ++i)
	{
		// Pobranie i-tej próbki (kolumny) z X
		matrix xi = X[i];
		double yi = Y(0, i); // Etykieta i-tej próbki

		// Obliczenie hipotezy h(x) = sigmoid(theta^T * x)
		double z = (trans(x) * xi)(0);
		double h = sigmoid(z);

		// Zabezpieczenie logarytmu przed NaN/Inf
		if (h < 1e-15) h = 1e-15;
		if (h > 1.0 - 1e-15) h = 1.0 - 1e-15;

		J += yi * log(h) + (1.0 - yi) * log(1.0 - h);
	}

	return matrix(1, 1, -J / m);
}

// Gradient funkcji kosztu logistycznego
matrix gf4R(matrix x, matrix ud1, matrix ud2)
{
	matrix X = ud1;
	matrix Y = ud2;
	int m = get_len(trans(Y));
	matrix grad(3, 1, 0.0); // Wektor gradientu o wymiarze takim jak theta (3x1)

	for (int i = 0; i < m; ++i)
	{
		matrix xi = X[i];
		double yi = Y(0, i);

		double z = (trans(x) * xi)(0);
		double h = sigmoid(z);

		// Sumowanie (h - y) * xi
		grad = grad + (h - yi) * xi;
	}

	return grad / (double)m;
}