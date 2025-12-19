#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	// Zmienne wej�ciowe:
	// ff - wska�nik do funkcji celu
	// N - liczba zmiennych funkcji celu
	// lb, ub - dolne i g�rne ograniczenie
	// epslion - zak��dana dok�adno�� rozwi�zania
	// Nmax - maksymalna liczba wywo�a� funkcji celu
	// ud1, ud2 - user data
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);									// losujemy macierz Nx1 stosuj�c rozk�ad jednostajny na przedziale [0,1]
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);// przeskalowywujemy rozwi�zanie do przedzia�u [lb, ub]
			Xopt.fit_fun(ff, ud1, ud2);							// obliczmy warto�� funkcji celu
			if (Xopt.y < epsilon)								// sprawdzmy 1. kryterium stopu
			{
				Xopt.flag = 1;									// flaga = 1 ozancza znalezienie rozwi�zanie z zadan� dok�adno�ci�
				break;
			}
			if (solution::f_calls > Nmax)						// sprawdzmy 2. kryterium stopu
			{
				Xopt.flag = 0;									// flaga = 0 ozancza przekroczenie maksymalne liczby wywo�a� funkcji celu
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2] { 0, 0 };
		solution::clear_calls();
		solution X0(x0);
		solution X1(x0 + d);
		X0.fit_fun(ff);
		X1.fit_fun(ff);
		vector<solution> x_vector;
		x_vector.push_back(X0);
		x_vector.push_back(X1);

		if (X0.y(0) == X1.y(0)) {
			p[0] = X0.x(0);
			p[1] = X1.x(0);
			return p;
		}

		if (X1.y(0) > X0.y(0)) {
			d = -d;
			X1.x(0) = X0.x(0) + d;
			X1.fit_fun(ff);
			x_vector[1] = X1;
			if (X1.y(0) >= X0.y(0)) {
				p[0] = X1.x(0);
				p[1] = X0.x(0) - d;
				return p;
			}
		}
		int i = 0;
		do
		{
			if (X0.f_calls > Nmax)
				break;
			i++;
			x_vector.push_back(x0 + (pow(alpha, i) * d));
			x_vector[i + 1].fit_fun(ff);
		} while (x_vector[i].y(0) >= x_vector[i + 1].y(0));
		if (d > 0)
		{
			p[0] = x_vector[i - 1].x(0);
			p[1] = x_vector[i + 1].x(0);
		}
		else {
			p[0] = x_vector[i + 1].x(0);
			p[1] = x_vector[i - 1].x(0);
		}
		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		Xopt.clear_calls();
		solution A(a), B(b), C, D;
		vector<double> ciag_fib;
		ciag_fib.push_back(1.0);
		ciag_fib.push_back(1.0);

		while (ciag_fib.back() <= (B.x(0) - A.x(0)) / epsilon) {
			ciag_fib.push_back(ciag_fib[ciag_fib.size() - 1] + ciag_fib[ciag_fib.size() - 2]);
		}

		int k = 1;
		while (true) {
			if (ciag_fib[k] > (B.x(0) - A.x(0)) / epsilon)
				break;
			k++;
		}
		C.x(0) = B.x(0) - (ciag_fib[k - 1] / ciag_fib[k]) * (B.x(0) - A.x(0));
		D.x(0) = A.x(0) + B.x(0) - C.x(0);
		C.fit_fun(ff);
		D.fit_fun(ff);

		for (int i = 0; i < k - 3; i++) {
			if (C.y(0) < D.y(0)) {
				B.x = D.x;
			}
			else {
				A.x = C.x;
			}
			C.x(0) = B.x(0) - (ciag_fib[k - i - 2] / ciag_fib[k - i - 1]) * (B.x(0) - A.x(0));
			D.x(0) = A.x(0) + B.x(0) - C.x(0);
			C.fit_fun(ff);
			D.fit_fun(ff);
		}

		Xopt.x = C.x(0);
		Xopt.y = C.y(0);
		Xopt.flag = 1;
		return Xopt;

	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}
}

solution fib_cout(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		Xopt.clear_calls(); // Resetuje global_f_calls_count
		solution A(a), B(b), C, D;
		std::vector<double> ciag_fib;
		ciag_fib.push_back(1.0);
		ciag_fib.push_back(1.0);
		std::ofstream fib_file("fib_tab2_output.csv");
		fib_file << "b-a\n";

		while (ciag_fib.back() <= (B.x(0) - A.x(0)) / epsilon) {
			ciag_fib.push_back(ciag_fib[ciag_fib.size() - 1] + ciag_fib[ciag_fib.size() - 2]);
		}

		int k = 1;
		while (true) {
			if (ciag_fib[k] > (B.x(0) - A.x(0)) / epsilon)
				break;
			k++;
		}

		// --- Inicjalizacja (2 ewaluacje) ---
		C.x(0) = B.x(0) - (ciag_fib[k - 1] / ciag_fib[k]) * (B.x(0) - A.x(0));
		D.x(0) = A.x(0) + B.x(0) - C.x(0);
		C.fit_fun(ff);
		D.fit_fun(ff);

		std::cout << "------------------------------------------" << std::endl;
		std::cout << "Początkowe wywołania funkcji: " << solution::f_calls << std::endl;
		std::cout << "a = " << A.x(0) << ", b = " << B.x(0) << std::endl;
		std::cout << "------------------------------------------" << std::endl;


		// Zgodnie z żądaniem, aby licznik iteracji zrównał się z f_calls,
		// inkrementujemy licznik iteracji o 2 w pętli i dodajemy początkowe 2 wywołania.
		for (int i = 0; i < k - 3; i++) {

			// Licznik f_calls (ewaluacji funkcji) po inicjalizacji
			int current_f_call = solution::f_calls;

			// Wypisanie numeru ewaluacji funkcji (f_calls) jako "iteracji"
			std::cout << "Ewaluacje f (f_calls) nr " << current_f_call + 1 << " i " << current_f_call + 2 << ":" << std::endl;
			std::cout << "Przed zawężeniem: a = " << A.x(0) << ", b = " << B.x(0) << std::endl;
			fib_file << B.x(0) - A.x(0) << "\n";

			// **Zawężanie przedziału**
			if (C.y(0) < D.y(0)) {
				B.x = D.x; // Nowe b staje się D
			}
			else {
				A.x = C.x; // Nowe a staje się C
			}

			// Wypisanie przedziału po zawężeniu
			std::cout << "Po zawężeniu:     a = " << A.x(0) << ", b = " << B.x(0) << std::endl;
			fib_file << B.x(0) - A.x(0) << "\n";

			// --- Ewaluacje w pętli (2 ewaluacje na krok) ---
			C.x(0) = B.x(0) - (ciag_fib[k - i - 2] / ciag_fib[k - i - 1]) * (B.x(0) - A.x(0));
			D.x(0) = A.x(0) + B.x(0) - C.x(0);
			C.fit_fun(ff);
			D.fit_fun(ff);

			std::cout << "Aktualna liczba f_calls: " << solution::f_calls << std::endl;
			std::cout << "------------------------------------------" << std::endl;
		}

		// Faza końcowa (ewaluacje dla ostatniego przedziału, 2 ewaluacje)
		if (C.y(0) < D.y(0)) {
			B.x = D.x;
		}
		else {
			A.x = C.x;
		}

		std::cout << "Ostateczny przedział: a = " << A.x(0) << ", b = " << B.x(0) << std::endl;
		std::cout << "Całkowita liczba wywołań funkcji (f_calls): " << solution::f_calls << std::endl;

		Xopt.x = C.x(0);
		Xopt.y = C.y(0);
		Xopt.flag = 1;
		fib_file.close();
		return Xopt;

	}
	catch (std::string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}
}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		Xopt.clear_calls();
		solution A(a), B(b), C, D, D0;
		C.x = (a + b) / 2;
		A.fit_fun(ff);
		B.fit_fun(ff);
		C.fit_fun(ff);
		double l, m;
		int i = 0;
		while (true) {
			l = A.y(0) * (pow(B.x(0), 2) - pow(C.x(0), 2)) + B.y(0) * (pow(C.x(0), 2) - pow(A.x(0), 2)) + C.y(0) * (pow(A.x(0), 2) - pow(B.x(0), 2));
			m = (A.y(0) * (B.x(0) - C.x(0))) + (B.y(0) * (C.x(0) - A.x(0))) + (C.y(0) * (A.x(0) - B.x(0)));
			if (m <= 0) {
				Xopt.flag = 0;
				break;
			}
			D0.x = D.x;
			D.x = l / (2 * m);
			D.fit_fun(ff);
			if (A.x(0) < D.x(0) && D.x(0) < C.x(0))
			{
				if (D.y(0) < C.y(0))
				{
					B.x = C.x;
					C.x = D.x;
					B.fit_fun(ff);
					C.fit_fun(ff);
				}
				else
					A.x = D.x;
				A.fit_fun(ff);
			}
			else if (C.x(0) < D.x(0) && D.x(0) < B.x(0))
			{
				if (D.y(0) < C.y(0))
				{
					A.x = C.x;
					C.x = D.x;
					A.fit_fun(ff);
					C.fit_fun(ff);
				}
				else
					B.x = D.x;
				B.fit_fun(ff);
			}
			else {
				Xopt.flag = 0;
				break;
			}
			if (i > Nmax) {
				Xopt.flag = 0;
				break;
			}
			if (B.x(0) - A.x(0) < epsilon || abs(D.x(0) - D0.x(0)) <= gamma)
			{
				break;
			}
		}
		D.fit_fun(ff);
		Xopt.x = D.x(0);
		Xopt.y = D.y(0);
		Xopt.flag = 1;

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution lag_cout(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		Xopt.clear_calls(); // Resetuje global_f_calls_count
		solution A(a), B(b), C, D, D0;
		C.x(0) = (a + b) / 2;
		std::ofstream lag_file("lag_tab2_output.csv");
		lag_file << "b-a\n";

		// --- Inicjalizacja (3 ewaluacje) ---
		A.fit_fun(ff);
		B.fit_fun(ff);
		C.fit_fun(ff);

		double l, m;
		int i = 0; // Licznik kroków iteracyjnych (nie f_calls)

		std::cout << "--------------------------------------------------------" << std::endl;
		std::cout << "Inicjalizacja (f_calls = 3):" << std::endl;
		std::cout << "Przedział: [a=" << A.x(0) << ", b=" << B.x(0) << "], c=" << C.x(0) << std::endl;
		std::cout << "--------------------------------------------------------" << std::endl;

		while (true) {
			i++;

			// Wypisanie numeru kroku iteracyjnego
			std::cout << "Krok iteracyjny nr " << i << " (Przed D.fit_fun):" << std::endl;

			l = A.y(0) * (pow(B.x(0), 2) - pow(C.x(0), 2)) + B.y(0) * (pow(C.x(0), 2) - pow(A.x(0), 2)) + C.y(0) * (pow(A.x(0), 2) - pow(B.x(0), 2));
			m = (A.y(0) * (B.x(0) - C.x(0))) + (B.y(0) * (C.x(0) - A.x(0))) + (C.y(0) * (A.x(0) - B.x(0)));

			if (m <= 0) {
				Xopt.flag = 0;
				std::cout << "Kryterium stopu: m <= 0. Niepowodzenie." << std::endl;
				break;
			}

			D0.x = D.x;
			D.x(0) = l / (2 * m); // Nowy punkt D (aproksymowane minimum)

			// 1. Ewaluacja D (1 ewaluacja)
			D.fit_fun(ff);

			std::cout << "  - Nowy D (aproksymowane minimum): " << D.x(0) << std::endl;
			std::cout << "  - f_calls po obliczeniu D: " << solution::f_calls << std::endl;
			std::cout << "  - Przed zawężeniem: [a=" << A.x(0) << ", b=" << B.x(0) << "], c=" << C.x(0) << std::endl;
			lag_file << B.x(0) - A.x(0) << "\n";

			if (A.x(0) < D.x(0) && D.x(0) < C.x(0)) // D jest w [A, C]
			{
				if (D.y(0) < C.y(0))
				{
					B.x = C.x;
					C.x = D.x;
					// 2. Ewaluacje: B.fit_fun jest zbędne, jeśli B=C (wartość y jest już znana),
					// ale oryginalny kod wywołuje fit_fun na B i C (czyli C=D)
					// Zgodnie z oryginalnym kodem:
					B.fit_fun(ff); // 2 ewaluacja w tym bloku
					C.fit_fun(ff); // 3 ewaluacja w tym bloku (ew. 2 jeśli C.x = D.x i D.y jest znane)
				}
				else
				{
					A.x = D.x;
					A.fit_fun(ff); // 2 ewaluacja w tym bloku
				}
				// Oryginalny kod, niestety, wywołuje A.fit_fun(ff) bez względu na warunek
				// w bloku, gdy D jest w [A, C]:
				// A.fit_fun(ff); // Ta linia w oryginalnym kodzie zepsuła by logikę!
				// Zamiast jej, usuwamy zbędne wywołanie A.fit_fun(ff) z końca, jeśli nie jest potrzebne:
			}
			else if (C.x(0) < D.x(0) && D.x(0) < B.x(0)) // D jest w [C, B]
			{
				if (D.y(0) < C.y(0))
				{
					A.x = C.x;
					C.x = D.x;
					A.fit_fun(ff); // 2 ewaluacja w tym bloku
					C.fit_fun(ff); // 3 ewaluacja w tym bloku
				}
				else
				{
					B.x = D.x;
					B.fit_fun(ff); // 2 ewaluacja w tym bloku
				}
				// B.fit_fun(ff); // Podobnie jak wyżej, to wywołanie zepsółoby logikę.
			}
			else {
				Xopt.flag = 0;
				std::cout << "Kryterium stopu: D poza [A, B] lub D=C. Niepowodzenie." << std::endl;
				break;
			}

			// Wypisanie przedziału po zawężeniu i f_calls
			std::cout << "  - Po zawężeniu:     [a=" << A.x(0) << ", b=" << B.x(0) << "], c=" << C.x(0) << std::endl;
			lag_file << B.x(0) - A.x(0) << "\n";
			std::cout << "  - Aktualna liczba f_calls: " << solution::f_calls << std::endl;


			if (i > Nmax) {
				Xopt.flag = 0;
				std::cout << "Kryterium stopu: Przekroczono maksymalną liczbę iteracji Nmax." << std::endl;
				break;
			}
			if (B.x(0) - A.x(0) < epsilon || std::abs(D.x(0) - D0.x(0)) <= gamma)
			{
				std::cout << "Kryterium stopu: Osiągnięto wymaganą dokładność (|b-a| < epsilon lub |D-D0| < gamma)." << std::endl;
				break;
			}
			std::cout << "--------------------------------------------------------" << std::endl;
		}

		// Finalna ewaluacja D (jest w oryginalnym kodzie poza pętlą, ale D.y jest już znane)
		D.fit_fun(ff);
		std::cout << "Całkowita liczba wywołań funkcji (f_calls): " << solution::f_calls << std::endl;

		Xopt.x = D.x(0);
		Xopt.y = D.y(0);
		Xopt.flag = 1;
		lag_file.close();

		return Xopt;
	}
	catch (std::string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		Xopt.clear_calls();
		solution XB, XB_temp, X(x0);
		do {
			XB.x = X.x;
			XB.fit_fun(ff);
			X = HJ_trial(ff, XB, s);
			if (X.y < XB.y) {
				do {
					XB_temp = XB;
					XB = X;
					X.x = 2 * XB.x - XB_temp.x;
					X.fit_fun(ff);
					X = HJ_trial(ff, X, s);
					if (X.f_calls > Nmax) {
						Xopt.x = XB.x;
						Xopt.y = XB.y;
						Xopt.flag = 0;
						return Xopt;
					}
				} while (X.y < XB.y);
				X = XB;
			}
			else {
				s = alpha * s;
			}
			if (X.f_calls > Nmax) {
				Xopt.x = XB.x;
				Xopt.y = XB.y;
				Xopt.flag = 0;
				return Xopt;
			}
		} while (s > epsilon);

		Xopt.x = XB.x;
		Xopt.y = XB.y;
		Xopt.flag = 1;

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		// Pobranie liczby wymiarów (zmiennych)
		int n = get_len(XB.x);
		// Utworzenie wektora bazowego e_j
		matrix e_j(n, 1, 0.0);
		// Stworzenie roboczej kopii rozwiązania
		solution X = XB;

		// Pętla po wszystkich wymiarach
		for (int j = 0; j < n; ++j)
		{
			// Ustawienie wektora bazowego e_j (1 na j-tej pozycji, reszta 0)
			e_j(j) = 1.0;

			// Utworzenie tymczasowego rozwiązania do testów
			solution X_test = X;

			// 1. Test w kierunku dodatnim: x + s * e^j
			X_test.x = X.x + s * e_j;
			X_test.fit_fun(ff, ud1, ud2);

			// Jeśli znaleziono lepszy punkt to zaktualizuj punkt bazowy dla tej próby
			if (X_test.y < X.y)
			{
				X = X_test;
			}
			else
			{
				// 2. Test w kierunku ujemnym: x - s * e^j
				X_test.x = X.x - s * e_j;
				X_test.fit_fun(ff, ud1, ud2);

				// Jeśli znaleziono lepszy punkt
				if (X_test.y < X.y)
				{
					// Zaktualizuj punkt bazowy dla tej próby
					X = X_test;
				}
			}
			// Jeśli żaden kierunek nie dał poprawy, X pozostaje bez zmian
		}
		// Zwróć najlepszy punkt znaleziony w tej próbie
		return X;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution HJ_with_trace(matrix(*ff)(matrix, matrix, matrix),
					   matrix x0, double s, double alpha,
					   double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		std::ofstream traceFile("trajektoria_HJ.csv");
		traceFile << "Iteracja;x1;x2;f\n";

		solution Xopt;
		Xopt.clear_calls();
		solution XB, XB_temp, X(x0);
		int iter = 0;

		do {
			XB.x = X.x;
			XB.fit_fun(ff, ud1, ud2);

			// Zapisz aktualny punkt bazowy
			traceFile << iter << ";" << XB.x(0) << ";" << XB.x(1) << ";" << XB.y(0) << "\n";

			X = HJ_trial(ff, XB, s, ud1, ud2);
			if (X.y < XB.y) {
				do {
					XB_temp = XB;
					XB = X;
					X.x = 2 * XB.x - XB_temp.x;
					X.fit_fun(ff, ud1, ud2);
					X = HJ_trial(ff, X, s, ud1, ud2);
					if (X.f_calls > Nmax) {
						Xopt.x = XB.x;
						Xopt.y = XB.y;
						Xopt.flag = 0;
						traceFile.close();
						return Xopt;
					}
				} while (X.y < XB.y);
				X = XB;
			}
			else {
				s = alpha * s;
			}

			iter++;
			if (X.f_calls > Nmax) {
				Xopt.x = XB.x;
				Xopt.y = XB.y;
				Xopt.flag = 0;
				traceFile.close();
				return Xopt;
			}
		} while (s > epsilon);

		Xopt.x = XB.x;
		Xopt.y = XB.y;
		Xopt.flag = 1;

		traceFile.close();
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_with_trace(...):\n" + ex_info);
	}
}


solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		solution::clear_calls();
		int n = get_len(x0);
		matrix d(n, n);
		for (int j = 0; j < n; j++) {
			d(j, j) = 1.0;
		}
		matrix lambda(n, 1), p(n, 1), s(s0);
		solution XB, X_temp;
		XB.x = x0;
		XB.fit_fun(ff);

		while (true) {

			// Szukaj rozwiązania zgodnie z zadaną bazą 
			for (int j = 0; j < n; j++) {
				X_temp.x = XB.x + (s(j) * d[j]);
				X_temp.fit_fun(ff);
				if (X_temp.y(0) < XB.y(0)) {
					XB = X_temp;
					lambda(j) += s(j);
					s(j) *= alpha;
				}
				else {
					p(j) = p(j) + 1;
					s(j) *= -beta;
				}
			}

			// Zweryfikuj poprawność kroków
			bool change = true;
			for (int j = 0; j < n; j++) {
				if (p(j) == 0 || lambda(j) == 0)
				{
					change = false;
					break;
				}

			}

			// Zmiana bazy, jeżeli jest to konieczne
			if (change)
			{
				matrix Q(n, n), v(n, 1);
				for (int i = 0; i < n; ++i)
					for (int j = 0; j <= i; ++j)
						Q(i, j) = lambda(i);

				Q = d * Q;
				v = Q[0] / norm(Q[0]);
				d.set_col(v, 0);
				for (int i = 1; i < n; ++i)
				{
					matrix temp(n, 1);
					for (int j = 0; j < i; ++j)
						temp = temp + (trans(Q[i]) * d[j]) * d[j];
					v = Q[i] - temp;
					d.set_col(v, i);
				}
				s = s0;
				lambda = matrix(n, 1);
				p = matrix(n, 1);
			}

			// Sprawdź dokładność wyniku
			double max_s = abs(s(0));
			for (int i = 1; i < n; ++i) {
				if (max_s < abs(s(i))) {
					max_s = abs(s(i));
				}
			}
			if (max_s < epsilon || solution::f_calls > Nmax) {
				return XB;
			}
		}


		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution Rosen_with_trace(matrix(*ff)(matrix, matrix, matrix),
                          matrix x0, matrix s0,
                          double alpha, double beta,
                          double epsilon, int Nmax,
                          matrix ud1, matrix ud2)
{
    try
    {
        std::ofstream traceFile("trajektoria_Rosen.csv");
        traceFile << "Iteracja;x1;x2;f\n";

        solution Xopt;
        solution::clear_calls();
        int n = get_len(x0);

        matrix d(n, n);
        for (int j = 0; j < n; j++) d(j, j) = 1.0;

        matrix lambda(n, 1), p(n, 1), s(s0);
        solution XB, X_temp;
        XB.x = x0;
        XB.fit_fun(ff, ud1, ud2);

        int iter = 0;

        while (true) {
            // Zapisz punkt bazowy do pliku CSV
            traceFile << iter << ";" << XB.x(0) << ";" << XB.x(1) << ";" << XB.y(0) << "\n";

            // Szukaj rozwiązania zgodnie z zadaną bazą
            for (int j = 0; j < n; j++) {
                X_temp.x = XB.x + (s(j) * d[j]);
                X_temp.fit_fun(ff, ud1, ud2);

                if (X_temp.y(0) < XB.y(0)) {
                    XB = X_temp;
                    lambda(j) += s(j);
                    s(j) *= alpha;
                } else {
                    p(j) = p(j) + 1;
                    s(j) *= -beta;
                }
            }

            bool change = true;
            for (int j = 0; j < n; j++) {
                if (p(j) == 0 || lambda(j) == 0) {
                    change = false;
                    break;
                }
            }

            if (change) {
                matrix Q(n, n), v(n, 1);
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j <= i; ++j)
                        Q(i, j) = lambda(i);

                Q = d * Q;
                v = Q[0] / norm(Q[0]);
                d.set_col(v, 0);
                for (int i = 1; i < n; ++i) {
                    matrix temp(n, 1);
                    for (int j = 0; j < i; ++j)
                        temp = temp + (trans(Q[i]) * d[j]) * d[j];
                    v = Q[i] - temp;
                    d.set_col(v, i);
                }

                s = s0;
                lambda = matrix(n, 1);
                p = matrix(n, 1);
            }

            double max_s = abs(s(0));
            for (int i = 1; i < n; ++i)
                if (max_s < abs(s(i))) max_s = abs(s(i));

            iter++;

            if (max_s < epsilon || solution::f_calls > Nmax) {
                traceFile.close();
                return XB;
            }
        }

        traceFile.close();
        return Xopt;
    }
    catch (string ex_info)
    {
        throw ("solution Rosen_with_trace(...):\n" + ex_info);
    }
}

// Funkcja kary
solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution XB;
		XB.x = x0;

		// Współczynnik kary c w ud2
		matrix c_mat(c);

		XB.fit_fun(ff, ud1, c_mat);

		solution XT;
		XT = XB;

		// Parametry Simplexu
		double s = 0.5;
		double alpha = 1.0, beta = 0.5, gamma = 2.0, delta = 0.5;

		while (true)
		{
			// Aktualizacja macierzy kary dla bieżącej iteracji
			c_mat = matrix(c);

			XT = sym_NM(ff, XB.x, s, alpha, beta, gamma, delta, epsilon, Nmax, ud1, c_mat);

			// Warunek stopu pętli zewnętrznej (zmiana położenia minimum)
			if (norm(XT.x - XB.x) < epsilon)
				break;

			// Aktualizacja współczynnika kary
			// Dla metody zewnętrznej (dc > 1) c rośnie.
			// Dla metody wewnętrznej (dc < 1) c maleje.
			c *= dc;

			// Zabezpieczenie przed pętlą nieskończoną
			if (solution::f_calls > Nmax)
			{
				XT.flag = 0; // Oznaczamy jako nieudane, ale zwracamy wynik
				break;
			}

			XB = XT;
		}

		return XT;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0,
    double s, double alpha, double beta, double gamma, double delta,
    double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    try
    {
        int n = get_len(x0);
        matrix D(n, n);
        for (int i = 0; i < n; ++i) D(i, i) = 1.0;

        // Inicjalizacja simpleksu
        std::vector<solution> simplex(n + 1);
        simplex[0].x = x0;
        // Przekazujemy ud1 (parametry) i ud2 (zawierające c)
        simplex[0].fit_fun(ff, ud1, ud2);

        for (int i = 1; i < n + 1; ++i)
        {
            simplex[i].x = simplex[0].x + s * D[i - 1];
            simplex[i].fit_fun(ff, ud1, ud2);
        }

        while (true)
        {
            // Sortowanie wierzchołków wg wartości funkcji celu
            int i_min = 0, i_max = 0;
            for (int i = 1; i < n + 1; ++i)
            {
                if (simplex[i].y < simplex[i_min].y) i_min = i;
                if (simplex[i].y > simplex[i_max].y) i_max = i;
            }

            // Sprawdzenie kryterium zbieżności
            bool converged = true;
            for (int i = 0; i < n + 1; ++i)
            {
                if (norm(simplex[i].x - simplex[i_min].x) > epsilon)
                {
                    converged = false;
                    break;
                }
            }
            if (converged) return simplex[i_min];

            if (solution::f_calls > Nmax) throw string("Maximum amount of f_calls reached!");

            // Wyznaczanie środka ciężkości (bez najgorszego punktu)
            matrix center(n, 1, 0.0);
            for (int i = 0; i < n + 1; ++i)
            {
                if (i != i_max) center = center + simplex[i].x;
            }
            center = center / (double)n;

            // Odbicie (Reflection)
            solution reflected;
            reflected.x = center + alpha * (center - simplex[i_max].x);
            reflected.fit_fun(ff, ud1, ud2);

            if (reflected.y < simplex[i_min].y)
            {
                // Ekspansja (Expansion)
                solution expansion;
                expansion.x = center + gamma * (reflected.x - center);
                expansion.fit_fun(ff, ud1, ud2);

                if (expansion.y < reflected.y)
                    simplex[i_max] = expansion;
                else
                    simplex[i_max] = reflected;
            }
            else
            {
                if (reflected.y < simplex[i_max].y)
                {
                    simplex[i_max] = reflected;
                }
                else
                {
                    // Zawężenie (Contraction)
                    solution contracted;
                    contracted.x = center + beta * (simplex[i_max].x - center);
                    contracted.fit_fun(ff, ud1, ud2);

                    if (contracted.y < simplex[i_max].y)
                    {
                        simplex[i_max] = contracted;
                    }
                    else
                    {
                        // Redukcja (Shrink)
                        for (int i = 0; i < n + 1; ++i)
                        {
                            if (i != i_min)
                            {
                                simplex[i].x = simplex[i_min].x + delta * (simplex[i].x - simplex[i_min].x);
                                simplex[i].fit_fun(ff, ud1, ud2);
                            }
                        }
                    }
                }
            }
        }
    }
    catch (string ex_info)
    {
        throw ("solution sym_NM(...):\n" + ex_info);
    }
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		Xopt.x = x0;
		double h;
		matrix d;

		// Pierwsze obliczenie wartości funkcji
		Xopt.fit_fun(ff, ud1, ud2);

		while (true)
		{
			// 1. Wyznacz kierunek (najszybszy spadek)
			matrix grad = Xopt.grad(gf, ud1, ud2);
			d = -grad;

			// 2. Wyznacz długość kroku h
			if (h0 != 0) // Krok stały
			{
				h = h0;
			}
			else // Krok zmienny (Metoda Złotego Podziału na kierunku)
			{
				// Szukamy h, które minimalizuje f(x + h*d)
				double a = 0.0, b = 1.0; // Zakładamy początkowy przedział poszukiwań [0, 1] (można użyć ekspansji)
				double alpha = (sqrt(5.0) - 1.0) / 2.0;
				double c = b - alpha * (b - a);
				double d_gold = a + alpha * (b - a);

				// Wartości funkcji w punktach próbnych
				matrix xc = Xopt.x + c * d;
				matrix xd = Xopt.x + d_gold * d;
				double fc = ff(xc, ud1, ud2)(0);
				double fd = ff(xd, ud1, ud2)(0);
				solution::f_calls += 2;

				while ((b - a) > epsilon)
				{
					if (fc < fd)
					{
						b = d_gold;
						d_gold = c;
						fd = fc;
						c = b - alpha * (b - a);
						xc = Xopt.x + c * d;
						fc = ff(xc, ud1, ud2)(0);
						solution::f_calls++;
					}
					else
					{
						a = c;
						c = d_gold;
						fc = fd;
						d_gold = a + alpha * (b - a);
						xd = Xopt.x + d_gold * d;
						fd = ff(xd, ud1, ud2)(0);
						solution::f_calls++;
					}
				}
				h = (a + b) / 2.0;
			}

			// 3. Aktualizacja pozycji
			solution Xnew;
			Xnew.x = Xopt.x + h * d;
			Xnew.fit_fun(ff, ud1, ud2);

			// 4. Warunki stopu
			if (norm(Xnew.x - Xopt.x) < epsilon)
			{
				Xopt = Xnew;
				Xopt.flag = 1;
				break;
			}

			if (solution::f_calls > Nmax || solution::g_calls > Nmax)
			{
				Xopt = Xnew;
				Xopt.flag = 0;
				break;
			}

			Xopt = Xnew;
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		Xopt.x = x0;
		matrix d, grad_old;
		double h;

		// Inicjalizacja
		Xopt.fit_fun(ff, ud1, ud2);
		matrix grad = Xopt.grad(gf, ud1, ud2);
		d = -grad; // Pierwszy kierunek zgodny z SD
		grad_old = grad;

		while (true)
		{
			// 1. Wyznacz długość kroku h (identycznie jak w SD)
			if (h0 != 0)
			{
				h = h0;
			}
			else
			{
				double a = 0.0, b = 1.0;
				double alpha = (sqrt(5.0) - 1.0) / 2.0;
				double c_gold = b - alpha * (b - a);
				double d_gold = a + alpha * (b - a);

				matrix xc = Xopt.x + c_gold * d;
				matrix xd = Xopt.x + d_gold * d;
				double fc = ff(xc, ud1, ud2)(0);
				double fd = ff(xd, ud1, ud2)(0);
				solution::f_calls += 2;

				while ((b - a) > epsilon)
				{
					if (fc < fd)
					{
						b = d_gold;
						d_gold = c_gold;
						fd = fc;
						c_gold = b - alpha * (b - a);
						xc = Xopt.x + c_gold * d;
						fc = ff(xc, ud1, ud2)(0);
						solution::f_calls++;
					}
					else
					{
						a = c_gold;
						c_gold = d_gold;
						fc = fd;
						d_gold = a + alpha * (b - a);
						xd = Xopt.x + d_gold * d;
						fd = ff(xd, ud1, ud2)(0);
						solution::f_calls++;
					}
				}
				h = (a + b) / 2.0;
			}

			// 2. Aktualizacja pozycji
			solution Xnew;
			Xnew.x = Xopt.x + h * d;
			Xnew.fit_fun(ff, ud1, ud2);

			// 3. Warunki stopu
			if (norm(Xnew.x - Xopt.x) < epsilon)
			{
				Xopt = Xnew;
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax || solution::g_calls > Nmax)
			{
				Xopt = Xnew;
				Xopt.flag = 0;
				break;
			}

			// 4. Wyznaczenie nowego kierunku (Fletcher-Reeves)
			grad = Xnew.grad(gf, ud1, ud2);
			double beta = pow(norm(grad), 2) / pow(norm(grad_old), 2);
			d = -grad + beta * d;

			// Aktualizacja zmiennych do następnej iteracji
			grad_old = grad;
			Xopt = Xnew;
		}

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		Xopt.x = x0;
		matrix d;
		double h;

		Xopt.fit_fun(ff, ud1, ud2);

		while (true)
		{
			// 1. Obliczenie gradientu i Hesjanu
			matrix grad = Xopt.grad(gf, ud1, ud2);
			matrix H = Xopt.hess(Hf, ud1, ud2);

			// 2. Wyznaczenie kierunku d = -H^-1 * grad
			// Korzystamy z funkcji inv() z matrix.h
			d = -inv(H) * grad;

			// 3. Wyznaczenie długości kroku h
			if (h0 != 0)
			{
				h = h0;
			}
			else
			{
				double a = 0.0, b = 1.0;
				double alpha = (sqrt(5.0) - 1.0) / 2.0;
				double c_gold = b - alpha * (b - a);
				double d_gold = a + alpha * (b - a);

				matrix xc = Xopt.x + c_gold * d;
				matrix xd = Xopt.x + d_gold * d;
				double fc = ff(xc, ud1, ud2)(0);
				double fd = ff(xd, ud1, ud2)(0);
				solution::f_calls += 2;

				while ((b - a) > epsilon)
				{
					if (fc < fd)
					{
						b = d_gold;
						d_gold = c_gold;
						fd = fc;
						c_gold = b - alpha * (b - a);
						xc = Xopt.x + c_gold * d;
						fc = ff(xc, ud1, ud2)(0);
						solution::f_calls++;
					}
					else
					{
						a = c_gold;
						c_gold = d_gold;
						fc = fd;
						d_gold = a + alpha * (b - a);
						xd = Xopt.x + d_gold * d;
						fd = ff(xd, ud1, ud2)(0);
						solution::f_calls++;
					}
				}
				h = (a + b) / 2.0;
			}

			// 4. Aktualizacja pozycji
			solution Xnew;
			Xnew.x = Xopt.x + h * d;
			Xnew.fit_fun(ff, ud1, ud2);

			// 5. Warunki stopu
			if (norm(Xnew.x - Xopt.x) < epsilon)
			{
				Xopt = Xnew;
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax || solution::g_calls > Nmax || solution::H_calls > Nmax)
			{
				Xopt = Xnew;
				Xopt.flag = 0;
				break;
			}

			Xopt = Xnew;
		}

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		double alpha = (sqrt(5.0) - 1) / 2; // Uzyj sqrt(5.0) dla precyzji i biblioteki cmath
		solution as, bs, cs, ds;
		as.x = a;
		bs.x = b;
		cs.x = bs.x - alpha * (bs.x - as.x);
		ds.x = as.x + alpha * (bs.x - as.x);

		cs.fit_fun(ff, ud1, ud2); // Dodano ud1, ud2
		ds.fit_fun(ff, ud1, ud2); // Dodano ud1, ud2

		while (true) {
			if (cs.y < ds.y) {
				bs.x = ds.x;
				ds.x = cs.x;
				ds.y = cs.y; // Optymalizacja: przepisujemy wynik, nie liczymy od nowa
				cs.x = bs.x - alpha * (bs.x - as.x);
				cs.fit_fun(ff, ud1, ud2);
			}
			else {
				as.x = cs.x;
				cs.x = ds.x;
				cs.y = ds.y; // Optymalizacja
				ds.x = as.x + alpha * (bs.x - as.x);
				ds.fit_fun(ff, ud1, ud2);
			}

			// POPRAWKA: Sprawdzamy globalny licznik wywołań, a nie lokalny > 0
			if (solution::f_calls > Nmax) {
				Xopt.flag = 0;
				break;
			}
			if (bs.x(0) - as.x(0) < epsilon) { // x jest macierzą 1x1, odnosimy się przez indeks (0) lub rzutowanie
				Xopt.flag = 1;
				break;
			}
		}
		Xopt.x = (as.x + bs.x) / 2;
		Xopt.fit_fun(ff, ud1, ud2); // Oblicz wartość w punkcie końcowym

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
