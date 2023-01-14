#include <cmath>
#include <cstdlib>
#include <iostream>
#include <time.h>

using namespace std;

#define CALKOWITE 8       // ilość bitów oznaczających liczby całkowite
#define ULAMEK 8          // ilość bitów oznaczających liczby po przecinku
#define ILOSC_PROBEK 1200 // ilość próbek na wykresie

// wzorcowe wartości dla wykresów
#define WZORCOWE_K 128
#define WZORCOWE_T 128
#define WZORCOWE_XI 0.5

typedef const long double &cldouble;

long double skokowa(cldouble time, cldouble K, cldouble T, cldouble XI) {
    long double e = exp(1);

    long double result = K * (1 - pow(e, (-XI * time / T)) / sqrt(1 - XI * XI) * sin(sqrt(1 - XI * XI) / T * time + atan(sqrt(1 - XI * XI) / XI)));
    return result;
}

long double impulsowa(cldouble time, cldouble K, cldouble T, cldouble XI) {
    long double e = exp(1);

    long double result = (K / (T * sqrt(1 - XI * XI)) * pow(e, (-XI * time / T)) * sin(sqrt(1 - XI * XI) / T * time));
    return result;
}
/*
    A - tablica binarna
    ilosc_calkowitych - ilosc liczb calkowitych w tablicy binarnej
    ilosc_w_ulamku - ilosc liczb po przecinku w tablicy binarnej
*/
long double bin_to_dec(const int A[], const int &ilosc_calkowitych, const int &ilosc_w_ulamku) {
    long double sum = 0;
    int dlugosc = ilosc_calkowitych + ilosc_w_ulamku - 1;
    if (ilosc_w_ulamku > 0) {
        for (int i = 0; i < ilosc_w_ulamku; i++) {
            if (A[dlugosc - i])
                sum += pow(2., -(ilosc_w_ulamku - i));
        }
    }

    if (ilosc_calkowitych > 0) {
        int p = 1;
        for (int i = 0; i < ilosc_calkowitych; ++i) {
            if (A[dlugosc - ilosc_w_ulamku - i])
                sum += p;
            p *= 2;
        }
    }
    return sum;
}
/*
    Funkcja sprawdza czy podana liczba znajduje się w przedziale poszukiwań
    A- tablica binarna
    calkowite - ilosc liczb calkowitych w liczbie
    ulamek - ilosc liczb w ulamku
*/
bool czy_w_przedziale(const int A[], const int &calkowite, const int &ulamek, cldouble start, cldouble end) {
    long double tmp = bin_to_dec(A, calkowite, ulamek);
    if (tmp < end && tmp > start)
        return true;
    return false;
}

/*
    Funkcja zeruje tablice
*/
void zeruj(int A[], const int &calkowite, const int &ulamek) {
    for (int i = 0; i < calkowite + ulamek; ++i)
        A[i] = 0;
}

/*
    Funkcja pomocnicza wyznacza początkową wartość liczby binarnej,
    podczas inicjalizacji
*/
void wyznacz_poczatek(int A[], const int &calkowite, const int &ulamek, cldouble start, cldouble end) {
    long double minimum = 0;
    long double tmp = 0;
    int *B = new int[calkowite + ulamek];
    zeruj(A, calkowite, ulamek);
    zeruj(B, calkowite, ulamek);

    for (int i = 0; i < calkowite + ulamek; ++i) {

        B[i] = 1;
        tmp = bin_to_dec(B, calkowite, ulamek);
        if ((start - minimum) >= tmp) {
            minimum += tmp;
            A[i] = 1;
        }
        B[i] = 0;
    }

    A[calkowite + ulamek - 1] = 1;
}

/*
    Funkcja losuje populacje początkową
*/
void inicjalizacja(int A[], const int &calkowite, const int &ulamek, cldouble start, cldouble end) {
    long double sum = 0;
    wyznacz_poczatek(A, calkowite, ulamek, start, end);
    for (int i = 0; i < calkowite + ulamek; i++) {
        A[i] = rand() % 2;
        if (!czy_w_przedziale(A, calkowite, ulamek, start, end)) {
            A[i] = 1 - A[i];
        }
    }
}

/*
    Funkcja obicza wartość wskaźnika jakości
*/
long double przystosowanie(cldouble K, cldouble T, cldouble XI) {
    long double sum_skok = 0;
    long double sum_impuls = 0;
    long double tmp;
    for (long double i = 0; i < ILOSC_PROBEK; i += 1) {

        tmp = skokowa(i, WZORCOWE_K, WZORCOWE_T, WZORCOWE_XI) - skokowa((long double)i, K, T, XI);
        tmp *= tmp;
        sum_skok += tmp;
        tmp = impulsowa(i, WZORCOWE_K, WZORCOWE_T, WZORCOWE_XI) - impulsowa((long double)i, K, T, XI);
        tmp *= tmp;
        sum_impuls += tmp;
    }

    long double f = sum_skok / (long double)ILOSC_PROBEK + sum_impuls / (long double)ILOSC_PROBEK;

    return f;
}

/*
    Funkcja pomocnicza do sortowania
*/
int compare(const void *a, const void *b) {
    const long double *x = (long double *)a;
    const long double *y = (long double *)b;

    if (*x < *y)
        return 1;
    else if (*x > *y)
        return -1;

    return 0;
}

/*
    Funkcja losuje osobników z koła ruletki, przy czym najlepszy osobnik o najniższej wartości
    wskaźnika jakości największą część koła
    przystosowanie - tablica z wartościa przystosowania osobników
    selekcja - tablica z wynikiem losowanie
*/
void selekcja_rankingowa(long double przystosowanie[], int selekcja[], const int &wielkosc_populacji) {

    // sortowanie tablicy przystosowań malejąco,
    // po to aby najlepszy osobnik zajmował największą część koła
    qsort(przystosowanie, wielkosc_populacji, sizeof(long double), compare);

    // tablica sum prefixowych
    long double *sum = new long double[wielkosc_populacji];
    sum[0] = 0;

    // sumowanie rankingów
    for (int i = 1; i < wielkosc_populacji; ++i) {
        sum[i] = i + sum[i - 1];
    }
    // liczenie dystrybuanty
    for (int i = 0; i < wielkosc_populacji; ++i)
        sum[i] /= sum[wielkosc_populacji - 1];

    for (int i = 0; i < wielkosc_populacji; ++i) {
        // r - losowanie z przedziału [0,1]
        long double r = rand() / (long double)(RAND_MAX);
        int j = 0;
        for (j = 0; sum[j] < r; ++j)
            ;
        selekcja[i] = j;
    }
    delete[] sum;
}

/*
    Funkcja krzyżuje 2 fenotypy
    KT1 - fenotyp 1
    KT2 - fenotyp 2
    calkowite - ilosc liczb calkowitych w tablicy fenotypu
    ulamek - ilosc liczb po przecinku w tablicy fenotypu
    start - lewy przedział zakresu poszukiwań
    end - prawy przedział zakresu poszukiwań
*/
void krzyzowanie(int KT1[], int KT2[], const int &calkowite, const int &ulamek, cldouble start, cldouble end) {
    int n = calkowite + ulamek;
    int pozycja = rand() % (n - 1) + 1;
    for (int i = pozycja; i < n; i++) {
        int p = KT1[i];
        KT1[i] = KT2[i];
        KT2[i] = p;
    }
    int a = czy_w_przedziale(KT1, calkowite, ulamek, start, end);
    int b = czy_w_przedziale(KT2, calkowite, ulamek, start, end);

    // jeżeli obie liczby znajdują się w przedziale losowania
    if (a && b)
        // wszystko OK
        return;
    else
        // przywróć fenotyp do poprzedniej wartośći
        for (int i = pozycja; i < n; i++) {
            int p = KT1[i];
            KT1[i] = KT2[i];
            KT2[i] = p;
        }
}
/*
    Funkcja mutuje jeden gen
    KT - fenotyp
    calkowite - ilosc liczb calkowitych w tablicy fenotypu
    ulamek - ilosc liczb po przecinku w tablicy fenotypu
    start - lewy przedział zakresu poszukiwań
    end - prawy przedział zakresu poszukiwań
*/
void mutacja(int KT[], int calkowite, int ulamek, cldouble start, cldouble end) {
    int pozycja = rand() % (calkowite + ulamek);
    KT[pozycja] = 1 - KT[pozycja];

    // jeżeli liczba zmutowana nie znajduje się w przedziale poszukiwań
    // przywróć jej poprzednią wartość
    if (!czy_w_przedziale(KT, calkowite, ulamek, start, end)) {
        KT[pozycja] = KT[pozycja] = 1 - KT[pozycja];
    }
}

/*
    Funkcja pomocnicza
*/
void kopiujDo(int dokad[], int skad[], int n) {
    for (int i = 0; i < n; i++) {
        dokad[i] = skad[i];
    }
}

int main() {
    srand(time(0));

    // inicjalizacja podstawowych parametrów zadania
    int wielkosc_populacji;
    int liczba_iteracji;
    long double prawdopodobienstwo_krzyzowania;
    long double prawdopodobienstwo_mutacji;
    long double start_K, end_K, start_T, end_T;

    long double tmp_K;
    long double tmp_T;
    long double tmp_XI;

    cout << "Podaj liczbe iteracji: ";
    cin >> liczba_iteracji;
    cout << "\nPodaj wielkosc populacji: ";
    cin >> wielkosc_populacji;
    cout << "\nPodaj zakres poszukiwań K <start,end>: ";
    cin >> start_K >> end_K;
    cout << "\nPodaj zakres poszukiwań T <start,end>: ";
    cin >> start_T >> end_T;
    cout << "\nPodaj przedzial prawdopodobienstwa krzyzowania: ";
    cin >> prawdopodobienstwo_krzyzowania;
    cout << "\nPodaj przedzial prawdopodobienstwa mutacji: ";
    cin >> prawdopodobienstwo_mutacji;

    // chromosom jest podzielony na 3 zmienne
    int K[wielkosc_populacji][CALKOWITE + ULAMEK];
    int T[wielkosc_populacji][CALKOWITE + ULAMEK];
    int XI[wielkosc_populacji][ULAMEK];

    for (int i = 0; i < wielkosc_populacji; ++i) {
        inicjalizacja(K[i], CALKOWITE, ULAMEK, start_K, end_K);
        inicjalizacja(T[i], CALKOWITE, ULAMEK, start_T, end_T);
        inicjalizacja(XI[i], 0, ULAMEK, 0, 1);
    }

    for (int i = 0; i < liczba_iteracji; i++) {
        // wartosci funkcji przystosowania osobnikow
        long double wartosc_przystosowania[wielkosc_populacji];

        // w tej zmiennej przetrzymywane są wyniki selekcji
        int selekcja[wielkosc_populacji];

        // Obliczanie
        for (int j = 0; j < wielkosc_populacji; j++) {
            tmp_K = bin_to_dec(K[j], CALKOWITE, ULAMEK);
            tmp_T = bin_to_dec(T[j], CALKOWITE, ULAMEK);
            tmp_XI = bin_to_dec(XI[j], 0, ULAMEK);

            wartosc_przystosowania[j] = przystosowanie(tmp_K, tmp_T, tmp_XI);
        }
        // Selekcja
        selekcja_rankingowa(wartosc_przystosowania, selekcja, wielkosc_populacji);

        // Kopiowanie wybranych osobnikow do tablicy pomocniczej
        int temp_K[wielkosc_populacji][CALKOWITE + ULAMEK];
        int temp_T[wielkosc_populacji][CALKOWITE + ULAMEK];
        int temp_XI[wielkosc_populacji][ULAMEK];

        // zapisywanie wylosowanych osobników do tablicy pomocniczej
        for (int j = 0; j < wielkosc_populacji; j++) {
            kopiujDo(temp_K[j], K[selekcja[j]], CALKOWITE + ULAMEK);
            kopiujDo(temp_T[j], T[selekcja[j]], CALKOWITE + ULAMEK);
            kopiujDo(temp_XI[j], T[selekcja[j]], ULAMEK);
        }

        // Kopiowanie zawartosci tablicy pomocniczej do podstawowej
        for (int j = 0; j < wielkosc_populacji; j++) {
            kopiujDo(K[j], temp_K[j], CALKOWITE + ULAMEK);
            kopiujDo(T[j], temp_T[j], CALKOWITE + ULAMEK);
            kopiujDo(XI[j], temp_XI[j], ULAMEK);
        }
        for (int j = 0; j < wielkosc_populacji; j++) {
            // krzyzowanie
            if (prawdopodobienstwo_krzyzowania < rand() / (long double)(RAND_MAX)) {
                int drugi = rand() % wielkosc_populacji;
                krzyzowanie(K[j], K[drugi], CALKOWITE, ULAMEK, start_T, end_T);
            }
            if (prawdopodobienstwo_krzyzowania < rand() / (long double)(RAND_MAX)) {
                int drugi = rand() % wielkosc_populacji;
                krzyzowanie(T[j], T[drugi], CALKOWITE, ULAMEK, start_T, end_T);
            }
            if (prawdopodobienstwo_krzyzowania < rand() / (long double)(RAND_MAX)) {
                int drugi = rand() % wielkosc_populacji;
                krzyzowanie(XI[j], XI[drugi], 0, ULAMEK, 0, 1);
            }
            // mutacja
            if (prawdopodobienstwo_mutacji < rand() / (long double)(RAND_MAX))
                mutacja(K[j], CALKOWITE, ULAMEK, start_K, end_K);
            if (prawdopodobienstwo_mutacji < rand() / (long double)(RAND_MAX))
                mutacja(T[j], CALKOWITE, ULAMEK, start_T, end_T);
            if (prawdopodobienstwo_mutacji < rand() / (long double)(RAND_MAX))
                mutacja(XI[j], 0, ULAMEK, 0, 1);
        }
    }

    long double najlepsze_przystosowanie;
    int najlepszy_osobnik;
    long double tmp_przystosownaie;
    cout << "\nPopulacja koncowa " << endl;
    for (int i = 0; i < wielkosc_populacji; i++) {
        tmp_K = bin_to_dec(K[i], CALKOWITE, ULAMEK);
        tmp_T = bin_to_dec(T[i], CALKOWITE, ULAMEK);
        tmp_XI = bin_to_dec(XI[i], 0, ULAMEK);

        // wyświetlanie ostatniej populacji
        cout << "K - ";
        for (int j = 0; j < CALKOWITE + ULAMEK; j++) {
            cout << K[i][j];
            if (j == 7)
                cout << ".";
        }
        cout << " = " << tmp_K;
        cout << "  T - ";
        for (int j = 0; j < CALKOWITE + ULAMEK; j++) {
            cout << T[i][j];
            if (j == 7)
                cout << ".";
        }
        cout << " = " << tmp_T;

        cout << "  XI - .";
        for (int j = 0; j < ULAMEK; j++)
            cout << XI[i][j];
        cout << " = " << tmp_XI;
        cout << endl;

        // wyszukiwanie najlepszego osobnika
        if (i == 0) {
            najlepsze_przystosowanie = przystosowanie(tmp_K, tmp_T, tmp_XI);
            najlepszy_osobnik = i;
        }
        tmp_przystosownaie = przystosowanie(tmp_K, tmp_T, tmp_XI);
        if (tmp_przystosownaie < najlepsze_przystosowanie) {
            najlepsze_przystosowanie = tmp_przystosownaie;
            najlepszy_osobnik = i;
        }
    }

    cout << "\n\nNajlepszy osobnik: " << najlepsze_przystosowanie
         << " " << bin_to_dec(K[najlepszy_osobnik], CALKOWITE, ULAMEK)
         << " " << bin_to_dec(T[najlepszy_osobnik], CALKOWITE, ULAMEK)
         << " " << bin_to_dec(XI[najlepszy_osobnik], 0, ULAMEK) << endl;

    // wyswietlanie parametrów
    long double sum_skok = 0;
    long double sum_impuls = 0;
    long double tmp_skok, tmp_impuls;
    long double tmp_g = 0;
    long double wzorzec_skok = 0, wzorzec_impuls = 0;
    long double sum_skok_wzorcowych = 0;
    long double sum_impuls_wzorcowych = 0;

    tmp_K = bin_to_dec(K[najlepszy_osobnik], CALKOWITE, ULAMEK);
    tmp_T = bin_to_dec(T[najlepszy_osobnik], CALKOWITE, ULAMEK);
    tmp_XI = bin_to_dec(XI[najlepszy_osobnik], 0, ULAMEK);

    for (long double i = 0; i < 3000; i += 1) {
        tmp_skok = skokowa(i, tmp_K, tmp_T, tmp_XI);
        tmp_impuls = impulsowa(i, tmp_K, tmp_T, tmp_XI);
        wzorzec_skok = skokowa(i, WZORCOWE_K, WZORCOWE_T, WZORCOWE_XI);
        wzorzec_impuls = impulsowa(i, WZORCOWE_K, WZORCOWE_T, WZORCOWE_XI);
        cout << tmp_skok << " " << wzorzec_skok << " " << tmp_impuls << " " << wzorzec_impuls << endl;
    }

    return 0;
}
