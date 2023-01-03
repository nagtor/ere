#include <iostream>
#include <cmath> // для pow(), cos()
#include <time.h> // для time()
#include <fstream>
#include <ap.h>
#include "linalg.h"
#include "solvers.h"

using namespace std;

const double a = -3.14; //
const double b = 3.14; //

const int K = 3; // кол-во конечных элементов (интервалов)
const int N = 3; // кол-во узлов конечного элемента (включая узлы на концах), многочлены степени (N-1)
const int L = 5*N; // L>N !!! - кол-во рандомных точек на конечном элементе




//для rand
const double eps = (b-a) / 10000; // минимальное расстояние между точками
const double pre = 6; // точность иксов, кол-во знаков после запятой
const int N0 = N+1;

double F(double x){
    return cos(x)*x ;
};


double RandomFloat(double min, double max, int precision){
    double x;
    x = rand() % (int)pow(10, precision); // рандомим число от 0 до 10^pre
    x = min + (x / pow(10, precision)) * (max - min); // масштабируем
    return x;
}


double L_k (int k, double X){ // k-ый/N многочлен Л степ (N-1) на интервале [-1; 1] в точке Х
    int j=0;
    double product = 1;
    double temp = 1;
    
    double *x;
    x = (double*) malloc (N*sizeof(double));
    
    // сетка из N узлов на интервале [-1; 1]
    for (j=0; j<N0; j++){
        x[j] = -1 + j * 2/double(N0-1);
        //cout << x[j] << endl;
    }
    
    for(j=0; j<N0; j++)
        if(k!=j)
            temp = temp * (X-x[j]) / (x[k]-x[j]);
    
    return temp;
};


// x -> x' = (x-a) / ( 1/2*(b-a) ) - 1
// y = y' = L_k (x')
//
double F_k (int k, double a0, double b0, double X){ // k-ая/N базиская  на [a0, b0]
    
    if ( a0 <= X and X <= b0 ){
        double y = L_k( k, (2*(X - a0) / (b0-a0)) - 1);
        return y;
    }
    else return 0;
};

double phi_j (int j, double X){ // j-ая базисная/N*K функция на [a, b], ГДЕ СТЫКИ НЕ СКЛЕЕНЫ, 0 <= j < N*K (!)
    //cout << j%N << " " << N << " " << a+(j/N)*(b-a)/K << " " << a+(j/N+1)*(b-a)/K << " " << X << endl;
    return F_k (j%(N0), a+(j/(N0))*(b-a)/K, a+(j/N0+1)*(b-a)/K, X);
};

double phi_i (int i, double X){ // i-ая/M базисная функция на [a, b], ГДЕ СТЫКИ СКЛЕЕНЫ, 0 <= i < M (!)
    int M = K * (N0-1) + 1;
    int j = i + i/(N0-1);
    double x_i=0;
    //cout << i%(N-1) << endl;
    //cout << endl << "i=" << i << "; j=" << j << "; N=" << N << "; x_i=" << x_i << endl;
    if (i%(N0-1)==0){
        if (i==M-1)
            return phi_j ( j-1, X);
        else{
            x_i = a + ( i/(N0-1) ) * (b-a)/K;
            //cout << "x_" << i << " = " << x_i << endl;
            if( X>x_i )
                return phi_j ( j, X);
            else
                return phi_j ( j-1, X);
        }
    }
    else
        return phi_j ( j, X);
};

double Approx (alglib::real_1d_array c, double X){
    int M =  K * (N0-1) + 1;
    double sum = 0;
    for (int i=0; i<M; i++)
        sum += c[i] * phi_i (i, X);
    return sum;
};

// функция считающая нормы
// С - кол-во точек в р/м сетке, на которой считаем нормы
// c[i] - коэф-ты аппроксимации по базисным функциям конечных элементов
int Norma(int C, alglib::real_1d_array c){
    ofstream fn("norma.txt");
    fn.setf( ios::scientific );
    int M = K * (N0-1) + K + 1;
    int flag;
    double temp;
    int k = 0;
    
    double max1 = 0;
    double max2 = F(a);
    
    double sum1_1 = 0;
    double sum1_2 = 0;
    
    double sum2_1 = 0;
    double sum2_2 = 0;
    
    // сетка из C иксов с шагом h/100 на отрезке [a;b]
    for(int i=1; i<C-1; i++){
        temp = a + i * (b-a)/(C-1);
        //cout << "temp = " << temp << endl;
        
        sum1_1 += abs (Approx(c, temp) - F(temp) );
        //cout << "|...-...| = " << abs (Approx(c, temp) - F(temp) ) << endl;
        //cout << "sum1_1 = " << sum1_1 << endl;
        sum1_2 += abs(F(temp));
        //cout << "F(temp) = " << abs ( F(temp) ) << endl;
        //cout << "sum1_2 = " << sum1_2 << endl << endl;
        
        
        sum2_1 += (Approx(c, temp) - F(temp))*(Approx(c, temp) - F(temp));
        sum2_2 += F(temp) * F(temp);
        
        if (max1 < abs(Approx(c, temp) - F(temp)))
            max1 = abs(Approx(c, temp) - F(temp));
        if (max2 < abs(F(temp)))
            max2 = abs(F(temp));
    }
    
    sum2_1 = sqrt(sum2_1);
    sum2_2 = sqrt(sum2_2);
    
    fn << endl << "Погрешности:        " << endl << endl;
    fn << "                             ||...||_1       ||...||_2       ||...||_inf" << endl;
    fn << "Абсолютная погрешность       " << sum1_1 << "    " << sum2_1 << "    " << max1 << endl;
    fn << "Относительная погрешность    " << sum1_1/sum1_2 << "    " << sum2_1/sum2_2 << "    " << max1/max2 << endl;
    
    fn.close();
    return 0;
};

//________________________________________________________________________________________________

int main(void){
    
    int i,j, k, l;
    //int N0 = N+1;
    double temp;
    int flag=0;
    
    ofstream fout("fout.txt");
    
    //if( L==N or L<N ){
        //cout << "Некорректные L и N :(" << endl;
        //return -1;
    //}
    
    int M = K * (N0-1) + 1; // общее кол-во узлов сетки
    
    double h = (b-a)/(M-1); // шаг для узлов
    
    fout << K << endl;
    fout << M << endl;
    fout << a << " " << b << endl;
    
    //рандомный вектор, на котором будем считать скалярное произведение
    double *X;
    X = (double*) malloc (L*K*sizeof(double));
    // заполняем
    srand(time(NULL));
    for(k=0; k<K; k++){
        //cout << "от " << a + k * (b-a)/K << " до " << a + (k+1) * (b-a)/K << ":" << endl;
        for(l=0; l<L; l++){
            X[l+k*L] = RandomFloat( a + k * (b-a)/K , a + (k+1) * (b-a)/K , pre);
            //cout << X[l+k*L] << endl;
        }
    }
  
    
    alglib::real_2d_array A("[[]]"); // создали пустой двумерный массив alglib
    A.setlength(M,M); // задали размер M*M (выделили память)
    // заполним alglib матрицу 
    for (i=0; i<M; i++){
        for(j=0; j<M; j++){
            temp = 0;
                for (k=0; k<K*L; k++)
                    temp += phi_i (i, X[k]) * phi_i (j, X[k]);
            A[i][j] = temp;
        }
    }
    
    alglib::real_1d_array B("[]"); // создали пустой вектор alglib
    B.setlength(M); // задали длину M (выделили память)

    for(i=0; i<M; i++){
        temp = 0;
        for(k=0; k<K*L; k++)
            temp += F(X[k]) * phi_i (i, X[k]);
        B[i] = temp;
    }
    
    // выведем массив
    /*
    cout << endl;
    for (i=0; i<M; i++){
        for(j=0; j<M; j++){
            cout << _A[i][j] << " ";
        }
        cout << endl;
    }
    */
    
    /*
    //выведем вектор
    for (i=0; i<M; i++)
        cout << _B[i] << endl;
    */
    
    /*  разложение A = LU, A имеет размер M*N
     LU decomposition of a general real matrix with row pivoting

     A is represented as A = P*L*U, where:
     * L is lower unitriangular matrix
     * U is upper triangular matrix
     * P = P0*P1*...*PK, K=min(M,N)-1,
       Pi - permutation matrix for I and Pivots[I]

     INPUT PARAMETERS:
         A       -   array[0..M-1, 0..N-1].
         M       -   number of rows in matrix A.
         N       -   number of columns in matrix A.


     OUTPUT PARAMETERS:
         A       -   matrices L and U in compact form:
                     * L is stored under main diagonal
                     * U is stored on and above main diagonal
         Pivots  -   permutation matrix in compact form.
                     array[0..Min(M-1,N-1)].
     
    void alglib::rmatrixlu(
         real_2d_array& a,
         ae_int_t m,
         ae_int_t n,
         integer_1d_array& pivots,
         const xparams _params = alglib::xdefault);
    */
    
    //cout << _A.rows() << " " << _A.cols() << endl;        методы для получения размеров 2-мерного массива (для вектора - length() )
    alglib::integer_1d_array pivots;
    const alglib::xparams _params = alglib::xdefault;
    alglib::rmatrixlu( A, A.rows(), A.cols(), pivots, _params ); // воспользовались LU-разложением
    
    /*
     Dense solver.

     This  subroutine  solves  a  system  A*x=b,  where A is NxN non-denegerate
     real matrix given by its LU decomposition, x and b are real vectors.  This
     is "fast-without-any-checks" version of the linear LU-based solver. Slower
     but more robust version is RMatrixLUSolve() function.

     Algorithm features:
     * O(N^2) complexity
     * fast algorithm without ANY additional checks, just triangular solver

     INPUT PARAMETERS
         LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
         P       -   array[0..N-1], pivots array, RMatrixLU result
         N       -   size of A
         B       -   array[0..N-1], right part

     OUTPUT PARAMETERS
         Info    -   return code:
                     * -3    matrix is exactly singular (ill conditioned matrices
                             are not recognized).
                             X is filled by zeros in such cases.
                     * -1    N<=0 was passed
                     *  1    task is solved
         B       -   array[N]:
                     * info>0    =>  overwritten by solution
                     * info=-3   =>  filled by zeros
     void alglib::rmatrixlusolvefast(
          real_2d_array lua,
          integer_1d_array p,
          ae_int_t n,
          real_1d_array& b,
          ae_int_t& info,
          const xparams _params = alglib::xdefault);
    */
    alglib::ae_int_t info;
    alglib::rmatrixlusolvefast(A, pivots, A.rows(), B, info, _params); //пытаемся решить систему (решенме Х вкладывается в вектор B)
    
    if (info != 1){
        cout << endl << "Проблемы с решением системы" << endl << endl;
        return -1;
    }
    
    
    //сетка из 1500 иксов на отрезке [a;b]
    // вывод в файл данных: x f(x) approximation(x)
    int C=1000;
    double sum = 0;
    for (i=0; i<C-1; i++ ){
        temp = a + i * (b-a)/(C-1);
        fout << temp << " " << F(temp) << " " << Approx(B, temp) << endl;
    }
    
    
    // сетка из 100*(M-1) иксов с шагом h/100 на отрезке [a;b]
    C = 100*(M-1);
    Norma(C, B);
    
    fout.close();
    
    free(X);
    system ("python3 3.py");
    return 0;
}
