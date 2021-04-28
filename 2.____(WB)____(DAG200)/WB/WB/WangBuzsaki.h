#include<iostream>                                                              //for cout
#include<vector>                                                                //for vector
#include <fstream>                                                              //infile
#include<numeric>                                                               //for accumulate
#include <sstream>																//to_string

//#include"Simplifier.h"

/*
#include <math.h>
#include <fstream>
#include <random>
#include <vector>
#include<numeric> //for accumulate
#include <time.h>
#include <string>
#include<iomanip>
#include <fstream>
#include <string>
#include<vector>                                                                //for vector
#include <algorithm>
#include <omp.h>
#include <sstream>																//to_string
#include <random>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>*/


using namespace std;
typedef std::vector<int> vector_1d_int;                                         //define type of vector integer one dimensional
typedef std::vector<vector_1d_int> vector_2d_int;                               //define type of vector integer two dimensional
typedef vector<double> vectotsaz_one_order;                                     //define type of vector double one dimensional
typedef vector <vectotsaz_one_order> vectotsaz_two_order;                       //define type of vector double two dimensional
using vec1i = vector<int>;                                                      //using of vector integer one dimensional
using vec2i = vector<vector<int>>;                                              //using of vector integer two dimensional
using vec1d = vector<double>;                                                   //using of vector double one dimensional
using vec2d = vector<vector<double>>;                                           //using of vector double two dimensional
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                                class WangBuzsaki                               @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
class WangBuzsaki {                                                              //create and define class 
private:                                                                        //private values
    const int N;                                                                //number of node
    const double dt;                                                            //step lenth
    double tfinal;                                                              //time final
    double phi;                                                                 // phi
    double I_mean;                                                              //mean I external
    double coupling;                                                            //stringh cupling
    string address_matrix;                                                      //name matrix file
    int num_steps;                                                              //number of spike

    vec1d I_app;                                                                //distribution corent
    vec1d IC;                                                                   //v & h & n & s
    vec2i column_connection;                                                    //number of column connection in [row][number of connection in one row]
    double couplingOverN = coupling / (N + 0.0);
    int THRESHOLD = -55;                                                        //threshold spike
    double Cm = 1.0;                                                            //Experimental values of the model
    double E_L = -65.0;                                                         //Experimental values of the model
    double E_Na = 55.0;                                                         //Experimental values of the model
    double E_K = -90.0;                                                         //Experimental values of the model
    double g_L = 0.1;                                                           //Experimental values of the model
    double g_Na = 35.0;                                                         //Experimental values of the model
    double g_K = 9.0;                                                           //Experimental values of the model
    double g_syn = 0.1;                                                         //Experimental values of the model
    double E_syn = -75;                                                         //E sync
    double alpha = 12;                                                          //alpha sync
    double beta = 0.1;                                                          //beta sync
public:                                                                         //public values
    WangBuzsaki
    (int N_, double dt_, double phi_, double I_, double t_final_, double coupling_, string address_matrix_) :
        N(N_), dt(dt_), phi(phi_), I_mean(I_), tfinal(t_final_), coupling(coupling_), address_matrix(address_matrix_)
    {
        num_steps = int(tfinal / dt);
    }



    void set_params( const vec2i&);

    //virtual ~WangBuzsaki() {}  //?

    void runge_kutta4_integrator(vec1d&);
    int Run();
    vec1d dydt(const vec1d&);

    vec1d rand_borger_200();
    void initial(const vec1d&);
    vec2i connection();
};
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                                  rand_borger_200                               @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
vec1d WangBuzsaki::rand_borger_200() {
    ifstream infile;
    infile.open("datarand200.txt");
    if (infile.fail())
    {
        cout << "Could not open file numbers." << "\n";
        //return 1;
    }
    else {
        vec1d phi_vec(N);
        double data;
        infile >> data;
        int i = 0;
        while (!infile.eof()) {
            phi_vec[i] = data;
            i++;
            infile >> data;
        }
        return phi_vec;
    }
}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                              initial conditions                                @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   n_inf_i  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void  n_inf_i(double* n_init, double* v_i, int size_matrix) {
    //double* alpha_n = new double[size_matrix];
    //double* beta_n = new double[size_matrix];
    for (int i = 0; i < size_matrix; i++) {
        double alpha_n = -0.01 * (v_i[i] + 34.0) / (exp(-0.1 * (v_i[i] + 34.0)) - 1.0);
        double beta_n = 0.125 * exp(-(v_i[i] + 44.0) / 80.0);
        n_init[i] = alpha_n / (alpha_n + beta_n);
    }
    //delete[] alpha_n;
    //delete[] beta_n;
}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  h_inf_i  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void  h_inf_i(double* h_init, double* v_i, int size_matrix) {
    //double* alpha_h = new double[size_matrix];
    //double* beta_h = new double[size_matrix];
    for (int i = 0; i < size_matrix; i++) {
        double alpha_h = 0.07 * exp(-(v_i[i] + 58.0) / 20.0);
        double beta_h = 1.0 / (exp(-0.1 * (v_i[i] + 28.0)) + 1.0);
        h_init[i] = alpha_h / (alpha_h + beta_h);
    }
    //delete[] alpha_h;
    //delete[] beta_h;
}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  m_inf_i  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void  m_inf_i(double* m_init, double* v_i, int size_matrix) {
    //double* alpha_m = new double[size_matrix];
    //double* beta_m = new double[size_matrix];
    for (int i = 0; i < size_matrix; i++) {
        double alpha_m = 0.1 * (v_i[i] + 35.0) / (1.0 - exp(-(v_i[i] + 35.0) / 10.0));
        double beta_m = 4.0 * exp(-(v_i[i] + 60.0) / 18.0);
        m_init[i] = alpha_m / (alpha_m + beta_m);
    }
    //delete[] alpha_m;
    //delete[] beta_m;
}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   dv/dt   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void  rk4_vinit(double* v, double* n, double* h, double* m, double I, double v_rev_i, int size_matrix, double dt) {
    for (int i = 0; i < size_matrix; i++) {
        int   g_Na_i = 35;
        int   E_Na_i = 55;
        int   g_K = 9;
        int   E_K = -90;
        double g_l = 0.1;
        int   E_l = -65;
        double  k1, k2, k3, k4;
        k1 = dt * (I - ((g_Na_i * h[i] * pow(m[i], 3) * (v[i] - E_Na_i)) + (g_K * pow(n[i], 4) * (v[i] - E_K)) + (g_l * (v[i] - E_l))));
        k2 = dt * (I - ((g_Na_i * h[i] * pow(m[i], 3) * ((v[i] + k1 / 2) - E_Na_i)) + (g_K * pow(n[i], 4) * ((v[i] + k1 / 2) - E_K)) + (g_l * ((v[i] + k1 / 2) - E_l))));
        k3 = dt * (I - ((g_Na_i * h[i] * pow(m[i], 3) * ((v[i] + k2 / 2) - E_Na_i)) + (g_K * pow(n[i], 4) * ((v[i] + k2 / 2) - E_K)) + (g_l * ((v[i] + k2 / 2) - E_l))));
        k4 = dt * (I - ((g_Na_i * h[i] * pow(m[i], 3) * ((v[i] + k3) - E_Na_i)) + (g_K * pow(n[i], 4) * ((v[i] + k3) - E_K)) + (g_l * ((v[i] + k3) - E_l))));
        v[i] = v[i] + double((1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4));
    }
}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   dn/dt   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void  rk4_n(double* v, double* n, double phi, int size_matrix, double dt) {

    double* tau_n = new double[size_matrix];
    double* n_inf = new double[size_matrix];
    for (int i = 0; i < size_matrix; i++) {
        double alpha_n = -0.01 * (v[i] + 34.0) / (exp(-0.1 * (v[i] + 34.0)) - 1.0);
        double beta_n = 0.125 * exp(-(v[i] + 44.0) / 80.0);
        double sum = alpha_n + beta_n;
        double tau_nn = (1.0 / sum);
        //int phi = 5;
        tau_n[i] = (tau_nn / phi);
        n_inf[i] = alpha_n / sum;

        double  k1, k2, k3, k4;

        k1 = dt * ((n_inf[i] - n[i]) / tau_n[i]);
        k2 = dt * ((n_inf[i] - (n[i] + k1 / 2)) / tau_n[i]);
        k3 = dt * ((n_inf[i] - (n[i] + k2 / 2)) / tau_n[i]);
        k4 = dt * ((n_inf[i] - (n[i] + k3)) / tau_n[i]);
        n[i] = n[i] + double((1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4));
    }

    delete[] tau_n;
    delete[] n_inf;
}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   dh/dt   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void  rk4_h(double* v, double* h, double phi, int size_matrix, double dt) {

    double* tau_h = new double[size_matrix];
    double* h_inf = new double[size_matrix];
    for (int i = 0; i < size_matrix; i++) {
        double alpha_h = 0.07 * exp(-(v[i] + 58.0) / 20.0);
        double beta_h = 1.0 / (exp(-0.1 * (v[i] + 28.0)) + 1.0);
        double sum = alpha_h + beta_h;
        double tau_hh = (1.0 / sum);
        //int phi = 5;
        tau_h[i] = (tau_hh / phi);
        h_inf[i] = alpha_h / sum;

        double  k1, k2, k3, k4;

        k1 = dt * ((h_inf[i] - h[i]) / tau_h[i]);
        k2 = dt * ((h_inf[i] - (h[i] + k1 / 2)) / tau_h[i]);
        k3 = dt * ((h_inf[i] - (h[i] + k2 / 2)) / tau_h[i]);
        k4 = dt * ((h_inf[i] - (h[i] + k3)) / tau_h[i]);
        h[i] = h[i] + double((1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4));
    }
    delete[] tau_h;
    delete[] h_inf;
}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                               initial conditions                                    @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void i1p(int size, double* a, double init) {

    for (int i = 0; i < size; i++) {
        a[i] = init;
    }
}
void i1p(int size, int* a, int init) {

    for (int i = 0; i < size; i++) {
        a[i] = init;
    }
}
void e1p(int size, double* a, double* b) {

    for (int i = 0; i < size; i++) {
        a[i] = b[i];
    }
}
void WangBuzsaki::initial(const vec1d& phi_vec) {                           //call initial in wangbuzsaki class
    double* V_ = new double[N];                                             //create pointer for calculate V
    double* N_ = new double[N];                                             //create pointer for calculate N
    double* H_ = new double[N];                                             //create pointer for calculate H
    double v_rev_i = -75.0;
    //__________________________define variables__________________________  //
    int max_spikes = 3;                                                     //maximum spike for each neurons
    double t_final_init = 2000,                                             //final time run for initial cundition
        t0_init = 0.0;                                                      //time conter and change with dt step 
    vector<int> done(N, 0);                                                 //define conter for finish loop
    int* num_spikes = new int[N];                                           //created 1d pointer for number of spike
    i1p(N, num_spikes, -1);                                                 //initional variable -1
    double* v_init = new double[N];                                         //created 1d pointer for voltage initial
    i1p(N, v_init, -70);                                                    //initional variable -70
    double* n_init = new double[N];                                         //created 1d pointer for N initial
    n_inf_i(n_init, v_init, N);                                             //initional variable  n_inf_i
    double* h_init = new double[N];                                         //created 1d pointer for H initial
    h_inf_i(h_init, v_init, N);                                             //initional variable  h_inf_i
    double* m_init = new double[N];                                         //created 1d pointer for M initial
    m_inf_i(m_init, v_init, N);                                             //initional variable  m_inf_i
    double** t_spikes = new double* [N];                                    //created 2d pointer
    for (int i = 0; i < N; i++) {                                           //
        t_spikes[i] = new double[3];                                        //
    }                                                                       //
    double* v_init_old = new double[N];                                     //saver last initial variable V
    double* n_init_old = new double[N];                                     //saver last initial variable N
    double* h_init_old = new double[N];                                     //saver last initial variable H
    //___________________calculate initial condition______________________  //
    while (accumulate(done.begin(), done.end(), 0) < N && t0_init < t_final_init) {
        //_________________________define variables___________________________  //
        double t0_init_old = t0_init;                                       //saver last initial variable time
        e1p(N, v_init_old, v_init);
        e1p(N, h_init_old, h_init);
        e1p(N, n_init_old, n_init);
        rk4_vinit(v_init, n_init, h_init, m_init, I_mean, v_rev_i, N, dt);
        rk4_n(v_init, n_init, phi, N, dt);
        rk4_h(v_init, h_init, phi, N, dt);
        m_inf_i(m_init, v_init, N);
        t0_init = t0_init + dt;
        //__________________________________________________________________find spikes_______________________________________

        vectotsaz_one_order ind(0, 0);
        for (int i = 0; i < N; i++) {
            if (v_init_old[i] >= -20 && v_init[i] < -20) {
                ind.push_back(i);
            }
        }
        int L = ind.size();

        for (int i = 0; i < L; i++) {
            int k = ind[i];
            num_spikes[k] = num_spikes[k] + 1;
            t_spikes[k][num_spikes[k]] = (t0_init_old * (-20 - v_init[k]) + t0_init * (v_init_old[k] - (-20))) / (v_init_old[k] - v_init[k]);
        }

        double* thr = new double[N];
        for (int i = 0; i < N; i++) {
            thr[i] = t_spikes[i][max_spikes - 1] + phi_vec[i] * (t_spikes[i][max_spikes - 1] - t_spikes[i][max_spikes - 2]);
        }

        vectotsaz_one_order ind2(0, 0);
        for (int i = 0; i < N; i++) {
            if (num_spikes[i] == (max_spikes - 1) && t0_init > thr[i] & t0_init_old <= thr[i]) {
                ind2.push_back(i);
            }
        }
        int L2 = ind2.size();

        //__________________________________________________________________add values in wb_init ____________________________
        for (int i = 0; i < L2; i++) {
            int k = ind2[i];
            V_[k] = (v_init_old[k] * (t0_init - thr[k]) + v_init[k] * (thr[k] - t0_init_old)) / dt;
            H_[k] = (h_init_old[k] * (t0_init - thr[k]) + h_init[k] * (thr[k] - t0_init_old)) / dt;
            N_[k] = (n_init_old[k] * (t0_init - thr[k]) + n_init[k] * (thr[k] - t0_init_old)) / dt;
            done[ind2[i]] = 1;
        }
        delete[] thr;
    }

    IC.resize(4 * N);

    for (int i = 0; i < N; i++) {
        IC[i] = V_[i];
        //cout << IC[i] << endl;
        IC[N + i] = H_[i];
        IC[2 * N + i] = N_[i];
    }
    done.clear();
    delete[] V_;
    delete[] H_;
    delete[] N_;
    delete[] num_spikes;
    delete[] v_init;
    delete[] n_init;
    delete[] h_init;
    delete[] m_init;
    delete[] t_spikes;
    delete[] v_init_old;
    delete[] n_init_old;
    delete[] h_init_old;
}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                          Read matrix connection                                @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void matrix(int size_matrix, const std::string address_matrix, double** a) {

    vector_2d_int matrix_2d;                                                    //create vector integer two dimensional(matrix_2d)
    const string files = address_matrix + ".txt";

    std::ifstream ifile(files);                      //read address of file .txt
    if (ifile.is_open()) {                                                      //if file was available 
        vector_1d_int matrix_1d;                                                //create vector integer one dimensional(matrix_1d)
        int num;                                                                //create integer number(num)
        while (ifile >> num) {                                                  //Set the read number of the file to the defined integer(num)
            matrix_1d.push_back(num);                                           //set defined integer(num) into the One after the last cell vector integer one dimensional(matrix_1d)
            if (matrix_1d.size() == size_matrix) {                              //if size of vector is full
                matrix_2d.push_back(matrix_1d);                                 //set vector integer one dimensional(matrix_1d) into the One after the last cell vector integer tow dimensional(matrix_2d)
                matrix_1d.clear();                                              //clean all cels of vector integer one dimensional(matrix_1d) and delete it
            }
        }
    }
    else {                                                                      //else if file wasn't available
        std::cout << "There was an error opening the input file!\n";            //print "..."
        exit(1);                                                                //means program(process) terminate normally unsuccessfully..
    }


    for (int i = 0; i < size_matrix; i++) {
        for (int j = 0; j < size_matrix; j++) {
            a[i][j] = matrix_2d[i][j];
        }
    }
    matrix_2d.clear();
}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                            column_connection                                   @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
vec2i WangBuzsaki::connection()
{
    //read matrix for example 200 dag[satr][soton] read matrix

    vec2i g_ii(N, vec1i(N));

    double** matrix_file = new double* [N];//created 2d pointer
    for (int i = 0; i < N; i++) {
        matrix_file[i] = new double[N];
    }
    matrix(N, address_matrix, matrix_file);//read file with address and size and puting up in 2d pointer//

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            g_ii[i][j] = matrix_file[i][j];
        }
    }
    delete[] matrix_file;
    //[satr][chandomin etesal dar satr]=soton
    int row = g_ii.size();
    int col = g_ii[0].size();
    vec2i column_connection;
    column_connection.resize(row);

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            if (g_ii[i][j] > 1.e-8) {
                column_connection[i].push_back(j);
            }
        }
    }

    IC.resize(4 * N);
    I_app.resize(N);
    for (int i = 0; i < N; i++)
    {
        IC[3 * N + i] = 0;
        I_app[i] = I_mean;
    }

    return column_connection;

}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                                     RUN                                        @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void WangBuzsaki::set_params(
    const vec2i& column_connection)

{

    this->column_connection = column_connection;

}

int WangBuzsaki::Run() {

    set_params(connection());//[satr][chandomin etesal dar satr]=soton

    vec1d y = IC;

    /*const string subName = address_matrix + "-" +//create name for save .txt_____dag2__K___phi___I____beta
        to_string(coupling) + "-" +
        to_string(phi) + "-" +
        to_string(I_mean) + "-" +
        to_string(beta) + ".txt";

    std::string synapseFileName = "syn-" + subName;//save synaps

    //cout << synapseFileName << endl;
    ofstream synapse(synapseFileName, ios::out | ios::trunc);*/

    double t = 0.0;
    double* V_old = new double[N];
    ofstream temp("temp.txt", ios::out | ios::trunc);
    for (int it = 1; it < num_steps; ++it)
    {
        for (int i = 0; i < N; i++) {
            V_old[i] = y[i];
        }
        //cout << it << '\t' << "t=" << t << '\t' << "y1=" << y[0] << endl;

        runge_kutta4_integrator(y);

        t = it * dt;

        for (int i = 0; i < N; i++) {
            if (V_old[i] < THRESHOLD & y[i] >= THRESHOLD  & t > tfinal - 500) {
            //& t0>t_final-500) {
                temp << t << '\t' << i << endl;
            }
        }

        /*synapse << t;
        for (int i = 0; i < N; i++)
            synapse << '\t' << y[i];
        synapse << endl;*/

    }

    return 0;
}


void WangBuzsaki::runge_kutta4_integrator(vec1d& y)
{

    int n = y.size();
    vec1d k1(n), k2(n), k3(n), k4(n);
    vec1d f(n);

    k1 = dydt(y);
    for (int i = 0; i < n; i++)
        f[i] = y[i] + 0.5 * dt * k1[i];

    k2 = dydt(f);
    for (int i = 0; i < n; i++)
        f[i] = y[i] + 0.5 * dt * k2[i];

    k3 = dydt(f);
    for (int i = 0; i < n; i++)
        f[i] = y[i] + dt * k3[i];

    k4 = dydt(f);
    for (int i = 0; i < n; i++)
        y[i] += (k1[i] + 2.0 * (k2[i] + k3[i]) + k4[i]) * dt / 6.0;
}

vec1d WangBuzsaki::dydt(const vec1d& x0)
{
    const int n2 = 2 * N;
    const int n3 = 3 * N;

    vec1d f(4 * N);

    for (int i = 0; i < N; i++)
    {
        double alpha_m = -0.1 * (x0[i] + 35.0) / (exp(-0.1 * (x0[i] + 35.0)) - 1.0);
        double alpha_h = 0.07 * exp(-(x0[i] + 58.0) / 20.0);
        double alpha_n = -0.01 * (x0[i] + 34.0) / (exp(-0.1 * (x0[i] + 34.0)) - 1.0);

        double beta_m = 4.0 * exp(-(x0[i] + 60.0) / 18.0);
        double beta_h = 1.0 / (exp(-0.1 * (x0[i] + 28.0)) + 1.0);
        double beta_n = 0.125 * exp(-(x0[i] + 44.0) / 80.0);

        double m = alpha_m / (alpha_m + beta_m);
        double F = 1.0 / (1.0 + exp(-0.5 * x0[i]));

        double I_Na = g_Na * m * m * m * x0[i + N] * (x0[i] - E_Na);
        double I_L = g_L * (x0[i] - E_L);
        double I_K = g_K * x0[i + n2] * x0[i + n2] * x0[i + n2] * x0[i + n2] * (x0[i] - E_K);

        double I_syn = 0.0;
        int counter = 0;

       // cout << column_connection_size[i].size() << "tdj" << endl;

        for (int j = 0; j < column_connection[i].size(); j++)
        {

            int k = column_connection[i][j];
            I_syn += (g_syn * x0[k + n3] * (x0[k] - E_syn));
        }

        I_syn = couplingOverN * I_syn;

        f[i] = -I_Na - I_K - I_L - I_syn + I_app[i];                             // dv/dt
        f[i + N] = phi * (alpha_h * (1 - x0[i + N]) - beta_h * x0[i + N]);       // dh/dt
        f[i + 2 * N] = phi * (alpha_n * (1 - x0[i + n2]) - beta_n * x0[i + n2]); // dn/dt
        f[i + 3 * N] = alpha * F * (1 - x0[i + n3]) - beta * x0[i + n3];         // ds/dt
    }

    return f;
}


