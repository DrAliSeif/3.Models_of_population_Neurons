#include <iostream>
#include"Simplifier.h"
#include"Header.h"
#include"initial.h"


#define size_matrix 100                                                         //size of matrix clom and row
#define dt 0.01

using namespace std;


int main() {
    clock_t t = clock();//time start record

    double
        sigma = 0,//standard deviation
        v_rev_i = -75.0,
        g_hat_ii = 0.5,
        p_ii = 1,
        g_hat_gap = 0,
        p_gap = 1,
        t_final = 500,
        dt2 = dt / 2,
        I = 1.5,
        tau_r_i = 0.5,//
        tau_peak_i = 0.5,//
        tau_d_i = 9;//


    double* V = new double[size_matrix];
    double* N = new double[size_matrix];
    double* H = new double[size_matrix];
    double* M = new double[size_matrix];
    double* Q = new double[size_matrix] {0};
    double* S = new double[size_matrix] {0};

    initial(V,H,N,M,I);

    //p1p(size_matrix, V);

//__________________________________________________________________function of tau dqi_______________________________
    double  tau_dq_i = tau_dq_function(tau_d_i, tau_r_i, tau_peak_i); //0.1163;
 //__________________________________________________________________Process network parameter_________________________
        //gii_phi_vec(g_ii,  phi_vec);//for 200

    double** g_ii = new double* [size_matrix];//created 2d pointer
    for (int i = 0; i < size_matrix; i++) {
        g_ii[i] = new double[size_matrix];
    }
    i2p(size_matrix, g_ii , ((g_hat_ii * p_ii) / (size_matrix * p_ii)));//initional variable g_ii*/



    ofstream temp("temp.txt", ios::out | ios::trunc);

    //__________________________________________________________________Solve the system using the midpoint method________
    double* V_old = new double[size_matrix];

    for (register double t0 = 0; t0 <= t_final - dt; t0 = t0 + dt) {


            rk4_run_step(g_ii, V, V_old, H, N, M, S, Q, I, v_rev_i, tau_dq_i, tau_r_i, tau_d_i);

            //rk4_v(g_ii, V, V_old, H, N, M, S, I, v_rev_i);
            //rk4_n(V, N);
           // rk4_h(V, H);
            //m_inf_i(M, V);
            //rk4_q(V, Q, tau_dq_i);
            //rk4_s(V, Q, S, tau_r_i, tau_d_i);






        for (int i = 0; i < size_matrix; i++) {
            if (V_old[i] < -20 & V[i] >= -20) {
                temp << t0 << '\t' << i << endl;
            }
        }
    }

    temp.close();//*/
    cout << "\nTime taken by program is :\t" << ((double)clock() - t) / CLOCKS_PER_SEC << " sec\n";
    cout << "finish\n";
    return 0;
}
