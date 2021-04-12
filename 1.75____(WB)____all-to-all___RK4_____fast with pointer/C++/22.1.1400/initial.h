#include<iostream>                                                              //for cout
#include <iostream>
#include <math.h>
#include <fstream>
#include <random>
#include <vector>
#include<numeric> //for accumulate
#include <time.h>
#include <string>
#include<iomanip>
#include"Simplifier.h"
#include"Header.h"



#define size_matrix 100                                                         //size of matrix clom and row
#define dt 0.01

using namespace std;
typedef vector<double> vectotsaz_one_order;
typedef vector <vectotsaz_one_order> vectotsaz_two_order;

void initial(double* V, double* H, double* N, double* M, double I) {

    double
        v_rev_i = -75.0;




//__________________________________________________________________random numbers of borger pic______________________
    double* phi_vec = new double[size_matrix];//created 1d pointer
    rand_borger(phi_vec);
    //__________________________________________________________________define variables__________________________________
    int
        max_spikes = 3;
    double
        t_final_init = 2000,
        t0_init = 0.0;

    vector<int> done(size_matrix, 0);
    int* num_spikes = new int[size_matrix];//created 1d pointer
    i1p(size_matrix, num_spikes, -1);//initional variable -1
    double* v_init = new double[size_matrix];
    i1p(size_matrix, v_init, -70);//initional variable -70
    double* n_init = new double[size_matrix];
    n_inf_i(n_init, v_init);//initional variable  n_inf_i
    double* h_init = new double[size_matrix];
    h_inf_i(h_init, v_init);//initional variable  h_inf_i
    double* m_init = new double[size_matrix];
    m_inf_i(m_init, v_init);//initional variable  m_inf_i
    //p1p(size_matrix, m_init);



    double** t_spikes = new double* [size_matrix];//created 2d pointer
    for (int i = 0; i < size_matrix; i++) {
        t_spikes[i] = new double[3];
    }

    int x = 0;

    double* v_init_old = new double[size_matrix];
    double* n_init_old = new double[size_matrix];
    double* h_init_old = new double[size_matrix];


    //__________________________________________________________________calculate initial condition_______________________
    while (accumulate(done.begin(), done.end(), 0) < size_matrix && t0_init < t_final_init) {
        x = x + 1;
        //__________________________________________________________________define variables__________________________________
        double
            t0_init_old = t0_init;

        e1p(size_matrix, v_init_old, v_init);
        e1p(size_matrix, h_init_old, h_init);
        e1p(size_matrix, n_init_old, n_init);


        rk4_vinit(v_init, n_init, h_init, m_init, I, v_rev_i);
        rk4_n(v_init, n_init);
        rk4_h(v_init, h_init);

        m_inf_i(m_init, v_init);

        t0_init = t0_init + dt;
        //__________________________________________________________________find spikes_______________________________________

        vectotsaz_one_order ind(0, 0);
        for (int i = 0; i < size_matrix; i++) {
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

        double* thr = new double[size_matrix];
        for (int i = 0; i < size_matrix; i++) {
            thr[i] = t_spikes[i][max_spikes - 1] + phi_vec[i] * (t_spikes[i][max_spikes - 1] - t_spikes[i][max_spikes - 2]);
        }

        vectotsaz_one_order ind2(0, 0);
        for (int i = 0; i < size_matrix; i++) {
            if (num_spikes[i] == (max_spikes - 1) && t0_init > thr[i] & t0_init_old <= thr[i]) {
                ind2.push_back(i);
            }
        }
        int L2 = ind2.size();

        //__________________________________________________________________add values in wb_init ____________________________
        for (int i = 0; i < L2; i++) {
            int k = ind2[i];
            V[k] = (v_init_old[k] * (t0_init - thr[k]) + v_init[k] * (thr[k] - t0_init_old)) / dt;
            H[k] = (h_init_old[k] * (t0_init - thr[k]) + h_init[k] * (thr[k] - t0_init_old)) / dt;
            N[k] = (n_init_old[k] * (t0_init - thr[k]) + n_init[k] * (thr[k] - t0_init_old)) / dt;
            done[ind2[i]] = 1;
        }
        delete[] thr;
    }//*/
    //p1p(size_matrix, V);
    m_inf_i(M, V);


    delete[] phi_vec;
    done.clear();
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
