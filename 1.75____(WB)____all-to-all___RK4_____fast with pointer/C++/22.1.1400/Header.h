#pragma once
#include<iostream>                                                              //for cout
#include"Matrixreader.h"

using namespace std;
#define size_matrix 100                                                         //size of matrix clom and row
#define  dt  0.01
#define address_matrix "matrix-SF-DAG1src-13.txt"


char address[]{ address_matrix };
//############################################################################################################################
void gii_phi_vec(double** g_ii, double* phi_vec) {

    double
        k_ii = 2.5,
        g = 0.1;

    double** matrix_file = new double* [size_matrix];//created 2d pointer
    for (int i = 0; i < size_matrix; i++) {
        matrix_file[i] = new double[size_matrix];
    }

    matrix(size_matrix, address, matrix_file);//read file with address and size and puting up in 2d pointer//

    srand((unsigned)time(NULL));
    double gkn = (g * k_ii) / (size_matrix);

    for (int i = 0; i < size_matrix; i++) {
        double data;
        data = (double)rand() / RAND_MAX;
        //phi_vec.push_back(data);
        phi_vec[i] = data;
        for (int j = 0; j < size_matrix; j++) {
            g_ii[i][j] = gkn * matrix_file[i][j];
        }
    }
    //p2p(5, 4, matrix_file);//show 2d pointer
    delete[] matrix_file;
}
//############################################################################################################################
double tau_peak_function(double tau_d, double tau_r, double tau_d_q) {

    //double dt = 0.01;
    double dt05 = dt / 2;
    double s = 0;
    double t = 0;
    double s_inc = exp(-t / tau_d_q) * (1 - s) / tau_r - s * tau_d;
    double t_old;
    double s_inc_old;
    double s_tmp;
    double s_inc_tmp;
    while (s_inc > 0) {
        t_old = t;
        s_inc_old = s_inc;
        s_tmp = s + dt05 * s_inc;
        s_inc_tmp = exp(-(t + dt05) / tau_d_q) * (1 - s_tmp) / tau_r - s_tmp / tau_d;
        s = s + dt * s_inc_tmp;
        t = t + dt;
        s_inc = exp(-t / tau_d_q) * (1 - s) / tau_r - s / tau_d;
    }
    return   (t_old * (-s_inc) + t * s_inc_old) / (s_inc_old - s_inc);
}
//##############################################################
double tau_dq_function(double tau_d, double tau_r, double tau_hat) {
    double tau_d_q_left = 1.0;
    while (tau_peak_function(tau_d, tau_r, tau_d_q_left) > tau_hat) {
        tau_d_q_left = tau_d_q_left / 2;
    }

    double tau_d_q_right = tau_r;
    while (tau_peak_function(tau_d, tau_r, tau_d_q_right) < tau_hat) {
        tau_d_q_right = tau_d_q_right * 2;
    }
    double tau_d_q_mid;

    while (tau_d_q_right - tau_d_q_left > pow(10, -12)) {
        tau_d_q_mid = (tau_d_q_left + tau_d_q_right) / 2;
        if (tau_peak_function(tau_d, tau_r, tau_d_q_mid) <= tau_hat) {
            tau_d_q_left = tau_d_q_mid;

        }
        else {
            tau_d_q_right = tau_d_q_mid;

        }
    }
    return (tau_d_q_left + tau_d_q_right) / 2;
}
//############################################################################################################################
void rand_borger(double* phi_vec) {

    ifstream infile;
    infile.open("datarand1.txt");
    if (infile.fail())
    {
        cout << "Could not open file numbers." << "\n";
        //return 1;
    }else{
        double data;
        infile >> data;
        int i = 0;
        while (!infile.eof()) {
            phi_vec[i] = data;
            i++;
            //phi_vec.push_back(data);
            infile >> data;
        }
    }


}
//############################################################################################################################
void  n_inf_i(double* n_init, double* v_i) {
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
//##############################################################

void  h_inf_i(double* h_init, double* v_i) {
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
//##############################################################

void  m_inf_i(double* m_init, double* v_i) {
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


//__________________________________Runge-Kutta calculations_________________________________//

void  rk4_vinit(double* v, double* n, double* h, double* m,double I, double v_rev_i) {


    for (int i = 0; i < size_matrix; i++) {

        int   g_Na_i = 35;
        int   E_Na_i = 55;
        int   g_K = 9;
        int   E_K = -90;
        double g_l = 0.1;
        int   E_l = -65;


        double  k1, k2, k3, k4;
        k1 = dt *  (I - ((g_Na_i * h[i] * pow(m[i], 3) * (v[i] - E_Na_i)) + (g_K * pow(n[i], 4) * (v[i] - E_K)) + (g_l * (v[i] - E_l))));
        k2 = dt * (I - ((g_Na_i * h[i] * pow(m[i], 3) * ((v[i] + k1 / 2) - E_Na_i)) + (g_K * pow(n[i], 4) * ((v[i] + k1 / 2) - E_K)) + (g_l * ((v[i] + k1 / 2) - E_l))));
        k3 = dt * (I - ((g_Na_i * h[i] * pow(m[i], 3) * ((v[i] + k2 / 2) - E_Na_i)) + (g_K * pow(n[i], 4) * ((v[i] + k2 / 2) - E_K)) + (g_l * ((v[i] + k2 / 2) - E_l))) );
        k4 = dt * (I - ((g_Na_i * h[i] * pow(m[i], 3) * ((v[i] + k3) - E_Na_i)) + (g_K * pow(n[i], 4) * ((v[i] + k3) - E_K)) + (g_l * ((v[i] + k3) - E_l))));
        v[i] = v[i] + double((1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4));
    }
}



void  rk4_n(double* v , double* n) {

    double* tau_n = new double[size_matrix];
    double* n_inf = new double[size_matrix];
    for (int i = 0; i < size_matrix; i++) {
        double alpha_n =-0.01 * (v[i] + 34.0) / (exp(-0.1 * (v[i] + 34.0)) - 1.0);
        double beta_n =0.125 * exp(-(v[i] + 44.0) / 80.0);
        double sum = alpha_n + beta_n;
        double tau_nn = ( 1.0 / sum);
        int phi = 5;
        tau_n[i] = (tau_nn / phi);
        n_inf[i]=alpha_n / sum;

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

void  rk4_h(double* v, double* h) {

    double* tau_h = new double[size_matrix];
    double* h_inf = new double[size_matrix];
    for (int i = 0; i < size_matrix; i++) {
        double alpha_h = 0.07 * exp(-(v[i] + 58.0) / 20.0);
        double beta_h = 1.0 / (exp(-0.1 * (v[i] + 28.0)) + 1.0);
        double sum = alpha_h + beta_h;
        double tau_hh = (1.0 / sum);
        int phi = 5;
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



void  rk4_q(double* v, double* q, double tau_dq_i) {

    for (int i = 0; i < size_matrix; i++) {
        double  k1, k2, k3, k4;
        k1 = dt * (0.5 * (1 + tanh(0.1 * v[i])) * (1.0 - q[i]) * 10.0 - (q[i] / tau_dq_i));
        k2 = dt * (0.5 * (1 + tanh(0.1 * v[i])) * (1.0 - (q[i] + k1 / 2)) * 10.0 - ((q[i] + k1 / 2) / tau_dq_i));
        k3 = dt * (0.5 * (1 + tanh(0.1 * v[i])) * (1.0 - (q[i] + k2 / 2)) * 10.0 - ((q[i] + k2 / 2) / tau_dq_i));
        k4 = dt * (0.5 * (1 + tanh(0.1 * v[i])) * (1.0 - (q[i] + k3)) * 10.0 - ((q[i] + k3) / tau_dq_i));
        q[i] = q[i] + double((1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4));
    }
}



void  rk4_s(double* v, double* q, double* s, double tau_r_i, double tau_d_i) {
    for (int i = 0; i < size_matrix; i++) {
        double  k1, k2, k3, k4;

        k1 = dt * (q[i] * (1.0 - s[i]) / tau_r_i - s[i] / tau_d_i);
        k2 = dt * (q[i] * (1.0 - (s[i] + k1 / 2)) / tau_r_i - (s[i] + k1 / 2) / tau_d_i);
        k3 = dt * (q[i] * (1.0 - (s[i] + k2 / 2)) / tau_r_i - (s[i] + k2 / 2) / tau_d_i);
        k4 = dt * (q[i] * (1.0 - (s[i] + k3)) / tau_r_i - (s[i] + k3) / tau_d_i);
        s[i] = s[i] + double((1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4));



    }
}


void rk4_v(double** g_ii, double* v, double* V_old, double* h, double* n, double* m, double* s, double I, double v_rev_i) {

    //double* gs = new double[size_matrix];




    for (int i = 0; i < size_matrix; i++) {

        V_old[i] = v[i];

        double gs=0;
        //gs[i] = 0;
        for (int j = 0; j < size_matrix; j++) {

            gs = gs + g_ii[i][j] * s[j];
        }

        int   g_Na_i = 35;
        int   E_Na_i = 55;
        int   g_K = 9;
        int   E_K = -90;
        double g_l = 0.1;
        int   E_l = -65;


        double  k1, k2, k3, k4;
        k1 = dt * (I - ((g_Na_i * h[i] * pow(m[i], 3) * (v[i] - E_Na_i)) + (g_K * pow(n[i], 4) * (v[i] - E_K)) + (g_l * (v[i] - E_l))) + (gs * (v_rev_i - v[i])));
        k2 = dt * (I - ((g_Na_i * h[i] * pow(m[i], 3) * ((v[i] + k1 / 2) - E_Na_i)) + (g_K * pow(n[i], 4) * ((v[i] + k1 / 2) - E_K)) + (g_l * ((v[i] + k1 / 2) - E_l))) + (gs * (v_rev_i - (v[i] + k1 / 2))));
        k3 = dt * (I - ((g_Na_i * h[i] * pow(m[i], 3) * ((v[i] + k2 / 2) - E_Na_i)) + (g_K * pow(n[i], 4) * ((v[i] + k2 / 2) - E_K)) + (g_l * ((v[i] + k2 / 2) - E_l))) + (gs * (v_rev_i - (v[i] + k2 / 2))));
        k4 = dt * (I - ((g_Na_i * h[i] * pow(m[i], 3) * ((v[i] + k3) - E_Na_i)) + (g_K * pow(n[i], 4) * ((v[i] + k3) - E_K)) + (g_l * ((v[i] + k3) - E_l))) + (gs * (v_rev_i - (v[i] + k3))));
        v[i] = v[i] + double((1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4));
    }

    //cout << "hi";
}





void rk4_run_step(double** g_ii, double* v, double* V_old, double* h, double* n, double* m, double* s, double* q, double I, double v_rev_i, double tau_dq_i, double tau_r_i, double tau_d_i) {

    double* tau_n = new double[size_matrix];
    double* n_inf = new double[size_matrix];
    double* tau_h = new double[size_matrix];
    double* h_inf = new double[size_matrix];

    for (int i = 0; i < size_matrix; i++) {

        V_old[i] = v[i];

        double gs = 0;
        //gs[i] = 0;
        for (int j = 0; j < size_matrix; j++) {

            gs = gs + g_ii[i][j] * s[j];
        }

        int   g_Na_i = 35;
        int   E_Na_i = 55;
        int   g_K = 9;
        int   E_K = -90;
        double g_l = 0.1;
        int   E_l = -65;


        double  k1, k2, k3, k4;
        k1 = dt * (I - ((g_Na_i * h[i] * pow(m[i], 3) * (v[i] - E_Na_i)) + (g_K * pow(n[i], 4) * (v[i] - E_K)) + (g_l * (v[i] - E_l))) + (gs * (v_rev_i - v[i])));
        k2 = dt * (I - ((g_Na_i * h[i] * pow(m[i], 3) * ((v[i] + k1 / 2) - E_Na_i)) + (g_K * pow(n[i], 4) * ((v[i] + k1 / 2) - E_K)) + (g_l * ((v[i] + k1 / 2) - E_l))) + (gs * (v_rev_i - (v[i] + k1 / 2))));
        k3 = dt * (I - ((g_Na_i * h[i] * pow(m[i], 3) * ((v[i] + k2 / 2) - E_Na_i)) + (g_K * pow(n[i], 4) * ((v[i] + k2 / 2) - E_K)) + (g_l * ((v[i] + k2 / 2) - E_l))) + (gs * (v_rev_i - (v[i] + k2 / 2))));
        k4 = dt * (I - ((g_Na_i * h[i] * pow(m[i], 3) * ((v[i] + k3) - E_Na_i)) + (g_K * pow(n[i], 4) * ((v[i] + k3) - E_K)) + (g_l * ((v[i] + k3) - E_l))) + (gs * (v_rev_i - (v[i] + k3))));
        v[i] = v[i] + double((1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4));



            double alpha_n = -0.01 * (v[i] + 34.0) / (exp(-0.1 * (v[i] + 34.0)) - 1.0);
            double beta_n = 0.125 * exp(-(v[i] + 44.0) / 80.0);
            double sum_n = alpha_n + beta_n;
            double tau_nn = (1.0 / sum_n);
            int phi = 5;
            tau_n[i] = (tau_nn / phi);
            n_inf[i] = alpha_n / sum_n;

           // double  k1, k2, k3, k4;

            k1 = dt * ((n_inf[i] - n[i]) / tau_n[i]);
            k2 = dt * ((n_inf[i] - (n[i] + k1 / 2)) / tau_n[i]);
            k3 = dt * ((n_inf[i] - (n[i] + k2 / 2)) / tau_n[i]);
            k4 = dt * ((n_inf[i] - (n[i] + k3)) / tau_n[i]);
            n[i] = n[i] + double((1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4));


            double alpha_h = 0.07 * exp(-(v[i] + 58.0) / 20.0);
            double beta_h = 1.0 / (exp(-0.1 * (v[i] + 28.0)) + 1.0);
            double sum_h = alpha_h + beta_h;
            double tau_hh = (1.0 / sum_h);
            //int phi = 5;
            tau_h[i] = (tau_hh / phi);
            h_inf[i] = alpha_h / sum_h;

           // double  k1, k2, k3, k4;

            k1 = dt * ((h_inf[i] - h[i]) / tau_h[i]);
            k2 = dt * ((h_inf[i] - (h[i] + k1 / 2)) / tau_h[i]);
            k3 = dt * ((h_inf[i] - (h[i] + k2 / 2)) / tau_h[i]);
            k4 = dt * ((h_inf[i] - (h[i] + k3)) / tau_h[i]);
            h[i] = h[i] + double((1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4));



            double alpha_m = 0.1 * (v[i] + 35.0) / (1.0 - exp(-(v[i] + 35.0) / 10.0));
            double beta_m = 4.0 * exp(-(v[i] + 60.0) / 18.0);
            m[i] = alpha_m / (alpha_m + beta_m);


            k1 = dt * (0.5 * (1 + tanh(0.1 * v[i])) * (1.0 - q[i]) * 10.0 - (q[i] / tau_dq_i));
            k2 = dt * (0.5 * (1 + tanh(0.1 * v[i])) * (1.0 - (q[i] + k1 / 2)) * 10.0 - ((q[i] + k1 / 2) / tau_dq_i));
            k3 = dt * (0.5 * (1 + tanh(0.1 * v[i])) * (1.0 - (q[i] + k2 / 2)) * 10.0 - ((q[i] + k2 / 2) / tau_dq_i));
            k4 = dt * (0.5 * (1 + tanh(0.1 * v[i])) * (1.0 - (q[i] + k3)) * 10.0 - ((q[i] + k3) / tau_dq_i));
            q[i] = q[i] + double((1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4));




            k1 = dt * (q[i] * (1.0 - s[i]) / tau_r_i - s[i] / tau_d_i);
            k2 = dt * (q[i] * (1.0 - (s[i] + k1 / 2)) / tau_r_i - (s[i] + k1 / 2) / tau_d_i);
            k3 = dt * (q[i] * (1.0 - (s[i] + k2 / 2)) / tau_r_i - (s[i] + k2 / 2) / tau_d_i);
            k4 = dt * (q[i] * (1.0 - (s[i] + k3)) / tau_r_i - (s[i] + k3) / tau_d_i);
            s[i] = s[i] + double((1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4));
    }

    delete[] tau_n;
    delete[] n_inf;
    delete[] tau_h;
    delete[] h_inf;

    //cout << "hi";
}