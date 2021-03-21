/************************************************************************************************/
/*** Topic: Voltage trace of a WB neuron with an inhibitory autapse.                          ***/
/*** Version Release 17.12 rev 11256                                                          ***/
/*** Date: 3/10/2021                                                               Ali-Seif   ***/
/*** Code implemented in Microsoft Visual Studio Enterprise 2019 C++ compiler                 ***/
/*** MSI: PX60 6QD/ DDR4                                                                      ***/
/*** Run under a Intel® Core™ i7-6700HQ CPU @ 2.60GHz × 64 based processor with 16 GB RAM     ***/
/************************************************************************************************/
#include <iostream>
#include <math.h>
#include <fstream>
#include <random>
#include <vector>
#include<numeric> //for accumulate
#include <time.h>
#include <string>

using namespace std;

typedef vector<double> vectotsaz_one_order;
typedef vector <vectotsaz_one_order> vectotsaz_two_order;
//##############################################################
//####                                                      ####
//####               Calculate alpha and betas              ####
//####                                                      ####
//##############################################################
//________________________________________________________________________________________________
//______________________________________I-cell____________________________________________________
//________________________________________________________________________________________________
double alpha_n_i(double v_i) {
    return   -0.01 * (v_i + 34.0) / (exp(-0.1 * (v_i + 34.0)) - 1.0);
}
double beta_n_i(double v_i) {
    return   0.125 * exp(-(v_i + 44.0) / 80.0);
}
double alpha_m_i(double v_i) {
    return  0.1 * (v_i + 35.0) / (1.0 - exp(-(v_i + 35.0) / 10.0));
}
double beta_m_i(double v_i) {
    return  4.0 * exp(-(v_i + 60.0) / 18.0);
}
double alpha_h_i(double v_i) {
    return   0.07 * exp(-(v_i + 58.0) / 20.0);
}
double beta_h_i(double v_i) {
    return   1.0 / (exp(-0.1 * (v_i + 28.0)) + 1.0);
}
//##############################################################
//####                                                      ####
//####      Calculate infinite activation variables         ####
//####                                                      ####
//##############################################################
//________________________________________________________________________________________________
//______________________________________I-cell____________________________________________________
//________________________________________________________________________________________________
double n_inf_i(double v_i) {
    return   alpha_n_i(v_i) / (alpha_n_i(v_i) + beta_n_i(v_i));
}
double h_inf_i(double v_i) {
    return   alpha_h_i(v_i) / (alpha_h_i(v_i) + beta_h_i(v_i));
}
double m_inf_i(double v_i) {
    return   alpha_m_i(v_i) / (alpha_m_i(v_i) + beta_m_i(v_i));
}
double tau_n_i(double v_i) {
    double tau_n = 1.0 / (alpha_n_i(v_i) + beta_n_i(v_i));
    int phi = 5;
    return   (tau_n / phi);
}
double tau_h_i(double v_i) {
    double tau_h = 1.0 / (alpha_h_i(v_i) + beta_h_i(v_i));
    int phi = 5;
    return    (tau_h / phi);
}
//##############################################################
//####                                                      ####
//####             Calculation of currents                  ####
//####                                                      ####
//##############################################################
//________________________________________________________________________________________________
//______________________________________I-cell____________________________________________________
//________________________________________________________________________________________________
double INa_i(double v_i, double h_i, double m_i) {
    int   g_Na_i = 35;
    int   E_Na_i = 55;
    return    g_Na_i * h_i * pow(m_i, 3) * (v_i - E_Na_i);
}
double IK_i(double v_i, double n_i) {
    int   g_K = 9;
    int   E_K = -90;
    return    g_K * pow(n_i, 4) * (v_i - E_K);
}
double Il_i(double v_i) {
    float g_l = 0.1;
    int   E_l = -65;
    return        g_l * (v_i - E_l);
}
double Isyn_i(double v_i, double g_iis_i, double v_rev_i) {
    int   E_rev_i = v_rev_i;
    return        g_iis_i * (E_rev_i - v_i);
}

//##############################################################
//####                                                      ####
//####            Differential Equations                    ####
//####                                                      ####
//##############################################################
//________________________________________________________________________________________________
//______________________________________I-cell____________________________________________________
//________________________________________________________________________________________________
double dvdt_i(double t, double v_i, double n_i, double h_i, double m_i, double g_iis_i , double I_i, double v_rev_i) {
    float I_app_i = I_i;
    float C_m = 1.0;
    return   (1 / C_m) * (I_app_i - (INa_i(v_i, h_i, m_i) + IK_i(v_i, n_i) + Il_i(v_i)) + Isyn_i(v_i, g_iis_i, v_rev_i));
}
double dndt_i(double t, double n_i, double v_i) {
    return   ((n_inf_i(v_i) - n_i) / tau_n_i(v_i));
}
double dhdt_i(double t, double h_i, double v_i) {
    return   ((h_inf_i(v_i) - h_i) / tau_h_i(v_i));
}
double dqdt_i(double t, double q_i, double v_i, double tau_dq_i) {
    return   0.5 * (1 + tanh(0.1 * v_i)) * (1.0 - q_i) * 10.0 - q_i / tau_dq_i;
}
double dsdt_i(double t, double s_i, double q_i, double v_i, double tau_r_i, double tau_d_i) {
    return   q_i * (1.0 - s_i) / tau_r_i - s_i / tau_d_i;
}
//##############################################################
//####                                                      ####
//####               tau_dq_function                        ####
//####                                                      ####
//##############################################################
//________________________________________________________________________________________________
//_____________________________________tau_peak___________________________________________________
//________________________________________________________________________________________________
double tau_peak_function(double tau_d, double tau_r, double tau_d_q) {

    double dt = 0.01;
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
//________________________________________________________________________________________________
//______________________________________tau_dq____________________________________________________
//________________________________________________________________________________________________
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


#include<iomanip>


//_______________________________________________________________________________________\\
//_____________              The principle of the program                   _____________\\
//_____________                                      @                      _____________\\
//_____________           @@       @@       @            @@     @           _____________\\
//_____________           @ @     @ @      @ @       @   @ @    @           _____________\\
//_____________           @  @   @  @     @   @      @   @  @   @           _____________\\
//_____________           @   @@@   @    @@@@@@@     @   @   @  @           _____________\\
//_____________           @    @    @   @       @    @   @    @ @           _____________\\
//_____________           @         @  @         @   @   @     @@           _____________\\
//_______________________________________________________________________________________

int main() {

//__________________________________________________________________variables_________________________________________
    double
        N = 100,//nunber of nodes
        sigma = 0,//standard deviation
        g_hat_ii = 0.5,
        p_ii = 1,
        g_hat_gap = 0,
        p_gap = 1,
        t_final = 500,
        dt = 0.01,
        dt2 = dt / 2,
        v_rev_i = -75.0,
        tau_r_i = 0.5,
        tau_peak_i = 0.5,
        tau_d_i = 9;

//__________________________________________________________________external corents__________________________________
    default_random_engine generator;
    normal_distribution<double> distribution(0, 1);
    vectotsaz_one_order I_app(N, 0);
    for (int i = 0; i < N; i++) {
        I_app[i] = 1.5 * (1 + distribution(generator) * sigma);
        //cout << I_app[i] << endl;
    }

//__________________________________________________________________Process network parameter_________________________
    vectotsaz_two_order g_ii(N, vectotsaz_one_order(N, 0));
    for (int i = 0; i < N; i++) {
        //cout << "i=" << i;
        for (int j = 0; j < N; j++) {
            g_ii[i][j] = (g_hat_ii * p_ii) / (N * p_ii);
            //cout <<"j="<<j <<g_ii[i][j] << '\t';
        }
        //cout << g_ii[i][i] << endl;
        //cout << endl;
    }

//__________________________________________________________________function of tau dqi_______________________________
    double  tau_dq_i = tau_dq_function(tau_d_i, tau_r_i, tau_peak_i); //0.1163;

//__________________________________________________________________created random numbers____________________________
    /*srand((unsigned)time(NULL));
    vectotsaz_one_order phi_vec(0, 0);
    for (int i = 0; i < N; i++) {
        double data;
        data = (double)rand() / RAND_MAX;
        phi_vec.push_back(data);
        //cout << phi_vec[i] << endl;
    }*/

//__________________________________________________________________random numbers of borger pic______________________
    ifstream infile;
    infile.open("datarand1.txt");
    vectotsaz_one_order phi_vec(0, 0);
    if (infile.fail())
    {
        cout << "Could not open file numbers." << "\n";
        return 1;
    }
    double data;
    infile >> data;
    while (!infile.eof()) {
        phi_vec.push_back(data);
        infile >> data;
    }

//MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
//                                  initial conditions (wb_init[i_ext_i][rand(num_i)])                              //
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

//__________________________________________________________________define variables__________________________________
    int 
        max_spikes = 3;

    double
        t_final_init = 2000,
        t0_init = 0.0;

    vector<int> 
        done(N, 0);

    vectotsaz_one_order
        v_init(N, -70),
        n_init(N, n_inf_i(-70)),
        h_init(N, h_inf_i(-70)),
        m_init(N, m_inf_i(-70)),
        num_spikes(N, -1);

    vectotsaz_two_order
        wb_init(N, vectotsaz_one_order(3, 0)),
        t_spikes(N, vectotsaz_one_order(3, 0));
    int x = 0;

//__________________________________________________________________calculate initial condition_______________________
    while (accumulate(done.begin(), done.end(), 0) < N && t0_init<t_final_init) {
        x = x + 1;
//__________________________________________________________________define variables__________________________________
        double
            t0_init_old = t0_init;

        vectotsaz_one_order 
            v_init_old(N, 0), n_init_old(N, 0), h_init_old(N, 0),
            v_inc(N, 0), n_inc(N, 0), h_inc(N, 0),
            v_tmp(N, 0), n_tmp(N, 0), h_tmp(N, 0), m_tmp(N, 0);

        v_init_old = v_init;
        h_init_old = h_init;
        n_init_old = n_init;

//__________________________________________________________________solve the system of Hodgkin-Huxley-like equations using the midpoint method
        for (int i = 0; i < N; i++) {
            v_inc[i] = dvdt_i(t0_init, v_init[i], n_init[i], h_init[i], m_init[i], 0, I_app[i], v_rev_i);
            n_inc[i] = dndt_i(t0_init, n_init[i], v_init[i]);
            h_inc[i] = dhdt_i(t0_init, h_init[i], v_init[i]);

            v_tmp[i] = v_init[i] + dt2 * v_inc[i];
            h_tmp[i] = h_init[i] + dt2 * h_inc[i];
            n_tmp[i] = n_init[i] + dt2 * n_inc[i];
            m_tmp[i] = m_inf_i(v_init[i]);

            v_inc[i] = dvdt_i(t0_init, v_tmp[i], n_tmp[i], h_tmp[i], m_tmp[i], 0, I_app[i], v_rev_i);
            n_inc[i] = dndt_i(t0_init, n_tmp[i], v_tmp[i]);
            h_inc[i] = dhdt_i(t0_init, h_tmp[i], v_tmp[i]);

            v_init[i] = v_init[i] + dt * v_inc[i];
            h_init[i] = h_init[i] + dt * h_inc[i];
            n_init[i] = n_init[i] + dt * n_inc[i];
            m_init[i] = m_inf_i(v_init[i]);
        }

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

        vectotsaz_one_order thr(N, 0);
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
            wb_init[k][0] = (v_init_old[k] * (t0_init - thr[k]) + v_init[k] * (thr[k] - t0_init_old)) / dt;
            wb_init[k][1] = (h_init_old[k] * (t0_init - thr[k]) + h_init[k] * (thr[k] - t0_init_old)) / dt;
            wb_init[k][2] = (n_init_old[k] * (t0_init - thr[k]) + n_init[k] * (thr[k] - t0_init_old)) / dt;
            done[ind2[i]] = 1;
        }
    }

    for (int i = 0; i < N; i++) {
        v_init[i] = wb_init[i][0];
        h_init[i] = wb_init[i][1];
        n_init[i] = wb_init[i][2];
        m_init[i] = m_inf_i(v_init[i]);
    }


    done.clear();
    num_spikes.clear();
    wb_init.clear();
    t_spikes.clear();
 //MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
 //                     solve the system of Hodgkin-Huxley-like equations using the midpoint method                  //
 //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

    vectotsaz_one_order
        v(N, 0), n(N, 0), h(N, 0), m(N, 0), q(N, 0), s(N, 0);

    v = v_init;
    n = n_init;
    h = h_init;
    m = m_init;

    v_init.clear();
    n_init.clear();
    h_init.clear();
    m_init.clear();

    ofstream temp("temp.txt", ios::out | ios::trunc);

//__________________________________________________________________Solve the system using the midpoint method________
    int num_spikes_i= 0;
    vectotsaz_one_order i_i_spikes(0, 0);
    vectotsaz_one_order t_i_spikes(0, 0);
    vectotsaz_one_order range(0, 0);
    int k = 0;

    for (register double t0 = 0; t0 <= t_final-dt; t0 = t0 + dt) {

        k = k + 1;
        vectotsaz_one_order
            v_inc(N, 0), n_inc(N, 0), h_inc(N, 0), q_inc(N, 0), s_inc(N, 0),
            v_tmp(N, 0), n_tmp(N, 0), h_tmp(N, 0), m_tmp(N, 0), q_tmp(N, 0), s_tmp(N, 0);
            
//__________________________________________________________________solve the system of Hodgkin-Huxley-like equations using the midpoint method
        double gs;
        for (int i = 0; i < N; i++) {
            gs = 0;
            for (int j = 0; j < N; j++) {
                gs=gs+g_ii[i][j] * s[j];
            }
            v_inc[i] = dvdt_i(t0, v[i], n[i], h[i], m[i], gs, I_app[i], v_rev_i);
            n_inc[i] = dndt_i(t0, n[i], v[i]);
            h_inc[i] = dhdt_i(t0, h[i], v[i]);
            q_inc[i] = dqdt_i(t0, q[i], v[i], tau_dq_i);
            s_inc[i] = dsdt_i(t0, s[i], q[i],v[i], tau_r_i, tau_d_i);
        }
        for (int i = 0; i < N; i++) {
            v_tmp[i] = v[i] + dt2 * v_inc[i];
            h_tmp[i] = h[i] + dt2 * h_inc[i];
            n_tmp[i] = n[i] + dt2 * n_inc[i];
            m_tmp[i] = m_inf_i(v_tmp[i]);
            q_tmp[i] = q[i] + dt2 * q_inc[i];
            s_tmp[i] = s[i] + dt2 * s_inc[i];
        }
        for (int i = 0; i < N; i++) {
            gs = 0;
            for (int j = 0; j < N; j++) {
                gs = gs + g_ii[i][j] * s_tmp[j];
            }
            v_inc[i] = dvdt_i(t0, v_tmp[i], n_tmp[i], h_tmp[i], m_tmp[i], gs, I_app[i], v_rev_i);
            n_inc[i] = dndt_i(t0, n_tmp[i], v_tmp[i]);
            h_inc[i] = dhdt_i(t0, h_tmp[i], v_tmp[i]);
            q_inc[i] = dqdt_i(t0, q_tmp[i], v_tmp[i], tau_dq_i);
            s_inc[i] = dsdt_i(t0, s_tmp[i], q_tmp[i], v_tmp[i], tau_r_i, tau_d_i);
        }
        vectotsaz_one_order v_old(N, 0);
        v_old = v;
        for (int i = 0; i < N; i++) {
            v[i] = v[i] + dt * v_inc[i];
            h[i] = h[i] + dt * h_inc[i];
            n[i] = n[i] + dt * n_inc[i];
            m[i] = m_inf_i(v[i]);
            q[i] = q[i] + dt * q_inc[i];
            s[i] = s[i] + dt * s_inc[i];
        }

//__________________________________________________________________find spikes_______________________________________
        vectotsaz_one_order which_i(0, 0);
        for (int i = 0; i < N; i++) {
            if (v_old[i] > -20 && v[i] <= -20) {
                which_i.push_back(i);
            }
            
        }
        int L = which_i.size();

//__________________________________________________________________find spikes (index node and time spikes)__________
        if (L > 0) {
            int conternodes = -1;
            for (int i = (num_spikes_i + 1) ; i <= (num_spikes_i + L); i++) {
                double timing = 0;
                conternodes = conternodes + 1;
                i_i_spikes.push_back(which_i[conternodes]+1);
                timing = (((-20 - v[which_i[conternodes]]) * (k - 1) * dt) + ((v_old[which_i[conternodes]] + 20) * k * dt)) / (-v[which_i[conternodes]] + v_old[which_i[conternodes]]);
                t_i_spikes.push_back(timing); 
            }
            num_spikes_i = num_spikes_i + L;
        }
    }

    for (int i = 0; i < i_i_spikes.size(); i++) {
    //cout << setiosflags(ios::fixed) << setprecision(55);
        temp << t_i_spikes[i]<<'\t'<< i_i_spikes[i] << endl;
    }

    temp.close();
    cout << "\nFinish" << endl;
    return 0;
}