/************************************************************************************************/
/*** Topic: Wang-Buzsaki model with Runge-Kutta 4th Order Method for 200 neurons              ***/
/*** in DAG and etc. network with cupling                                           Ali-Seif  ***/
/*** Version Release 17.12 rev 11256                                                          ***/
/*** Date: 5/3/2021                                                                           ***/
/*** Code implemented in Microsoft Visual Studio Enterprise 2019 C++ compiler                 ***/
/*** MSI: PX60 6QD/ DDR4                                                                      ***/
/*** Run under a Intel® Core™ i7-6700HQ CPU @ 2.60GHz × 64 based processor with 16 GB RAM     ***/
/************************************************************************************************/
#include <iostream>                                                             //import library for cout and endl and etc.
#include"WangBuzsaki.h"                                                         //import library WangBuzsaki
                                                                                //
using namespace std;                                                            //using standard 
                                                                                //
int main() {                                                                    //
    clock_t t = clock();                                                        //time start record
                                                                                //
    const string address_matrix = "matrix-SF-DAG1src-2";                        //name matrix file
    const int size_matrix = 200;                                                //number of node size of matrix clom and row
    const double dt = 0.05;                                                     //step lenth
    const double t_final = 2500;                                                //time final
    const double coupling = 2.5;                                                //stringh cupling
    double phi = 8.2;                                                           // phi
    double I = 0.3;                                                             //mean I external
                                                                                //
    ofstream syncrony_spike("syncrony spike.txt", ios::out | ios::trunc);       //create file .txt for save sync spike
    ofstream syncrony_voltage("syncrony voltage.txt", ios::out | ios::trunc);   //create file .txt for save sync volt
                                                                                //
    WangBuzsaki WB(size_matrix, dt, phi, I, t_final, coupling, address_matrix); //call class and input initial variable
    WB.rand_borger_200();                                                       //read random constant from .txt
    WB.connection();                                                            //read matrix DAG or etc. network 
    //for (WB.phi = 4.5; WB.phi < 9; WB.phi += 0.1) {                           //change phi paramiter from 4.5 to 9 with lenth step 0.1
    //  for (WB.I_mean = 0; WB.I_mean < 3; WB.I_mean += 0.1) {                  //change I paramiter from 0 to 3 with lenth step 0.1
            WB.initial();                                                       //reset initial paramiter after each run
            WB.Run();                                                           //run wang modle and calculate sync
            syncrony_spike << WB.spikeSynchrony << '\t';                        //print spike sync in file
            syncrony_voltage << WB.voltageSynchrony << '\t';                    //print volt sync in file
            cout << "-------------------------------->\tphi=" +                 //print information on the comand window
                to_string(WB.phi) + "\t-\tI=" +                                 //
                to_string(WB.I_mean)+ "\tspike = " +                            //
                to_string(WB.spikeSynchrony) + "\t - \tvolt = " +               //
                to_string(WB.voltageSynchrony) << endl;                         //
    //    }                                                                     //
    //    syncrony_spike << endl;                                               //press inter in file .txt
    //    syncrony_voltage << endl;                                             //press inter in file .txt
    //}                                                                         //
    cout << "\nTime taken by program is :\t" << ((double)clock() - t) / CLOCKS_PER_SEC << " sec\nfinish\n";
    return 0;                                                                   //run program was correct
}
