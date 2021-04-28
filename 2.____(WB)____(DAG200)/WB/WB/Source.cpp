#include <iostream>
#include"WangBuzsaki.h"

using namespace std;

int main() {
    clock_t t = clock();                                                        //time start record

    const int size_matrix = 200;                                                //number of node size of matrix clom and row
    const double dt = 0.05;                                                     //step lenth
    const double phi = 6.9;                                                     // phi
    const double I = 0.35;                                                      //mean I external
    const double t_final = 2500;                                                //time final
    const double coupling = 2.5;                                                //stringh cupling
    const string address_matrix = "matrix-SF-DAG1src-2";                        //name matrix file

    WangBuzsaki WB(size_matrix, dt, phi, I, t_final, coupling, address_matrix); //call class and input initial variable
    WB.initial(WB.rand_borger_200());                                           //run 2000ms without coupling for neurons for initial cundition
    WB.Run();

    std::cout << "\nTime taken by program is :\t" << ((double)clock() - t) / CLOCKS_PER_SEC << " sec\n";
    std::cout << "finish\n";
    return 0;
}
