# 3.Models_of_population_Neurons
Basic Network Simulations

__________________________________________________________________
--------------------------------------------------------------
# 3.Wang-Buzsáki (Sync volt & spike) (DAG & etc.)
## Wang-Buzsáki model with Runge-Kutta 4th Order Method for 200 neurons in DAG and etc. network with cupling
<p align="center">
 <img src="https://github.com/aliseif321/3.Models_of_population_Neurons/blob/main/3.Wang-Buzs%C3%A1ki_____(Sync%20volt%20&%20spike)___(DAG%20&%20etc.)/Pic/phi=8.2%20%20-%20%20%20I=0.3%20%20-%20%20%20K=1.25%20%20%20-%20%20%20dt=0.05.png?raw=true" >
 </p>
 
_______________________________________________________
```ruby
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


``` 
__________________________________________________________________
--------------------------------------------------------------

# 1.all-to-all 100 WB neurons
## Spike rastergram of 100 synaptically coupled WB neurons.

Figure 31.3 shows an ING rhythm generated by 100 WB neurons coupled with inhibitory synapses. Note that it takes longer to reach synchrony than in, for instance, the PING network of Fig. 30.4. In Fig. 31.3, conditions for synchronization are ideal: There is no heterogeneity in external drives, synaptic connectivity is allto-all, and all synapses have the same strength. This is why perfect synchrony is reached in the limit as t → ∞. With a modest level of drive heterogeneity (σI = 0.03), it takes even longer to reach (approximate) synchrony; see Fig. 31.4. A similar effect is seen when 15% of the synaptic connections in Fig. 31.3 are omitted at random (pII = 0.85), and the remaining ones strengthened by the factor 100/85; see Fig. 31.5. As in Chapter 30, this is not primarily an effect of sparseness and randomness per se, but of variations in the total amount of synaptic input per cell; see Fig. 31.6. With greater heterogeneity, the rhythm disappears; see Fig. 31.7, and also exercise 3. For an analysis of the sensitivity of ING rhythms to heterogeneity, see.

<p align="center">
 <img src="https://github.com/aliseif321/3.Models_of_population_Neurons/blob/main/1.____(WB)____all-to-all/Book/Untitled.png?raw=true" >
 </p>
 
 ## C++
 
 <p align="center">
 <img src="https://github.com/aliseif321/3.Models_of_population_Neurons/blob/main/1.____(WB)____all-to-all/C++/Picture/Untitled.png?raw=true" >
 </p>


 ## MATLAB
 
  <p align="center">
 <img src="https://github.com/aliseif321/3.Models_of_population_Neurons/blob/main/1.____(WB)____all-to-all/MATLAB/Picture/Untitled.png?raw=true" >
 </p>



------------------------------------------------------
_______________________________________________________

__________________________________________________________________
--------------------------------------------------------------

# 1.5.all-to-all 100 WB neurons with Runge-kutaa 4th order method
## Spike rastergram of 100 synaptically coupled WB neurons with Runge-kutaa 4th order method.

 
  <p align="center">
 <img src="https://github.com/aliseif321/3.Models_of_population_Neurons/blob/main/1.5.____(WB)____all-to-all___RK4/C++/Picture/Untitled.png?raw=true" >
 </p>
 
 ## voltage trace of 100 synaptically coupled WB neurons with Runge-kutaa 4th order method.

 
  <p align="center">
 <img src="https://github.com/aliseif321/3.Models_of_population_Neurons/blob/main/1.5.____(WB)____all-to-all___RK4/C++/Picture/Untitled%20(Recovered).png?raw=true" >
 </p>
