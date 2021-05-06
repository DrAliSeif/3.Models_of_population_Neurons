/************************************************************************************************/
/*** Topic: Wang-Buzsaki model with Runge-Kutta 4th Order Method for 200 neurons              ***/
/*** in DAG and etc. network with cupling                                           Ali-Seif  ***/
/*** Version Release 17.12 rev 11256                                                          ***/
/*** Date: 5/3/2021                                                                           ***/
/*** Code implemented in Microsoft Visual Studio Enterprise 2019 C++ compiler                 ***/
/*** MSI: PX60 6QD/ DDR4                                                                      ***/
/*** Run under a Intel® Core™ i7-6700HQ CPU @ 2.60GHz × 64 based processor with 16 GB RAM     ***/
/************************************************************************************************/
#include<iostream>                                                              //for cout
#include<vector>                                                                //for vector
#include <fstream>                                                              //infile
#include<numeric>                                                               //for accumulate
#include <sstream>																//to_string
#include <random>                                                               //sort
#include <assert.h>                                                             //assert
using namespace std;                                                            //using standard 
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
class WangBuzsaki {                                                             //create and define class 
private:                                                                        //private values
                                                                                //
    const int N;                                                                //number of node
    const double dt;                                                            //step lenth
    double tfinal;                                                              //time final
    double coupling;                                                            //stringh cupling
    string address_matrix;                                                      //name matrix file
    int num_steps;                                                              //number of spike
    double* randborg = new double[N];                                           //creat pointer1D for save random in .txt file
    double* I_app = new double[N];                                              //distribution corent
    double* IC =new double[4*N];                                                //v & h & n & s
    vec2i column_connection;                                                    //number of column connection in [row][number of connection in one row]
    double couplingOverN = coupling / (N + 0.0);                                //wight copling string in each neuron
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
                                                                                //
    WangBuzsaki                                                                 //
    (int N_, double dt_, double phi_, double I_, double t_final_, double coupling_, string address_matrix_) ://input data in class
        N(N_), dt(dt_), phi(phi_), I_mean(I_), tfinal(t_final_), coupling(coupling_), address_matrix(address_matrix_)//change name input to privet
    {num_steps = int(tfinal / dt);}                                             //number of steps
    double phi;                                                                 // phi
    double I_mean;                                                              //mean I external
    double spikeSynchrony = 0.0;                                                //define variable spike sync
    double voltageSynchrony = 0.0;                                              //define variable volt sync
                                                                                //
    void rand_borger_200();                                                     //read rand rumber from .txt to variable
    void h_inf(int);                                                            //calculate h inf 
    void n_inf(int);                                                            //calculate n inf 
    void initial();                                                             //set initial condition
    void connection();                                                          //read matrix from .txt and calculate cunections
    void runge_kutta4_integrator(double*);                                      //runge kutta 4
    void dydt(double*,const double*);                                           //one step wang for all paramiter
    void Run();                                                                 //run wang for all time 
    double calculate_spike_synchrony(vec1d&);                                   //calculate spike sync
    double calculate_voltage_synchrony(double**);                               //calculate volt sync
    void calculate_synchrony_measures(double**, const vec2d&, double&, double&);//chek cundition that calculate synchrony
};                                                                              //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                                  rand_borger_200                               @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void WangBuzsaki::rand_borger_200() {                                           //propertis of read random generator from .txt file
    ifstream infile;                                                            //create space for read file
    infile.open("datarand200.txt");                                             //read address and file .txt
    if (infile.fail())                                                          //if could not find file
    {                                                                           //
        cout << "Could not open file numbers." << "\n";                         //type this paragraph
    }                                                                           //
    else {                                                                      //if find file
        double data;                                                            //define double for transport .txt file
        infile >> data;                                                         //transport first data from file to variable 
        int i = 0;                                                              //create conter for cont 
        while (!infile.eof()) {                                                 //while ine the file was data do this  
            randborg[i] = data;                                                 //save variable in to the pointer1D randborg
            i++;                                                                //i=i+1
            infile >> data;                                                     //transport next data from file to variable 
        }                                                                       //
    }                                                                           //
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                                    h_inf_i                                     @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void  WangBuzsaki::h_inf(int i) {                                               //propertis of calculate h_inf by V
    double alpha_h = 0.07 * exp(-(IC[i] + 58.0) / 20.0);                        //calculate alpha_h by V
    double beta_h = 1.0 / (exp(-0.1 * (IC[i] + 28.0)) + 1.0);                   //calculate betha_h by V
    IC[N + i] = alpha_h / (alpha_h + beta_h);                                   //calculate initial h for i th node
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                                    n_inf_i                                     @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void WangBuzsaki::n_inf(int i) {                                                //propertis of calculate n_inf by V
        double alpha_n = -0.01 * (IC[i] + 34.0) / (exp(-0.1 * (IC[i] + 34.0)) - 1.0);//calculate alpha_n by V
        double beta_n = 0.125 * exp(-(IC[i] + 44.0) / 80.0);                    //calculate betha_n by V
        IC[2 * N + i] = alpha_n / (alpha_n + beta_n);                           //calculate initial n for i th node
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                                    initial                                     @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void WangBuzsaki::initial() {                                                   //call initial in wangbuzsaki class
    for (int i = 0; i < N; i++) {                                               //loop for from 0 to N 
        IC[i] = (randborg[i] * 20.0 - 70.0);                                    //change random number to between -70 to -50
        h_inf(i);                                                               //call Function h_inf
        n_inf(i);                                                               //call Function n_inf
        IC[3 * N + i] = 0;                                                      //initial condition for S
        I_app[i] = I_mean;                                                      //all corrent was I
    }                                                                           //
    cout << "initial conditions done--------->"<< endl;                         //type this sentence
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                          Read matrix connection                                @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void matrix(int size_matrix, const std::string address_matrix, double** a) {    //read matrix connection 
    vector_2d_int matrix_2d;                                                    //create vector integer two dimensional(matrix_2d)
    std::ifstream ifile(address_matrix + ".txt");                               //read address of file .txt
    if (ifile.is_open()) {                                                      //if file was available 
        vector_1d_int matrix_1d;                                                //create vector integer one dimensional(matrix_1d)
        int num;                                                                //create integer number(num)
        while (ifile >> num) {                                                  //Set the read number of the file to the defined integer(num)
            matrix_1d.push_back(num);                                           //set defined integer(num) into the One after the last cell vector integer one dimensional(matrix_1d)
            if (matrix_1d.size() == size_matrix) {                              //if size of vector is full
                matrix_2d.push_back(matrix_1d);                                 //set vector integer one dimensional(matrix_1d) into the One after the last cell vector integer tow dimensional(matrix_2d)
                matrix_1d.clear();                                              //clean all cels of vector integer one dimensional(matrix_1d) and delete it
            }                                                                   //
        }                                                                       //
    }                                                                           //
    else {                                                                      //else if file wasn't available
        std::cout << "There was an error opening the input file!\n";            //print "..."
        exit(1);                                                                //means program(process) terminate normally unsuccessfully..
    }                                                                           //
    for (int i = 0; i < size_matrix; i++) {                                     //loop for from 0 to 200 for row
        for (int j = 0; j < size_matrix; j++) {                                 //loop for from 0 to 200 for clumn
            a[i][j] = matrix_2d[i][j];                                          //move data matrix in to the  matrix_file pointer 2D
        }                                                                       //
    }                                                                           //
    matrix_2d.clear();                                                          //remove vectore
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                            column_connection                                   @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void WangBuzsaki::connection()                                                  //propertis of read matrix connection from .txt file and calculate number of connection in each row
{                                                                               //
    vec2i g_ii(N, vec1i(N));                                                    //created 2d vectore
    double** matrix_file = new double* [N];                                     //created 2d pointer
    for (int i = 0; i < N; i++) {                                               //
        matrix_file[i] = new double[N];                                         //
    }                                                                           //
    matrix(N, address_matrix, matrix_file);                                     //read matrix for example 200 dag[satr][soton] with address and size and puting up in 2d pointer//
    for (int i = 0; i < N; i++) {                                               //loop for from 0 to 200 for row
        for (int j = 0; j < N; j++) {                                           //loop for from 0 to 200 for clumn
            g_ii[i][j] = matrix_file[i][j];                                     //move data matrix in to the  g_ii vector 2D
        }                                                                       //
    }                                                                           //
    delete[] matrix_file;                                                       //remove pointer
    column_connection.resize(N);                                                //change size vector to 200 
    for (int i = 0; i < N; i++){                                                //loop for from 0 to 200 for row
        for (int j = 0; j < N; j++){                                            //loop for from 0 to 200 for clumn
            if (g_ii[i][j] > 1.e-8) {                                           //if for this connection was 1
                column_connection[i].push_back(j);                              //[satr][chandomin etesal dar satr]=soton 
            }                                                                   //
        }                                                                       //
    }                                                                           //
    g_ii.clear();                                                               //remove vector
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                            spike_synchrony                                     @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
vec1d calculate_isi(const vec1d& v)                                             //function that input is  vectore of time spike and output is vectore of inter spike interval
{                                                                               //
    int n = v.size();                                                           //n is number of spikes
    vec1d v_isi(n - 1);                                                         //creat vectore 1d with size of n-1
    for (int i = 0; i < (n - 1); i++)                                           //loop for from 0 to n-1
        v_isi[i] = v[i + 1] - v[i];                                             //calculate interval
    return v_isi;                                                               //return vector inter spike interval ISI
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                                  mean                                          @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
inline double mean(const vec1d& v, const int id)                                //average the vector from element "id" to end of the vector
{                                                                               //
    assert(v.size() > id);                                                      //Provided that the denominator is not zero
    return accumulate(v.begin() + id, v.end(), 0.0) / (v.size() - id);          //accumulate(first vector,last vector,number of sum since external=0) and sum all number of vector and devided number thats means average
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                         calculate_spike_synchrony                              @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
double WangBuzsaki::calculate_spike_synchrony(vec1d& spiketrains)               //
{                                                                               //
    int n = spiketrains.size();                                                 //number of timespike in interval time
    sort(spiketrains.begin(), spiketrains.end());                               //sort all timespike from The least to the most
    auto tau = calculate_isi(spiketrains);                                      //calculate inter spike interval between each tow timespike(auto means i dont know type of this variable and hey computer Recognize type this)
    vec1d tau2(tau.size());                                                     //define vectore 1d with size number inter spike interval
    for (int i = 0; i < tau.size(); i++)                                        //loop  for from 0 to size of inter spike
        tau2[i] = tau[i] * tau[i];                                              //power 2
    double tau_m = mean(tau, 0);                                                //<tau> 
    double tau2_m = mean(tau2, 0);                                              //<tau2> 
    double tau_m2 = tau_m * tau_m;                                              //<tau>2
    double burst = ((sqrt(tau2_m - tau_m2)) / (tau_m + 0.0) - 1.0) / (sqrt(N) + 0.0);//calculate sync value
    return burst;                                                               //return spike sync number
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                       calculate_voltage_synchrony                              @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
double WangBuzsaki::calculate_voltage_synchrony(double** voltages)              //calculate volt sync
{                                                                               //
    int multiplier = 0.25 / dt;                                                 //for decreas calculation volt sync each 0.25ms calculate volt sync
    int index_transition = (int)(2000 / dt);                                    //number of step that start save and calculate sync (2000ms or 40000th)
    int nstep = int((num_steps - index_transition) / double(multiplier)) - 1;   //time interval of step to save and calculate sync of step that 500ms or 2000th that because first nod was 0 numbe rof steps was 1999
    vec1d vg(nstep);                                                            //create vg vectore 1d with size step 1999 
    for (int j = 0; j < nstep; j++)                                             //for loop from 0 to last step
    {                                                                           //
        for (int i = 0; i < N; i++)                                             //for loop for all node
            vg[j] += voltages[i][j];                                            //sum all voltage for each node
        vg[j] /= (N + 0.0);                                                     //average each node
    }                                                                           //
    double vg_m = mean(vg, 0);                                                  //average for all node average
    double O = 0.0;                                                             //define o=0 initial variable
    for (int j = 0; j < nstep; j++)                                             //count 0 to 1999 step for save data to voltage vector
    {                                                                           //
        double tmp = (vg[j] - vg_m);                                            //diferents between each average node and average all node
        O += (tmp * tmp);                                                       //add squer diferent to O variable 
    }                                                                           //
    O /= (nstep + 0.0);                                                         //average since diferents squer
    vec1d v(nstep);                                                             //define vectore v with size of steps 1999 
    double denom = 0.0;                                                         //define new variable
    for (int i = 0; i < N; i++)                                                 //for loop that cont all node N=200
    {                                                                           //
        for (int j = 0; j < nstep; j++)                                         //for loop that cont step since 0 to 1999
            v[j] = voltages[i][j];                                              //sum all voltage in each node in all time saved
        double m = mean(v, 0);                                                  //average voltage vector for all node
        double sigma = 0;                                                       //define variable sigma
        for (int j = 0; j < nstep; j++)                                         //for loop that cont step since 0 to 1999
        {                                                                       //
            double tmp = (v[j] - m);                                            //diferents between each average node and average all node
            sigma += (tmp * tmp);                                               //add squer diferent to O variable 
        }                                                                       //
        sigma /= (nstep + 0.0);                                                 //average since diferents squer
        denom += sigma;                                                         //sum average for all node
    }                                                                           //
    denom /= (N + 0.0);                                                         //average for all node
    double xi = O / (denom + 0.0);                                              //calculate volt sync
    return xi;                                                                  //return sync volt
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                                  flatten                                       @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
template <typename T>                                                           //define new template that use typename
vector<T> flatten(const vector<vector<T>>& v)                                   //call 2d vector
{                                                                               //
    size_t total_size = 0;                                                      //create initial total size
    for (const auto& sub : v)                                                   //cont all size of 2d vector
        total_size += sub.size();                                               //sum all size in each row
    vector<T> result;                                                           //define vector 1d
    result.reserve(total_size);                                                 //resize vector 1d with total size
    for (const auto& sub : v)                                                   //cont all size of 2d vector
        result.insert(result.end(), sub.begin(), sub.end());                    //add all row in 1d vector
    return result;                                                              //return vectore T result
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                        calculate_synchrony_measures                            @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void WangBuzsaki::calculate_synchrony_measures(double** voltages,const vec2d& spiketrains,double& vol_syn,double& spk_syn)//call variable that save for calculate sync
{                                                                               //
    vec1d flatten_vec = flatten(spiketrains);                                   //compress 2d vector in 1d vector
    if (flatten_vec.size() < 2)                                                 //if in this interval 2000ms to 2500ms had not spike 
    {                                                                           //
        spk_syn = 0.0;                                                          //print and type non sync
        vol_syn = 0.0;                                                          //print and type non sync
    }                                                                           //
    else                                                                        //else if in this interval 2000ms to 2500ms had spike 
    {                                                                           //
        spk_syn = abs(calculate_spike_synchrony(flatten_vec));                  //calculate sync spike and save this number in the spk_sync to print in the file
        vol_syn = calculate_voltage_synchrony(voltages);                        //calculate sync volt and save this number in the vol_sync to print in the file
    }                                                                           //
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                                     RUN                                        @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void WangBuzsaki::Run() {                                                       //call Run function
    double* aux = new double[N];                                                //created 1d pointer for chek spike
    int multiplier = 0.25 / dt;                                                 //for decreas calculation volt sync each 0.25ms calculate volt sync
    int index_transition = (int)(2000 / dt);                                    //number of step that start save and calculate sync (2000ms or 40000th)
    int m = int((num_steps - index_transition) / double(multiplier)) - 1;       //time interval of step to save and calculate sync of step that 500ms or 2000th that because first nod was 0 numbe rof steps was 1999
    double** voltages = new double* [N];                                        //created 2d pointer to put voltage for calculate sync volt
    for (int i = 0; i < N; i++) {                                               //N row
        voltages[i] = new double[m]{ NULL };                                    //m clumn
    }                                                                           //
    cout << address_matrix + " Run start--->"  << endl;                         //type address and Start run
    ofstream temp("temp.txt", ios::out | ios::trunc);                           //create file .txt for save data spike
    vec2d spikes(N);                                                            //create 2d vector for put time spike in this vectore that size was N*1
    double t = 0.0;                                                             //initial time for start run
    int counter = 0;                                                            //conter of spike node
    for (int it = 1; it < num_steps; ++it)                                      //loop fo 1 to 2500*2 or (2500/0.5)
    {                                                                           //
        runge_kutta4_integrator(IC);                                            //calculate runge kutta for all v n h s with k1 to k4 step by step
        t = it * dt;                                                            //calculate time now
        for (int i = 0; i < N; i++)                                             //loop for all nod 200
        {                                                                       //
            if ((IC[i] > THRESHOLD) && aux[i] == 0)                             //if cross threshold and fier spike
            {                                                                   //
                aux[i] = 1;                                                     //for i neuron tell this neuron was upper threshold
                if (t > 2000) {                                                 //if time was upper 2000 ms
                    temp << t << '\t' << i << endl;                             //print this time and nod in file
                    spikes[i].push_back(t);                                     //add time t in vectore spike for neuron i to one row
                }                                                               //
            }                                                                   //
            else if ((IC[i] < THRESHOLD))                                       //if voltage was under threshold
                aux[i] = 0;                                                     //for i neuron tell this neuron was under threshold
            if ((it % multiplier == 0) && (it > index_transition))              //if time t was multiplier of 0.25ms and time was after 2000ms or(after 4000it)
            {                                                                   //
               voltages[i][counter] = IC[i];                                    //add voltage to calculate sync volt [number of node][number of volt]
            }                                                                   //
        }                                                                       //
        if ((it % multiplier == 0) && (it > index_transition))                  //if time t was multiplier of 0.25ms and time was after 2000ms or(after 4000it)
        {                                                                       //after loop for to cont
            counter++;                                                          //move next step for save voltage 
        }                                                                       //
    }                                                                           //
    delete aux;                                                                 //remove chek spike
    spikeSynchrony = 0.0;                                                       //initial valiu befor calculate sync
    voltageSynchrony = 0.0;                                                     //initial valiu befor calculate sync
    calculate_synchrony_measures(voltages,spikes,voltageSynchrony,spikeSynchrony);//calculate sync spike and volt
    delete voltages;                                                            //remove pointer that save voltage to calculate sync volt
    spikes.clear();                                                             //remove vectore that save spike to calculate sync spike
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                            runge_kutta4_integrator                             @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void WangBuzsaki::runge_kutta4_integrator(double* y)                            //with have IC(V h n s) we calculate runge kutta for all parameter
{                                                                               //
    int n4 =4* N;                                                               //size of IC
    double* k1 = new double[n4];                                                //created 1d pointer for k1
    double* k2 = new double[n4];                                                //created 1d pointer
    double* k3 = new double[n4];                                                //created 1d pointer
    double* k4 = new double[n4];                                                //created 1d pointer
    double* f = new double[n4];                                                 //created 1d pointer for change to runge
    for (int i = 0; i < n4; i++) {                                              //for loop from 0 to 4*200 V h n s
        f[i] = y[i];                                                            //each IC[i] save to f[i]
    }                                                                           //
    dydt(k1,f);                                                                 //calculate for all V h n s one step k1
    for (int i = 0; i < n4; i++)                                                //for loop from 0 to 4*200 V h n s
        f[i] = y[i] + 0.5 * dt * k1[i];                                         //calculate gradient k1
    dydt(k2, f);                                                                //calculate for all V h n s one step k2
    for (int i = 0; i < n4; i++)                                                //for loop from 0 to 4*200 V h n s
        f[i] = y[i] + 0.5 * dt * k2[i];                                         //calculate gradient k2
    dydt(k3, f);                                                                //calculate for all V h n s one step k3
    for (int i = 0; i < n4; i++)                                                //for loop from 0 to 4*200 V h n s
        f[i] = y[i] + dt * k3[i];                                               //calculate gradient k3
    dydt(k4, f);                                                                //calculate for all V h n s one step k4
    for (int i = 0; i < n4; i++)                                                //for loop from 0 to 4*200 V h n s
        y[i] += (k1[i] + 2.0 * (k2[i] + k3[i]) + k4[i]) * dt / 6.0;             //calculate gradient total
    delete f;                                                                   //remove pointer
    delete k1;                                                                  //remove pointer
    delete k2;                                                                  //remove pointer
    delete k3;                                                                  //remove pointer
    delete k4;                                                                  //remove pointer
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                                     dydt                                       @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void WangBuzsaki::dydt(double* ff,const double* x0)                             //calculate one step for all parameter v h n s
{                                                                               //
    const int n2 = 2 * N;                                                       //first n
    const int n3 = 3 * N;                                                       //first s
    for (int i = 0; i < N; i++)                                                 //loop for from all node
    {                                                                           //
        double alpha_m = -0.1 * (x0[i] + 35.0) / (exp(-0.1 * (x0[i] + 35.0)) - 1.0);//alpha_m
        double alpha_h = 0.07 * exp(-(x0[i] + 58.0) / 20.0);                    //alpha_h
        double alpha_n = -0.01 * (x0[i] + 34.0) / (exp(-0.1 * (x0[i] + 34.0)) - 1.0);//alpha_n
        double beta_m = 4.0 * exp(-(x0[i] + 60.0) / 18.0);                      //beta_m
        double beta_h = 1.0 / (exp(-0.1 * (x0[i] + 28.0)) + 1.0);               //beta_h
        double beta_n = 0.125 * exp(-(x0[i] + 44.0) / 80.0);                    //beta_n
        double m = alpha_m / (alpha_m + beta_m);                                //m_inf
        double F = 1.0 / (1.0 + exp(-0.5 * x0[i]));                             //F function
        double I_Na = g_Na * m * m * m * x0[i + N] * (x0[i] - E_Na);            //I_Na
        double I_L = g_L * (x0[i] - E_L);                                       //I_L
        double I_K = g_K * x0[i + n2] * x0[i + n2] * x0[i + n2] * x0[i + n2] * (x0[i] - E_K);//I_K
        double I_syn = 0.0;                                                     //I_syn befor calculate reset
        for (int j = 0; j < column_connection[i].size(); j++)                   //for loop for each row with size of number of connection in to the each row
        {                                                                       //
            int k = column_connection[i][j];                                    //number of column connection that means i=output j=input neuron
            I_syn += (g_syn * x0[k + n3] * (x0[k] - E_syn));                    //sum calculation all nodes that connect to i ---> k means node that connect to i
        }                                                                       //
        I_syn = couplingOverN * I_syn;                                          //add Weight Connection
        ff[i] = -I_Na - I_K - I_L - I_syn + I_app[i];                           // dv/dt
        ff[i + N] = phi * (alpha_h * (1 - x0[i + N]) - beta_h * x0[i + N]);     // dh/dt
        ff[i + n2] = phi * (alpha_n * (1 - x0[i + n2]) - beta_n * x0[i + n2]);  // dn/dt
        ff[i + n3] = alpha * F * (1 - x0[i + n3]) - beta * x0[i + n3];          // ds/dt
    }                                                                           //
}                                                                               //