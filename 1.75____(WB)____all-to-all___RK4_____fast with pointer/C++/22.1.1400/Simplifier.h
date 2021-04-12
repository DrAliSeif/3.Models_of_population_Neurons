#pragma once
#include<iostream>
using namespace std;

void endl() {
    cout << endl;
}

void p1p(int size, double* a) {
    for (int i = 0; i < size; i++) {
        cout << i<<'\t'<<a[i];
        endl();
    }
}
void p1p(int size, int* a) {
    for (int i = 0; i < size; i++) {
        cout << i << '\t' << a[i];
        endl();
    }
}

void i1p(int size, double* a, double init) {

    for (int i = 0; i < size; i++) {
        a[i]= init;
    }
}
void i1p(int size, int* a, int init) {

    for (int i = 0; i < size; i++) {
        a[i] = init;
    }
}


void i2p(int size, double** a, double init) {

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            a[i][j] = init;
        }
    }

}


void p2p(int size_x, int size_y,double** a) {

    for (int i = 0; i < size_x; i++) {
         for (int j = 0; j < size_y; j++) {
               cout << a[i][j] << '\t';
         }
         endl();
    }

}


void e1p(int size, double* a, double* b) {

    for (int i = 0; i < size; i++) {
        a[i] = b[i];
    }
}
