#ifndef MOL2_H
#define MOL2_H

#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<string.h>
#include<stdlib.h>
#include<stdio.h>
#include<cmath>

using namespace std;

class Mol2 {
public:


    string molname;
    int N ;
    int Natomtypes;
    int Nbonds;
    int Nres;
    vector<double>charges;
    vector<double>masses;
    vector<string>gaffatoms;
    vector<string> atomnames;
    vector<vector<double> > xyz;
    vector<string> resnames;
    vector< vector< string> > bonds;
    Mol2(string molfile);


};

#endif // MOL2_H
