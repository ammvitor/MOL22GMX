#ifndef GAFFPARSER_H
#define GAFFPARSER_H
#include <iostream>
#include<vector>
#include<iostream>
#include<fstream>
#include<vector>
#include<stdio.h>
#include<stdlib.h>
#include<string>
#include<time.h>
#include <stdio.h>
#include <string.h>
#include <memory>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdio.h>
#include <string.h>


using namespace std;

class gaffparser
{
public:
    gaffparser(string gaff);
    vector< string> atommass;
    vector< double> masses;
    vector< string>  atom;
    vector < vector< string> > dihedral;
    vector < vector < string> > angles;
    vector < vector < string> > bonds;
    vector < vector < string> > improper;


    vector < string > atoms_id;

    vector < vector < double> > atom_parm;

    vector < vector < double> > dihedral_parm;
    vector < vector < double> > angles_parm;
    vector < vector < double> > bond_parm;
    vector < vector < double> > improper_parm;
    vector < double> GetParm_bonds(string atom_a, string atom_b);
    vector < double> GetParm_angles(string atom_a, string atom_b, string atom_c);
    vector < double> GetParm_dihedrals(string atom_a, string atom_b, string atom_c, string atom_d);
    vector < double> GetParm_impropers( string atom_a, string atom_b, string atom_c, string atom_d);
    vector < double> GetParm_atom(string atom_a);
    double GetParm_masses(string atom_a);
     vector < int > reorg_impropers( string atom_a, string atom_b, string atom_c, string atom_d, int index_a, int index_b,int index_c,int index_d);
};

#endif // GAFFPARSER_H
