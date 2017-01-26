#include "gaffparser.h"

gaffparser::gaffparser(string gaff)
{

    string buff1, buff2, buff3, buff4;

    vector < string> temp_str;
    vector < double> temp_int;

    float k0, ang0 , period, mass;
    string line ;
    ifstream myfile (gaff.c_str());
    if(myfile){
        while ( line != "[bonds]" )
        {
          getline (myfile,line);
          if(line != "[bonds]"){
             buff1 = line.substr(0,2);

             buff1.erase(remove (buff1.begin(), buff1.end(), ' '), buff1.end());
             this->atommass.push_back(buff1);
             mass = atof(line.substr(3,5).c_str());
             this->masses.push_back(mass);
          }

        }

          //PArsing the bonds
        while ( line != "[angles]" )
        {
          getline (myfile,line);
          if(line != "[bonds]" && line != "[angles]"){
             buff1 = line.substr(0,2);
             ;
             buff1.erase(remove (buff1.begin(), buff1.end(), ' '), buff1.end());
             buff2 = line.substr(3,2);
             buff2.erase(remove(buff2.begin(), buff2.end(), ' '), buff2.end());
             temp_str.push_back(buff1);
             temp_str.push_back(buff2);
             this->bonds.push_back(temp_str);
             temp_str.clear();
             k0 = atof(line.substr(7,6).c_str());
             ang0= atof(line.substr(16,5).c_str());
             temp_int.push_back(k0);
             temp_int.push_back(ang0);
             this->bond_parm.push_back(temp_int);
             temp_int.clear();
          }


        }

        while ( line != "[dihedral]" )
        {
          getline (myfile,line);
          if(line != "[dihedral]" && line != "[angles]"){
              buff1 = line.substr(0,2);
              buff1.erase(remove(buff1.begin(), buff1.end(), ' '), buff1.end());
              buff2 = line.substr(3,2);
              buff2.erase(remove(buff2.begin(), buff2.end(), ' '), buff2.end());
              buff3 = line.substr(6,2);
              buff3.erase(remove(buff3.begin(), buff3.end(), ' '), buff3.end());
             temp_str.push_back(buff1);
             temp_str.push_back(buff2);
             temp_str.push_back(buff3);
             this->angles.push_back(temp_str);
             temp_str.clear();
             k0 = atof(line.substr(11,11).c_str());
             ang0= atof(line.substr(20).c_str());
             temp_int.push_back(k0);
             temp_int.push_back(ang0);
             this->angles_parm.push_back(temp_int);
             temp_int.clear();
          }


        }


           //Parsing dihedrals


         while ( line != "[improper]" )
          {
             getline (myfile,line);
             if(line != "[dihedral]" && line != "[improper]"){
                 buff1 = line.substr(0,2);
                 buff1.erase(remove(buff1.begin(), buff1.end(), ' '), buff1.end());
                 buff2 = line.substr(3,2);
                 buff2.erase(remove(buff2.begin(), buff2.end(), ' '), buff2.end());
                 buff3 = line.substr(6,2);
                 buff3.erase(remove(buff3.begin(), buff3.end(), ' '), buff3.end());
                 buff4 = line.substr(9,2);
                 buff4.erase(remove(buff4.begin(), buff4.end(), ' '), buff4.end());
                 temp_str.push_back(buff1);
                 temp_str.push_back(buff2);
                 temp_str.push_back(buff3);
                 temp_str.push_back(buff4);


                 this->dihedral.push_back(temp_str);
                 temp_str.clear();
                 period = atof(line.substr(14,1).c_str());
                 k0 = atof(line.substr(19,5).c_str());
                 ang0= atof(line.substr(31,7).c_str());

                 temp_int.push_back(period);
                 temp_int.push_back(k0);
                 temp_int.push_back(ang0);
                 this->dihedral_parm.push_back(temp_int);
                 temp_int.clear();
          }
         }


          //improper dihedral
          while ( line != "[end]" )
          {
            getline (myfile,line);
              if(line != "[end]" && line != "[improper]"){
                  buff1 = line.substr(0,2);
                  buff1.erase(remove(buff1.begin(), buff1.end(), ' '), buff1.end());
                  buff2 = line.substr(3,2);
                  buff2.erase(remove(buff2.begin(), buff2.end(), ' '), buff2.end());
                  buff3 = line.substr(6,2);
                  buff3.erase(remove(buff3.begin(), buff3.end(), ' '), buff3.end());
                  buff4 = line.substr(9,2);
                  buff4.erase(remove(buff4.begin(), buff4.end(), ' '), buff4.end());
                  temp_str.push_back(buff1);
                  temp_str.push_back(buff2);
                  temp_str.push_back(buff3);
                  temp_str.push_back(buff4);
                  this->improper.push_back(temp_str);
                  temp_str.clear();
                  k0 = atof(line.substr(19,5).c_str());
                  ang0= atof(line.substr(32,9).c_str());
                  period = atof(line.substr(48,5).c_str());

                  temp_int.push_back(k0);
                  temp_int.push_back(ang0);
                  temp_int.push_back(period);
                  this->improper_parm.push_back(temp_int);
                  temp_int.clear();
              }
          }


          while ( line != "[aftermol]" )
          {
            getline (myfile,line);

              if(line != "[end]" && line != "[aftermol]"){
                  buff1 = line.substr(2,2);
                  buff1.erase(remove(buff1.begin(), buff1.end(), ' '), buff1.end());

                  this->atoms_id.push_back(buff1);
                  k0 = atof(line.substr(14,5).c_str());
                  ang0= atof(line.substr(22,6).c_str());
                  temp_int.push_back(k0);
                  temp_int.push_back(ang0);
                  this->atom_parm.push_back(temp_int);

                  temp_int.clear();
              }
          }

        }
        else{
        printf("Cant Find gaff Parameters files.");
    }

    }

vector < double> gaffparser::GetParm_bonds(string atom_a, string atom_b){
    for(int i =0; i < this->bonds.size(); i++){
        if(bonds[i][0] == atom_a && bonds[i][1] == atom_b){

            return this->bond_parm[i];
        }
        else if(bonds[i][0] == atom_b && bonds[i][1] == atom_a){
            return this->bond_parm[i];
        }

    }

    cout << "Missing bond parameters for :"  << atom_b << " " << atom_a <<  endl;
    exit(0);

}

vector < double> gaffparser::GetParm_angles(string atom_a, string atom_b, string atom_c){
    for(int i =0; i < this->angles.size(); i++){
        if(angles[i][0] == atom_a && angles[i][1] == atom_b && angles[i][2] == atom_c){

            return this->angles_parm[i];
        }
        else if(angles[i][0] == atom_c && angles[i][1] == atom_b && angles[i][2] == atom_a){

            return this->angles_parm[i];

        }
    }
    cout << atom_a << " " << atom_b << " " << atom_c <<  endl;
    cout << "Missing angles parameters for :"  << atom_b << " " << atom_a << " " << atom_c <<  endl;
    exit(0);


}
vector < double> gaffparser::GetParm_dihedrals(string atom_a, string atom_b, string atom_c, string atom_d){


    for(int i =0; i < this->dihedral.size(); i++){

        if(this->dihedral[i][0] == atom_a && this->dihedral[i][1] == atom_b && this->dihedral[i][2] == atom_c && this->dihedral[i][3] == atom_d){

            return this->dihedral_parm[i];
        }

        if(this->dihedral[i][0] == atom_d && this->dihedral[i][1] == atom_c && this->dihedral[i][2] == atom_b && this->dihedral[i][3] == atom_a){

            return this->dihedral_parm[i];
        }
        if(this->dihedral[i][0] == "X" && this->dihedral[i][1] == atom_b && this->dihedral[i][2] == atom_c && this->dihedral[i][3] == "X"){
            return this->dihedral_parm[i];
        }
        if(this->dihedral[i][0] == "X" && this->dihedral[i][1] == atom_c && this->dihedral[i][2] == atom_b && this->dihedral[i][3] == "X"){

            return this->dihedral_parm[i];
        }

   }
    cout << "Missing propers dihedral parameters for :"  << atom_b << " " << atom_a << " " << atom_c << " " << atom_d <<  endl;
    exit(0);


}
vector < double> gaffparser::GetParm_impropers( string atom_a, string atom_b, string atom_c, string atom_d){



    for(int i =0; i < this->improper.size(); i++){


          if(this->improper[i][1] == atom_b && this->improper[i][2] == atom_c){

              if(this->improper[i][0] == atom_a && this->improper[i][3] == atom_d){
                  return this->improper_parm[i];

              }
              if(this->improper[i][0] == atom_d && this->improper[i][3] == atom_a){
                  return this->improper_parm[i];

              }

              if(this->improper[i][0] == "X" && this->improper[i][3] == atom_a){
                  return this->improper_parm[i];

              }
              if(this->improper[i][0] == "X" && this->improper[i][3] == atom_d){
                  return this->improper_parm[i];

              }

          }

          if(this->improper[i][1] == atom_c && this->improper[i][2] == atom_b){

              if(this->improper[i][0] == atom_a && this->improper[i][3] == atom_d){
                  return this->improper_parm[i];

              }
              if(this->improper[i][0] == atom_d && this->improper[i][3] == atom_a){
                  return this->improper_parm[i];

              }
              if(this->improper[i][0] == "X" && this->improper[i][3] == atom_a){
                  return this->improper_parm[i];

              }
              if(this->improper[i][0] == "X" && this->improper[i][3] == atom_d){
                  return this->improper_parm[i];

              }

          }

          if(this->improper[i][1] == "X" && this->improper[i][2] == atom_c){

              if(this->improper[i][3] == atom_d){

                  return this->improper_parm[i];

              }

              if(this->improper[i][3] == atom_a){

                  return this->improper_parm[i];

              }
          }

          if(this->improper[i][1] == "X" && this->improper[i][2] == atom_b){

              if(this->improper[i][3] == atom_a){
                  return this->improper_parm[i];

              }
              if(this->improper[i][3] == atom_d){
                  return this->improper_parm[i];

              }

          }

      }
      vector < double> nulo;
      nulo.push_back(9999);
      return nulo;


}


vector < double> gaffparser::GetParm_atom(string atom_a){
    for(int i =0; i < this->atoms_id.size(); i++){
        if(atoms_id[i] == atom_a){
            return this->atom_parm[i];
        }

    }
}

double gaffparser::GetParm_masses(string atom_a){
    for(int i =0; i < this->atommass.size(); i++){
        if(atommass[i] == atom_a){
            return this->masses[i];
        }

    }
}

vector < int> gaffparser::reorg_impropers( string atom_a, string atom_b, string atom_c, string atom_d, int index_a, int index_b,int index_c,int index_d){

    vector< int> organized;

    for(int i =0; i < this->improper.size(); i++){


          if(this->improper[i][1] == atom_b && this->improper[i][2] == atom_c){

              if(this->improper[i][0] == atom_a && this->improper[i][3] == atom_d){
                  organized.push_back(index_a);
                  organized.push_back(index_b);
                  organized.push_back(index_c);
                  organized.push_back(index_d);

                  return organized;

              }
              if(this->improper[i][0] == atom_d && this->improper[i][3] == atom_a){
                  organized.push_back(index_d);
                  organized.push_back(index_b);
                  organized.push_back(index_c);
                  organized.push_back(index_a);

                  return organized;


              }

              if(this->improper[i][0] == "X" && this->improper[i][3] == atom_a){
                  organized.push_back(index_d);
                  organized.push_back(index_b);
                  organized.push_back(index_c);
                  organized.push_back(index_a);

                  return organized;


              }
              if(this->improper[i][0] == "X" && this->improper[i][3] == atom_d){
                  organized.push_back(index_a);
                  organized.push_back(index_b);
                  organized.push_back(index_c);
                  organized.push_back(index_d);

                  return organized;


              }

          }

          if(this->improper[i][1] == atom_c && this->improper[i][2] == atom_b){

              if(this->improper[i][0] == atom_a && this->improper[i][3] == atom_d){
                  organized.push_back(index_a);
                  organized.push_back(index_c);
                  organized.push_back(index_b);
                  organized.push_back(index_d);

                  return organized;


              }
              if(this->improper[i][0] == atom_d && this->improper[i][3] == atom_a){
                  organized.push_back(index_d);
                  organized.push_back(index_c);
                  organized.push_back(index_b);
                  organized.push_back(index_a);

                  return organized;
              }
              if(this->improper[i][0] == "X" && this->improper[i][3] == atom_a){
                  organized.push_back(index_d);
                  organized.push_back(index_c);
                  organized.push_back(index_b);
                  organized.push_back(index_a);

                  return organized;

              }
              if(this->improper[i][0] == "X" && this->improper[i][3] == atom_d){
                  organized.push_back(index_a);
                  organized.push_back(index_c);
                  organized.push_back(index_b);
                  organized.push_back(index_d);

                  return organized;


              }

          }

          if(this->improper[i][1] == "X" && this->improper[i][2] == atom_c){

              if(this->improper[i][3] == atom_d){

                  organized.push_back(index_a);
                  organized.push_back(index_b);
                  organized.push_back(index_c);
                  organized.push_back(index_d);

                  return organized;

              }

              if(this->improper[i][3] == atom_a){

                  organized.push_back(index_d);
                  organized.push_back(index_b);
                  organized.push_back(index_c);
                  organized.push_back(index_a);

                  return organized;

              }
          }

          if(this->improper[i][1] == "X" && this->improper[i][2] == atom_b){

              if(this->improper[i][3] == atom_a){
                  organized.push_back(index_d);
                  organized.push_back(index_c);
                  organized.push_back(index_b);
                  organized.push_back(index_a);

                  return organized;

              }
              if(this->improper[i][3] == atom_d){
                  organized.push_back(index_a);
                  organized.push_back(index_c);
                  organized.push_back(index_b);
                  organized.push_back(index_d);

                  return organized;

              }

          }

      }

}
