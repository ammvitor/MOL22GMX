#include <iostream>
#include<vector>
#include<iostream>
#include<fstream>
#include<vector>
#include<stdio.h>
#include<stdlib.h>
#include<string>
#include<time.h>
#include <tr1/memory>
#include <memory>
#include <iostream>
#include <iomanip>
#include<algorithm>
#include<gaffparser.h>
#include <stdio.h>
#include <string>
#include<mol2.h>
#include <getopt.h>

//-i /home/jvscunha/Models/Lyz_ligands/3HTB_lig.mol2 -x /home/jvscunha/Models/Lyz_ligands/3htb_lig.itp -g /home/jvscunha/Models/Lyz_ligands/3htb.gro -p /home/jvscunha/Models/Lyz_ligands/3HTB_rec.gro -c /home/jvscunha/Models/Lyz_ligands/3htb_complex.gro -t /home/jvscunha/Models/Lyz_ligands/topol.top
//-i /home/jvscunha/Models/ag_ligands/FAD.mol2 -x /home/jvscunha/Models/ag_ligands/fad_lig.itp -g /home/jvscunha/Models/ag_ligands/FAD.gro -p /home/jvscunha/Models/ag_ligands/FAD_rec.gro -c /home/jvscunha/Models/ag_ligands/FAD_complex.gro -t /home/jvscunha/Models/Lyz_ligands/topol.top
using namespace std;

//*reorganize the topology file adding the atomtype//
void toppy_top_rewriter(string topolfile, string ligand_atomtypes_itp, string ligand_itp, string mol_name){

    string include_atomtypes = "#include \"" + ligand_atomtypes_itp + "\" ";
    string include_parmitp = "#include \"" + ligand_itp + "\" ";

        ifstream myfile;
        myfile.open(topolfile.c_str());
        string line ="";
        int count = 0;
        string out_name = topolfile + ".new";
        FILE* out_lig = fopen(out_name.c_str(),"w");

        while (myfile)

        {
            getline (myfile,line);


            if(line[0] != ';'){
                if(count == 3 ){
                    fprintf(out_lig,"%s \n", include_atomtypes.c_str());
                    count++;
                }
                if(count == 4 ){
                    fprintf(out_lig,"%s \n", include_parmitp.c_str());
                    count++;
                }
                fprintf(out_lig,"%s \n", line.c_str());
                    count++;


            }

            else{
               fprintf(out_lig,"%s \n", line.c_str());
              //printf("%s \n", line.c_str());
            }


        }
            fprintf(out_lig,"%s              1 \n", mol_name.c_str());


}

bool equal_vector(std::vector<int> a, std::vector<int> b) {
  std::sort(a.begin(), a.end());
  std::sort(b.begin(), b.end());

  return std::equal(a.begin(), a.end(), b.begin());
}

static bool IsEqual(vector<int> a, vector<int> b)
{
    sort(a.begin(), a.end());
    sort(b.begin(), b.end());
    return (a == b);
}
static bool IsEqual(vector<string> a, vector<string> b)
{
    sort(a.begin(), a.end());
    sort(b.begin(), b.end());
    return (a == b);
}

/*
"-i /home/jvscunha/Models/Lyz_ligands/3HTB_lig.mol2 "
"-x /home/jvscunha/Models/Lyz_ligands/3htb_lig.itp "
"-g /home/jvscunha/Models/Lyz_ligands/3htb.gro "
"-p /home/jvscunha/Models/Lyz_ligands/3HTB_REC2.gro "
"-c /home/jvscunha/Models/Lyz_ligands/3htb_complex.gro"
*/

int main(int argc, char **argv)
{

    cout << " local " << argv[0] << endl;


  if (argc < 2) {

      cout << "Usage %s -i <inputfile_lig_file> \n "
              "-p <input_receptor_pdb -optional->  \n"
              "-c <output_complex_gro - only if -p -> \n "
              "-g <output_ligand_gro> \n "
              "-x <output_ligand_itp> [-h]\n" << endl;

    return 1;
  }
   int c;
   string  inputfile_lig_file, input_receptor_gro, input_topol, output_complex_gro, output_ligand_gro, output_ligand_itp,gaff_local;
   input_receptor_gro = "";


  while ((c = getopt(argc, argv, "i:p:c:g:t:h:x:")) != -1)



          switch (c){
          case 'i':
              inputfile_lig_file = string(optarg);
              break;
          case 'h':
              printf("Usage %s -i <inputfile_lig_file> \n -p <input_receptor_pdb -optional->  \n -c <output_complex_gro - only if -p -> \n -g <output_ligand_gro> \n -x <output_ligand_itp> [-h]\n");
              break;
              exit(1);
          case 'g':
              output_ligand_gro = string(optarg);
              break;
          case 'p':
              input_receptor_gro = string(optarg);
              break;
          case 'c':
              output_complex_gro = string(optarg);
              break;
          case 't':
              input_topol = string(optarg);
              break;
          case 'x':
              output_ligand_itp = string(optarg);
              break;

         }

  string atomtypes_out = output_ligand_itp.substr(0, output_ligand_itp.size()-4)+ "_atomtypes.itp";



  cout << "******************************************************************" << endl;
  cout << "                         MOL22GMX                                " << endl;
  cout << "               A GAFF parametrizer for GROMACS                   " << endl;
  cout << "                    Newcastle University                         " << endl;
  cout << "                 Cunha, JVS and Bronowska, AK                    " << endl;
  cout << "******************************************************************" << endl;



      // print the topo

   // Read the file
  cout << "* Reading Ligand .mol2 : " << inputfile_lig_file <<endl;
  cout << "******************************************************************" << endl;
   Mol2* mol_mol2 = new Mol2(inputfile_lig_file);

    cout << "* Molecule Name : " << mol_mol2->resnames[0] << endl;
    cout << "* Total Atoms : " << mol_mol2->N << endl;

    toppy_top_rewriter(input_topol,atomtypes_out,output_ligand_itp,mol_mol2->resnames[0]);

    std::string temp(argv[1]);
       // print out the molecule.
       cout << "* Printing Ligand gro : " << output_ligand_gro << endl;
       cout << "******************************************************************" << endl;

       vector< string> atomtypes;
       vector< vector < string> > grotypes_rec;
       vector< vector < string> > grotypes_lig;

       vector < string> gro_temp;

         string gro_name, gro_type;
         string gro_res_num, total_atom;
         string gro_x, gro_y, gro_z;
         string box_vec_1, box_vec_2, box_vec_3, box_vec_4,box_vec_5,box_vec_6,box_vec_7,box_vec_8,box_vec_9;
         int count =0;



         FILE* out_lig = fopen(output_ligand_gro.c_str(),"w");
         fprintf(out_lig,"Ligand\n");
         fprintf(out_lig," %i\n", mol_mol2->N);

         for(int i =0; i < mol_mol2->N; i++){



         fprintf(out_lig,"%5i%3s%7s%5i%8.3f%8.3f%8.3f\n",
                 1,
                 mol_mol2->resnames[0].c_str(),
                 mol_mol2->atomnames[i].c_str(),
                 i+1,
                 mol_mol2->xyz[i][0]/10,
                 mol_mol2->xyz[i][1]/10,
                 mol_mol2->xyz[i][2]/10);
         }

         fprintf(out_lig,"   0.00000   0.00000   0.00000");
         fclose(out_lig);
        string _lol;


       if(input_receptor_gro != ""){
        count =0;
       ifstream myfile2;
        myfile2.open(input_receptor_gro.c_str());
        cout << "* Reading Receptor gro : " << input_receptor_gro << endl;
        cout << "******************************************************************" << endl;


         if (myfile2.is_open())
         {
             getline(myfile2,_lol);
             getline(myfile2,total_atom);

             //myfile2 >> total_atom;
           while (!myfile2.eof())
           {
             if(count < atoi(total_atom.c_str())){
             myfile2 >> gro_name >> gro_type >> gro_res_num >> gro_x >> gro_y >> gro_z;
             gro_temp.push_back(gro_name);
             gro_temp.push_back(gro_type);
             gro_temp.push_back(gro_res_num);
             gro_temp.push_back(gro_x);
             gro_temp.push_back(gro_y);
             gro_temp.push_back(gro_z);
             grotypes_rec.push_back(gro_temp);
             gro_temp.clear();
             count++;
             }
             else{
                 myfile2 >> box_vec_1 >> box_vec_2 >> box_vec_3 >> box_vec_4 >> box_vec_5 >> box_vec_6 >> box_vec_7 >> box_vec_8 >> box_vec_9;
             }
           }
           myfile2.close();
         }
         else {
             cout << "Unable to open receptor .gro file";
               exit(1);
         }

         cout << "* Printing Complex .gro : " << output_complex_gro << endl;
         cout << "******************************************************************" << endl;

         string final_resnun= gro_name.substr(0, gro_name.size()-3);

         string buff = mol_mol2->molname;
         FILE* out_complex = fopen(output_complex_gro.c_str(),"w");
         fprintf(out_complex,"COMPLEX\n");
         fprintf(out_complex," %i\n", grotypes_rec.size()+mol_mol2->N);
         for(int i =0; i < grotypes_rec.size(); i++){
         fprintf(out_complex,"%8s%7s%5s%8s%8s%8s\n",grotypes_rec[i][0].c_str(),grotypes_rec[i][1].c_str(),grotypes_rec[i][2].c_str(),grotypes_rec[i][3].c_str(),grotypes_rec[i][4].c_str(),grotypes_rec[i][5].c_str());
         }
         for(int i =0; i < mol_mol2->N; i++){



             fprintf(out_complex,"%5i%3s%7s%5i%8.3f%8.3f%8.3f\n",
                     atoi(final_resnun.c_str())+1,
                     mol_mol2->resnames[0].c_str(),
                     mol_mol2->atomnames[i].c_str(),
                     i+1+grotypes_rec.size(),
                     mol_mol2->xyz[i][0]/10,
                     mol_mol2->xyz[i][1]/10,
                     mol_mol2->xyz[i][2]/10);
             }

         fprintf(out_complex,"%10s%10s%10s%10s%10s%10s%10s%10s%10s",box_vec_1.c_str(),box_vec_2.c_str(),box_vec_3.c_str(),box_vec_4.c_str(),box_vec_5.c_str(),box_vec_6.c_str(),box_vec_7.c_str(),box_vec_8.c_str(),box_vec_9.c_str());
         fclose(out_complex);

    }



   vector < vector < string> > atom_id_bonds;
   vector < vector < string> > id_double_bonds;
   vector < vector < string> > id_aromatic_bonds;
   vector < vector < string> > id_amine_bonds;

   vector < string> vec_temp;

   FILE* out_itp = fopen(output_ligand_itp.c_str(),"w");

   FILE* out_itp_atomtypes = fopen(atomtypes_out.c_str(),"w");


   for(int i = 0; i < mol_mol2->bonds.size(); i++){
       if(mol_mol2->bonds[i][2] == "ar"){
           vec_temp.push_back((mol_mol2->bonds[i][0].c_str()));
           vec_temp.push_back((mol_mol2->bonds[i][1].c_str()));

           id_aromatic_bonds.push_back(vec_temp);

       }if(mol_mol2->bonds[i][2] == "2"){
           vec_temp.push_back((mol_mol2->bonds[i][0].c_str()));
           vec_temp.push_back((mol_mol2->bonds[i][1].c_str()));
           id_double_bonds.push_back(vec_temp);

       }if(mol_mol2->bonds[i][2] == "am"){
           vec_temp.push_back((mol_mol2->bonds[i][0].c_str()));
           vec_temp.push_back((mol_mol2->bonds[i][1].c_str()));
           id_amine_bonds.push_back(vec_temp);

       }
       vec_temp.clear();

   }





   vector< string> atom_id_atoms = mol_mol2->gaffatoms;

   // unique atom paramters.


   sort( atom_id_atoms.begin(), atom_id_atoms.end() );
   atom_id_atoms.erase(unique( atom_id_atoms.begin(), atom_id_atoms.end() ), atom_id_atoms.end() );





   // angle finder.
   vector <vector < string> > atom_id_angles;
   vector < string> ang_Atom_triad;
     //Angles Finder
     for(int i =0; i < mol_mol2->bonds.size(); i++){
         for(int j = 0 ; j < mol_mol2->bonds.size(); j++){
             if(i!=j){

                      if(mol_mol2->bonds[i][0] == mol_mol2->bonds[j][0]){
                          ang_Atom_triad.push_back(mol_mol2->bonds[i][1]);
                          ang_Atom_triad.push_back(mol_mol2->bonds[i][0]);
                          ang_Atom_triad.push_back(mol_mol2->bonds[j][1]);
                          atom_id_angles.push_back(ang_Atom_triad);
                          ang_Atom_triad.clear();

                      }
                      if(mol_mol2->bonds[i][0] == mol_mol2->bonds[j][1]){
                          ang_Atom_triad.push_back(mol_mol2->bonds[j][0]);
                          ang_Atom_triad.push_back(mol_mol2->bonds[i][0]);
                          ang_Atom_triad.push_back(mol_mol2->bonds[i][1]);
                          atom_id_angles.push_back(ang_Atom_triad);
                          ang_Atom_triad.clear();

                      }
                      if(mol_mol2->bonds[i][1] == mol_mol2->bonds[j][0]){
                          ang_Atom_triad.push_back(mol_mol2->bonds[j][1]);
                          ang_Atom_triad.push_back(mol_mol2->bonds[j][0]);
                          ang_Atom_triad.push_back(mol_mol2->bonds[i][0]);
                          atom_id_angles.push_back(ang_Atom_triad);
                          ang_Atom_triad.clear();

                      }
                      if(mol_mol2->bonds[i][1] == mol_mol2->bonds[j][1]){
                          ang_Atom_triad.push_back(mol_mol2->bonds[i][0]);
                          ang_Atom_triad.push_back(mol_mol2->bonds[i][1]);
                          ang_Atom_triad.push_back(mol_mol2->bonds[j][0]);
                          atom_id_angles.push_back(ang_Atom_triad);
                          ang_Atom_triad.clear();

                      }


             }

         }
     }

     vector <vector < string> > atom_id_angles_2 = atom_id_angles;
     vector <vector < string> > atom_id_angles_3;


     //Angles Degenerescency killer
     while(atom_id_angles.size() > 0){

         for(int k = 1; k < atom_id_angles_2.size(); k++){


            if(IsEqual(atom_id_angles[0],atom_id_angles_2[k])){

            atom_id_angles_3.push_back(atom_id_angles[0]);
            atom_id_angles.erase(atom_id_angles.begin());
            atom_id_angles_2.erase(atom_id_angles_2.begin());

            atom_id_angles_2.erase(atom_id_angles_2.begin()+k-1);
            atom_id_angles.erase(atom_id_angles.begin()+k-1);
            break;

        }

      }

    }


     vector <vector < string> > atom_id_dihedral;
     vector < string> ang_Atom_quad;



     //Dihedral Finder
     for(int i =0; i < mol_mol2->bonds.size(); i++){
         for(int j = 0 ; j < mol_mol2->bonds.size(); j++){
             for(int k = 0; k < mol_mol2->bonds.size(); k++){


             if(i!=k && k!=j && j!=i){

                      if(mol_mol2->bonds[i][0] == mol_mol2->bonds[j][0]){

                          if(mol_mol2->bonds[i][1] == mol_mol2->bonds[k][0]){
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[j][1]);
                              atom_id_dihedral.push_back(ang_Atom_quad);
                              ang_Atom_quad.clear();

                          }
                          if(mol_mol2->bonds[i][1] == mol_mol2->bonds[k][1]){
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[j][1]);
                              atom_id_dihedral.push_back(ang_Atom_quad);
                              ang_Atom_quad.clear();

                          }
                          if(mol_mol2->bonds[j][1] == mol_mol2->bonds[k][0]){
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[j][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][1]);
                              atom_id_dihedral.push_back(ang_Atom_quad);
                              ang_Atom_quad.clear();

                          }
                          if(mol_mol2->bonds[j][1] == mol_mol2->bonds[k][1]){
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[j][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][1]);
                              atom_id_dihedral.push_back(ang_Atom_quad);
                              ang_Atom_quad.clear();

                          }


                      }


                      if(mol_mol2->bonds[i][0] == mol_mol2->bonds[j][1]){

                          if(mol_mol2->bonds[i][1] == mol_mol2->bonds[k][0]){
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[j][0]);
                              atom_id_dihedral.push_back(ang_Atom_quad);
                              ang_Atom_quad.clear();

                          }
                          if(mol_mol2->bonds[i][1] == mol_mol2->bonds[k][1]){
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[j][0]);
                              atom_id_dihedral.push_back(ang_Atom_quad);
                              ang_Atom_quad.clear();

                          }
                          if(mol_mol2->bonds[j][0] == mol_mol2->bonds[k][0]){
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[j][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][1]);
                              atom_id_dihedral.push_back(ang_Atom_quad);
                              ang_Atom_quad.clear();

                          }
                          if(mol_mol2->bonds[j][0] == mol_mol2->bonds[k][1]){
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[j][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][1]);
                              atom_id_dihedral.push_back(ang_Atom_quad);
                              ang_Atom_quad.clear();

                          }

                      }


                      if(mol_mol2->bonds[i][1] == mol_mol2->bonds[j][0]){

                          if(mol_mol2->bonds[i][0] == mol_mol2->bonds[k][0]){

                              ang_Atom_quad.push_back(mol_mol2->bonds[k][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[j][1]);
                              atom_id_dihedral.push_back(ang_Atom_quad);
                              ang_Atom_quad.clear();

                          }
                          if(mol_mol2->bonds[i][0] == mol_mol2->bonds[k][1]){
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[j][1]);
                              atom_id_dihedral.push_back(ang_Atom_quad);
                              ang_Atom_quad.clear();

                          }
                          if(mol_mol2->bonds[j][1] == mol_mol2->bonds[k][0]){

                              ang_Atom_quad.push_back(mol_mol2->bonds[k][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[j][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][0]);
                              atom_id_dihedral.push_back(ang_Atom_quad);
                              ang_Atom_quad.clear();

                          }
                          if(mol_mol2->bonds[j][1] == mol_mol2->bonds[k][1]){

                              ang_Atom_quad.push_back(mol_mol2->bonds[k][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[j][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][0]);
                              atom_id_dihedral.push_back(ang_Atom_quad);
                              ang_Atom_quad.clear();

                          }


                      }


                      if(mol_mol2->bonds[i][1] == mol_mol2->bonds[j][1]){

                          if(mol_mol2->bonds[i][0] == mol_mol2->bonds[k][0]){

                              ang_Atom_quad.push_back(mol_mol2->bonds[k][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[j][0]);
                              atom_id_dihedral.push_back(ang_Atom_quad);
                              ang_Atom_quad.clear();

                          }
                          if(mol_mol2->bonds[i][0] == mol_mol2->bonds[k][1]){

                              ang_Atom_quad.push_back(mol_mol2->bonds[k][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[j][0]);
                              atom_id_dihedral.push_back(ang_Atom_quad);
                              ang_Atom_quad.clear();

                          }
                          if(mol_mol2->bonds[j][0] == mol_mol2->bonds[k][0]){

                              ang_Atom_quad.push_back(mol_mol2->bonds[k][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[j][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][0]);
                              atom_id_dihedral.push_back(ang_Atom_quad);
                              ang_Atom_quad.clear();

                          }
                          if(mol_mol2->bonds[j][0] == mol_mol2->bonds[k][1]){

                              ang_Atom_quad.push_back(mol_mol2->bonds[k][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[j][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][0]);
                              atom_id_dihedral.push_back(ang_Atom_quad);
                              ang_Atom_quad.clear();
                      }
             }



             }

         }
       }
     }

     vector <vector < string> > atom_id_dihedral_2 = atom_id_dihedral;
     vector <vector < string> > atom_id_dihedral_3;
     vector <int> del_index;

     //Dihedral Degenerescency killer


     while(atom_id_dihedral.size() > 0){

         for(int k = 1; k < atom_id_dihedral.size(); k++){


             if(IsEqual(atom_id_dihedral_2[0],atom_id_dihedral[k])){
              del_index.push_back(k);

             }

        }

          if(del_index.size() > 0){
              atom_id_dihedral_3.push_back(atom_id_dihedral[0]);
              atom_id_dihedral_2.erase(atom_id_dihedral_2.begin());
              atom_id_dihedral.erase(atom_id_dihedral.begin());
          for(int i =0; i < del_index.size(); i++){
             atom_id_dihedral.erase(atom_id_dihedral.begin()+del_index[i]-i-1);
              atom_id_dihedral_2.erase(atom_id_dihedral_2.begin()+del_index[i]-i-1);
         }
          }
          else{
              atom_id_dihedral_3.push_back(atom_id_dihedral[0]);
              atom_id_dihedral_2.erase(atom_id_dihedral_2.begin());
              atom_id_dihedral.erase(atom_id_dihedral.begin());
          }

          del_index.clear();


      }




     //atom_id_dihedral.clear();
     //ang_Atom_quad.clear();

     vector <vector < string> > atom_id_pairs;
     vector < string> pair_temp;

     //pairs finder

     for(int i =0; i < atom_id_dihedral_3.size(); i++){
             pair_temp.push_back(atom_id_dihedral_3[i][0]);
             pair_temp.push_back(atom_id_dihedral_3[i][3]);
             atom_id_pairs.push_back(pair_temp);
             pair_temp.clear();
        }

     vector <vector < string> > atom_id_pairs_2 = atom_id_pairs;
     vector <vector < string> > atom_id_pairs_3;


     //pair Degenerescency killer

     std::sort(atom_id_pairs_2.begin(), atom_id_pairs_2.end());
     atom_id_pairs_2.erase(std::unique(atom_id_pairs_2.begin(), atom_id_pairs_2.end()), atom_id_pairs_2.end());



     vector < vector < string> > atom_id_impropers;

     //improper dihedral finder

     for(int i =0; i < mol_mol2->bonds.size(); i++){
         for(int j = 0 ; j < mol_mol2->bonds.size(); j++){
             for(int k = 0; k < mol_mol2->bonds.size(); k++){

             if(i!=k && k!=j && j!=i){
              if(mol_mol2->bonds[i][2] == "am" || mol_mol2->bonds[i][2] == "ar" || mol_mol2->bonds[i][2] == "2" ){

                      if(mol_mol2->bonds[i][0] == mol_mol2->bonds[j][0]){

                          if(mol_mol2->bonds[i][0] == mol_mol2->bonds[k][0]){

                              ang_Atom_quad.push_back(mol_mol2->bonds[j][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][1]);
                              atom_id_impropers.push_back(ang_Atom_quad);
                              ang_Atom_quad.clear();

                          }


                          if(mol_mol2->bonds[i][0] == mol_mol2->bonds[k][1]){

                              ang_Atom_quad.push_back(mol_mol2->bonds[j][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][0]);
                              atom_id_impropers.push_back(ang_Atom_quad);
                              ang_Atom_quad.clear();

                          }


                      }


                      if(mol_mol2->bonds[i][0] == mol_mol2->bonds[j][1]){

                          if(mol_mol2->bonds[i][0] == mol_mol2->bonds[k][0]){

                              ang_Atom_quad.push_back(mol_mol2->bonds[j][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][1]);
                              atom_id_impropers.push_back(ang_Atom_quad);
                              ang_Atom_quad.clear();

                          }
                          if(mol_mol2->bonds[i][0] == mol_mol2->bonds[k][1]){

                              ang_Atom_quad.push_back(mol_mol2->bonds[j][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][0]);
                              atom_id_impropers.push_back(ang_Atom_quad);
                              ang_Atom_quad.clear();

                          }


                      }


                      if(mol_mol2->bonds[i][1] == mol_mol2->bonds[j][0]){

                          if(mol_mol2->bonds[i][1] == mol_mol2->bonds[k][0]){

                              ang_Atom_quad.push_back(mol_mol2->bonds[j][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][1]);
                              atom_id_impropers.push_back(ang_Atom_quad);
                              ang_Atom_quad.clear();

                          }
                          if(mol_mol2->bonds[i][1] == mol_mol2->bonds[k][1]){

                              ang_Atom_quad.push_back(mol_mol2->bonds[j][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][0]);
                              atom_id_impropers.push_back(ang_Atom_quad);
                              ang_Atom_quad.clear();

                          }

                      }


                      if(mol_mol2->bonds[i][1] == mol_mol2->bonds[j][1]){

                          if(mol_mol2->bonds[i][1] == mol_mol2->bonds[k][0]){

                              ang_Atom_quad.push_back(mol_mol2->bonds[j][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][1]);
                              atom_id_impropers.push_back(ang_Atom_quad);
                              ang_Atom_quad.clear();

                          }
                          if(mol_mol2->bonds[i][1] == mol_mol2->bonds[k][1]){

                              ang_Atom_quad.push_back(mol_mol2->bonds[j][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][0]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[i][1]);
                              ang_Atom_quad.push_back(mol_mol2->bonds[k][0]);
                              atom_id_impropers.push_back(ang_Atom_quad);
                              ang_Atom_quad.clear();

                          }

                 }


                }
             }


         }

       }
     }


    string program = argv[0];
     string position = program.substr(0, program.size()-9);
     gaffparser *GP = new gaffparser(position+"/dat/gaff2.dat") ;
     vector < vector < double> >  atom_parm_final;
     vector < double >  atommass;

     vector < vector < double> >  bonds_parm_final;
     vector < vector < double> >  angles_parm_final;
     vector < vector < double> >  dihedral_parm_final;
     vector < vector < double> >  improper_parm_final;

     vector < double> temporary_vec;

     for(int i = 0; i < mol_mol2->N; i++){

        atommass.push_back(GP->GetParm_masses(mol_mol2->gaffatoms[i]));

    }



     //scaffold for atom parametrization
     for(int i = 0; i < atom_id_atoms.size(); i++){

        atom_parm_final.push_back(GP->GetParm_atom(atom_id_atoms[i]));

    }



     //scaffold for bonds parametrization
     for(int i = 0; i < mol_mol2->bonds.size(); i++){
          temporary_vec = GP->GetParm_bonds(mol_mol2->gaffatoms[(atoi(mol_mol2->bonds[i][0].c_str())-1)] ,
                  mol_mol2->gaffatoms[(atoi(mol_mol2->bonds[i][1].c_str())-1)]);
        bonds_parm_final.push_back(temporary_vec);

    }


    //scaffold for the angle parametrization

    for(int i = 0; i < atom_id_angles_3.size(); i++){
        temporary_vec = GP->GetParm_angles(mol_mol2->gaffatoms[atoi(atom_id_angles_3[i][0].c_str())-1] ,mol_mol2->gaffatoms[atoi(atom_id_angles_3[i][1].c_str())-1],mol_mol2->gaffatoms[atoi(atom_id_angles_3[i][2].c_str())-1]);
        angles_parm_final.push_back(temporary_vec);
    }


    //scaffold for the dihedral parametrization
    for(int i = 0; i < atom_id_dihedral_3.size(); i++){

        temporary_vec = GP->GetParm_dihedrals(mol_mol2->gaffatoms[atoi(atom_id_dihedral_3[i][0].c_str())-1] ,mol_mol2->gaffatoms[atoi(atom_id_dihedral_3[i][1].c_str())-1],mol_mol2->gaffatoms[atoi(atom_id_dihedral_3[i][2].c_str())-1],mol_mol2->gaffatoms[atoi(atom_id_dihedral_3[i][3].c_str())-1]);
        dihedral_parm_final.push_back(temporary_vec);
    }



    vector <vector < string> > atom_id_impropers_2;

    for(int i = 0; i < atom_id_impropers.size(); i++){



       temporary_vec = GP->GetParm_impropers(mol_mol2->gaffatoms[atoi(atom_id_impropers[i][0].c_str())-1], mol_mol2->gaffatoms[atoi(atom_id_impropers[i][1].c_str())-1],mol_mol2->gaffatoms[atoi(atom_id_impropers[i][2].c_str())-1], mol_mol2->gaffatoms[atoi(atom_id_impropers[i][3].c_str())-1]);


       if(temporary_vec[0]!= 9999){


            atom_id_impropers_2.push_back(atom_id_impropers[i]);

        }
    }


    temporary_vec.clear();
    atom_id_impropers.clear();






    vector <vector < string> > atom_id_impropers_4 = atom_id_impropers_2;
    vector <vector < string> > atom_id_impropers_3;


    while(atom_id_impropers_2.size() > 0){

        for(int k = 1; k < atom_id_impropers_4.size(); k++){


           if(IsEqual(atom_id_impropers_2[0],atom_id_impropers_4[k])){
            del_index.push_back(k);

           }

      }

        if(del_index.size() > 0){
            atom_id_impropers_3.push_back(atom_id_impropers_4[0]);
            atom_id_impropers_2.erase(atom_id_impropers_2.begin());
            atom_id_impropers_4.erase(atom_id_impropers_4.begin());
        for(int i =0; i < del_index.size(); i++){
           atom_id_impropers_4.erase(atom_id_impropers_4.begin()+del_index[i]-i-1);
            atom_id_impropers_2.erase(atom_id_impropers_2.begin()+del_index[i]-i-1);
       }
        }
        else{
            atom_id_impropers_3.push_back(atom_id_impropers_4[0]);
            atom_id_impropers_2.erase(atom_id_impropers_2.begin());
            atom_id_impropers_4.erase(atom_id_impropers_4.begin());
        }

        del_index.clear();

   }



    atom_id_impropers_4 = atom_id_impropers_3;
    vector <vector < int> > atom_id_impropers_final_int;

    vector < int> temporary_vec_int;
    for (int i =0; i < atom_id_impropers_3.size(); i++){
             temporary_vec_int =  GP->reorg_impropers(
                      mol_mol2->gaffatoms[atoi(atom_id_impropers_3[i][0].c_str())-1],
                      mol_mol2->gaffatoms[atoi(atom_id_impropers_3[i][1].c_str())-1],
                      mol_mol2->gaffatoms[atoi(atom_id_impropers_3[i][2].c_str())-1],
                      mol_mol2->gaffatoms[atoi(atom_id_impropers_3[i][3].c_str())-1],
                      atoi(atom_id_impropers_3[i][0].c_str()),
                      atoi(atom_id_impropers_3[i][1].c_str()),
                      atoi(atom_id_impropers_3[i][2].c_str()),
                      atoi(atom_id_impropers_3[i][3].c_str()));
    atom_id_impropers_final_int.push_back(temporary_vec_int);
    temporary_vec_int.clear();
    }

    for(int i = 0; i < atom_id_impropers_3.size(); i++){


        temporary_vec = GP->GetParm_impropers(mol_mol2->gaffatoms[atom_id_impropers_final_int[i][0]-1], mol_mol2->gaffatoms[atom_id_impropers_final_int[i][1]-1],mol_mol2->gaffatoms[atom_id_impropers_final_int[i][2]-1],mol_mol2->gaffatoms[atom_id_impropers_final_int[i][3]-1]);


        improper_parm_final.push_back(temporary_vec);

    }

    cout << "* Printing Ligand .itp : " << output_ligand_itp << endl;

    //test printer to itp
        fprintf(out_itp_atomtypes,"[ atomtypes ]\n");
        fprintf(out_itp_atomtypes,";name   bond_type     mass     charge   ptype   sigma         epsilon \n");
     for(int i = 0; i < atom_id_atoms.size(); i++){
        fprintf(out_itp_atomtypes," %-2s       %-2s          0.00000  0.00000   A%16.5e%14.5e\n", atom_id_atoms[i].c_str() , atom_id_atoms[i].c_str() , atom_parm_final[i][0]*0.0891*2 ,atom_parm_final[i][1]* 4.184);

     }
        fclose(out_itp_atomtypes);

     fprintf(out_itp,"[ moleculetype ]\n");
     fprintf(out_itp,"%s         3\n",mol_mol2->resnames[0].c_str());

     fprintf(out_itp,"\n");

        fprintf(out_itp,"[ atoms ]\n");
       fprintf(out_itp,";   nr  type  resnr resid  atom  cgnr   charge     mass\n");


            cout << "* Printing atoms in the .itp file " << endl;

     for(int i = 0; i < mol_mol2->N; i++)
     {
         fprintf(out_itp,"%6i%5s%6i%6s%6s%5i%13.6lf%14.6lf\n", i+1 , mol_mol2->gaffatoms[i].c_str() , 1 , mol_mol2->molname.c_str() , mol_mol2->atomnames[i].c_str() , i+1 , mol_mol2->charges[i] ,atommass[i]);

     }

     cout << "* Printing bonds in the .itp file " << endl;

     fprintf(out_itp,"\n");

    fprintf(out_itp,"[ bonds ]\n");

     for(int i = 0; i < mol_mol2->bonds.size(); i++)
     {

          fprintf(out_itp,"%6s%7s%4i%14.4e%14.4e\n", mol_mol2->bonds[i][0].c_str() ,
                  mol_mol2->bonds[i][1].c_str() , 1 ,
                  bonds_parm_final[i][1]*0.1 , bonds_parm_final[i][0]*200*4.186);

     }
     cout << "* Printing Pairs in the .itp file " << endl;


     fprintf(out_itp,"\n");

    fprintf(out_itp,"[ pairs ]\n");
     for(int i = 0; i < atom_id_pairs_2.size(); i++)
     {
        fprintf(out_itp,"%6s%7s%7i\n", atom_id_pairs_2[i][0].c_str() , atom_id_pairs_2[i][1].c_str(), 1);
     }
     cout << "* Printing Angles in the .itp file " << endl;


     fprintf(out_itp,"\n");

     fprintf(out_itp,"[ angles ]\n");

     for(int i = 0; i < atom_id_angles_3.size(); i++)
     {
       fprintf(out_itp,"%6s%7s%7s%7i%14.4e%14.4e\n", atom_id_angles_3[i][0].c_str() , atom_id_angles_3[i][1].c_str() , atom_id_angles_3[i][2].c_str() , 1 , angles_parm_final[i][1] , 2*4.186*angles_parm_final[i][0]);

     }
     cout << "* Printing Proper Dihedrals in the .itp file " << endl;


     fprintf(out_itp,"\n");

    fprintf(out_itp,"[ dihedrals ]\n");
     for(int i = 0; i < atom_id_dihedral_3.size(); i++)
     {
        fprintf(out_itp,"%6s%7s%7s%7s%7i%9.2f%10.5f%4.0f\n", atom_id_dihedral_3[i][0].c_str() , atom_id_dihedral_3[i][1].c_str() , atom_id_dihedral_3[i][2].c_str() , atom_id_dihedral_3[i][3].c_str() , 9 , dihedral_parm_final[i][2] , dihedral_parm_final[i][1]*1.046 , dihedral_parm_final[i][0]);

     }
     cout << "* Printing Impropers Dihedrals in the .itp file " << endl;
     cout << "******************************************************************" << endl;

     fprintf(out_itp,"\n");

     fprintf(out_itp,"[ dihedrals ];Impropers\n");
    for(int i = 0; i < atom_id_impropers_final_int.size(); i++)
    {

        fprintf(out_itp,"%6i%7i%7i%7i%7i%9.2f%10.5f%4.0f\n", atom_id_impropers_final_int[i][0] , atom_id_impropers_final_int[i][1] , atom_id_impropers_final_int[i][2] , atom_id_impropers_final_int[i][3] , 4 , improper_parm_final[i][1] , improper_parm_final[i][0]* 4.184 , improper_parm_final[i][2]);

    }
       fclose(out_itp);

       string include_atomtypes = "#include \"" + atomtypes_out + "\" ";
       string include_parmitp = "#include \"" + output_ligand_itp + "\" ";

       ifstream myfile;
       myfile.open(input_topol.c_str());
       string line ="";
       count = 0;
       string out_name_top = "/home/jvscunha/Models/Lyz_ligands/teste_newtop.top";
       FILE* out_top = fopen(out_name_top.c_str(),"w");

       while (myfile)

       {
           getline (myfile,line);

           if(line[0] != ';'){
               if(count == 3 ){
                   fprintf(out_top,"%s \n", include_atomtypes.c_str());
                   count++;
               }
               if(count == 4 ){
                   fprintf(out_top,"%s \n", include_parmitp.c_str());
                   count++;
               }
               fprintf(out_top,"%s \n", line.c_str());
                   count++;


           }
           else{
              fprintf(out_top,"%s \n", line.c_str());
             //printf("%s \n", line.c_str());
           }

       }

   }




//-i /home/jvscunha/Models/Lyz_ligands/3jzb2.pdb  -x /home/jvscunha/Models/Lyz_ligands/3htb.itp -g /home/jvscunha/Models/Lyz_ligands/3htb.gro -p /home/jvscunha/Models/Lyz_ligands/3htb_rec.gro -c /home/jvscunha/Models/Lyz_ligands/3htb_complex.gro
//-i /home/jvscunha/Models/ag_ligand/FAD.pdb  -x /home/jvscunha/Models/ag_ligand/FAD.itp -g /home/jvscunha/Models/ag_ligand/FAD.gro -p /home/jvscunha/Models/Lyz_ligands/3htb_rec.gro -c /home/jvscunha/Models/Lyz_ligands/3htb_complex.gro
