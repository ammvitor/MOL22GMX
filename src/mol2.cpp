#include "mol2.h"

Mol2::Mol2(string molfile)
{

            FILE *mol2file;
            mol2file = fopen(molfile.c_str(), "r");

            int tint;
            float tx, ty, tz;
            vector<double> txyz;
            int tres;
            float tcharge;
            int count=0;
            char tatomtype[10];
            char resname[10];
            string cpstr;

            char str[80];
            if (mol2file !=NULL){
                str[0]='#';
                while(str[0] !='@'){
                    fgets(str, 80, mol2file);

                }

                fgets(str, 80, mol2file);
                this->molname = str;
                this->molname = this->molname.substr(0,this->molname.size()-5);
                fscanf(mol2file, "%d %d %d %d %d", &this->N, &this->Nbonds, &this->Nres, &tint, &tint);


                cpstr = string(str);
                while (cpstr.substr(0,13) != "@<TRIPOS>ATOM"){
                    fgets(str, 80, mol2file);
                    cpstr = string(str);
                }

                for (int i=0; i<this->N; i++){
                    fscanf(mol2file, "%d %s %f %f %f %5s%d %s %f\n", &tint, str, &tx, &ty, &tz, tatomtype, &tres, &resname, &tcharge);
                    txyz.push_back(tx);
                    txyz.push_back(ty);
                    txyz.push_back(tz);
                    this->xyz.push_back(txyz);
                    txyz.clear();

                    this->charges.push_back(tcharge);
                    this->atomnames.push_back(str);
                    this->resnames.push_back(resname);
                    this->gaffatoms.push_back(tatomtype);

                }

                fscanf(mol2file, "%s\n", str);
                if (str[0] != '@'){
                    while (str[0] != '@'){
                        fgets(str, 80, mol2file);
                    }
                }

                vector<string> bond;
                char s1[6], s2[6], s3[5];
                for (int i=0; i<this->Nbonds; i++){
                    fscanf(mol2file, "%d%s%s%s\n", &tint, s1, s2, s3);

                    bond.push_back(string(s1));
                    bond.push_back(string(s2));
                    bond.push_back(string(s3));
                    this->bonds.push_back(bond);
                    bond.clear();
                }
            }
            else {
                printf("Mol2 file could not be opened! Please check!\n");
                exit(1);
            }

            fclose(mol2file);
        }


