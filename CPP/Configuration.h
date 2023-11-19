#pragma once
#include "Water.h"
#include <iostream>

/**
 * Version: August 2023
 */

/**
 * This class creates a Configuration object, with an array of Water objects
 */
class Configuration {
	protected:
		Water** molecs; //The pointer to the pointer to the first element of the array
		int N_MOLEC= 0; //The number of Water objects in the array
		Vector* bounds; //The bounds of the system

		/**
		 * Transforms a string to a Vector that corresponds to Coordinates
		 * @param l Line that contains the coordinates
		 * @return A Vector to the coordinate read
		 */
		Vector* getCoordFromLine(const string l) {
			//   57SOL    HW1  226   0.175   2.143   2.375
			//0         1         2 ______  3_____  __4___
			float x= stof(l.substr(22,6))*10;
			float y= stof(l.substr(30,6))*10;
			float z= stof(l.substr(38,6))*10;
			return new Vector(x,y,z);
		}

		/**
		 * Transforms a string to a Vector that corresponds to Bounds (last line of the .gro file)
		 * @param l Line that contains the bounds
		 * @return A Vector that indicates the position of the last point in the box
		 */
		Vector* setBoundsFromLine(const string l) {
			float x= stof(l.substr(2,9))*10;
			float y= stof(l.substr(12,9))*10;
			float z= stof(l.substr(22,9))*10;

			return new Vector(x,y,z);
		}


	public:
        //Getters
		int getNMolec() { return N_MOLEC; }
		Water* getMolec(int id) { return molecs[id-1]; }
		Vector* getBounds() { return bounds; }

        /**
         * Basic constructor without velocity. You must use delete(*Configuration) after using it
         * @param file_dif Path to the .gro file or similar
         * @param water_model Use defined values WMODEL_
         */
		Configuration(const string file_dir, const int water_model) {
			string line;
			const float O_MASS= 15.99940;
			const float H_MASS= 1.00800;

			int n_atoms= 0;
			if(water_model == WMODEL_SPCE) n_atoms= 3;
			else if(water_model == WMODEL_TIP3P) n_atoms= 3;
			else if(water_model == WMODEL_TIP4P) n_atoms= 4;
			else if(water_model == WMODEL_TIP5P) n_atoms= 5;

			ifstream f(file_dir);
			getline(f, line);
			//cout << "Reading " << line << endl; //Title
			getline(f, line);
			N_MOLEC= stoi(line)/n_atoms;

			molecs= new Water*[N_MOLEC]; //It must be deleted because of the new; See Destructor

			Atom* o;
			Atom* h1;
			Atom* h2;
			Atom* l1;
			Atom* l2;

			for(int m= 1; m <= N_MOLEC; m++) {
				getline(f, line);
				o= new Atom(getCoordFromLine(line), O_MASS);
				getline(f, line);
				h1= new Atom(getCoordFromLine(line), H_MASS);
				getline(f, line);
				h2= new Atom(getCoordFromLine(line), H_MASS);

				if(n_atoms >= 4) {
					getline(f, line);
					l1= new Atom(getCoordFromLine(line));
				}
				if(n_atoms >= 5) {
					getline(f, line);
					l2= new Atom(getCoordFromLine(line));
				}

                if(water_model == WMODEL_SPCE) molecs[m-1]= new SPCEWater(m,o,h1,h2);
                else if(water_model == WMODEL_TIP3P) molecs[m-1]= new TIP3PWater(m,o,h1,h2);
                else if(water_model == WMODEL_TIP4P) molecs[m-1]= new TIP4PWater(m,o,h1,h2,l1);
                else if(water_model == WMODEL_TIP5P) molecs[m-1]= new TIP5PWater(m,o,h1,h2,l1,l2);
			}

			getline(f, line);
			bounds= setBoundsFromLine(line);
			f.close();
		}

        /**
         * Destructor. It destroys the molecule array and bounds
         */
		~Configuration() {
			for(int m= 0; m < N_MOLEC; m++) {
				delete(molecs[m]);
            }
			delete(molecs);
			delete(bounds);
		}

};
