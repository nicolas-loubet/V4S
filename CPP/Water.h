#pragma once
#include "Particle.h"
#include "Atom.h"

#define WMODEL_TIP3P 30
#define WMODEL_SPCE 3
#define WMODEL_TIP4P 4
#define WMODEL_TIP5P 5

#define D_MOLECULE 0
#define T0_MOLECULE 1
#define T1_MOLECULE 2
#define T2_MOLECULE 3

/**
 * Version: August 2023
 */

/**
 * This class creates a Water molecule object, with the 3 atoms: H2O
 */
class Water : public Particle {
	protected:
		Atom* oxygen;
		Atom* hydrogen_1;
		Atom* hydrogen_2;
		int classif= ERROR_VALUE; //It's used with the defines

	public:

        //Setters and getters
		void setAtom_O(Atom* a) {
			delete(oxygen);
			oxygen= a;
			pos= a->getPosition();
		}
		void setAtom_H1(Atom* a) { delete(hydrogen_1); hydrogen_1= a; }
		void setAtom_H2(Atom* a) { delete(hydrogen_2); hydrogen_2= a; }
		void setClassification(int c) { classif= c; }
		Atom* getAtom_O() { return oxygen; }
		Atom* getAtom_H1() { return hydrogen_1; }
		Atom* getAtom_H2() { return hydrogen_2; }
		int getClassification() { return classif; }

        /**
         * Basic constructor. You must use delete(*Water) after using it
         * @param id Number of ID to identify the molecule in the configuration
         * @param o *Atom for oxygen
         * @param h1 *Atom for one hydrogen
         * @param h2 *Atom for the other hydrogen
         */
		Water(int id, Atom* o, Atom* h1, Atom* h2): Particle(id, o->getPosition()), oxygen(o), hydrogen_1(h1), hydrogen_2(h2) {}

        /**
         * Constructor for assigment
         * @param w Other Water molecule already defined
         */
        Water(const Water& w): Water(w.ID, w.oxygen, w.hydrogen_1, w.hydrogen_2) {}

        /**
         * Destructor. It destroys the atoms and the position of the molecule
         * Virtual because each model redefines it but at the end calls implicitly this one
         */
		virtual ~Water() {
			if(*pos == *(oxygen->getPosition())) //If the Vectors are the same, it only deletes one of it avoinding free() exception
                pos= NULL;
			delete(oxygen);
			delete(hydrogen_1);
			delete(hydrogen_2);
		}

		/**
		 * Returns an array of *Atom of all the atoms of the Water molecule
		 * Virtual, must be redefined explicitly
		 * @return The array of *Atom
		 */
		virtual Atom** getArrAtoms() = 0;

		/**
		 * Calculates the potential of this molecule with another specified
		 * Virtual, must be redefined explicitly
		 * @param m *Water to the other molecule
		 * @param bounds The coordinate of the last point, so the components are the width, height and length
		 * @return the potential energy in kJ/mol of the interaction between this two molecules
		 */
		virtual float potentialWith(Water* m, Vector* bounds) = 0;

        /**
         * Calculates the Lennard-Jones potential energy between two atoms
         * @param m *Particle to the other molecule to calculate the LJPotential with
         * @param s Parameter of the Lennard-Jones potential ; s (Anstrong)
         * @param e Parameter of the Lennard-Jones potential ; e (kJ/mol)
         * @param bounds The coordinate of the last point, so the components are the width, height and length
         * @return The L-J potential energy in kJ/mol between the two molecules
         */
        float getLJPotential(Particle* m, const float s, const float e, Vector* bounds) {
            const float R= dist(this->pos, m->getPosition(), bounds);
            return 4*e*(pow(s/R,12)-pow(s/R,6));
        }

};

/**
 * This class creates a Water molecule object, with 3 atoms: H2O (SPC/E model https://doi.org/10.1021/j100308a038)
 */
class SPCEWater : public Water {
	public:
        /**
         * Basic constructor. You must use delete(*SPCEWater) after using it
         * @param id Number of ID to identify the molecule in the configuration
         * @param o *Atom for oxygen
         * @param h1 *Atom for one hydrogen
         * @param h2 *Atom for the other hydrogen
         */
		SPCEWater(int id, Atom* o, Atom* h1, Atom* h2): Water(id, o, h1, h2) {
            o->setCharge(-0.8476);
            h1->setCharge(0.4238);
            h2->setCharge(0.4238);
		}

        /**
         * Constructor for assigment
         * @param w Other SPCEWater molecule already defined
         */
        SPCEWater(const SPCEWater& w): SPCEWater(w.ID, w.oxygen, w.hydrogen_1, w.hydrogen_2) {}

        /**
         * Destructor. It destroys the atoms and the position of the molecule
         */
		~SPCEWater() {/*It only calls Water destructor*/}

		/**
		 * Returns an array of *Atom of all the atoms of the SPCEWater molecule
		 * The array must be deleted
		 * @return The array of *Atom
		 */
		Atom** getArrAtoms() override {
            Atom** output= new Atom*[3];
            output[0]= this->oxygen;
			output[1]= this->hydrogen_1;
			output[2]= this->hydrogen_2;
			return output;
		}

		/**
		 * Calculates the potential of this molecule with another specified
		 * @param m *SPCEWater to the other molecule
		 * @param bounds The coordinate of the last point, so the components are the width, height and length
		 * @return the potential energy in kJ/mol of the interaction between this two molecules
		 */
		float potentialWith(Water* m, Vector* bounds) override {
			SPCEWater* molec= (SPCEWater*) m;

			float Vtot= 0;

			Atom** arr_this= this->getArrAtoms();
			Atom** arr_other= m->getArrAtoms();

			for(int i= 0; i < 3; i++)
                for(int j= 0; j < 3; j++)
                    Vtot+= arr_this[i]->getCoulombPotential(arr_other[j], bounds);

            delete(arr_this);
            delete(arr_other);

			Vtot+= getLJPotential(m, 3.166, 0.650, bounds);

			return Vtot;
		}
};

/**
 * This class creates a Water molecule object, with 3 atoms: H2O (TIP3P https://doi.org/10.1063/1.445869 )
 */
class TIP3PWater : public SPCEWater {
	public:
        /**
         * Basic constructor. You must use delete(*SPCEWater) after using it
         * @param id Number of ID to identify the molecule in the configuration
         * @param o *Atom for oxygen
         * @param h1 *Atom for one hydrogen
         * @param h2 *Atom for the other hydrogen
         */
		TIP3PWater(int id, Atom* o, Atom* h1, Atom* h2): SPCEWater(id, o, h1, h2) {
            o->setCharge(-0.8340);
            h1->setCharge(0.4170);
            h2->setCharge(0.4170);
		}

		/**
		 * Calculates the potential of this molecule with another specified
		 * @param m *SPCEWater to the other molecule
		 * @param bounds The coordinate of the last point, so the components are the width, height and length
		 * @return the potential energy in kJ/mol of the interaction between this two molecules
		 */
		float potentialWith(Water* m, Vector* bounds) override {
			TIP3PWater* molec= (TIP3PWater*) m;

			float Vtot= 0;

			Atom** arr_this= this->getArrAtoms();
			Atom** arr_other= m->getArrAtoms();

			for(int i= 0; i < 3; i++)
                for(int j= 0; j < 3; j++)
                    Vtot+= arr_this[i]->getCoulombPotential(arr_other[j], bounds);

            delete(arr_this);
            delete(arr_other);

			Vtot+= getLJPotential(m, 3.15061, 0.6364, bounds);

			return Vtot;
		}
};

/**
 * This class creates a Water molecule object, with 4 atoms: H2O+M (TIP4P/2005 model https://doi.org/10.1063/1.2121687)
 */
class TIP4PWater : public Water {
	protected:
		Atom* electrons;
	public:
        //Setters and getters
		void setAtom_L1(Atom* a) { delete(electrons); electrons= a; }
		Atom* getAtom_L1() { return electrons; }

        /**
         * Basic constructor. You must use delete(*TIP4PWater) after using it
         * @param id Number of ID to identify the molecule in the configuration
         * @param o *Atom for oxygen
         * @param h1 *Atom for one hydrogen
         * @param h2 *Atom for the other hydrogen
         * @param l1 *Atom for one "electrons" atom
         */
		TIP4PWater(int id, Atom* o, Atom* h1, Atom* h2, Atom* l1): Water(id, o, h1, h2), electrons(l1) {
            o->setCharge(0);
            h1->setCharge(0.5564);
            h2->setCharge(0.5564);
            l1->setCharge(-1.1128);
		}

        /**
         * Constructor for assigment
         * @param w Other TIP4PWater molecule already defined
         */
        TIP4PWater(const TIP4PWater& w): TIP4PWater(w.ID, w.oxygen, w.hydrogen_1, w.hydrogen_2, w.electrons) {}

        /**
         * Destructor. It destroys the atoms and the position of the molecule
         */
		~TIP4PWater() {
			delete(electrons);
		}

		/**
		 * Returns an array of *Atom of all the atoms of the TIP4PWater molecule
		 * The array must be deleted
		 * @return The array of *Atom
		 */
		Atom** getArrAtoms() override {
            Atom** output= new Atom*[4];
            output[0]= this->oxygen;
			output[1]= this->hydrogen_1;
			output[2]= this->hydrogen_2;
			output[3]= this->electrons;
			return output;
		}

		/**
		 * Calculates the potential of this molecule with another specified
		 * @param m *TIP4PWater to the other molecule
		 * @param bounds The coordinate of the last point, so the components are the width, height and length
		 * @return the potential energy in kJ/mol of the interaction between this two molecules
		 */
		float potentialWith(Water* m, Vector* bounds) override {
			TIP4PWater* molec= (TIP4PWater*) m;

			float Vtot= 0.;

			Atom** arr_this= this->getArrAtoms();
			Atom** arr_other= m->getArrAtoms();

			for(int i= 0; i < 4; i++)
                for(int j= 0; j < 4; j++)
                    Vtot+= arr_this[i]->getCoulombPotential(arr_other[j], bounds);

            delete(arr_this);
            delete(arr_other);

			Vtot+= getLJPotential(m, 3.1589, 0.774907, bounds);

			return Vtot;
		}
};

/**
 * This class creates a Water molecule object, with 5 atoms: H2O+L2 (TIP5P model https://doi.org/10.1063/1.481505) [There is a 2018 update]
 */
class TIP5PWater : public Water {
	protected:
		Atom* electrons_1;
		Atom* electrons_2;
	public:
        //Setters and getters
		void setAtom_L1(Atom* a) { delete(electrons_1); electrons_1= a; }
		Atom* getAtom_L1() { return electrons_1; }
		void setAtom_L2(Atom* a) { delete(electrons_2); electrons_2= a; }
		Atom* getAtom_L2() { return electrons_2; }

        /**
         * Basic constructor. You must use delete(*TIP5PWater) after using it
         * @param id Number of ID to identify the molecule in the configuration
         * @param o *Atom for oxygen
         * @param h1 *Atom for one hydrogen
         * @param h2 *Atom for the other hydrogen
         * @param l1 *Atom for one "electrons" atom
         * @param l2 *Atom for the other "electrons" atom
         */
		TIP5PWater(int id, Atom* o, Atom* h1, Atom* h2, Atom* l1, Atom* l2): Water(id, o, h1, h2), electrons_1(l1), electrons_2(l2) {
            o->setCharge(0);
            h1->setCharge(0.241);
            h2->setCharge(0.241);
            l1->setCharge(-0.241);
            l2->setCharge(-0.241);
		}

        /**
         * Constructor for assigment
         * @param w Other TIP5PWater molecule already defined
         */
        TIP5PWater(const TIP5PWater& w): TIP5PWater(w.ID, w.oxygen, w.hydrogen_1, w.hydrogen_2, w.electrons_1, w.electrons_2) {}

        /**
         * Destructor. It destroys the atoms and the position of the molecule
         */
		~TIP5PWater() {
			delete(electrons_1);
			delete(electrons_2);
		}

		/**
		 * Returns an array of *Atom of all the atoms of the TIP5PWater molecule
		 * The array must be deleted
		 * @return The array of *Atom
		 */
		Atom** getArrAtoms() override {
            Atom** output= new Atom*[5];
            output[0]= this->oxygen;
			output[1]= this->hydrogen_1;
			output[2]= this->hydrogen_2;
			output[3]= this->electrons_1;
			output[4]= this->electrons_2;
			return output;
		}

		/**
		 * Calculates the potential of this molecule with another specified
		 * @param m *TIP5PWater to the other molecule
		 * @param bounds The coordinate of the last point, so the components are the width, height and length
		 * @return the potential energy in kJ/mol of the interaction between this two molecules
		 */
		float potentialWith(Water* m, Vector* bounds) override {
			TIP5PWater* molec= (TIP5PWater*) m;

			float Vtot= 0.;

			Atom** arr_this= this->getArrAtoms();
			Atom** arr_other= m->getArrAtoms();

			for(int i= 0; i < 5; i++)
                for(int j= 0; j < 5; j++)
                    Vtot+= arr_this[i]->getCoulombPotential(arr_other[j], bounds);

            delete(arr_this);
            delete(arr_other);

			Vtot+= getLJPotential(m, 3.0970, 0.74480, bounds);

			return Vtot;
		}
};
