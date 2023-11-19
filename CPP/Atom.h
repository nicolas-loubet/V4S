#pragma once
#include "Vector.h"
#include <iostream>

/**
 * Version: August 2023
 */

/**
 * This class creates an Atom object, with a position in a 3D system and optionally a velocity
 */
class Atom : public Particle {
	private:
		Vector* vel;
		float charge;
		float mass;
		string type;

	public:
        //Setters and getters
		void setVelocity(Vector* velocity) { vel= velocity; }
		Vector* getVelocity() { return vel; }
		void setCharge(const float c) { charge= c; }
		float getCharge() { return charge; }
		void setMass(const float m) { mass= m; }
		float getMass() { return mass; }
		void setType(const string t) { type= t; }
		string getType() { return type; }

        /**
         * Basic constructor without velocity

         * @param position Is the *Vector of the position
         * @param mass Is the mass of the Atom
         * @param ID Is the ID of the Atom
         */
		Atom(Vector* position, const float mass= 0, const int ID= 0): Particle(ID, position), mass(mass), vel(NULL) {}

        /**
         * Basic constructor with velocity
         * @param position Is the *Vector of the position
         * @param velocity Is the *Vector of the velocity
         * @param mass Is the mass of the Atom
         */
		Atom(Vector* position, Vector* velocity, const float mass= 0, const int ID= 0): Particle(ID, position), mass(mass), vel(velocity) {}

        /**
         * Constructor for assigment
         * @param a Other Atom already defined
         */
        Atom(const Atom& a): Particle(a.ID, a.pos), mass(a.mass), charge(a.charge), vel(a.vel) {}

        /**
         * Destructor. It destroys the position and velocity if not already destroyed.
         */
		~Atom() {
			if(vel != NULL) {
                delete(vel);
                vel= NULL;
			}
		}

        /**
         * Return the momentum. Must be deleted after used
         * @return The vel vector * the mass of the atom
         */
		Vector* getMomentum() {
            return *vel*mass;
		}


        /**
         * Calculates the electrostatic potential energy between two atoms
         * @param a *Atom to the other atom
         * @param bounds The coordinate of the last point, so the components are the width, height and length
         * @return The potential energy in kJ/mol between the two atoms
         */
        float getCoulombPotential(Atom* a, Vector* bounds) {
            const float K=  1389.35458; //kJ/mol to e-2
            return (K*this->getCharge()*a->getCharge()) / dist(this->getPosition(), a->getPosition(), bounds);
        }


};
