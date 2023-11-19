#pragma once
#include "Util.h"

/**
 * Version: August 2023
 */

/**
 * This class creates a Particle object, with an ID to identify it and a position in a 3D system
 */
class Particle {
	protected:
		Vector* pos;
		int ID= -1;

	public:
        //Setters and getters
		void setPosition(Vector* position) { pos= position; }
		Vector* getPosition() { return pos; }
		void setID(int i) { ID= i; }
		int getID() { return ID; }

        /**
         * Basic constructor
         * @param id Is the ID to identify the Particle
         * @param position Is the *Vector of the position
         */
		Particle(int id, Vector* position) {
			pos= position;
			ID= id;
		}

        /**
         * Constructor for assigment
         * @param p Other Particle already defined
         */
        Particle(const Particle& p): ID(p.ID), pos(p.pos) {}

        /**
         * Destructor. It destroys the position if not already destroyed.
         */
		virtual ~Particle() {
            if(pos != NULL) { //It could be assign as NULL if it was deleted in "Atom.h" because the position of the Water is the same of the Oxygen
                delete(pos);
                pos= NULL;
            }
		}

		/**
         * It measures the distance to another Particle
         * @param p *Particle to other particle
         * @param bounds The coordinate of the last point, so the components are the width, height and length
         */
		float distanceTo(Particle* p, Vector* bounds) {
			return dist(this->pos, p->pos, bounds);
		}

};
