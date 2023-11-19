
/**
 * Version: August 2023
 */

struct Tetrahedron {
    Vector* H1;
    Vector* H2;
    Vector* L1;
    Vector* L2;
};

Vector** checkPBC(Vector* O, Vector* H1, Vector* H2, Vector* bounds) {
    Vector** output= new Vector*[3];

    output[0]= new Vector(O->x,O->y,O->z);
    Vector* H_i[2]= {H1, H2};
    for(int i= 0; i < 2; i++) {
        Vector* dif= *H_i[i]-*O;
        Vector* b_1= *bounds/2;
        Vector* b_2= *bounds/(-2);

        if(dif->x > b_1->x) dif->x-= bounds->x;
        if(dif->x < b_2->x) dif->x+= bounds->x;
        if(dif->y > b_1->y) dif->y-= bounds->y;
        if(dif->y < b_2->y) dif->y+= bounds->y;
        if(dif->z > b_1->z) dif->z-= bounds->z;
        if(dif->z < b_2->z) dif->z+= bounds->z;

        output[i+1]= *O+*dif;

        delete(dif);
        delete(b_1);
        delete(b_2);
    }
    return output;
}

Tetrahedron getPerfectTetrahedron(Vector* O_real, Vector* H1_real, Vector* H2_real, float phy, Vector* bounds) {
    const float R= 1; //Anstrongs between O and perfect vertices
    const float theta= acos(-1./3.)/2; //Perfect angle for tetrahedron
    phy/= 2; //I need the angle between OH and b

    Vector** pbc_vectors= checkPBC(O_real, H1_real, H2_real, bounds);
    Vector* O= pbc_vectors[0];
    Vector* H1= pbc_vectors[1];
    Vector* H2= pbc_vectors[2];
    delete(pbc_vectors);

    //I step
    Vector* OH1= (*H1-*O);
    Vector* OH2= (*H2-*O);
    Vector* OH1_norm= *OH1*(1./OH1->magnitude());
    Vector* OH2_norm= *OH2*(1./OH2->magnitude());

    //II step
    Vector* b= *OH1_norm + *OH2_norm;
    delete(OH1_norm);
    delete(OH2_norm);

    //III step
    Vector* nu_OH_sin_norm= *OH1 % *OH2;
    Vector* nu_OH= *nu_OH_sin_norm*(1./nu_OH_sin_norm->magnitude());
    delete(nu_OH_sin_norm);

    //IV step
    Vector* hidrogenos[2]= {OH1, OH2};
    Vector* h[2];
    for(int i= 0; i < 2; i++) {
        Vector* OHi= hidrogenos[i];
        float k1= R*b->magnitude()*cos(theta);
        float k2= R*OHi->magnitude()*cos(theta-phy);

        float A_matriz[3][3]= {
            {   b->x   ,   b->y   ,   b->z   },
            {  OHi->x  ,  OHi->y  ,  OHi->z  },
            { nu_OH->x , nu_OH->y , nu_OH->z }
        };
        delete(OHi); //delete OH1,2
        float m_independientes[3]= { k1, k2, 0. };

        Vector* oh= CramersRule(A_matriz, m_independientes);
        h[i]= *oh + *O;
        delete(oh);
    }
    delete(b);

    //V step
    Vector* suma_h= *h[0]+*h[1];
    Vector* m_H= *suma_h*(0.5);
    delete(suma_h);
    Vector* mH_H= *h[0] - *m_H;
    float delta= mH_H->magnitude();
    delete(mH_H);

    //VI step
    Vector* O_2= *O*2.;
    Vector* m_L= *O_2 - *m_H;
    delete(O_2);
    delete(m_H);

    //VII step
    Vector* mL_L1= *nu_OH*delta;
    Vector* mL_L2= *nu_OH*(-delta);
    delete(nu_OH);
    Vector* L1= *mL_L1 + *m_L;
    Vector* L2= *mL_L2 + *m_L;
    delete(mL_L1);
    delete(mL_L2);
    delete(m_L);

    Tetrahedron output;
    output.H1= h[0];
    output.H2= h[1];
    output.L1= L1;
    output.L2= L2;

    delete(O);
    delete(H1);
    delete(H2);

    return output;
}
