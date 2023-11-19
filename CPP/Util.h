#pragma once
#include "Vector.h"
#include <math.h>
#include <string>
#include <fstream>
using namespace std;

#define ERROR_VALUE -3

/**
 * Version: August 2023
 */

/**
 * Calculates the distance between two coordinates (vectors) given the bounds to use periodic boundary conditions
 * @param c1 One of the coordinates
 * @param c2 The other coordinate
 * @param bounds The coordinate of the last point, so the components are the width, height and length
 * @return The distance between the two coordinates
 */
float dist(Vector* c1, Vector* c2, Vector* bounds) {
    Vector* dif= *c1-*c2;
    Vector* b_1= *bounds/2;
    Vector* b_2= *bounds/(-2);

    if(dif->x > b_1->x) dif->x-= bounds->x;
    if(dif->x < b_2->x) dif->x+= bounds->x;
    if(dif->y > b_1->y) dif->y-= bounds->y;
    if(dif->y < b_2->y) dif->y+= bounds->y;
    if(dif->z > b_1->z) dif->z-= bounds->z;
    if(dif->z < b_2->z) dif->z+= bounds->z;

    float output= dif->magnitude();

    delete(dif);
    delete(b_1);
    delete(b_2);

    return output;
}

/**
 * Calculates the angle that 3 points form in a 3D system
 * @param c1 One of the point in the edges
 * @param c2 The center point
 * @param c3 The other point in an edge
 * @param bounds The coordinate of the last point, so the components are the width, height and length
 * @return The angle in radians formed by the 3 points
 */
float getAngle(Vector* c1, Vector* c2, Vector* c3, Vector* bounds) {
    float a= dist(c1,c3,bounds); //Opposite to the angle
    float b= dist(c1,c2,bounds);
    float c= dist(c2,c3,bounds);
    //Law of cosines
    return abs(acos((pow(b,2)+pow(c,2)-pow(a,2))/(2*b*c)));
}

/**
 * Function that combines two sorted arrays creating a sorted array
 * @param arr A *float that corresponds to the first member of the array
 * @param p_inic The position of the first element of the first subarray
 * @param p_medium The position of the last element of the first subarray (and the first element of the second subarray)
 * @param p_fin The position of the last element of the second subarray
 * @return true if it succeded, false if error
 */
bool merge(float* arr, const int p_inic, const int p_medium, const int p_fin) {
    if(p_inic > p_medium || p_medium > p_fin)
        return -1;

    float arr_aux[p_fin-p_inic+1];
    int p_first_half, p_second_half, p_aux;

    p_first_half= p_inic;
    p_second_half= p_medium+1;
    p_aux= 0;

    //Sorts the elements of the two halves in the auxiliary array
    while(p_first_half <= p_medium && p_second_half <= p_fin) {
        if(arr[p_first_half] <= arr[p_second_half])
            arr_aux[p_aux++]= arr[p_first_half++];
        else
            arr_aux[p_aux++] = arr[p_second_half++];
    }

    while(p_first_half <= p_medium)
        arr_aux[p_aux++] = arr[p_first_half++];
    while(p_second_half <= p_fin)
        arr_aux[p_aux++] = arr[p_second_half++];

    //Copies the elements of the auxiliary array in the original array
    for(int i= 0; i < p_fin-p_inic+1; i++)
        arr[p_inic+i]= arr_aux[i];

    return true;
}

/**
 * Function that divides recursively an array in 2. Then, changes positions so they are sorted. Finally, combines the halves.
 * @param arr A *float that corresponds to the first member of the array
 * @param p_inic The position of the first element of the subarray
 * @param p_fin The position of the last element of the subarray
 * @return true if it succeded, false if error
 */
bool splitMerge(float* arr, const int p_inic, const int p_fin) {
    if(p_inic > p_fin)
        return false;
    const int p_medium= (p_inic + p_fin)/2;

    if(p_inic != p_fin) {
        splitMerge(arr, p_inic, p_medium); //First half
        splitMerge(arr, p_medium+1, p_fin); //Second half
        return merge(arr, p_inic, p_medium, p_fin); //Merge halves
    }
    return true;
}

/**
 * Function that sorts a floats array using merge sort. O(n log(n))
 * @param arr A *float that corresponds to the first member of the array
 * @param size The number of elements of the array
 * @return true if it succeded, false if error
 */
bool sort(float* arr, int size){
    if(size < 0)
        return false;
    return splitMerge(arr, 0, size-1);
}

/**
 * Calculates the determinant of a 3x3 matrix with the Sarrus Rule
 * @param matrix A 3x3 float array
 * @return The matrix's determinant
 */
float determinant_3x3(float matrix[3][3]) {
    return matrix[0][0]*matrix[1][1]*matrix[2][2] + matrix[0][1]*matrix[1][2]*matrix[2][0] + matrix[0][2]*matrix[1][0]*matrix[2][1]
         - matrix[0][2]*matrix[1][1]*matrix[2][0] - matrix[0][0]*matrix[1][2]*matrix[2][1] - matrix[0][1]*matrix[1][0]*matrix[2][2];
}

/**
 * Calculates the x,y,z values (a Vector) obteined from the Cramer's Rule
 * @param matrix A 3x3 float array of the base matrix
 * @param m_independents A float array of 3 values corresponding to the independt values
 * @return The Vector with the x, y and z values obtained. (The Vector must be removed with "delete()")
 */
Vector* CramersRule(float matrix[3][3], float m_independents[3]) {
        float A= determinant_3x3(matrix);
        float Ax_matrix[3][3]= {
            { m_independents[0] , matrix[0][1] , matrix[0][2] },
            { m_independents[1] , matrix[1][1] , matrix[1][2] },
            { m_independents[2] , matrix[2][1] , matrix[2][2] }
        };
        float Ax= determinant_3x3(Ax_matrix);
        float Ay_matrix[3][3]= {
            { matrix[0][0] , m_independents[0] , matrix[0][2] },
            { matrix[1][0] , m_independents[1] , matrix[1][2] },
            { matrix[2][0] , m_independents[2] , matrix[2][2] }
        };
        float Ay= determinant_3x3(Ay_matrix);
        float Az_matrix[3][3]= {
            { matrix[0][0] , matrix[0][1] , m_independents[0] },
            { matrix[1][0] , matrix[1][1] , m_independents[1] },
            { matrix[2][0] , matrix[2][1] , m_independents[2] }
        };
        float Az= determinant_3x3(Az_matrix);
        return new Vector(Ax/A, Ay/A, Az/A);
}
