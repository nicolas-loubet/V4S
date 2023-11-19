# V4S
This is a guide related to the new V4S parameter, a structural indicator introduced in https://doi.org/10.48550/arXiv.2311.08087
---
abstract: |
  This document has for aim explaining how we compute the
  *V*<sub>4*S*</sub> values for any molecule in a system, providing the
  code in FORTRAN and C++. Those should be adapted for each particular
  case. Firstly, we should explain the calculation of the tetrahedron.
  Secondly, the two algorithm should be explained in a generic way, so
  the same idea could be used in other cases.
author:
- Nicolás A. Loubet<sup>1</sup>
- Alejandro R. Verde<sup>1</sup>
- Gustavo A. Appignanesi<sup>1\*</sup>
title: README
---

# The perfect tetrahedron

In this first section, we will explain how to calculate the four points
of the perfect tetrahedron for a water molecule. This is intended for
water models such as TIP3P, TIP4P, TIP5P, SPC, SPC/E, among others. The
model must have 3 coordinates for each molecule: The two hydrogens and
the oxygen. This last one will be the center of the tetrahedron. The
necessary steps can be summarized in two objectives:

1.  Open the angle H-O-H (not necessary if the angle is already the
    desirable *a**r**c**c**o**s*(−1/3), the perfect tetrahedron angle
    109.471<sup>∘</sup> = 1, 9106*r**a**d*).

2.  Locate the opposite two points that complete the tetrahedron.

In this sense, from the three initial coordinates we can find the two
vectors that conform the molecule, that we should name
$\\overrightarrow{OH_1}$ and $\\overrightarrow{OH_2}$.

$$\\overrightarrow{OH_i}=\\overrightarrow{H_i}-\\overrightarrow{O}$$

From this two vectors we can find the bisector. This one is parallel to
the dipole moment.

$$\\overrightarrow{b}=\\frac{\\overrightarrow{OH_1}} {\\left\|\\overrightarrow{OH_1}\\right\|} + \\frac{\\overrightarrow{OH_2}} {\\left\|\\overrightarrow{OH_2}\\right\|}$$

Also, we can find a normal vector perpendicular to the two
*O**H*<sub>*i*</sub> vectors, using the cross product.

$$\\overrightarrow{\\eta}= \\frac{\\overrightarrow{OH_1}\\times\\overrightarrow{OH_2}} {\\left\|\\overrightarrow{OH_1}\\times\\overrightarrow{OH_2}\\right\|}$$

<span id="fig1" label="fig1">\[fig1\]</span>

Now, we can find the parameters of the plane
(*A**x* + *B**y* + *C**z* + *D* = 0) that contains the oxygen and the
two hydrogens, and, therefore, the two points corresponding to the
hydrogens when the angle is open to the perfect tetrahedron angle. As we
are using the vectors as if oxygen was the center of coordinates, then
we must set
*D* =  − *A**x*<sub>0</sub> − *B**y*<sub>0</sub> − *C**z*<sub>0</sub> = 0.
Note that *η⃗* = *A**x̂* + *B**ŷ* + *C**ẑ*

For computing the point i of the tetrahedron, the nearest to the
*H*<sub>*i*</sub>, we should note that the dot products are:

$$\\begin{aligned}
    \\overrightarrow{oh_i} \\cdot \\overrightarrow{b} &= R \\left\| \\overrightarrow{b} \\right\| cos(\\theta)\\\\
    \\overrightarrow{oh_i} \\cdot \\overrightarrow{OH_i} &= R \\left\| \\overrightarrow{OH_i} \\right\| cos(\\theta-\\phi)\\end{aligned}$$
Being R the distance between the oxygen and each point of the
tetrahedron, that we defined as 1Å.

From this data (and the plane definition), we can find the coordinates
of the *o**h*<sub>*i*</sub> vectors (for example, with the Cramer’s
rule). Then, we can find the perfect point (*h*<sub>1, 2</sub>) for each
one.

$$M=
    \\begin{bmatrix}
        x\_{\\overrightarrow{b}} & y\_{\\overrightarrow{b}} & z\_{\\overrightarrow{b}} \\\\
        x\_{\\overrightarrow{OH_i}} & y\_{\\overrightarrow{OH_i}} & z\_{\\overrightarrow{OH_i}} \\\\
        A & B & C \\\\
    \\end{bmatrix}$$
$$M_x=
    \\begin{bmatrix}
        R \\left\| \\overrightarrow{b} \\right\| cos(\\theta) & y\_{\\overrightarrow{b}} & z\_{\\overrightarrow{b}} \\\\
        R \\left\| \\overrightarrow{OH_i} \\right\| cos(\\theta-\\phi) & y\_{\\overrightarrow{OH_i}} & z\_{\\overrightarrow{OH_i}} \\\\
        0 & B & C \\\\
    \\end{bmatrix}$$
$$M_y=
    \\begin{bmatrix}
        x\_{\\overrightarrow{b}} & R \\left\| \\overrightarrow{b} \\right\| cos(\\theta) & z\_{\\overrightarrow{b}} \\\\
        x\_{\\overrightarrow{OH_i}} & R \\left\| \\overrightarrow{OH_i} \\right\| cos(\\theta-\\phi) & z\_{\\overrightarrow{OH_i}} \\\\
        A & 0 & C \\\\
    \\end{bmatrix}$$
$$M_z=
    \\begin{bmatrix}
        x\_{\\overrightarrow{b}} & y\_{\\overrightarrow{b}} & R \\left\| \\overrightarrow{b} \\right\| cos(\\theta) \\\\
        x\_{\\overrightarrow{OH_i}} & y\_{\\overrightarrow{OH_i}} & R \\left\| \\overrightarrow{OH_i} \\right\| cos(\\theta-\\phi) \\\\
        A & B & 0
    \\end{bmatrix}$$
$$\\begin{aligned}
    \\overrightarrow{oh_i} &=  \\frac{\|M_x\|}{\|M\|} \\hat{x}+ \\frac{\|M_y\|}{\|M\|} \\hat{y}+  \\frac{\|M_z\|}{\|M\|} \\hat{z}\\\\
    \\overrightarrow{h_i} &= \\overrightarrow{oh_i} + \\overrightarrow{O}\\end{aligned}$$

<span id="fig2" label="fig2">\[fig2\]</span>

At this point we finished the first part, finding two points. If the
water model has already the right angle (such as SPC, SPC/E, or OPC3),
you can start from here.

Using the image as a guide, we should consider the fact that the middle
points of the segments between the *h*<sub>*i*</sub> and
*e*<sub>*i*</sub> define a line that contains the oxygen point. Also,
the distances to this one are equal. So, we can consider that
$\\overrightarrow{m_h} - \\overrightarrow{O} = - \\left( \\overrightarrow{m_e} - \\overrightarrow{O} \\right)$
from this, we have
$\\overrightarrow{m_e} = 2\\overrightarrow{O} - \\overrightarrow{m_h}$.
Also, we can define the distance from any *h*<sub>*i*</sub> to the
*m*<sub>*h*</sub> (that is equal to the distance between the
*e*<sub>*i*</sub> and the *m*<sub>*e*</sub>) as *δ*.

Also, we need to find a vector that is parallel to the line formed by
the two *e*<sub>*i*</sub> points, and the *m*<sub>*e*</sub>. As it must
be perpendicular to the plane formed by the oxygen and the two
*h*<sub>*i*</sub> points, then we have already found *η*. With this
considerations, then:

$$\\begin{aligned}
    \\overrightarrow{m_h} &=  \\frac{1}{2} \\left( \\overrightarrow{h_1} + \\overrightarrow{h_2} \\right)\\\\
    \\delta &= \\left\| \\overrightarrow{h_i} - \\overrightarrow{m_h} \\right\|\\\\
    \\overrightarrow{m_e} &= 2\\overrightarrow{O} - \\overrightarrow{m_h}\\\\
    \\overrightarrow{e\_{1,2}} &= \\pm \\delta \\overrightarrow{\\eta} + \\overrightarrow{m_e}\\end{aligned}$$

This way, we have found the four points from the water molecule.

# Generic algorithm

In this second section, we should explain how to compute the
*V*<sub>4*S*</sub> of a molecule i, in a system of *N*<sub>*w*</sub>
water molecules and *N*<sub>*a*</sub> atoms that are not water.

1.  We have to compute the four points of the perfect tetrahedron, using
    the equations explained before.

2.  We create 4 variables where we will calculate the total potential
    energy for each one of the four points, initialized as 0.

3.  In a loop, we consider each water molecule and atom in the system,
    different from the central molecule.

4.  We calculate the distance between this atom/molecule and the water
    molecule i. If it is more than 6Å, we skip it.

5.  Then, we iterate over each tetrahedron point, and calculate the
    distance between the atom/molecule and the tetrahedron point.
    Comparing the four distances, we find the one that is closest to the
    atom/molecule.

6.  If the distance from this closest point to the atom/molecule is 5Åor
    less, we calculate the potential energy interaction between the
    central molecule and this atom/molecule, considering Coulombic and
    Lennard-Jones potentials. This value is added to the variable
    created in the 2nd step corresponding to to this point.

7.  We repeat from step 4 with the next molecule/atom.

8.  When we have all the contributions, we find the greatest value
    stored in the variables. This value is *V*<sub>4*S*</sub>. The other
    values are reffered as *V*<sub>1*S*</sub>, *V*<sub>2*S*</sub> and
    *V*<sub>3*S*</sub>, that could be also useful in some cases.

We show in this GitHub some examples coded in FORTRAN and C++, as a
guide. Those files should be adapted to your particular case.

