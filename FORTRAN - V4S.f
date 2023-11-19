       module parameters
			implicit none
			integer,parameter :: nO=216 !number of TIP3P molecules
			real dx,dy,dz,lxbox,lybox,lzbox
			real p1(3),p2(3),p3(3),p4(3)
			real v_4S
			real m_det
			real xyz_cramer(3)
			real R, vtot

			real n0
			real xOi(nO),yOi(nO),zOi(nO)
			real xH1i(nO),yH1i(nO),zH1i(nO)
			real xH2i(nO),yH2i(nO),zH2i(nO)
		end module parameters

		       program main
		use parameters
		implicit none

		lxbox= 50
		lybox= 50
		lzbox= 50

		call tetrahedron(16.04,12.27,17.50,
     &			15.69,13.07,17.91,15.66,11.56,18.02)

		print*, p1
		print*, p2
		print*, p3
		print*, p4

       end program main

       subroutine da(x1,y1,z1,x2,y2,z2)
			use parameters
			implicit none
			real x1,y1,z1,x2,y2,z2

			dx = x1 - x2
			dy = y1 - y2
			dz = z1 - z2

			if( dx .gt.  lxbox/2 ) dx = dx - lxbox
			if( dx .lt. -lxbox/2 ) dx = dx + lxbox
			if( dy .gt.  lybox/2 ) dy = dy - lybox
			if( dy .lt. -lybox/2 ) dy = dy + lybox
			if( dz .gt.  lzbox/2 ) dz = dz - lzbox
			if( dz .lt. -lzbox/2 ) dz = dz + lzbox
       end

   	   subroutine matrixDet(m)
   	   		use parameters
   	   		implicit none
   	   		real m(3,3)
			m_det= m(1,1)*m(2,2)*m(3,3) + m(2,1)*m(3,2)*m(1,3) + 
     &			m(3,1)*m(1,2)*m(2,3) - m(1,3)*m(2,2)*m(3,1) - m(2,3)*
     & 			m(3,2)*m(1,1) - m(3,3)*m(1,2)*m(2,1)
   	   end

   	   subroutine cramerRule(m,a)
   	   		use parameters
   	   		implicit none
   	   		real m(3,3),a(3)
   	   		real copy_m(3,3)
   	   		real Ma
   	   		integer*1 i,j

   	   		call matrixDet(m)
   	   		Ma= m_det

   	   		do i=1,3
	   	   		copy_m= m
	   	   		do j=1,3
					copy_m(i,j)= a(j)
				end do
				call matrixDet(copy_m)
				xyz_cramer(i)= m_det/Ma
			end do
   	   end

       subroutine v4S(ID,xO,yO,zO,xH1,yH1,zH1,xH2,yH2,zH2)
			use parameters
			implicit none
			integer ID
			real xO,yO,zO, xH1,yH1,zH1, xH2,yH2,zH2
			real x,y,z
			real dist1,dist2,dist3,dist4
			real sum(4)
			integer p_close,j,p_4S

			call tetrahedron(xO,yO,zO,xH1,yH1,zH1,xH2,yH2,zH2)

			sum(1)= 0.
			sum(2)= 0.
			sum(3)= 0.
			sum(4)= 0.

			do j=1,nO
				if(j.eq.ID) cycle

				x= xOi(j)
				y= yOi(j)
				z= zOi(j)

				call da(x,y,z,xOi(j),yOi(j),zOi(j))
				dist1= sqrt(dx**2 + dy**2 + dz**2)
				if(dist1.gt.6.) cycle

				call da(x,y,z,p1(1),p1(2),p1(3))
				dist1= sqrt(dx**2 + dy**2 + dz**2)
				call da(x,y,z,p2(1),p2(2),p2(3))
				dist2= sqrt(dx**2 + dy**2 + dz**2)
				call da(x,y,z,p3(1),p3(2),p3(3))
				dist3= sqrt(dx**2 + dy**2 + dz**2)
				call da(x,y,z,p4(1),p4(2),p4(3))
				dist4= sqrt(dx**2 + dy**2 + dz**2)

				p_close= 0
				if((dist1.lt.5.).and.(((dist1.lt.dist2)
     & .and.(dist1.lt.dist3)).and.(dist1.lt.dist4))) then
					p_close= 1
				end if
				if((dist2.lt.5.).and.(((dist2.lt.dist1)
     & .and.(dist2.lt.dist3)).and.(dist2.lt.dist4))) then
					p_close= 2
				end if
				if((dist3.lt.5.).and.(((dist3.lt.dist1)
     & .and.(dist3.lt.dist2)).and.(dist3.lt.dist4))) then
					p_close= 3
				end if
				if((dist4.lt.5.).and.(((dist4.lt.dist1)
     & .and.(dist4.lt.dist2)).and.(dist4.lt.dist3))) then
					p_close= 4
				end if

				call potWat(ID,j)
				sum(p_close)= sum(p_close)+vtot

			end do

			p_4S= 1
			do j=2,4
				if(sum(p_4S).lt.sum(j)) then
					p_4S= j
				end if
			end do

			v_4S= sum(p_4S)

       end

       subroutine tetrahedron(xO,yO,zO,xH1,yH1,zH1,xH2,yH2,zH2)
			use parameters
			implicit none
			real xO,yO,zO, xH1,yH1,zH1, xH2,yH2,zH2

			real,parameter :: R_perfect = 1.
			real,parameter :: theta= acos(-1./3.)/2

			real dOH1, dOH2, dH1H2
			real phy, magn_nu, magn_B
			real nu_A,nu_B,nu_C,nu_OH(3)
			real k1,k2
			real xB,yB,zB
			real m(3,3), a(3)

			integer*1 j
			real arrO(3)

			real mH(3), delta, mL(3)
			!END OF PARAMETERS

			!Calculate the distances
			call da(xO,yO,zO,xH1,yH1,zH1)
			dOH1= sqrt(dx**2 + dy**2 + dz**2)
			call da(xO,yO,zO,xH2,yH2,zH2)
			dOH2= sqrt(dx**2 + dy**2 + dz**2)
			call da(xH1,yH1,zH1,xH2,yH2,zH2)
			dH1H2= sqrt(dx**2 + dy**2 + dz**2)

			!Get the angle of the molecule (the half)
			phy= acos((dOH1**2 + dOH2**2 - dH1H2**2)/(2*dOH1*dOH2)) /2

			!Get the normal vector to the plane of the three atoms
			nu_A= ((yH1-yO)*(zH2-zO) - (yH2-yO)*(zH1-zO))
			nu_B= ((xH2-xO)*(zH1-zO) - (xH1-xO)*(zH2-zO))
			nu_C= ((xH1-xO)*(yH2-yO) - (xH2-xO)*(yH1-yO))
			magn_nu= sqrt(nu_A**2 + nu_B**2 + nu_C**2)
			nu_A= nu_A/magn_nu
			nu_B= nu_B/magn_nu
			nu_C= nu_C/magn_nu
			nu_OH= reshape((/ nu_A, nu_B, nu_C /),shape(nu_OH))

			!Get the bisector of the angle H-O-H
			xB= (xH1-xO)/dOH1 + (xH2-xO)/dOH2
			yB= (yH1-yO)/dOH1 + (yH2-yO)/dOH2
			zB= (zH1-zO)/dOH1 + (zH2-zO)/dOH2
			magn_B= sqrt(xB**2 + yB**2 + zB**2)

			k1= R_perfect*magn_B*cos(theta)
			k2= R_perfect*dOH1*cos(theta-phy)
			a= reshape((/ k1, k2, 0. /),shape(a))
			arrO= reshape((/ xO, yO, zO /),shape(arrO))

			!Calculate the h1 vector -> p1
			m= reshape((/ xB,yB,zB,xH1-xO,yH1-yO,
     &			zH1-zO,nu_A,nu_B,nu_C /),shape(m))

			!m= reshape((/ xB,xH1-xO,nu_A,yB,yH1-yO,
    ! &			nu_B,zB,zH1-zO,nu_C /),shape(m))
			call cramerRule(m,a)
			do j=1,3
				p1(j)= xyz_cramer(j)+arrO(j)
			end do

			!Calculate the h2 vector -> p2
			a(2)= R_perfect*dOH2*cos(theta-phy)
			m(1,2)= xH2-xO
			m(2,2)= yH2-yO
			m(3,2)= zH2-zO
			call cramerRule(m,a)
			do j=1,3
				p2(j)= xyz_cramer(j)+arrO(j)
			end do

			!Get the midpoint and mid-distance
			do j=1,3
				mH(j)= (p1(j)+p2(j) ) / 2
			end do
			call da(p1(1),p1(2),p1(3),mH(1),mH(2),mH(3))
			delta= sqrt(dx**2 + dy**2 + dz**2)

			!Get the midpoint between the other two points
			do j=1,3
				mL(j)= 2*arrO(j) - mH(j)
			end do

			!Calculate the other two points
			do j=1,3
				p3(j)=  delta*nu_OH(j)+mL(j)
				p4(j)= -delta*nu_OH(j)+mL(j)
			end do

       end

       subroutine potWat (x,y)
			use parameters
			implicit none
			integer x,y,xx,yy
			real vc, vlj

			real,parameter  :: qe= 1.
			real,parameter  :: kcoul=1389.35458 !kJ mol-1 A e-2

c   	   !TIP3P parameters
			real,parameter ::Os = 3.15061  !sigma/A
			real,parameter ::Oe = 0.6364   !Oxygen epsilon L-J/kJ mol-1
			real,parameter ::Oc = -0.8340  !charge oxygen
			real,parameter ::Hc = 0.4170   !charge hydrogen

			vc=0.
			VLJ=0.
			vtot=0.

			xx=x
			yy=y

c   	  !H1j contributions  (only coulomb interactions cause eH=0)
			call da(xH1i(xx),yH1i(xx),zH1i(xx),xH1i(yy),yH1i(yy),zH1i(yy))
			R = sqrt(dx**2 + dy**2 + dz**2)
			Vc=(Kcoul*qe*Hc**2)/R !coulomb H1-H1     V = q2 / 4 pi e0 r
			Vtot=Vtot+Vc

			call da(xH2i(xx),yH2i(xx),zH2i(xx),xH1i(yy),yH1i(yy),zH1i(yy))
			R = sqrt(dx**2 + dy**2 + dz**2)
			Vc=(Kcoul*qe*Hc**2)/R !coulomb H1-H2
			Vtot=Vtot+Vc

			call da(xOi(xx),yOi(xx),zOi(xx),xH1i(yy),yH1i(yy),zH1i(yy))
			R = sqrt(dx**2 + dy**2 + dz**2)
			Vc=(Kcoul*qe*Hc*Oc)/R !coulomb H1-O
			Vtot=Vtot+Vc


c   	  !H2j contributions  (only coulomb interactions cause eH=0)
			call da(xH1i(xx),yH1i(xx),zH1i(xx),xH2i(yy),yH2i(yy),zH2i(yy))
			R = sqrt(dx**2 + dy**2 + dz**2)
			Vc=(Kcoul*qe*Hc**2)/R !coulomb H2-H1
			Vtot=Vtot+Vc

			call da(xH2i(xx),yH2i(xx),zH2i(xx),xH2i(yy),yH2i(yy),zH2i(yy))
			R = sqrt(dx**2 + dy**2 + dz**2)
			Vc=(Kcoul*qe*Hc**2)/R !coulomb H2-H2
			Vtot=Vtot+Vc

			call da(xOi(xx),yOi(xx),zOi(xx),xH2i(yy),yH2i(yy),zH2i(yy))
			R = sqrt(dx**2 + dy**2 + dz**2)
			Vc=(Kcoul*qe*Hc*Oc)/R !coulomb H2-O
			Vtot=Vtot+Vc


c   	  !Oj contributions  (coulomb+L-J)
c   	  !Coulomb
			call da(xH1i(xx),yH1i(xx),zH1i(xx),xOi(yy),yOi(yy),zOi(yy))
			R = sqrt(dx**2 + dy**2 + dz**2)
			Vc=(Kcoul*qe*Oc*Hc)/R !coulomb O-H1
			Vtot=Vtot+Vc

			call da(xH2i(xx),yH2i(xx),zH2i(xx),xOi(yy),yOi(yy),zOi(yy))
			R = sqrt(dx**2 + dy**2 + dz**2)
			Vc=(Kcoul*qe*Hc*Oc)/R !coulomb O-H2
			Vtot=Vtot+Vc

			call da(xOi(xx),yOi(xx),zOi(xx),xOi(yy),yOi(yy),zOi(yy))
			R = sqrt(dx**2 + dy**2 + dz**2)
			Vc=(Kcoul*qe*Oc**2)/R !coulomb O-O
			Vtot=Vtot+Vc

c   	  !Lennard - Jones
			VLJ=VlJ+4*Oe*((Os/R)**12-(Os/R)**6)
			vtot=vtot+VLJ
		end

       subroutine potOthers(x,y,sigma,epsilon,charge)
			use parameters
			implicit none
			real sigma,epsilon,charge,Os,Oe
			real vc, vlj
			integer x,y,xx,yy

			real,parameter  :: qe= 1.
			real,parameter  :: kcoul=1389.35458 !kJ mol-1 A e-2

c   	   !TIP3P parameters
			real,parameter ::Oc = -0.8340  !charge oxygen
			real,parameter ::Hc = 0.4170   !charge hydrogen
			Os = sigma
			Oe = epsilon

			vc=0.
			VLJ=0.
			vtot=0.

			xx=x
			yy=y

c   	  !H1j contributions  (only coulomb interactions cause eH=0)
			call da(xH1i(xx),yH1i(xx),zH1i(xx),xH1i(yy),yH1i(yy),zH1i(yy))
			R = sqrt(dx**2 + dy**2 + dz**2)
			Vc=(Kcoul*qe*Hc*charge)/R !coulomb H1-H1     V = q2 / 4 pi e0 r
			Vtot=Vtot+Vc


c   	  !H2j contributions  (only coulomb interactions cause eH=0)
			call da(xH1i(xx),yH1i(xx),zH1i(xx),xH2i(yy),yH2i(yy),zH2i(yy))
			R = sqrt(dx**2 + dy**2 + dz**2)
			Vc=(Kcoul*qe*Hc*charge)/R !coulomb H2-H1
			Vtot=Vtot+Vc

c   	  !Oj contributions  (coulomb+L-J)
c   	  !Coulomb
			call da(xH1i(xx),yH1i(xx),zH1i(xx),xOi(yy),yOi(yy),zOi(yy))
			R = sqrt(dx**2 + dy**2 + dz**2)
			Vc=(Kcoul*qe*Oc*charge)/R !coulomb O-H1
			Vtot=Vtot+Vc

c   	  !Lennard - Jones
			VLJ=VlJ+4*Oe*((Os/R)**12-(Os/R)**6)
			vtot=vtot+VLJ
       end
