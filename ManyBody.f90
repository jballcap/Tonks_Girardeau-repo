


!**************************************************************************************
!				MAIN PROGRAM
!**************************************************************************************
PROGRAM sequation 
IMPLICIT NONE
  INTEGER, PARAMETER 	::  n=2*1000, npart=240, nstates=300, Nt=2*100  !grid points, # of particles, # of states
  INTEGER 		:: i, j, np
  REAL(8)		:: xmax, summ,dx, kmax, pi, tmax, lav
  REAL(8), DIMENSION(n,n)	:: Ham0, Hamev  !Hamiltonian
  REAL(8), DIMENSION(n, nstates):: u0, uev
  REAL(8), DIMENSION(n)		:: dens0, densev  !single particle density
  REAL(8), DIMENSION(n)		:: xarr, karr !k and x in the mesh
  REAL(8), DIMENSION(npart)     :: groundfidel
  REAL(8), DIMENSION(Nt)	:: t
  REAL(8), DIMENSION(Nt,npart)  :: fidpart
 

!Defining variables for LAPACK
  INTEGER, PARAMETER 	:: LDA=2*1000  !Leading dimension 
  INTEGER :: INFO  		 !variable that gets info from LAPACK
  INTEGER :: LWORK 		 !will be the size for the workspace
  INTEGER, PARAMETER :: LWMAX=10000 	 !maximun size of the workspace 
  REAL(8), DIMENSION(LWMAX) :: WORK, WORK1 !workspace of LAPACK subroutine
  REAL(8), DIMENSION(n) :: W0, Wev	 !vector of eigenvalues
  CHARACTER(LEN=1) :: COMPE, UT     !characters to interact with LAPACK
  
 
  tmax= 25.0 !maximun time
  xmax=15.0 !size of the grid
  pi=4.0*ATAN(1.0)
  dx=2.0*(xmax/(1.0*n))

  kmax=pi/dx


OPEN(UNIT=7, FILE='eigenvalues0.txt', STATUS='replace', ACTION='write')
OPEN(UNIT=8, FILE='eigenvectors0.txt', STATUS='replace', ACTION='write')
OPEN(UNIT=9, FILE='eigenvalues_ev.txt', STATUS='replace', ACTION='write')
OPEN(UNIT=10, FILE='eigenvectors_ev.txt', STATUS='replace', ACTION='write')
OPEN(UNIT=11, FILE='sing_part_dens0.txt', STATUS='replace', ACTION='write')
OPEN(UNIT=12, FILE='sing_part_densev.txt', STATUS='replace', ACTION='write')
OPEN(UNIT=13, FILE='ground_fidel.txt', STATUS='replace', ACTION='write')
OPEN(UNIT=14, FILE='fidelity.txt', STATUS='replace', ACTION='write')


!creating the mesh and constructing the matrix elements of H=K+V
  CALL fgrid(n, xmax, xarr, karr, Ham0, Hamev)


  COMPE="V"  !V to tell lapack to compute eigenvalues and eigenvectors
  UT="U"
 
  LWORK=-1
  CALL DSYEV(COMPE, UT,  n, Ham0, LDA, W0, WORK, LWORK, INFO  )
  LWORK= MIN(LWMAX, INT(WORK(1)))

!solving eigenvalue problem
  CALL DSYEV(COMPE, UT, n, Ham0, LDA, W0, WORK, LWORK, INFO)

   IF( INFO.GT.0 ) THEN
    WRITE(7, *) 'The algorithm failed to compute eigenvalues.'
     STOP
   ELSE IF (INFO.EQ.0) THEN 
!Exporting eigenvalues('eigenvalues0.txt')
    DO i=1,nstates
     WRITE(7,'(I5, F11.4)') i, W0(i)
    END DO
   END IF  
   CLOSE(7)

!Choosing the quantity of eigenvectors(nstates)  INITIAL eigenvectors
  summ=0.0
  DO i=1, nstates
   DO j=1, n 
    summ=summ+(Ham0(j,i))**2
   END DO
  u0(:, i)=Ham0(:,i)/((SQRT(summ*dx)))  !check: product with dx, doesnt normlize
  summ=0.0
  END DO
!Exporting eigenvectors('eigenvectors0.txt')
   DO i=1,nstates
    WRITE(8,'(240(F11.4))') u0(:,i)
   END DO
   CLOSE(8)





 LWORK=-1
 CALL DSYEV(COMPE, UT,  n, Hamev, LDA, Wev, WORK1, LWORK, INFO  )
  LWORK= MIN(LWMAX, INT(WORK1(1)))

!solving eigenvalue problem
  CALL DSYEV(COMPE, UT, n, Hamev, LDA, Wev, WORK1, LWORK, INFO)

   IF( INFO.GT.0 ) THEN
    WRITE(9, *) 'The algorithm failed to compute eigenvalues.'
     STOP
   ELSE IF (INFO.EQ.0) THEN 
!Exporting eigenvalues('eigenvalues_ev.txt')
    DO i=1,nstates
     WRITE(9,'(I5, F11.4)') i, Wev(i)
    END DO
   END IF  
   CLOSE(9)

!Choosing the quantity of eigenvectors(nstates) EVOLVED eigenvectors  
  summ=0.0
  DO i=1, nstates
   DO j=1, n 
    summ=summ+(Hamev(j,i))**2
   END DO
  uev(:, i)=Hamev(:,i)/((SQRT(summ*dx)))  !check: product with dx, doesnt normlize
  summ=0.0
  END DO
!Exporting eigenvectors('eigenvectors_ev.txt')
   DO i=1,nstates
    WRITE(10,'(240(F11.4))') uev(:,i)
   END DO
   CLOSE(10)




!Computing the single particle density
! CALL onepartdens(n, nstates, npart, u0, dens0)
! CALL onepartdens(n, nstates, npart, uev, densev)
!Exporting the values of the single particle density('sing_part_dens0.txt')
!   DO i=1, n
!    WRITE(11,'(2(F11.5))'), xarr(i), dens0(i)
!   END DO
!   CLOSE(11)
!Exporting the values of the single particle density('sing_part_densev.txt')
!   DO i=1, n
!    WRITE(12,'(2(F11.5))'), xarr(i), densev(i)
!   END DO
!   CLOSE(12)




!Computing the ground state fidelity
!  CALL groundfid(n, npart, dx, u0(1:n, 1:npart), uev(1:n, 1:npart), groundfidel)
!exporting groundfidel(i)
!  DO i=1, npart
!  WRITE(13, *) i, groundfidel(i)
!  END DO
!  CLOSE(13)


!Computing the fidelity(Loschmidt echo)
  CALL fidelity(Nt, n, npart, nstates, dx, tmax, u0(1:n, 1:nstates), uev(1:n, 1:nstates), W0(1:nstates), Wev(1:nstates), t, fidpart)
!exporting fidelitypart(time, Npart)
  DO i=1, Nt
    WRITE(14, *) t(i), fidpart(i,:) 
  END DO
  CLOSE(14)


!Calculating the average fidelity
  np=240
  CALL loschav(Nt,np, npart, tmax, fidpart, lav )
  PRINT*, lav


  
END PROGRAM


!************************************************************
!    calculate the matrix elements of the Hamiltonian
!*************************************************************
SUBROUTINE fgrid(n, xmax, xarr, karr, Ham0, Hamev)
IMPLICIT NONE
INTEGER, INTENT(IN)			:: n
REAL(8), INTENT(IN)			:: xmax
INTEGER					:: i, j
REAL(8)					:: kmax, dx, dk, pi, hbar, mass
REAL(8), DIMENSION(n), INTENT(OUT)	:: xarr, karr
REAL(8), DIMENSION(n,n)			:: T,vpot, vpot1, vpot2
REAL(8), DIMENSION(n,n), INTENT(OUT)	:: Ham0, Hamev


hbar=1 !1.0546e-34
mass=1 !87*1.673e-27
pi= 4.0*ATAN(1.0)

!creating the grid
!n !number of points
!xmax !size of the grid

dx=2.0*(xmax/(1.0*n))

kmax=pi/dx
dk=2.0*pi/(n*1.0*dx)


!constructing the fourier grid 
!First need to allocate space for the kinetic energy matrix
  DO i=0,n-1
   xarr(i+1)=i*dx !-xmax
   karr(i+1)=i*dk !-kmax
  END DO


!Constructing the kinetic energy matrix
  DO i=1,n
   DO j=1,n
   IF (i.eq.j) THEN 
    T(i,j)=((hbar**2)/(2*mass))*((kmax**2)/3)*(1+2/(n**2))
   ELSE
    T(i,j)=((hbar**2)/(2*mass))*(2*(kmax**2)/((1.0*n)**2))*((-1)**(j-i))/((SIN(pi*(j-i)/(1.0*n)))**2)
   END IF
   END DO
  END DO
  

  
  
!Calling the potential subroutine to generate v-matrix

  CALL potlat(n,xmax,xarr,karr,vpot1)
  CALL potwell(n,vpot2)

  vpot=vpot2 +vpot1

  
 Ham0=T+vpot !box+lattice 
 Hamev=T+vpot2 !box

END SUBROUTINE fgrid


!**************************************************************
!		Time average fidelity
!**************************************************************
SUBROUTINE  loschav(Nt, np, npart, tmax, fidel, lav ) !considered time, fidelity vector f(dt), losch. ave
IMPLICIT NONE
INTEGER, INTENT(IN)				:: Nt, np, npart
REAL(8), INTENT(IN)				:: tmax
REAL(8), DIMENSION(Nt, npart), INTENT(IN)       :: fidel
REAL(8), INTENT(OUT)				:: lav
INTEGER						:: i
REAL(8)						:: dt, time, summ, eps, lavdummy,dd !eps is the precision 

OPEN(UNIT=15, FILE='time_av.txt', STATUS='replace', ACTION='write')

time=0.0
summ=0.0
eps=1e-6
lavdummy=0.0
lav=0.0
dd=0.0

dt=tmax/(1.0*Nt)

  DO i=1, Nt 
   lavdummy=lav     
   summ=summ+fidel(i, np)!*dt    
   lav= (summ/(1.0*i))
   dd=ABS(lav-lavdummy)
   WRITE(15, *) i, lav , dd
    IF (dd.le.eps) THEN
     PRINT*, 'The time average of L is', lav, 'at', i ,'for', np
    END IF	
  END DO
    CLOSE(15)

END SUBROUTINE loschav



!*************************************************************
!		      Fidelity
!************************************************************
SUBROUTINE fidelity(Nt, n, npart, nstates, dx, tmax, h0, hev, eval0, evalev, t, fid)
IMPLICIT NONE
INTEGER, INTENT(IN)			  :: Nt, n, npart, nstates
REAL(8), INTENT(IN)			  :: dx, tmax
REAL(8), DIMENSION(nstates), INTENT(IN)	  :: eval0, evalev !eigenvalues matrices
REAL(8), DIMENSION(n, nstates), INTENT(IN):: h0, hev
REAL(8), DIMENSION(Nt,npart), INTENT(OUT) :: fid
REAL(8), DIMENSION(Nt), INTENT(OUT)	  :: t !vector that keeps values of t
COMPLEX*16, DIMENSION(nstates)		  :: eE0, eEev !time dependent exponential
REAL(8), DIMENSION(nstates, n)		  :: h0t, hevt 
REAL(8), DIMENSION(nstates, nstates)	  :: c0, cev 
COMPLEX*16, DIMENSION(nstates, npart)     :: c0_time, cev_time 
COMPLEX*16, DIMENSION(n, npart)		  :: psi_time, psiev_time
COMPLEX*16, DIMENSION(npart, n)		  :: psi_timeT
COMPLEX*16, DIMENSION(npart, npart)       :: prodmat
INTEGER					  :: i, j, itime
REAL(8)					  :: dt, hbar
COMPLEX*16				  :: im
COMPLEX*16, EXTERNAL			  :: determinant




!Defining some time parameters
im=(0.0, 1.0)
hbar=1.0
dt= (1.0*tmax)/(1.0*Nt)


 

  DO i=1, Nt
   t(i)=(i-1)*dt
  END DO

!constructing the overlaps between eigenfunctions
  DO i=1, n
   DO j=1, nstates
   h0t(j,i)=h0(i,j)  !trasposed matrices
   hevt(j,i)=hev(i,j)
   END DO
  END DO
 c0=dx*(MATMUL(h0t, h0))
 cev=dx*(MATMUL(hevt,h0))

!Constructing the time dependent fidelity

 DO itime=1, Nt
   !constructing the time dependent coeficients  
   
    eE0=exp(-im*eval0*t(itime)/hbar)
    eEev=exp(-im*evalev*t(itime)/hbar)
   
 

  
   DO i=1, npart
    DO j=1, nstates
    cev_time(j,i)=eEev(j)*cev(j,i)
    c0_time(j,i)=eE0(j)*c0(j,i)
    END DO
   END DO
   
   

  
   !Computing the time dependent wavefunction
   psi_time=MATMUL(h0,c0_time)
   psiev_time=MATMUL(hev, cev_time)

   DO i=1,n
    DO j=1,npart
     psi_timeT(j,i)=CONJG(psi_time(i,j))
    END DO
   END DO

   prodmat=dx*MATMUL(psi_timeT, psiev_time)
  
   DO i=1, npart
     fid(itime, i)=(ZABS(determinant(i,prodmat(1:i, 1:i))))**2
   END DO
 END DO  

END SUBROUTINE fidelity
 






!**************************************************************
!		Ground State Fidelity
!**************************************************************
SUBROUTINE groundfid(n, npart, dx, h0, hev, groundfidel)
IMPLICIT NONE
INTEGER, INTENT(IN)			:: n, npart
REAL(8), INTENT(IN)			:: dx
REAL(8), DIMENSION(npart), INTENT(OUT) 	:: groundfidel
REAL(8), DIMENSION(n,npart), INTENT(IN)	:: h0, hev ! matrix eigenvectors
REAL(8), DIMENSION(npart,n)		:: hevt, h0t 
COMPLEX*16, DIMENSION(npart, npart)	:: cevol
COMPLEX*16, EXTERNAL			:: determinant
REAL(8)					::  f
INTEGER					:: i, j





!Obtaining the trasposed matrix for hev
  DO i=1,n
   DO j=1,npart
    hevt(j,i)=hev(i,j)
    h0t(j,i)=h0(i,j)
   END DO
  END DO

  cevol=dx*(MATMUL(hevt, h0))


   


   DO i=1, npart  
   groundfidel(i)=(ZABS(determinant(i, cevol(1:i, 1:i))))**2   
   END DO
    
  
 


END SUBROUTINE groundfid



!****************************************************************
!		Determinant of a Matrix
!****************************************************************
COMPLEX*16 FUNCTION determinant(N, mat) RESULT(det)
IMPLICIT NONE
INTEGER, INTENT(IN)	:: N
COMPLEX*16, DIMENSION(N,N), INTENT(IN)	:: mat
COMPLEX*16, DIMENSION(N,N)  :: mataux
INTEGER :: i, j, info 
INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
REAL(8)	:: sgn


  DO i=1, N
   DO j=1, N
    mataux(i,j)=mat(i,j)
   END DO
  END DO
  
  ALLOCATE(ipiv(N))

  ipiv= 0

  CALL zgetrf(N, N, mataux, N, ipiv, info)

  det= 1.0

  DO i=1, N
   det=det*mataux(i,i)
  END DO

  sgn= 1.0

  DO i=1, N
    IF (ipiv(i) /= i) THEN
     sgn=-sgn
    END IF
  END DO

  det=sgn*det


  DEALLOCATE(ipiv) 
END FUNCTION determinant


!*************************************************************
!		single particle density subroutine
!**************************************************************
SUBROUTINE onepartdens(n, nstates, npart, u, dens)
IMPLICIT NONE
INTEGER, INTENT(IN)			    ::n, nstates, npart !# of points, # ptcles
REAL(8), DIMENSION(n,nstates), INTENT(IN)   :: u !eigenvectors of the Ham.
INTEGER					    :: i, j
REAL(8), DIMENSION(1,npart)		    :: d
REAL(8), DIMENSION(npart, n)		    :: psi !wave function matrix
REAL(8), DIMENSION(1,n), INTENT(OUT)	    :: dens

!copying elements from H to storage the eigenfunctions for the needed # of particles, taking the traspose

  DO i=1,n
   DO j=1, npart
  psi(j,i)=u(i,j)
   END DO
  END DO
  
!Defining an auxiliary matrix
  d=1.0

!calculating the amplitude density
  DO i=1,n
   DO j=1,npart 
    psi(j,i)=(psi(j,i))**2
   END DO
  END DO

!computing the density(notice: array of 1xn) 
  dens=MATMUL(d, psi)

END SUBROUTINE onepartdens

!**************************************************
!		lattice potential 
!**************************************************
SUBROUTINE potlat(n, xmax, xarr, karr, vpot)
IMPLICIT NONE
INTEGER, INTENT(IN)			:: n
REAL(8), INTENT(IN)			:: xmax
REAL(8), DIMENSION(n), INTENT(IN)	:: xarr, karr
REAL(8), DIMENSION(n,n), INTENT(OUT)	:: vpot
INTEGER					:: i, j, msite  !site number
REAL(8)					:: m, hbar, xunit, tunit, eunit, 						   enrec, pi, v0 

 pi=4.0*ATAN(1.0)
 m=87*1.673e-27
 hbar=1.0546e-34
 xunit=1e-6
 tunit=(2*m*xunit**2)/hbar
 eunit=(hbar**2)/(2*m*(xunit**2))
 msite=240

!Defining some parameters
 enrec=(hbar**2/(2.0*m))*(msite*pi/(2.0*xmax*xunit))**2 !recoil energy
 v0=0.5*enrec/eunit

  DO i=1,n
   DO j=1,n
    IF (i.eq.j) THEN
     vpot(i,j)= v0*(COS((pi*1.0*msite/(2*xmax))*xarr(i)))**2
    ELSE 
     vpot(i,j)=0.0
    END IF
   END DO
  END DO

END SUBROUTINE potlat



!******************************************************
!		potential well
!******************************************************
SUBROUTINE potwell(n,vpot)
IMPLICIT NONE
INTEGER, INTENT(IN)::n
REAL(8), DIMENSION(n,n), INTENT(OUT)::vpot

 vpot=0.0
 vpot(1,1)=100000
 vpot(n,n)=100000


END SUBROUTINE potwell


!**************************************************
!	 Potential(harmonic oscillator) 
!**************************************************
SUBROUTINE pothar(n,xarr,karr,vpot)
IMPLICIT NONE
INTEGER, INTENT(IN)::n
INTEGER:: i, j
REAL(8),DIMENSION(n), INTENT(IN):: xarr, karr
REAL(8):: m, w
REAL(8), DIMENSION(n,n), INTENT(OUT):: vpot

m=1
w=1

  DO i=1,n
   DO j=1,n
    IF (i.eq.j) THEN
     vpot(i,j)=0.5*m*(w**2)*(xarr(i))**2
    ELSE
     vpot(i,j)=0.0 
    END IF
   END DO
  END DO


END SUBROUTINE pothar





