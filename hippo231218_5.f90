module parameter
INTEGER,parameter :: n_t=200000, m_t=248, nlit=n_t/m_t, mlit=n_t/5000  !total num. of particle,timestep n_t=33000
INTEGER,parameter :: INTEGMX=2147483647,INTEGST=584287,INTEG=48828125 ! 2^31-1, x0, 5^11 Random number seeds 
INTEGER,parameter :: ldim=3, mdim=2 ! Dimensions of the cellular environment and the mapping space
double precision, PARAMETER :: pi=4.D0*datan(1.D0) ! the ratio of the circumference of a circle to its diameter
double precision, PARAMETER :: cbrt2=2.D0**(1.D0/3.D0) ! 
double precision, PARAMETER :: dt=3.75D0, & ! time step [sec]
m0 = 2.0D-12,   & ! mass of the initial cell [kg] 
k_fric  = 5.D-13, & !  frictional coefficient [kg/sec]
epsilon = 2.D-15, &  !energy [kg um^2/sec^2]
al_LJ=0.3D0, lambda=0.2D0, &  !energy [kg um^2/sec^2]
gamma = cbrt2 - 1.D0  ! initial arrangement ratio of the M period (trial 0<gamma<1 )
double precision, PARAMETER :: twop1r6=(2.D0 - al_LJ*(1.D0 - lambda)**(2.D0) )**( 1.D0/6.D0) ! softcore correction of twop1r6=(2.D0)**(1.D0/6.D0),
double precision, PARAMETER :: twom1r6=(2.D0 - al_LJ*(1.D0 - lambda)**(2.D0) )**(-1.D0/6.D0) ! softcore correction of twom1r6=(2.D0)**(-1.D0/6.D0)
double precision, PARAMETER :: sgm0d=5.0D0, sgm0 = sgm0d*twom1r6 !  sgm0d and sgm0 [um]
! Reaction constants
double precision, PARAMETER :: a1=5.D-2, a2=5.D-4, b1=1.D0, a3=1.D-4 ! reaction rates
double precision, PARAMETER :: ttc=20.D0*3.6D3, dsgm=(cbrt2 - 1.D0)*dt*sgm0/ttc ! period and growth rate of cells 
double precision, PARAMETER :: ttm=1.D0*3.6D3 ! M-periods [sec] 
INTEGER,parameter :: nnc=int(ttc/dt), mmt=int(dble(n_t)/dble(nnc)), nnt=dint(2.D0**(dble(mmt+1))) ! total cell count
INTEGER,parameter :: ndim=3 ! Dimensions of reaction(ndim) and rec_* 
double precision, PARAMETER :: Xth = 0.012D0 ! a threshold level of the cell growth or cellcycle arrest  [uM]
double precision Xinit(ndim) ; data Xinit/0.1D0, 0.4D0, 0.D0/ !initial concentration
INTEGER,parameter :: nrho=4 ! number of points for density rho
double precision rhon(nrho) ; data rhon/1.D-90,0.5D0, 1.D0,1.5D0/ 
double precision d_am(0:3) ; data d_am/2.0D0, 3.5D0, 1.0D0, 4.0D0/ !colors of cyan, orange, red, blue  
complex(kind(0d0)) :: ci=dcmplx(1.D0, 0.0D0), cj=dcmplx(0.D0, 1.D0), c0=dcmplx(0.D0, 0.D0)
end module parameter
!
program HIPPO
    call sHIPPO   
end program HIPPO
!
!
subroutine sHIPPO
    USE parameter
    INTEGER :: i, j, k, l, n, i_t, j_t, nn(n_t), ix, id_am(nnt,m_t), cm(nnt), cmtemp(nnt)
    double precision :: t, r(nnt,ldim),v(nnt,ldim),rb(ldim),rtemp(nnt,ldim),vtemp(nnt,ldim)
    double precision :: rmat(nnt, m_t, ldim), sgmmm(nnt, m_t), sgm(nnt), mc(nnt,ldim)
    double precision :: sij, s_inv,ss_inv, ss_inv2, ss_inv3, r_cut,u_LJ,u_inv
    double precision :: f(nnt,ldim), g(nnt,ldim), Xd(ndim,nnt),al(nnt)
    double precision :: r2,u_rn(2),theta,phi, rij
    double precision :: X0, X( ndim,nnt), Xb( ndim,nnt), rhom(3,0:n_t)
    double precision :: Xtemp( ndim,nnt) 
    double precision :: X1(nnt), X2(nnt), X3(nnt), surv_ratio, Xm(ndim,3,0:n_t)
    double precision :: Xd1(ndim,nnt), Xd2(ndim,nnt), Xd3(ndim,nnt), Xd4(ndim,nnt),r_i(ldim),r_j(ldim)
    double precision :: Y( ndim), Yd(ndim),Yrho(nrho,ndim,0:n_t)
    double precision :: Yd1(ndim), Yd2(ndim), Yd3(ndim), Yd4(ndim)
    double precision :: costh, sinth, cosph, sinph, e_r(ldim), rho(nnt), sum
    character filename*128
!
    OPEN(UNIT=10, file='surat_t231218_5.dat',status='replace')
    OPEN(UNIT=20, file='Xyzw_t231218_5.dat',status='replace')
!    
    write(6,*)'m_t=',m_tnlitdat
    write(6,*)'epsilon[J],k_fric[kg/s]=',epsilon*1.D-12,k_fric
    nn = 2 ! setting initial number of cells   
  ! Initial positions and velocities of i-th cells 
    ix=0 
    call suniform_rn(2, ix, u_rn) 
    theta=2.D0*pi*u_rn(1) ! theta: uniform [0, 2 pi) random number 
    phi  =     pi*u_rn(2) ! phi:   uniform [0,   pi) random number 
    sgm = sgm0 ; mc = m0 ! setting the initial size and mass of cells
    costh=dcos(theta) ; sinth=dsin(theta) ; cosph=dcos(phi) ; sinph=dsin(phi)
    e_r(1) = +sinph*costh ; e_r(2) = +sinph*sinth ; e_r(3) =  cosph  
    r(1,:)=-gamma*sgm0d*e_r(:) ; v(1,:)=-(1.D0-gamma)*sgm0d*e_r(:)/ttm
    r(2,:)=+gamma*sgm0d*e_r(:) ; v(2,:)=+(1.D0-gamma)*sgm0d*e_r(:)/ttm       
    al = 0.D0  ! initializing ages of cells
  ! Reaction initial values 
  ! Reference from the material of HIPPO path mathematical model 
    Xm = 0.D0 ; X  = 0.D0 ; X0 =0.15D0 
     do i = 1, nn(1) ; X(:, i) = Xinit(:) ;  enddo    !initiallizing X [uM]
    X1 = 0.D0 ; X2 = 0.D0 ; X3 = 0.D0    !!initiallizing X1, X2, X3
    Xm(:, 3, 0)= X(:,1) 
    do l=1,ndim ; Xm(l,:,0)=X(l,1) ; enddo  ! initiallizing output Xm at t=0
    rhom = 0.D0 ! initiallizing packing fraction 
    cm=1 ;id_am=1
!
    DO i_t = 1, n_t ; j_t = i_t / nlit ; t = i_t * dt ; if(i_t .ge. 2)nn(i_t)=nn(i_t-1)
! (2)    
   DO i = 1, nn(i_t) 
    if(X(2, i) .gt. Xth)then ; al(i) = al(i) + dt ; sgm(i) = sgm(i) + dsgm ; mc(i,:)=m0*(sgm(i)/sgm0)**(3.D0) ; endif ! aging, growth radius and mass [kg] of cells
    if(X(2, i) .gt. Xth .and. al(i) .le. ttm)cm(i) = 1  !  orange; Cell growth and M  phase 
    if(X(2, i) .gt. Xth .and. al(i) .gt. ttm)cm(i) = 3  !  red;    Cell growth and G1 phase 
    if(X(2, i) .le. Xth .and. al(i) .le. ttm)cm(i) = 0  !  cyan;   Cellcycle arrest and M  phase 
    if(X(2, i) .le. Xth .and. al(i) .gt. ttm)cm(i) = 2  !  blue;   Cellcycle arrest and G1 phase
    if(mod(i_t, nlit) .eq. 0)id_am(i, j_t) = cm(i)  !  for output color index 
   enddo !do i=1,nn(i_t)  
! (3)   
   Xtemp = X ; rtemp = r ; vtemp = v ; cmtemp = cm
   DO i = 1, nn(i_t)
    if(sgm(i) .ge. cbrt2*sgm0 )then 
       nn(i_t) = nn(i_t) + 1 ;  n = nn(i_t) ! Proliferation of Cell division  
       sgm(i) = sgm0 ; sgm(n) = sgm0
       al( i) = 0.D0 ; al( n) = 0.D0
       call suniform_rn(2, ix, u_rn) 
       theta=2.D0*pi*u_rn(1) ! theta: uniform [0, 2 pi) random number 
       phi  =     pi*u_rn(2) ! phi:   uniform [0,   pi) random number 
       costh=dcos(theta) ; sinth=dsin(theta) ; cosph=dcos(phi) ; sinph=dsin(phi)
       e_r(1) = +sinph*costh ; e_r(2) = +sinph*sinth ; e_r(3) =  cosph  
       r(i,:)=rtemp(i,:)-gamma*sgm0d*e_r(:) ; v(i,:)=-(1.D0-gamma)*sgm0d*e_r(:)/ttm
       r(n,:)=rtemp(i,:)+gamma*sgm0d*e_r(:) ; v(n,:)= (1.D0-gamma)*sgm0d*e_r(:)/ttm
       X(:,i)=Xtemp(:,i) ; X(:,n)=Xtemp(:,i) 
       cm(i) = 1 ; cm(n) = 1
       if(mod(i_t, nlit) .eq. 0)id_am(i, j_t) = cm(i)
       if(mod(i_t, nlit) .eq. 0)id_am(n, j_t) = cm(n)
    endif ! Cell position, velocity, age, concentration X and color index cm updates
   enddo !do i=1,nn(i_t) 
   n = nn(i_t)
  ! (4)
   f = 0.D0  !    Lennard-Jones potential force
   do i = 1, n
    do j = 1, n ; if(j .eq. i)goto 100  
      rb(:) = r(i,:) - r(j,:) ! vector rb_l of i-j-th intercellular distance  
      r2  =0.D0 ; do l = 1, ldim ; r2 = r2  + rb(l)**(2.D0) ; enddo ; rij = dsqrt(r2)
      sij = (sgm(i) + sgm(j)) ; r_cut = 3.D0 * sij   
      if(rij .eq. 0.D0 .or. rij .gt. r_cut )goto 100
       s_inv = 1.D0 / sij ; ss_inv=s_inv**(2.D0) ; ss_inv2=ss_inv**(2.D0) ; ss_inv3=ss_inv2*ss_inv
       u_LJ = al_LJ*(1.D0 - lambda)**(2.D0) + r2**(3.D0) * ss_inv3 ; u_inv = 1.D0 / u_LJ
       f(i,:)=f(i,:)+(rb(:)/rij)*lambda*24.D0*epsilon*r2**(2.D0)*rij*ss_inv3*u_inv*(2.D0*u_inv**(2.D0)-u_inv)
  100  continue
    enddo !do j = 1, n
   enddo !do i = 1, n    
! id_amc
  !    Frictional force of cells
     do i = 1, n ; g(i,:) = -k_fric * v(i,:)  ; enddo   
  ! 
  !   cell positions and velocities evolved by Verlet method - 1/2
      v  = v + 5.D-1 * dt * ( f + g ) / mc
      r    = r + dt * v 
  !
  !    Lennard-Jones potential force
      f = 0.D0 
        do i = 1, n
          do j = 1, n ; if(j .eq. i)goto 300  
            rb(:) = r(i,:) - r(j,:) ! vector rb_l of i-j-th intercellular distance  
            r2  =0.D0 ; do l = 1, ldim ; r2   = r2  + rb(l)**(2.D0) ; enddo ; rij = dsqrt(r2)
            sij = (sgm(i) + sgm(j)) ; r_cut = 3.D0 * sij   
            if(rij .eq. 0.D0 .or. rij .gt. r_cut )goto 300
             s_inv = 1.D0 / sij ; ss_inv=s_inv**(2.D0) ; ss_inv2=ss_inv**(2.D0) ; ss_inv3=ss_inv2*ss_inv
             u_LJ = al_LJ*(1.D0 - lambda)**(2.D0) + r2**(3.D0) * ss_inv3 ; u_inv = 1.D0 / u_LJ 
            f(i,:)=f(i,:)+(rb(:)/rij)*lambda*24.D0*epsilon*r2**(2.D0)*rij*ss_inv3*u_inv*(2.D0*u_inv**(2.D0)-u_inv) 
        300  continue
          enddo !do j = 1, n
        enddo !do i = 1, n    
        ! id_amc
          !    Frictional force of cells
             do i = 1, n ; g(i,:) = -k_fric * v(i,:) ;enddo   
  !
  !   cell positions and velocities evolved by Verlet method - 2/2
        v  = v + 5.D-1 * dt * ( f + g ) / mc 
  !       
  !        
      if(mod(i_t, nlit) .eq. 0)then     
        do i=1,n
          do j=1,ldim  
            rmat(i, j_t, j) = r(i,j)
          enddo !j=1,ldim
          sgmmm(i, j_t)  = sgm(i) 
        enddo !  i=1,n        
      endif
! 
!   reaction equations solved by Runge Kutta method   
!   
!
!  packing fraction
!             
    rho= 0.D0 
    DO i = 1, n
     sum = 0.D0     
     do j = 1, n  ; if(j .eq. i)goto 400  
      rb(:) = r(i,:) - r(j,:) 
      rij = 0.D0 ; do l = 1, ldim ; rij = rij + rb(l)**(2.D0) ; enddo ; rij = dsqrt(rij) 
      sum = sum + dexp(2.D0-rij/sgm0d)/12.D0  !dexp(2.D0)/12.D0
400   continue 
     enddo ! j = 1, n
     rho(i) = sum       
    enddo ! i = 1, n
!
!     Solving reaction equations by Eigen value problem
!        
            Xb = X ; call sEVP(rho, Xb, X, X0, n)
!            
            DO j = 1, n       
              X1(j)  = X(1, j) ; X2(j)  = X(2, j) ; X3(j)  = X(3, j)
            ENDDO !  j = 1, n
!               
            Xm(1,1,i_t) = MAXVAL(X1)   ; Xm(1,2,i_t) = MINVAL(  X1, MASK=X1   .GT. 0.0)
            Xm(2,1,i_t) = MAXVAL(X2)   ; Xm(2,2,i_t) = MINVAL(  X2, MASK=X2   .GT. 0.0)
            Xm(3,1,i_t) = MAXVAL(X3)   ; Xm(3,2,i_t) = MINVAL(  X3, MASK=X3   .GT. 0.0)           
            DO l = 1, ndim 
              sum=0.D0 ; DO i = 1, n ; sum = sum + X(l, i) ; ENDDO
              Xm(l, 3, i_t)= sum/dfloat(n) 
            ENDDO !  l = 1, ndim
!
            rhom(1,i_t) = MAXVAL(rho) ; rhom(2,i_t) = MINVAL(rho, MASK=rho .GT. 0.0)
            sum=0.D0 ; DO i = 1, n ; sum = sum + rho(i) ; ENDDO !  i = 1, n
            rhom(3, i_t) = sum/dfloat(n)   ! cell-mean of rhom
!        
    END DO ! i_t = 1, n_t
!
!   X_rho=0 and X_rho=1 
!        
    Y  = 0.D0 ; Y(1)=Xinit(1) ; Y(2)=Xinit(2) ; Y(3)=Xinit(3)    !initiallizing Y [uM]
    do l = 1, nrho ; Yrho(l,:,0) = Y(:) ; enddo
!
    DO j = 1, nrho   !; write(6,*)'rho=',rho
        Y = Xinit
        DO i_t = 1, n_t ; j_t = i_t / nlit !; if(mod(i_t, nlit) .eq. 0)write(6,*)'j_t=',j_t
        t     = i_t * dt 
! 
!     Solving reaction equations by Eigen value problem
!   
        rho(1)=rhon(j) ; Xb(:,1)=Y(:) ; call sEVP(rho, Xb, X, X0, 1) ; Y(:) = X(:,1)
   !
            Yrho(j,:,i_t) = Y(:)
       END DO ! i_t = 1, n_t 
    Enddo ! j = 1, nrho
!
!   outputs
    DO j_t = 1, m_t ; i_t = j_t * nlit
        t     = i_t * dt ; t = t / 3.6D3 ! [sec] to [hour]
        write (filename, '("xyz", i4.4, ".csv")') j_t
        open (17, file=filename, status='replace')
        write(17,'(A)')'Data Format, 23'
        write(17,'(A)')'memo1'
        write(17,'(A)')'x, X, X, R, data'
        surv_ratio = 0.D0
        DO i = 1, nn(i_t) 
          write(17,'(E13.6,A1,E13.6,A1,E13.6,A1,E13.6,A1,F6.2)')rmat(i, j_t, 1),",",rmat(i, j_t, 2),",",rmat(i, j_t, 3), &
                                                              ",",sgmmm(i, j_t)*twop1r6,",", d_am(id_am(i, j_t))
          if(d_am(id_am(i, j_t)) .ge. 3.D0)surv_ratio = surv_ratio + 1.D0 
        END DO !  i = 1, nn(i_t)
        surv_ratio = surv_ratio/dfloat(nn(i_t))
        write(10,'(1000(E13.6,2X))')t, dfloat(nn(i_t)), surv_ratio
      END DO ! j_t = 1, m_t
!!!
      DO i_t = 0, n_t !, mlit 
        t     = i_t * dt 
        if(dmod(t, ttc) .lt. 360.D0  .or. mod(i_t, mlit) .eq. 0)then 
        write(20,'(1000(E13.6,2X))')t / 3.6D3, (Xm(i,1,i_t), Xm(i,3,i_t), Xm(i,2,i_t), i=1, ndim), & ! t [sec] to [hour]
                                      Xth, ((Yrho(j,l,i_t), j=1,nrho), l=1,ndim),(rhom(i,i_t), i=1,3)
        endif    
       END DO !  i_t=0, n_t 
      CLOSE(10) ; CLOSE(20)
end subroutine sHIPPO
!
subroutine sEVP(rho, Xb, X, X0, n)
    USE parameter
    INTEGER :: i, j
    double precision :: b, c, Xb(ndim,nnt), X(ndim,nnt) 
    complex(kind(0d0)) :: cX(ndim), Y(ndim), e(ndim,ndim), p(ndim,ndim), q(ndim,ndim)
    complex(kind(0d0)) :: Yb(ndim),cXb(ndim),lam(ndim)
    double precision :: X0, rho(nnt)
    character filename*128
!
     do i = 1, n ; cXb(:) = ci * Xb(:, i) ; if(dabs(rho(i)) .le. 1.D-90)rho(i) = 1.D-90
      b   = a1 + a2 + a3 + b1 * rho(i) * X0 ; c = (a1+a2)*a3 + a2*b1*rho(i)*X0
      lam(1) = 0.5D0*(-b*ci + cdsqrt((b**(2.D0) - 4.D0*c)*ci))  ! Eigenfrequency 1 
      lam(2) = 0.5D0*(-b*ci - cdsqrt((b**(2.D0) - 4.D0*c)*ci))  ! Eigenfrequency 2
      lam(3) = c0                                          ! Eigenfrequency 3 (0 rad/s)
!     
      p(1, 1)=lam(1)+ci*a3; p(2,1)=a1*(lam(1)+ci*a3)/(lam(1)+ci*a2) ; p(3,1)=ci*b1*rho(i)*X0
      p(1, 2)=lam(2)+ci*a3; p(2,2)=a1*(lam(2)+ci*a3)/(lam(2)+ci*a2) ; p(3,2)=ci*b1*rho(i)*X0
      p(1, 3)=lam(3)+ci*a3; p(2,3)=a1*(lam(3)+ci*a3)/(lam(3)+ci*a2) ; p(3,3)=ci*b1*rho(i)*X0
!  
      call sinvA(p,q)
      call sAv(q, cXb, Yb)
!
      do j = 1, ndim ; Y(j) = cdexp(lam(j)*dt)*Yb(j) ; enddo !j = 1, ndim 
      call sAv(p, Y, cX) 
      X(:, i) = dreal(cX(:)) 
     enddo ! i = 1, n 
    end subroutine sEVP
!
! X v -> Y
        subroutine sAv(X,v,w)
            USE parameter
            IMPLICIT NONE 
            INTEGER :: i, j, k
            complex(kind(0d0)) :: X(ndim,ndim),v(ndim),w(ndim) 
        !
            w = c0
            do i=1,ndim         
              do k=1,ndim ; w(i) = w(i) + X(i,k)*v(k) ; enddo
            enddo 
        end subroutine sAv  
! X^-1 -> Y
        subroutine sinvA(X,Y)
            USE parameter
            IMPLICIT NONE 
            complex(kind(0d0)) :: X(ndim,ndim),Y(ndim,ndim),detA 
        !      
         detA=X(1,1)*X(2,2)*X(3,3)+X(1,2)*X(2,3)*X(3,1)+X(1,3)*X(2,1)*X(3,2) &
             -X(1,3)*X(2,2)*X(3,1)-X(1,1)*X(2,3)*X(3,2)-X(1,2)*X(2,1)*X(3,3)
         Y(1,1)= X(2,2)*X(3,3)-X(2,3)*X(3,2);Y(1,2)=-X(1,2)*X(3,3)+X(1,3)*X(3,2);Y(1,3)= X(1,2)*X(2,3)-X(1,3)*X(2,2)
         Y(2,1)=-X(2,1)*X(3,3)+X(2,3)*X(3,1);Y(2,2)= X(1,1)*X(3,3)-X(1,3)*X(3,1);Y(2,3)=-X(1,1)*X(2,3)+X(1,3)*X(2,1)
         Y(3,1)= X(2,1)*X(3,2)-X(2,2)*X(3,1);Y(3,2)=-X(1,1)*X(3,2)+X(1,2)*X(3,1);Y(3,3)= X(1,1)*X(2,2)-X(1,2)*X(2,1)
         Y = Y/detA     
        end subroutine sinvA              
!    
!> pseudo-randomized and uniform [0,1) random number by Linear congruential method
  SUBROUTINE suniform_rn(ni, ix, res)
    USE parameter
    IMPLICIT NONE
    double precision :: res(ni), AINTEGMX 
    integer i, ni, ix
!
    AINTEGMX = REAL ( INTEGMX )
! 
   ! IF ( ix .LT. 0 ) PAUSE
    IF ( ix .EQ. 0 ) ix = INTEGST
    DO  i = 1, ni 
        ix = ix * INTEG
        IF (ix) 10, 20, 20
10      ix   = (ix + INTEGMX) + 1
20      res(I) = REAL(ix)/AINTEGMX 
    ENDDO
  END SUBROUTINE suniform_rn
