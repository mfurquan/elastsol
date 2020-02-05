module hexa_lagr
! LAGRANGIAN SHAPE FUNCTIONS FOR HEXAHEDRA
!======================================================================
   !use iso_fortran_env, only: rp => REAL64
   implicit none


!            A D J U S T A B L E   P A R A M E T E R S
!----------------------------------------------------------------------
   integer,parameter ::     &
      rp   = 8,             & ! real precision
      deg  = 2,             & ! degree of interpolation
      nqd1 = 3                ! no. of gauss-quadrature points / dimension

!     F I X E D  ( H A R D - C O D E D )  P A R A M E T E R S
!----------------------------------------------------------------------
   integer,parameter :: nsd = 3, ntd = 6

!  1D Gaussian quadrature
   integer,parameter,private :: maxint = 6
   real(kind=rp),parameter,private ::               &
      wg(maxint) = [                                & ! weights
         2._rp,                                     & ! 1 point
         1._rp, 1._rp,                              & ! 2 point
         5._rp/9._rp,8._rp/9._rp,5._rp/9._rp],      & ! 3 point
!         
      xg(maxint) = [                                & ! coordinates
         0._rp,                                     & ! 1 point
         -1._rp/sqrt(3._rp),1._rp/sqrt(3._rp),      & ! 2 point
         -sqrt(3._rp/5._rp),0._rp,sqrt(3._rp/5._rp)]  ! 3 point

!      S E L F - E V A L U A T I N G    P A R A M E T E R S
!----------------------------------------------------------------------
   integer,parameter ::                                                  &
      qi   = (nqd1 - 1)*nqd1/2 + 1,                                      &
      qf   = (nqd1 + 1)*nqd1/2,                                          &
      nen1 =   deg + 1,                                                  &
      nen  = nen1**nsd,                                                  &
      nqd  = nqd1**nsd


   integer,private :: i, j, k, p, q, r

   logical,parameter,private,dimension(nqd1,nen1,nen1) ::                &
      mq = spread(reshape([((i/=j, i = 1,nen1), j = 1,nen1)],            &
                            [nen1,nen1]),1,nqd1)

   logical,parameter,private,dimension(nqd1,nen1,nen1,nen1-1) ::         &
      cq = spread(reshape([(((i/=j.AND.(i-j/=k.AND.i-j/=k-nen1),         &
                              i = 1,nen1), j = 1,nen1), k = 1,nen1-1)],  &
                          [nen1,nen1,nen1-1]),1,nqd1)
   
   real(kind=rp),parameter,private,dimension(nen1) ::                 &
      xm    = [(-1+2*i/deg, i=0,deg)]                 

   real(kind=rp),parameter,private,dimension(nqd1,nen1,nen1) ::       &
      xq    = spread(spread(xg(qi:qf),2,nen1),3,nen1),                &
      xa    = spread(spread(xm,1,nen1),1,nqd1),                       &
      xb    = spread(spread(xm,1,nqd1),3,nen1)        


!   1 D  S H A P E  F U N C T I O N S  A N D   D E R I V A T I V E S
!----------------------------------------------------------------------
   real(kind=rp),parameter,private,dimension(nqd1,nen1) ::            &
      phi1    = product(xq - xa,3,mq)/product(xb - xa,3,mq),          &
      phi1_xi = sum(product(spread(xq - xa,4,nen1-1),3,cq),3)         &
              / product(xb - xa,3,mq)

!  3 D   G A U S S   Q U A D R A T U R E   W E I G H T S 
!----------------------------------------------------------------------
   real(kind=rp),parameter,dimension(nqd) ::                          &
      wq    = [(((wg(qi+i)*wg(qi+j)*wg(qi+k),                         &
                  i = 0,nqd1-1), j = 0,nqd1-1), k = 0,nqd1-1)]

!  3 D   S H A P E   F U N C T I O N S   A N D   D E R I V A T I V E S
!----------------------------------------------------------------------
   real(kind=rp),parameter,dimension(nqd,nen) ::                      &
      phi     = reshape([((((((phi1(i,p)*phi1(j,q)*phi1(k,r),         &
                        i = 1,nqd1), j = 1,nqd1), k = 1,nqd1),        &
                        p = 1,nen1), q = 1,nen1), r = 1,nen1)],       &
                        [nqd,nen])

   real(kind=rp),parameter,dimension(nsd,nen,nqd) ::                  &
      phi_xi = reshape([((((((phi1_xi(i,p)*phi1(j,q)*phi1(k,r),       &
                              phi1(i,p)*phi1_xi(j,q)*phi1(k,r),       &
                              phi1(i,p)*phi1(j,q)*phi1_xi(k,r),       &
                        p = 1,nen1), q = 1,nen1), r = 1,nen1),        &
                        i = 1,nqd1), j = 1,nqd1), k = 1,nqd1)],       &
                        [nsd,nen,nqd])
end module hexa_lagr
