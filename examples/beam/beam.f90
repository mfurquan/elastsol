program beam
   ! use iso_fortran_env, only: rp => REAL64
   use elastsol
   use vtk_legacy
   implicit none

   integer :: i, j, k
   integer,parameter :: nenv = 8 
   ! Adjustable parameters
   real(kind=rp),parameter ::                                            &
      L(nsd)    = +1._rp*[0.02_rp, 0.02_rp, 0.35_rp],                    &
      density   = 1.E3_rp,                                               &
      stiffness = 1.4E6_rp,                                              &
      poisson   = 0.4_rp,                                                &
      g         = 2.0_rp,                                                &
      dt        = 0.01_rp
   integer, parameter ::                                                 &
      N(nsd)    = [4, 1, 50],                                            &
      nstep     = 300,                                                   &
      iout      = 5

   ! Self evaluating and fixed parameters
   integer,parameter ::                                                  &
      nelem     = PRODUCT(N),                                            &
      nelemv    = 8*nelem,                                               &
      nnod      = PRODUCT(2*N+1),                                        &
      nfnod     = (2*N(1)+1)*(2*N(2)+1),                                 &
      ntip      = nfnod*(2*N(3))+(nfnod+1)/2 
   real(kind=rp),parameter ::                                            &
      dx(nsd)   = L/(2*N),                                               &
      pi        = acos(-1._rp)

   integer :: ien(nen,nelem), ienv(nenv,nelemv), dirich_nods(nfnod*nsd), &
              neuman_nods(nnod), z(2*N(1)+1,2*N(2)+1,2*N(3)+1)
   character(len=11)                  :: fname
   real(kind=rp),dimension(nsd,nnod)  :: x, u, v, a
   real(kind=rp) :: force(nsd,nnod), EI, rhoA

   x   = RESHAPE([(((dx*[i,j,k], i=0,2*N(1)), j=0,2*N(2)), k=0,2*N(3))], &
                 [nsd,nnod])
   
   z   = RESHAPE([(i,i=1,nnod)],2*N+1)
   ien = RESHAPE([(((RESHAPE(z(2*i-1:2*i+1,2*j-1:2*j+1,2*k-1:2*k+1),     &
            [nen]),i=1,N(1)),j=1,N(2)),k=1,N(3))],[nen,nelem])

   dirich_nods = [([(i,i=1,           nfnod)],j=1,nsd)]
   neuman_nods = [(  i,i=1,nnod)]

   EI    = stiffness*L(1)**3*L(2)/12._rp
   rhoA  = density*L(1)*L(2)
   force = g*rhoA*L(3)*SPREAD([1._rp,0._rp,0._rp],2,nnod)/nnod
   write(*,*) 'fL^4/8EI =',SUM(force)*L(3)**3/(8._rp*EI)
   write(*,*) 'f_n1 =',1.875_rp**2*sqrt(EI/rhoA)/(2._rp*pi*L(3)**2)

   ! form ienv for visualization
   ienv = split_HEXA(ien)

   call init_elastsol(density,x,ien,dirich_nods,[nfnod,nfnod,nfnod],     &
                      nnod,nelem)

   u = 0._rp
   call deform(u,force,neuman_nods,nnod,stiffness,poisson)
   write(*,*) 'steady tip displacement =',u(2,ntip)
   !force = 0._rp

   u = 0._rp; v = 0._rp; a = 0._rp
   open(10,file='tip_disp.dat')
   do i=1,nstep
      call deform(u,force,neuman_nods,nnod,stiffness,poisson,v,a,dt)
      write(10,*) i*dt,u(:,ntip)
      if(mod(i,iout)==0) then
         write(fname,'(a4i3.3a4)') 'beam',i/iout,'.vtk'
         call init_vtk(fname,'elastsol_beamtest')
         call vtkwrite_mesh(x,ienv,'HEXA')
         call vtkwrite_vector(u,'displacement')
         call exit_vtk
      end if
      write(logfile,*) 'Completed time step:',i
   end do
   call exit_elastsol
   close(10)
end program beam
