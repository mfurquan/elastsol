program bar
   ! use iso_fortran_env, only: rp => REAL64
   use elastsol
   use vtk_legacy
   implicit none

   integer :: i, j, k
   integer,parameter :: nenv = 8 
   ! Adjustable parameters
   logical,      parameter :: steady = .TRUE.
   real(kind=rp),parameter ::                                            &
      L(nsd)    = +2._rp*[1._rp, 1._rp, 20._rp],                         &
      x0(nsd)   = -1._rp*[1._rp, 1._rp, 1._rp],                          &
      density   = 1._rp,                                                 &
      stiffness = 100._rp,                                               &
      poisson   = 0.35_rp,                                               &
      dt        = 0.05_rp
   integer, parameter ::                                                 &
      N(nsd)    = [4, 4, 10],                                            &
      nstep     = 1,                                                     &
      iout      = 1

   ! Self evaluating and fixed parameters
   integer,parameter ::                                                  &
      nelem     = PRODUCT(N),                                            &
      nelemv    = 8*nelem,                                               &
      nnod      = PRODUCT(2*N+1),                                        &
      nfnod     = (2*N(1)+1)*(2*N(2)+1),                                 &
      ntip      = nfnod*(2*N(3))+(nfnod+1)/2 
   real(kind=rp),parameter ::                                            &
      dx(nsd)   = L/(2*N)

   integer :: ien(nen,nelem), ienv(nenv,nelemv), dirich_nods(nfnod*nsd), &
              neuman_nods(nfnod), z(2*N(1)+1,2*N(2)+1,2*N(3)+1)
   character(len=11)                  :: fname
   real(kind=rp),dimension(nsd,nnod)  :: x, u, v, a
   real(kind=rp),dimension(nsd,nfnod) :: force

   x   = RESHAPE([(((x0+dx*[i,j,k],                                      &
           i=0,2*N(1)), j=0,2*N(2)), k=0,2*N(3))],[nsd,nnod])
   
   z   = RESHAPE([(i,i=1,nnod)],2*N+1)
   ien = RESHAPE([(((RESHAPE(z(2*i-1:2*i+1,2*j-1:2*j+1,2*k-1:2*k+1),     &
            [nen]),i=1,N(1)),j=1,N(2)),k=1,N(3))],[nen,nelem])

   dirich_nods = [([(i,i=1,           nfnod)],j=1,nsd)]
   neuman_nods = [(  i,i=nnod-nfnod+1,nnod )]

   force = SPREAD([0._rp,0._rp,10._rp],2,nfnod)                          &
         / nfnod!* SPREAD([1,2,1,2,4,2,1,2,1]/16._rp,1,nsd)

   write(*,*) 'FL/EA=',SUM(force)*L(3)/(L(1)*L(2)*stiffness)
   !write(*,*) 'FL/GA=',SUM(force)*L(3)*2._rp*(1._rp+poisson)             &
   !                   /(L(1)*L(2)*stiffness)

   ! form ienv for visualization
   ienv = split_HEXA(ien)

   call init_elastsol(density,x,ien,dirich_nods,[nfnod,nfnod,nfnod],     &
                      nnod,nelem)

   u = 0._rp; v = 0._rp; a = 0._rp
   open(10,file='tip_disp.dat')
   !call deform(u,force,neuman_nods,nfnod,stiffness,poisson,v,a,dt)
   !force = 0._rp
   do i=1,nstep
      call deform(u,force,neuman_nods,nfnod,stiffness,poisson)!,v,a,dt)
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
end program bar
