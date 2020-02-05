program beam
   ! use iso_fortran_env, only: rp => REAL64
   use elastsol
   use vtk_legacy
   implicit none

   integer :: i, j, k
   integer,parameter :: nenv = 8 
   ! Adjustable parameters
   logical,      parameter :: steady = .TRUE.
   real(kind=rp),parameter ::                                            &
      L(nsd)    = [0.2_rp, 0.2_rp, 1._rp],                               &
      density   = 1._rp,                                                 &
      stiffness = 100._rp,                                               &
      poisson   = 0.35_rp,                                               &
      dt        = 0.05_rp
   integer, parameter ::                                                 &
      N(nsd)    = [4, 4, 1],                                            &
      nstep     = 1,                                                   &
      iout      = 1

   ! Self evaluating and fixed parameters
   integer,parameter ::                                                  &
      nelem     = PRODUCT(N),                                            &
      nelemv    = 8*nelem,                                               &
      nnod      = PRODUCT(2*N+1),                                        &
      nfnod     = (2*N(2)+1)*(2*N(3)+1),                                 &
      ntip      = nfnod*(2*N(1))+(nfnod+1)/2 
   real(kind=rp),parameter ::                                            &
      dx(nsd)   = L/(2*N)

   integer :: ien(nen,nelem), ienv(nenv,nelemv), dirich_nods(nfnod*nsd), &
              neuman_nods(nfnod), z(2*N(3)+1,2*N(2)+1,2*N(1)+1)
   character(len=11)                  :: fname
   real(kind=rp),dimension(nsd,nnod)  :: x, u, v, a
   real(kind=rp),dimension(nsd*nfnod) :: force

   x   = RESHAPE([(((dx(1)*i,dx(2)*j,dx(3)*k,                            &
            k=0,2*N(3)), j=0,2*N(2)), i=0,2*N(1))],[nsd,nnod])
   
   z   = RESHAPE([(i,i=1,nnod)],2*N(3:1:-1)+1)
   ien = RESHAPE([(((RESHAPE(z(2*k-1:2*k+1,2*j-1:2*j+1,2*i-1:2*i+1),     &
            [nen]),i=1,N(1)),j=1,N(2)),k=1,N(3))],[nen,nelem])

   dirich_nods = [([(i,i=1,           nfnod)],j=1,nsd)]
   neuman_nods = [(  i,i=nnod-nfnod+1,nnod)]

   force       = [([100.0_rp,0._rp,0._rp],i=1,nfnod)]

   ! form ienv for visualization
   do concurrent (i = 1:nelem)
      ienv(:,(i-1)*nenv+1) = [ien( 1,i),ien( 2,i),ien( 5,i),ien( 4,i),   &
                              ien(10,i),ien(11,i),ien(14,i),ien(13,i)]
      ienv(:,(i-1)*nenv+2) = [ien( 2,i),ien( 3,i),ien( 6,i),ien( 5,i),   &
                              ien(11,i),ien(12,i),ien(15,i),ien(14,i)]
      ienv(:,(i-1)*nenv+3) = [ien( 4,i),ien( 5,i),ien( 8,i),ien( 7,i),   &
                              ien(13,i),ien(14,i),ien(17,i),ien(16,i)]
      ienv(:,(i-1)*nenv+4) = [ien( 5,i),ien( 6,i),ien( 9,i),ien( 8,i),   &
                              ien(14,i),ien(15,i),ien(18,i),ien(17,i)]
      ienv(:,(i-1)*nenv+5) = [ien(10,i),ien(11,i),ien(14,i),ien(13,i),   &
                              ien(19,i),ien(20,i),ien(23,i),ien(22,i)]
      ienv(:,(i-1)*nenv+6) = [ien(11,i),ien(12,i),ien(15,i),ien(14,i),   &
                              ien(20,i),ien(21,i),ien(24,i),ien(23,i)]
      ienv(:,(i-1)*nenv+7) = [ien(13,i),ien(14,i),ien(17,i),ien(16,i),   &
                              ien(22,i),ien(23,i),ien(26,i),ien(25,i)]
      ienv(:,(i-1)*nenv+8) = [ien(14,i),ien(15,i),ien(18,i),ien(17,i),   &
                              ien(23,i),ien(24,i),ien(27,i),ien(26,i)]
   end do

   call init_elastsol(density,x,ien,dirich_nods,[nfnod,nfnod,nfnod],     &
                      nnod,nelem)

   u = x; v = 0._rp; a = 0._rp
   open(10,file='tip_disp.dat')
   !call deform(u,v,a,force,neuman_nods,nfnod,stiffness,poisson,dt)
   !force = 0._rp
   do i=1,nstep
      call deform(u,v,a,force,neuman_nods,nfnod,stiffness,poisson,       &
                  dt,steady)
      write(10,*) i*dt,u(:,ntip)-x(:,ntip)
      if(mod(i,iout)==0) then
         write(fname,'(a4i3.3a4)') 'beam',i/iout,'.vtk'
         call init_vtk(fname,'elastsol_beamtest')
         call vtkwrite_mesh(x,ienv,'HEXA')
         call vtkwrite_vector(u-x,'displacement')
         call exit_vtk
      end if
      write(logfile,*) 'Completed time step:',i
   end do
   call exit_elastsol
   close(10)
end program beam
