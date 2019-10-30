program beam
   use iso_fortran_env, only: rp => REAL64
   use elastsol
   use vtk_legacy
   implicit none

   integer,parameter :: nsd = 3, nen = 27
   ! Adjustable parameters
   real(kind=rp),parameter ::                                            &
      L(nsd)    = [3.5_rp, 0.2_rp, 1._rp],                               &
      density   = 1._rp,                                                 &
      stiffness = 100._rp,                                               &
      poisson   = 0.35_rp,                                               &
      dt        = 0.05_rp
   integer, parameter ::                                                 &
      N(nsd)    = [10, 4, 1],                                            &
      nstep     = 100,                                                   &
      iout      = 10

   ! Self evaluating and fixed parameters
   integer,parameter :: nelem = PRODUCT(N), nnod  = PRODUCT(2*N+1),      &
                        nfnod = 3*(2*N(2)+1),                            &
                        ntip = nfnod*(2*N(1))+(nfnod+1)/2 
   real,   parameter :: dx(nsd) = L/(2*N)


   integer                           :: ien(nen,nelem),                  &
                                        a(2*N(3)+1,2*N(2)+1,2*N(1)+1),   &
                                        i, j, k, dirich_nods(nfnod*nsd)
   real(kind=rp),dimension(nsd,nnod) :: x, u, v, a
   character(len=11)                 :: fname

   x   = RESHAPE([(((dx(3)*i,dx(2)*j,dx(1)*k,                            &
            k=0,2*N(3)), j=0,2*N(2)), i=0,2*N(1))],[nsd,nnod])
   
   a   = RESHAPE([(i,i=1,nnod)],2*N(3:1:-1)+1)
   ien = RESHAPE([(((RESHAPE(a(2*k-1:2*k+1,2*j-1:2*j+1,2*i-1:2*i+1),     &
            [nen]),i=1,N(1)),j=1,N(2)),k=1,N(3))],[nen,nelem])
   dirich_nods = [([(ien(1:9,i),i=1,(N(2)-1)*N(1)+1,N(1))]),j=1,nsd]

   call init_elastsol(density,x,ien,dirich_nods,[nfnod,nfnod,nfnod],     &
                      nnod,nelem)

   u = x; v = 0._rp; a = 0._rp
   open(10,file='tip_disp.dat')
   do i=1,nstep
      call deform(u,v,a,load,nnod,stiffness,poisson,dt)
      write(10,*) i*dt,u(:,ntip)
      if(mod(i,iout)==0) then
         write(fname,'(a4i3.3a4)') 'beam',i/iout,'.vtk'
         call init_vtk(fname,'elastsol_beamtest')
         call vtkwrite_mesh(x,ien,'HEXA')
         call vtkwrite_vector(u-x,'displacement')
         call exit_vtk
      end if
      write(logfile,*) 'Completed time step:',i
   end do
   call exit_elastsol
   close(10)
end program beam
