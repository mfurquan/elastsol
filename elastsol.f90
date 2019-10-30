module elastsol
   use iso_fortran_env, only: rp => REAL64
   use hexa_lagr
   use csr3s, only: matrix_csr3 => matrix
   implicit none

   logical,      parameter,private :: debug = .TRUE.,                    &
                                      sparse_matrix = .TRUE.
   real(kind=rp),parameter,private :: tol = 1.E-7, beta = 0.25_rp,       &                                        au = 0.5_rp

   integer,parameter                      :: ned = nsd*nen, logfile=7
   integer,                  private,save :: nnod, nele, neqn
   integer,      allocatable,private,save :: id(:,:), conn(:,:)
   real(kind=rp),allocatable,private,save :: X(:,:), M(:), a_tp(:,:)
   type(matrix),             private,save :: mat

contains
   subroutine init_elastsol(rho0,coor,con,fix_nod,nfix,numnp,numel)
      real(kind=rp),intent(in) :: rho0, coord(nsd,numnp)
      integer,intent(in)  :: numnp, numel, nfix(nsd),                    &
                             fix_nod(SUM(nfix)), con(nen,numel)
      logical,allocatable :: dmat(:,:)

      nnod = numnp; nele = numel
      ! save coordinates and connectivity array
      X    = coord; conn = con

      ! generate ID array
      allocate(id(nsd,nnod))
      call gen_id(id,neqn)

      ! find sparse matrix strcuture
      if(sparse_matrix)
         allocate(dmat(neqn,neqn))
         call get_matrix_struct(dmat)
         call mat%create(dmat)
         deallocate(dmat)
      end if

      ! form mass matrix
      allocate(M(mat%nval))
      do concurrent (iele = 1:nele)
         i_ele = iloc(iele,id)
         X_ele = floc(iele, X)
         call assemble_mat(M,i_ele,m_ele(X_ele))
      end do
   contains
      pure subroutine gen_id(idn,neq)
         integer,intent(inout) :: idn, neq
         integer :: isd, nstart, nend, inod

         ! tag Dirichlet nodes
         id = 0; nend = 0
         do isd = 1,nsd
            nstart = nend + 1
            nend   = nend + nfix(isd)
            do inod = nstart,nend
               idn(isd,fix_nod(inod)) = -1
            end do
         end do

         ! count equations/unknowns
         neq = 0
         do inod = 1,nnod
            do isd = 1,nsd
               if(idn(isd,inod)==0) then
                  neq = neq + 1
                  idn(isd,inod) = neq
               end if
            end do
         end do
      end subroutine gen_id

      pure subroutine get_matrix_struct(A)
         logical,intent(inout) :: A(neqn,neqn)
         integer :: isd, jsd, ien, jen, ied, jed, iele 

         A = .FALSE.
         do concurrent (iele = 1:nele, ien = 1:nen, jen = 1:nen          &
                                       isd = 1:nsd, jsd = 1:nsd)
            ied = id(isd,conn(ien,iele))
            jed = id(jsd,conn(jen,iele))
            A(ied,jed) = .TRUE.
         end do 
      end subroutine get_matrix_struct

      pure function m_ele(Xe)
         real(kind=rp),intent(in) :: Xe(nsd,nen)
         real(kind=rp)            :: m_ele(ned,ned), intgrl
         integer                  :: isd, ien, jen

         do concurrent (ien = 1:nen, jen = 1:nen)
            intgrl = rho0*integ(phi(:,ien)*phi(:,jen),Xe)
            do concurrent (isd = 1:nsd)
               m_ele((ien-1)*nsd+isd,(jen-1)*nsd+isd) = intgrl
            end do
         end do
      end function m_ele
   end subroutine init_elastsol

   subroutine exit_elastsol
      deallocate(conn)
      deallocate(id)
      deallocate(X)
      deallocate(M)
      deallocate(a_tp)
      deallocate(du)
      call mat%free
   end subroutine exit_elastsol

   subroutine deform(u,v,a,F,Fnod,nF,E,nu,dt)
      real(kind=rp),dimension(nsd,nnod),intent(inout) :: u, v, a
      integer,      intent(in) :: nF, Fnod(nF)
      real(kind=rp),intent(in) :: F(nsd,nF), E, nu, dt
      real(kind=rp) :: resnorm, F_ele(nqd,nsd,nsd), lambda, mu
      integer :: iter, iele

      lambda = nu*E/((1._rp+nu)(1._rp-2._rp*nu))
      mu     =    E/((1._rp+nu)*2._rp)

      ! initialize Newton-Raphson iterations
      a_tp = a
      a    = -(v/(dt*beta) + (0.5_rp-beta)*a_tp/beta)
      v    = v + dt*((1._rp-tau)*a_tp + tau*a)

      newton : do iter =1,maxiter
         res = 0._rp; Kaug = 0._rp
         do concurrent (iele = 1:nele)
            i_ele = iloc(iele,id)
            X_ele = floc(iele,X)
            F_ele = grad(floc(iele,u),X_ele)
            ! assembling internal force vector
            call assemble_vec(res, i_ele,fint(F_ele,X_ele))
            ! assembling internal tangent stiffness matrix
            call assemble_mat(Kaug,i_ele,ktan(F_ele,X_ele))
         end do
         
         res = Fext(F,Fnod,nF) & ! External force vector
             - Finrt(a)        & ! Inertial force vector
             - res               ! Internal force vector

         !check convergence
         resnorm = NORM2(res)/neqn
         if(resnorm > tol) then
            write(logfile,*) 'newton iter:',iter,'||res||=',resnorm
         else
            exit newton
         end if


         ! finding displacement increment
         Kaug =  Kaug + M/(dt*dt*beta)  ! augmented stiffness matrix
         call mat%solve(Kaug,du,res)

         call update_statvec(u,du)
         call update_statvec(v,du*tau/(beta*dt))
         call update_statvec(a,du/(beta*dt*dt))
      end do newton

      if(iter >= maxiter) erro stop 'In deform: failed to converge in &
              & Newtom iterations!'
   contains
      pure function Finrt(acc)
         real(kind=rp),intent(in) :: acc(nsd,nnod)
         real(kind=rp)            :: Finrt(neqn), ac(neqn)

         call renew_solvec(ac,acc)
         call mkl_dcsrgemv('N',neqn,M,mat%row_i,mat%col,ac,Finrt)
      end function Finrt

      pure function Fext(frc,fr_nod,nfn)
         real(kind=rp),intent(in) :: frc(nsd,nfn)
         integer,      intent(in) :: fr_nod(nfn), nfn
         real(kind=rp) :: Fext(neqn)

         Fext = 0._rp
         do concurrent (isd = 1:nsd, ifn = 1:nfn, id(isd,fr_nod(ifn))>0)
            Fext(id(isd,fr_nod(ifn))) = frc(isd,ifn)
         end do
      end function Fext

      pure function fint(Fe,Xe)
         real(kind=rp),intent(in) :: Fe(nqd,nsd,nsd), Xe(nsd,nen)
         real(kind=rp) :: fint(ned), phi_X(nqd,nsd,nen), S(nqd,ntd)
         integer :: ien

         S = Stress(Fe)
         phi_X = dphi_d(Xe)
         do concurrent (ien = 1:nen)
            fint((ien-1)*nsd+1:ien*nsd)                                  &
               = integn(MATMULQ(BT(Fe,phi_X(:,:,ien)),S),Xe)
         end do
      end function fint

      pure function ktan(Fe,Xe)
         real(kind=rp),intent(in) :: Fe(nqd,nsd,nsd), Xe(nsd,nen)
         real(kind=rp) :: ktan(ned,ned), BL(nqd,nsd,ntd,nen)

         phi_X = dphi_d(Xe)
         do concurrent (ien = 1:nen)
            BL(:,:,:,ien) = BT(Fe,phi_X(:,:,ien))
         end do

         do concurrent (ien = 1:nen, jen = 1:nen)
            ktan((ien-1)*nsd+1:ien*nsd,(jen-1)*nsd+1:jen*nsd)
               = inregn(MATMULQ(BL(:,:,:,ien),DT(BL(:,:,:,jen))),Xe)
         end do
      contains
         pure function DT(T)
            real(kind=rp),intent(in) :: T(nqd,nsd,ntd)
            real(kind=rp) :: DT(nqd,ntd,nsd), trace(nqd,nsd)

            trace = T(:,:,1)+T(:,:,2)+T(:,:,3)
            DT(:,1,:) = 2._rp*mu*T(:,:,1) + lambda*trace
            DT(:,2,:) = 2._rp*mu*T(:,:,2) + lambda*trace
            DT(:,3,:) = 2._rp*mu*T(:,:,3) + lambda*trace
            DT(:,4,:) =       mu*T(:,:,4)
            DT(:,5,:) =       mu*T(:,:,5)
            DT(:,6,:) =       mu*T(:,:,6)
         end function DT
      end function ktan

      pure function Stress(F)
         real(kind=rp),intent(in) :: F(nqd,nsd,nsd)
         real(kind=rp) :: Stress(nqd,ntd), GLE(nqd,nsd,nsd), lambda, mu

         ! calc: F^T.F
         GLE = 0._rp
         do concurrent (iqd = 1:nqd,isd = 1:nsd, jsd = 1:nsd)
            GLE(iqd,isd,jsd) = GLE(iqd,isd,jsd)
                             + F(iqd,jsd,isd)*F(iqd,isd,jsd)
         end do
         GLE = MATMUL(TRANSPOSE(F),F)/2._rp
         GLE(:,1,1) = GLE(:,1,1)-0.5_rp
         GLE(:,2,2) = GLE(:,2,2)-0.5_rp
         GLE(:,3,3) = GLE(:,3,3)-0.5_rp
         Stress = D(GLE)
      contains
         pure function D(T)
            real(kind=rp),intent(in) :: T(nqd,nsd,nsd)
            real(kind=rp) :: D(nqd,ntd), trace(nqd,nsd,nsd)

            trace  = T(:,1,1) + T(:,2,2) + T(:,3,3)
            D(:,1) = 2._rp*mu*T(:,1,1) + lambda*trace
            D(:,2) = 2._rp*mu*T(:,2,2) + lambda*trace
            D(:,3) = 2._rp*mu*T(:,3,3) + lambda*trace
            D(:,4) =       mu*(T(:,1,2)+T(:,2,1))
            D(:,5) =       mu*(T(:,2,3)+T(:,3,2))
            D(:,6) =       mu*(T(:,1,3)+T(:,3,1))
         end function D
      end function Stress

      pure function BT(F,p_X)
         real(kind=rp),intent(in) :: F(nqd,nsd,nsd), p_X(nqd,nsd)
         integer,intent(in) :: ien
         real(kind=rp) :: BT(nqd,nsd,ntd)

         do concurrent (iqd = 1:nqd)
            BT(iqd,:,1) = [F(iqd,1,1)*p_X(iqd,1),                        &
                           F(iqd,2,1)*p_X(iqd,1),                        &
                           F(iqd,3,1)*p_X(iqd,1)]
            BT(iqd,:,2) = [F(iqd,1,2)*p_X(iqd,2),                        &
                           F(iqd,2,2)*p_X(iqd,2),                        &
                           F(iqd,3,2)*p_X(iqd,2)]
            BT(iqd,:,3) = [F(iqd,1,3)*p_X(iqd,3),                        &
                           F(iqd,2,3)*p_X(iqd,3),                        &
                           F(iqd,3,3)*p_X(iqd,3)]
            BT(iqd,:,4) = [F(iqd,1,1)*p_X(iqd,2)+F(iqd,1,2)*p_X(iqd,1),  &
                           F(iqd,2,1)*p_X(iqd,2)+F(iqd,2,2)*p_X(iqd,1),  &
                           F(iqd,3,1)*p_X(iqd,2)+F(iqd,3,2)*p_X(iqd,1)]
            BT(iqd,:,5) = [F(iqd,1,2)*p_X(iqd,3)+F(iqd,1,3)*p_X(iqd,2),  &
                           F(iqd,2,2)*p_X(iqd,3)+F(iqd,2,3)*p_X(iqd,2),  &
                           F(iqd,3,2)*p_X(iqd,3)+F(iqd,3,3)*p_X(iqd,2)]
            BT(iqd,:,6) = [F(iqd,1,1)*p_X(iqd,3)+F(iqd,1,3)*p_X(iqd,1),  &
                           F(iqd,2,1)*p_X(iqd,3)+F(iqd,2,3)*p_X(iqd,1),  &
                           F(iqd,3,1)*p_X(iqd,3)+F(iqd,3,3)*p_X(iqd,1)]
         end do
      end function BT
   end subroutine deform

!================ M A T H E M A T I C A L   T O O L S ====================
   pure function integ(intgrnd,xe)
      real(kind=rp),intent(in) :: intgrnd(nqd), xe(nsd,nen)
      real(kind=rp)            :: integ, jac(nqd)
      integer                  :: iqd

      do concurrent (iqd = 1:nqd)
         jac(iqd)   = det(d_dxi(xe,iqd))
      end do
      integ = SUM(intgrnd*wq*jac)
   end function integ

   pure function integn(intgrnd,xe)
      real(kind=rp),intent(in) :: intgrnd(nqd,nsd), xe(nsd,nen)
      real(kind=rp) :: integn(nsd)
      integer :: iqd, isd

      integn = 0._rp
      do concurrent (iqd = 1:nqd, isd = 1:nsd)
         integn(isd) = integn(isd)                                       &
                     + intgrnd(iqd,isd)*wq(iqd)*det(d_dxi(xe,iqd))
      end do
   end function integn

   pure function dphi_d(Xe)
      real(kind=rp),intent(in) :: Xe(nsd,nen)
      real(kind=rp)            :: dphi_d(nqd,nsd,nen)

      do concurrent (iqd = 1:nqd, ien = 1:nen)
         dphi_d(iqd,:,ien) = MATMUL(inv(d_dxi(Xe,iqd)),phi_xi(:,ien,iqd))
      end do
   contains
      pure function inv(A)
         real(kind=rp),intent(in) :: A(nsd,nsd)
         real(kind=rp)            :: inv(nsd,nsd)

         inv(:,1) = [A(2,2)*A(3,3)-A(3,2)*A(2,3),                        &
                     A(2,3)*A(3,1)-A(3,3)*A(2,1),                        &
                     A(2,1)*A(3,2)-A(3,1)*A(2,2)]
         inv(:,2) = [A(1,3)*A(3,2)-A(3,3)*A(1,2),                        &
                     A(1,1)*A(3,3)-A(3,1)*A(1,3),                        &
                     A(1,2)*A(3,1)-A(3,2)*A(1,1)]
         inv(:,3) = [A(1,2)*A(2,3)-A(2,2)*A(1,3),                        &
                     A(1,3)*A(2,1)-A(2,3)*A(1,1),                        &
                     A(1,1)*A(2,2)-A(2,1)*A(1,2)]
         inv = inv/det(A)
      end function inv
   end function dphi_d

   pure function grad(ue,Xe)
      real(kind=rp),intent(in) :: ue(nsd,nen), Xe(nsd,nen)
      real(kind=rp) :: grad(nqd,nsd,nsd)

      do concurrent (iqd = 1:nqd,isd = 1:nen, jsd = 1:nsd, ken = 1:nen)
         grad(iqd,isd,jsd) = grad(iqd,isd,jsd)                           &
                           + ue(isd,ken)*dphi_dX(iqd,jsd,ken)
      end do
   end function grad

   pure function d_dxi(ue,iqd)
      real(kind=rp),intent(in) :: ue(nsd,nen)
      integer,      intent(in) :: iqd
      real(kind=rp)            :: d_dxi(nsd,nsd)
      integer                  :: isd, jsd, ken

      d_dxi = 0._rp
      do concurrent (isd = 1:nsd, jsd = 1:nsd, ken = 1:nen)
         d_dxi(isd,jsd) = d_dxi(isd,jsd) + ue(isd,ken)*phi_xi(jsd,ken,iqd)
      end do
   end function d_dxi

   pure function det(A)
      real(kind=rp),intent(in) :: A(nsd,nsd)
      real(kind=rp)            :: det

      det = A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))                      &
          + A(1,2)*(A(3,1)*A(2,3)-A(2,1)*A(3,3))                      &
          + A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
   end function det

   pure function MATMULQ(A,B)
      real(kind=rp),intent(in) :: A(nqd,nsd,ntd), B(nqd,ntd)
      real(kind=rp) :: MATMULQ(nqd,nsd)
      integer :: iqd, isd, itd

      MATMULQ = 0._rp
      do concurrent (iqd = 1:nqd, isd = 1:nsd, itd = 1:ntd)
         MATMULQ(iqd,isd) = MATMULQ(iqd,isd) + A(iqd,isd,itd)*B(iqd,itd)
      end do
   end function MATMULQ

!================ A S S E M B L Y    P R O C E D U R E S =================
   pure subroutine assemble_mat(A,ide,ae)
      real(kind=rp),intent(inout) :: A(neqn,neqn)
      real(kind=rp),intent(in)    :: ae(ned,ned)
      integer,      intent(in)    :: ide(ned)
      integer                     :: ied, jed

      do concurrent (ied = 1:ned, jed = 1:ned, ide(ied)>.AND.ide(jed)>0)
         call mat%add2entry(A,ide(ied),ide(jed),ae(ied,jed))
      end do
   end subroutine assemble_mat

   pure subroutine assemble_vec(V,ide,ve)
      real(kind=rp),intent(inout) :: V(neqn)
      real(kind=rp),intent(in)    :: ve(ned) 
      integer,      intent(in)    :: ide(ned)
      integer                     :: ied

      do concurrent (ied = 1:ned)
         V(ide(ied)) = V(ide(ied)) + ve(ied)
      end do
   end subroutine assemble_vec

   pure function iloc(iele,a)
      implicit none
      integer,intent(in) :: iele, a(nsd,nnod)
      integer            :: iloc(nsd,nen), conn_ele(nen), isd, ien

      conn_ele = conn(:,iele)
      do concurrent (isd = 1:nsd, ien = 1:nen)
         iloc(isd,ien) = a(isd,conn_ele(ien))
      end do
   end function iloc

   pure function floc(iele,a)
      implicit none
      integer,      intent(in) :: iele
      real(kind=rp),intent(in) :: a(nsd,nnod)
      integer                  :: iloc(nsd,nen), conn_ele(nen), isd, ien

      conn_ele = conn(:,iele)
      do concurrent (isd = 1:nsd, ien = 1:nen)
         iloc(isd,ien) = a(isd,conn_ele(ien))
      end do
   end function floc

   pure function update_statvec(q,dq)
      real(kind=rp),intent(in) :: q(nsd,nnod), dq(neqn)
      integer                  :: isd, inod

      do concurrent (isd = 1:nsd, inod = 1:nnod, id(isd,inod)>0)
         q(isd,inod) = q(isd,inod) + dq(id(isd,inod))
      end do
   end function update_statvec

   pure subroutine renew_solvec(x,y)
      real(kind=rp),intent(inout) :: x(neqn)
      real(kind=rp),intent(in)    :: y(nsd,nnod)
      integer                     :: isd, inod

      do concurrent (isd = 1:nsd, inod = 1:nnod, id(isd,inod)>0)
         x(id(isd,inod)) = y(isd,inod)
      end do
   end subroutine renew_solvec
end module elastsol
