module elastsol
   !use iso_fortran_env, only: rp => REAL64
   use hexa_lagr
   use csr3s, only: matrix => matrix_csr3
   implicit none

   logical,      parameter,private :: debug = .TRUE.,                    &
                                      sparse_matrix = .TRUE.
   real(kind=rp),parameter,private :: tol = 1.E-7, beta = 0.25_rp,       &
                                      tau = 0.5_rp
   integer,      parameter,private :: maxiter = 50

   integer,parameter                      :: ned = nsd*nen, logfile=7
   integer,                  private,save :: nnod, nele, neqn
   integer,      allocatable,private,save :: id(:,:), conn(:,:)
   real(kind=rp),allocatable,private,save :: X(:,:), M(:), Kaug(:),      &
                                             a_tp(:,:), du(:)
   type(matrix),             private,save :: mat

!======================= I N T E R F A C E S =============================
   interface operator (.times.)    ! UV
      module procedure A_times_v
   end interface operator (.times.)

   interface operator (.dot.)      ! U^T V
      module procedure A_dot_v, u_dot_v, A_dot_B
   end interface operator (.dot.)

   interface integ      ! integrate
      module procedure integ_scal, integ_vec, integ_mat
   end interface integ
!-------------------------------------------------------------------------

contains
   subroutine init_elastsol(rho0,coord,con,fix_nod,nfix,numnp,numel)
      integer,      intent(in) :: numnp, numel, nfix(nsd),               &
                                  fix_nod(SUM(nfix)), con(nen,numel)
      real(kind=rp),intent(in) :: rho0, coord(nsd,numnp)
      logical,allocatable      :: dmat(:,:)
      integer       :: i_ele(nsd,nen), iele
      real(kind=rp) :: X_ele(nsd,nen)

      nnod = numnp; nele = numel
      ! save coordinates and connectivity array
      X    = coord; conn = con

      ! generate ID array
      allocate(id(nsd,nnod))
      call gen_id(id,neqn)

      ! find sparse matrix structure
      if(sparse_matrix) then
         allocate(dmat(neqn,neqn))
         call get_matrix_struct(dmat)
         call mat%create(dmat)
         deallocate(dmat)
      end if

      ! matrices and solution vector
      allocate(M(mat%nval))
      allocate(Kaug(mat%nval))
      allocate(du(neqn))

      ! form mass matrix
      do concurrent (iele = 1:nele)
         i_ele = iloc(iele,id)
         X_ele = floc(iele, X)
         call assemble_mat(M,i_ele,m_ele(X_ele))
      end do
   contains
      pure subroutine gen_id(idn,neq)
         integer,intent(inout) :: idn(nsd,nnod), neq
         integer :: isd, nstart, nend, inod

         ! tag Dirichlet nodes
         idn = 0; nend = 0
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
         do concurrent (iele = 1:nele, ien = 1:nen, jen = 1:nen,         &
                                       isd = 1:nsd, jsd = 1:nsd)
            ied = id(isd,conn(ien,iele))
            jed = id(jsd,conn(jen,iele))
            if(ied>0 .AND. jed>0) A(ied,jed) = .TRUE.
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
      deallocate(Kaug)
      deallocate(a_tp)
      deallocate(du)
      call mat%free
   end subroutine exit_elastsol

   subroutine deform(u,v,a,F,Fnod,nF,E,nu,dt,steady)
      real(kind=rp),dimension(nsd,nnod),intent(inout) :: u, v, a
      integer,      intent(in) :: nF, Fnod(nF)
      real(kind=rp),intent(in) :: F(nsd,nF), E, nu, dt
      logical,      intent(in),optional :: steady
      real(kind=rp) :: resnorm, F_ele(nqd,nsd,nsd), X_ele(nsd,nen),      &
                       res(neqn), stdy_switch, lambda, mu
      integer :: iter, iele, i_ele(nsd,nen)

      stdy_switch = 1._rp
      if(present(steady)) then
         if(steady) stdy_switch = 0._rp
      end if

      lambda = nu*E/((1._rp+nu)*(1._rp-2._rp*nu))
      mu     =    E/((1._rp+nu)*2._rp)

      ! initialize Newton-Raphson iterations
      a_tp = a
      a    = -(v/(dt*beta) + (0.5_rp-beta)*a_tp/beta)
      v    = v + dt*((1._rp-tau)*a_tp + tau*a)

      newton : do iter = 1,maxiter
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
         
         res = Fext(F,Fnod,nF)       & ! External force vector
             - Finrt(a)*stdy_switch  & ! Inertial force vector
             - res                     ! Internal force vector

         !check convergence
         resnorm = NORM2(res)/neqn
         if(resnorm > tol) then
            write(logfile,*) 'newton iter:',iter,'||res||=',resnorm
         else
            exit newton
         end if


         ! finding displacement increment
         Kaug = Kaug + M/(dt*dt*beta)*stdy_switch !augment stiffness matrix
         du   = mat%solve(Kaug,res)

         call update_statvec(u,du)
         call update_statvec(v,du*tau/(beta*dt))
         call update_statvec(a,du/(beta*dt*dt))
      end do newton

      if(iter >= maxiter) error stop 'In deform: failed to converge in &
              & Newton iterations!'
   contains
      function Finrt(acc)
         real(kind=rp),intent(in) :: acc(nsd,nnod)
         real(kind=rp)            :: Finrt(neqn), ac(neqn)

         call renew_solvec(ac,acc)
         Finrt = mat%times(M,ac)
      end function Finrt

      pure function Fext(frc,fr_nod,nfn)
         real(kind=rp),intent(in) :: frc(nsd,nfn)
         integer,      intent(in) :: fr_nod(nfn), nfn
         real(kind=rp) :: Fext(neqn)
         integer       :: isd, ifn

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
               = integ(B(Fe,phi_X(:,:,ien)).dot.S,Xe)
         end do
      end function fint

      pure function ktan(Fe,Xe)
         real(kind=rp),intent(in) :: Fe(nqd,nsd,nsd), Xe(nsd,nen)
         real(kind=rp) :: ktan(ned,ned), BL(nqd,ntd,nsd,nen),            &
                          S(nqd,nsd,nsd), phi_X(nqd,nsd,nen)
         integer :: ien, jen
         real(kind=rp),parameter :: I(nsd,nsd) = RESHAPE(                &
            [1._rp,0._rp,0._rp,                                          &
             0._rp,1._rp,0._rp,                                          &
             0._rp,0._rp,1._rp],[nsd,nsd])

         phi_X = dphi_d(Xe)
         do concurrent (ien = 1:nen)
            BL(:,:,:,ien) = B(Fe,phi_X(:,:,ien))
         end do

         S = Voigt2tensor(Stress(Fe))
         
         do concurrent (ien = 1:nen, jen = 1:nen)
            ktan((ien-1)*nsd+1:ien*nsd,(jen-1)*nsd+1:jen*nsd)            &
               = integ(phi_X(:,:,ien).dot.(S.times.phi_X(:,:,jen)),Xe)*I &
               + integ(BL(:,:,:,ien).dot.C(BL(:,:,:,jen)),Xe)
         end do
      end function ktan

      pure function Stress(F)
         real(kind=rp),intent(in) :: F(nqd,nsd,nsd)
         real(kind=rp) :: Stress(nqd,ntd), L(nqd,nsd,nsd), E(nqd,ntd)

         L      = (F.dot.F)/2._rp
         E(:,1) =     L(:,1,1) - 0.5_rp
         E(:,2) =     L(:,2,2) - 0.5_rp
         E(:,3) =     L(:,3,3) - 0.5_rp
         E(:,4) = 2.0*L(:,1,2)
         E(:,5) = 2.0*L(:,2,3)
         E(:,6) = 2.0*L(:,1,3)
         Stress = ConstitutiveLaw(E)
      end function Stress

      pure function Voigt2tensor(S) result (T)
         real(kind=rp),intent(in) :: S(nqd,ntd)
         real(kind=rp) :: T(nqd,nsd,nsd)

         T(:,1,1) = S(:,1); T(:,1,2) = S(:,4); T(:,1,3) = S(:,6)
         T(:,2,1) = S(:,4); T(:,2,2) = S(:,3); T(:,2,3) = S(:,5)
         T(:,3,1) = S(:,6); T(:,3,2) = S(:,5); T(:,3,3) = S(:,3)
      end function Voigt2tensor

      pure function ConstitutiveLaw(T) Result (D)
        real(kind=rp),intent(in) :: T(nqd,ntd)
        real(kind=rp) :: D(nqd,ntd), trace(nqd)

        ! St. Venant-Kirchhoff's material
        trace = T(:,1)+T(:,2)+T(:,3)
        D(:,1) = 2._rp*mu*T(:,1) + lambda*trace
        D(:,2) = 2._rp*mu*T(:,2) + lambda*trace
        D(:,3) = 2._rp*mu*T(:,3) + lambda*trace
        D(:,4) =       mu*T(:,4)
        D(:,5) =       mu*T(:,5)
        D(:,6) =       mu*T(:,6)
     end function ConstitutiveLaw

      pure function C(T) ! Linearized elasticity law
        real(kind=rp),intent(in) :: T(nqd,ntd,nsd)
        real(kind=rp) :: C(nqd,ntd,nsd), trace(nqd,nsd)

        trace = T(:,1,:)+T(:,2,:)+T(:,3,:)
        C(:,1,:) = 2._rp*mu*T(:,1,:) + lambda*trace
        C(:,2,:) = 2._rp*mu*T(:,2,:) + lambda*trace
        C(:,3,:) = 2._rp*mu*T(:,3,:) + lambda*trace
        C(:,4,:) =       mu*T(:,4,:)
        C(:,5,:) =       mu*T(:,5,:)
        C(:,6,:) =       mu*T(:,6,:)
     end function C

      pure function B(F,p_X)
         real(kind=rp),intent(in) :: F(nqd,nsd,nsd), p_X(nqd,nsd)
         real(kind=rp) :: B(nqd,ntd,nsd)
         integer       :: iqd

         do concurrent (iqd = 1:nqd)
            B(iqd,:,1) = [F(iqd,1,1)*p_X(iqd,1),                         &
                          F(iqd,1,2)*p_X(iqd,2),                         &
                          F(iqd,1,3)*p_X(iqd,3),                         &
                          F(iqd,1,1)*p_X(iqd,2)+F(iqd,1,2)*p_X(iqd,1),   &
                          F(iqd,1,2)*p_X(iqd,3)+F(iqd,1,3)*p_X(iqd,2),   &
                          F(iqd,1,1)*p_X(iqd,3)+F(iqd,1,3)*p_X(iqd,1)]
            B(iqd,:,2) = [F(iqd,2,1)*p_X(iqd,1),                         &
                          F(iqd,2,2)*p_X(iqd,2),                         &
                          F(iqd,2,3)*p_X(iqd,3),                         &
                          F(iqd,2,1)*p_X(iqd,2)+F(iqd,2,2)*p_X(iqd,1),   &
                          F(iqd,2,2)*p_X(iqd,3)+F(iqd,2,3)*p_X(iqd,2),   &
                          F(iqd,2,1)*p_X(iqd,3)+F(iqd,2,3)*p_X(iqd,1)]
            B(iqd,:,3) = [F(iqd,3,1)*p_X(iqd,1),                         &
                          F(iqd,3,2)*p_X(iqd,2),                         &
                          F(iqd,3,3)*p_X(iqd,3),                         &
                          F(iqd,3,1)*p_X(iqd,2)+F(iqd,3,2)*p_X(iqd,1),   &
                          F(iqd,3,2)*p_X(iqd,3)+F(iqd,3,3)*p_X(iqd,2),   &
                          F(iqd,3,1)*p_X(iqd,3)+F(iqd,3,3)*p_X(iqd,1)]
         end do
      end function B
   end subroutine deform

!================ M A T H E M A T I C A L   T O O L S ====================
   pure function integ_scal(intgrnd,xe)
      real(kind=rp),intent(in) :: intgrnd(nqd), xe(nsd,nen)
      real(kind=rp)            :: integ_scal, jac(nqd)
      integer                  :: iqd

      do concurrent (iqd = 1:nqd)
         jac(iqd)   = det(d_dxi(xe,iqd))
      end do
      integ_scal = SUM(intgrnd*wq*jac)
   end function integ_scal

   pure function integ_vec(intgrnd,xe)
      real(kind=rp),intent(in) :: intgrnd(nqd,nsd), xe(nsd,nen)
      real(kind=rp) :: integ_vec(nsd)
      integer :: iqd, isd

      integ_vec = 0._rp
      do concurrent (iqd = 1:nqd, isd = 1:nsd)
         integ_vec(isd) = integ_vec(isd)                                 &
                        + intgrnd(iqd,isd)*wq(iqd)*det(d_dxi(xe,iqd))
      end do
   end function integ_vec

   pure function integ_mat(intgrnd,xe)
      real(kind=rp),intent(in) :: intgrnd(nqd,nsd,nsd), xe(nsd,nen)
      real(kind=rp) :: integ_mat(nsd,nsd)
      integer :: iqd, isd, jsd

      integ_mat = 0._rp
      do concurrent (iqd = 1:nqd, isd = 1:nsd, jsd = 1:nsd)
         integ_mat(isd,jsd) = integ_mat(isd,jsd)                         &
                     + intgrnd(iqd,isd,jsd)*wq(iqd)*det(d_dxi(xe,iqd))
      end do
   end function integ_mat

   pure function dphi_d(Xe)
      real(kind=rp),intent(in) :: Xe(nsd,nen)
      real(kind=rp)            :: dphi_d(nqd,nsd,nen)
      integer                  :: iqd, ien

      do concurrent (iqd = 1:nqd, ien = 1:nen)
         dphi_d(iqd,:,ien) = MATMUL(phi_xi(:,ien,iqd),inv(d_dxi(Xe,iqd)))
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
      real(kind=rp) :: grad(nqd,nsd,nsd), dphi_dX(nqd,nsd,nen)
      integer       :: iqd, isd, jsd, ken

      dphi_dX = dphi_d(Xe)
      grad = 0._rp
      do concurrent (iqd = 1:nqd, isd = 1:nsd, jsd = 1:nsd, ken = 1:nen)
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

      det = A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))                         &
          + A(1,2)*(A(3,1)*A(2,3)-A(2,1)*A(3,3))                         &
          + A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
   end function det

   pure function A_times_v(A,v) ! calc matrix-vector product: A v
      real(kind=rp), intent(in) :: A(:,:,:), v(:,:)
      real(kind=rp) :: A_times_v(nqd,SIZE(v,2))
      integer       :: m, n, iqd, i, j

      m = SIZE(A,2); n = SIZE(A,3)
      A_times_v = 0._rp
      do concurrent (iqd = 1:nqd, i = 1:m, j = 1:n)
         A_times_v(iqd,i) = A_times_v(iqd,i) + A(iqd,i,j)*v(iqd,j)
      end do
   end function A_times_v

   pure function A_dot_v(A,v) ! calc matrix-vector product: A.v = A^T v
      real(kind=rp), intent(in) :: A(:,:,:), v(:,:)
      real(kind=rp) :: A_dot_v(nqd,SIZE(v,2))
      integer       :: m, n, iqd, i, j

      m = SIZE(A,2); n = SIZE(A,3)
      A_dot_v = 0._rp
      do concurrent (iqd = 1:nqd, i = 1:m, j = 1:n)
         A_dot_v(iqd,j) = A_dot_v(iqd,j) + A(iqd,i,j)*v(iqd,i)
      end do
   end function A_dot_v

   pure function u_dot_v(u,v) ! calc matrix-vector product: u.v = u^T v
      real(kind=rp), intent(in) :: u(:,:), v(:,:)
      real(kind=rp) :: u_dot_v(nqd)
      integer       :: n, iqd, i

      n = SIZE(u,2)
      u_dot_v = 0._rp
      do concurrent (iqd = 1:nqd, i = 1:n)
         u_dot_v(iqd) = u_dot_v(iqd) + u(iqd,i)*v(iqd,i)
      end do
   end function u_dot_v

   pure function A_dot_B(A,B) ! calc matrix-matrix product: A.B = A^T B
      real(kind=rp), intent(in) :: A(:,:,:), B(:,:,:)
      real(kind=rp) :: A_dot_B(nqd,SIZE(A,3),SIZE(B,3))
      integer       :: l, m, n, iqd, i, j, k

      l = SIZE(A,2); m = SIZE(A,3); n = SIZE(B,3)
      A_dot_B = 0._rp
      do concurrent (iqd = 1:nqd, i = 1:l, j = 1:m, k = 1:n)
         A_dot_B(iqd,j,k) = A_dot_B(iqd,j,k) + A(iqd,i,j)*B(iqd,i,k)
      end do
   end function A_dot_B

!================ A S S E M B L Y    P R O C E D U R E S =================
   pure subroutine assemble_mat(A,ide,ae)
      real(kind=rp),intent(inout) :: A(neqn,neqn)
      real(kind=rp),intent(in)    :: ae(ned,ned)
      integer,      intent(in)    :: ide(ned)
      integer                     :: ied, jed

      do concurrent (ied = 1:ned, jed = 1:ned, ide(ied)>0.AND.ide(jed)>0)
         call mat%add2entry(A,ide(ied),ide(jed),ae(ied,jed))
      end do
   end subroutine assemble_mat

   pure subroutine assemble_vec(V,ide,ve)
      real(kind=rp),intent(inout) :: V(neqn)
      real(kind=rp),intent(in)    :: ve(ned) 
      integer,      intent(in)    :: ide(ned)
      integer                     :: ied

      do concurrent (ied = 1:ned, ide(ied)>0)
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
      real(kind=rp)            :: floc(nsd,nen)
      integer                  :: conn_ele(nen), isd, ien

      conn_ele = conn(:,iele)
      do concurrent (isd = 1:nsd, ien = 1:nen)
         floc(isd,ien) = a(isd,conn_ele(ien))
      end do
   end function floc

   pure subroutine update_statvec(q,dq)
      real(kind=rp),intent(inout) :: q(nsd,nnod)
      real(kind=rp),intent(in)    :: dq(neqn)
      integer                  :: isd, inod

      do concurrent (isd = 1:nsd, inod = 1:nnod, id(isd,inod)>0)
         q(isd,inod) = q(isd,inod) + dq(id(isd,inod))
      end do
   end subroutine update_statvec

   pure subroutine renew_solvec(x,y)
      real(kind=rp),intent(inout) :: x(neqn)
      real(kind=rp),intent(in)    :: y(nsd,nnod)
      integer                     :: isd, inod

      do concurrent (isd = 1:nsd, inod = 1:nnod, id(isd,inod)>0)
         x(id(isd,inod)) = y(isd,inod)
      end do
   end subroutine renew_solvec
end module elastsol
