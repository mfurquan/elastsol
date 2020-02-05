!=======================================================================
! This module defines data structure for
! Compressed Sparse Row, 3 array variant (CSR3) matrix structure 
! Values of entries are to be supplied separately
!-----------------------------------------------------------------------
! M. Furquan | Oct 2017
!=======================================================================
include "mkl_dss.f90"
module csr3s
   !use iso_fortran_env
   use mkl_dss
   implicit none

   integer,parameter,private :: rp = 8!REAL64

   type :: matrix_csr3
      integer :: nval, nrow
      integer,      allocatable :: col(:), row_i(:)
      type(mkl_dss_handle)      :: handle
      contains
        procedure :: create
        procedure :: add2entry, replace_entry
        procedure :: times
        procedure :: solve
        procedure :: free
   end type matrix_csr3

   contains

      subroutine free(this)
         implicit none
         class(matrix_csr3) :: this
         integer           :: ierr
         deallocate(this%col)
         deallocate(this%row_i)
         ierr = dss_delete(this%handle,MKL_DSS_DEFAULTS)
      end subroutine free

      pure subroutine add2entry(this,val,irow,icol,a)
         implicit none
         class(matrix_csr3),intent(in)    :: this
         integer,           intent(in)    :: irow, icol
         real(kind=rp),     intent(in)    :: a
         real(kind=rp),     intent(inout) :: val(this%nval)
         integer                          :: ival

         do ival = this%row_i(irow),this%row_i(irow+1)
            if(this%col(ival)==icol) exit
         end do

         val(ival) = val(ival) + a
      end subroutine add2entry

      pure subroutine replace_entry(this,val,irow,icol,a)
         implicit none
         class(matrix_csr3),intent(in)    :: this
         integer,           intent(in)    :: irow, icol
         real(kind=rp),     intent(in)    :: a
         real(kind=rp),     intent(inout) :: val(this%nval)
         integer                          :: ival

         do ival = this%row_i(irow),this%row_i(irow+1)
            if(this%col(ival)==icol) exit
         end do

         val(ival) = a
      end subroutine replace_entry

      function solve(this,val,rhs)
         implicit none
         class(matrix_csr3),intent(inout) :: this
         real(kind=rp),     intent(in)    :: val(this%nval), rhs(this%nrow)
         real(kind=rp)                    :: solve(this%nrow)
         integer                          :: ierr

         ierr = dss_factor(this%handle,MKL_DSS_INDEFINITE,val)
         ierr = dss_solve(this%handle,MKL_DSS_DEFAULTS,rhs,1,solve)
      end function solve

      function times(this,A,v)
         class(matrix_csr3),intent(in) :: this
         real(kind=rp),     intent(in) :: A(this%nval), v(this%nrow)
         real(kind=rp)                 :: times(this%nrow)
         
         call mkl_dcsrgemv('N',this%nrow,A,this%row_i,this%col,v,times)
      end function times

      subroutine create(this,struct_matT)
         implicit none
         class(matrix_csr3),intent(inout) :: this
         logical,           intent(in)    :: struct_matT(:,:)
         integer                          :: i, j, ival, rval, &
                                             ierr, perm(1)

         this%nval = COUNT(struct_matT)
         this%nrow = SIZE(struct_matT,1)
         allocate(this%col(this%nval))
         allocate(this%row_i(this%nrow+1))

         ival = 0
         do i = 1,this%nrow
            rval = 0
            do j = 1,this%nrow
               if(struct_matT(j,i)) then
                  ival = ival + 1
                  rval = rval + 1
                  this%col(ival) = j
               end if
               if(rval==1) this%row_i(i) = ival
            end do
         end do
         this%row_i(this%nrow+1) = ival + 1

         ierr = dss_create(this%handle,MKL_DSS_DEFAULTS)
         ierr = dss_define_structure(this%handle,                     &
                MKL_DSS_NON_SYMMETRIC,this%row_i,                     &
                this%nrow,this%nrow,this%col,this%nval)
         ierr = dss_reorder(this%handle,MKL_DSS_DEFAULTS,perm)
      end subroutine create
end module csr3s
