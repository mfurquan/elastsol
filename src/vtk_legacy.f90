module vtk_legacy
   implicit none
   integer,parameter,private :: rp = 8, nsd = 3
   integer,private,save :: nn, ne, nen, funit

   contains
       subroutine init_vtk(filename,title)
         implicit none
         character(len=*),intent(in) :: filename, title

         open(newunit=funit,file=filename)
         write(funit,'(a26)') "# vtk DataFile Version 2.0"
         write(funit,'(a15)') title
         write(funit,'(a5)') "ASCII"
         write(funit,*)
         write(funit,'(a25)') "DATASET UNSTRUCTURED_GRID"
      end subroutine init_vtk

      subroutine vtkwrite_mesh(x,conn,elemtype)
         implicit none
         character(len=*),intent(in) :: elemtype
         integer,intent(in) :: conn(:,:)
         real(kind=rp),intent(in) :: x(:,:)
         integer :: i, j
         character(len=13) :: fstr

         nn  = size(x,2)
         nen = size(conn,1)
         ne  = size(conn,2)
         write(funit,'(a7,i0,a6)') 'POINTS ',nn,' float'
         do i = 1,nn
            !write(funit,*) x(1,i),x(2,i),x(3,i)
            write(funit,'(e16.9,a,e16.9,a,e16.9)') x(1,i),' ',x(2,i),' ',x(3,i)
         end do

         write(funit,*)
         write(fstr,'(a4,i0,a7)') '(i0,',nen,'(a,i0))'
         write(funit,'(a6,i0,a,i0)') "CELLS ",ne,' ',(nen+1)*ne
         do i = 1,ne
            write(funit,trim(fstr)) nen,(' ',conn(j,i)-1,j=1,nen)
         end do

         write(funit,*)
         write(funit,'(a11,i0)') "CELL_TYPES ",ne
         select case(elemtype)
         case('TRI')
            nen = 5
         case('QUAD')
            nen = 9
         case('TETRA')
            nen = 10
         case('HEXA')
            nen = 12
         case('WEDGE')
            nen = 13
         end select
         do i = 1,ne
            write(funit,'(i0)') nen
         end do
         write(funit,*)

         write(funit,'(a11i0)') "POINT_DATA ",nn
      end subroutine vtkwrite_mesh

      subroutine vtkwrite_scalar(q,variable)
         implicit none
         real(kind=rp),intent(in) :: q(nn)
         character(len=*),intent(in) :: variable
         character(len=5) :: fstr
         integer :: i

         write(fstr,'(a2,i2,a)') '(a',14+len(variable),')'
         write(funit,fstr) 'SCALARS '//variable//' float'
         write(funit,'(a20)') 'LOOKUP_TABLE default'
         do i = 1,nn
            write(funit,'(f0.10)') q(i)
         end do

         write(funit,*)
      end subroutine vtkwrite_scalar

      subroutine vtkwrite_vector(v,variable)
         implicit none
         real(kind=rp),intent(in) :: v(nsd,nn)
         character(len=*),intent(in) :: variable
         character(len=5) :: fstr
         integer :: i

         write(fstr,'(a2,i2,a)') '(a',14+len(variable),')'
         write(funit,fstr) 'VECTORS '//variable//' float'
         do i = 1,nn
            !write(funit,*) v(1,i),v(2,i),v(3,i)
            write(funit,'(e16.9,a,e16.9,a,e16.9)') v(1,i),' ',v(2,i),' ',v(3,i)
         end do

         write(funit,*)
      end subroutine vtkwrite_vector

      subroutine exit_vtk
         close(funit)
      end subroutine exit_vtk
end module vtk_legacy
