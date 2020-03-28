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
            write(funit,'(e16.9,a,e16.9,a,e16.9)')                       &
                    v(1,i),' ',v(2,i),' ',v(3,i)
         end do

         write(funit,*)
      end subroutine vtkwrite_vector

      subroutine exit_vtk
         close(funit)
      end subroutine exit_vtk

      pure function split_HEXA(con)
         integer,intent(in) :: con(:,:)
         integer :: split_HEXA(8,8*SIZE(con,2)), i, nelem

         nelem = SIZE(con,2)
         do concurrent (i = 1:nelem)
            split_HEXA(:,(i-1)*8+1) = [con( 1,i),con( 2,i),con( 5,i),    &
                    con( 4,i),con(10,i),con(11,i),con(14,i),con(13,i)]
            split_HEXA(:,(i-1)*8+2) = [con( 2,i),con( 3,i),con( 6,i),    &
                    con( 5,i), con(11,i),con(12,i),con(15,i),con(14,i)]
            split_HEXA(:,(i-1)*8+3) = [con( 4,i),con( 5,i),con( 8,i),    &
                    con( 7,i),con(13,i),con(14,i),con(17,i),con(16,i)]
            split_HEXA(:,(i-1)*8+4) = [con( 5,i),con( 6,i),con( 9,i),    &
                    con( 8,i),con(14,i),con(15,i),con(18,i),con(17,i)]
            split_HEXA(:,(i-1)*8+5) = [con(10,i),con(11,i),con(14,i),    &
                    con(13,i),con(19,i),con(20,i),con(23,i),con(22,i)]
            split_HEXA(:,(i-1)*8+6) = [con(11,i),con(12,i),con(15,i),    &
                    con(14,i),con(20,i),con(21,i),con(24,i),con(23,i)]
            split_HEXA(:,(i-1)*8+7) = [con(13,i),con(14,i),con(17,i),    &
                    con(16,i),con(22,i),con(23,i),con(26,i),con(25,i)]
            split_HEXA(:,(i-1)*8+8) = [con(14,i),con(15,i),con(18,i),    &
                    con(17,i),con(23,i),con(24,i),con(27,i),con(26,i)]
               end do
      end function split_HEXA
end module vtk_legacy
