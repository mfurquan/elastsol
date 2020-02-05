module hexa20
   use iso_fortran_env, only : rp => REAL64
   use tensor
   implicit none

   integer,parameter :: nsd = 3, nen = 20, neq = 14, ner = 6
   real(kind=rp),parameter,private :: a = sqrt(19._rp/30._rp),          &
                                      b = sqrt(19._rp/33._rp),          &
                                      c = 1._rp/sqrt(3._rp)

   real(kind=rp),parameter ::                                           &
      wr(ner) = [SPREAD(1._rp,1,8)],                                    &
      wq(neq) = [SPREAD(320._rp/361._rp,1,6),                           &
                 SPREAD(121._rp/361._rp,1,8)],                          &
!
      xr(nsd,ner) = RESHAPE([-c, -c, -c,  c, -c, -c, -c,  c, -c,        &
                              c,  c, -c, -c, -c,  c,  c, -c,  c,        &
                             -c,  c,  c,  c,  c, c],                    &
                            [nsd,ner]),                                 &
      xq(nsd,neq) = RESHAPE([a, 0._rp, 0._rp, -a, 0._rp, 0._rp,         &
                             0._rp, a, 0._rp, 0._rp, -a, 0._rp,         &
                             0._rp, 0._rp, a, 0._rp, 0._rp, -a,         &
                             b, b, b, -b, b, b, b, -b, b, -b, -b, b, b  &
                             b, -b, b, -b, -b, -b, -b, -b, -b, b, -b],  &
                            [nsd,neq])

   contains

      pure function shapefn(xi)
         implicit none
         real(kind=rp),intent(in) :: xi(2)
         real(kind=rp)            :: shapefn(8)

         m1 = 1.0_rp - xi(1)
         m2 = 1.0_rp - xi(2)

         p1 = 1.0_rp + xi(1)
         p2 = 1.0_rp + xi(2)
         shapefn = [                                                    &
            m1*m2*(m1+m2+1.0_rp)/4.0_rp,                         &
            p1*m2*(p1+m2+1.0_rp)/4.0_rp,                         &
            p1*p2*(p1+p2+1.0_rp)/4.0_rp,                         &
            m1*p2*(m1+p2+1.0_rp)/4.0_rp,                         &
            m1*m2*p2/2.0_rp,                                        &
            p1*m2*p2/2.0_rp,                                        &
            m1*p1*m2/2.0_rp,                                        &
            m1*p1*p2/2.0_rp,                                        &
      end function shapefn

      pure function deriv_shapefn(xi)
         implicit none
         real(kind=rp),intent(in) :: xi(2)
         type(vector)             :: deriv_shapefn(8)

         m1 = 1.0 - xi(1)
         m2 = 1.0 - xi(2)

         p1 = 1.0 + xi(1)
         p2 = 1.0 + xi(2)
         deriv_shapefn = [                                              &
            vector((-m2*(m1+m2+1.0_rp)-m1*m2)/4.0_rp,          &
                   (-m1*(m1+m2+1.0_rp)-m1*m2)/4.0_rp, 0._rp),          &
            vector(( m2*(p1+m2+1.0_rp)+p1*m2)/4.0_rp,          &
                   (-p1*(p1+m2+1.0_rp)-p1*m2)/4.0_rp, 0._rp),         &
            vector((-p2*m3*(m1+m2+p3-1.0_rp)+p1*p2*m3)/8.0_rp,          &
                   (-p1*m3*(m1+m2+p3-1.0_rp)+p1*p2*m3)/8.0_rp,          &
                   ( p1*p2*(m1+m2+p3-1.0_rp)-p1*p2*m3)/8.0_rp),         &
            vector(( p2*m3*(p1+m2+p3-1.0_rp)-m1*p2*m3)/8.0_rp,          &
                   (-m1*m3*(p1+m2+p3-1.0_rp)+m1*p2*m3)/8.0_rp,          &
                   ( m1*p2*(p1+m2+p3-1.0_rp)-m1*p2*m3)/8.0_rp),         &
            vector(( m2*p3*(p1+p2+m3-1.0_rp)-m1*m2*p3)/8.0_rp,          &
                   ( m1*p3*(p1+p2+m3-1.0_rp)-m1*m2*p3)/8.0_rp,          &
                   (-m1*m2*(p1+p2+m3-1.0_rp)+m1*m2*p3)/8.0_rp),         &
            vector((-m2*p3*(m1+p2+m3-1.0_rp)+p1*m2*p3)/8.0_rp,          &
                   ( p1*p3*(m1+p2+m3-1.0_rp)-p1*m2*p3)/8.0_rp,          &
                   (-p1*m2*(m1+p2+m3-1.0_rp)+p1*m2*p3)/8.0_rp),         &
            vector((-p2*p3*(m1+m2+m3-1.0_rp)+p1*p2*p3)/8.0_rp,          &
                   (-p1*p3*(m1+m2+m3-1.0_rp)+p1*p2*p3)/8.0_rp,          &
                   (-p1*p2*(m1+m2+m3-1.0_rp)+p1*p2*p3)/8.0_rp),         &
            vector(( p2*p3*(p1+m2+m3-1.0_rp)-m1*p2*p3)/8.0_rp,          &
                   (-m1*p3*(p1+m2+m3-1.0_rp)+m1*p2*p3)/8.0_rp,          &
                   (-m1*p2*(p1+m2+m3-1.0_rp)+m1*p2*p3)/8.0_rp),         &
            vector((-p1*m2*m3+m1*m2*m3)/4.0_rp,                         &
                   (-m1*p1*m3         )/4.0_rp,                         &
                   (-m1*p1*m2         )/4.0_rp),                        &
            vector(( m2*p2*m3         )/4.0_rp,                         &
                   (-p1*p2*m3+p1*m2*m3)/4.0_rp,                         &
                   (-p1*m2*p2         )/4.0_rp),                        &
            vector((-p1*p2*m3+m1*p2*m3)/4.0_rp,                         &
                   ( m1*p1*m3         )/4.0_rp,                         &
                   (-m1*p1*p2         )/4.0_rp),                        &
            vector((-m2*p2*m3         )/4.0_rp,                         &
                   (-m1*p2*m3+m1*m2*m3)/4.0_rp,                         &
                   (-m1*m2*p2         )/4.0_rp),                        &
            vector((-p1*m2*p3+m1*m2*p3)/4.0_rp,                         &
                   (-m1*p1*p3         )/4.0_rp,                         &
                   ( m1*p1*m2         )/4.0_rp),                        &
            vector(( m2*p2*p3         )/4.0_rp,                         &
                   (-p1*p2*p3+p1*m2*p3)/4.0_rp,                         &
                   ( p1*m2*p2         )/4.0_rp),                        &
            vector((-p1*p2*p3+m1*p2*p3)/4.0_rp,                         &
                   ( m1*p1*p3         )/4.0_rp,                         &
                   ( m1*p1*p2         )/4.0_rp),                        &
            vector((-m2*p2*p3         )/4.0_rp,                         &
                   (-m1*p2*p3+m1*m2*p3)/4.0_rp,                         &
                   ( m1*m2*p2         )/4.0_rp),                        &
            vector((-m2*m3*p3         )/4.0_rp,                         &
                   (-m1*m3*p3         )/4.0_rp,                         &
                   (-m1*m2*p3+m1*m2*m3)/4.0_rp),                        &
            vector(( m2*m3*p3         )/4.0_rp,                         &
                   (-p1*m3*p3         )/4.0_rp,                         &
                   (-p1*m2*p3+p1*m2*m3)/4.0_rp),                        &
            vector(( p2*m3*p3         )/4.0_rp,                         &
                   ( p1*m3*p3         )/4.0_rp,                         &
                   (-p1*p2*p3+p1*p2*m3)/4.0_rp),                        &
            vector((-p2*m3*p3         )/4.0_rp,                         &
                   ( m1*m3*p3         )/4.0_rp,                         &
                   (-m1*p2*p3+m1*p2*m3)/4.0_rp)]
      end function deriv_shapefn
end module hexa20
