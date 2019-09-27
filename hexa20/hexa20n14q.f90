module hexa20n14q
   use iso_fortran_env, only: rp => REAL64
   implicit none

   integer,parameter :: nsd = 3, ntd = 6, nen = 20, nqd = 14
   real(kind=rp),parameter ::
      phi(nqd,nen) = RESHAPE([
     a**2/8._rp+a/8._rp+(-1._rp)/4._rp
     a**2/8._rp-a/8._rp+(-1._rp)/4._rp
     a**2/8._rp-a/8._rp+(-1._rp)/4._rp
     a**2/8._rp+a/8._rp+(-1._rp)/4._rp
     a**2/8._rp+a/8._rp+(-1._rp)/4._rp
     a**2/8._rp-a/8._rp+(-1._rp)/4._rp
     a**2/8._rp-a/8._rp+(-1._rp)/4._rp
     a**2/8._rp+a/8._rp+(-1._rp)/4._rp
     1._rp/4._rp-a**2/4._rp
     a/4._rp+1._rp/4._rp
     1._rp/4._rp-a**2/4._rp
     1._rp/4._rp-a/4._rp
     1._rp/4._rp-a**2/4._rp
     a/4._rp+1._rp/4._rp
     1._rp/4._rp-a**2/4._rp
     1._rp/4._rp-a/4._rp
     1._rp/4._rp-a/4._rp
     a/4._rp+1._rp/4._rp
     a/4._rp+1._rp/4._rp
     1._rp/4._rp-a/4._rp

      ],[nqd,nen])
      phi_xi(nsd,nen,nqd) =
end module hexa20n14q
