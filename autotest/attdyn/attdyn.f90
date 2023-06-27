SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, IJAC, ICP(*)
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
  DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
  DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM),DFDP(NDIM,*)

  DOUBLE PRECISION x, y, A, a1, b, c, d, alpha, epsilon, theta, sigma, sigma_c, D1, D2

  x = U(1)
  y = U(2)
  A = U(3)

  a1 = PAR(1)
  b = PAR(2)
  c = PAR(3)
  d = PAR(4)
  alpha = PAR(5)
  epsilon = PAR(6)
  theta = PAR(7)
  
  sigma_c = (a1*d - 2*b*c + 2*sqrt(b**2*c**2 - a1*b*c*d)) / (a1**2)
  sigma = theta*sigma_c

  D1 = (sigma*a1 + d) / (4*sigma)
  D2 = sigma*D1

  F(1) = a1*x + b*y - x*y**2 - 2*D1*A*x
  F(2) = c*x + d*y + x*y**2 - 2*D2*A*y
  F(3) = alpha*A*(1 - A)*(epsilon - x)*(epsilon + x)

END SUBROUTINE FUNC

SUBROUTINE STPNT(NDIM,U,PAR,T)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
  DOUBLE PRECISION, INTENT(IN) :: T

  PAR(1:7) = (/ -1.1d0, -2.0d0, 1.0d0, 1.0d0, 1.0d0, 0.01d0, 0.3d0 /)

  U(1) = 0
  U(2) = 0
  U(3) = 1
   
END SUBROUTINE STPNT

SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
  DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
  DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
  DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)

!X FB(1)=
!X FB(2)=

END SUBROUTINE BCND

SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
  DOUBLE PRECISION, INTENT(IN) :: PAR(*)
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
  DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
  DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)

!X FI(1)=

END SUBROUTINE ICND

SUBROUTINE FOPT(NDIM,U,ICP,PAR,IJAC,FS,DFDU,DFDP)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
  DOUBLE PRECISION, INTENT(OUT) :: FS
  DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM),DFDP(*)

!X FS=

END SUBROUTINE FOPT

SUBROUTINE PVLS(NDIM,U,PAR)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
  DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)

END SUBROUTINE PVLS