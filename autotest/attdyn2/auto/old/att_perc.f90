      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

! Evaluates the algebraic equations or ODE right hand side

! Input arguments :
!      NDIM   :   Dimension of the ODE system 
!      U      :   State variables
!      ICP    :   Array indicating the free parameter(s)
!      PAR    :   Equation parameters

! Values to be returned :
!      F      :   ODE right hand side values

! Normally unused Jacobian arguments : IJAC, DFDU, DFDP (see manual)

      IMPLICIT NONE
      INTEGER NDIM, IJAC, ICP(*)
      DOUBLE PRECISION U(NDIM), PAR(*), F(NDIM), DFDU(*), DFDP(*)
      DOUBLE PRECISION U1,U2,U3
      DOUBLE PRECISION A,B,C,D
      DOUBLE PRECISION D1,D2,SIGMA

      A=-1
      B=-2
      C=1
      D=1

      SIGMA=PAR(1)*((sqrt(A*D-B*C)-sqrt(-B*C))/A)**2
      D1=( (SIGMA*A+D)/(4*SIGMA) )
      D2=D1*SIGMA
      
      U1=U(1)
      U2=U(2)
      U3=U(3)

      F(1) = A*U1 +B*U2 -U1*U2**2 -2*D1*U1*U3
      F(2) = C*U1 +D*U2 +U1*U2**2 -2*D2*U2*U3
      F(3) = PAR(2)*U3*(1-U3)*(PAR(3)-U1**2)

      END SUBROUTINE FUNC
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

! Input arguments :
!      NDIM   :   Dimension of the ODE system 

! Values to be returned :
!      U      :   A starting solution vector
!      PAR    :   The corresponding equation-parameter values
!      T      :	  Not used here

      IMPLICIT NONE
      INTEGER NDIM
      DOUBLE PRECISION U(NDIM), PAR(*), T

! Initialize the equation parameters
       PAR(1)=0.95
       PAR(2)=4.0
       PAR(3)=0.01

! Initialize the solution
!       U(1)=-0.100000000000000
!       U(2)=+0.141691283930165
!       U(3)=+0.662432168104518

       U(1)=-0.100000000000000
       U(2)=+0.153904428850466
       U(3)=+0.800127522126962

      END SUBROUTINE STPNT
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! The following subroutines are not used here,
! but they must be supplied as dummy routines

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
!----------------------------------------------------------------------
!----------------------------------------------------------------------
