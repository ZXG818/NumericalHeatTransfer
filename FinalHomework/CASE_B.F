CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      MODULE USER_L
C******************************************************************
       INTEGER*4 I,J
       real(8) :: TIN       ! air inlet temperature
       real(8) :: TW        ! wall temperature
       real(8) :: UIN       ! air inlet velocity
       real(8) :: UOUT      ! air outlet velocity
       real(8) :: AMU       ! kinematic viscosity
       real(8) :: COND      ! conductivity
       real(8) :: FLOWIN    ! inlet flow rate calculation
       real(8) :: FL        ! outlet flow rate calculation
       real(8) :: FACTOR
       real(8) :: PR
C******************************************************************
      END MODULE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE USER
C*******************************************************************
      USE START_L
      USE USER_L
      IMPLICIT NONE
C******************************************************************
C                   CASE-A CALCULATION
C                ALL UNITS ARE MINIMETER
C******************************************************************

C
      ENTRY GRID
      TITLE(1)='  VEL_U'
      TITLE(2)='  VEL_V'
      TITLE(3)=' STR_FN'
      TITLE(4)=' TEMP. '
      TITLE(11)='PRESSURE'
      TITLE(12)='DENSITY'
      RELAX(1)=0.8             ! relaxtion
      RELAX(2)=0.8             ! relaxtion
      LSOLVE(1)=.TRUE.
      LSOLVE(4)=.TRUE.
      LPRINT(1)=.TRUE.
      LPRINT(2)=.TRUE.
      LPRINT(3)=.TRUE.
      LPRINT(4)=.TRUE.
      LPRINT(11)=.TRUE.
      LPRINT(12)=.TRUE.
      LAST=500
      XL=25.0E-3                  ! 20mm
      YL=0.65E-3                  ! half of the tube in radial direction, use the symmetry boundary.
      L1=100
      M1=65
      CALL UGRID
      RETURN
C
      ENTRY START
      TIN=20.0          ! inlet temperature
      TW=40.0           ! wall temperature
      UIN=30.0          ! the velocity inlet
      UOUT=UIN*YCV(M2)/Y(M1)*TW/TIN        ! estimate the outlet U velocity
      DO J=1,M1
          DO I=1,L1
              V(I,J)=0.0   ! initailize the V velocity field
              U(I,J)=UOUT  ! initialize the U velocity field
              U(I,M1)=0.0  ! upper wall, the U is zero
              T(I,J)=TIN    ! initialize the temperature field
          ENDDO
      ENDDO
      ! THE INLET FIN
      DO J=1,M1
          DO I=1,L1
              ! the length of the upper fin is 4mm
              IF((I>=1.AND.I<=16).AND.(J>=55))THEN
                  U(I,J)=0.
                  V(I,J)=0.
                  T(I,M1)=TW
              ENDIF
          ENDDO
      ENDDO
      ! THE FIRST BOTTOM FIN 
      DO J=1, M1
          DO I=1,L1
              ! the length of the bottom fin is 2mm
              IF((I>=16.AND.I<=24).AND.(J<=10))THEN
                  U(I,J)=0.
                  V(I,J)=0.
                  T(I,1)=TW
              ENDIF
          ENDDO
      ENDDO
      ! THE SECOND UPPER FIN
      DO J=1,M1
          DO I=1,L1
              ! the length of the upper fin is 4mm
              IF((I>=24.AND.I<=40).AND.(J>=55))THEN
                  U(I,J)=0.
                  V(I,J)=0.
                  T(I,M1)=TW
              ENDIF
          ENDDO
      ENDDO
      ! THE SECOND BOTTOM FIN
      DO J=1,M1
          DO I=1,L1
              ! the length of the bottom fin is 2mm
              IF((I>=40.AND.I<=48).AND.(J<=10))THEN
                  U(I,J)=0.
                  V(I,J)=0.
                  T(I,1)=TW
              ENDIF
          ENDDO
      ENDDO
      ! THE THIRD UPPER FIN
      DO J=1,M1
          DO I=1,L1
              ! the length of the upper fin is 4mm
              IF((I>=48.AND.I<=64).AND.(J>=55))THEN
                  U(I,J)=0.
                  V(I,J)=0.
                  T(I,M1)=TW
              ENDIF
          ENDDO
      ENDDO
      ! inlet velocity
      DO J=1,M1
          DO I=1,L1
              IF(J<=55)THEN
                  U(2,J)=UIN
                  T(1,J)=TIN
              ENDIF
          ENDDO
      ENDDO
      U(2,1:M1) = UIN      ! U velocity inlet boundary
      T(1,2:M2) = TIN      ! inlet temperature
      PR=0.700             ! Prandtl number
      AMU=1.7894E-5
      CPCON=1006.43E-3
      COND=AMU*CPCON/PR    ! coefficient of the thermal conductivity
      RETURN
      
C     Density
      ENTRY DENSE
      DO J=1,M1
          DO I=1, L1
              RHO(I,J)=1.225
          ENDDO
      ENDDO
      RETURN
      
C     Boundary condition
      ENTRY BOUND          ! Needs to be met the MASS CONSERVATION.
      ! THE INLET FIN
      DO J=1,M1
          DO I=1,L1
              ! the length of the upper fin is 4mm
              IF((I>=1.AND.I<=16).AND.(J>=55))THEN
                  U(I,J)=0.
                  V(I,J)=0.
                  T(I,1)=TW
              ENDIF
          ENDDO
      ENDDO
      ! THE FIRST BOTTOM FIN
      DO J=1, M1
          DO I=1,L1
              ! the length of the bottom fin is 2mm
              IF((I>=16.AND.I<=24).AND.(J<=10))THEN
                  U(I,J)=0.
                  V(I,J)=0.
                  T(I,M1)=TW
              ENDIF
          ENDDO
      ENDDO
      ! THE SECOND UPPER FIN
      DO J=1,M1
          DO I=1,L1
              ! the length of the upper fin is 4mm
              IF((I>=24.AND.I<=40).AND.(J>=55))THEN
                  U(I,J)=0.
                  V(I,J)=0.
                  T(I,M1)=TW
              ENDIF
          ENDDO
      ENDDO
      ! THE SECOND BOTTOM FIN
      DO J=1,M1
          DO I=1,L1
              ! the length of the bottom fin is 2mm
              IF((I>=40.AND.I<=48).AND.(J<=10))THEN
                  U(I,J)=0.
                  V(I,J)=0.
                  T(I,1)=TW
              ENDIF
          ENDDO
      ENDDO
      ! THE THIRD UPPER FIN
      DO J=1,M1
          DO I=1,L1
              ! the length of the upper fin is 4mm
              IF((I>=48.AND.I<=64).AND.(J>=55))THEN
                  U(I,J)=0.
                  V(I,J)=0.
                  T(I,M1)=TW
              ENDIF
          ENDDO
      ENDDO
      ! GLOBAL MASS CONSERVATION
      IF(ITER==0) THEN
          FLOWIN=0.0 ! inlet flow rate calculation
          DO J=2,M2
              IF(J<=55)THEN
                  FLOWIN=FLOWIN+RHO(1,1)*U(2,J)*YCV(M2)  
              ENDIF
          ENDDO
      ENDIF
      FL=0.0 ! outlet flow rate calculation
      DO J=2,M2
          FL=FL+RHO(1,1)*U(L2,J)*YCV(M2)
      ENDDO
      
      ! get the factor
      FACTOR=FLOWIN/FL
      DO J=2,M2
          U(L1,J)=U(L2,J)*FACTOR     ! local one way assumption.
          T(L1,J)=T(L2,J)            ! gradient Temperature equal 0.
      ENDDO
      RETURN
C
      ENTRY OUTPUT
      IF(ITER==0) THEN
      PRINT 401
      WRITE(8,401)
  401 FORMAT('   ITER',7X,'SMAX',11X,'SSUM',10X,'V(4,7)',
     1 9X,'T(4,7)')
      ELSE
      PRINT 403,ITER,SMAX,SSUM,V(4,7),T(4,7)
      WRITE(8,403) ITER,SMAX,SSUM,V(4,7),T(4,7)
  403 FORMAT(I6,1P4E15.3)
      ENDIF
      IF(ITER==LAST) THEN
C     only for output, no computing use.
      DO J=1,M1
          V(L1,J)=V(L2,J)
          T(L1,J)=T(L2,J)
          RHO(L1,J)=RHO(L2,J)
      ENDDO
      CALL PRINT
      ENDIF
      RETURN
C
      ENTRY GAMSOR
      DO J=1,M1
          DO I=1,L1
              ! MOMENTUM EQUATION
              IF(NF/=4)THEN
                  GAM(I,J)=AMU
                  ! THE INLET FIN
                  IF((I>=1.AND.I<=16).AND.(J>=55))THEN
                      GAM(I,J)=1.0E20
                  ENDIF
                  ! THE FIRST BOTTOM FIN
                  IF((I>=16.AND.I<=24).AND.(J<=10))THEN
                      GAM(I,J)=1.0E20
                  ENDIF
                  ! THE SECOND UPPER FIN
                  IF((I>=24.AND.I<=40).AND.(J>=55))THEN
                      GAM(I,J)=1.0E20
                  ENDIF
                  ! THE SECOND BOTTOM FIN
                  IF((I>=40.AND.I<=48).AND.(J<=10))THEN
                      GAM(I,J)=1.0E20
                  ENDIF
                  ! THE THIRD UPPER FIN
                  IF((I>=48.AND.I<=64).AND.(J>=55))THEN
                      GAM(I,J)=1.0E20
                  ENDIF
                  GAM(L1,J)=0.                      ! LOCAL ONE WAY
                  IF(NF/=2) GAM(I,1)=0.             ! BOTTOM IS SYMMETRY
                  IF(NF/=2) GAM(I,M1)=0.            ! UPPER IS SUMMETRY
              ENDIF
              ! ENERGY EQUATION
              IF(NF==4)THEN
                  GAM(I,J)=COND
                  GAM(I,1)=0.                       ! BOTTOM IS SYMMETRY
                  GAM(I,M1)=0.                      ! UPPER IS SYMMETRY
                  ! THE INLET FIN
                  IF((I>=1.AND.I<=16).AND.(J>=55))THEN
                      GAM(I,J)=1.0E20
                  ENDIF
                  ! THE FIRST BOTTOM FIN
                  IF((I>=16.AND.I<=24).AND.(J<=10))THEN
                      GAM(I,J)=1.0E20
                  ENDIF
                  ! THE SECOND UPPER FIN
                  IF((I>=24.AND.I<=40).AND.(J>=55))THEN
                      GAM(I,J)=1.0E20
                  ENDIF
                  ! THE SECOND BOTTOM FIN
                  IF((I>=40.AND.I<=48).AND.(J<=10))THEN
                      GAM(I,J)=1.0E20
                  ENDIF
                  ! THE THIRD UPPER FIN
                  IF((I>=48.AND.I<=64).AND.(J>=55))THEN
                      GAM(I,J)=1.0E20
                  ENDIF
                  GAM(L1,J)=0.
              ENDIF
          ENDDO
      ENDDO
      RETURN
      END