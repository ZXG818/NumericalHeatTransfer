CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      MODULE USER_L
C******************************************************************
       INTEGER*4 I,J
       REAL(8) :: TIN       ! air inlet temperature
       REAL(8) :: TW        ! wall temperature
       REAL(8) :: UIN       ! air inlet velocity
       REAL(8) :: UOUT      ! air outlet velocity
       REAL(8) :: AMU       ! kinematic viscosity
       REAL(8) :: COND      ! conductivity
       REAL(8) :: FLOWIN    ! inlet flow rate calculation
       REAL(8) :: FL        ! outlet flow rate calculation
       REAL(8) :: FACTOR
       REAL(8) :: PR
       REAL(8) :: PPin      ! calculate the pumping power 
       REAL(8) :: PPout
       CHARACTER(LEN=20) :: FILENAME='PUMPINGPOWER' ! PP FILENAME
       INTEGER(4) :: FILEID=100
       REAL(8) :: Tout
       REAL(8) :: Tave
       REAL(8) :: Hin
       REAL(8) :: Hout
       REAL(8) :: Have ! the average convective heat transfer coefficient
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
      UIN=1.5          ! the velocity inlet
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
      CPCON=1006.43
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
      ! write down the pumping power
      IF(ITER==LAST)THEN
          PPin=0.
          PPout=0.
          ! calculate the inner PP
          DO J=1,M1
              IF(J<55)THEN
                  PPin = PPin + RHO(1,J)*U(2,J)*YCV(J)*P(1,J)
              ENDIF
          ENDDO
          ! calculate the outter PP
          DO J=1,M1
              PPout = PPout + RHO(L1,J)*U(L1,J)*YCV(J)*P(L1,J)
          ENDDO
          ! write down to the file
          ! open the file
          OPEN(UNIT=FILEID, FILE=FILENAME)
          ! write the title
          WRITE(FILEID, *)
     1    " INLET P    INLET U    OUTLET P    OUTLET U    "
          DO J=1,M1
              IF(J<55)THEN
                  WRITE(FILEID, "(F8.4,2X,F8.4,4X,F8.2,2X,F8.4)")
     1            P(1,J), U(2,J), P(L1,J), U(L1,J)
              ELSE
                  WRITE(FILEID, "(22X, F8.2,2X,F8.4)") P(L1,J),U(L1,J)
              ENDIF
          ENDDO
          ! write the pumping power
          WRITE(FILEID, *)" INLET PP  OUTLET PP  PUMPING POWER"
          WRITE(FILEID,"(F8.4,4X,F8.4,4X,F8.4)") PPin, PPout, PPin-PPout
          ! calculate the heat transfer rate
          ! first we should calculate the outlet enthalpy
          WRITE(FILEID, *)" **************************************** "
          WRITE(FILEID, *)" *         HEAT TRANSFER RATE           * "
          WRITE(FILEID, *)" **************************************** "
          Hin=0.
          Hout=0.
          Tout=0.
          DO J=1,M1
              IF(J<55) THEN
                  Hin=Hin+CPCON*RHO(1,J)*T(1,J)*RHO(1,J)*U(2,J)*YCV(J)
              ENDIF
              Hout=Hout+CPCON*RHO(L1,J)*T(L1,J)*RHO(L1,J)*U(L1,J)*YCV(J)
              Tout=Tout+T(L1,J)
          ENDDO
          Tave=Tout/M1
          Have=(Hout-Hin)/((4.2*3 + 2.2*2)*1E-3)/(TAVE-TIN)
          WRITE(FILEID, *)"ENTHALPY_IN ENTHALPY_OUT H_AVERAGE T_AVERAGE"
          WRITE(FILEID, "(4F10.4)")HIN,HOUT,HAVE,TAVE
          WRITE(FILEID, *)"THE FINAL HEAT TRANSFER RATE IS:",
     1                    HAVE*(TAVE-TIN)*17.0E-3
          WRITE(FILEID, *)"THE FINAL HEAT FLUX IS :", HAVE*(TAVE-TIN)
          ! close the file
          CLOSE(FILEID)
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
                  GAM(L1,J)=0.            ! LOCAL ONE WAY
                  IF(NF/=2) GAM(I,1)=0.   ! BOTTOM IS SYMMETRY, EXCEPT V
                  IF(NF/=2) GAM(I,M1)=0.  ! UPPER IS SUMMETRY, EXCEPT V
              ENDIF
              ! ENERGY EQUATION
              IF(NF==4)THEN
                  GAM(I,J)=COND/CPCON
                  GAM(I,1)=0.             ! BOTTOM IS SYMMETRY
                  GAM(I,M1)=0.            ! UPPER IS SYMMETRY
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