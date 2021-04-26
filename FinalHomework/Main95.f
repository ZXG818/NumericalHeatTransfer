CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C This computer program was copied from the graduate student course
C program of the University of Minnesota.  Part of it was re-formu-
C lated to meet the microcomputer environment.  Some inappropriate
C expressions were also corrected. The program is used only for the
C teaching purpose.  No part of it may be published. You may use it
C as a frame to re-develop your own code for research purpose.
C --------Instructor of Numerical Heat Transfer, XJTU,1998.12-------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C The current version of the programm was updated from Fortran 77 to Fortran 95  
C by Dr.YuTong Mu, Dr. Li Chen,  Dr. Kong Lin of NHT group of XJTU in 2013.01-04, 
C contact us with muyutong@stu.xjtu.edu.cn  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C*******************************************************************
      MODULE START_L
       PARAMETER(NI=100,NJ=200,NIJ=NI,NFMAX=10,NFX4=NFMAX+4)
C******************************************************************
      CHARACTER*8 TITLE(NFX4)
      LOGICAL LSOLVE(NFX4),LPRINT(NFX4),LBLK(NFX4),LSTOP
      REAL*8,DIMENSION(NI,NJ,NFX4)::F
      REAL*8,DIMENSION(NI,NJ,6)::COF,COFU,COFV,COFP
      REAL*8,DIMENSION(NI,NJ)::P,RHO,GAM,CP,CON,AIP,AIM,AJP,AJM,AP
      REAL*8,DIMENSION(NI,NJ)::U,V,PC,T,DU,DV,UHAT,VHAT 
      REAL*8,DIMENSION(NI)::X,XU,XDIF,XCV,XCVS,XCVI,XCVIP
      REAL*8,DIMENSION(NJ)::Y,YV,YDIF,YCV,YCVS,YCVR,YCVRS,ARX,ARXJ,
     1    ARXJP,R,RMN,SX,SXMN
      REAL*8,DIMENSION(NI)::FV,FVP,FX,FXM
      REAL*8,DIMENSION(NJ)::FY,FYM
      REAL*8,DIMENSION(NIJ)::PT,QT
      REAL*8 RELAX(NFX4),TIME,DT,XL,YL,RHOCON,CPCON
      INTEGER*4 NF,NP,NRHO,NGAM,NCP,L1,L2,L3,M1,M2,M3,
     1    IST,JST,ITER,LAST,MODE,NTIMES(NFX4),IPREF,JPREF
      REAL*8 SMAX,SSUM
      REAL*8 FLOW,DIFF,ACOF
C******************************************  
      EQUIVALENCE(F(1,1,1),U(1,1)),(F(1,1,2),V(1,1)),(F(1,1,3),PC(1,1))
     1,(F(1,1,4),T(1,1)) 
      EQUIVALENCE(F(1,1,11),P(1,1)),(F(1,1,12),RHO(1,1)),(F(1,1,13)
     1,GAM(1,1)),(F(1,1,14),CP(1,1))
      EQUIVALENCE(COF(1,1,1),CON(1,1)),(COF(1,1,2),AIP(1,1)),
     1(COF(1,1,3),AIM(1,1)),(COF(1,1,4),AJP(1,1)),
     2(COF(1,1,5),AJM(1,1)),(COF(1,1,6),AP(1,1))
      REAL*8,DIMENSION(NI)::TH,THU,THDIF,THCV,THCVS
      REAL*8 THL
      EQUIVALENCE(X,TH),(XU,THU),(XDIF,THDIF),(XCV,THCV),
     1(XCVS,THCVS),(XL,THL)
      DATA LSTOP,LSOLVE,LPRINT/.FALSE.,NFX4*.FALSE.,NFX4*.FALSE./
      DATA LBLK/NFX4*.TRUE./
      DATA MODE,LAST,TIME,ITER/1,5,0.,0/
      DATA RELAX,NTIMES/NFX4*1.,NFX4*1/
      DATA DT,IPREF,JPREF,RHOCON,CPCON/1.E+30,1,1,1.,1./
      
      REAL*8,DIMENSION(NI,NJ)::W
      EQUIVALENCE (F(1,1,5),W(1,1))
  
      END MODULE
C -------------------------MAIN PROGRAM----------------------------
C*******************************************************************
      PROGRAM MAIN
      USE START_L
      IMPLICIT NONE
C********************************************************************
      OPEN(8,FILE='RESULT.TXT')
      CALL GRID
      CALL SETUP1
      CALL START
      DO WHILE(.NOT.LSTOP)
      CALL DENSE
      CALL BOUND
      CALL OUTPUT
      CALL SETUP2
      ENDDO
      CALL OUTPUT
      CLOSE(8)
      STOP
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DIFLOW
C********************************************************************
      USE START_L
      IMPLICIT NONE
      REAL*8 TEMP
C********************************************************************
C******************************************************************
      ACOF=DIFF
      IF(FLOW==0.) RETURN
      TEMP=DIFF-ABS(FLOW)*0.1
      ACOF=0.
      IF(TEMP<=0.) RETURN
      TEMP=TEMP/DIFF
      ACOF=DIFF*TEMP**5
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SOLVE
      USE START_L
      IMPLICIT NONE
      INTEGER*4 ISTF,JSTF,IT1,IT2,JT1,JT2,NT,N,I,J,II,JJ
      REAL*8 BL,BLP,BLM,BLC,DENOM,TEMP
C******************************************************************
      ISTF=IST-1
      JSTF=JST-1
      IT1=L2+IST
      IT2=L3+IST
      JT1=M2+JST
      JT2=M3+JST
C******************************************************************
      DO 999 NT=1,NTIMES(NF)
      N=NF
C-------------------------------------------------------------------
      IF(LBLK(NF)) THEN
      PT(ISTF)=0.
      QT(ISTF)=0.
      DO 11 I=IST,L2
      BL=0.
      BLP=0.
      BLM=0.
      BLC=0.
      DO  12 J=JST,M2
      BL=BL+AP(I,J)
      IF(J/=M2)  BL=BL-AJP(I,J)
      IF(J/=JST) BL=BL-AJM(I,J)
      BLP=BLP+AIP(I,J)
      BLM=BLM+AIM(I,J)
      BLC=BLC+CON(I,J)+AIP(I,J)*F(I+1,J,N)+AIM(I,J)*F(I-1,J,N)
     1   +AJP(I,J)*F(I,J+1,N)+AJM(I,J)*F(I,J-1,N)-AP(I,J)*F(I,J,N)
   12 ENDDO
      DENOM=BL-PT(I-1)*BLM
      IF(ABS(DENOM/BL)<1.E-10) DENOM=1.E25
      PT(I)=BLP/DENOM
      QT(I)=(BLC+BLM*QT(I-1))/DENOM
   11 ENDDO
      BL=0.
      DO 13 II=IST,L2
      I=IT1-II
      BL=BL*PT(I)+QT(I)
      DO 14 J=JST,M2
      F(I,J,N)=F(I,J,N)+BL
   14 ENDDO
   13 ENDDO
C------------------------------------------------------------------------
      PT(JSTF)=0.
      QT(JSTF)=0.
      DO 21 J=JST,M2
      BL=0.
      BLP=0.
      BLM=0.
      BLC=0.
      DO 22  I=IST,L2
      BL=BL+AP(I,J)
      IF(I/=L2) BL=BL-AIP(I,J)
      IF(I/=IST) BL=BL-AIM(I,J)
      BLP=BLP+AJP(I,J)
      BLM=BLM+AJM(I,J)
      BLC=BLC+CON(I,J)+AIP(I,J)*F(I+1,J,N)+AIM(I,J)*F(I-1,J,N)
     1   +AJP(I,J)*F(I,J+1,N)+AJM(I,J)*F(I,J-1,N)-AP(I,J)*F(I,J,N)
   22 ENDDO
      DENOM=BL-PT(J-1)*BLM
      IF(ABS(DENOM/BL)<1.E-10) DENOM=1.E25
      PT(J)=BLP/DENOM
      QT(J)=(BLC+BLM*QT(J-1))/DENOM
   21 ENDDO
      BL=0.
      DO  23 JJ=JST,M2
      J=JT1-JJ
      BL=BL*PT(J)+QT(J)
      DO 24 I=IST,L2
      F(I,J,N)=F(I,J,N)+BL
   24 ENDDO
   23 ENDDO
   10 ENDIF
C--------------------------------------------------------------
      DO 90 J=JST,M2
      PT(ISTF)=0.
      QT(ISTF)=F(ISTF,J,N)
      DO 70 I=IST,L2
      DENOM=AP(I,J)-PT(I-1)*AIM(I,J)
      PT(I)=AIP(I,J)/DENOM
      TEMP=CON(I,J)+AJP(I,J)*F(I,J+1,N)+AJM(I,J)*F(I,J-1,N)
      QT(I)=(TEMP+AIM(I,J)*QT(I-1))/DENOM
   70 ENDDO
      DO 80 II=IST,L2
      I=IT1-II
      F(I,J,N)=F(I+1,J,N)*PT(I)+QT(I)
   80 ENDDO
   90 ENDDO
C------------------------------------------------------------------
      DO 190 JJ=JST,M3
      J=JT2-JJ
      PT(ISTF)=0.
      QT(ISTF)=F(ISTF,J,N)
      DO 170 I=IST,L2
      DENOM=AP(I,J)-PT(I-1)*AIM(I,J)
      PT(I)=AIP(I,J)/DENOM
      TEMP=CON(I,J)+AJP(I,J)*F(I,J+1,N)+AJM(I,J)*F(I,J-1,N)
      QT(I)=(TEMP+AIM(I,J)*QT(I-1))/DENOM
  170 ENDDO
      DO 180 II=IST,L2
      I=IT1-II
      F(I,J,N)=F(I+1,J,N)*PT(I)+QT(I)
  180 ENDDO
  190 ENDDO
C------------------------------------------------------------------
      DO 290 I=IST,L2
      PT(JSTF)=0.
      QT(JSTF)=F(I,JSTF,N)
      DO 270 J=JST,M2
      DENOM=AP(I,J)-PT(J-1)*AJM(I,J)
      PT(J)=AJP(I,J)/DENOM
      TEMP=CON(I,J)+AIP(I,J)*F(I+1,J,N)+AIM(I,J)*F(I-1,J,N)
      QT(J)=(TEMP+AJM(I,J)*QT(J-1))/DENOM
  270 ENDDO
      DO 280 JJ=JST,M2
      J=JT1-JJ
      F(I,J,N)=F(I,J+1,N)*PT(J)+QT(J)
  280 ENDDO
  290 ENDDO
C-------------------------------------------------------------------
      DO 390 II=IST,L3
      I=IT2-II
      PT(JSTF)=0.
      QT(JSTF)=F(I,JSTF,N)
      DO 370 J=JST,M2
      DENOM=AP(I,J)-PT(J-1)*AJM(I,J)
      PT(J)=AJP(I,J)/DENOM
      TEMP=CON(I,J)+AIP(I,J)*F(I+1,J,N)+AIM(I,J)*F(I-1,J,N)
      QT(J)=(TEMP+AJM(I,J)*QT(J-1))/DENOM
  370 ENDDO
      DO 380 JJ=JST,M2
      J=JT1-JJ
      F(I,J,N)=F(I,J+1,N)*PT(J)+QT(J)
  380 ENDDO
  390 ENDDO
C*******************************************
  999 ENDDO
      ENTRY RESET
      DO 400 J=2,M2
      DO 401 I=2,L2
      CON(I,J)=0.
      AP(I,J)=0.
  401 ENDDO
  400 ENDDO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SETUP
C******************************************************************
      USE START_L
      IMPLICIT NONE 
      INTEGER*4 I,J,K,N
      REAL*8 REL,FL,FLM,FLP,GM,GMM,VOL,APT,AREA,SXT,SXB,ARHO
C******************************************************************
    1 FORMAT(//15X,'COMPUTATION  IN  CARTESIAN  COORDINATES')
    2 FORMAT(//15X,'COMPUTATION FOR AXISYMMETRIC SITUATION')
    3 FORMAT(//15X,'COMPUTATION   IN   POLAR   COORDINATES')
    4 FORMAT(14X,40(1H*),//)
C-------------------------------------------------------------------
      ENTRY SETUP1
      NP=NFMAX+1
      NRHO=NP+1
      NGAM=NRHO+1
      NCP=NGAM+1
      L2=L1-1
      L3=L2-1
      M2=M1-1
      M3=M2-1
      X(1)=XU(2)
      DO 5 I=2,L2
      X(I)=0.5*(XU(I+1)+XU(I))
    5 ENDDO
      X(L1)=XU(L1)
      Y(1)=YV(2)
      DO 10 J=2,M2
      Y(J)=0.5*(YV(J+1)+YV(J))
   10 ENDDO
      Y(M1)=YV(M1)
      DO 15 I=2,L1
      XDIF(I)=X(I)-X(I-1)
   15 ENDDO
      DO 18 I=2,L2
      XCV(I)=XU(I+1)-XU(I)
   18 ENDDO
      DO 20 I=3,L2
      XCVS(I)=XDIF(I)
   20 ENDDO
      XCVS(3)=XCVS(3)+XDIF(2)
      XCVS(L2)=XCVS(L2)+XDIF(L1)
      DO 22 I=3,L3
      XCVI(I)=0.5*XCV(I)
      XCVIP(I)=XCVI(I)
   22 ENDDO
      XCVIP(2)=XCV(2)
      XCVI(L2)=XCV(L2)
      DO 35 J=2,M1
      YDIF(J)=Y(J)-Y(J-1)
   35 ENDDO
      DO 40 J=2,M2
      YCV(J)=YV(J+1)-YV(J)
   40 ENDDO
      DO 45 J=3,M2
      YCVS(J)=YDIF(J)
   45 ENDDO
      YCVS(3)=YCVS(3)+YDIF(2)
      YCVS(M2)=YCVS(M2)+YDIF(M1)
      IF(MODE==1) THEN
      DO 52 J=1,M1
      RMN(J)=1.0
      R(J)=1.0
   52 ENDDO
      ELSE
      DO 50 J=2,M1
      R(J)=R(J-1)+YDIF(J)
   50 ENDDO
      RMN(2)=R(1)
      DO 60 J=3,M2
      RMN(J)=RMN(J-1)+YCV(J-1)
   60 ENDDO
      RMN(M1)=R(M1)
      ENDIF
      DO 57 J=1,M1
      SX(J)=1.
      SXMN(J)=1.
      IF(MODE==3) THEN
      SX(J)=R(J)
      IF(J/=1) SXMN(J)=RMN(J)
      ENDIF
   57 ENDDO
      DO 62 J=2,M2
      YCVR(J)=R(J)*YCV(J)
      ARX(J)=YCVR(J)
      IF(MODE==3) ARX(J)=YCV(J)
   62 ENDDO
      DO 64 J=4,M3
      YCVRS(J)=0.5*(R(J)+R(J-1))*YDIF(J)
   64 ENDDO
      YCVRS(3)=0.5*(R(3)+R(1))*YCVS(3)
      YCVRS(M2)=0.5*(R(M1)+R(M3))*YCVS(M2)
      IF(MODE==2) THEN
      DO 65 J=3,M3
      ARXJ(J)=0.25*(1.+RMN(J)/R(J))*ARX(J)
      ARXJP(J)=ARX(J)-ARXJ(J)
   65 ENDDO
      ELSE
      DO 66 J=3,M3
      ARXJ(J)=0.5*ARX(J)
      ARXJP(J)=ARXJ(J)
   66 ENDDO
      ENDIF
      ARXJP(2)=ARX(2)
      ARXJ(M2)=ARX(M2)
      DO 70 J=3,M3
      FV(J)=ARXJP(J)/ARX(J)
      FVP(J)=1.-FV(J)
   70 ENDDO
      DO 85 I=3,L2
      FX(I)=0.5*XCV(I-1)/XDIF(I)
      FXM(I)=1.-FX(I)
   85 ENDDO
      FX(2)=0.
      FXM(2)=1.
      FX(L1)=1.
      FXM(L1)=0.
      DO 90 J=3,M2
      FY(J)=0.5*YCV(J-1)/YDIF(J)
      FYM(J)=1.-FY(J)
   90 ENDDO
      FY(2)=0.
      FYM(2)=1.
      FY(M1)=1.
      FYM(M1)=0.
CON,AP,U,V,RHO,CP,PC AND P ARRAYS ARE INITIALIZED HERE
      DO 96 J=1,M1
      DO 95 I=1,L1
      PC(I,J)=0.
      U(I,J)=0.
      V(I,J)=0.
      CON(I,J)=0.
      AP(I,J)=0.
      RHO(I,J)=RHOCON
      CP(I,J)=CPCON
      P(I,J)=0.
   95 ENDDO
   96 ENDDO
      IF(MODE==1) PRINT 1
      IF(MODE==1) WRITE(8,1)
      IF(MODE==2) PRINT 2
      IF(MODE==2) WRITE(8,2)
      IF(MODE==3) PRINT 3
      IF(MODE==3) WRITE(8,3)
      PRINT 4
      WRITE(8,4)
      RETURN
C---------------------------------------------------------------
      ENTRY SETUP2
COEFFICIENTS FOR THE U EQUATION
      NF=1
      IF(LSOLVE(NF)) THEN
      IST=3
      JST=2
      CALL GAMSOR
      REL=1.-RELAX(NF)
      DO 102 I=3,L2
      FL=XCVI(I)*V(I,2)*RHO(I,1)
      FLM=XCVIP(I-1)*V(I-1,2)*RHO(I-1,1)
      FLOW=R(1)*(FL+FLM)
      DIFF=R(1)*(XCVI(I)*GAM(I,1)+XCVIP(I-1)*GAM(I-1,1))/YDIF(2)
      CALL DIFLOW
      AJM(I,2)=ACOF+MAX(0.,FLOW)
  102 ENDDO
      DO 103 J=2,M2
      FLOW=ARX(J)*U(2,J)*RHO(1,J)
      DIFF=ARX(J)*GAM(1,J)/(XCV(2)*SX(J))
      CALL DIFLOW
      AIM(3,J)=ACOF+MAX(0.,FLOW)
      DO 104 I=3,L2
      IF(I==L2) THEN
      FLOW=ARX(J)*U(L1,J)*RHO(L1,J)
      DIFF=ARX(J)*GAM(L1,J)/(XCV(L2)*SX(J))
      ELSE     
      FL=U(I,J)*(FX(I)*RHO(I,J)+FXM(I)*RHO(I-1,J))
      FLP=U(I+1,J)*(FX(I+1)*RHO(I+1,J)+FXM(I+1)*RHO(I,J))
      FLOW=ARX(J)*0.5*(FL+FLP)
      DIFF=ARX(J)*GAM(I,J)/(XCV(I)*SX(J))
      ENDIF
      CALL DIFLOW
      AIM(I+1,J)=ACOF+MAX(0.,FLOW)
      AIP(I,J)=AIM(I+1,J)-FLOW
      IF(J==M2) THEN
      FL=XCVI(I)*V(I,M1)*RHO(I,M1)
      FLM=XCVIP(I-1)*V(I-1,M1)*RHO(I-1,M1)
      DIFF=R(M1)*(XCVI(I)*GAM(I,M1)+XCVIP(I-1)*GAM(I-1,M1))/YDIF(M1)
      ELSE
      FL=XCVI(I)*V(I,J+1)*(FY(J+1)*RHO(I,J+1)+FYM(J+1)*RHO(I,J))
      FLM=XCVIP(I-1)*V(I-1,J+1)*(FY(J+1)*RHO(I-1,J+1)+FYM(J+1)*
     1 RHO(I-1,J))
      GM=GAM(I,J)*GAM(I,J+1)/(YCV(J)*GAM(I,J+1)+YCV(J+1)*GAM(I,J)+
     1 1.0E-30)*XCVI(I)
      GMM=GAM(I-1,J)*GAM(I-1,J+1)/(YCV(J)*GAM(I-1,J+1)+YCV(J+1)*
     1 GAM(I-1,J)+1.E-30)*XCVIP(I-1)
      DIFF=RMN(J+1)*2.*(GM+GMM)
      ENDIF
      FLOW=RMN(J+1)*(FL+FLM)
      CALL DIFLOW
      AJM(I,J+1)=ACOF+MAX(0.,FLOW)
      AJP(I,J)=AJM(I,J+1)-FLOW
      VOL=YCVR(J)*XCVS(I)
      APT=(RHO(I,J)*XCVI(I)+RHO(I-1,J)*XCVIP(I-1))
     1/(XCVS(I)*DT)
      AP(I,J)=AP(I,J)-APT
      CON(I,J)=CON(I,J)+APT*U(I,J)
      AP(I,J)=(-AP(I,J)*VOL+AIP(I,J)+AIM(I,J)+AJP(I,J)+AJM(I,J))
     1/RELAX(NF)
      CON(I,J)=CON(I,J)*VOL+REL*AP(I,J)*U(I,J)
      DU(I,J)=VOL/(XDIF(I)*SX(J))
      DU(I,J)=DU(I,J)/AP(I,J)
  104 ENDDO
  103 ENDDO
      COFU(IST:L2,JST:M2,1:6)=COF(IST:L2,JST:M2,1:6)
COEFFICIENTS FOR THE  V  EQUATION----------------------------------
      NF=2
      CALL RESET
      IST=2
      JST=3
      CALL GAMSOR
      REL=1.-RELAX(NF)
      DO 202 I=2,L2
      AREA=R(1)*XCV(I)
      FLOW=AREA*V(I,2)*RHO(I,1)
      DIFF=AREA*GAM(I,1)/YCV(2)
      CALL DIFLOW
      AJM(I,3)=ACOF+MAX(0.,FLOW)
  202 ENDDO
      DO 203 J=3,M2
      FL=ARXJ(J)*U(2,J)*RHO(1,J)
      FLM=ARXJP(J-1)*U(2,J-1)*RHO(1,J-1)
      FLOW=FL+FLM
      DIFF=(ARXJ(J)*GAM(1,J)+ARXJP(J-1)*GAM(1,J-1))/(XDIF(2)*SXMN(J))
      CALL DIFLOW
      AIM(2,J)=ACOF+MAX(0.,FLOW)
      DO 204 I=2,L2
      IF(I==L2) THEN
      FL=ARXJ(J)*U(L1,J)*RHO(L1,J)
      FLM=ARXJP(J-1)*U(L1,J-1)*RHO(L1,J-1)
      DIFF=(ARXJ(J)*GAM(L1,J)+ARXJP(J-1)*GAM(L1,J-1))/(XDIF(L1)*SXMN(J))
      ELSE
      FL=ARXJ(J)*U(I+1,J)*(FX(I+1)*RHO(I+1,J)+FXM(I+1)*RHO(I,J))
      FLM=ARXJP(J-1)*U(I+1,J-1)*(FX(I+1)*RHO(I+1,J-1)+FXM(I+1)*
     1 RHO(I,J-1))
      GM=GAM(I,J)*GAM(I+1,J)/(XCV(I)*GAM(I+1,J)+XCV(I+1)*GAM(I,J)+
     1 1.E-30)*ARXJ(J)
      GMM=GAM(I,J-1)*GAM(I+1,J-1)/(XCV(I)*GAM(I+1,J-1)+XCV(I+1)*
     1 GAM(I,J-1)+1.0E-30)*ARXJP(J-1)
      DIFF=2.*(GM+GMM)/SXMN(J)
      ENDIF
      FLOW=FL+FLM
      CALL DIFLOW
      AIM(I+1,J)=ACOF+MAX(0.,FLOW)
      AIP(I,J)=AIM(I+1,J)-FLOW
      IF(J==M2) THEN
      AREA=R(M1)*XCV(I)
      FLOW=AREA*V(I,M1)*RHO(I,M1)
      DIFF=AREA*GAM(I,M1)/YCV(M2)
      ELSE
      AREA=R(J)*XCV(I)
      FL=V(I,J)*(FY(J)*RHO(I,J)+FYM(J)*RHO(I,J-1))*RMN(J)
      FLP=V(I,J+1)*(FY(J+1)*RHO(I,J+1)+FYM(J+1)*RHO(I,J))*RMN(J+1)
      FLOW=(FV(J)*FL+FVP(J)*FLP)*XCV(I)
      DIFF=AREA*GAM(I,J)/YCV(J)
      ENDIF
      CALL DIFLOW
      AJM(I,J+1)=ACOF+MAX(0.,FLOW)
      AJP(I,J)=AJM(I,J+1)-FLOW
      VOL=YCVRS(J)*XCV(I)
      SXT=SX(J)
      IF(J==M2) SXT=SX(M1)
      SXB=SX(J-1)
      IF(J==3) SXB=SX(1)
      APT=(ARXJ(J)*RHO(I,J)*0.5*(SXT+SXMN(J))+ARXJP(J-1)*RHO(I,J-1)*
     10.5*(SXB+SXMN(J)))/(YCVRS(J)*DT)
      AP(I,J)=AP(I,J)-APT
      CON(I,J)=CON(I,J)+APT*V(I,J)
      AP(I,J)=(-AP(I,J)*VOL+AIP(I,J)+AIM(I,J)+AJP(I,J)+AJM(I,J))
     1/RELAX(NF)
      CON(I,J)=CON(I,J)*VOL+REL*AP(I,J)*V(I,J)
      DV(I,J)=VOL/YDIF(J)
      DV(I,J)=DV(I,J)/AP(I,J)
  204 ENDDO
  203 ENDDO
      COFV(IST:L2,JST:M2,1:6)=COF(IST:L2,JST:M2,1:6)
CALCULATE UHAT AND VHAT
      DO 150 J=2,M2
      DO 151 I=3,L2
      UHAT(I,J)=(COFU(I,J,2)*U(I+1,J)+COFU(I,J,3)*U(I-1,J)+COFU(I,J,4)
     1 *U(I,J+1)+COFU(I,J,5)*U(I,J-1)+COFU(I,J,1))/COFU(I,J,6)
  151 ENDDO
  150 ENDDO
      DO 250 J=3,M2
      DO 251 I=2,L2
      VHAT(I,J)=(COFV(I,J,2)*V(I+1,J)+COFV(I,J,3)*V(I-1,J)+COFV(I,J,4)
     1 *V(I,J+1)+COFV(I,J,5)*V(I,J-1)+COFV(I,J,1))/COFV(I,J,6)
  251 ENDDO
  250 ENDDO
COEFFICIENTS FOR THE PRESSURE EQUATION-------------------
      NF=3
      CALL RESET
      IST=2
      JST=2
      CALL GAMSOR
      DO 410 J=2,M2
      DO 411 I=2,L2
      VOL=YCVR(J)*XCV(I)
      CON(I,J)=CON(I,J)*VOL
  411 ENDDO
  410 ENDDO
      DO 402 I=2,L2
      ARHO=R(1)*XCV(I)*RHO(I,1)
      CON(I,2)=CON(I,2)+ARHO*V(I,2)
      AJM(I,2)=0.
  402 ENDDO
      DO 403 J=2,M2
      ARHO=ARX(J)*RHO(1,J)
      CON(2,J)=CON(2,J)+ARHO*U(2,J)
      AIM(2,J)=0.
      DO 404 I=2,L2
      IF(I==L2) THEN
      ARHO=ARX(J)*RHO(L1,J)
      CON(I,J)=CON(I,J)-ARHO*U(L1,J)
      AIP(I,J)=0.
      ELSE
      ARHO=ARX(J)*(FX(I+1)*RHO(I+1,J)+FXM(I+1)*RHO(I,J))
      FLOW=ARHO*UHAT(I+1,J)
      CON(I,J)=CON(I,J)-FLOW
      CON(I+1,J)=CON(I+1,J)+FLOW
      AIP(I,J)=ARHO*DU(I+1,J)
      AIM(I+1,J)=AIP(I,J)
      ENDIF
      IF(J==M2) THEN
      ARHO=RMN(M1)*XCV(I)*RHO(I,M1)
      CON(I,J)=CON(I,J)-ARHO*V(I,M1)
      AJP(I,J)=0.
      ELSE
      ARHO=RMN(J+1)*XCV(I)*(FY(J+1)*RHO(I,J+1)+FYM(J+1)*RHO(I,J))
      FLOW=ARHO*VHAT(I,J+1)
      CON(I,J)=CON(I,J)-FLOW
      CON(I,J+1)=CON(I,J+1)+FLOW
      AJP(I,J)=ARHO*DV(I,J+1)
      AJM(I,J+1)=AJP(I,J)
      ENDIF
      AP(I,J)=AIP(I,J)+AIM(I,J)+AJP(I,J)+AJM(I,J)
  404 ENDDO
  403 ENDDO
      DO 421 J=2,M2
      DO 422 I=2,L2
      AP(I,J)=AP(I,J)/RELAX(NP)
      CON(I,J)=CON(I,J)+(1.0-RELAX(NP))*AP(I,J)*P(I,J)
  422 ENDDO
  421 ENDDO
      COFP(IST:L2,JST:M2,2:5)=COF(IST:L2,JST:M2,2:5)
      NF=NP
      CALL SOLVE
C
COMPUTE U AND V
      NF=1
      IST=3
      JST=2
      COF(IST:L2,JST:M2,1:6)=COFU(IST:L2,JST:M2,1:6)
      DO 551 J=JST,M2
      DO 552 I=IST,L2
      CON(I,J)=CON(I,J)+DU(I,J)*AP(I,J)*(P(I-1,J)-P(I,J))
  552 ENDDO
  551 ENDDO
      CALL SOLVE
C
      NF=2
      IST=2
      JST=3
      COF(IST:L2,JST:M2,1:6)=COFV(IST:L2,JST:M2,1:6)
      DO 553 J=JST,M2
      DO 554 I=IST,L2
      CON(I,J)=CON(I,J)+DV(I,J)*AP(I,J)*(P(I,J-1)-P(I,J))
  554 ENDDO
  553 ENDDO
      CALL SOLVE
COEFFICIENTS FOR THE PRESSURE CORRECTION EQUATION
      NF=3
      CALL RESET
      IST=2
      JST=2
      COF(IST:L2,JST:M2,2:5)=COFP(IST:L2,JST:M2,2:5)
      CALL GAMSOR
      SMAX=0.
      SSUM=0.
      DO 510 J=2,M2
      DO 511 I=2,L2
      VOL=YCVR(J)*XCV(I)
      CON(I,J)=CON(I,J)*VOL
  511 ENDDO
  510 ENDDO
      DO 502 I=2,L2
      ARHO=R(1)*XCV(I)*RHO(I,1)
      CON(I,2)=CON(I,2)+ARHO*V(I,2)
  502 ENDDO
      DO 503 J=2,M2
      ARHO=ARX(J)*RHO(1,J)
      CON(2,J)=CON(2,J)+ARHO*U(2,J)
      DO 504 I=2,L2
      IF(I==L2) THEN
      ARHO=ARX(J)*RHO(L1,J)
      CON(I,J)=CON(I,J)-ARHO*U(L1,J)
      ELSE
      ARHO=ARX(J)*(FX(I+1)*RHO(I+1,J)+FXM(I+1)*RHO(I,J))
      FLOW=ARHO*U(I+1,J)
      CON(I,J)=CON(I,J)-FLOW
      CON(I+1,J)=CON(I+1,J)+FLOW
      ENDIF
      IF(J==M2) THEN
      ARHO=RMN(M1)*XCV(I)*RHO(I,M1)
      CON(I,J)=CON(I,J)-ARHO*V(I,M1)
      ELSE
      ARHO=RMN(J+1)*XCV(I)*(FY(J+1)*RHO(I,J+1)+FYM(J+1)*RHO(I,J))
      FLOW=ARHO*V(I,J+1)
      CON(I,J)=CON(I,J)-FLOW
      CON(I,J+1)=CON(I,J+1)+FLOW
      ENDIF
      AP(I,J)=AIP(I,J)+AIM(I,J)+AJP(I,J)+AJM(I,J)
      PC(I,J)=0.
      SMAX=MAX(SMAX,ABS(CON(I,J)))
      SSUM=SSUM+CON(I,J)
  504 ENDDO
  503 ENDDO
      CALL SOLVE
COME HERE TO CORRECT THE VELOCITIES
      DO 521 J=2,M2
      DO 522 I=2,L2
      IF(I/=2) U(I,J)=U(I,J)+DU(I,J)*(PC(I-1,J)-PC(I,J))
      IF(J/=2) V(I,J)=V(I,J)+DV(I,J)*(PC(I,J-1)-PC(I,J))
  522 ENDDO
  521 ENDDO
  500 ENDIF
COEFFICIENTS FOR OTHER EQUATIONS-----------------------------------
      IST=2
      JST=2
      DO 600 N=4,NFMAX
      NF=N
      IF(LSOLVE(NF)) THEN
      CALL GAMSOR
      IF(LSOLVE(4)) THEN
      DO I=1,L1
      DO J=1,M1
      RHO(I,J)=RHO(I,J)*CP(I,J)
      ENDDO
      ENDDO
      ENDIF
      REL=1.-RELAX(NF)
      DO 602 I=2,L2
      AREA=R(1)*XCV(I)
      FLOW=AREA*V(I,2)*RHO(I,1)
      DIFF=AREA*GAM(I,1)/YDIF(2)
      CALL DIFLOW
      AJM(I,2)=ACOF+MAX(0.,FLOW)
  602 ENDDO
      DO 603 J=2,M2
      FLOW=ARX(J)*U(2,J)*RHO(1,J)
      DIFF=ARX(J)*GAM(1,J)/(XDIF(2)*SX(J))
      CALL DIFLOW
      AIM(2,J)=ACOF+MAX(0.,FLOW)
      DO 604 I=2,L2
      IF(I==L2) THEN
      FLOW=ARX(J)*U(L1,J)*RHO(L1,J)
      DIFF=ARX(J)*GAM(L1,J)/(XDIF(L1)*SX(J))
      ELSE
      FLOW=ARX(J)*U(I+1,J)*(FX(I+1)*RHO(I+1,J)+FXM(I+1)*RHO(I,J))
      DIFF=ARX(J)*2.*GAM(I,J)*GAM(I+1,J)/((XCV(I)*GAM(I+1,J)+
     1 XCV(I+1)*GAM(I,J)+1.0E-30)*SX(J))
      ENDIF
      CALL DIFLOW
      AIM(I+1,J)=ACOF+MAX(0.,FLOW)
      AIP(I,J)=AIM(I+1,J)-FLOW
      AREA=RMN(J+1)*XCV(I)
      IF(J==M2) THEN
      FLOW=AREA*V(I,M1)*RHO(I,M1)
      DIFF=AREA*GAM(I,M1)/YDIF(M1)
      ELSE
      FLOW=AREA*V(I,J+1)*(FY(J+1)*RHO(I,J+1)+FYM(J+1)*RHO(I,J))
      DIFF=AREA*2.*GAM(I,J)*GAM(I,J+1)/(YCV(J)*GAM(I,J+1)+
     1 YCV(J+1)*GAM(I,J)+1.0E-30)
      ENDIF
      CALL DIFLOW
      AJM(I,J+1)=ACOF+MAX(0.,FLOW)
      AJP(I,J)=AJM(I,J+1)-FLOW
      VOL=YCVR(J)*XCV(I)
      APT=RHO(I,J)/DT
      AP(I,J)=AP(I,J)-APT
      CON(I,J)=CON(I,J)+APT*F(I,J,NF)
      AP(I,J)=(-AP(I,J)*VOL+AIP(I,J)+AIM(I,J)+AJP(I,J)+AJM(I,J))
     1/RELAX(NF)
      CON(I,J)=CON(I,J)*VOL+REL*AP(I,J)*F(I,J,NF)
  604 ENDDO
  603 ENDDO
      CALL SOLVE
      IF(LSOLVE(4)) THEN
      DO I=1,L1
      DO J=1,M1
      RHO(I,J)=RHO(I,J)/CP(I,J)
      ENDDO
      ENDDO
      ENDIF
      ENDIF
  600 ENDDO
      TIME=TIME+DT
      ITER=ITER+1
      IF(ITER>=LAST) LSTOP=.TRUE.
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SUPPLY
C******************************************************************
      USE START_L
      IMPLICIT NONE
      REAL*8 DX,DY,RHOM,PREF,FSHOW
      INTEGER*4 I,J,N,JJ,IEND,JEND,IBEG,JBEG,IFST,JFST,JFL
C******************************************************************
   10 FORMAT(1X,26(1H*),3X,A10,3X,26(1H*))
   20 FORMAT(1X,4H I =,I6,6I10)
   30 FORMAT(1X,' J')
   40 FORMAT(1X,I3,2X,1P7E10.3)
   50 FORMAT(1X,1H )
   51 FORMAT(2X,'I =',2X,7(I4,5X))
   52 FORMAT(2X,'X =',1P7E10.3)
   53 FORMAT(1X,' TH =',1P7E10.3)
   54 FORMAT(2X,'J =',2X,7(I4,5X))
   55 FORMAT(2X,'Y =',1P7E10.3)
C******************************************************************
      ENTRY UGRID
      XU(2)=0.
      DX=XL/FLOAT(L1-2)
      DO 1 I=3,L1
      XU(I)=XU(I-1)+DX
    1 ENDDO
      YV(2)=0.
      DY=YL/FLOAT(M1-2)
      DO 2 J=3,M1
      YV(J)=YV(J-1)+DY
    2 ENDDO
      RETURN
C************************************************
      ENTRY PRINT
      IF(LPRINT(3)) THEN
CALCULATE THE STREAM FUNCTION---------------------------------------
      F(2,2,3)=0.
      DO 81 I=2,L1
      IF(I/=2) F(I,2,3)=F(I-1,2,3)-RHO(I-1,1)*V(I-1,2)
     1*R(1)*XCV(I-1)
      DO 82 J=3,M1
      RHOM=FX(I)*RHO(I,J-1)+FXM(I)*RHO(I-1,J-1)
      F(I,J,3)=F(I,J-1,3)+RHOM*U(I,J-1)*ARX(J-1)
   82 ENDDO
   81 ENDDO
      ENDIF
C
      IF(LPRINT(NP)) THEN
C
CONSTRUCT BOUNDARY PRESSURES BY EXTRAPOLATION
      DO 91 J=2,M2
      P(1,J)=(P(2,J)*XCVS(3)-P(3,J)*XDIF(2))/XDIF(3)
      P(L1,J)=(P(L2,J)*XCVS(L2)-P(L3,J)*XDIF(L1))/XDIF(L2)
   91 ENDDO
      DO 92 I=2,L2
      P(I,1)=(P(I,2)*YCVS(3)-P(I,3)*YDIF(2))/YDIF(3)
      P(I,M1)=(P(I,M2)*YCVS(M2)-P(I,M3)*YDIF(M1))/YDIF(M2)
   92 ENDDO
      P(1,1)=P(2,1)+P(1,2)-P(2,2)
      P(L1,1)=P(L2,1)+P(L1,2)-P(L2,2)
      P(1,M1)=P(2,M1)+P(1,M2)-P(2,M2)
      P(L1,M1)=P(L2,M1)+P(L1,M2)-P(L2,M2)
      PREF=P(IPREF,JPREF)
      DO 93 J=1,M1
      DO 94 I=1,L1
      P(I,J)=P(I,J)-PREF
   94 ENDDO
   93 ENDDO
      ENDIF
C
      PRINT 50
      WRITE(8,50)
      IEND=0
      DO WHILE(IEND/=L1) 
      IBEG=IEND+1
      IEND=IEND+7
      IEND=MIN0(IEND,L1)
      PRINT 50
      WRITE(8,50)
      PRINT 51,(I,I=IBEG,IEND)
      WRITE(8,51) (I,I=IBEG,IEND)
      IF(MODE/=3) THEN
      PRINT 52,(X(I),I=IBEG,IEND)
      WRITE(8,52) (X(I),I=IBEG,IEND)
      ELSE
      PRINT 53,(X(I),I=IBEG,IEND)
      WRITE(8,53) (X(I),I=IBEG,IEND)
      ENDIF
      ENDDO
      IF(IEND==L1) THEN
      JEND=0
      PRINT 50
      WRITE(8,50)
      DO WHILE(JEND/=M1)
      JBEG=JEND+1
      JEND=JEND+7
      JEND=MIN0(JEND,M1)
      PRINT 50
      WRITE(8,50)
      PRINT 54,(J,J=JBEG,JEND)
      WRITE(8,54) (J,J=JBEG,JEND)
      PRINT 55,(Y(J),J=JBEG,JEND)
      WRITE(8,55) (Y(J),J=JBEG,JEND)
      ENDDO  
      ENDIF
C
      DO 999 N=1,NCP
      NF=N
      IF(LPRINT(NF)) THEN
      PRINT 50
      WRITE(8,50)
      PRINT 10,TITLE(NF)
      WRITE(8,10) TITLE(NF)
      IFST=1
      JFST=1
      IF(NF==1.OR.NF==3) IFST=2
      IF(NF==2.OR.NF==3) JFST=2
      IBEG=IFST-7
      DO WHILE(IEND<L1.OR.IBEG==-5.OR.IBEG==-6)
      IBEG=IBEG+7
      IEND=IBEG+6
      IEND=MIN0(IEND,L1)
      PRINT 50
      WRITE(8,50)
      PRINT 20,(I,I=IBEG,IEND)
      WRITE(8,20) (I,I=IBEG,IEND)
      PRINT 30
      WRITE(8,30)
      JFL=JFST+M1
      DO 115 JJ=JFST,M1
      J=JFL-JJ
      PRINT 40,J,(F(I,J,NF),I=IBEG,IEND)
      WRITE(8,40) J,(F(I,J,NF),I=IBEG,IEND)
  115 ENDDO
      ENDDO
      ENDIF
  999 ENDDO
      OPEN(9,FILE="RESULT.DAT")
      WRITE(9,'("VARIABLES=X,Y",$)')
	DO NF=1,NCP
	IF(LPRINT(NF)) WRITE(9,'(",",A7,$)') TITLE(NF)
	ENDDO
      WRITE(9,'(/,"ZONE I=",I4,",J=",I4,",T=T",$)') L1,M1
      DO J=1,M1
      DO I=1,L1
	WRITE(9,'(/,E11.3,E11.3,$)') X(I),Y(J)
      DO NF=1,NCP
	IF(LPRINT(NF)) THEN
      FSHOW=F(I,J,NF)
      IF(NF==1) THEN
        IF(I==1) FSHOW=U(2,J)
        IF(I>=2.AND.I<=L2)  FSHOW=(U(I,J)+U(I+1,J))/2
        IF(I==L1) FSHOW=U(L1,J)
      ENDIF
      IF(NF==2) THEN
        IF(J==1) FSHOW=V(I,2)
        IF(J>=2.AND.J<=M2) FSHOW=(V(I,J)+V(I,J+1))/2
        IF(J==M1) FSHOW=V(I,M1)
      ENDIF
      WRITE(9,'(E11.3,$)') FSHOW
      ENDIF
      ENDDO
	ENDDO
      ENDDO
      CLOSE(9)
      RETURN
      END