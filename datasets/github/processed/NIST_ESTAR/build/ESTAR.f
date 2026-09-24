      PROGRAM ESTAR
C
C     8 September 1992.
C
C     Version 1.2
C
C     Written by Martin J. Berger,   
C     National Institute of Standards and Technology,
C     Gaithersburg, MD 20899.
C
C     Minor changes made in October of 1998 for web version 1.0
C     These included: a line too long, and file open routines that
C     were not supported by our compiler (Johnathan S. Coursey)
C
C     Change made June of 1999 for web version 1.2
C     This included: rounding data to 4 significant figures
C     (it is only accurate to 3 anyway) and formatting data into
C     exponential format (Johnathan S. Coursey)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ATB(100),MZ(14),WT(14),G(14),AT(14),NC(26),BD(26),
     1 T(1001),TL(1001),CLOSS(1001),RLOSS(1001),RLOSSL(1001),
     2 TLOSS(1001),TLOSSL(1001),AST(1001),BST(1001),CST(1001),DST(1001),
     4 ART(1001),BRT(1001),CRT(1001),DRT(1001),DLT(1001),ER(113),
     5 ERL(113),RLOS(113),RLOST(113),RLOSTL(113),ARL(113),BRL(113),
     6 CRL(113),DRL(113),ALF(1000),EPS(1000),EN(1000),F(1000),Q(1200),
     7 YQ(1200),YQL(1200),D(1200),ADEL(1200),BDEL(1200),CDEL(1200),
     8 DDEL(1200),GRAND(81),GRAND1(81),RG(113),RAD(113)
      CHARACTER ENGIN*30,MAT*72,FF*1,OUTPUT*30,
     1 DATFIL*30,HEAD1*13,BL*61,HEAD2*120,HEAD3*120
      DATA ATB/
     1 1.00794D0,4.002602D0,6.941D0,9.012182D0,10.811D0,12.011D0,
     2 14.00674D0,15.9994D0,18.9984032D0,20.1797D0,22.989768D0,
     3 24.3050D0,26.981539D0,28.0855D0,30.973762D0,32.066D0,35.4527D0,
     4 39.948D0,39.0983D0,40.078D0,44.955910D0,47.88D0,50.9415D0,
     5 51.9961D0,54.93805D0,55.847D0,58.93320D0,58.69D0,63.546D0,
     6 65.39D0,69.723D0,72.61D0,74.92159D0,78.96D0,79.904D0,83.80D0,
     7 85.4678D0,87.62D0,88.90585D0,91.224D0,92.90638D0,95.94D0,
     8 97.9072D0,101.07D0,102.9055D0,106.42D0,107.8682D0,112.411D0,
     9 114.82D0,118.710D0,121.75D0,127.60D0,126.90447D0,131.29D0,
     1 132.90543D0,137.327D0,138.9055D0,140.115D0,140.90765D0,144.24D0,
     2 144.9127D0,150.36D0,151.965D0,157.25D0,158.92534D0,162.50D0,
     3 164.93032D0,167.26D0,168.93421D0,173.04D0,174.967D0,178.49D0,
     4 180.9479D0,183.85D0,186.207D0,190.2D0,192.22D0,195.08D0,
     5 196.96654D0,200.59D0,204.3833D0,207.2D0,208.98037D0,208.9824D0,
     6 209.9871D0,222.0176D0,223.0197D0,226.0254D0,227.0278D0,
     7 232.0381D0,231.03588D0,238.0289D0,237.0482D0,239.0522D0,
     8 243.0614D0,247.0703D0,247.0703D0,251.0796D0,252.083D0,
     9 257.0951D0/
      DATA LKMAX/113/
      DATA ER/1.00E-03,1.25E-03,1.50E-03,1.75E-03,2.00E-03,2.50E-03,
     1        3.00E-03,3.50E-03,4.00E-03,4.50E-03,5.00E-03,5.50E-03,
     2        6.00E-03,7.00E-03,8.00E-03,9.00E-03,1.00E-02,1.25E-02,
     3        1.50E-02,1.75E-02,2.00E-02,2.50E-02,3.00E-02,3.50E-02,
     4        4.00E-02,4.50E-02,5.00E-02,5.50E-02,6.00E-02,7.00E-02,
     5        8.00E-02,9.00E-02,1.00E-01,1.25E-01,1.50E-01,1.75E-01,
     6        2.00E-01,2.50E-01,3.00E-01,3.50E-01,4.00E-01,4.50E-01,
     7        5.00E-01,5.50E-01,6.00E-01,7.00E-01,8.00E-01,9.00E-01,
     8        1.00E+00,1.25E+00,1.50E+00,1.75E+00,2.00E+00,2.50E+00,
     9        3.00E+00,3.50E+00,4.00E+00,4.50E+00,5.00E+00,5.50E+00,
     1        6.00E+00,7.00E+00,8.00E+00,9.00E+00,1.00E+01,1.25E+01,
     2        1.50E+01,1.75E+01,2.00E+01,2.50E+01,3.00E+01,3.50E+01,
     3        4.00E+01,4.50E+01,5.00E+01,5.50E+01,6.00E+01,7.00E+01,
     4        8.00E+01,9.00E+01,1.00E+02,1.25E+02,1.50E+02,1.75E+02,
     5        2.00E+02,2.50E+02,3.00E+02,3.50E+02,4.00E+02,4.50E+02,
     6        5.00E+02,5.50E+02,6.00E+02,7.00E+02,8.00E+02,9.00E+02,
     7        1.00E+03,1.25E+03,1.50E+03,1.75E+03,2.00E+03,2.50E+03,
     8        3.00E+03,3.50E+03,4.00E+03,4.50E+03,5.00E+03,5.50E+03,
     9        6.00E+03,7.00E+03,8.00E+03,9.00E+03,1.00E+04/
      DATA QBEG/1.0D-04/,NUMQ/50/,LMAX/1101/,MDAUX/0/,ISTORE/2/
      DATA COFF/0.307072D0/,RMASS/0.510999906D0/,MGRD/21/
      DATA LP/73/
    1 FORMAT(1H )
    5 FORMAT(A)
      OPEN (3,FILE='UEDAT',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=1224)
      OPEN (4,FILE='UCOMP',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=268)
      BL='                                                             ' 
      FF=CHAR(12)
      PRINT *,' Output options: '
      PRINT *,'    1) Stopping powers, ranges and radiation yields,'
      PRINT *,'       for standard energy grid'
      PRINT *,'    2) Stopping powers only,'
      PRINT *,'       for user-selected energy grid'
      PRINT *,' Choose 1 or 2: '
      READ *, MOUT
      PRINT *,' Options for entering properties of stopping material:'
      PRINT *,'    1) Use input from composition file UCOMP'
      PRINT *,'    2) Enter composition data from keyboard, as prompted'
      PRINT *,' Choose 1 or 2: '
      READ *, INMAT
      GO TO (15,10), INMAT
   10 CALL CPREP (MAT,KMAT,POT,RHO,MMAX,MZ,WT,ISTORE,DATFIL)
      POTL=LOG(POT*1.0D-06)
      GO TO 35
   15 PRINT *,' Enter ID number of material: '
      READ *, JPIC
      IF(JPIC-278)20,20,16
   16 IF(JPIC-906)18,17,18
   17 JPIC=279 
      GO TO 20
   18 PRINT 19, JPIC
   19 FORMAT(I6,' is not an allowed ID number.')
      STOP 1 
   20 READ (4,REC=JPIC) MAT,MMAX,ZAG,POT,RHO,MZ,WT
      PRINT 25,POT
   25 FORMAT(' I-value from UCOMP file is = ',F6.1,' eV.')
      PRINT *, ' Is this value acceptable (1=yes,2=no): '
      READ *, IPOT
      IF(IPOT.EQ.1) GO TO 30 
      PRINT *,' Enter desired I-value (eV): '
      READ *, POT
   30 POTL=LOG(POT*1.0E-06)
   35 IF(MOUT.EQ.1) GO TO 80
      PRINT *,' Options for entering energy list:'
      PRINT *,'    1) Use default file ENG.ELE'
      PRINT *,'    2) Use prepared file'
      PRINT *,'    3) Entry from keyboard'
      PRINT *,' Choose 1, 2, or 3: '
      READ *, INEN
      GO TO (40,50,70), INEN
   40 ENGIN='ENG.ELE'
      GO TO 55  
   50 PRINT *,' Enter name of energy-list file: '
      READ 5, ENGIN
   55 OPEN (7,FILE=ENGIN)
      READ (7,*) IMAX
      READ (7,*) (T(I),I=1,IMAX)
      CLOSE (7)
   60 TMIN=1.0E12
      DO 61 I=1,IMAX
      IF(T(I).LT.TMIN) TMIN=T(I)
   61 CONTINUE
      IF(TMIN-9.9999D-04)62,64,64
   62 PRINT 63
   63 FORMAT(' At least one of the specified energies is below 0.001 MeV
     1, out or range.')
      STOP 
   64 IF(TMIN-0.01D0)65,80,80
   65 PRINT 66
   66 FORMAT('  Warning: at energies below 0.01 MeV, accuracy of')
      PRINT 67
   67 FORMAT('           collision stopping may be poor.')      
      GO TO 80 
   70 PRINT *,' Specify the number of energies in the list: '
      READ *, IMAX
      PRINT *,' Enter all energies (in MeV): '
      READ *,(T(I),I=1,IMAX)
      GO TO 60
   80 PRINT *,' Enter name of output file: '
      READ 5, OUTPUT
      OPEN (UNIT=8,FILE=OUTPUT)
      QFAC=10.0D0**(1.0D0/DBLE(NUMQ))
      Q(1)=QBEG
      DO 90 L=2,LMAX
   90 Q(L)=Q(L-1)*QFAC
      DO 100 LK=1,LKMAX
  100 ERL(LK)=LOG(ER(LK))
      IF(MOUT.EQ.2) GO TO 106
      IMAX=LKMAX
      DO 105 I=1,IMAX
      T(I)=ER(I)
  105 TL(I)=ERL(I)
  106 GTOT=0.0
      DO 110 M=1,MMAX
      JZ=MZ(M)
      AT(M)=ATB(JZ) 
      Z=DBLE(MZ(M))
      A=AT(M) 
      G(M)=WT(M)*Z/A 
  110 GTOT=GTOT+G(M)
      ZAV=GTOT
      HOM=28.81593D0*SQRT(RHO*ZAV)
      PHIL=2.0D0*LOG(POT/HOM)
      CBAR=PHIL+1.0D0 
      DO 120 M=1,MMAX
  120 G(M)=G(M)/GTOT
      NBAS=0
      DO 130 LK=1,LKMAX
  130 RLOST(LK)=0.0 
      DO 200 M=1,MMAX
      IZ=MZ(M)
      READ (3,REC=IZ) NMAX,LKIN,NC,BD,RLOS 
      IF(LKIN.NE.LKMAX) STOP 2
      DO 140 LK=1,LKMAX
  140 RLOST(LK)=RLOST(LK)+WT(M)*RLOS(LK)
      IF(NC(NMAX))150,170,170 
  150 NC(NMAX)=-NC(NMAX)
      IF(MMAX-1)160,160,170
  160 BD(NMAX)=0.0
  170 NSUM=0
      DO 180 N=1,NMAX
  180 NSUM=NSUM+NC(N)
      SUM=DBLE(NSUM)
      DO 190 N=1,NMAX
      NN=N+NBAS
      F(NN)=NC(N)*G(M)/SUM
  190 EN(NN)=BD(N)
  200 NBAS=NBAS+NMAX
      DO 210 LK=1,LKMAX
  210 RLOSTL(LK)=LOG(RLOST(LK))
      CALL SCOF(ERL,RLOSTL,LKMAX,ARL,BRL,CRL,DRL)
      NMAX=NBAS
      DO 390 N=1,NMAX
  390 ALF(N)=2.0/3.0
      IF(EN(NMAX))400,400,410
  400 ALF(NMAX)=1.0 
  410 DO 420 N=1,NMAX
  420 EPS(N)=(EN(N)/HOM)**2
      ROOT=1.0
  430 FUN=-PHIL
      DER=0.0
      DO 440 N=1,NMAX
      TRM=ROOT*EPS(N)+ALF(N)*F(N)
      FUN=FUN+F(N)*LOG(TRM)
  440 DER=DER+F(N)*EPS(N)/TRM 
      DROOT=FUN/DER 
      ROOT=ROOT-DROOT
      IF(ABS(DROOT)-0.00001)450,450,430 
  450 FACTOR=SQRT(ROOT)
      DO 460 N=1,NMAX
  460 EPS(N)=ROOT*EPS(N)
      IF(EN(NMAX))470,470,480
  470 YCUT=0.0
      GO TO 500
  480 SUM=0.0
      DO 490 N=1,NMAX
  490 SUM=SUM+F(N)/EPS(N)
      YCUT=1.0/SUM
  500 DO 530 L=1,LMAX
      SUM=0.0
      DO 510 N=1,NMAX
  510 SUM=SUM+F(N)/(EPS(N)+Q(L))
      YQ(L)=1.0/SUM 
      YQL(L)=LOG(YQ(L))
      SUM=0.0
      DO 520 N=1,NMAX
      ARG=1.0+Q(L)/(EPS(N)+ALF(N)*F(N)) 
  520 SUM=SUM+F(N)*LOG(ARG)
      D(L)=SUM-Q(L)/(YQ(L)+1.0)
  530 CONTINUE
      TCUT=RMASS*(SQRT(YCUT+1.0)-1.0)
      IF(MDAUX.EQ.0) GO TO 560
      OPEN (UNIT=10,FILE='DAUX')
      WRITE (10,540)
  540 FORMAT(9X,'ZAV',9X,'RHO',9X,'HOM',9X,'POT',6X,'FACTOR',
     1 8X,'CBAR',8X,'TCUT',2X,'MMAX',2X,'NMAX')
      WRITE (10,545) ZAV,RHO,HOM,POT,FACTOR,CBAR,TCUT,MMAX,NMAX
  545 FORMAT(7F12.6,2I6)
      WRITE (10,1)
      WRITE (10,550) (YQ(L),L=1,LMAX,10)
      WRITE (10,550) (D(L),L=1,LMAX,10)
  550 FORMAT(1P4E17.9)
      CLOSE (10)
  560 CALL SCOF(YQL,D,LMAX,ADEL,BDEL,CDEL,DDEL)
      DO 595 I=1,IMAX
      E=T(I)
      TAU=T(I)/RMASS
      Y=TAU*(TAU+2.0)
      BETQ=Y/((TAU+1.0)**2)
      DELTA=0.0D0
      IF(Y.LT.YQ(1)) GO TO 590
      IF(Y.LE.YCUT) GO TO 590
      IF(Y-YQ(LMAX))585,585,580
  580 PRINT *,' YQ(LMAX) out of range'
      STOP 3
  585 YL=LOG(Y)
      CALL BSPOL(YL,YQL,ADEL,BDEL,CDEL,DDEL,LMAX,DELTA)
  590 DENCOR=0.5D0*DELTA
      DLT(I)=DELTA
      SPART=LOG(T(I))-POTL+0.5*LOG(1.0+0.5*TAU)-DENCOR
      TERM=(1.0-BETQ)*(1.0+(TAU**2)/8.0-(2.0*TAU+1.0)*LOG(2.0))
      STNUM=SPART+0.5*TERM
      CLOSS(I)=COFF*ZAV*STNUM/BETQ
      EL=LOG(E)
      CALL BSPOL(EL,ERL,ARL,BRL,CRL,DRL,LKMAX,RES)
      RLOSS(I)=EXP(RES)
      TLOSS(I)=CLOSS(I)+RLOSS(I)
  595 CONTINUE
      IF(MOUT.EQ.2) GO TO 700
      DO 600 I=1,IMAX
      RLOSSL(I)=LOG(RLOSS(I))
  600 TLOSSL(I)=LOG(TLOSS(I))
      CALL SCOF(TL,TLOSSL,IMAX,AST,BST,CST,DST)
      CALL SCOF(TL,RLOSSL,IMAX,ART,BRT,CRT,DRT)
      RG(1)=0.5D0*T(1)/TLOSS(1)
      RAD(1)=0.5D0*T(1)*RLOSS(1)/TLOSS(1)
      DO 690 I=2,IMAX
      ETMAX=T(I)
      ETMIN=T(I-1)
      EDIFF=(ETMAX-ETMIN)/DBLE(MGRD-1)
      DET=EDIFF/3.0D0
      DO 685 M=1,MGRD
      ETL=LOG(ETMAX-EDIFF*DBLE(M-1))
      CALL BSPOL(ETL,TL,AST,BST,CST,DST,IMAX,RES)
      GRAND(M)=1.0D0*EXP(-RES)
      CALL BSPOL(ETL,TL,ART,BRT,CRT,DRT,IMAX,RESR)
  685 GRAND1(M)=EXP(RESR)*GRAND(M)
      CALL GRAL(DET,GRAND,MGRD,STEP)
      CALL GRAL(DET,GRAND1,MGRD,DRAD)
      RG(I)=RG(I-1)+STEP
      RAD(I)=RAD(I-1)+DRAD
  690 CONTINUE
      DO 695 I=1,IMAX
  695 RAD(I)=RAD(I)/T(I)
      GO TO 790
  700 WRITE (8,705) MAT
  705 FORMAT(6X,'Electrons in ',A)
      WRITE (8,1)
      WRITE (8,706) 
  706 FORMAT(9X,'Z/A',2X,'Density, g/cm3',4X,'I, eV')
      WRITE (8,707) ZAG,RHO,POT
  707 FORMAT(0PF12.6,1PE16.4,0PF9.1)
      WRITE (8,1)
      WRITE (8,710)
  710 FORMAT('     T = kinetic energy, MeV')
      WRITE (8,720)
  720 FORMAT(' CLOSS = collision stopping power, MeV cm2/g')
      WRITE (8,730) 
  730 FORMAT(' RLOSS = radiative stopping power, MeV cm2/g')
      WRITE (8,740)
  740 FORMAT(' TLOSS = total stopping power, MeV cm2/g')
      WRITE (8,750)
  750 FORMAT(' DELTA = density effect parameter delta') 
      WRITE (8,1)
      WRITE (8,760)
  760 FORMAT(11X,'T',7X,'CLOSS',7X,'RLOSS',7X,'TLOSS',7X,'DELTA')
      DO 780 I=1,IMAX
      WRITE (8,770) T(I),CLOSS(I),RLOSS(I),TLOSS(I),DLT(I)
C changed from FORMAT(1PE12.4,1P3E12.4,0PF12.5)
C we wanted all exponential and only 4 sig figs
C Jack Coursey 6/1999
  770 FORMAT(1P5E12.3)
  780 CONTINUE
      GO TO 875
  790 HEAD1='Electrons in '
      HEAD2=HEAD1//MAT
      DO 791 L=120,1,-1
      IF(HEAD2(L:L).NE.' ') GO TO 792
  791 CONTINUE
  792 LENGTH=L
      LSHIFT=6+(LP-LENGTH)/2
      HEAD3=BL(1:LSHIFT)//HEAD2
      WRITE (8,5) HEAD3
      WRITE (8,1)
      WRITE (8,793) 
  793 FORMAT(30X,'Density, g/cm3',5x,'I, eV')
      WRITE (8,794) RHO,POT
  794 FORMAT(30X,1PE12.4,0PF12.1)
      WRITE (8,1)
      WRITE (8,1)
      WRITE (8,1)
      WRITE (8,800)
  800 FORMAT ( 2X,'    KINETIC         STOPPING POWER             CSDA
c 800 FORMAT ( 2X,'    ENERGY          STOPPING POWER             CSDA
     1   RADIATION   DENSITY')
      WRITE (8,810)
  810 FORMAT ( 2X,'    ENERGY  COLLISION  RADIATIVE    TOTAL      RANGE
c 810 FORMAT ( 2X,'            COLLISION  RADIATIVE    TOTAL      RANGE
     1     YIELD      EFFECT')
      WRITE (8,820)
  820 FORMAT ( 2X,'
     1                DELTA')
      WRITE (8,830)
  830 FORMAT ( 2X,'     MeV    MeV cm2/g  MeV cm2/g  MeV cm2/g    g/cm2'
     1)
      DO 855 I=17,56
      IF (REAL((I-1)/8).NE.REAL(I-1)/8.0) GO TO 840
      WRITE (8,1)
  840 WRITE (8,850) T(I),CLOSS(I),RLOSS(I),TLOSS(I),RG(I),RAD(I),DLT(I)
C changed from FORMAT (2X,0PF10.4,1P611.3)
C we wanted all exponential: Jack Coursey 6/1999
  850 FORMAT (2X,1P7E11.3)
  855 CONTINUE
      WRITE (8,5) FF
      WRITE (8,5) HEAD3
      WRITE (8,1)
      WRITE (8,793) 
      WRITE (8,794) RHO,POT
      WRITE (8,1)
      WRITE (8,1)
      WRITE (8,1)
      WRITE (8,800)
      WRITE (8,810)
      WRITE (8,820)
      WRITE (8,830)
      DO 870 I=57,97
      IF (REAL((I-1)/8).NE.REAL(I-1)/8.0) GO TO 860
      WRITE (8,1)
  860 WRITE (8,850) T(I),CLOSS(I),RLOSS(I),TLOSS(I),RG(I),RAD(I),DLT(I)
  870 CONTINUE
  875 CLOSE (8)
      GO TO (880,900), ISTORE
  880 PRINT 890, DATFIL
  890 FORMAT('  Composition file ',A)
      PRINT *,' has been generated.'
  900 STOP  
      END
      SUBROUTINE CPREP (MAT,KMAT,POT,RHO,MMAX,JZ,WT,ISTOR,OUTFIL)
C     30 Jul 89. Prepares input data for compound or mixture. 
C                Derived in part from SPEC in XCOM, 2 Apr 87.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION JZ(14),WT(140),JZ1(14),WT1(14),LH(100),WATE(100),
     1 FRAC(100),JZIP(100),ZAV1(14),POT1(14)
      CHARACTER FORMLA*72,FRM(100)*30,MAT*72,MAT1*72,
     1 OUTFIL*30,IDNO*3 
      DATA RHOCUT/0.1/
    1 FORMAT(1H )
    5 FORMAT(A)
      PRINT *,' Enter name of material: '
      READ 5,  MAT
      PRINT *,' Enter density of material (g/cm3): '
      READ *, RHO
      NPHAS=2
      IF(RHO.LT.RHOCUT) NPHAS=1 
      PRINT *,' Should the composition data be stored in a file?'
      PRINT *,'     (1 = yes; 2 = no): '
      READ *, ISTOR
      IF(ISTOR.EQ.2) GO TO 8
      PRINT *,' Enter name of this data file: '
      READ 5, OUTFIL
    8 PRINT *,' Options for type of material:'
      PRINT *,'    1) Element'
      PRINT *,'    2) Compound'
      PRINT *,'    3) Mixture of elements and/or compounds'
      PRINT *,' Choose 1, 2 or 3: '
      READ *, KMAT
      GO TO (10,12,20), KMAT
   10 PRINT *,' Enter chemical symbol for element: '
      GO TO 15
   12 PRINT *,' Enter chemical formula for compound: '
   15 READ 5,FORMLA
      CALL FORMEL(FORMLA,KMAT,NPHAS,MMAX,JZ,WT,ZAV,POT)
      GO TO 280 
   20 PRINT *,' How many components in mixture? Enter number: '  
      READ *,NCOMP
      DO 90 N=1,NCOMP
      PRINT 30, N
   30 FORMAT(' Choice for component',I3,': ')
      PRINT *,'    1) Use input from composition file UCOMP'
      PRINT *,'    2) Enter chemical formula'
      PRINT *,' Choose 1 or 2: '
      READ *,JINP
      JZIP(N)=0
      GO TO (40,50) JINP
   40 PRINT 45, N
   45 FORMAT(' Enter ID number for component',I3,': ')
      READ *, JZIP(N)
      CALL INDEX(JZIP(N),IDNO)
      FRM(N)=IDNO
      GO TO 70
   50 PRINT 60, N
   60 FORMAT(' Enter chemical symbol or formula for component',
     1 I3,': ')
      READ 5, FRM(N)
   70 PRINT 80, N
   80 FORMAT(' Enter fraction by weight for component',I3,': ')
      READ *, FRAC(N)
   90 CONTINUE
      PRINT 1
      SUMF=0.0
      DO 100 N=1,NCOMP
  100 SUMF=SUMF+FRAC(N)
      PRINT 110
  110 FORMAT('   Component    Fraction')
      PRINT 120
  120 FORMAT('               by Weight')
      DO 140 N=1,NCOMP
      PRINT 130,N,FRAC(N),FRM(N)
  130 FORMAT(I12,F12.6,3X,A)
  140 CONTINUE
      PRINT 150, SUMF
  150 FORMAT(6X,'Sum = ',F12.6)
      PRINT 1
      PRINT *,' Options for accepting or rejecting composition data:'
      PRINT *,'     1. Accept, but let program normalize fractions'
      PRINT *,'        by weight so that their sum is unity'
      PRINT *,'     2. Reject, and enter different set of fractions'
      PRINT *,' Choose 1 or 2: '
      READ *, MSUMGO
      GO TO (160,20), MSUMGO
  160 DO 170 N=1,NCOMP
  170 FRAC(N)=FRAC(N)/SUMF
      DO 180 L=1,100
  180 LH(L)=0
      DO 250 N=1,NCOMP
      IF(JZIP(N))190,190,200
  190 CALL FORMEL (FRM(N),KMAT,NPHAS,MAX,JZ1,WT1,ZAV1(N),POT1(N))
      GO TO 220
  200 JPIC=JZIP(N)
      IF(JPIC-278)210,210,201
  201 IF(JPIC-906)203,202,203
  202 JPIC=279 
      GO TO 210
  203 PRINT 204, JPIC
  204 FORMAT(I6,' is not an allowed ID number.')
      STOP 4 
  210 READ (4,REC=JPIC) MAT1,MAX,ZAVIN1,POTIN1,RHO1,JZ1,WT1
      ZAV1(N)=ZAVIN1
      POT1(N)=POTIN1
  220 DO 250 M=1,MAX
      IN=JZ1(M)
      IF(LH(IN))230,230,240
  230 LH(IN)=1
      WATE(IN)=FRAC(N)*WT1(M)
      GO TO 250
  240 WATE(IN)=WATE(IN)+FRAC(N)*WT1(M)
  250 CONTINUE
      LL=0
      DO 260 L=1,100
      IF(LH(L))260,260,255
  255 LL=LL+1
      JZ(LL)=L
      WT(LL)=WATE(L)
  260 CONTINUE
      MMAX=LL
      ZAV=0.0
      POTL=0.0
      DO 270 N=1,NCOMP
      ZAV=ZAV+FRAC(N)*ZAV1(N)
  270 POTL=POTL+FRAC(N)*ZAV1(N)*LOG(POT1(N))
      POT=EXP(POTL/ZAV)
  280 PRINT 290, POT
  290 FORMAT('  I-value computed by program is = ',F6.1,' eV.')
      PRINT *, ' Is this value acceptable (1=yes,2=no): '
      READ *, IPOT
      IF(IPOT.EQ.1) GO TO 295 
      PRINT *,' Enter desired I-value (eV): '
      READ *, POT
  295 GO TO (300,340), ISTOR
  300 OPEN (UNIT=8,FILE=OUTFIL)
      WRITE (8,310) MAT
  310 FORMAT(1X,A)
      WRITE (8,320) MMAX,ZAV,POT,RHO
  320 FORMAT(I6,0P2F12.6,1PE12.5)
      WRITE (8,330) (JZ(M),WT(M),M=1,MMAX)
  330 FORMAT(6(I3,F9.6))
      CLOSE (8)
  340 RETURN
      END
      SUBROUTINE FORMEL (W,KMAT,NPHAS,MMAX,JZ,WT,ZAV,POT)
C     23 Jul 89. Interprets chemical formula.  
C                Derived from FORM in XCOM, 24 Mar 87
C
C                W: chemical formula
C                KMAT: Phase index
C                MMAX: number of atomic constituents
C                JZ(M),M=1,MMAX: atomic numbers of constituents
C                WT(M),M-1,MMAX: fractions by weight
C                ZAV:  <Z/A>
C                POT: mean excitation energy of compound (eV)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION MASH1(26),MASH2(418),IC(72),K(72),NZ(100),MS(100),
     1 ATWTD(100),POTGAS(9),POTCON(9),POTH(100),JZ(100),WT(100)
      CHARACTER W*72
      DATA MASH1/0,5,6,0,0,9,0,1,53,0,19,0,0,7,8,15,0,0,16,
     1 0,92,23,74,0,39,0/
      DATA (MASH2(I),I=1,111)/70,36*0,54,3*0,73,8*0,65,8*0,43,51,88,
     1 7*0,21,37,6*0,52,3*0,91,5*0,34,2*0,82,6*0,75,30,2*0,11,2*0,
     2 90,3*0,46,0,41/
      DATA (MASH2(I),I=112,205)/2*0,22,7*0,57,0,14,45,3*0,60,5*0,40,
     1 2*0,10,2*0,81,8*0,69,9*0,62,5*0,12,2*0,50,2*0,31,0,28,4*0,
     2 86,10*0,61,3*0,3,3*0,2,64,5*0,38/
      DATA (MASH2(I),I=206,288)/0,72,84,3*0,20,3*0,80,0,26,3*0,56,6*0,
     1 25,5*0,59,0,93,42,48,2*0,44,5*0,58,0,89,2*0,78,76,2*0,98,
     2 4,3*0,94,6*0,49,15*0,36,47,0,67/
      DATA (MASH2(I),I=289,418)/0,100,3*0,83,7*0,71,2*0,77,5*0,17,97,
     1 7*0,96,10*0,13,3*0,87,2*0,27,0,95,4*0,68,8*0,99,10*0,24,
     2 6*0,63,0,55,35,9*0,18,6*0,29,0,33,8*0,85,8*0,79,5*0,66/
      DATA (ATWTD(K1),K1=1,60)/
     1    1.00794D0,       4.002602D0,      6.941D0,         9.012182D0,
     2   10.811D0,        12.011D0,        14.00674D0,      15.9994D0,
     3   18.9984032D0,    20.1797D0,       22.989768D0,     24.3050D0,
     4   26.981539D0,     28.0855D0,       30.973762D0,     32.066D0,
     5   35.4527D0,       39.948D0,        39.0983D0,       40.078D0,
     6   44.955910D0,     47.88D0,         50.9415D0,       51.9961D0,
     7   54.93805D0,      55.847D0,        58.93320D0,      58.69D0,
     8   63.546D0,        65.39D0,         69.723D0,        72.61D0,
     9   74.92159D0,      78.96D0,         79.904D0,        83.80D0,
     1   85.4678D0,       87.62D0,         88.90585D0,      91.224D0,
     2   92.90638D0,      95.94D0,         97.9072D0,      101.07D0,  
     3  102.9055D0,      106.42D0,        107.8682D0,      112.411D0,
     4  114.82D0,        118.710D0,       121.75D0,        127.60D0,
     5  126.90447D0,     131.29D0,        132.90543D0,     137.327D0,
     6  138.9055D0,      140.115D0,       140.90765D0,     144.24D0/
      DATA (ATWTD(K1),K1=61,100)/
     1  144.9127D0,      150.36D0,        151.965D0,       157.25D0,
     2  158.92534D0,     162.50D0,        164.93032D0,     167.26D0,
     3  168.93421D0,     173.04D0,        174.967D0,       178.49D0,
     4  180.9479D0,      183.85D0,        186.207D0,       190.2D0,
     5  192.22D0,        195.08D0,        196.96654D0,     200.59D0,
     6  204.3833D0,      207.2D0,         208.98037D0,     208.9824D0,
     7  209.9871D0,      222.0176D0,      223.0197D0,      226.0254D0,
     8  227.0278D0,      232.0381D0,      231.03588D0,     238.0289D0,
     9  237.0482D0,      239.0522D0,      243.0614D0,      247.0703D0,
     1  247.0703D0,      251.0796D0,      252.083D0,       257.0951D0/
      DATA POTH/19.2,41.8,40.0,63.7,76.0,78.0,82.0,95.0,
     1  115.,137.,149.,156.,166.,173.,173.,180.,159.29,
     2  188.,190.,191.,216.,233.,245.,257.,272.,286.,
     3  297.,311.,322.,330.,334.,350.,347.,348.,357.,
     4  352.,363.,366.,379.,393.,417.,424.,428.,441.,
     5  449.,470.,470.,469.,488.,488.,487.,485.,491.,
     6  482.,488.,491.,501.,523.,535.,546.,560.,574.,
     7  580.,591.,614.,628.,650.,658.,674.,684.,694.,
     8  705.,718.,727.,736.,746.,757.,790.,790.,800.,
     9  810.,823.,823.,830.,825.,794.,827.,826.,841.,
     1  847.,878.,890.,902.,921.,934.,939.,952.,966.,
     2  980.,994./
      DATA POTGAS/19.2,41.8,34.0,38.6,49.0,70.0,82.0,
     3  97.0,115.0/
      DATA POTCON/19.2,41.8,45.2,72.0,85.9,81.0,82.0,
     4  106.0,112.0/
      DO 116 L=1,72
      IC(L)=ICHAR(W(L:L))
      IF(IC(L)-32)101,102,103
  101 K(L)=1
      GO TO 116
  102 K(L)=2
      GO TO 116
  103 IF(IC(L)-48)104,105,105
  104 K(L)=1
      GO TO 116
  105 IF(IC(L)-58)106,107,107
  106 K(L)=3
      GO TO 116
  107 IF(IC(L)-65)108,109,109
  108 K(L)=1
      GO TO 116
  109 IF(IC(L)-91)110,111,111
  110 K(L)=4
      GO TO 116
  111 IF(IC(L)-97)112,113,113
  112 K(L)=1
      GO TO 116
  113 IF(IC(L)-123)114,115,115
  114 K(L)=5
      GO TO 116
  115 K(L)=1
  116 CONTINUE
      L=1
      M=0
  117 IF(K(L)-2)118,118,119
  118 L=L+1
      GO TO 117
  119 LMIN=L
  120 KG=K(L)
      IF(L-LMIN)130,130,140
  130 GO TO (150,150,150,160,150), KG
  140 GO TO (150,470,150,160,150), KG
  150 STOP 5
  160 KG1=K(L+1)
      GO TO (170,180,180,180,240), KG1
  170 STOP 6
  180 ICC=IC(L)-64
      JT=MASH1(ICC)
      IF(JT)190,190,200
  190 STOP 7
  200 M=M+1
      JZ(M)=JT
      GO TO (170,210,230,220,240), KG1
  210 NZ(M)=1
      GO TO 470
  220 NZ(M)=1
      L=L+1
      GO TO 120
  230 IN=L+1
      GO TO 390
  240 ICC=9*IC(L+1)-10*IC(L)+9
      IF(ICC-1)310,250,250
  250 IF(ICC-418)260,260,310
  260 IF(ICC-208)300,270,300
  270 M=M+1
      IF(IC(L)-71)290,280,290
  280 JZ(M)=32
      GO TO 330
  290 JZ(M)=84
      GO TO 330
  300 JT=MASH2(ICC)
      IF(JT)310,310,320
  310 STOP 8
  320 M=M+1
      JZ(M)=JT
  330 KG2=K(L+2)
      GO TO (340,350,380,360,370), KG2
  340 STOP 9
  350 NZ(M)=1
      GO TO 470
  360 NZ(M)=1
      L=L+2
      GO TO 120
  370 STOP 10
  380 IN=L+2
  390 INN=IN
      IS=0
      NZ(M)=0
  400 IF(K(INN)-3)420,410,420
  410 IS=IS+1
      MS(IS)=IC(INN)-48
      INN=INN+1
      GO TO 400
  420 ISM=IS
      KFAC=1
  430 NZ(M)=NZ(M)+KFAC*MS(IS)
      KFAC=10*KFAC
      IS=IS-1
      IF(IS)440,440,430
  440 IF(NZ(M))450,450,460
  450 STOP 12
  460 L=IN+ISM
      GO TO 120
  470 MMAX=M
      ASUM=0.0
      DO 480 M=1,MMAX
      JM=JZ(M)
  480 ASUM=ASUM+ATWTD(JM)*DBLE(NZ(M))
      DO 490 M=1,MMAX
      JM=JZ(M)
  490 WT(M)=ATWTD(JM)*DBLE(NZ(M))/ASUM
      ZAV=0.0
      POTL=0.0
      DO 540 M=1,MMAX
      JM=JZ(M)
      ZA=DBLE(JZ(M))/ATWTD(JM)
      ZAV=ZAV+WT(M)*ZA
      IF(KMAT-1)492,492,495
  492 POTM=POTH(JM)
      GO TO 540
  495 IF(JM-10)500,530,530
  500 GO TO (510,520), NPHAS
  510 POTM=POTGAS(JM)
      GO TO 540
  520 POTM=POTCON(JM)
      GO TO 540
  530 POTM=1.13*POTH(JM)
  540 POTL=POTL+WT(M)*ZA*LOG(POTM)
      POT=EXP(POTL/ZAV)
      RETURN
      END
      SUBROUTINE INDEX (L,TAG)
      CHARACTER TAG*3,IND(0:9)*1
      DATA IND/'0','1','2','3','4','5','6','7','8','9'/
      I1=L/100
      M=L-100*I1
      I2=M/10
      I3=M-10*I2
      TAG=IND(I1)//IND(I2)//IND(I3)
      RETURN
      END
      SUBROUTINE SCOF(X,F,NMAX,A,B,C,D) 
C     22 Feb 83
C  IF S LIES BETWEEN X(M) AND X(M+1), THEN
C  F(S)=((D(M)*S+C(M))*S+B(M))*S+A(M) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(1000),F(1000),A(1000),B(1000),C(1000),D(1000)
      M1=2
      M2=NMAX-1 
      S=0.0 
      DO 10 M=1,M2
      D(M)=X(M+1)-X(M)
      R=(F(M+1)-F(M))/D(M)
      C(M)=R-S
   10 S=R 
      S=0.0 
      R=0.0 
      C(1)=0.0
      C(NMAX)=0.0 
      DO 20 M=M1,M2 
      C(M)=C(M)+R*C(M-1)
      B(M)=(X(M-1)-X(M+1))*2.0-R*S
      S=D(M)
   20 R=S/B(M)
      MR=M2 
      DO 30 M=M1,M2 
      C(MR)=(D(MR)*C(MR+1)-C(MR))/B(MR) 
   30 MR=MR-1 
      DO 40 M=1,M2
      S=D(M)
      R=C(M+1)-C(M) 
      D(M)=R/S
      C(M)=C(M)*3.0 
      B(M)=(F(M+1)-F(M))/S-(C(M)+R)*S 
   40 A(M)=F(M) 
      RETURN
      END 
      SUBROUTINE BSPOL(S,X,A,B,C,D,N,G)
C     22 FEB 83
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(1000),A(1000),B(1000),C(1000),D(1000)
      IF (X(1).GT.X(N)) GO TO 10
      IDIR=0
      MLB=0 
      MUB=N 
      GO TO 20
   10 IDIR=1
      MLB=N 
      MUB=0 
   20 IF (S.GE.X(MUB+IDIR)) GO TO 60
      IF (S.LE.X(MLB+1-IDIR)) GO TO 70
      ML=MLB
      MU=MUB
      GO TO 40
   30 IF (IABS(MU-ML).LE.1) GO TO 80
   40 MAV=(ML+MU)/2 
      IF (S.LT.X(MAV)) GO TO 50 
      ML=MAV
      GO TO 30
   50 MU=MAV
      GO TO 30
   60 MU=MUB+2*IDIR-1 
      GO TO 90
   70 MU=MLB-2*IDIR+1 
      GO TO 90
   80 MU=MU+IDIR-1
   90 Q=S-X(MU) 
      G=((D(MU)*Q+C(MU))*Q+B(MU))*Q+A(MU) 
      RETURN
      END 
      SUBROUTINE GRAL(DELTA,G,N,RESULT)
C     SUBROUTINE GRAL(DELTA,G,N,RESULT); 5 MAY 86; FORMERLY INT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION G(10001)
      NL1=N-1 
      NL2=N-2 
      IF (FLOAT (N) -2.0*FLOAT (N/2)) 100,100,10
C     IF N IS ODD, GO TO 10 - IF N IS EVEN, GO TO 100 
   10  IF (N-1) 15,15,20
   15 SIGMA=0.0 
      GO TO 70
   20 IF(N-3) 30,30,40
   30 SIGMA=G(1)+4.0*G(2)+G(3)
      GO TO 70
   40 SUM4=0.0
      DO 50 K=2,NL1,2 
   50 SUM4=SUM4+G(K)
      SUM2=0.0
      DO 60 K=3,NL2,2 
   60 SUM2=SUM2+G(K)
      SIGMA=G(1)+4.0*SUM4+2.0*SUM2+G(N) 
   70 RESULT=DELTA*SIGMA
      RETURN
  100 IF(N-2)110,110,120
  110 SIGMA=1.5*(G(1)+G(2)) 
      GO TO 70
  120 IF(N-4)130,130,140
  130 SIGMA=1.125*(G(1)+3.0*G(2)+3.0*G(3)+G(4)) 
      GO TO 70
  140 IF(N-6)150,150,160
  150 SIGMA=G(1)+3.875*G(2)+2.625*G(3)+2.625*G(4)+3.875*G(5)+G(6) 
      GO TO 70
  160 IF (N-8)170,170,180 
  170 SIGMA=G(1)+3.875*G(2)+2.625*G(3)+2.625*G(4)+3.875*G(5)+2.0*G(6) 
     1+4.0*G(7)+G(8)
      GO TO 70
  180 SIG6=G(1)+3.875*G(2)+2.625*G(3)+2.625*G(4)+3.875*G(5)+G(6)
      SUM4=0.0
      DO 190 K=7,NL1,2
  190 SUM4=SUM4+G(K)
      SUM2=0.0
      DO 200 K=8,NL2,2
  200 SUM2=SUM2+G(K)
      SIGMA=SIG6+G(6)+4.0*SUM4+2.0*SUM2+G(N)
      GO TO 70
      END
