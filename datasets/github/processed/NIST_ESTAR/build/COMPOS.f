      PROGRAM COMP
C
C     8 September 1992.
C
C     Written by Martin J. Berger,   
C     National Institute of Standards and Technology,
C     Gaithersburg, MD 20899.
C
C     Minor changes made in June of 1999 for web version 1.2
C     These included: file open routines that
C     were not supported by our compiler (Johnathan S. Coursey)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION MZ(14),WT(14)
      CHARACTER MAT*72,OUTPUT*30
    1 FORMAT(1H )
    5 FORMAT(A)
    6 FORMAT(1X,A)
      PRINT *,' Enter ID number: '
      READ *, ID
C Option to print to the screen was removed because it wasn't working
      PRINT *,' Enter name of output file: '
      READ 5, OUTPUT
      OPEN (8,FILE=OUTPUT)
      IF(ID-278)50,50,10
   10 IF(ID-906)30,20,30
   20 ID=279 
      GO TO 50
   30 PRINT 40, ID
   40 FORMAT(I6,' is not an allowed ID number.')
      STOP 
   50 OPEN (7,FILE='UCOMP',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=268)
      READ (7,REC=ID) MAT,MMAX,ZAG,POT,RHO,MZ,WT
      WRITE (8,55) ID
   55 FORMAT(' ID number ',I3)
      WRITE (8,1)
      WRITE (8,6) MAT
      WRITE (8,1)
      WRITE (8,60) RHO
   60 FORMAT('             Density (g/cm3) = ',F12.5) 
      WRITE (8,1)
      WRITE (8,70) POT
   70 FORMAT(' Mean Excitation Energy (eV) = ',F8.1)
      WRITE (8,1)
      WRITE (8,80)
   80 FORMAT(' COMPOSITION:')
      WRITE (8,90)
   90 FORMAT('        Z = Atomic number')
      WRITE (8,100)
  100 FORMAT('       WT = fraction by weight')
      WRITE (8,1)
      WRITE (8,110)
  110 FORMAT(15X,'Z',8X,'WT')  
      DO 130 M=1,MMAX
      WRITE (8,120) MZ(M),WT(M)
  120 FORMAT(10X,I6,F10.6)
  130 CONTINUE
      STOP
      END
      
