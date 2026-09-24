      PROGRAM CONVERT
C
C     8 September 1992.
C
C     Written by Martin J. Berger,   
C     National Institute of Standards and Technology,
C     Gaithersburg, MD 20899.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION MZ(14),WT(14)
      CHARACTER MAT*72
      DATA KMAX/279/
    5 FORMAT(1X,A)
      OPEN (7,FILE='FCOMP')
      OPEN (8,FILE='UCOMP',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=268)
      DO 20 K=1,KMAX
      READ (7,5) MAT
      READ (7,*) MMAX,ZAG,POT,RHO
      READ (7,10) (MZ(M),WT(M),M=1,MMAX)
   10 FORMAT(6(I3,F9.6))
      WRITE (8,REC=K) MAT,MMAX,ZAG,POT,RHO,MZ,WT
   20 CONTINUE
      STOP
      END
