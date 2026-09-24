      PROGRAM EDCONV
C
C     8 September 1992.
C
C     Written by Martin J. Berger,   
C     National Institute of Standards and Technology,
C     Gaithersburg, MD 20899.
C
C     Converts formatted electron data file FEDAT into
C     binary direct-access file UEDAT.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION NC(26),BD(26),RLOS(113)
      OPEN (7,FILE='FEDAT')
      OPEN (8,FILE='UEDAT',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=1224)
      DO 10 K=1,100
      READ (7,*) NMAX,LKMAX
      READ (7,*) (NC(N),N=1,NMAX)
      READ (7,*) (BD(N),N=1,NMAX)
      READ (7,*) (RLOS(L),L=1,LKMAX)
      WRITE (8,REC=K) NMAX,LKMAX,NC,BD,RLOS
   10 CONTINUE
      STOP
      END
