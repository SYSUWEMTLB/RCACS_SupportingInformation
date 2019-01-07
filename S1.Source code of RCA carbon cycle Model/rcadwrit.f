      SUBROUTINE  RCAWBUF(ISYSTM,A1,A2,A3,A4,A5)
C
C        RCAWBUF/RCAWRIT IS USED TO WRITE RCA DUMPS TO DISK
C                        IT REDUCES DISK I/O BY PROVIDING AN INTERNAL
C                        BUFFER TO HOLD INTERMEDIATE VALUES
      SAVE
      INCLUDE  'RCACM'      
      
      REAL  BUF1(5000),BUF2(10000),BUF3(20000),BUF4(30000),BUF5(50000),
     .BUF6(80000),BUF7(100000),BUF8(130000),BUF9(150000),BUFFER(200000)
      
      EQUIVALENCE 
     .   (BUFFER(1),BUF1(1)) , (BUFFER(1),BUF2(1)) , (BUFFER(1),BUF3(1))
     . , (BUFFER(1),BUF4(1)) , (BUFFER(1),BUF5(1)) , (BUFFER(1),BUF6(1))
     . , (BUFFER(1),BUF7(1)) , (BUFFER(1),BUF8(1)) , (BUFFER(1),BUF9(1))

      DATA  INITBUF/0/,ICNT/0/
      IF(INITBUF.EQ.0) THEN
       DO IB=1,200000
        BUFFER(IB)=0.0
       ENDDO
       INITBUF=1
      ENDIF

      BUFFER(ICNT+1) = A1
      BUFFER(ICNT+2) = A2
      BUFFER(ICNT+3) = A3
      BUFFER(ICNT+4) = A4
      BUFFER(ICNT+5) = A5
      ICNT = ICNT + 5
      IF(ICNT.GT.200000)  THEN
         WRITE(OUT,1200)   ICNT
 1200    FORMAT(//10X,'RCADWRIT ERROR...','ICNT = ',I10/
     .       10X,'MAXIMUM BUFFER SIZE IS EXCEEDED'/
     .       10X,'(PRODUCT OF NDMPS*NOSYS*5 CANNOT EXCEED 200000'/
     .       10X,'RCA PROGRAM TERMINATED'//)
       CALL EXIT
      ENDIF
      RETURN

      ENTRY RCAWRIT

C        BUFFER IS COMPLETE - DUMP TO DISK

      IF(IDDOPT.EQ.0) THEN
        WRITE(12)   TIME
      ELSE
        IF(INITB.EQ.0)  WRITE(12)   TIME
        IF(INITB.GE.1)  WRITE(12)   TIME-(FLOAT(IPRNTDSECS)/86400.)/2.
      ENDIF 
      IF(ICNT.LE.5000)   WRITE(13) BUF1
      IF(ICNT.GT.5000 .AND. ICNT.LE.10000)   WRITE(13) BUF2   
      IF(ICNT.GT.10000 .AND. ICNT.LE.20000)   WRITE(13) BUF3
      IF(ICNT.GT.20000 .AND. ICNT.LE.30000)   WRITE(13) BUF4
      IF(ICNT.GT.30000 .AND. ICNT.LE.50000)   WRITE(13) BUF5
      IF(ICNT.GT.50000 .AND. ICNT.LE.80000)   WRITE(13) BUF6
      IF(ICNT.GT.80000 .AND. ICNT.LE.100000)   WRITE(13) BUF7
      IF(ICNT.GT.100000 .AND. ICNT.LE.130000)   WRITE(13) BUF8
      IF(ICNT.GT.130000 .AND. ICNT.LE.150000)   WRITE(13) BUF9
      IF(ICNT.GT.150000)   WRITE(13) BUFFER
C        RESET ICNT
      RECSIZ = ICNT
      ICNT = 0
      RETURN
      END
