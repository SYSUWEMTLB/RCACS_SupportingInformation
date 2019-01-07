      SUBROUTINE RCA10(SVALS,BVALS,NDIM,NO,NXT,IFILE,IPWLOPT,SCALE)
C 
C        RCA10 UPDATES THE FORCING FUNCTIONS
C 

      SAVE
      INCLUDE 'RCACM' 
      REAL      NXT,SVALS(NDIM,NOSYS),BVALS(NDIM,NOSYS),SCALE(NOSYS)
      INTEGER   NO(NOSYS),IBINOPT

      WRITE(OUT,8000)
 8000 FORMAT(//)
      IF (IFILE.EQ.33) THEN
        IBINOPT = IBNRYRDOPTS(2)
        ISCALE = ISCALPS
      ELSEIF (IFILE.EQ.34) THEN
        IBINOPT = IBNRYRDOPTS(3)
        ISCALE = ISCALNPS
      ELSEIF (IFILE.EQ.35) THEN
        IBINOPT = IBNRYRDOPTS(4)
        ISCALE = ISCALFL
      ELSEIF (IFILE.EQ.36) THEN
        IBINOPT = IBNRYRDOPTS(5)
        ISCALE = ISCALATM
      ELSE
        PRINT *,'*** Error in rca10.f - Invalid file number: IFILE = ',
     .           IFILE
        STOP
      ENDIF

       IF(IPWLOPT.EQ.1)  THEN
        OLDTIME=NXT
        IF(IBINOPT.EQ.0)  THEN
         READ(IFILE,2700,ERR=950,END=970)  NXTSECS
 2700    FORMAT(10X,I10)
        ELSE
         READ(IFILE,ERR=960,END=970)  NXTSECS
        ENDIF
        NXTSECS=ISCALE*NXTSECS
        NXT=NXTSECS/86400.
       ENDIF

      DO 100 ISYS=1,NOSYS 
       IF(NO(ISYS).EQ.0)   GO TO 100
       IF(IFILE.EQ.33)  THEN
        IF(IPWLOPT.EQ.0) THEN
          WRITE(OUT,2601)   ISYS,NXT
 2601     FORMAT(30X,'POINT SOURCE LOADINGS FOR SYSTEM',I4,' AT TIME ='
     .       ,F8.2,' DAYS')
          IF(LIST(3).EQ.1) WRITE(OUT,2701)
 2701     FORMAT(/12X,8(4X,'PS(T)',4X))
        ELSE
          IF(INITB.NE.0) THEN
            WRITE(OUT,2601)   ISYS,OLDTIME
            IF(LIST(3).EQ.1) WRITE(OUT,2701)
          ENDIF
        ENDIF
       ENDIF
       IF(IFILE.EQ.34)  THEN
        IF(IPWLOPT.EQ.0) THEN
          WRITE(OUT,2602)   ISYS,NXT
 2602     FORMAT(30X,'NONPOINT SOURCE LOADINGS FOR SYSTEM',I4,
     .       ' AT TIME =',F8.2,' DAYS')
          IF(LIST(3).EQ.1) WRITE(OUT,2702)
 2702     FORMAT(/12X,8(4X,'NPS(T)',3X))
        ELSE
          IF(INITB.NE.0) THEN
            WRITE(OUT,2602)   ISYS,OLDTIME
            IF(LIST(3).EQ.1) WRITE(OUT,2702)
          ENDIF
        ENDIF
       ENDIF
       IF(IFILE.EQ.35)  THEN
        IF(IPWLOPT.EQ.0) THEN
          WRITE(OUT,2603)   ISYS,NXT
 2603     FORMAT(30X,'FALL-LINE LOADINGS FOR SYSTEM',I4,' AT TIME ='
     .       ,F8.2,' DAYS')
          IF(LIST(3).EQ.1) WRITE(OUT,2703)
 2703     FORMAT(/12X,8(4X,'FL(T)',4X))
        ELSE
          IF(INITB.NE.0) THEN
            WRITE(OUT,2603)   ISYS,OLDTIME
            IF(LIST(3).EQ.1) WRITE(OUT,2703)
          ENDIF
        ENDIF
       ENDIF
       IF(IFILE.EQ.36)  THEN
        IF(IPWLOPT.EQ.0) THEN
          WRITE(OUT,2604)   ISYS,NXT
 2604     FORMAT(30X,'ATMOSPHERIC LOADINGS FOR SYSTEM',I4,' AT TIME ='
     .       ,F8.2,' DAYS')
        ELSE
          IF(INITB.NE.0) THEN
            WRITE(OUT,2604)   ISYS,OLDTIME
          ENDIF
        ENDIF
       ENDIF

      IF(IPWLOPT.EQ.1)  THEN
       IF(IFILE.NE.36)  THEN 
         DO I=1,NO(ISYS)
          SVALS(I,ISYS)=BVALS(I,ISYS)
         ENDDO
       ELSE
         SVALC=0.0
         DO IY=1,NY
          DO IX=1,NX
           IPOS = NX*(IY-1)+IX
           SVALS(IPOS,ISYS) = BVALS(IPOS,ISYS)
           IF(BVALS(IPOS,ISYS).NE.0.0) SVALC=BVALS(IPOS,ISYS)
          ENDDO
         ENDDO
       ENDIF
      ENDIF

      IF (IBINOPT .EQ. 0) THEN
        IF(IFILE.NE.36) THEN
         READ(IFILE,2750,ERR=980,END=970) (BVALS(I,ISYS),I=1,NO(ISYS))
 2750    FORMAT(10X,7F10.0)
        ELSE
         IF(NO(ISYS).EQ.1)  THEN
          READ(IFILE,2750,ERR=960,END=970)  BVALC
          DO IY=1,NY
           DO IX=1,NX
            IPOS = NX*(IY-1)+IX
            BVALS(IPOS,ISYS) = BVALC
           ENDDO
          ENDDO
         ELSE
          DO IY=1,NY
           IPOS = NX*(IY-1)
           READ(IFILE,2750,ERR=990,END=970)  (BVALS(IPOS+IX,ISYS),
     .          IX=1,NX)
           DO IX=1,NX
            IPOS = NX*(IY-1)+IX
            IF(FSM(IX,IY).NE.1.)  BVALS(IPOS,ISYS)=0.0
           ENDDO
          ENDDO
         ENDIF
        ENDIF

      ELSE

        IF(IFILE.NE.36) THEN
         READ(IFILE,ERR=960,END=970) (BVALS(I,ISYS),I=1,NO(ISYS))
        ELSE
         IF(NO(ISYS).EQ.1)  THEN
          READ(IFILE,ERR=960,END=970)  BVALC
          DO IY=1,NY
           DO IX=1,NX
            IPOS = NX*(IY-1)+IX
            BVALS(IPOS,ISYS) = BVALC
           ENDDO
          ENDDO
         ELSE
          DO IY=1,NY
           IPOS = NX*(IY-1)
           READ(IFILE,ERR=960,END=970)  (BVALS(IPOS+IX,ISYS),IX=1,NX)
           DO IX=1,NX
            IPOS = NX*(IY-1)+IX
            IF(FSM(IX,IY).NE.1.)  BVALS(IPOS,ISYS)=0.0
           ENDDO
          ENDDO
         ENDIF
        ENDIF

      ENDIF

      IF(LIST(3).EQ.1)  THEN
        IF(IFILE.NE.36)  THEN
          IF(IPWLOPT.EQ.0) THEN
           WRITE(OUT,2800) (BVALS(I,ISYS),I=1,NO(ISYS))
 2800      FORMAT(10X,8E13.3)
          ELSE
           IF(INITB.NE.0)
     .         WRITE(OUT,2800) ((SVALS(I,ISYS)/1000.),I=1,NO(ISYS))
          ENDIF
        ELSE
          IF(NO(ISYS).EQ.1) THEN
           IF(IPWLOPT.EQ.0) WRITE(OUT,2800)  BVALC
           IF(IPWLOPT.EQ.1 .AND. INITB.NE.0) WRITE(OUT,2800) SVALC/1000.
          ELSE
           DO IY=1,NY
            IPOS = NX*(IY-1)
            IF(IPWLOPT.EQ.0)
     .          WRITE(OUT,2850)  IY,(BVALS(IPOS+IX,ISYS),IX=1,NX)
            IF(IPWLOPT.EQ.1 .AND. INITB.NE.0)
     .          WRITE(OUT,2850) IY,((SVALS(IPOS+IX,ISYS)/1000.),IX=1,NX)
 2850       FORMAT(I5,8E13.3/(5X,8E13.3))
           ENDDO
          ENDIF
        ENDIF
      ENDIF

C     Multiply by scale factor and convert to MG*M^3/DAY-L

      TOTLD = 0.0
      IF(IFILE.NE.36)  THEN
        DO 40 I=1,NO(ISYS)
         BVALS(I,ISYS) = SCALE(ISYS)*BVALS(I,ISYS)*1000.
         IF(IPWLOPT.EQ.0) THEN
           TOTLD = TOTLD + BVALS(I,ISYS)
         ELSE
           TOTLD = TOTLD + SVALS(I,ISYS)
         ENDIF
   40   CONTINUE
      ELSE
        IF(IFILE.EQ.36)  THEN
         TOTLD = 0.0
         DO 50 IY=1,NY
          DO 50 IX=1,NX
           IPOS = NX*(IY-1)+IX
           BVALS(IPOS,ISYS) = 1000.*SCALE(ISYS)*BVALS(IPOS,ISYS)
           IF(IPWLOPT.EQ.0) THEN
             TOTLD = TOTLD + XAZ(IX,IY)*BVALS(IPOS,ISYS)
           ELSE
             TOTLD = TOTLD + XAZ(IX,IY)*SVALS(IPOS,ISYS)
           ENDIF
   50    CONTINUE
        ENDIF
      ENDIF
  
      IF(IPWLOPT.EQ.0 .OR. INITB.NE.0) WRITE(OUT,2900)  ISYS,TOTLD/1000.
 2900 FORMAT(35X,'TOTAL LOADING FOR SYSTEM',I3,' =',E13.4,' KG/DAY'/)

      IF(IPWLOPT.EQ.1)  THEN
        IF(IFILE.NE.36)  THEN
         DO I=1,NO(ISYS)
          SVALS(I,ISYS)=(SVALS(I,ISYS)-BVALS(I,ISYS))/(OLDTIME-NXT)
         ENDDO
        ELSE
         DO IY=1,NY
          DO IX=1,NX
           IPOS = NX*(IY-1)+IX
           SVALS(IPOS,ISYS)=(SVALS(IPOS,ISYS)-BVALS(IPOS,ISYS))/
     .                           (OLDTIME-NXT)
          ENDDO
         ENDDO
        ENDIF
      ENDIF

  100 CONTINUE

      IF(IPWLOPT.EQ.0)  THEN
       IF(IBINOPT.EQ.0)  THEN
        READ(IFILE,2700,ERR=950,END=970)  NXTSECS
       ELSE
        READ(IFILE,ERR=960,END=970)  NXTSECS
       ENDIF
       NXTSECS=ISCALE*NXTSECS
       NXT=NXTSECS/86400.
      ENDIF

      RETURN

  950 IN=IFILE
      WRITE(OUT,9500)
 9500 FORMAT(///20X,'ERROR ... INPUT ERROR ENCOUNTERED WHILE ATTEMPTING 
     .TO READ NEXT TIME FOR LOADING FILE'/20X,'RCA TERMINATED'//)
      CALL FMTER
      CALL EXIT
  960 WRITE(OUT,9600)  IFILE,ISYS,ITIMESECS
 9600 FORMAT(///20X,'ERROR READING BINARY FILE NUMBER ',I3,' FOR SYSTEM'
     .   I3,' AT TIME =',I10,' SECONDS (',E10.3,' DAYS)')
      CALL EXIT
  970 WRITE(OUT,9700)
 9700 FORMAT(//5X,
     .  'END OF FILE ENCOUNTERED WHILE READING LOADING FILE'//)
      INPCHCK=2
      RETURN
  980 WRITE(OUT,9800)  ISYS
 9800 FORMAT(///20X,'ERROR ... INPUT ERROR ENCOUNTERED WHILE ATTEMPTING
     .TO READ LOADS FOR SYSTEM',I3/20X,'RCA TERMINATED')
      CALL FMTER
      CALL EXIT
  990 WRITE(OUT,9900)  ISYS
 9900 FORMAT(///20X,'ERROR ... INPUT ERROR ENCOUNTERED WHILE ATTEMPTING
     .TO READ ATMOSPHERIC LOADS FOR SYSTEM',I3/20X,'RCA TERMINATED')
      CALL FMTER
      CALL EXIT
      END 
