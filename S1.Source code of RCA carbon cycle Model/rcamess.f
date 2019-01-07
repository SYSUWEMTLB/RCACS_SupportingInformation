      SUBROUTINE RCAMESS(I,DATA)
C 
C        RCAMESS IS USED TO INFORM USER CONCERNING CHANGES IN INTEGRATION
C                PROCEDURE AND/OR STABILITY DURING THE SIMULATION 
C 
      SAVE
      INCLUDE 'RCACM' 

      GO TO (10,20,30),I

   10 WRITE(OUT,1000)  TIME,ISYS,IX,IY,IZ,CARAY(IX,IY,IZ,ISYS),DATA
     .   ,CDARAY(IX,IY,IZ,ISYS) 
 1000 FORMAT(/////12X,'STABILITY CRITERIA VIOLATED AT TIME',F9.3,' IN SY
     .STEM',I4,' SEGMENT',3I4 / 4X,'C =',E15.6,' EXCEEDED CMAX =',E15.6, 
     .  5X,'DERIVATIVE AT PREVIOUS STEP =',E15.6// 19X,'COMPUTATION DISC
     .ONTINUED, DUMPS UNTIL LAST PRINTOUT WILL FOLLOW') 
      RETURN

c  20 WRITE(OUT,2000)   ISYS,TIME,DATA
   20 CONTINUE
c     print *,'   IN,OUT,MXACTS,ISYS,IX,IY,IZ,LIST(5)'
c     print *,   IN,OUT,MXACTS,ISYS,IX,IY,IZ,LIST
      
      IDAY=INT(ITIMESECS/86400)
      IHOUR=INT(ITIMESECS/3600)-IDAY*24
      IMINUTE=INT(ITIMESECS/60)-IDAY*1440-IHOUR*60
      ISECOND=INT(ITIMESECS)-IDAY*86400-IHOUR*3600-IMINUTE*60
      WRITE(OUT,2000)   ISYS,IDAY,IHOUR,IMINUTE,ISECOND,CMIN(ISYS),DATA
 2000 FORMAT(/39X,'NEGATIVE CONCS COMPUTED IN SYSTEM',I3,' AT TIME =',
     .   I3,':',I2,':',I2,':',I2,' Days:Hrs:Mins:Secs'
     .   /22X,'CONC ADJUSTED TO CMIN=',E10.3,' WHICH RESULTED IN',
     .   E10.3,' KGS BEING ADDED TO THE SYSTEM'/30X,'THE FOLLOWING SEGME
     .NTS ADJUSTED TO CMIN...')
C2000 FORMAT(/39X,'NEGATIVE CONCS COMPUTED IN SYSTEM',I3,' AT TIME =',
C    .   F7.3/22X,'QUARTER CONC ADJUSTMENTS MADE WHICH RESULTED IN',
C    .   E10.3,' KGS BEING ADDED TO THE SYSTEM'/30X,'THE FOLLOWING SEGME
C    .NTS RECEIVED THE QUARTER CONC ADJUSTMENTS...')
      RETURN

   30 WRITE(OUT,3000)  TIME,DATA
 3000 FORMAT(/20X,'ERROR AT TIME = ',F8.3,' DAYS  INVALID INTEGRATION ST
     .EPSIZE SPECIFIED, DT =',E13.4/20X,'RCA WILL BE TERMINATED')
      RETURN
      END 
