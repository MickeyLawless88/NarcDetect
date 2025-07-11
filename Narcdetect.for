C.##....##....###....########...######..########..########.########.########..######..########
C.###...##...##.##...##.....##.##....##.##.....##.##..........##....##.......##....##....##...
C.####..##..##...##..##.....##.##.......##.....##.##..........##....##.......##..........##...
C.##.##.##.##.....##.########..##.......##.....##.######......##....######...##..........##...
C.##..####.#########.##...##...##.......##.....##.##..........##....##.......##..........##...
C.##...###.##.....##.##....##..##....##.##.....##.##..........##....##.......##....##....##...
C.##....##.##.....##.##.....##..######..########..########....##....########..######.....##...
C
C     DRUG DETECTION TIME CALCULATOR FOR ORAL FLUID TESTS
C     ANSI FORTRAN IV PROGRAM
C     CALCULATES DETECTION WINDOW IN SECONDS FOR VARIOUS DRUGS
C     INCLUDES ACCUMULATION EFFECTS FROM CHRONIC USE
C     ENHANCED WITH ROUTES OF ADMINISTRATION VARIABLES
C     
      PROGRAM DRUGDT
      IMPLICIT REAL (A-H,O-Z)
      INTEGER DRUG, DOSAGE, WEIGHT, AGE, METAB, ROUTE
      REAL HALFLF, DETECT, ELIMRT, FACTOR, FENTDS, USEDUR, ACCUMF
      REAL STEADY, BUILDUP, TOTCONC, BIOAVAIL, ABSORPT
      CHARACTER*20 DNAME, DINPUT, RNAME, RINPUT
      
C     SET FENTANYL DOSE CONSTANT (1.0 = 1 GRAM = 1000MG)
      FENTDS = 1.0
      
C     BANNER AND INTRODUCTION
      WRITE(6,100)
100   FORMAT('1','DRUG DETECTION TIME CALCULATOR',/,
     1       ' FOR ORAL FLUID (SALIVA) TESTING',//,
     2       ' ESTIMATES TIME UNTIL NON-DETECTABLE',/,
     3       ' BASED ON PHARMACOKINETIC PARAMETERS',/,
     4       ' INCLUDES CHRONIC USE ACCUMULATION',/,
     5       ' AND ROUTES OF ADMINISTRATION',//)
      
C     DRUG SELECTION MENU
      WRITE(6,101)
101   FORMAT(' AVAILABLE DRUGS',/,
     1       ' FENTANYL, NITAZENES, AMPHETAMINE,',/,
     2       ' METHAMPHETAMINE, DEXTROAMPHETAMINE,',/,
     3       ' HYDROMORPHONE, OXYCODONE, MORPHINE,',/,
     4       ' HYDROCODONE, CODEINE, PETHIDINE,',/,
     5       ' BARBITURATES, BENZODIAZEPINES, ALCOHOL,',/,
     6       ' LSD, KETAMINE, MESCALINE, PSILOCYBIN,',/,
     7       ' DMT, GHB, METHAQUALONE, METHADONE,',/,
     8       ' DEXTROPROPOXYPHENE',//,
     9       ' ENTER DRUG NAME ',$)
      READ(5,'(A)') DINPUT
      
C     CONVERT INPUT TO UPPERCASE FOR COMPARISON
      CALL UPPER(DINPUT)
      
C     DETERMINE DRUG NUMBER FROM NAME
      DRUG = 0
      IF (DINPUT.EQ.'FENTANYL') DRUG = 1
      IF (DINPUT.EQ.'NITAZENES') DRUG = 2
      IF (DINPUT.EQ.'AMPHETAMINE') DRUG = 3
      IF (DINPUT.EQ.'METHAMPHETAMINE') DRUG = 4
      IF (DINPUT.EQ.'DEXTROAMPHETAMINE') DRUG = 5
      IF (DINPUT.EQ.'HYDROMORPHONE') DRUG = 6
      IF (DINPUT.EQ.'OXYCODONE') DRUG = 7
      IF (DINPUT.EQ.'MORPHINE') DRUG = 8
      IF (DINPUT.EQ.'HYDROCODONE') DRUG = 9
      IF (DINPUT.EQ.'CODEINE') DRUG = 10
      IF (DINPUT.EQ.'PETHIDINE') DRUG = 11
      IF (DINPUT.EQ.'MEPERIDINE') DRUG = 11
      IF (DINPUT.EQ.'BARBITURATES') DRUG = 12
      IF (DINPUT.EQ.'BENZODIAZEPINES') DRUG = 13
      IF (DINPUT.EQ.'ALCOHOL') DRUG = 14
      IF (DINPUT.EQ.'ETHANOL') DRUG = 14
      IF (DINPUT.EQ.'LSD') DRUG = 15
      IF (DINPUT.EQ.'KETAMINE') DRUG = 16
      IF (DINPUT.EQ.'MESCALINE') DRUG = 17
      IF (DINPUT.EQ.'PSILOCYBIN') DRUG = 18
      IF (DINPUT.EQ.'DMT') DRUG = 19
      IF (DINPUT.EQ.'GHB') DRUG = 20
      IF (DINPUT.EQ.'METHAQUALONE') DRUG = 21
      IF (DINPUT.EQ.'METHADONE') DRUG = 22
      IF (DINPUT.EQ.'DEXTROPROPOXYPHENE') DRUG = 23
      IF (DINPUT.EQ.'PROPOXYPHENE') DRUG = 23
      
      IF (DRUG.EQ.0) THEN
        WRITE(6,*) 'INVALID DRUG NAME ', DINPUT
        WRITE(6,*) 'PLEASE CHECK SPELLING AND TRY AGAIN'
        STOP
      
C     ROUTE OF ADMINISTRATION SELECTION
C     ---------------------------------

      WRITE(6,108)
108   FORMAT(/,' AVAILABLE ROUTES OF ADMINISTRATION',/,
     1       ' ORAL, INTRAVENOUS, INTRAMUSCULAR,',/,
     2       ' SUBCUTANEOUS, INTRANASAL, INHALATION,',/,
     3       ' SUBLINGUAL, TRANSDERMAL, RECTAL,',/,
     4       ' BUCCAL, TOPICAL',//,
     5       ' ENTER ROUTE OF ADMINISTRATION ',$)
      READ(5,'(A)') RINPUT
      
C     CONVERT ROUTE INPUT TO UPPERCASE
C     --------------------------------

      CALL UPPER(RINPUT)
      
C     DETERMINE ROUTE NUMBER FROM NAME
C     ----------------------

      ROUTE = 0
      IF (RINPUT.EQ.'ORAL') ROUTE = 1
      IF (RINPUT.EQ.'INTRAVENOUS') ROUTE = 2
      IF (RINPUT.EQ.'IV') ROUTE = 2
      IF (RINPUT.EQ.'INTRAMUSCULAR') ROUTE = 3
      IF (RINPUT.EQ.'IM') ROUTE = 3
      IF (RINPUT.EQ.'SUBCUTANEOUS') ROUTE = 4
      IF (RINPUT.EQ.'SC') ROUTE = 4
      IF (RINPUT.EQ.'INTRANASAL') ROUTE = 5
      IF (RINPUT.EQ.'NASAL') ROUTE = 5
      IF (RINPUT.EQ.'INHALATION') ROUTE = 6
      IF (RINPUT.EQ.'INHALED') ROUTE = 6
      IF (RINPUT.EQ.'SMOKING') ROUTE = 6
      IF (RINPUT.EQ.'SUBLINGUAL') ROUTE = 7
      IF (RINPUT.EQ.'SL') ROUTE = 7
      IF (RINPUT.EQ.'TRANSDERMAL') ROUTE = 8
      IF (RINPUT.EQ.'PATCH') ROUTE = 8
      IF (RINPUT.EQ.'RECTAL') ROUTE = 9
      IF (RINPUT.EQ.'PR') ROUTE = 9
      IF (RINPUT.EQ.'BUCCAL') ROUTE = 10
      IF (RINPUT.EQ.'TOPICAL') ROUTE = 11
      
      IF (ROUTE.EQ.0) THEN
        WRITE(6,*) 'INVALID ROUTE NAME ', RINPUT
        WRITE(6,*) 'PLEASE CHECK SPELLING AND TRY AGAIN'
        STOP
      ENDIF
      
C     INPUT PARAMETERS
C     ----------------

      WRITE(6,102)
102   FORMAT(/,' ENTER DOSAGE IN MG ',$)
      READ(5,*) DOSAGE
      
      WRITE(6,103)
103   FORMAT(' ENTER BODY WEIGHT IN KG ',$)
      READ(5,*) WEIGHT
      
      WRITE(6,104)
104   FORMAT(' ENTER AGE IN YEARS ',$)
      READ(5,*) AGE
      
      WRITE(6,105)
105   FORMAT(' METABOLISM RATE (1=SLOW, 2=NORMAL, 3=FAST) ',$)
      READ(5,*) METAB
      
      WRITE(6,107)
107   FORMAT(' DURATION OF USE IN HOURS (24.0=1 DAY) ',$)
      READ(5,*) USEDUR
      
C     SET DRUG-SPECIFIC PARAMETERS USING COMPUTED GOTO
      GOTO (10,20,30,40,50,60,70,80,90,91,92,93,94,95,96,97,98,99,
     1      100,101,102,103,104), DRUG
      WRITE(6,*) 'DRUG SELECTION ERROR'
      STOP
      
C     FENTANYL PARAMETERS (FENTDS = 1.0 = 1000MG)

10    DNAME = 'FENTANYL'
      HALFLF = 2.0
      CUTOFF = 1.0
      DOSINT = 4.0
      WRITE(6,106) FENTDS*1000.0
106   FORMAT(' FENTANYL DOSE CONSTANT ',F6.0,' MG')
      GOTO 200
      
C     NITAZENES (GENERAL ESTIMATE)

20    DNAME = 'NITAZENES'
      HALFLF = 3.5
      CUTOFF = 0.5
      DOSINT = 6.0
      GOTO 200
      
C     AMPHETAMINE

30    DNAME = 'AMPHETAMINE'
      HALFLF = 12.0
      CUTOFF = 25.0
      DOSINT = 12.0
      GOTO 200
      
C     METHAMPHETAMINE  

40    DNAME = 'METHAMPHETAMINE'
      HALFLF = 10.0
      CUTOFF = 25.0
      DOSINT = 8.0
      GOTO 200
      
C     DEXTROAMPHETAMINE

50    DNAME = 'DEXTROAMPHETAMINE'
      HALFLF = 11.0
      CUTOFF = 25.0
      DOSINT = 12.0
      GOTO 200
      
C     HYDROMORPHONE

60    DNAME = 'HYDROMORPHONE'
      HALFLF = 2.3
      CUTOFF = 1.0
      DOSINT = 4.0
      GOTO 200
      
C     OXYCODONE

70    DNAME = 'OXYCODONE'
      HALFLF = 3.2
      CUTOFF = 5.0
      DOSINT = 6.0
      GOTO 200
      
C     MORPHINE

80    DNAME = 'MORPHINE'
      HALFLF = 2.0
      CUTOFF = 10.0
      DOSINT = 4.0
      GOTO 200
      
C     HYDROCODONE

90    DNAME = 'HYDROCODONE'
      HALFLF = 3.8
      CUTOFF = 5.0
      DOSINT = 6.0
      GOTO 200
      
C     CODEINE

91    DNAME = 'CODEINE'
      HALFLF = 2.9
      CUTOFF = 10.0
      DOSINT = 6.0
      GOTO 200
      
C     PETHIDINE/MEPERIDINE

92    DNAME = 'PETHIDINE/MEPERIDINE'
      HALFLF = 3.2
      CUTOFF = 25.0
      DOSINT = 6.0
      GOTO 200
      
C     BARBITURATES (PHENOBARBITAL AS REFERENCE)

93    DNAME = 'BARBITURATES'
      HALFLF = 72.0
      CUTOFF = 50.0
      DOSINT = 24.0
      GOTO 200
      
C     BENZODIAZEPINES (DIAZEPAM AS REFERENCE)

94    DNAME = 'BENZODIAZEPINES'
      HALFLF = 43.0
      CUTOFF = 10.0
      DOSINT = 24.0
      GOTO 200
      
C     ALCOHOL (ETHANOL)
95    DNAME = 'ALCOHOL'
      HALFLF = 1.0
      CUTOFF = 5.0
      DOSINT = 2.0
      GOTO 200
      
C     LSD (LYSERGIC ACID DIETHYLAMIDE)

96    DNAME = 'LSD'
      HALFLF = 3.6
      CUTOFF = 0.5
      DOSINT = 12.0
      GOTO 200
      
C     KETAMINE

97    DNAME = 'KETAMINE'
      HALFLF = 2.5
      CUTOFF = 25.0
      DOSINT = 4.0
      GOTO 200
      
C     MESCALINE

98    DNAME = 'MESCALINE'
      HALFLF = 6.0
      CUTOFF = 25.0
      DOSINT = 12.0
      GOTO 200
      
C     PSILOCYBIN (AS PSILOCIN)

99    DNAME = 'PSILOCYBIN'
      HALFLF = 2.5
      CUTOFF = 1.0
      DOSINT = 8.0
      GOTO 200
      
C     DMT (N,N-DIMETHYLTRYPTAMINE)

100   DNAME = 'DMT'
      HALFLF = 0.25
      CUTOFF = 1.0
      DOSINT = 1.0
      GOTO 200
      
C     GHB (GAMMA-HYDROXYBUTYRIC ACID)

101   DNAME = 'GHB'
      HALFLF = 0.5
      CUTOFF = 5.0
      DOSINT = 2.0
      GOTO 200
      
C     METHAQUALONE

102   DNAME = 'METHAQUALONE'
      HALFLF = 24.0
      CUTOFF = 25.0
      DOSINT = 12.0
      GOTO 200
      
C     METHADONE

103   DNAME = 'METHADONE'
      HALFLF = 22.0
      CUTOFF = 25.0
      DOSINT = 24.0
      GOTO 200
      
C     DEXTROPROPOXYPHENE (PROPOXYPHENE)

104   DNAME = 'DEXTROPROPOXYPHENE'
      HALFLF = 14.0
      CUTOFF = 10.0
      DOSINT = 8.0
      
C     CONTINUE WITH ROUTE-SPECIFIC CALCULATIONS
200   CONTINUE
      
C     SET ROUTE-SPECIFIC PARAMETERS
C     BIOAVAILABILITY AND ABSORPTION RATE CONSTANTS
      GOTO (300,310,320,330,340,350,360,370,380,390,400), ROUTE
      
C     ORAL ADMINISTRATION

300   RNAME = 'ORAL'
      BIOAVAIL = 0.7
      ABSORPT = 1.5
      ORALFAC = 0.01
      GOTO 500
      
C     INTRAVENOUS ADMINISTRATION

310   RNAME = 'INTRAVENOUS'
      BIOAVAIL = 1.0
      ABSORPT = 0.1
      ORALFAC = 0.05
      GOTO 500
      
C     INTRAMUSCULAR ADMINISTRATION

320   RNAME = 'INTRAMUSCULAR'
      BIOAVAIL = 0.9
      ABSORPT = 0.5
      ORALFAC = 0.03
      GOTO 500
      
C     SUBCUTANEOUS ADMINISTRATION

330   RNAME = 'SUBCUTANEOUS'
      BIOAVAIL = 0.8
      ABSORPT = 0.8
      ORALFAC = 0.025
      GOTO 500
      
C     INTRANASAL ADMINISTRATION

340   RNAME = 'INTRANASAL'
      BIOAVAIL = 0.6
      ABSORPT = 0.3
      ORALFAC = 0.02
      GOTO 500
      
C     INHALATION ADMINISTRATION

350   RNAME = 'INHALATION'
      BIOAVAIL = 0.9
      ABSORPT = 0.1
      ORALFAC = 0.04
      GOTO 500
      
C     SUBLINGUAL ADMINISTRATION

360   RNAME = 'SUBLINGUAL'
      BIOAVAIL = 0.8
      ABSORPT = 0.5
      ORALFAC = 0.02
      GOTO 500
      
C     TRANSDERMAL ADMINISTRATION

370   RNAME = 'TRANSDERMAL'
      BIOAVAIL = 0.9
      ABSORPT = 4.0
      ORALFAC = 0.015
      GOTO 500
      
C     RECTAL ADMINISTRATION

380   RNAME = 'RECTAL'
      BIOAVAIL = 0.7
      ABSORPT = 1.0
      ORALFAC = 0.015
      GOTO 500
      
C     BUCCAL ADMINISTRATION

390   RNAME = 'BUCCAL'
      BIOAVAIL = 0.75
      ABSORPT = 0.8
      ORALFAC = 0.025
      GOTO 500
      
C     TOPICAL ADMINISTRATION

400   RNAME = 'TOPICAL'
      BIOAVAIL = 0.1
      ABSORPT = 8.0
      ORALFAC = 0.005
      
500   CONTINUE
      
C     CALCULATE ROUTE-ADJUSTED SINGLE DOSE CONCENTRATION
C     ADJUST FOR DRUG-SPECIFIC ROUTE COMPATIBILITY
      CALL ROUTEADJ(DRUG, ROUTE, BIOAVAIL, ORALFAC, ABSORPT)
      
C     CALCULATE SINGLE DOSE CONCENTRATION
C     FOR FENTANYL, USE FENTDS CONSTANT IF SELECTED
C     FOR ALCOHOL, USE SPECIAL CALCULATION (MG TO ML CONVERSION)
      GOTO (510,520,520,520,520,520,520,520,520,520,520,520,520,
     1      521,520,520,520,520,520,520,520,520,520), DRUG
510   SINGLE = FENTDS * 1000.0 * ORALFAC * BIOAVAIL / REAL(WEIGHT)
      GOTO 530
521   SINGLE = REAL(DOSAGE) * ORALFAC * BIOAVAIL * 0.5 / REAL(WEIGHT)
      GOTO 530
520   SINGLE = REAL(DOSAGE) * ORALFAC * BIOAVAIL / REAL(WEIGHT)
      
530   CONTINUE
      
C     ADJUST HALF-LIFE FOR ABSORPTION RATE (FLIP-FLOP KINETICS)
C     WHEN ABSORPTION IS SLOWER THAN ELIMINATION
     
      IF (ABSORPT.GT.HALFLF*0.693) THEN
        HALFLF = HALFLF * (1.0 + ABSORPT/(HALFLF*0.693))
      
C     ADJUST FOR AGE (SLOWER ELIMINATION WITH AGE)
      AGEFAC = 1.0
      GOTO (540,541,542,543), MIN0(4,MAX0(1,(AGE-20)/15+1))
540   AGEFAC = 1.15
      GOTO 550
541   AGEFAC = 1.0
      GOTO 550
542   AGEFAC = 0.85
      GOTO 550
543   AGEFAC = 0.7
      
550   CONTINUE
      
C     ADJUST FOR METABOLISM RATE
      GOTO (560,570,580), METAB
560   METFAC = 0.7
      GOTO 590
570   METFAC = 1.0
      GOTO 590
580   METFAC = 1.4
      
590   CONTINUE
      
C     CALCULATE EFFECTIVE ELIMINATION RATE
      ELIMRT = 0.693 / HALFLF
      ELIMRT = ELIMRT * AGEFAC * METFAC
      
C     CALCULATE ACCUMULATION FACTOR BASED ON DURATION OF USE
C     ASSUMES REGULAR DOSING AT DOSINT INTERVALS
      NDOSES = INT(USEDUR / DOSINT) + 1
      
C     CALCULATE STEADY-STATE ACCUMULATION FACTOR
C     R = 1 - EXP(-k*tau) WHERE tau = DOSING INTERVAL
      R = 1.0 - EXP(-ELIMRT * DOSINT)
      
C     ACCUMULATION FACTOR FOR N DOSES
C     AF = (1 - R^N) / (1 - R) FOR R != 1
      GOTO (600,601), MIN0(2,MAX0(1,INT(SIGN(1.0,ABS(R-1.0)-0.001))+2))
600   ACCUMF = REAL(NDOSES)
      GOTO 610
601   ACCUMF = (1.0 - R**NDOSES) / (1.0 - R)
      
610   CONTINUE
      
C     CALCULATE TOTAL ACCUMULATED CONCENTRATION
      TOTCONC = SINGLE * ACCUMF
      
C     CALCULATE STEADY-STATE CONCENTRATION (IF REACHED)
      STEADY = SINGLE / (1.0 - EXP(-ELIMRT * DOSINT))
      
C     DETERMINE BUILDUP PERCENTAGE TO STEADY STATE
      BUILDUP = (TOTCONC / STEADY) * 100.0
      GOTO (620,621), MIN0(2,MAX0(1,INT(SIGN(1.0,BUILDUP-100.0))+2))
620   GOTO 630
621   BUILDUP = 100.0
      
630   CONTINUE
      
C     CALCULATE TIME TO REACH CUTOFF CONCENTRATION
C     USING EXPONENTIAL DECAY C(t) = C0 * EXP(-kt)
C     SOLVE FOR t WHEN C(t) = CUTOFF
      DETECT = 0.0
      GOTO (640,641), MIN0(2,MAX0(1,INT(SIGN(1.0,TOTCONC-CUTOFF))+2))
640   GOTO 650
641   DETECT = ALOG(TOTCONC / CUTOFF) / ELIMRT
      
650   CONTINUE
      
C     CONVERT TO SECONDS
      DETECT = DETECT * 3600.0
      
C     CALL GRAPHING SUBROUTINE
      CALL PLOTXY(TOTCONC, ELIMRT, CUTOFF, HALFLF, USEDUR)
      
C     OUTPUT RESULTS
      WRITE(6,700) DNAME
700   FORMAT(//,' DETECTION TIME CALCULATION FOR ',A20,/)
      
      WRITE(6,701) DOSAGE, WEIGHT, AGE
701   FORMAT(' INPUT PARAMETERS',/,
     1       '   DOSAGE ',I4,' MG',/,
     2       '   WEIGHT ',I3,' KG',/,
     3       '   AGE ',I3,' YEARS')
      
      GOTO (710,720,730), METAB
710   WRITE(6,702) 'SLOW'
      GOTO 740
720   WRITE(6,702) 'NORMAL'
      GOTO 740
730   WRITE(6,702) 'FAST'
702   FORMAT('   METABOLISM ',A6)
      
740   CONTINUE
      
      WRITE(6,708) USEDUR, USEDUR/24.0
708   FORMAT('   DURATION OF USE ',F6.1,' HOURS (',F5.2,' DAYS)')
      
      WRITE(6,709) RNAME, BIOAVAIL*100.0, ABSORPT
709   FORMAT('   ROUTE ',A12,' (BIOAVAIL ',F5.1,'%,',
     1       ' ABS RATE ',F5.2,' HR)')
      
      GOTO (750,760), MIN0(2,MAX0(1,DRUG))
750   WRITE(6,707) FENTDS*1000.0
707   FORMAT('   FENTANYL DOSE ',F6.0,' MG (CONSTANT)')
760   CONTINUE
      
      WRITE(6,703) HALFLF, CUTOFF, DOSINT
703   FORMAT(/,' PHARMACOKINETIC DATA',/,
     1       '   HALF-LIFE ',F5.1,' HOURS',/,
     2       '   CUTOFF ',F6.1,' NG/ML',/,
     3       '   DOSING INTERVAL ',F5.1,' HOURS')
      
      WRITE(6,710) NDOSES, ACCUMF
710   FORMAT('   NUMBER OF DOSES ',I3,/,
     1       '   ACCUMULATION FACTOR ',F6.2)
      
      WRITE(6,704) SINGLE, TOTCONC, ELIMRT*3600.0
704   FORMAT('   SINGLE DOSE CONC ',F8.2,' NG/ML',/,
     1       '   TOTAL ACCUM CONC ',F8.2,' NG/ML',/,
     2       '   ELIM RATE ',F8.4,' /HOUR')
      
      WRITE(6,711) STEADY, BUILDUP
711   FORMAT('   STEADY-STATE CONC ',F8.2,' NG/ML',/,
     1       '   BUILDUP TO SS ',F5.1,' %')
      
C     CONVERT SECONDS TO READABLE FORMAT
      ISEC = INT(DETECT)
      IHOUR = ISEC / 3600
      IMIN = (ISEC - IHOUR*3600) / 60
      ISEC = ISEC - IHOUR*3600 - IMIN*60
      
      WRITE(6,705) DETECT
705   FORMAT(/,' DETECTION TIME ',F10.0,' SECONDS')
      
      WRITE(6,706) IHOUR, IMIN, ISEC
706   FORMAT(' EQUIVALENT TO ',I3,' HOURS, ',I2,' MINUTES, ',
     1       I2,' SECONDS')
      
C     ADDITIONAL TIME CALCULATIONS
      IDAYS = IHOUR / 24
      IHOUR = IHOUR - IDAYS*24
      WRITE(6,712) IDAYS, IHOUR, IMIN, ISEC
712   FORMAT(' FULL FORMAT ',I2,' DAYS, ',I2,' HOURS, ',
     1       I2,' MINUTES, ',I2,' SECONDS')
      
C     WARNINGS AND DISCLAIMERS
      WRITE(6,800)
800   FORMAT(//,' ** IMPORTANT DISCLAIMERS **',/,
     1       ' - ESTIMATES BASED ON POPULATION AVERAGES',/,
     2       ' - INDIVIDUAL VARIATION CAN BE SIGNIFICANT',/,
     3       ' - CHRONIC USE CALCULATIONS ARE SIMPLIFIED',/,
     4       ' - ASSUMES REGULAR DOSING INTERVALS',/,
     5       ' - ROUTE-SPECIFIC PARAMETERS ARE ESTIMATES',/,
     6       ' - DRUG-ROUTE COMBINATIONS MAY VARY',/,
     7       ' - FOR RESEARCH/EDUCATIONAL USE ONLY',/,
     8       ' - CONSULT TOXICOLOGY REFERENCES FOR',/,
     9       '   CLINICAL APPLICATIONS')
      
      STOP
      END
      
C     SUBROUTINE TO ADJUST ROUTE PARAMETERS FOR DRUG COMPATIBILITY
      SUBROUTINE ROUTEADJ(DRUG, ROUTE, BIOAVAIL, ORALFAC, ABSORPT)
      IMPLICIT REAL (A-H,O-Z)
      INTEGER DRUG, ROUTE
      REAL BIOAVAIL, ORALFAC, ABSORPT
      
C     ADJUST PARAMETERS BASED ON DRUG-ROUTE COMPATIBILITY
C     SOME DRUGS ARE NOT SUITABLE FOR CERTAIN ROUTES
      
C     ALCOHOL SPECIFIC ADJUSTMENTS

      IF (DRUG.EQ.14) THEN
        IF (ROUTE.EQ.2.OR.ROUTE.EQ.3.OR.ROUTE.EQ.4) THEN
C         IV/IM/SC NOT TYPICAL FOR ALCOHOL
          BIOAVAIL = BIOAVAIL * 0.1
        
        IF (ROUTE.EQ.6) THEN

C         INHALATION (VAPOR) HIGHLY BIOAVAILABLE
          BIOAVAIL = 0.95
          ABSORPT = 0.05
      
C     FENTANYL SPECIFIC ADJUSTMENTS
     
      IF (DRUG.EQ.1) THEN
        IF (ROUTE.EQ.8) THEN

C         TRANSDERMAL PATCH - SUSTAINED RELEASE
         
          ABSORPT = 12.0
          BIOAVAIL = 0.92
        
        IF (ROUTE.EQ.7) THEN

C         SUBLINGUAL - HIGH BIOAVAILABILITY
        
          BIOAVAIL = 0.8

C     STIMULANTS (AMPHETAMINES) ADJUSTMENTS
     
      IF (DRUG.GE.3.AND.DRUG.LE.5) THEN
        IF (ROUTE.EQ.5) THEN

C         INTRANASAL - COMMON ROUTE
  
          BIOAVAIL = 0.8
          ABSORPT = 0.2
        
        IF (ROUTE.EQ.6) THEN

C         SMOKING/INHALATION
      
          BIOAVAIL = 0.7
          ABSORPT = 0.08
        
C     OPIOIDS GENERAL ADJUSTMENTS

      IF ((DRUG.GE.6.AND.DRUG.LE.11).OR.DRUG.EQ.22) THEN
        IF (ROUTE.EQ.2) THEN

C         IV - FULL BIOAVAILABILITY
     
          BIOAVAIL = 1.0
          ORALFAC = 0.08
        
        IF (ROUTE.EQ.5) THEN

C         INTRANASAL - MODERATE BIOAVAILABILITY
          BIOAVAIL = 0.65
        
C     PSYCHEDELICS ADJUSTMENTS
 
      IF (DRUG.GE.15.AND.DRUG.LE.19) THEN
        IF (ROUTE.EQ.6) THEN

C         INHALATION NOT TYPICAL FOR MOST PSYCHEDELICS
         
          BIOAVAIL = BIOAVAIL * 0.3
        
        IF (DRUG.EQ.19) THEN

C         DMT - TYPICALLY SMOKED/VAPORIZED
          
          IF (ROUTE.EQ.6) THEN
  
            BIOAVAIL = 0.8
            ABSORPT = 0.02
          

C     BENZODIAZEPINES ADJUSTMENTS
    
      IF (DRUG.EQ.13) THEN
        IF (ROUTE.EQ.7) THEN

C         SUBLINGUAL - GOOD ABSORPTION
         
          BIOAVAIL = 0.9
          ABSORPT = 0.3
        
        IF (ROUTE.EQ.9) THEN

C         RECTAL - GOOD ABSORPTION FOR SEIZURE CONTROL
          
          BIOAVAIL = 0.8
          ABSORPT = 0.5
        
C     KETAMINE ADJUSTMENTS
   
      IF (DRUG.EQ.16) THEN
        IF (ROUTE.EQ.5) THEN

C         INTRANASAL - COMMON ROUTE
      
          BIOAVAIL = 0.5
          ABSORPT = 0.3
        
        IF (ROUTE.EQ.3) THEN

C         INTRAMUSCULAR - CLINICAL USE
         
          BIOAVAIL = 0.93
          ABSORPT = 0.3
        
C     GHB ADJUSTMENTS
    
      IF (DRUG.EQ.20) THEN
        IF (ROUTE.NE.1) THEN

C         GHB PRIMARILY ORAL
    
          BIOAVAIL = BIOAVAIL * 0.5
        
C     TOPICAL ROUTE RESTRICTIONS
   
      IF (ROUTE.EQ.11) THEN

C         TOPICAL GENERALLY LOW SYSTEMIC ABSORPTION
       
        IF (DRUG.NE.1.AND.DRUG.NE.22) THEN
C         EXCEPT FOR FENTANYL AND METHADONE PATCHES
          
          BIOAVAIL = 0.05
          ORALFAC = 0.002
        
      RETURN
      END
      
C     SUBROUTINE PLOTXY - ENHANCED FOR ACCUMULATION DISPLAY
      SUBROUTINE PLOTXY(C0, KELIM, CUTOFF, THALF, DURATION)
      IMPLICIT REAL (A-H,O-Z)
      DIMENSION CONC(61), TIME(61), A(119)
      CHARACTER*1 CHAR(10), BLANK, PLUS
      DATA BLANK,PLUS,CHAR(1),CHAR(2),CHAR(3),CHAR(4),CHAR(5),
     1CHAR(6),CHAR(7),CHAR(8),CHAR(9),CHAR(10)
     1/1H ,1H+,1HA,1HB,1HC,1HD,1HE,1HF,1HG,1HH,1HI,1HJ/
      
C     CALCULATE TIME POINTS (0 TO 8 HALF-LIVES OR DURATION*2)
      
      TMAX = AMAX1(8.0 * THALF, DURATION * 2.0)
      DT = TMAX / 60.0
      
      WRITE(6,600)
600   FORMAT(//,
     1' ====================================================',/,
     2'   PLASMA CONCENTRATION vs TIME WITH ACCUMULATION   ',/,
     3'        (INCLUDES CHRONIC USE BUILD-UP EFFECTS)     ',/,
     4'        (ADJUSTED FOR ROUTE OF ADMINISTRATION)      ',/,
     5' ====================================================',//)
      
C     CALCULATE CONCENTRATION AT EACH TIME POINT
     
      DO 10 I = 1, 61
        TIME(I) = REAL(I-1) * DT
        CONC(I) = C0 * EXP(-KELIM * TIME(I))
10    CONTINUE
      
C     FIND SCALING FACTORS
      
      CMAX = C0
      
C     WRITE HEADER WITH SCALE INFO
     
      WRITE(6,601) 50, 118, 200E-03
601   FORMAT(' ',31X,'NO. OF POINTS ON ABSCISSA = ',I5/32X,
     1       'NO. OF INTERVALS ON ORDINATE = ',I3,
     2       ', INCREMENT = ',E10.3)
      
      WRITE(6,602) CMAX, CMAX*0.0
602   FORMAT(' ',E16.6,E101.6/10X,'(',99X,'1'/' 1',8X,'1',
     1       2(49X,'['),8X,'I')
      
C     INITIALIZE PLOT ARRAY
    
      DO 15 I = 1, 119
15      A(I) = BLANK
      
C     PLOT CONCENTRATION DECAY CURVE
      
      DO 40 I = 1, 61

C       CLEAR ARRAY WITH BLANKS
        
        DO 20 J = 1, 119
20        A(J) = BLANK
        
C       ADD GRID MARKERS EVERY 10 POSITIONS
       
        IF (MOD(I,6).EQ.1) THEN
          DO 25 J = 10, 119, 10
25          A(J) = PLUS
      
C       PLOT CONCENTRATION POINT
        
        IPOS = INT(CONC(I) * 118.0 / CMAX) + 1
        IPOS = MIN0(119, MAX0(1, IPOS))
        A(IPOS) = CHAR(1)
        
C       PLOT CUTOFF LINE IF VISIBLE
       
        ICUT = INT(CUTOFF * 118.0 / CMAX) + 1
        IF (ICUT.GE.1.AND.ICUT.LE.119) THEN
          A(ICUT) = CHAR(3)

C       OUTPUT LINE
      
        WRITE(6,603) (A(J), J=1,119)
603     FORMAT(' ',119A1)
        
C       CLEAR GRID MARKERS IF THEY WERE SET
       
        IF (MOD(I,6).EQ.1) THEN
          DO 30 J = 10, 119, 10
30          A(J) = BLANK
40    CONTINUE
      
C     CALCULATE FINAL STATISTICS
     
      TFINAL = 0.0
      IF (C0.GT.CUTOFF) THEN
        TFINAL = ALOG(C0/CUTOFF) / KELIM
        
      WRITE(6,605) TFINAL, C0, CUTOFF, DURATION
605   FORMAT(/,' LEGEND A = CONCENTRATION CURVE',/,
     1       '         C = DETECTION CUTOFF THRESHOLD',/,
     2       '         + = TIME GRID MARKERS',//,
     3       ' ANALYSIS TIME TO NON-DETECTION = ',F6.2,' HOURS',/,
     4       '          ACCUMULATED PEAK CONC = ',F8.2,' NG/ML',/,
     5       '          CUTOFF LEVEL = ',F8.2,' NG/ML',/,
     6       '          DURATION OF USE = ',F6.2,' HOURS',//)
      
      RETURN
      END
      
C     SUBROUTINE TO CONVERT STRING TO UPPERCASE
      SUBROUTINE UPPER(STRING)
      CHARACTER*20 STRING
      INTEGER I, ICHAR
      
      DO 10 I = 1, 20
        ICHAR = ICHAR(STRING(I:I))
        IF (ICHAR.GE.97.AND.ICHAR.LE.122) THEN
          STRING(I:I) = CHAR(ICHAR - 32)
        
10    CONTINUE
      
      RETURN
      END
      
C     SUBROUTINE FOR ROUTE-SPECIFIC ORAL FLUID TRANSFER CALCULATIONS
      SUBROUTINE ORFLUID(DRUG, ROUTE, PLASMA, SALIVA, TRANSFER)
      IMPLICIT REAL (A-H,O-Z)
      INTEGER DRUG, ROUTE
      REAL PLASMA, SALIVA, TRANSFER
      
C     CALCULATE SALIVA/PLASMA RATIO BASED ON DRUG PROPERTIES
C     AND ROUTE OF ADMINISTRATION
      
C     BASE TRANSFER RATIOS (SALIVA/PLASMA)
      
      GOTO (100,200,300,400,500,600,700,800,900,910,920,
     1      930,940,950,960,970,980,990,1000,1010,1020,
     2      1030,1040), DRUG
      
C     FENTANYL

100   TRANSFER = 0.8
      GOTO 2000
      
C     NITAZENES

200   TRANSFER = 0.6
      GOTO 2000
      
C     AMPHETAMINE

300   TRANSFER = 2.8
      GOTO 2000
      
C     METHAMPHETAMINE

400   TRANSFER = 4.2
      GOTO 2000
      
C     DEXTROAMPHETAMINE

500   TRANSFER = 2.9
      GOTO 2000
      
C     HYDROMORPHONE

600   TRANSFER = 0.7
      GOTO 2000
      
C     OXYCODONE

700   TRANSFER = 1.1
      GOTO 2000
      
C     MORPHINE

800   TRANSFER = 1.2
      GOTO 2000
      
C     HYDROCODONE

900   TRANSFER = 1.0
      GOTO 2000
      
C     CODEINE

910   TRANSFER = 1.4
      GOTO 2000
      
C     PETHIDINE

920   TRANSFER = 3.5
      GOTO 2000
      
C     BARBITURATES

930   TRANSFER = 0.3
      GOTO 2000
      
C     BENZODIAZEPINES

940   TRANSFER = 0.2
      GOTO 2000
      
C     ALCOHOL

950   TRANSFER = 1.0
      GOTO 2000
      
C     LSD

960   TRANSFER = 0.1
      GOTO 2000
      
C     KETAMINE

970   TRANSFER = 0.8
      GOTO 2000
      
C     MESCALINE

980   TRANSFER = 1.5
      GOTO 2000
      
C     PSILOCYBIN

990   TRANSFER = 0.3
      GOTO 2000
      
C     DMT

1000  TRANSFER = 0.2
      GOTO 2000
      
C     GHB

1010  TRANSFER = 0.8
      GOTO 2000
      
C     METHAQUALONE

1020  TRANSFER = 0.6
      GOTO 2000
      
C     METHADONE

1030  TRANSFER = 0.9
      GOTO 2000
      
C     DEXTROPROPOXYPHENE

1040  TRANSFER = 1.8
      
2000  CONTINUE
      
C     ADJUST TRANSFER RATIO BASED ON ROUTE
C     SOME ROUTES AFFECT SALIVA CONTAMINATION
     
      IF (ROUTE.EQ.1) THEN

C       ORAL - POTENTIAL CONTAMINATION
       
        TRANSFER = TRANSFER * 1.5
      IF (ROUTE.EQ.7) THEN

C       SUBLINGUAL - HIGH CONTAMINATION
        
        TRANSFER = TRANSFER * 3.0
      
      IF (ROUTE.EQ.10) THEN

C       BUCCAL - HIGH CONTAMINATION
        
        TRANSFER = TRANSFER * 2.5
      IF (ROUTE.EQ.5) THEN

C       INTRANASAL - SOME CONTAMINATION VIA POSTNASAL DRIP
       
        TRANSFER = TRANSFER * 1.2
      IF (ROUTE.EQ.6) THEN

C       INHALATION - CONTAMINATION VIA ORAL CAVITY
        
        TRANSFER = TRANSFER * 1.8

C     CALCULATE SALIVA CONCENTRATION
      
      SALIVA = PLASMA * TRANSFER
      
      RETURN
      END
