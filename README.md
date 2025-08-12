# NARCDETECT – Drug Detection Time Calculator v2.0

**For Oral Fluid (Saliva) Testing**  
**Includes NMR Spectrum Simulation**

---

## Overview

**NARCDETECT v2.0** is a pharmacokinetics-based calculator for estimating the time until a given drug becomes non-detectable in oral fluid (saliva) testing.  
It supports a wide range of drug classes and accounts for:

- **Pharmacokinetic parameters**  
- **Chronic use accumulation effects**  
- **Route-specific bioavailability and absorption rates**  
- **Metabolism speed variations**  
- **NMR spectrum simulation for forensic identification**

---

## Supported Drug Classes

**Synthetic Opioids:**
- Fentanyl
- Nitazenes
- Hydromorphone
- Oxycodone
- Hydrocodone
- Dextropropoxyphene

**Natural Opioids:**
- Morphine
- Codeine
- Diamorphine (Heroin)
- Pethidine
- Methadone

**Stimulants:**
- Amphetamine
- Methamphetamine
- Dextroamphetamine

**Depressants:**
- Barbiturates
- Benzodiazepines
- Methaqualone

**Psychedelics:**
- LSD
- Mescaline
- Psilocybin
- DMT

**Other:**
- Alcohol
- Ketamine
- GHB

---

## Input Parameters

The program prompts for:

1. **Drug name** (from available list)  
2. **Route of administration**  
   - Oral, IV, IM, SC, Intranasal, Inhalation, Sublingual, Transdermal, Rectal, Buccal, Topical
3. **Dosage** (mg)  
4. **Body weight** (kg)  
5. **Age** (years)  
6. **Metabolism rate** (1 = Slow, 2 = Normal, 3 = Fast)  
7. **Duration of use** (hours)  

---

## Output

### 1. Pharmacokinetic Analysis
- **Plasma concentration vs time** chart  
  - Includes accumulation from chronic use  
  - Adjusted for route-specific parameters  
  - Detection cutoff line shown
- **Key statistics:**
  - Time to non-detection
  - Peak concentration
  - Duration of use
  - Dosing interval
  - Accumulation factor
  - Steady-state concentration

### 2. Detection Time
- Presented in **seconds**, **hours/minutes**, and **full date-time format**

### 3. Disclaimers
- Estimates are based on population averages  
- Individual variation is significant  
- Chronic use model simplified  
- Not for medical or legal decision-making

---

## NMR Spectrum Simulation

When enabled, the program generates a synthetic **¹H NMR spectrum** for the chosen drug:

- **Chemical shift range:** 0.0 – 12.0 ppm  
- **Graphical spectrum output** with peak markers  
- **Peak assignments** with shift, intensity, width, and likely chemical group

Example peak table:

| Shift (ppm) | Intensity | Width  | Assignment         |
|-------------|-----------|--------|--------------------|
| 6.80        | 50.0      | 0.09   | O–CH, N–CH         |
| 6.50        | 50.0      | 0.09   | O–CH, N–CH         |
| 4.20        | 60.0      | 0.08   | O–CH, N–CH         |
| 3.00        | 90.0      | 0.09   | CH₂, CH₃ Alpha     |
| 2.10        | 100.0     | 0.10   | CH₂, CH₃ Alpha     |

---

## Example Run

====================================================================
NARCDETECT - Drug Detection Time Calculator v2.0
FOR ORAL FLUID (SALIVA) TESTING
====================================================================

ESTIMATES TIME UNTIL NON-DETECTABLE
BASED ON PHARMACOKINETIC PARAMETERS
INCLUDES CHRONIC USE ACCUMULATION
AND ROUTES OF ADMINISTRATION
WITH NMR SPECTRUM SIMULATION

AVAILABLE DRUGS BY TYPE:
====================================================================
SYNTHETIC OPIOIDS:     NATURAL OPIOIDS:       STIMULANTS:

  FENTANYL               MORPHINE               AMPHETAMINE
  NITAZENES              CODEINE                METHAMPHETAMINE
  HYDROMORPHONE          DIAMORPHINE (HEROIN)   DEXTROAMPHETAMINE
  OXYCODONE              PETHIDINE
  HYDROCODONE            METHADONE
  DEXTROPROPOXYPHENE

DEPRESSANTS:           PSYCHEDELICS:          OTHER:

  BARBITURATES           LSD                    ALCOHOL
  BENZODIAZEPINES        MESCALINE              KETAMINE
  METHAQUALONE           PSILOCYBIN             GHB
                         DMT
====================================================================

Enter drug name: HEROIN
AVAILABLE ROUTES OF ADMINISTRATION:
ORAL, INTRAVENOUS, INTRAMUSCULAR,
SUBCUTANEOUS, INTRANASAL, INHALATION,
SUBLINGUAL, TRANSDERMAL, RECTAL,
BUCCAL, TOPICAL

Enter route of administration: IV
Enter dosage in mg: 1000 
Enter body weight in kg: 76
Enter age in years: 28
Metabolism rate (1=SLOW, 2=NORMAL, 3=FAST): 3
Duration of use in hours (24.0=1 day): 48 
====================================================================
  PLASMA CONCENTRATION vs TIME WITH ACCUMULATION
       (INCLUDES CHRONIC USE BUILD-UP EFFECTS)
       (ADJUSTED FOR ROUTE OF ADMINISTRATION)
====================================================================

Time range: 0 to 96.0 hours
Maximum concentration: 13.68 ng/mL
Cutoff level: 10.00 ng/mL

         +         +         +         +         +         +         +         +     C   +         +         +       A 
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A        +         +         +         +         +         +         +         +     C   +         +         +         
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A        +         +         +         +         +         +         +         +     C   +         +         +         
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A        +         +         +         +         +         +         +         +     C   +         +         +         
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A        +         +         +         +         +         +         +         +     C   +         +         +         
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A        +         +         +         +         +         +         +         +     C   +         +         +         
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A        +         +         +         +         +         +         +         +     C   +         +         +         
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A        +         +         +         +         +         +         +         +     C   +         +         +         
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A        +         +         +         +         +         +         +         +     C   +         +         +         
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A        +         +         +         +         +         +         +         +     C   +         +         +         
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A                                                                                    C                                 
A        +         +         +         +         +         +         +         +     C   +         +         +         

LEGEND: A = CONCENTRATION CURVE
        C = DETECTION CUTOFF THRESHOLD
        + = TIME GRID MARKERS

ANALYSIS: Time to non-detection = 0.05 hours
          Peak concentration = 13.68 ng/mL
          Duration of use = 48.00 hours


====================================================================
DETECTION TIME CALCULATION FOR DIAMORPHINE
====================================================================

INPUT PARAMETERS:
  Dosage: 1000 mg
  Weight: 76 kg
  Age: 28 years
  Metabolism: FAST
  Duration of use: 48.0 hours (2.00 days)
  Route: INTRAVENOUS (Bioavail 100.0%, Abs rate 0.10 hr)

PHARMACOKINETIC DATA:
  Half-life: 0.2 hours
  Cutoff: 10.0 ng/mL
  Dosing interval: 4.0 hours
  Number of doses: 13
  Accumulation factor: 13.00
  Single dose conc: 1.05 ng/mL
  Total accum conc: 13.68 ng/mL
  Elim rate: 5.7423 /hour
  Steady-state conc: 1.05 ng/mL
  Buildup to SS: 100.0%

DETECTION TIME: 197 seconds
EQUIVALENT TO: 0 hours, 3 minutes, 16 seconds
FULL FORMAT: 0 days, 0 hours, 3 minutes, 16 seconds

** IMPORTANT DISCLAIMERS **
- Estimates based on population averages
- Individual variation can be significant
- Chronic use calculations are simplified
- Assumes regular dosing intervals
- Route-specific parameters are estimates
- For research/educational use only

Generate NMR spectrum simulation? (Y/N): Y
====================================================================
          1H NMR SPECTRUM SIMULATION FOR DIAMORPHINE
       CONCENTRATION: 1000.00 NG/ML IN SAMPLE
       CHEMICAL SHIFT RANGE: 0.0 - 12.0 PPM
       SYNTHETIC SPECTRUM FOR IDENTIFICATION
====================================================================

Maximum intensity = 1010.20 (relative)
Chemical shift scale: 12.0 to 0.0 PPM

      12        10         8         6         4         2         0
       |         |         |         |         |         |         |
+                   +                   +                   +                   +                  *+                   +
+                   +                   +                   +                   +                  *+                   +
+                   +                   +                   +                   +                  *+                   +
+                   +                   +                   +                   +                  *+                   +
+                   +                   +                   +                   +                  *+                   +
+                   +                   +                   +                   +         *        *+                   +
+                   +                   +                   +                   +         *        *+                   +
+                   +                   +                   +                   +         *        *+                   +
+                   +                   +                   +                   +         *        *+                   +
+                   +                   +                   +                   +         *        *+                   +
+                   +                   +                   +                   +         *        *+                   +
+                   +                   +                   +                   +         *        *+                   +
+                   +                   +                   +                   +         *        *+                   +
+                   +                   +                   +                   +         *        *+                   +
+                   +                   +                   +                   +         *        *+                   +
+                   +                   +                   +                   +         *        *+                   +
+                   +                   +                   +                   +         *        *+                   +
+                   +                   +                   +                   +         *        *+                   +
+                   +                   +                   +                   +         *        *+                   +
+                   +                   +                   +                   +         *        *+                   +
+                   +                   +                   +                 * +         *        *+                   +
+                   +                   +                   +                 * +         *        *+                   +
+                   +                   +                   +                 * +         *        *+                   +
+                   +                   +                   +                 * +         *        *+                   +
+                   +                   +           *  *    +                 * +         *        *+                   +
+                   +                   +           *  *    +                 * +         *        *+                   +
+                   +                   +           *  *    +                 * +         *       ***                   +
+                   +                   +           *  *    +                 * +         *       ***                   +
+                   +                   +           *  *    +                 * +         *       ***                   +
+                   +                   +           *  *    +                 * +         *       ***                   +
+                   +                   +           *  *    +                 * +        ***      ***                   +
+                   +                   +           *  *    +                 * +        ***      ***                   +
+                   +                   +           *  *    +                 * +        ***      ***                   +
+                   +                   +           *  *    +                 * +        ***      ***                   +
+                   +                   +           *  *    +                 * +        ***      ***                   +
+                   +                   +           * **    +                 * +        ***      ***                   +
+                   +                   +           ****    +                 * +        ***      ***                   +
+                   +                   +           ****    +                 * +        ***      ***                   +
+                   +                   +           *****   +                ***+        ***      ***                   +
+                   +                   +          ******   +                ***+        ***      ***                   +
+                   +                   +          ******   +                ***+        ***     ****                   +
+                   +                   +          ******   +                ***+        ***     *****                  +
+                   +                   +          ******   +                ***+       *****    *****                  +
+                   +                   +          ******   +                ***+       *****    *****                  +
+                   +                   +          ******   +                ***+       *****    *****                  +
+                   +                   +          ******   +                ***+       ******  ******                  +
+                   +                   +         ********  +               *****      ****************                 +
+                   +                   +         ********  +               *****      ****************                 +
+                   +                   +        ********** +              *******   *******************                +
+                   +                   +      **************            **********************************             +

PEAK ASSIGNMENTS:
SHIFT(PPM)  INTENSITY  WIDTH   ASSIGNMENT
----------  ---------  -----   ----------
    6.80       50.0     0.09   O-CH, N-CH
    6.50       50.0     0.09   O-CH, N-CH
    4.20       60.0     0.08   O-CH, N-CH
    3.00       90.0     0.09   CH2, CH3 ALPHA
    2.10      100.0     0.10   CH2, CH3 ALPHA

SPECTRUM ANALYSIS:
NUMBER OF PEAKS DETECTED: 5
MAXIMUM PEAK INTENSITY:   1010.20
SAMPLE CONCENTRATION:     1000.00 NG/ML
INTEGRATION COMPLETE

* = SPECTRAL PEAK    + = PPM GRID LINES

----------------------------------------------------------------------------------------------------------------------------

**Example Output Highlights:**
- Peak concentration: 13.68 ng/mL  
- Detection cutoff: 10.00 ng/mL  
- Time to non-detection: **197 seconds** (~3 minutes, 16 seconds)  
- NMR simulation with 5 peaks detected

---

## Author Information

**Name:** Mickey W. Lawless  
**System:** KAYPRO XL  
**Host CPU:** NEC V20  
**Host OS:** MS-DOS v6.22  
**Compiler:** Turbo C Version 1.0  
**Compiler Copyright:** © 1987 Borland International Inc.

**Source Files:**
- `NARCDTC.C` – 1,041 lines
- `NARCII.C` – 1,116 lines

---

## Development Notes

This program was written and compiled entirely on original vintage hardware:

- **Machine:** Kaypro XL — a classic CP/M- and MS-DOS-capable system from the mid-1980s.  
- **Processor:** NEC V20 — binary-compatible with the Intel 8088, but with additional instruction set support and improved performance.  
- **Operating System:** MS-DOS v6.22 — a late-generation DOS release, running on legacy hardware for development authenticity.  
- **Compiler:** Turbo C 1.0 — Borland’s first C compiler for DOS, prized for its compact size and fast compile times.

Source code was written with constraints in mind, including:
- Minimal reliance on non-standard libraries.
- Compatibility with small memory models.
- Optimized for execution speed on sub-10 MHz CPUs.
- Fully self-contained build process with no external dependencies.

This retrocomputing approach preserves the feel of 1980s/early 1990s software development while producing a functional modern-use calculation tool.

---

## Disclaimer

This tool is for **educational and research purposes only**.  
It is **not** intended for clinical, diagnostic, or legal use.



