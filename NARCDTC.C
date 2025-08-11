/*
 * NARCDETECT - Drug Detection Time Calculator for Oral Fluid Tests
 * Turbo C Version - Compatible with Turbo C 1.0/2.0
 * 
 * By: Mickey W. Lawless
 * Date: August 11, 2025
 *
 * Calculates detection window for various drugs including:
 * - Accumulation effects from chronic use
 * - Route of administration variables
 * - NMR spectrum simulation for identification
 * 
 * Based on pharmacokinetic parameters and oral fluid testing
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/* Maximum constants */
#define MAX_PEAKS 20
#define MAX_DRUG_NAME 25
#define MAX_ROUTE_NAME 15
#define SPECTRUM_WIDTH 121
#define PLOT_HEIGHT 50
#define PLOT_WIDTH 119

/* Drug constants */
#define NUM_DRUGS 23
#define NUM_ROUTES 11

/* Drug types */
enum {
    DRUG_FENTANYL = 1,
    DRUG_NITAZENES = 2,
    DRUG_AMPHETAMINE = 3,
    DRUG_METHAMPHETAMINE = 4,
    DRUG_DEXTROAMPHETAMINE = 5,
    DRUG_HYDROMORPHONE = 6,
    DRUG_OXYCODONE = 7,
    DRUG_MORPHINE = 8,
    DRUG_HYDROCODONE = 9,
    DRUG_CODEINE = 10,
    DRUG_PETHIDINE = 11,
    DRUG_BARBITURATES = 12,
    DRUG_BENZODIAZEPINES = 13,
    DRUG_ALCOHOL = 14,
    DRUG_LSD = 15,
    DRUG_KETAMINE = 16,
    DRUG_MESCALINE = 17,
    DRUG_PSILOCYBIN = 18,
    DRUG_DMT = 19,
    DRUG_GHB = 20,
    DRUG_METHAQUALONE = 21,
    DRUG_METHADONE = 22,
    DRUG_DEXTROPROPOXYPHENE = 23
};

/* Route types */
enum {
    ROUTE_ORAL = 1,
    ROUTE_INTRAVENOUS = 2,
    ROUTE_INTRAMUSCULAR = 3,
    ROUTE_SUBCUTANEOUS = 4,
    ROUTE_INTRANASAL = 5,
    ROUTE_INHALATION = 6,
    ROUTE_SUBLINGUAL = 7,
    ROUTE_TRANSDERMAL = 8,
    ROUTE_RECTAL = 9,
    ROUTE_BUCCAL = 10,
    ROUTE_TOPICAL = 11
};

/* Structure definitions */
typedef struct {
    char name[MAX_DRUG_NAME];
    float halflife;
    float cutoff;
    float dosing_interval;
} DrugData;

typedef struct {
    char name[MAX_ROUTE_NAME];
    float bioavailability;
    float absorption_rate;
    float oral_factor;
} RouteData;

typedef struct {
    float shifts[MAX_PEAKS];
    float intensities[MAX_PEAKS];
    float widths[MAX_PEAKS];
    int num_peaks;
} NMRData;

/* Global variables */
static DrugData drugs[NUM_DRUGS + 1];
static RouteData routes[NUM_ROUTES + 1];
static float fentanyl_dose_constant = 1.0f;

/* Function prototypes */
void initialize_drug_data(void);
void initialize_route_data(void);
void print_banner(void);
void print_drug_menu(void);
void print_route_menu(void);
int get_drug_selection(void);
int get_route_selection(void);
void get_input_parameters(int *dosage, int *weight, int *age, int *metab, float *duration);
void adjust_route_parameters(int drug, int route, float *bioavail, float *oral_fac, float *absorpt);
void calculate_detection_time(int drug, int route, int dosage, int weight, int age, 
                            int metab, float duration);
void plot_concentration_curve(float c0, float kelim, float cutoff, float thalf, float duration);
void nmr_plot(int drug, float concentration, NMRData *nmr_data);
void get_peak_label(int drug, int peak_no, float shift, char *label);
void generate_nmr_data(int drug, NMRData *nmr_data);
void str_upper(char *str);
int str_compare_upper(const char *str1, const char *str2);
float max_float(float a, float b);
float min_float(float a, float b);
int max_int(int a, int b);
int min_int(int a, int b);

/* Main program */
int main(void)
{
    int drug, route;
    int dosage, weight, age, metab;
    float duration;
    char answer;

    /* Initialize data tables */
    initialize_drug_data();
    initialize_route_data();

    /* Print program banner */
    print_banner();

    /* Get drug selection */
    drug = get_drug_selection();
    if (drug == 0) {
        printf("Invalid drug selection. Exiting.\n");
        return 1;
    }

    /* Get route selection */
    route = get_route_selection();
    if (route == 0) {
        printf("Invalid route selection. Exiting.\n");
        return 1;
    }

    /* Get input parameters */
    get_input_parameters(&dosage, &weight, &age, &metab, &duration);

    /* Calculate and display results */
    calculate_detection_time(drug, route, dosage, weight, age, metab, duration);

    /* Ask if user wants NMR spectrum */
    printf("\nGenerate NMR spectrum simulation? (Y/N): ");
    scanf(" %c", &answer);
    if (toupper(answer) == 'Y') {
        NMRData nmr_data;
        generate_nmr_data(drug, &nmr_data);
        nmr_plot(drug, (float)dosage, &nmr_data);
    }

    return 0;
}

void initialize_drug_data(void)
{
    /* Initialize drug database */
    strcpy(drugs[DRUG_FENTANYL].name, "FENTANYL");
    drugs[DRUG_FENTANYL].halflife = 2.0f;
    drugs[DRUG_FENTANYL].cutoff = 1.0f;
    drugs[DRUG_FENTANYL].dosing_interval = 4.0f;

    strcpy(drugs[DRUG_NITAZENES].name, "NITAZENES");
    drugs[DRUG_NITAZENES].halflife = 3.5f;
    drugs[DRUG_NITAZENES].cutoff = 0.5f;
    drugs[DRUG_NITAZENES].dosing_interval = 6.0f;

    strcpy(drugs[DRUG_AMPHETAMINE].name, "AMPHETAMINE");
    drugs[DRUG_AMPHETAMINE].halflife = 12.0f;
    drugs[DRUG_AMPHETAMINE].cutoff = 25.0f;
    drugs[DRUG_AMPHETAMINE].dosing_interval = 12.0f;

    strcpy(drugs[DRUG_METHAMPHETAMINE].name, "METHAMPHETAMINE");
    drugs[DRUG_METHAMPHETAMINE].halflife = 10.0f;
    drugs[DRUG_METHAMPHETAMINE].cutoff = 25.0f;
    drugs[DRUG_METHAMPHETAMINE].dosing_interval = 8.0f;

    strcpy(drugs[DRUG_DEXTROAMPHETAMINE].name, "DEXTROAMPHETAMINE");
    drugs[DRUG_DEXTROAMPHETAMINE].halflife = 11.0f;
    drugs[DRUG_DEXTROAMPHETAMINE].cutoff = 25.0f;
    drugs[DRUG_DEXTROAMPHETAMINE].dosing_interval = 12.0f;

    strcpy(drugs[DRUG_HYDROMORPHONE].name, "HYDROMORPHONE");
    drugs[DRUG_HYDROMORPHONE].halflife = 2.3f;
    drugs[DRUG_HYDROMORPHONE].cutoff = 1.0f;
    drugs[DRUG_HYDROMORPHONE].dosing_interval = 4.0f;

    strcpy(drugs[DRUG_OXYCODONE].name, "OXYCODONE");
    drugs[DRUG_OXYCODONE].halflife = 3.2f;
    drugs[DRUG_OXYCODONE].cutoff = 5.0f;
    drugs[DRUG_OXYCODONE].dosing_interval = 6.0f;

    strcpy(drugs[DRUG_MORPHINE].name, "MORPHINE");
    drugs[DRUG_MORPHINE].halflife = 2.0f;
    drugs[DRUG_MORPHINE].cutoff = 10.0f;
    drugs[DRUG_MORPHINE].dosing_interval = 4.0f;

    strcpy(drugs[DRUG_HYDROCODONE].name, "HYDROCODONE");
    drugs[DRUG_HYDROCODONE].halflife = 3.8f;
    drugs[DRUG_HYDROCODONE].cutoff = 5.0f;
    drugs[DRUG_HYDROCODONE].dosing_interval = 6.0f;

    strcpy(drugs[DRUG_CODEINE].name, "CODEINE");
    drugs[DRUG_CODEINE].halflife = 2.9f;
    drugs[DRUG_CODEINE].cutoff = 10.0f;
    drugs[DRUG_CODEINE].dosing_interval = 6.0f;

    strcpy(drugs[DRUG_PETHIDINE].name, "PETHIDINE");
    drugs[DRUG_PETHIDINE].halflife = 3.2f;
    drugs[DRUG_PETHIDINE].cutoff = 25.0f;
    drugs[DRUG_PETHIDINE].dosing_interval = 6.0f;

    strcpy(drugs[DRUG_BARBITURATES].name, "BARBITURATES");
    drugs[DRUG_BARBITURATES].halflife = 72.0f;
    drugs[DRUG_BARBITURATES].cutoff = 50.0f;
    drugs[DRUG_BARBITURATES].dosing_interval = 24.0f;

    strcpy(drugs[DRUG_BENZODIAZEPINES].name, "BENZODIAZEPINES");
    drugs[DRUG_BENZODIAZEPINES].halflife = 43.0f;
    drugs[DRUG_BENZODIAZEPINES].cutoff = 10.0f;
    drugs[DRUG_BENZODIAZEPINES].dosing_interval = 24.0f;

    strcpy(drugs[DRUG_ALCOHOL].name, "ALCOHOL");
    drugs[DRUG_ALCOHOL].halflife = 1.0f;
    drugs[DRUG_ALCOHOL].cutoff = 5.0f;
    drugs[DRUG_ALCOHOL].dosing_interval = 2.0f;

    strcpy(drugs[DRUG_LSD].name, "LSD");
    drugs[DRUG_LSD].halflife = 3.6f;
    drugs[DRUG_LSD].cutoff = 0.5f;
    drugs[DRUG_LSD].dosing_interval = 12.0f;

    strcpy(drugs[DRUG_KETAMINE].name, "KETAMINE");
    drugs[DRUG_KETAMINE].halflife = 2.5f;
    drugs[DRUG_KETAMINE].cutoff = 25.0f;
    drugs[DRUG_KETAMINE].dosing_interval = 4.0f;

    strcpy(drugs[DRUG_MESCALINE].name, "MESCALINE");
    drugs[DRUG_MESCALINE].halflife = 6.0f;
    drugs[DRUG_MESCALINE].cutoff = 25.0f;
    drugs[DRUG_MESCALINE].dosing_interval = 12.0f;

    strcpy(drugs[DRUG_PSILOCYBIN].name, "PSILOCYBIN");
    drugs[DRUG_PSILOCYBIN].halflife = 2.5f;
    drugs[DRUG_PSILOCYBIN].cutoff = 1.0f;
    drugs[DRUG_PSILOCYBIN].dosing_interval = 8.0f;

    strcpy(drugs[DRUG_DMT].name, "DMT");
    drugs[DRUG_DMT].halflife = 0.25f;
    drugs[DRUG_DMT].cutoff = 1.0f;
    drugs[DRUG_DMT].dosing_interval = 1.0f;

    strcpy(drugs[DRUG_GHB].name, "GHB");
    drugs[DRUG_GHB].halflife = 0.5f;
    drugs[DRUG_GHB].cutoff = 5.0f;
    drugs[DRUG_GHB].dosing_interval = 2.0f;

    strcpy(drugs[DRUG_METHAQUALONE].name, "METHAQUALONE");
    drugs[DRUG_METHAQUALONE].halflife = 24.0f;
    drugs[DRUG_METHAQUALONE].cutoff = 25.0f;
    drugs[DRUG_METHAQUALONE].dosing_interval = 12.0f;

    strcpy(drugs[DRUG_METHADONE].name, "METHADONE");
    drugs[DRUG_METHADONE].halflife = 22.0f;
    drugs[DRUG_METHADONE].cutoff = 25.0f;
    drugs[DRUG_METHADONE].dosing_interval = 24.0f;

    strcpy(drugs[DRUG_DEXTROPROPOXYPHENE].name, "DEXTROPROPOXYPHENE");
    drugs[DRUG_DEXTROPROPOXYPHENE].halflife = 14.0f;
    drugs[DRUG_DEXTROPROPOXYPHENE].cutoff = 10.0f;
    drugs[DRUG_DEXTROPROPOXYPHENE].dosing_interval = 8.0f;
}

void initialize_route_data(void)
{
    /* Initialize route database */
    strcpy(routes[ROUTE_ORAL].name, "ORAL");
    routes[ROUTE_ORAL].bioavailability = 0.7f;
    routes[ROUTE_ORAL].absorption_rate = 1.5f;
    routes[ROUTE_ORAL].oral_factor = 0.01f;

    strcpy(routes[ROUTE_INTRAVENOUS].name, "INTRAVENOUS");
    routes[ROUTE_INTRAVENOUS].bioavailability = 1.0f;
    routes[ROUTE_INTRAVENOUS].absorption_rate = 0.1f;
    routes[ROUTE_INTRAVENOUS].oral_factor = 0.05f;

    strcpy(routes[ROUTE_INTRAMUSCULAR].name, "INTRAMUSCULAR");
    routes[ROUTE_INTRAMUSCULAR].bioavailability = 0.9f;
    routes[ROUTE_INTRAMUSCULAR].absorption_rate = 0.5f;
    routes[ROUTE_INTRAMUSCULAR].oral_factor = 0.03f;

    strcpy(routes[ROUTE_SUBCUTANEOUS].name, "SUBCUTANEOUS");
    routes[ROUTE_SUBCUTANEOUS].bioavailability = 0.8f;
    routes[ROUTE_SUBCUTANEOUS].absorption_rate = 0.8f;
    routes[ROUTE_SUBCUTANEOUS].oral_factor = 0.025f;

    strcpy(routes[ROUTE_INTRANASAL].name, "INTRANASAL");
    routes[ROUTE_INTRANASAL].bioavailability = 0.6f;
    routes[ROUTE_INTRANASAL].absorption_rate = 0.3f;
    routes[ROUTE_INTRANASAL].oral_factor = 0.02f;

    strcpy(routes[ROUTE_INHALATION].name, "INHALATION");
    routes[ROUTE_INHALATION].bioavailability = 0.9f;
    routes[ROUTE_INHALATION].absorption_rate = 0.1f;
    routes[ROUTE_INHALATION].oral_factor = 0.04f;

    strcpy(routes[ROUTE_SUBLINGUAL].name, "SUBLINGUAL");
    routes[ROUTE_SUBLINGUAL].bioavailability = 0.8f;
    routes[ROUTE_SUBLINGUAL].absorption_rate = 0.5f;
    routes[ROUTE_SUBLINGUAL].oral_factor = 0.02f;

    strcpy(routes[ROUTE_TRANSDERMAL].name, "TRANSDERMAL");
    routes[ROUTE_TRANSDERMAL].bioavailability = 0.9f;
    routes[ROUTE_TRANSDERMAL].absorption_rate = 4.0f;
    routes[ROUTE_TRANSDERMAL].oral_factor = 0.015f;

    strcpy(routes[ROUTE_RECTAL].name, "RECTAL");
    routes[ROUTE_RECTAL].bioavailability = 0.7f;
    routes[ROUTE_RECTAL].absorption_rate = 1.0f;
    routes[ROUTE_RECTAL].oral_factor = 0.015f;

    strcpy(routes[ROUTE_BUCCAL].name, "BUCCAL");
    routes[ROUTE_BUCCAL].bioavailability = 0.75f;
    routes[ROUTE_BUCCAL].absorption_rate = 0.8f;
    routes[ROUTE_BUCCAL].oral_factor = 0.025f;

    strcpy(routes[ROUTE_TOPICAL].name, "TOPICAL");
    routes[ROUTE_TOPICAL].bioavailability = 0.1f;
    routes[ROUTE_TOPICAL].absorption_rate = 8.0f;
    routes[ROUTE_TOPICAL].oral_factor = 0.005f;
}

void print_banner(void)
{
    printf("====================================================================\n");
    printf("NARCDETECT - Drug Detection Time Calculator v2.0\n");
    printf("FOR ORAL FLUID (SALIVA) TESTING\n");
    printf("====================================================================\n\n");
    printf("ESTIMATES TIME UNTIL NON-DETECTABLE\n");
    printf("BASED ON PHARMACOKINETIC PARAMETERS\n");
    printf("INCLUDES CHRONIC USE ACCUMULATION\n");
    printf("AND ROUTES OF ADMINISTRATION\n");
    printf("WITH NMR SPECTRUM SIMULATION\n\n");
}

void print_drug_menu(void)
{
    printf("AVAILABLE DRUGS:\n");
    printf("FENTANYL, NITAZENES, AMPHETAMINE,\n");
    printf("METHAMPHETAMINE, DEXTROAMPHETAMINE,\n");
    printf("HYDROMORPHONE, OXYCODONE, MORPHINE,\n");
    printf("HYDROCODONE, CODEINE, PETHIDINE,\n");
    printf("BARBITURATES, BENZODIAZEPINES, ALCOHOL,\n");
    printf("LSD, KETAMINE, MESCALINE, PSILOCYBIN,\n");
    printf("DMT, GHB, METHAQUALONE, METHADONE,\n");
    printf("DEXTROPROPOXYPHENE\n\n");
    printf("Enter drug name: ");
}

void print_route_menu(void)
{
    printf("\nAVAILABLE ROUTES OF ADMINISTRATION:\n");
    printf("ORAL, INTRAVENOUS, INTRAMUSCULAR,\n");
    printf("SUBCUTANEOUS, INTRANASAL, INHALATION,\n");
    printf("SUBLINGUAL, TRANSDERMAL, RECTAL,\n");
    printf("BUCCAL, TOPICAL\n\n");
    printf("Enter route of administration: ");
}

int get_drug_selection(void)
{
    char input[50];
    int i;

    print_drug_menu();
    fgets(input, sizeof(input), stdin);
    
    /* Remove newline */
    input[strcspn(input, "\n")] = 0;
    str_upper(input);

    /* Search for drug */
    for (i = 1; i <= NUM_DRUGS; i++) {
        if (str_compare_upper(input, drugs[i].name) == 0) {
            return i;
        }
    }

    /* Check common alternatives */
    if (str_compare_upper(input, "MEPERIDINE") == 0) return DRUG_PETHIDINE;
    if (str_compare_upper(input, "ETHANOL") == 0) return DRUG_ALCOHOL;
    if (str_compare_upper(input, "PROPOXYPHENE") == 0) return DRUG_DEXTROPROPOXYPHENE;

    return 0; /* Not found */
}

int get_route_selection(void)
{
    char input[50];
    int i;

    print_route_menu();
    fgets(input, sizeof(input), stdin);
    
    /* Remove newline */
    input[strcspn(input, "\n")] = 0;
    str_upper(input);

    /* Search for route */
    for (i = 1; i <= NUM_ROUTES; i++) {
        if (str_compare_upper(input, routes[i].name) == 0) {
            return i;
        }
    }

    /* Check common alternatives */
    if (str_compare_upper(input, "IV") == 0) return ROUTE_INTRAVENOUS;
    if (str_compare_upper(input, "IM") == 0) return ROUTE_INTRAMUSCULAR;
    if (str_compare_upper(input, "SC") == 0) return ROUTE_SUBCUTANEOUS;
    if (str_compare_upper(input, "NASAL") == 0) return ROUTE_INTRANASAL;
    if (str_compare_upper(input, "INHALED") == 0) return ROUTE_INHALATION;
    if (str_compare_upper(input, "SMOKING") == 0) return ROUTE_INHALATION;
    if (str_compare_upper(input, "SL") == 0) return ROUTE_SUBLINGUAL;
    if (str_compare_upper(input, "PATCH") == 0) return ROUTE_TRANSDERMAL;
    if (str_compare_upper(input, "PR") == 0) return ROUTE_RECTAL;

    return 0; /* Not found */
}

void get_input_parameters(int *dosage, int *weight, int *age, int *metab, float *duration)
{
    printf("\nEnter dosage in mg: ");
    scanf("%d", dosage);

    printf("Enter body weight in kg: ");
    scanf("%d", weight);

    printf("Enter age in years: ");
    scanf("%d", age);

    printf("Metabolism rate (1=SLOW, 2=NORMAL, 3=FAST): ");
    scanf("%d", metab);

    printf("Duration of use in hours (24.0=1 day): ");
    scanf("%f", duration);
}

void adjust_route_parameters(int drug, int route, float *bioavail, float *oral_fac, float *absorpt)
{
    /* Drug-specific route adjustments */
    
    /* Alcohol adjustments */
    if (drug == DRUG_ALCOHOL) {
        if (route == ROUTE_INTRAVENOUS || route == ROUTE_INTRAMUSCULAR || 
            route == ROUTE_SUBCUTANEOUS) {
            *bioavail *= 0.1f; /* IV/IM/SC not typical for alcohol */
        }
        if (route == ROUTE_INHALATION) {
            *bioavail = 0.95f; /* Vapor highly bioavailable */
            *absorpt = 0.05f;
        }
    }

    /* Fentanyl adjustments */
    if (drug == DRUG_FENTANYL) {
        if (route == ROUTE_TRANSDERMAL) {
            *absorpt = 12.0f; /* Sustained release */
            *bioavail = 0.92f;
        }
        if (route == ROUTE_SUBLINGUAL) {
            *bioavail = 0.8f; /* High bioavailability */
        }
    }

    /* Stimulants (amphetamines) adjustments */
    if (drug >= DRUG_AMPHETAMINE && drug <= DRUG_DEXTROAMPHETAMINE) {
        if (route == ROUTE_INTRANASAL) {
            *bioavail = 0.8f; /* Common route */
            *absorpt = 0.2f;
        }
        if (route == ROUTE_INHALATION) {
            *bioavail = 0.7f; /* Smoking */
            *absorpt = 0.08f;
        }
    }

    /* Opioids adjustments */
    if ((drug >= DRUG_HYDROMORPHONE && drug <= DRUG_PETHIDINE) || drug == DRUG_METHADONE) {
        if (route == ROUTE_INTRAVENOUS) {
            *bioavail = 1.0f; /* Full bioavailability */
            *oral_fac = 0.08f;
        }
        if (route == ROUTE_INTRANASAL) {
            *bioavail = 0.65f; /* Moderate bioavailability */
        }
    }

    /* Psychedelics adjustments */
    if (drug >= DRUG_LSD && drug <= DRUG_DMT) {
        if (route == ROUTE_INHALATION && drug != DRUG_DMT) {
            *bioavail *= 0.3f; /* Not typical for most psychedelics */
        }
        if (drug == DRUG_DMT && route == ROUTE_INHALATION) {
            *bioavail = 0.8f; /* DMT typically smoked */
            *absorpt = 0.02f;
        }
    }

    /* Benzodiazepines adjustments */
    if (drug == DRUG_BENZODIAZEPINES) {
        if (route == ROUTE_SUBLINGUAL) {
            *bioavail = 0.9f; /* Good absorption */
            *absorpt = 0.3f;
        }
        if (route == ROUTE_RECTAL) {
            *bioavail = 0.8f; /* Good for seizure control */
            *absorpt = 0.5f;
        }
    }

    /* Ketamine adjustments */
    if (drug == DRUG_KETAMINE) {
        if (route == ROUTE_INTRANASAL) {
            *bioavail = 0.5f; /* Common route */
            *absorpt = 0.3f;
        }
        if (route == ROUTE_INTRAMUSCULAR) {
            *bioavail = 0.93f; /* Clinical use */
            *absorpt = 0.3f;
        }
    }

    /* GHB adjustments */
    if (drug == DRUG_GHB && route != ROUTE_ORAL) {
        *bioavail *= 0.5f; /* Primarily oral */
    }

    /* Topical route restrictions */
    if (route == ROUTE_TOPICAL) {
        if (drug != DRUG_FENTANYL && drug != DRUG_METHADONE) {
            *bioavail = 0.05f; /* Low systemic absorption */
            *oral_fac = 0.002f;
        }
    }
}

void calculate_detection_time(int drug, int route, int dosage, int weight, int age, 
                            int metab, float duration)
{
    float halflife, cutoff, dosing_interval;
    float bioavail, absorpt, oral_fac;
    float single_conc, total_conc, steady_conc;
    float elim_rate, age_factor, metab_factor;
    float accumulation_factor, r_factor;
    float detection_time, buildup;
    int num_doses;
    int hours, minutes, seconds, days;

    /* Get drug parameters */
    halflife = drugs[drug].halflife;
    cutoff = drugs[drug].cutoff;
    dosing_interval = drugs[drug].dosing_interval;

    /* Get route parameters */
    bioavail = routes[route].bioavailability;
    absorpt = routes[route].absorption_rate;
    oral_fac = routes[route].oral_factor;

    /* Apply route adjustments */
    adjust_route_parameters(drug, route, &bioavail, &oral_fac, &absorpt);

    /* Calculate single dose concentration */
    if (drug == DRUG_FENTANYL) {
        single_conc = fentanyl_dose_constant * 1000.0f * oral_fac * bioavail / (float)weight;
    } else if (drug == DRUG_ALCOHOL) {
        single_conc = (float)dosage * oral_fac * bioavail * 0.5f / (float)weight;
    } else {
        single_conc = (float)dosage * oral_fac * bioavail / (float)weight;
    }

    /* Adjust half-life for absorption rate (flip-flop kinetics) */
    if (absorpt > halflife * 0.693f) {
        halflife = halflife * (1.0f + absorpt / (halflife * 0.693f));
    }

    /* Age factor adjustment */
    if (age < 35) age_factor = 1.15f;
    else if (age < 50) age_factor = 1.0f;
    else if (age < 65) age_factor = 0.85f;
    else age_factor = 0.7f;

    /* Metabolism factor */
    if (metab == 1) metab_factor = 0.7f;      /* Slow */
    else if (metab == 2) metab_factor = 1.0f; /* Normal */
    else metab_factor = 1.4f;                 /* Fast */

    /* Calculate elimination rate */
    elim_rate = 0.693f / halflife;
    elim_rate = elim_rate * age_factor * metab_factor;

    /* Calculate accumulation */
    num_doses = (int)(duration / dosing_interval) + 1;
    r_factor = 1.0f - exp(-elim_rate * dosing_interval);

    if (fabs(r_factor - 1.0f) < 0.001f) {
        accumulation_factor = (float)num_doses;
    } else {
        accumulation_factor = (1.0f - pow(r_factor, (float)num_doses)) / (1.0f - r_factor);
    }

    /* Calculate concentrations */
    total_conc = single_conc * accumulation_factor;
    steady_conc = single_conc / (1.0f - exp(-elim_rate * dosing_interval));
    buildup = (total_conc / steady_conc) * 100.0f;
    if (buildup > 100.0f) buildup = 100.0f;

    /* Calculate detection time */
    if (total_conc > cutoff) {
        detection_time = log(total_conc / cutoff) / elim_rate;
    } else {
        detection_time = 0.0f;
    }

    /* Plot concentration curve */
    plot_concentration_curve(total_conc, elim_rate, cutoff, halflife, duration);

    /* Display results */
    printf("\n====================================================================\n");
    printf("DETECTION TIME CALCULATION FOR %s\n", drugs[drug].name);
    printf("====================================================================\n\n");

    printf("INPUT PARAMETERS:\n");
    printf("  Dosage: %d mg\n", dosage);
    printf("  Weight: %d kg\n", weight);
    printf("  Age: %d years\n", age);
    printf("  Metabolism: %s\n", (metab == 1) ? "SLOW" : (metab == 2) ? "NORMAL" : "FAST");
    printf("  Duration of use: %.1f hours (%.2f days)\n", duration, duration / 24.0f);
    printf("  Route: %s (Bioavail %.1f%%, Abs rate %.2f hr)\n", 
           routes[route].name, bioavail * 100.0f, absorpt);

    if (drug == DRUG_FENTANYL) {
        printf("  Fentanyl dose: %.0f mg (constant)\n", fentanyl_dose_constant * 1000.0f);
    }

    printf("\nPHARMACOKINETIC DATA:\n");
    printf("  Half-life: %.1f hours\n", halflife);
    printf("  Cutoff: %.1f ng/mL\n", cutoff);
    printf("  Dosing interval: %.1f hours\n", dosing_interval);
    printf("  Number of doses: %d\n", num_doses);
    printf("  Accumulation factor: %.2f\n", accumulation_factor);
    printf("  Single dose conc: %.2f ng/mL\n", single_conc);
    printf("  Total accum conc: %.2f ng/mL\n", total_conc);
    printf("  Elim rate: %.4f /hour\n", elim_rate);
    printf("  Steady-state conc: %.2f ng/mL\n", steady_conc);
    printf("  Buildup to SS: %.1f%%\n", buildup);

    /* Convert detection time to readable format */
    seconds = (int)(detection_time * 3600.0f);
    hours = seconds / 3600;
    minutes = (seconds - hours * 3600) / 60;
    seconds = seconds - hours * 3600 - minutes * 60;
    days = hours / 24;
    hours = hours - days * 24;

    printf("\nDETECTION TIME: %.0f seconds\n", detection_time * 3600.0f);
    printf("EQUIVALENT TO: %d hours, %d minutes, %d seconds\n", 
           (int)(detection_time * 3600.0f) / 3600, minutes, seconds);
    printf("FULL FORMAT: %d days, %d hours, %d minutes, %d seconds\n", 
           days, hours, minutes, seconds);

    /* Disclaimers */
    printf("\n** IMPORTANT DISCLAIMERS **\n");
    printf("- Estimates based on population averages\n");
    printf("- Individual variation can be significant\n");
    printf("- Chronic use calculations are simplified\n");
    printf("- Assumes regular dosing intervals\n");
    printf("- Route-specific parameters are estimates\n");
    printf("- For research/educational use only\n");
}

void plot_concentration_curve(float c0, float kelim, float cutoff, float thalf, float duration)
{
    float conc[61], time[61];
    char plot_line[PLOT_WIDTH + 1];
    float tmax, dt, cmax;
    int i, j, pos, cutoff_pos;

    printf("\n====================================================================\n");
    printf("  PLASMA CONCENTRATION vs TIME WITH ACCUMULATION\n");
    printf("       (INCLUDES CHRONIC USE BUILD-UP EFFECTS)\n");
    printf("       (ADJUSTED FOR ROUTE OF ADMINISTRATION)\n");
    printf("====================================================================\n\n");

    /* Calculate time points */
    tmax = max_float(8.0f * thalf, duration * 2.0f);
    dt = tmax / 60.0f;

    /* Calculate concentrations */
    cmax = c0;
    for (i = 0; i < 61; i++) {
        time[i] = (float)i * dt;
        conc[i] = c0 * exp(-kelim * time[i]);
    }

    printf("Time range: 0 to %.1f hours\n", tmax);
    printf("Maximum concentration: %.2f ng/mL\n", cmax);
    printf("Cutoff level: %.2f ng/mL\n\n", cutoff);

    /* Plot the curve */
    for (i = 0; i < 61; i++) {
        /* Clear plot line */
        for (j = 0; j < PLOT_WIDTH; j++) {
            plot_line[j] = ' ';
        }
        plot_line[PLOT_WIDTH] = '\0';

        /* Add grid markers every 10 positions */
        if (i % 6 == 0) {
            for (j = 9; j < PLOT_WIDTH; j += 10) {
                plot_line[j] = '+';
            }
        }

        /* Plot concentration point */
        pos = (int)(conc[i] * (PLOT_WIDTH - 2) / cmax);
        if (pos >= 0 && pos < PLOT_WIDTH) {
            plot_line[pos] = 'A';
        }

        /* Plot cutoff line */
        cutoff_pos = (int)(cutoff * (PLOT_WIDTH - 2) / cmax);
        if (cutoff_pos >= 0 && cutoff_pos < PLOT_WIDTH && cutoff_pos != pos) {
            plot_line[cutoff_pos] = 'C';
        }

        printf("%s\n", plot_line);
    }

    printf("\nLEGEND: A = CONCENTRATION CURVE\n");
    printf("        C = DETECTION CUTOFF THRESHOLD\n");
    printf("        + = TIME GRID MARKERS\n\n");

    if (c0 > cutoff) {
        float final_time = log(c0 / cutoff) / kelim;
        printf("ANALYSIS: Time to non-detection = %.2f hours\n", final_time);
        printf("          Peak concentration = %.2f ng/mL\n", c0);
        printf("          Duration of use = %.2f hours\n\n", duration);
    }
}

void generate_nmr_data(int drug, NMRData *nmr_data)
{
    int i;

    /* Initialize */
    nmr_data->num_peaks = 0;
    for (i = 0; i < MAX_PEAKS; i++) {
        nmr_data->shifts[i] = 0.0f;
        nmr_data->intensities[i] = 0.0f;
        nmr_data->widths[i] = 0.1f;
    }

    /* Generate drug-specific NMR data */
    switch (drug) {
        case DRUG_FENTANYL:
            nmr_data->num_peaks = 4;
            nmr_data->shifts[0] = 7.2f;  /* Phenyl protons */
            nmr_data->shifts[1] = 3.8f;  /* N-CH3 protons */
            nmr_data->shifts[2] = 2.4f;  /* Piperidine protons */
            nmr_data->shifts[3] = 1.2f;  /* Ethyl protons */
            nmr_data->intensities[0] = 100.0f;
            nmr_data->intensities[1] = 150.0f;
            nmr_data->intensities[2] = 200.0f;
            nmr_data->intensities[3] = 120.0f;
            break;

        case DRUG_AMPHETAMINE:
        case DRUG_METHAMPHETAMINE:
        case DRUG_DEXTROAMPHETAMINE:
            nmr_data->num_peaks = 4;
            nmr_data->shifts[0] = 7.3f;  /* Phenyl protons */
            nmr_data->shifts[1] = 2.8f;  /* CH2-Phenyl */
            nmr_data->shifts[2] = 3.1f;  /* CH-NH2 */
            nmr_data->shifts[3] = 1.1f;  /* CH3 */
            nmr_data->intensities[0] = 100.0f;
            nmr_data->intensities[1] = 80.0f;
            nmr_data->intensities[2] = 60.0f;
            nmr_data->intensities[3] = (drug == DRUG_METHAMPHETAMINE) ? 90.0f : 0.0f;
            if (drug != DRUG_METHAMPHETAMINE) nmr_data->num_peaks = 3;
            break;

        case DRUG_MORPHINE:
        case DRUG_HYDROMORPHONE:
        case DRUG_OXYCODONE:
        case DRUG_HYDROCODONE:
        case DRUG_CODEINE:
            nmr_data->num_peaks = 5;
            nmr_data->shifts[0] = 6.8f;  /* Aromatic */
            nmr_data->shifts[1] = 6.5f;  /* Aromatic */
            nmr_data->shifts[2] = 4.2f;  /* CH-OH */
            nmr_data->shifts[3] = 3.0f;  /* N-CH3 */
            nmr_data->shifts[4] = 2.1f;  /* CH2 */
            nmr_data->intensities[0] = 50.0f;
            nmr_data->intensities[1] = 50.0f;
            nmr_data->intensities[2] = 60.0f;
            nmr_data->intensities[3] = 90.0f;
            nmr_data->intensities[4] = 100.0f;
            break;

        case DRUG_KETAMINE:
            nmr_data->num_peaks = 3;
            nmr_data->shifts[0] = 7.5f;  /* Aromatic */
            nmr_data->shifts[1] = 4.1f;  /* CH-N */
            nmr_data->shifts[2] = 2.5f;  /* CH2 */
            nmr_data->intensities[0] = 80.0f;
            nmr_data->intensities[1] = 60.0f;
            nmr_data->intensities[2] = 100.0f;
            break;

        case DRUG_LSD:
            nmr_data->num_peaks = 6;
            nmr_data->shifts[0] = 8.1f;  /* Indole NH */
            nmr_data->shifts[1] = 7.4f;  /* Aromatic */
            nmr_data->shifts[2] = 7.0f;  /* Aromatic */
            nmr_data->shifts[3] = 6.8f;  /* Aromatic */
            nmr_data->shifts[4] = 4.0f;  /* CH-N */
            nmr_data->shifts[5] = 1.3f;  /* CH3 */
            nmr_data->intensities[0] = 30.0f;
            nmr_data->intensities[1] = 50.0f;
            nmr_data->intensities[2] = 50.0f;
            nmr_data->intensities[3] = 50.0f;
            nmr_data->intensities[4] = 60.0f;
            nmr_data->intensities[5] = 90.0f;
            break;

        default:
            /* Generic spectrum */
            nmr_data->num_peaks = 3;
            nmr_data->shifts[0] = 7.0f;  /* Aromatic region */
            nmr_data->shifts[1] = 3.5f;  /* Aliphatic CH */
            nmr_data->shifts[2] = 1.5f;  /* Methyl region */
            nmr_data->intensities[0] = 100.0f;
            nmr_data->intensities[1] = 80.0f;
            nmr_data->intensities[2] = 120.0f;
            break;
    }

    /* Set default widths */
    for (i = 0; i < nmr_data->num_peaks; i++) {
        nmr_data->widths[i] = 0.08f + (float)(rand() % 20) / 1000.0f; /* 0.08-0.10 */
    }
}

void nmr_plot(int drug, float concentration, NMRData *nmr_data)
{
    float spectrum[SPECTRUM_WIDTH];
    float freq[SPECTRUM_WIDTH];
    char plot_line[SPECTRUM_WIDTH + 1];
    float spec_max, thresh, delta, peak_val, width, intensity;
    int i, j, line, pos;
    char peak_labels[MAX_PEAKS][25];

    printf("\n====================================================================\n");
    printf("          1H NMR SPECTRUM SIMULATION FOR %s\n", drugs[drug].name);
    printf("       CONCENTRATION: %.2f NG/ML IN SAMPLE\n", concentration);
    printf("       CHEMICAL SHIFT RANGE: 0.0 - 12.0 PPM\n");
    printf("       SYNTHETIC SPECTRUM FOR IDENTIFICATION\n");
    printf("====================================================================\n\n");

    /* Initialize spectrum array (12.0 to 0.0 PPM) */
    for (i = 0; i < SPECTRUM_WIDTH; i++) {
        freq[i] = 12.0f - (float)i * 0.1f;
        spectrum[i] = 0.0f;
    }

    /* Add peaks to spectrum */
    for (j = 0; j < nmr_data->num_peaks; j++) {
        if (nmr_data->shifts[j] >= 0.0f && nmr_data->shifts[j] <= 12.0f) {
            width = max_float(0.05f, nmr_data->widths[j]);
            intensity = nmr_data->intensities[j] * concentration / 100.0f;

            for (i = 0; i < SPECTRUM_WIDTH; i++) {
                delta = fabs(freq[i] - nmr_data->shifts[j]);
                peak_val = intensity / (1.0f + pow(delta / width, 2.0f));
                spectrum[i] += peak_val;
            }
        }
    }

    /* Find maximum for scaling */
    spec_max = 0.0f;
    for (i = 0; i < SPECTRUM_WIDTH; i++) {
        spec_max = max_float(spec_max, spectrum[i]);
    }
    if (spec_max <= 0.0f) spec_max = 1.0f;

    /* Print scale information */
    printf("Maximum intensity = %.2f (relative)\n", spec_max);
    printf("Chemical shift scale: 12.0 to 0.0 PPM\n\n");
    printf("      12        10         8         6         4         2         0\n");
    printf("       |         |         |         |         |         |         |\n");

    /* Plot spectrum (50 lines, top to bottom) */
    for (line = PLOT_HEIGHT; line >= 1; line--) {
        thresh = spec_max * (float)line / (float)PLOT_HEIGHT;

        /* Clear plot line */
        for (i = 0; i < SPECTRUM_WIDTH; i++) {
            plot_line[i] = ' ';
        }
        plot_line[SPECTRUM_WIDTH] = '\0';

        /* Plot spectrum points above threshold */
        for (i = 0; i < SPECTRUM_WIDTH; i++) {
            if (spectrum[i] >= thresh) {
                plot_line[i] = '*';
            }
            /* Add grid lines at PPM marks (every 20 points = 2 PPM) */
            if (i % 20 == 0 && plot_line[i] == ' ') {
                plot_line[i] = '+';
            }
        }

        printf("%s\n", plot_line);
    }

    /* Print peak assignments */
    if (nmr_data->num_peaks > 0) {
        printf("\nPEAK ASSIGNMENTS:\n");
        printf("SHIFT(PPM)  INTENSITY  WIDTH   ASSIGNMENT\n");
        printf("----------  ---------  -----   ----------\n");

        for (j = 0; j < nmr_data->num_peaks; j++) {
            if (nmr_data->shifts[j] >= 0.0f && nmr_data->shifts[j] <= 12.0f) {
                get_peak_label(drug, j + 1, nmr_data->shifts[j], peak_labels[j]);
                printf("%8.2f    %7.1f    %5.2f   %s\n",
                       nmr_data->shifts[j], nmr_data->intensities[j],
                       nmr_data->widths[j], peak_labels[j]);
            }
        }
    }

    printf("\nSPECTRUM ANALYSIS:\n");
    printf("NUMBER OF PEAKS DETECTED: %d\n", nmr_data->num_peaks);
    printf("MAXIMUM PEAK INTENSITY:   %.2f\n", spec_max);
    printf("SAMPLE CONCENTRATION:     %.2f NG/ML\n", concentration);
    printf("INTEGRATION COMPLETE\n\n");
    printf("* = SPECTRAL PEAK    + = PPM GRID LINES\n\n");
}

void get_peak_label(int drug, int peak_no, float shift, char *label)
{
    /* Generic labels based on chemical shift regions */
    if (shift >= 10.0f) {
        strcpy(label, "AROMATIC H");
    } else if (shift >= 7.0f) {
        strcpy(label, "AROMATIC/VINYL H");
    } else if (shift >= 4.0f) {
        strcpy(label, "O-CH, N-CH");
    } else if (shift >= 2.0f) {
        strcpy(label, "CH2, CH3 ALPHA");
    } else if (shift >= 1.0f) {
        strcpy(label, "CH2, CH3 BETA");
    } else {
        strcpy(label, "CH3 ALIPHATIC");
    }

    /* Drug-specific assignments */
    if (drug == DRUG_FENTANYL) {
        switch (peak_no) {
            case 1: strcpy(label, "PHENYL H"); break;
            case 2: strcpy(label, "FENTANYL N-CH3"); break;
            case 3: strcpy(label, "PIPERIDINE H"); break;
            case 4: strcpy(label, "ETHYL H"); break;
        }
    } else if (drug >= DRUG_AMPHETAMINE && drug <= DRUG_DEXTROAMPHETAMINE) {
        switch (peak_no) {
            case 1: strcpy(label, "PHENYL H"); break;
            case 2: strcpy(label, "CH2-PHENYL"); break;
            case 3: strcpy(label, "CH-NH2"); break;
            case 4: strcpy(label, "CH3 (IF METH)"); break;
        }
    }
}

/* Utility functions */
void str_upper(char *str)
{
    int i;
    for (i = 0; str[i]; i++) {
        str[i] = toupper(str[i]);
    }
}

int str_compare_upper(const char *str1, const char *str2)
{
    char s1[100], s2[100];
    strcpy(s1, str1);
    strcpy(s2, str2);
    str_upper(s1);
    str_upper(s2);
    return strcmp(s1, s2);
}

float max_float(float a, float b)
{
    return (a > b) ? a : b;
}

float min_float(float a, float b)
{
    return (a < b) ? a : b;
}

int max_int(int a, int b)
{
    return (a > b) ? a : b;
}

int min_int(int a, int b)
{
    return (a < b) ? a : b;
}
