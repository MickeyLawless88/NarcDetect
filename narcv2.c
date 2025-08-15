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
#define NUM_DRUGS 24
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
    DRUG_DEXTROPROPOXYPHENE = 23,
    DRUG_DIAMORPHINE = 24
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
    float halflife_saliva;      /* Half-life in oral fluid (hours) */
    float halflife_urine;       /* Half-life in urine (hours) */
    float cutoff_saliva;        /* Cutoff concentration in ng/mL (saliva) */
    float cutoff_urine;         /* Cutoff concentration in ng/mL (urine) */
    float dosing_interval;
    char metabolite_info[100];  /* Primary metabolites detected */
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
void plot_concentration_curve(float c0, float kelim, float cutoff, float thalf, float duration, float dosing_interval, float single_dose_conc, float absorption_rate);
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
    /* Initialize drug database with realistic pharmacokinetic parameters for both matrices */
    
    /* FENTANYL - Highly potent synthetic opioid */
    strcpy(drugs[DRUG_FENTANYL].name, "FENTANYL");
    drugs[DRUG_FENTANYL].halflife_saliva = 7.0f;   /* Saliva detection: 8-24 hours */
    drugs[DRUG_FENTANYL].halflife_urine = 20.0f;   /* Urine detection: 1-3 days */
    drugs[DRUG_FENTANYL].cutoff_saliva = 1.0f;     /* 1 ng/mL saliva */
    drugs[DRUG_FENTANYL].cutoff_urine = 2.0f;      /* 2 ng/mL urine */
    drugs[DRUG_FENTANYL].dosing_interval = 4.0f;
    strcpy(drugs[DRUG_FENTANYL].metabolite_info, "Parent drug + norfentanyl");

    /* NITAZENES - New synthetic opioids, similar to fentanyl */
    strcpy(drugs[DRUG_NITAZENES].name, "NITAZENES");
    drugs[DRUG_NITAZENES].halflife_saliva = 8.0f;
    drugs[DRUG_NITAZENES].halflife_urine = 24.0f;
    drugs[DRUG_NITAZENES].cutoff_saliva = 0.5f;    /* More potent */
    drugs[DRUG_NITAZENES].cutoff_urine = 1.0f;
    drugs[DRUG_NITAZENES].dosing_interval = 6.0f;
    strcpy(drugs[DRUG_NITAZENES].metabolite_info, "Parent drug + hydroxy metabolites");

    /* AMPHETAMINE - Classic stimulant, longer detection */
    strcpy(drugs[DRUG_AMPHETAMINE].name, "AMPHETAMINE");
    drugs[DRUG_AMPHETAMINE].halflife_saliva = 8.0f;  /* 1-3 days saliva */
    drugs[DRUG_AMPHETAMINE].halflife_urine = 30.0f;  /* 1-4 days urine */
    drugs[DRUG_AMPHETAMINE].cutoff_saliva = 50.0f;
    drugs[DRUG_AMPHETAMINE].cutoff_urine = 500.0f;   /* Higher urine cutoff */
    drugs[DRUG_AMPHETAMINE].dosing_interval = 12.0f;
    strcpy(drugs[DRUG_AMPHETAMINE].metabolite_info, "Unchanged drug (80%) + metabolites");

    /* METHAMPHETAMINE - Longer detection than amphetamine */
    strcpy(drugs[DRUG_METHAMPHETAMINE].name, "METHAMPHETAMINE");
    drugs[DRUG_METHAMPHETAMINE].halflife_saliva = 12.0f; /* 1-4 days saliva */
    drugs[DRUG_METHAMPHETAMINE].halflife_urine = 36.0f;  /* 3-6 days urine */
    drugs[DRUG_METHAMPHETAMINE].cutoff_saliva = 50.0f;
    drugs[DRUG_METHAMPHETAMINE].cutoff_urine = 500.0f;
    drugs[DRUG_METHAMPHETAMINE].dosing_interval = 8.0f;
    strcpy(drugs[DRUG_METHAMPHETAMINE].metabolite_info, "Parent drug + amphetamine metabolite");

    /* DEXTROAMPHETAMINE - Medical amphetamine */
    strcpy(drugs[DRUG_DEXTROAMPHETAMINE].name, "DEXTROAMPHETAMINE");
    drugs[DRUG_DEXTROAMPHETAMINE].halflife_saliva = 9.0f;
    drugs[DRUG_DEXTROAMPHETAMINE].halflife_urine = 32.0f;
    drugs[DRUG_DEXTROAMPHETAMINE].cutoff_saliva = 50.0f;
    drugs[DRUG_DEXTROAMPHETAMINE].cutoff_urine = 500.0f;
    drugs[DRUG_DEXTROAMPHETAMINE].dosing_interval = 12.0f;
    strcpy(drugs[DRUG_DEXTROAMPHETAMINE].metabolite_info, "Unchanged drug + hydroxylated metabolites");

    /* HYDROMORPHONE - Semi-synthetic opioid */
    strcpy(drugs[DRUG_HYDROMORPHONE].name, "HYDROMORPHONE");
    drugs[DRUG_HYDROMORPHONE].halflife_saliva = 3.0f;  /* 12-36 hours */
    drugs[DRUG_HYDROMORPHONE].halflife_urine = 11.0f;  /* 1-3 days */
    drugs[DRUG_HYDROMORPHONE].cutoff_saliva = 1.0f;
    drugs[DRUG_HYDROMORPHONE].cutoff_urine = 10.0f;
    drugs[DRUG_HYDROMORPHONE].dosing_interval = 4.0f;
    strcpy(drugs[DRUG_HYDROMORPHONE].metabolite_info, "Parent drug + hydromorphone-3-glucuronide");

    /* OXYCODONE - Semi-synthetic opioid */
    strcpy(drugs[DRUG_OXYCODONE].name, "OXYCODONE");
    drugs[DRUG_OXYCODONE].halflife_saliva = 4.5f;  /* 1-2 days */
    drugs[DRUG_OXYCODONE].halflife_urine = 19.0f;  /* 1-4 days */
    drugs[DRUG_OXYCODONE].cutoff_saliva = 5.0f;
    drugs[DRUG_OXYCODONE].cutoff_urine = 100.0f;
    drugs[DRUG_OXYCODONE].dosing_interval = 6.0f;
    strcpy(drugs[DRUG_OXYCODONE].metabolite_info, "Parent drug + oxymorphone + glucuronides");

    /* MORPHINE - Natural opioid, main heroin metabolite */
    strcpy(drugs[DRUG_MORPHINE].name, "MORPHINE");
    drugs[DRUG_MORPHINE].halflife_saliva = 3.5f;   /* 1-3 days */
    drugs[DRUG_MORPHINE].halflife_urine = 15.0f;   /* 1-4 days */
    drugs[DRUG_MORPHINE].cutoff_saliva = 10.0f;
    drugs[DRUG_MORPHINE].cutoff_urine = 300.0f;    /* Higher screening cutoff */
    drugs[DRUG_MORPHINE].dosing_interval = 4.0f;
    strcpy(drugs[DRUG_MORPHINE].metabolite_info, "Parent drug + morphine-3-glucuronide + M6G");

    /* HYDROCODONE - Semi-synthetic opioid */
    strcpy(drugs[DRUG_HYDROCODONE].name, "HYDROCODONE");
    drugs[DRUG_HYDROCODONE].halflife_saliva = 4.0f;
    drugs[DRUG_HYDROCODONE].halflife_urine = 18.0f; /* 1-4 days */
    drugs[DRUG_HYDROCODONE].cutoff_saliva = 5.0f;
    drugs[DRUG_HYDROCODONE].cutoff_urine = 100.0f;
    drugs[DRUG_HYDROCODONE].dosing_interval = 6.0f;
    strcpy(drugs[DRUG_HYDROCODONE].metabolite_info, "Parent drug + hydromorphone + glucuronides");

    /* CODEINE - Natural opioid, metabolizes to morphine */
    strcpy(drugs[DRUG_CODEINE].name, "CODEINE");
    drugs[DRUG_CODEINE].halflife_saliva = 3.0f;
    drugs[DRUG_CODEINE].halflife_urine = 12.0f;    /* 1-2 days */
    drugs[DRUG_CODEINE].cutoff_saliva = 10.0f;
    drugs[DRUG_CODEINE].cutoff_urine = 300.0f;
    drugs[DRUG_CODEINE].dosing_interval = 6.0f;
    strcpy(drugs[DRUG_CODEINE].metabolite_info, "Parent drug + morphine + norcodeine");

    /* PETHIDINE/MEPERIDINE - Synthetic opioid */
    strcpy(drugs[DRUG_PETHIDINE].name, "PETHIDINE");
    drugs[DRUG_PETHIDINE].halflife_saliva = 4.0f;
    drugs[DRUG_PETHIDINE].halflife_urine = 16.0f;  /* 1-4 days */
    drugs[DRUG_PETHIDINE].cutoff_saliva = 25.0f;
    drugs[DRUG_PETHIDINE].cutoff_urine = 200.0f;
    drugs[DRUG_PETHIDINE].dosing_interval = 6.0f;
    strcpy(drugs[DRUG_PETHIDINE].metabolite_info, "Parent drug + norpethidine");

    /* BARBITURATES - Long-acting CNS depressants */
    strcpy(drugs[DRUG_BARBITURATES].name, "BARBITURATES");
    drugs[DRUG_BARBITURATES].halflife_saliva = 120.0f; /* 1-15+ days */
    drugs[DRUG_BARBITURATES].halflife_urine = 240.0f;  /* 2-30+ days */
    drugs[DRUG_BARBITURATES].cutoff_saliva = 50.0f;
    drugs[DRUG_BARBITURATES].cutoff_urine = 200.0f;
    drugs[DRUG_BARBITURATES].dosing_interval = 24.0f;
    strcpy(drugs[DRUG_BARBITURATES].metabolite_info, "Parent drugs + hydroxylated metabolites");

    /* BENZODIAZEPINES - Variable detection depending on specific drug */
    strcpy(drugs[DRUG_BENZODIAZEPINES].name, "BENZODIAZEPINES");
    drugs[DRUG_BENZODIAZEPINES].halflife_saliva = 72.0f; /* 1-10+ days */
    drugs[DRUG_BENZODIAZEPINES].halflife_urine = 168.0f; /* 3-30+ days */
    drugs[DRUG_BENZODIAZEPINES].cutoff_saliva = 10.0f;
    drugs[DRUG_BENZODIAZEPINES].cutoff_urine = 200.0f;
    drugs[DRUG_BENZODIAZEPINES].dosing_interval = 24.0f;
    strcpy(drugs[DRUG_BENZODIAZEPINES].metabolite_info, "Parent drugs + oxazepam + glucuronides");

    /* ALCOHOL - Short detection window */
    strcpy(drugs[DRUG_ALCOHOL].name, "ALCOHOL");
    drugs[DRUG_ALCOHOL].halflife_saliva = 1.0f;    /* 6-12 hours direct */
    drugs[DRUG_ALCOHOL].halflife_urine = 2.0f;     /* 6-24 hours direct */
    drugs[DRUG_ALCOHOL].cutoff_saliva = 25.0f;     /* 25 mg/dL */
    drugs[DRUG_ALCOHOL].cutoff_urine = 100.0f;     /* 100 mg/dL */
    drugs[DRUG_ALCOHOL].dosing_interval = 2.0f;
    strcpy(drugs[DRUG_ALCOHOL].metabolite_info, "Ethanol + EtG (up to 80 hours urine)");

    /* LSD - Very low concentrations, short window */
    strcpy(drugs[DRUG_LSD].name, "LSD");
    drugs[DRUG_LSD].halflife_saliva = 5.0f;        /* 6-24 hours */
    drugs[DRUG_LSD].halflife_urine = 8.0f;         /* 1-5 days */
    drugs[DRUG_LSD].cutoff_saliva = 0.5f;          /* Ultra-low */
    drugs[DRUG_LSD].cutoff_urine = 0.5f;
    drugs[DRUG_LSD].dosing_interval = 12.0f;
    strcpy(drugs[DRUG_LSD].metabolite_info, "Parent drug + iso-LSD + nor-LSD");

    /* KETAMINE - Dissociative anesthetic */
    strcpy(drugs[DRUG_KETAMINE].name, "KETAMINE");
    drugs[DRUG_KETAMINE].halflife_saliva = 3.5f;   /* 24-48 hours */
    drugs[DRUG_KETAMINE].halflife_urine = 14.0f;   /* 2-4 days */
    drugs[DRUG_KETAMINE].cutoff_saliva = 25.0f;
    drugs[DRUG_KETAMINE].cutoff_urine = 100.0f;
    drugs[DRUG_KETAMINE].dosing_interval = 4.0f;
    strcpy(drugs[DRUG_KETAMINE].metabolite_info, "Parent drug + norketamine + dehydronorketamine");

    /* MESCALINE - Psychedelic phenethylamine */
    strcpy(drugs[DRUG_MESCALINE].name, "MESCALINE");
    drugs[DRUG_MESCALINE].halflife_saliva = 8.0f;  /* 1-3 days */
    drugs[DRUG_MESCALINE].halflife_urine = 36.0f;  /* 2-7 days */
    drugs[DRUG_MESCALINE].cutoff_saliva = 25.0f;
    drugs[DRUG_MESCALINE].cutoff_urine = 100.0f;
    drugs[DRUG_MESCALINE].dosing_interval = 12.0f;
    strcpy(drugs[DRUG_MESCALINE].metabolite_info, "Parent drug + 3,4,5-trimethoxyphenylacetic acid");

    /* PSILOCYBIN - Detected as psilocin */
    strcpy(drugs[DRUG_PSILOCYBIN].name, "PSILOCYBIN");
    drugs[DRUG_PSILOCYBIN].halflife_saliva = 3.0f; /* 6-24 hours */
    drugs[DRUG_PSILOCYBIN].halflife_urine = 13.0f; /* 1-3 days */
    drugs[DRUG_PSILOCYBIN].cutoff_saliva = 1.0f;
    drugs[DRUG_PSILOCYBIN].cutoff_urine = 10.0f;
    drugs[DRUG_PSILOCYBIN].dosing_interval = 8.0f;
    strcpy(drugs[DRUG_PSILOCYBIN].metabolite_info, "Psilocin (active metabolite) + glucuronide");

    /* DMT - Very short detection window */
    strcpy(drugs[DRUG_DMT].name, "DMT");
    drugs[DRUG_DMT].halflife_saliva = 0.5f;        /* 15-60 minutes */
    drugs[DRUG_DMT].halflife_urine = 2.0f;         /* 2-24 hours */
    drugs[DRUG_DMT].cutoff_saliva = 1.0f;
    drugs[DRUG_DMT].cutoff_urine = 10.0f;
    drugs[DRUG_DMT].dosing_interval = 1.0f;
    strcpy(drugs[DRUG_DMT].metabolite_info, "Indole-3-acetic acid + 6-hydroxyindole-3-acetic acid");

    /* GHB - Short detection window */
    strcpy(drugs[DRUG_GHB].name, "GHB");
    drugs[DRUG_GHB].halflife_saliva = 1.0f;        /* 4-8 hours */
    drugs[DRUG_GHB].halflife_urine = 6.0f;         /* 12-24 hours */
    drugs[DRUG_GHB].cutoff_saliva = 5.0f;
    drugs[DRUG_GHB].cutoff_urine = 10.0f;
    drugs[DRUG_GHB].dosing_interval = 2.0f;
    strcpy(drugs[DRUG_GHB].metabolite_info, "Parent drug (endogenous levels present)");

    /* METHAQUALONE - Sedative-hypnotic */
    strcpy(drugs[DRUG_METHAQUALONE].name, "METHAQUALONE");
    drugs[DRUG_METHAQUALONE].halflife_saliva = 36.0f; /* 1-3 days */
    drugs[DRUG_METHAQUALONE].halflife_urine = 72.0f;  /* 7-14 days */
    drugs[DRUG_METHAQUALONE].cutoff_saliva = 25.0f;
    drugs[DRUG_METHAQUALONE].cutoff_urine = 200.0f;
    drugs[DRUG_METHAQUALONE].dosing_interval = 12.0f;
    strcpy(drugs[DRUG_METHAQUALONE].metabolite_info, "Parent drug + hydroxylated metabolites");

    /* METHADONE - Long-acting opioid agonist */
    strcpy(drugs[DRUG_METHADONE].name, "METHADONE");
    drugs[DRUG_METHADONE].halflife_saliva = 48.0f; /* 1-10+ days */
    drugs[DRUG_METHADONE].halflife_urine = 86.0f;  /* 3-14+ days */
    drugs[DRUG_METHADONE].cutoff_saliva = 25.0f;
    drugs[DRUG_METHADONE].cutoff_urine = 200.0f;
    drugs[DRUG_METHADONE].dosing_interval = 24.0f;
    strcpy(drugs[DRUG_METHADONE].metabolite_info, "Parent drug + EDDP + EMDP metabolites");

    /* PROPOXYPHENE - Synthetic opioid */
    strcpy(drugs[DRUG_DEXTROPROPOXYPHENE].name, "DEXTROPROPOXYPHENE");
    drugs[DRUG_DEXTROPROPOXYPHENE].halflife_saliva = 18.0f; /* 6-48 hours */
    drugs[DRUG_DEXTROPROPOXYPHENE].halflife_urine = 48.0f;  /* 1-2 days */
    drugs[DRUG_DEXTROPROPOXYPHENE].cutoff_saliva = 10.0f;
    drugs[DRUG_DEXTROPROPOXYPHENE].cutoff_urine = 300.0f;
    drugs[DRUG_DEXTROPROPOXYPHENE].dosing_interval = 8.0f;
    strcpy(drugs[DRUG_DEXTROPROPOXYPHENE].metabolite_info, "Parent drug + norpropoxyphene");

    /* HEROIN - Detected primarily via metabolites */
    strcpy(drugs[DRUG_DIAMORPHINE].name, "DIAMORPHINE");
    drugs[DRUG_DIAMORPHINE].halflife_saliva = 8.0f;  /* Via 6-MAM: 6-24 hours */
    drugs[DRUG_DIAMORPHINE].halflife_urine = 24.0f;  /* Via morphine: 1-4 days */
    drugs[DRUG_DIAMORPHINE].cutoff_saliva = 2.0f;    /* 6-MAM cutoff */
    drugs[DRUG_DIAMORPHINE].cutoff_urine = 10.0f;    /* 6-MAM cutoff */
    drugs[DRUG_DIAMORPHINE].dosing_interval = 4.0f;
    strcpy(drugs[DRUG_DIAMORPHINE].metabolite_info, "6-MAM (specific) + morphine + morphine glucuronides");
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
    /* Clear screen and display ASCII Art Title first */
    printf("\n\n\n");
    printf(".##....##....###....########...######..########..########.########.########..######..########\n");
    printf(".###...##...##.##...##.....##.##....##.##.....##.##..........##....##.......##....##....##...\n");
    printf(".####..##..##...##..##.....##.##.......##.....##.##..........##....##.......##..........##...\n");
    printf(".##.##.##.##.....##.########..##.......##.....##.######......##....######...##..........##...\n");
    printf(".##..####.#########.##...##...##.......##.....##.##..........##....##.......##..........##...\n");
    printf(".##...###.##.....##.##....##..##....##.##.....##.##..........##....##.......##....##....##...\n");
    printf(".##....##.##.....##.##.....##..######..########..########....##....########..######.....##...\n");
    printf("\n");
    printf("====================================================================\n");
    printf("NARCDETECT - Drug Detection Time Calculator v2.0\n");
    printf("FOR ORAL FLUID (SALIVA) TESTING\n");
    printf("====================================================================\n\n");
    printf("ESTIMATES TIME UNTIL NON-DETECTABLE\n");
    printf("BASED ON PHARMACOKINETIC PARAMETERS\n");
    printf("INCLUDES CHRONIC USE ACCUMULATION\n");
    printf("AND ROUTES OF ADMINISTRATION\n");
    printf("WITH NMR SPECTRUM SIMULATION\n\n");
    fflush(stdout);  /* Force immediate display */
}

void print_drug_menu(void)
{
    printf("AVAILABLE DRUGS BY TYPE:\n");
    printf("====================================================================\n");
    printf("SYNTHETIC OPIOIDS:     NATURAL OPIOIDS:       STIMULANTS:\n");
    printf("\n");
    printf("  FENTANYL               MORPHINE               AMPHETAMINE\n");
    printf("  NITAZENES              CODEINE                METHAMPHETAMINE\n");
    printf("  HYDROMORPHONE          DIAMORPHINE (HEROIN)   DEXTROAMPHETAMINE\n");
    printf("  OXYCODONE              PETHIDINE\n");
    printf("  HYDROCODONE            METHADONE\n");
    printf("  DEXTROPROPOXYPHENE\n");
    printf("\n");
    printf("DEPRESSANTS:           PSYCHEDELICS:          OTHER:\n");
    printf("\n");
    printf("  BARBITURATES           LSD                    ALCOHOL\n");
    printf("  BENZODIAZEPINES        MESCALINE              KETAMINE\n");
    printf("  METHAQUALONE           PSILOCYBIN             GHB\n");
    printf("                         DMT\n");
    printf("====================================================================\n\n");
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
    if (str_compare_upper(input, "HEROIN") == 0) return DRUG_DIAMORPHINE;

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

    /* Check common alternatives and abbreviations */
    /* Intravenous */
    if (str_compare_upper(input, "IV") == 0) return ROUTE_INTRAVENOUS;
    if (str_compare_upper(input, "I.V.") == 0) return ROUTE_INTRAVENOUS;
    if (str_compare_upper(input, "I.V") == 0) return ROUTE_INTRAVENOUS;
    if (str_compare_upper(input, "INJECTION") == 0) return ROUTE_INTRAVENOUS;
    
    /* Intramuscular */
    if (str_compare_upper(input, "IM") == 0) return ROUTE_INTRAMUSCULAR;
    if (str_compare_upper(input, "I.M.") == 0) return ROUTE_INTRAMUSCULAR;
    if (str_compare_upper(input, "I.M") == 0) return ROUTE_INTRAMUSCULAR;
    if (str_compare_upper(input, "MUSCLE") == 0) return ROUTE_INTRAMUSCULAR;
    
    /* Subcutaneous */
    if (str_compare_upper(input, "SC") == 0) return ROUTE_SUBCUTANEOUS;
    if (str_compare_upper(input, "SQ") == 0) return ROUTE_SUBCUTANEOUS;
    if (str_compare_upper(input, "SUBQ") == 0) return ROUTE_SUBCUTANEOUS;
    if (str_compare_upper(input, "S.C.") == 0) return ROUTE_SUBCUTANEOUS;
    if (str_compare_upper(input, "SUB-Q") == 0) return ROUTE_SUBCUTANEOUS;
    
    /* Intranasal */
    if (str_compare_upper(input, "IN") == 0) return ROUTE_INTRANASAL;
    if (str_compare_upper(input, "NASAL") == 0) return ROUTE_INTRANASAL;
    if (str_compare_upper(input, "SNORT") == 0) return ROUTE_INTRANASAL;
    if (str_compare_upper(input, "SNORTING") == 0) return ROUTE_INTRANASAL;
    if (str_compare_upper(input, "NOSE") == 0) return ROUTE_INTRANASAL;
    
    /* Inhalation */
    if (str_compare_upper(input, "INH") == 0) return ROUTE_INHALATION;
    if (str_compare_upper(input, "INHALED") == 0) return ROUTE_INHALATION;
    if (str_compare_upper(input, "SMOKING") == 0) return ROUTE_INHALATION;
    if (str_compare_upper(input, "SMOKE") == 0) return ROUTE_INHALATION;
    if (str_compare_upper(input, "VAPING") == 0) return ROUTE_INHALATION;
    if (str_compare_upper(input, "VAPE") == 0) return ROUTE_INHALATION;
    
    /* Oral */
    if (str_compare_upper(input, "PO") == 0) return ROUTE_ORAL;
    if (str_compare_upper(input, "P.O.") == 0) return ROUTE_ORAL;
    if (str_compare_upper(input, "MOUTH") == 0) return ROUTE_ORAL;
    if (str_compare_upper(input, "SWALLOW") == 0) return ROUTE_ORAL;
    if (str_compare_upper(input, "PILL") == 0) return ROUTE_ORAL;
    if (str_compare_upper(input, "TABLET") == 0) return ROUTE_ORAL;
    
    /* Sublingual */
    if (str_compare_upper(input, "SL") == 0) return ROUTE_SUBLINGUAL;
    if (str_compare_upper(input, "S.L.") == 0) return ROUTE_SUBLINGUAL;
    if (str_compare_upper(input, "UNDER TONGUE") == 0) return ROUTE_SUBLINGUAL;
    if (str_compare_upper(input, "SUB") == 0) return ROUTE_SUBLINGUAL;
    
    /* Transdermal */
    if (str_compare_upper(input, "TD") == 0) return ROUTE_TRANSDERMAL;
    if (str_compare_upper(input, "PATCH") == 0) return ROUTE_TRANSDERMAL;
    if (str_compare_upper(input, "SKIN") == 0) return ROUTE_TRANSDERMAL;
    
    /* Rectal */
    if (str_compare_upper(input, "PR") == 0) return ROUTE_RECTAL;
    if (str_compare_upper(input, "P.R.") == 0) return ROUTE_RECTAL;
    if (str_compare_upper(input, "RECTAL") == 0) return ROUTE_RECTAL;
    if (str_compare_upper(input, "SUPPOSITORY") == 0) return ROUTE_RECTAL;
    
    /* Buccal */
    if (str_compare_upper(input, "BUC") == 0) return ROUTE_BUCCAL;
    if (str_compare_upper(input, "CHEEK") == 0) return ROUTE_BUCCAL;
    
    /* Topical */
    if (str_compare_upper(input, "TOP") == 0) return ROUTE_TOPICAL;
    if (str_compare_upper(input, "CREAM") == 0) return ROUTE_TOPICAL;
    if (str_compare_upper(input, "GEL") == 0) return ROUTE_TOPICAL;

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
    if ((drug >= DRUG_HYDROMORPHONE && drug <= DRUG_PETHIDINE) || drug == DRUG_METHADONE || drug == DRUG_DIAMORPHINE) {
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
    float halflife_saliva, halflife_urine, cutoff_saliva, cutoff_urine, dosing_interval;
    float bioavail, absorpt, oral_fac;
    float single_conc_saliva, single_conc_urine, total_conc_saliva, total_conc_urine;
    float steady_conc_saliva, steady_conc_urine;
    float elim_rate_saliva, elim_rate_urine, age_factor, metab_factor;
    float accumulation_factor, r_factor_saliva, r_factor_urine;
    float detection_time_saliva, detection_time_urine, buildup_saliva, buildup_urine;
    int num_doses;
    int hours_s, minutes_s, seconds_s, days_s;
    int hours_u, minutes_u, seconds_u, days_u;

    /* Get drug parameters for both matrices */
    halflife_saliva = drugs[drug].halflife_saliva;
    halflife_urine = drugs[drug].halflife_urine;
    cutoff_saliva = drugs[drug].cutoff_saliva;
    cutoff_urine = drugs[drug].cutoff_urine;
    dosing_interval = drugs[drug].dosing_interval;

    /* Get route parameters */
    bioavail = routes[route].bioavailability;
    absorpt = routes[route].absorption_rate;
    oral_fac = routes[route].oral_factor;

    /* Apply route adjustments */
    adjust_route_parameters(drug, route, &bioavail, &oral_fac, &absorpt);

    /* Calculate single dose concentration for both matrices */
    if (drug == DRUG_FENTANYL) {
        single_conc_saliva = fentanyl_dose_constant * 1000.0f * oral_fac * bioavail / (float)weight;
        single_conc_urine = single_conc_saliva;
    } else if (drug == DRUG_ALCOHOL) {
        single_conc_saliva = (float)dosage * oral_fac * bioavail * 0.5f / (float)weight;
        single_conc_urine = single_conc_saliva;
    } else {
        single_conc_saliva = (float)dosage * oral_fac * bioavail / (float)weight;
        single_conc_urine = single_conc_saliva;
    }

    /* Adjust half-life for absorption rate (flip-flop kinetics) */
    if (absorpt > halflife_saliva * 0.693f) {
        halflife_saliva = halflife_saliva * (1.0f + absorpt / (halflife_saliva * 0.693f));
    }
    if (absorpt > halflife_urine * 0.693f) {
        halflife_urine = halflife_urine * (1.0f + absorpt / (halflife_urine * 0.693f));
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

    /* Calculate elimination rates for both matrices */
    elim_rate_saliva = 0.693f / halflife_saliva;
    elim_rate_saliva = elim_rate_saliva * age_factor * metab_factor;
    elim_rate_urine = 0.693f / halflife_urine;
    elim_rate_urine = elim_rate_urine * age_factor * metab_factor;

    /* Calculate accumulation for saliva */
    num_doses = (int)(duration / dosing_interval) + 1;
    r_factor_saliva = 1.0f - exp(-elim_rate_saliva * dosing_interval);
    r_factor_urine = 1.0f - exp(-elim_rate_urine * dosing_interval);

    if (fabs(r_factor_saliva - 1.0f) < 0.001f) {
        accumulation_factor = (float)num_doses;
    } else {
        accumulation_factor = (1.0f - pow(r_factor_saliva, (float)num_doses)) / (1.0f - r_factor_saliva);
    }

    /* Calculate concentrations for saliva (primary matrix) */
    total_conc_saliva = single_conc_saliva * accumulation_factor;
    steady_conc_saliva = single_conc_saliva / (1.0f - exp(-elim_rate_saliva * dosing_interval));
    buildup_saliva = (total_conc_saliva / steady_conc_saliva) * 100.0f;
    if (buildup_saliva > 100.0f) buildup_saliva = 100.0f;

    /* Calculate concentrations for urine */
    if (fabs(r_factor_urine - 1.0f) < 0.001f) {
        accumulation_factor = (float)num_doses;
    } else {
        accumulation_factor = (1.0f - pow(r_factor_urine, (float)num_doses)) / (1.0f - r_factor_urine);
    }
    total_conc_urine = single_conc_urine * accumulation_factor;
    steady_conc_urine = single_conc_urine / (1.0f - exp(-elim_rate_urine * dosing_interval));
    buildup_urine = (total_conc_urine / steady_conc_urine) * 100.0f;
    if (buildup_urine > 100.0f) buildup_urine = 100.0f;

    /* Calculate detection times for both matrices */
    if (total_conc_saliva > cutoff_saliva) {
        detection_time_saliva = log(total_conc_saliva / cutoff_saliva) / elim_rate_saliva;
    } else {
        detection_time_saliva = 0.0f;
    }
    
    if (total_conc_urine > cutoff_urine) {
        detection_time_urine = log(total_conc_urine / cutoff_urine) / elim_rate_urine;
    } else {
        detection_time_urine = 0.0f;
    }

    /* Plot concentration curve for saliva (primary) */
    plot_concentration_curve(total_conc_saliva, elim_rate_saliva, cutoff_saliva, halflife_saliva, duration, dosing_interval, single_conc_saliva, absorpt);

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

    printf("\nPHARMACOKINETIC DATA (SALIVA):\n");
    printf("  Half-life: %.1f hours\n", halflife_saliva);
    printf("  Cutoff: %.1f ng/mL\n", cutoff_saliva);
    printf("  Dosing interval: %.1f hours\n", dosing_interval);
    printf("  Number of doses: %d\n", num_doses);
    printf("  Single dose conc: %.2f ng/mL\n", single_conc_saliva);
    printf("  Total accum conc: %.2f ng/mL\n", total_conc_saliva);
    printf("  Elim rate: %.4f /hour\n", elim_rate_saliva);
    printf("  Steady-state conc: %.2f ng/mL\n", steady_conc_saliva);
    printf("  Buildup to SS: %.1f%%\n", buildup_saliva);

    printf("\nPHARMACOKINETIC DATA (URINE):\n");
    printf("  Half-life: %.1f hours\n", halflife_urine);
    printf("  Cutoff: %.1f ng/mL\n", cutoff_urine);
    printf("  Single dose conc: %.2f ng/mL\n", single_conc_urine);
    printf("  Total accum conc: %.2f ng/mL\n", total_conc_urine);
    printf("  Elim rate: %.4f /hour\n", elim_rate_urine);
    printf("  Steady-state conc: %.2f ng/mL\n", steady_conc_urine);
    printf("  Buildup to SS: %.1f%%\n", buildup_urine);

    /* Convert detection times to readable format */
    seconds_s = (int)(detection_time_saliva * 3600.0f);
    hours_s = seconds_s / 3600;
    minutes_s = (seconds_s - hours_s * 3600) / 60;
    seconds_s = seconds_s - hours_s * 3600 - minutes_s * 60;
    days_s = hours_s / 24;
    hours_s = hours_s - days_s * 24;

    seconds_u = (int)(detection_time_urine * 3600.0f);
    hours_u = seconds_u / 3600;
    minutes_u = (seconds_u - hours_u * 3600) / 60;
    seconds_u = seconds_u - hours_u * 3600 - minutes_u * 60;
    days_u = hours_u / 24;
    hours_u = hours_u - days_u * 24;

    printf("\nDETECTION TIME (SALIVA): %.0f seconds\n", detection_time_saliva * 3600.0f);
    printf("EQUIVALENT TO: %d hours, %d minutes, %d seconds\n", 
           (int)(detection_time_saliva * 3600.0f) / 3600, minutes_s, seconds_s);
    printf("FULL FORMAT: %d days, %d hours, %d minutes, %d seconds\n", 
           days_s, hours_s, minutes_s, seconds_s);

    printf("\nDETECTION TIME (URINE): %.0f seconds\n", detection_time_urine * 3600.0f);
    printf("EQUIVALENT TO: %d hours, %d minutes, %d seconds\n", 
           (int)(detection_time_urine * 3600.0f) / 3600, minutes_u, seconds_u);
    printf("FULL FORMAT: %d days, %d hours, %d minutes, %d seconds\n", 
           days_u, hours_u, minutes_u, seconds_u);

    printf("\nMETABOLITE INFO: %s\n", drugs[drug].metabolite_info);

    /* Disclaimers */
    printf("\n** IMPORTANT DISCLAIMERS **\n");
    printf("- Estimates based on population averages\n");
    printf("- Individual variation can be significant\n");
    printf("- Chronic use calculations are simplified\n");
    printf("- Assumes regular dosing intervals\n");
    printf("- Route-specific parameters are estimates\n");
    printf("- For research/educational use only\n");
}

void plot_concentration_curve(float c0, float kelim, float cutoff, float thalf, float duration, float dosing_interval, float single_dose_conc, float absorption_rate)
{
    float conc[61], time[61];
    char plot_line[PLOT_WIDTH + 1];
    float tmax, dt, cmax;
    float ka, absorbed_conc, eliminated_conc;
    float time_since_last_dose, dose_time;
    int i, j, pos, cutoff_pos, num_doses, dose_num;

    printf("\n====================================================================\n");
    printf("  SALIVA CONCENTRATION vs TIME WITH ACCUMULATION\n");
    printf("       (INCLUDES CHRONIC USE BUILD-UP EFFECTS)\n");
    printf("       (ADJUSTED FOR ROUTE OF ADMINISTRATION)\n");
    printf("====================================================================\n\n");

    /* Calculate absorption rate constant */
    ka = 0.693f / absorption_rate;  /* Absorption half-life to rate constant */
    if (ka < 0.1f) ka = 0.1f; /* Minimum absorption rate for IV/fast routes */
    
    /* Calculate time points - extend to show full elimination */
    tmax = max_float(duration + 8.0f * thalf, 24.0f);
    dt = tmax / 60.0f;
    num_doses = (int)(duration / dosing_interval) + 1;

    /* Calculate concentrations with realistic pharmacokinetics */
    cmax = 0.0f;
    for (i = 0; i < 61; i++) {
        time[i] = (float)i * dt;
        conc[i] = 0.0f;
        
        /* Add contribution from each dose during the dosing period */
        for (dose_num = 0; dose_num < num_doses; dose_num++) {
            dose_time = (float)dose_num * dosing_interval;
            
            if (time[i] >= dose_time) {
                time_since_last_dose = time[i] - dose_time;
                
                /* Two-compartment model with absorption and elimination */
                if (absorption_rate < 0.5f) {
                    /* Fast absorption (IV, inhalation, intranasal) */
                    absorbed_conc = single_dose_conc;
                } else {
                    /* Slower absorption with flip-flop kinetics */
                    absorbed_conc = single_dose_conc * (1.0f - exp(-ka * time_since_last_dose));
                }
                
                /* Elimination from time of absorption */
                eliminated_conc = absorbed_conc * exp(-kelim * time_since_last_dose);
                conc[i] += eliminated_conc;
            }
        }
        
        /* Continue elimination after dosing stops */
        if (time[i] > duration) {
            float time_since_end = time[i] - duration;
            conc[i] *= exp(-kelim * time_since_end);
        }
        
        /* Track maximum for scaling */
        if (conc[i] > cmax) cmax = conc[i];
    }
    
    /* Ensure reasonable scaling */
    if (cmax < cutoff * 2.0f) cmax = cutoff * 2.0f;
    if (cmax < 1.0f) cmax = 1.0f;

    printf("Time range: 0 to %.1f hours\n", tmax);
    printf("Maximum concentration: %.2f ng/mL\n", cmax);
    printf("Cutoff level: %.2f ng/mL\n", cutoff);
    printf("Dosing period: %.1f hours (%d doses)\n\n", duration, num_doses);

    /* Plot the curve */
    for (i = 0; i < 61; i++) {
        /* Clear plot line */
        for (j = 0; j < PLOT_WIDTH; j++) {
            plot_line[j] = ' ';
        }
        plot_line[PLOT_WIDTH] = '\0';

        /* Add time grid markers */
        if (i % 6 == 0) {
            for (j = 9; j < PLOT_WIDTH; j += 10) {
                if (plot_line[j] == ' ') plot_line[j] = '+';
            }
        }
        
        /* Mark end of dosing period */
        if (duration > 0 && fabs(time[i] - duration) < dt) {
            int end_pos = (int)(PLOT_WIDTH * 0.1f);
            if (end_pos >= 0 && end_pos < PLOT_WIDTH && plot_line[end_pos] == ' ') {
                plot_line[end_pos] = '|';
            }
        }

        /* Plot concentration point */
        if (conc[i] > 0.001f) {
            pos = (int)(conc[i] * (PLOT_WIDTH - 2) / cmax);
            if (pos >= 0 && pos < PLOT_WIDTH) {
                plot_line[pos] = '*';
            }
        }

        /* Plot cutoff line */
        cutoff_pos = (int)(cutoff * (PLOT_WIDTH - 2) / cmax);
        if (cutoff_pos >= 0 && cutoff_pos < PLOT_WIDTH && plot_line[cutoff_pos] != '*') {
            plot_line[cutoff_pos] = '-';
        }

        printf("%s\n", plot_line);
    }

    printf("\nLEGEND: * = CONCENTRATION CURVE\n");
    printf("        - = DETECTION CUTOFF THRESHOLD\n");
    printf("        + = TIME GRID MARKERS (every %.1f hrs)\n", tmax/10.0f);
    printf("        | = END OF DOSING PERIOD\n\n");

    /* Analysis */
    if (cmax > cutoff) {
        float detection_time = 0.0f;
        /* Find when concentration drops below cutoff */
        for (i = 60; i >= 0; i--) {
            if (conc[i] > cutoff) {
                detection_time = time[i];
                break;
            }
        }
        printf("ANALYSIS: Time to non-detection = %.1f hours (%.1f days)\n", detection_time, detection_time/24.0f);
        printf("          Peak concentration = %.2f ng/mL\n", cmax);
        printf("          Dosing duration = %.1f hours (%.1f days)\n", duration, duration/24.0f);
        printf("          Absorption rate = %.2f hours\n", absorption_rate);
        printf("          Elimination half-life = %.1f hours\n\n", thalf);
    } else {
        printf("ANALYSIS: Peak concentration (%.2f ng/mL) below cutoff\n", cmax);
        printf("          No detection expected with these parameters\n\n");
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
        case DRUG_DIAMORPHINE:
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
    int i, j, line;
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
    
    /* Create properly aligned scale for 121 characters (0-120 indices) */
    /* Major grid lines at positions: 0, 20, 40, 60, 80, 100, 120 */
    /* Each label needs to be centered above its grid marker */
    printf("12.0                10.0                8.0                 6.0                 4.0                 2.0                 0.0\n");
    printf("|                   |                   |                   |                   |                   |                   |\n");

    /* Plot spectrum (50 lines, top to bottom) */
    for (line = PLOT_HEIGHT; line >= 1; line--) {
        thresh = spec_max * (float)line / (float)PLOT_HEIGHT;

        /* Clear plot line */
        for (i = 0; i < SPECTRUM_WIDTH; i++) {
            plot_line[i] = ' ';
        }
        plot_line[SPECTRUM_WIDTH] = '\0';

        /* Add horizontal grid lines at intensity intervals */
        if (line % 10 == 0) {
            /* Major horizontal grid lines every 20% intensity (lines 10, 20, 30, 40, 50) */
            for (i = 0; i < SPECTRUM_WIDTH; i++) {
                if (plot_line[i] == ' ') {
                    plot_line[i] = '-';
                }
            }
        } else if (line % 5 == 0) {
            /* Minor horizontal grid lines every 10% intensity (lines 5, 15, 25, 35, 45) */
            for (i = 0; i < SPECTRUM_WIDTH; i++) {
                if (plot_line[i] == ' ') {
                    plot_line[i] = '.';
                }
            }
        }

        /* Plot spectrum points above threshold */
        for (i = 0; i < SPECTRUM_WIDTH; i++) {
            if (spectrum[i] >= thresh) {
                plot_line[i] = '*';
            }
        }
        
        /* Add vertical grid lines at major PPM marks (every 20 points = 2 PPM) */
        /* Grid at: 12.0, 10.0, 8.0, 6.0, 4.0, 2.0, 0.0 PPM */
        /* Positions: 0, 20, 40, 60, 80, 100, 120 */
        for (i = 0; i < SPECTRUM_WIDTH; i += 20) {
            if (i < SPECTRUM_WIDTH && plot_line[i] != '*') {
                plot_line[i] = '|';
            }
        }
        
        /* Add minor vertical grid lines at 1 PPM intervals (every 10 points) */
        for (i = 10; i < SPECTRUM_WIDTH; i += 10) {
            if ((i % 20) != 0 && plot_line[i] != '*' && plot_line[i] != '|') {
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
    printf("* = SPECTRAL PEAK    | = MAJOR PPM GRID (2 PPM)    + = MINOR PPM GRID (1 PPM)\n\n");
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
    int i;
    
    /* Copy and convert to uppercase in one step */
    for (i = 0; str1[i] && i < 99; i++) {
        s1[i] = toupper(str1[i]);
    }
    s1[i] = '\0';
    
    for (i = 0; str2[i] && i < 99; i++) {
        s2[i] = toupper(str2[i]);
    }
    s2[i] = '\0';
    
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