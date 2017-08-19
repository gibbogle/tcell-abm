#include <qstring.h>
#include "params.h"

Params::Params()
{
	PARAM_SET params[] = {
{"TC_AVIDITY_MEDIAN", 0.5, 0.1, 10.0,
"TCR avidity median parameter",
"TCR avidity has a lognormal distribution, described by the median and shape parameters. \n\
(TCR stimulation rate is proportional to the product of TC avidity and DC antigen density.)"},

{"TC_AVIDITY_SHAPE", 1.1, 1.01, 3.0,
"TCR avidity shape parameter",
"TCR avidity has a lognormal distribution, described by the median and shape parameters.\n\
The shape value must be greater than 1, and values close to 1 give distributions that are close to normal."},

{"TC_CD8_FRACTION", 0.33, 0, 0,
"T cell CD8 fraction",
"The fraction of T cells that are CD8."},

{"TC_COGNATE_FRACTION_CD4", 0.0001, 0, 0,
"CD4 T cell cognate fraction",
"The fraction of CD4 T cells that are cognate, i.e. recognize and respond to the antigen on DCs."},

{"TC_COGNATE_FRACTION_CD8", 0.0001, 0, 0,
"CD8 T cell cognate fraction",
"The fraction of CD8 T cells that are cognate, i.e. recognize and respond to the antigen on DCs."},

{"TC_STIM_RATE_CONSTANT", 1, 0, 0,
"TCR stimulation rate constant",
"Rate constant Ks for TCR stimulation, where:\n\
rate of TCR stimulation = Ks*(TCR avidity)*(DC antigen density)\n\
[molecules/min]"},
	
{"TC_STIM_HALFLIFE", 24.0, 0, 0,
"TCR stimulation halflife",
"Integrated TCR stimulation decays with a specified halflife. \n\
[hours]"},

{"DIVIDE1_MEDIAN", 6.0, 1, 10,
"1st division time median parameter",
"The time taken for the first T cell division, after full activation, has a lognormal distribution, described by the median and shape parameters. \n\
[hours]"},

{"DIVIDE1_SHAPE", 1.1, 1.01, 3.0,
"1st division time shape parameter",
"The time taken for the first T cell division, after full activation, has a lognormal distribution, described by the median and shape parameters."},

{"DIVIDE2_MEDIAN", 5.0, 1, 10,
"Later division time median parameter",
"The time taken for later T cell divisions has a lognormal distribution, described by the median and shape parameters.\n\
[hours]"},

{"DIVIDE2_SHAPE", 1.1, 1.01, 3.0,
"Later division time shape parameter",
"The time taken for later T cell divisions has a lognormal distribution, described by the median and shape parameters."},

{"MOTILITY_BETA", 0.35, 0, 0,
"Motility speed parameter",
"T cell motility is described by speed and persistence parameters, each in the range 0 - 1. \n\
 Median T cell speed is roughly proportional to MOTILITY_BETA.  Range: 0.25 - 0.5."},

{"MOTILITY_RHO", 0.76, 0, 0,
"Motility persistence parameter",
"T cell motility is described by speed and persistence parameters, each in the range 0 - 1. \n\
 MOTILITY_RHO determines the extent to which motion is in the same direction from one time step to the next.  Range: 0.5 - 0.9"},

{"DC_ANTIGEN_MEDIAN", 150, 1, 1000,
"DC antigen density median parameter",
"Antigen density is the number of pMHC per DC.  It has a lognormal distribution, described by the median and shape parameters.\n\
(TCR stimulation rate is proportional to the product of TC avidity and DC antigen density.)"},

{"DC_ANTIGEN_SHAPE", 1.15, 1.01, 4.0,
"DC antigen density shape parameter",
"Antigen density has a lognormal distribution, described by the median and shape parameters."},

{"DC_LIFETIME_MEDIAN", 10.0, 0, 0,
"DC lifetime median parameter",
"DC lifetime has a lognormal distribution, described by the median and shape parameters.\n\
[days]"},

{"DC_LIFETIME_SHAPE", 1.2, 1.01, 4.0,
"DC lifetime shape parameter",
"DC lifetime has a lognormal distribution, described by the median and shape parameters."},

{"DC_BIND_DELAY", 3.0, 0, 0,
"DC binding delay",
"After a T cell unbinds from a DC there is a delay before it can bind to a DC again.\n\
[mins]"},

{"DC_DENS_HALFLIFE", 48.0, 0, 0,
"Antigen density half-life",
"Antigen density on a DC decays with a specified half-life.\n\
[hours]"},

{"MAX_TC_BIND", 50, 0, 0,
"Max T cell binding/DC",
"The maximum number of T cells that can be in simultaneous contact with a DC."},

{"MAX_COG_BIND", 15, 0, 0,
"Max cognate T cell binding/DC",
"The maximum number of cognate T cells that can bind simultaneously to a DC."},

{"IL2_THRESHOLD", 100, 0, 0,
"Stimulation threshold for IL-2",
"Integrated TCR stimulation needed to initiate IL-2/CD5 production."},

{"ACTIVATION_THRESHOLD", 100, 0, 0,
"Stimulation threshold for activation",
"Integrated TCR stimulation level needed for full activation.\n\
 If after passing this threshold the cell fails to achieve sufficient activation for the first division it will die."},

{"FIRST_DIVISION_THRESHOLD", 300, 0, 0,
"Stimulation threshold for first division",
"Integrated TCR stimulation level needed for first division.\n\
 (Also used for UNSTAGED mode)"},

{"DIVISION_THRESHOLD", 100, 0, 0,
"Stimulation threshold for subsequent divisions",
"Integrated TCR stimulation level needed for subsequent divisions.\n\
(Also used for UNSTAGED mode)"},

{"EXIT_THRESHOLD", 100, 0, 0,
"Stimulation threshold for exit",
"Integrated TCR stimulation level below which exit is permitted (using Exit Rule #2)."},

{"STIMULATION_LIMIT", 500, 0, 0,
"Maximum stimulation level",
"Maximum integrated TCR stimulation level (saturation level).\n\
(Also used for UNSTAGED mode)"},

{"THRESHOLD_FACTOR", 1, 0, 0,
"Threshold factor",
"All stimulation threshold values are scaled by this factor."},

{"STAGED_CONTACT_RULE", 0, 0, 0,
"Contact duration rule",
"Rule for contact durations: \n\
    Constant = fixed (using specified median values) \n\
    Hill = linear variation with stimulation rate, which is a Hill function \n\
    Henrickson = lognormal distn., median depends linearly on total stimulation."},

{"USE_OVERSTIMULATION", 0, 0, 1,
"Overstimulation effect?",
"Suppressive effect of TCR overstimulation effect turned on"},

{"OVERSTIM_F1", 0.5, 0, 0,
"Overstimulation f1",
"Rule for proliferation suppression effect of excessive TCR stimulation: \n\
     S = total TCR stimulation, Sdiv1 = threshold for 1st division, P = probability that DC binding and cell division will be permanently suppressed \n\
     With f = (S - Sdiv1)/Sdiv1: \n\
     If f < f1, P = 0   if f > f2, P = 1,  else P = (f-f1)/(f2-f1)"},

{"OVERSTIM_F2", 1.0, 0, 0,
"Overstimulation f2",
"Rule for proliferation suppression effect of excessive TCR stimulation: \n\
     S = total TCR stimulation, Sdiv1 = threshold for 1st division, P = probability that DC binding and cell division will be permanently suppressed \n\
     With f = (S - Sdiv1)/Sdiv1: \n\
     If f < f1, P = 0   if f > f2, P = 1,  else P = (f-f1)/(f2-f1)"},

{"STIM_HILL_THRESHOLD", 0.01, 0, 1,
"Stimulation rate threshold",
"If (normalized avidity)*(normalized pMHC) < this threshold there is no TCR signal."},

{"STIM_HILL_N", 1, 0, 0,
"Stimulation rate Hill N",
"Stimulation rate is a Hill function of x = (normalized avidity)*(normalized pMHC), with parameters N and C. \n\
 H(x;N,c) = (1 + C^N).x^N/(x^N + C^N)"},

{"STIM_HILL_C", 0.3, 0, 0,
"Stimulation rate Hill C",
"Stimulation rate is a Hill function of x = (normalized avidity)*(normalized pMHC), with parameters N and C. \n\
 H(x;N,c) = (1 + C^N).x^N/(x^N + C^N)"},

{"ACTIVATION_MODE", 1, 0, 0,
"Activation mode",
"The activation mode is either STAGED or UNSTAGED."},

{"BINDTIME_HILL_THRESHOLD", 0.01, 0, 1,
"Signalling threshold",
"If normalized signal strength is below this threshold there will be no binding interaction,\n\
 effectively the cognate cell will behave like a non-cognate cell."},

{"BINDTIME_HILL_N", 3, 0, 0,
"Bind duration Hill N",
"The duration of the T cell-DC interaction is given by a scaled Hill function of normalized rate of TCR stimulation.\n\
 The normalized stimulation rate dS/dt is the product of normalized TCR avidity and normalized DC antigent density.\n\
 The Hill function is defined by two parameters, N and C: H(x) = (1 + C^N).x^N/(x^N + C^N) where x (=dS/dt) ranges 0-1.\n\
 The binding duration lies between the specified min and max values: Bind duration = Tmin + (Tmax - Tmin).H(x)"},

{"BINDTIME_HILL_C", 0.5, 0, 1,
"Bind duration Hill C",
"The duration of the T cell-DC interaction is given by a scaled Hill function of normalized rate of TCR stimulation.\n\
 The normalized stimulation rate dS/dt is the product of normalized TCR avidity and normalized DC antigent density.\n\
 The Hill function is defined by two parameters, N and C: H(x) = (1 + C^N).x^N/(x^N + C^N) where x (=dS/dt) ranges 0-1.\n\
 The binding duration lies between the specified min and max values: Bind duration = Tmin + (Tmax - Tmin).H(x)"},

{"BINDTIME_HILL_MIN", 10, 0, 0,
"Minimum bind duration (mins)",
"The binding duration lies between the specified min and max values.  The min value is the shortest kinapse interaction time."},

{"BINDTIME_HILL_MAX", 12, 0, 0,
"Maximum bind duration (hrs)",
"The binding duration lies between the specified min and max values.  The max value is the longest synapse interaction time."},

{"BINDTIME_ACTIVATED", 10, 0, 0,
"Activated cell bind duration (mins)",
"Activated cell bind duration."},

{"UNSTAGED_MIN_DIVIDE_T", 24, 0, 0,
"Min time to start 1st division (hrs)",
"Regardless of the level of integrated signal, this is the minimum time that must have elapsed since the cell first encountered \n\
 antigen for the first division to be initiated.  (Note that the first division time applies at this point.)"},

{"MAXIMUM_AVIDITY", 1, 0, 0,
"Avidity upper limit",
"T cell TCR avidity levels are normalized by dividing by the maximum possible value, to give values in the range 0 - 1."},

{"MAXIMUM_ANTIGEN", 250, 0, 0,
"DC antigen upper limit",
"DC antigen density levels are normalized by dividing by the maximum possible value, to give values in the range 0 - 1."},

{"K1_CD69", 0.4, 0, 0,
"CD69 K1 parameter",
"K1 parameter of the CD69 ODE (K1_CD69 = K1C, K2_CD69 = K2C, K1_S1PR1 = K1S, K2_S1PR1 = K2S)\n\
 dS/dt = normalised stimulation rate \n\
 dCD69/dt = K1C*(1-CD69)*dS/dt - K2C*CD69 \n\
 dS1PR1/dt = K1S*(1-S1PR1) - K2S*CD69*S1PR1"},

{"K2_CD69", 0.01, 0, 0,
"CD69 K2 parameter",
"K2 parameter of the CD69 ODE"},

{"K1_S1PR1", 0.01, 0, 0,
"S1PR1 K1 parameter",
"K1 parameter of the S1PR1 ODE"},

{"K2_S1PR1", 0.05, 0, 0,
"S1PR1 K2 parameter",
"K2 parameter of the S1PR1 ODE"},

{"S1PR1_EXIT_THRESHOLD", 0.6, 0, 0,
"S1PR1 exit threshold",
"S1PR1 level required for egress permission with EXIT_RULE 3"},

{"NX", 100, 0, 0,
"Lattice size",
"Dimension of the lattice (number of sites in X, Y and Z directions).  Typically 4*BLOB_RADIUS is OK."},

{"BLOB_RADIUS", 23.1, 0, 0,
"Initial blob size",
"The radius of the initial spherical blob of T cells, as number of sites.  (18.38, 23.1, 29.1, 36.7 -> 25k, 50k, 100k, 200k sites)"},

{"TC_FRACTION", 0.6, 0, 0,
"T cell fraction",
"Fraction of the paracortical volume occupied by T cells."},

{"FLUID_FRACTION", 0.1, 0, 0,
"Fluid fraction",
"Fraction of the paracortical volume occupied by fluid."},

{"DC_RADIUS", 19.0, 0, 0,
"DC SOI radius",
"Radius of DC sphere of influence.\n\
[um]"},

{"TC_TO_DC",2000.0, 0, 0,
"T cells per DC",
"Ratio of T cell count to DC count in the paracortex initially."},

{"DCrate_100k", 0.1, 0, 0,
"DC influx rate parameter",
"DC influx rate corresponding to an initial T cell count of 100k, and an inflammation level of 1.0.\n\
[/100k/min]"},

{"T_DC1", 3.5, 0, 0,
"DC constant influx duration",
"Duration of constant DC influx (T_DC1).\n\
[days]"},

{"T_DC2", 4.5, 0, 0,
"DC influx end time",
"Time of cessation of DC influx (T_DC2).  Influx reduces linearly to zero between T_DC1 and T_DC2.\n\
[days]"},

{"DC_INJECTION", 0, 0, 1,
"DCs injected?",
"DCs were injected into experimental animals"},

{"DC_INJECTION_FILE", 0, 0, 0,
"DC_injection.dat",
"File with schedule of injection of DCs."},

{"USE_TRAFFIC", 1, 0, 1,
"T cell trafficking?",
"T cell trafficking is simulated (ingress and egress)"},

{"USE_HEV_PORTALS", 1, 0, 1,
"Use HEV influx portals?",
"T cells enter through a discrete number of HEVs"},

{"USE_EXIT_CHEMOTAXIS", 0, 0, 1,
"T cell exit chemotaxis?",
"S1P-modulated T cell chemotaxis towards exit portals is simulated"},

{"USE_DC_CHEMOTAXIS", 0, 0, 1,
"T cell DC chemotaxis?",
"T cell chemotaxis towards DCs is simulated"},

{"COGNATE_ONLY", 1, 0, 1,
"Simulate only cognate cells?",
"For fast simulation, no non-cognate cells are simulated."},

{"HALVE_CD69", 0, 0, 1,
"Halve CD69 on cell division?",
"When a cell divides, the level of CD69 expression is halved."},

{"CD8_EFFECTOR_PROB", 0, 0, 0,
"Prob of effector switch on division",
"When a CD8 cell divides, it can switch to effector function with DC killing capability."},

{"USE_DESENSITISATION", 0, 0, 1,
"Use desensitisation?",
"Probability of forming a stable interaction with a DC is a function of T cell activation level."},

{"DESENS_STIM_THRESH", 0.5, 0, 0,
"Desensitisation stimulation threshold",
"Binding probability factor = 1 for S < stimulation threshold, ramping down to 0 when S = stimulation limit"},

{"DESENS_STIM_LIMIT", 1.0, 0, 0,
"Desensitisation stimulation limit",
"Binding probability factor = 1 for S < stimulation threshold, ramping down to 0 when S = stimulation limit"},

{"EXIT_RULE", 2, 0, 0,
"Cell egress control rule",
"Cognate cell egress control rules: \n\
 Generation threshold = egress permitted when cell generation exceeds a threshold \n\
 Stimulation threshold = egress permitted when cell stimulation level is below a threshold \n\
 S1PR1 threshold = egress permitted when cell S1PR1 expression level exceeds a threshold \n\
 Unlimited = egress permitted when a cell is adjacent to an exit portal"},

{"RESIDENCE_TIME_CD4", 12.0, 0, 0,
"CD4 T cell residence time",
"CD4 T cell residence time (based on no DCs).\n\
[hours]"},

{"RESIDENCE_TIME_CD8", 24.0, 0, 0,
"CD8 T cell residence time",
"CD8 T cell residence time (based on no DCs).\n\
[hours]"},

{"EXIT_PROB_CD4", 0.05, 0, 0,
"CD4 T cell exit probability",
"Exit probability for a CD4 T cell at an exit portal"},

{"EXIT_PROB_CD8", 0.05, 0, 0,
"CD8 T cell exit probability",
"Exit probability for a CD8 T cell at an exit portal"},

{"SIMULATE_PERIPHERY", 0, 0, 1,
"Simulate proliferation in the periphery?",
"Simulate proliferation in the periphery"},

{"INFLAMM_DAYS1", 3.5, 0, 0,
"Inflammation plateau duration",
"Period over which the level of inflammation signal from the periphery is constant.\n\
[days]"},

{"INFLAMM_DAYS2", 4.5, 0, 0,
"Inflammation cessation time",
"Time at which the level of inflammation signal from the periphery goes to zero.\n\
[days]"},

{"INFLAMM_LEVEL", 0.5, 0, 0,
"Inflammation level",
"The plateau inflammation signal level."},

{"CCR1_STRENGTH", 1.0, 0, 0,
"CCR1_STRENGTH",
"Strength factor for CCR1 receptor"},

{"CCR1_1", 1, 0, 0,
"CCR1_1",
"Level for naive cell"},

{"CCR1_2", 1, 0, 0,
"CCR1_2",
"Level for cell that has met antigen"},

{"CCR1_3", 1, 0, 0,
"CCR1_3",
"Level for activated cell"},

{"CCR1_SAT_THRESHOLD", 1, 0, 0,
"CCR1_SAT_THRESHOLD",
"CCL3 signal level to saturate CCR1 receptor"},

{"CCR1_REFRACTORY_TIME", 5, 0, 0,
"CCR1_REFRACTORY_TIME",
"Time for CCR1 receptor to recover sensitivity after desensitization (min)"},

{"CCL3_DIFF_COEFF", 0.01, 0, 0,
"CCL3 diffusion coeff",
"CCL3 diffusion coefficient"},

{"CCL3_HALFLIFE", 11.55, 0, 0,
"CCL3 half life",
"CCL3 half life"},

{"DCCCL3_0", 0, 0, 0,
"Use CCL3 secretion?",
"Use CCL3 secretion? (otherwise use concentration)"},

{"DCCCL3_1", 1, 0, 0,
"Use CCL3 concentration?",
"Use CCL3 concentration? (otherwise use secretion)"},

{"DC_CCL3_RATE", 1, 0, 0,
"Rate of CCL3 secretion",
"DC rate of secretion of CCL3"},

{"DC_CCL3_CONC", 1, 0, 0,
"CCL3 concentration",
"DC CCL3 concentration"},

{"DC_CCL3_RADIUS", 2.025, 0, 0,
"Radius of CCL3 conc",
"Sites within this radius of the DC receive CCL3 concentration (+ all DC sites)"},

{"NDAYS", 10.0, 0, 0,
"Number of days",
"Length of the simulation.\n\
[days]"},

{"SEED1", 1234, 0, 0,
"First RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"SEED2", 5678, 0, 0,
"Second RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"NCPU", 4, 1, 8,
"Number of CPUs",
"Number of CPUs to use for the simulation."},

{"NT_ANIMATION", 20, 0, 0,
 "Animation interval (timesteps)",
 "Interval between animation screen updates (timesteps).  One timestep = 15 sec."},

{"FACS_INTERVAL", 0, 0, 0,
"FACS plot output interval (h)",
"If > 0, this is the interval (in hours) at which FACS plots will be generated.\n\
 The file naming is: FACS_h#### where #### is the hour"},

{"SPECIES", 1, 0, 1,
"Species",
"Animal species (mouse or human)"},

{"IN_VITRO", 0, 0, 1,
"In vitro simulation",
"Select in vitro simulation, i.e. in two dimensions"},

{"IV_WELL_DIAMETER", 0.4, 0.1, 5.0,
"In vitro well diameter",
"Diameter of the well for the simulated experiment"},

{"IV_NTCELLS", 4000, 1000, 100000,
"Number of T cells",
"Initial number of T cells for the selected well size, for the simulated experiment"},

{"IV_COGNATE_FRACTION", 0.1, 0.01, 1.0,
"Fraction of cognate T cells",
"The fraction of the initial T cell population that can recognize antigen on the DCs"},

{"IV_SHOW_NONCOGNATE", 0, 0, 1,
"Display non-cognate cells",
"Display both non-cognate and cognate T cells"},

{"SPECIAL_CASE", 0, 0, 0,
"Special case simulation",
"Select one of the hard-coded special cases (> 0)"},

{"SPECIAL_CASE_FILE", 0, 0, 0,
"",
"Input file required by the selected hard-coded special case"},

{"INPUT_FILE", 0, 0, 0,
"fixed.inpdata",
"The auxiliary input file contains data that (almost!) never changes"},

{"DUMMY_IMAGE_BASENAME", 0, 0, 0,
"3D image base file name",
"The base name, including path, for saved 3D images. For example, movie/frame will give image files movie/frame00000.jpg, movie/frame00001.jpg, ... \n\
 NOTE the use of the forward slash /."},

{"DUMMY_STIMULATION_PLOT", 0, 0, 0,
"Normalized (pMHC*avidity)",
"In the STAGED case binding durations depend on stage, and the stimulation rate is a Hill function of x = (normalized pMHC)*(normalised avidity). \n\
 Note that the normalised value is the actual value divided by the specified upper limit."},

{"DUMMY_BINDTIME_PLOT", 0, 0, 0,
"Normalized (pMHC*avidity)",
"In the UNSTAGED case the stimulation rate is simply x = (normalized pMHC)*(normalised avidity), and the binding duration is a Hill function of x. \n\
 Note that the normalised value is the actual value divided by the specified upper limit."},

// Time-series plots
    {"nDC",                     1, 0,1,"","No. of DCs"},
    {"ntot_LN",                 1, 0,1,"","No of cells in the LN"},
    {"ncogseed",                1, 0,1,"",""},
    {"ncog_LN",                 1, 0,1,"",""},
    {"ncog_PER",                1, 0,1,"",""},
    {"nbnd",                    0, 0,1,"",""},
    {"nexits",                  0, 0,1,"",""},
    {"nteffgen0",               0, 0,1,"",""},
    {"nteffgen",                0, 0,1,"",""},
    {"act",                     1, 0,1,"",""},
    {"stim_LN",                 1, 0,1,"",""},
    {"stim_PER",                0, 0,1,"",""},
    {"stimrate_LN",             1, 0,1,"",""},
    {"DCcontact_time",          0, 0,1,"",""},
    {"DCtravel_time",           0, 0,1,"",""},
    {"DCbind_time",             0, 0,1,"",""},
    {"Bound_fraction",          0, 0,1,"",""},
    {"nDC_SOI",                 0, 0,1,"",""},
    {"noDC_contact",            0, 0,1,"",""},
    {"noDC_contacttime",        0, 0,1,"",""},
    {"totDC_contacttime_LN",    0, 0,1,"",""},
    {"totDC_contacttime_PER",   0, 0,1,"",""},
// Profile plots
    {"CD69",                    0, 0,1,"",""},
    {"S1PR1",                   1, 0,1,"",""},
    {"Stimulation",             1, 0,1,"",""},
    {"Stimulation rate",        1, 0,1,"",""},
    {"Avidity LN",              1, 0,1,"",""},
    {"Avidity PER",             0, 0,1,"",""},
    {"Generation LN",           0, 0,1,"",""},
    {"DC contact time (min)",   0, 0,1,"",""},
    {"DC bind time (min)",      0, 0,1,"",""}


};
    nParams = sizeof(params)/sizeof(PARAM_SET);
    workingParameterList = new PARAM_SET[nParams];
    for (int i=0; i<nParams; i++) {
        workingParameterList[i] = params[i];
    }
}


PARAM_SET Params::get_param(int k)
{
	return workingParameterList[k];
}

void Params::set_value(int k, double v)
{
	workingParameterList[k].value = v;
}

void Params::set_label(int k, QString str)
{
	workingParameterList[k].label = str;
}
