#include <qstring.h>
#include "params.h"

Params::Params()
{
	PARAM_SET params[] = {
{"TC_AVIDITY_MEDIAN", 1.0, 0.1, 10.0,
"TCR avidity median parameter",
"TCR avidity has a lognormal distribution, described by the median and shape parameters.\n\
(TCR stimulation rate is proportional to the product of TC avidity and DC antigen density.)"},

{"TC_AVIDITY_SHAPE", 1.1, 1.01, 3.0,
"TCR avidity shape parameter",
"TCR avidity has a lognormal distribution, described by the median and shape parameters.\n\
The shape value must be greater than 1, and values close to 1 give distributions that are close to normal."},

{"TC_COGNATE_FRACTION", 0.0001, 0.0, 0.002,
"T cell cognate fraction",
"The fraction of T cells that are cognate, i.e. recognize and respond to the antigen on DCs."},

{"TC_STIM_RATE_CONSTANT", 1, 0.0, 0.0,
"TCR stimulation rate constant",
"Rate constant Ks for TCR stimulation, where:\n\
rate of TCR stimulation = Ks*(TCR avidity)*(DC antigen density)\n\
[molecules/min]"},
	
{"TC_STIM_HALFLIFE", 24.0, 0.0, 100.0,
"TCR stimulation halflife",
"Integrated TCR stimulation decays with a specified halflife. \n\
[hours]"},

{"DIVIDE1_MEDIAN", 6.0, 0.0, 0.0,
"1st division time median parameter",
"The time taken for the first T cell division, after full activation, has a lognormal distribution, described by the median and shape parameters. \n\
[hours]"},

{"DIVIDE1_SHAPE", 1.1, 1.01, 3.0,
"1st division time shape parameter",
"The time taken for the first T cell division, after full activation, has a lognormal distribution, described by the median and shape parameters."},

{"DIVIDE2_MEDIAN", 5.0, 0, .0,
"Later division time median parameter",
"The time taken for later T cell divisions has a lognormal distribution, described by the median and shape parameters.\n\
[hours]"},

{"DIVIDE2_SHAPE", 1.1, 1.01, 3.0,
"Later division time shape parameter",
"The time taken for later T cell divisions has a lognormal distribution, described by the median and shape parameters."},

{"MOTILITY_BETA", 0.35, 0.25, 0.5,
"Motility speed parameter",
"T cell motility is described by speed and persistence parameters, each in the range 0 - 1. Median T cell speed is roughly proportional to MOTILITY_BETA."},

{"MOTILITY_RHO", 0.76, 0.5, 0.9,
"Motility persistence parameter",
"T cell motility is described by speed and persistence parameters, each in the range 0 - 1. MOTILITY_RHO determines the extent to which motion is in the same direction from one time step to the next."},

{"DC_ANTIGEN_MEDIAN", 150, 50.0, 1000.0,
"DC antigen density median parameter",
"Antigen density is the number of pMHC per DC.  It has a lognormal distribution, described by the median and shape parameters.\n\
(TCR stimulation rate is proportional to the product of TC avidity and DC antigen density.)"},

{"DC_ANTIGEN_SHAPE", 1.5, 1.01, 3.0,
"DC antigen density shape parameter",
"Antigen density has a lognormal distribution, described by the median and shape parameters."},

{"DC_LIFETIME_MEDIAN", 3.0, 1.0, 20.0,
"DC lifetime median parameter",
"DC lifetime has a lognormal distribution, described by the median and shape parameters.\n\
[days]"},

{"DC_LIFETIME_SHAPE", 1.2, 1.01, 3.0,
"DC lifetime shape parameter",
"DC lifetime has a lognormal distribution, described by the median and shape parameters."},

{"DC_BIND_DELAY", 3.0, 0.0, 20.0,
"DC binding delay",
"After a T cell unbinds from a DC there is a delay before it can bind to a DC again.\n\
[mins]"},

{"DC_DENS_HALFLIFE", 6.01, 0.1, 100.0,
"Antigen density half-life",
"Antigen density on a DC decays with a specified half-life.\n\
[hours]"},

{"MAX_TC_BIND", 100, 50, 200,
"Max T cell binding/DC",
"The maximum number of T cells that can be in simultaneous contact with a DC."},

{"MAX_COG_BIND", 15, 5, 20,
"Max cognate T cell binding/DC",
"The maximum number of cognate T cells that can bind simultaneously to a DC."},

{"IL2_THRESHOLD", 150, 10.0, 500.0,
"Stimulation threshold for IL-2",
"Integrated TCR stimulation needed to initiate IL-2/CD5 production."},

{"ACTIVATION_THRESHOLD", 150, 10.0, 500.0,
"Stimulation threshold for activation",
"Integrated TCR stimulation level needed for full activation."},

{"FIRST_DIVISION_THRESHOLD", 300, 10.0, 1000.0,
"Stimulation threshold for first division",
"Integrated TCR stimulation level needed for first division."},

{"DIVISION_THRESHOLD", 80, 10.0, 1000.0,
"Stimulation threshold for subsequent divisions",
"Integrated TCR stimulation level needed for subsequent divisions."},

{"EXIT_THRESHOLD", 480, 10.0, 1000.0,
"Stimulation threshold for exit",
"Integrated TCR stimulation level below which exit is permitted (using Exit Rule #2)."},

{"STIMULATION_LIMIT", 1000, 0.0, 0.0,
"Maximum stimulation level",
"Maximum integrated TCR stimulation level (saturation level)."},

{"NX", 100, 100, 300,
"Lattice size",
"Dimension of the lattice (number of sites in X, Y and Z directions).  Typically 4*BLOB_RADIUS is OK."},

{"BLOB_RADIUS", 23.1, 15.0, 50.0,
"Initial blob size",
"The radius of the initial spherical blob of T cells, as number of sites.  (18.38, 23.1, 29.1, 36.7 -> 25k, 50k, 100k, 200k sites)"},

{"TC_FRACTION", 0.6, 0.4, 0.8,
"T cell fraction",
"Fraction of the paracortical volume occupied by T cells."},

{"FLUID_FRACTION", 0.1, 0.05, 0.2,
"Fluid fraction",
"Fraction of the paracortical volume occupied by fluid."},

{"DC_RADIUS", 19.0, 10.0, 30.0,
"DC radius",
"Radius of DC sphere of influence.\n\
[um]"},

{"TC_TO_DC",2000.0, 50, 10000,
"T cells per DC",
"Ratio of T cell count to DC count in the paracortex initially."},

{"DCrate_100k", 0.1, 0.0, 5.0,
"DC influx rate parameter",
"DC influx rate corresponding to an initial T cell count of 100k, and an inflammation level of 1.0.\n\
[/100k/min]"},

{"T_DC1", 3.5, 0.0, 10.0,
"DC constant influx duration",
"Duration of constant DC influx (T_DC1).\n\
[days]"},

{"T_DC2", 4.5, 0.0, 10.0,
"DC influx end time",
"Time of cessation of DC influx (T_DC2).  Influx reduces linearly to zero between T_DC1 and T_DC2.\n\
[days]"},

{"DC_INJECTION", 0, 0, 1,
"DCs injected?",
"DCs were injected into experimental animals"},

{"T_DC_INJECTION", 0, 0.0, 0.0,
"DC injection time",
"Time of injection of DCs (T_DC_INJECTION).  In an experiment simulation antigen density begins to decay at this time.\n\
[hours]"},

{"USE_TRAFFIC", 1, 0, 1,
"T cell trafficking?",
"T cell trafficking is simulated (ingress and egress)"},

{"RESIDENCE_TIME", 24.0, 12.0, 36.0,
"T cell residence time",
"T cell residence time (based on no DCs).\n\
[hours]"},

{"INFLAMM_DAYS1", 3.5, 0.0, 10.0,
"Inflammation plateau duration",
"Period over which the level of inflammation signal from the periphery is constant.\n\
[days]"},

{"INFLAMM_DAYS2", 4.5, 0.0, 10.0,
"Inflammation cessation time",
"Time at which the level of inflammation signal from the periphery goes to zero.\n\
[days]"},

{"INFLAMM_LEVEL", 1.0, 0.0, 10.0,
"Inflammation level",
"The plateau inflammation signal level."},

{"EXIT_RULE", 3, 1, 3,
"Exit rule",
"T cell exit rule.  1 = use NGEN_EXIT, 2 = use EXIT_THRESHOLD, 3 = no restriction."},

{"EXIT_REGION", 1, 1, 3,
"Exit region",
"Determines blob region for cell exits: 1 = everywhere, 2 = lower half of blob, 3 = by chemotaxis, via discrete exits."},

{"CHEMO_RADIUS", 30.0, 10.0, 200.0,
"Radius of chemotactic influence",
"Range of chemotactic influence of an exit site or DC on T cell motion.  At this distance the influence is reduced to 5% of its maximum value.\n\
[um]"},

{"CHEMO_K_EXIT", 0.5, 0.0, 1.0,
"Exit chemotaxis influence parameter",
"Strength of chemotactic influence on T cell motion towards exits."},

{"CHEMO_K_DC", 0.0, 0.0, 10.0,
"DC chemotaxis influence parameter",
"Strength of chemotactic influence on T cell motion towards DCs."},

{"NDAYS", 1.0, 0.0, 30.0,
"Number of days",
"Length of the simulation.\n\
[days]"},

{"SEED1", 1234, 0, 0,
"First RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"SEED2", 5678, 0, 0,
"Second RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"NCPU", 1, 1, 8,
"Number of CPUs",
"Number of CPUs to use for the simulation."},

{"NT_ANIMATION", 20, 0, 0,
 "Animation interval (timesteps)",
 "Interval between animation screen updates (timesteps).  One timestep = 15 sec."},

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
"Display both non-cognate and cognate T cells"}

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
