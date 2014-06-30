//ABM
#include <qstring.h>
#include "graphs.h"
#include "log.h"

LOG_USE();

Graphs::Graphs()
{
    GRAPH_SET tsGraphSet[] = {

//        summaryData(1:22) = (/ int(tnow/60), istep, NDCalive, ntot_LN, nseed, ncog(1), ncog(2), ndead, &
//            nbnd, int(InflowTotal), Nexits, nteffgen0, nteffgen, nact, navestim(1), navestim(2), navestimrate(1), &
//            naveDCtime, naveDCtraveltime, naveDCbindtime, nbndfraction, nDCSOI /)

//{"dummy",
//"",
//"",
//0, false, 0, 1, 0, true},

{"nDC",
"Antigen Presenting Cells",
"",  //No. of cells ",
2, true, 0, 1, 0, true},

{"ntot_LN",
"Total LN Population",
"",  //No. of cells",
3, true, 0, 1, 0, true},

{"ncogseed",
"Seed Cognate Cells",
"",  //No. of cells",
4, true, 0, 1, 0, true},

{"ncog_LN",
"Cognate T Cells (LN)",
"",  //No. of cells",
5, true, 0, 1, 0, true},

{"ncog_PER",
"Activated Cells (PER)",
"",  //No. of cells",
6, true, 0, 1, 0, true},

{"nbnd",
"Bound Cognate Cells",
"",  //No. of cells",
8, false, 0, 1, 0, true},

{"nexits",
"No. of Exit Portals",
"",  //No. of portals",
10, false, 0, 1, 0, true},

{"nteffgen0",
"Efferent Unactivated",
"",  //No. of cells",
11, false, 0, 1, 0, true},

{"nteffgen",
"Efferent Activated",
"",  //No. of cells",
12, true, 0, 1, 0, true},

{"act",
"Total DC Antigen",
"",
13, true, 0, .001, 0, true},

{"stim_LN",
"Average Stimulation (LN)",
"",
14, true, 0, .001, 1.0, true},

{"stim_PER",
"Average Stimulation (PER)",
"",
15, false, 0, .001, 1.0, true},

{"stimrate_LN",
"Stimulation Rate (LN)",
"",
16, false, 0, .001, 1.0, true},

{"DCcontact_time",
"1st DC Contact (min)",
"",
17, false, 0, .001, 0, true},

{"DCtravel_time",
"Inter-DC Travel (min)",
"",
18, false, 0, .001, 0, true},

{"DCbind_time",
"DC Bind Time (min)",
"",
19, false, 0, .001, 0, true},

{"Bound_fraction",
"Bound Fraction",
"",
20, true, 0, .001, 1.0, true},

{"nDC_SOI",
"Cells in DC SOI",
"",
21, false, 0, .001, 0, true},

{"CD69",
"CD69 Profile",
"",  //Fraction",
PROFILE_CD69, false, 0, 1, 1.0, false},

{"S1PR1",
"S1PR1 Profile",
"",  //Fraction",
PROFILE_S1PR1, true, 0, 1, 1.0, false},

{"Stimulation",
"Stimulation Profile (LN)",
"",  //Fraction",
PROFILE_STIM, true, 0, 1, 1.0, false},

{"Stimulation Rate",
"Stim Rate Profile (LN)",
"",  //Fraction",
PROFILE_STIMRATE, false, 0, 1, 1.0, false},

{"Avidity LN",
"Avidity Profile (LN)",
"",  //Fraction",
PROFILE_AVIDITY_LN, true, 0, 1, 1.0, false},

{"Avidity PER",
"Avidity Profile (PER)",
"",  //Fraction",
PROFILE_AVIDITY_PER, false, 0, 1, 1.0, false},

{"Generation LN",
"Generation Profile (LN)",
"",  //Fraction",
PROFILE_GENERATION_LN, false, 0, 1, 1.0, false},

{"DC Contact Time (min)",
"First DC Contact Time Profile",
"",  //Fraction",
PROFILE_FIRSTDCCONTACTTIME, false, 0, 1, 1.0, false},

{"DC Bind Time (min)",
"DC Bind Time Profile",
"",  //Fraction",
PROFILE_DCBINDTIME, false, 0, 1, 1.0, false}

};
    // Note: tsGraphs[] is constant = tsGraphSet[]
    // .active is either true or false, depending on checkBox_selectgraph
    // Note: the "ts" prefix is misleading because the list can include profile plots.

    n_tsGraphs = sizeof(tsGraphSet)/sizeof(GRAPH_SET);
    tsGraphs = new GRAPH_SET[n_tsGraphs];
    for (int i=0; i<n_tsGraphs; i++) {
        tsGraphs[i] = tsGraphSet[i];
    }
    graphList = new GRAPH_SET[maxGraphs];
    nGraphs = maxGraphs;
}


GRAPH_SET Graphs::get_graph(int k)
{
	return graphList[k];
}

int Graphs::get_dataIndex(int k)
{
	return graphList[k].dataIndex;
}

QString Graphs::get_tag(int k)
{
	return graphList[k].tag;
}

QString Graphs::get_title(int k)
{
//    QString title;
//    if (!graphList[k].active)
//        title = "";
//    else
//        title = graphList[k].title;
    return graphList[k].title;
}

QString Graphs::get_yAxisTitle(int k)
{
//    QString title;
//    if (!graphList[k].active)
//        title = "";
//    else
//        title = graphList[k].yAxisTitle;
    return graphList[k].yAxisTitle;
}

double Graphs::get_maxValue(int k) {
	return graphList[k].maxValue;
}

double Graphs::get_scaling(int k) {
	return graphList[k].scaling;
}

double Graphs::get_yscale(int k) {
    return graphList[k].yscale;
}

double Graphs::get_xscale(double xmax) {
    int n = 1;
    for (;;) {
        if (xmax <= n) break;
        n++;
    }
    return double(n);
}

bool Graphs::isActive(int k)
{
	return graphList[k].active;
}

bool Graphs::isTimeseries(int k)
{
    return graphList[k].ts;
}

bool Graphs::isProfile(int k)
{
    return !graphList[k].ts;
}

void Graphs::set_maxValue(int k, double v)
{
	graphList[k].maxValue = v;
}

void Graphs::makeGraphList()
{
    char msg[128];
    int k = maxGraphs;
    int nts = 0;
    for (int i=0; i<n_tsGraphs; i++) {
        if (tsGraphs[i].active) {
            k--;
            graphList[k] = tsGraphs[i];
            nts++;
//            if (nts == maxGraphs - non_ts) break;
            if (nts == maxGraphs) break;
        }
    }
//    int ndummy = maxGraphs - nts - non_ts;
    int ndummy = maxGraphs - nts;
//    sprintf(msg,"nts: %d  ndummy: %d",nts,ndummy);
//    LOG_MSG(msg);
    for (k=0; k<ndummy; k++) {
        graphList[k].active = false;
        graphList[k].ts = true;
        graphList[k].tag = "dummy";
        graphList[k].scaling = 1;
        graphList[k].title = "";
        graphList[k].yAxisTitle = "";
    }
//    for (k=ndummy; k<ndummy + non_ts; k++) {
//        graphList[k].tag = "non_ts";
//        graphList[k].active = true;
//        graphList[k].ts = false;
//        graphList[k].scaling = 1;
//    }
    nGraphs = maxGraphs;
//    sprintf(msg,"nGraphs: %d",nGraphs);
//    LOG_MSG(msg);
//    for (k=0; k<nGraphs; k++) {
//        LOG_QMSG(graphList[k].tag);
//        sprintf(msg,"k: %d  active: %d  ts: %d",k,graphList[k].active,graphList[k].ts);
//        LOG_MSG(msg);
//    }
}
