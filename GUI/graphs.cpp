//ABM
#include <qstring.h>
#include "graphs.h"
#include "log.h"

LOG_USE();

Graphs::Graphs()
{
    GRAPH_SET tsGraphSet[] = {

//{"dummy",
//"",
//"",
//0, false, 0, 1, 0, true},

{"ntot_LN",
"Total T Cell Population in LN",
"No. of cells",
4, true, 0, 1, 0, true},

{"ncog_LN",
"Cognate T Cells in LN",
"No. of cells",
6, true, 0, 1, 0, true},

{"ncogseed",
"Seed Cognate Cells",
"No. of cells",
5, true, 0, 1, 0, true},

{"nbnd",
"Bound Cognate Cells",
"No. of cells",
9, true, 0, 1, 0, true},

{"teffgen",
"Efferent Activated Cells",
"No. of cells",
12, true, 0, 1, 0, true},

{"nDC",
"Antigen Presenting Cells",
"No. of cells",
2, true, 0, 1, 0, true},

{"act",
"Total DC Antigen Activity",
"",
3, true, 0, .01, 0, true},

{"stim",
"Average Stimulation",
"",
13, true, 0, .01, 1.0, true},

{"stimrate",
"Average Stimulation Rate",
"",
14, true, 0, .01, 1.0, true},

{"CD69",
"CD69 Profile",
"Fraction",
PROFILE_CD69, true, 0, 1, 1.0, false},

{"S1PR1",
"S1PR1 Profile",
"Fraction",
PROFILE_S1PR1, true, 0, 1, 1.0, false},

{"Stimulation",
"Stimulation Profile",
"Fraction",
PROFILE_STIM, true, 0, 1, 1.0, false},

{"Stimulation Rate",
"Stimulation Rate Profile",
"Fraction",
PROFILE_STIMRATE, true, 0, 1, 1.0, false},

{"ncog_PER",
"Activated T Cells in Periphery",
"No. of cells",
7, false, 0, 1, 0, true},

{"nexits",
"Exit Portals",
"No. of portals",
11, false, 0, 1, 0, true}

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
