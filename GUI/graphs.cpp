#include <qstring.h>
#include "graphs.h"

Graphs::Graphs()
{
	GRAPH_SET graphs[] = {

{"dummy",
"",
"",
1, false, 0, 1},

{"teffgen",
"Efferent Activated Cognate Cells",
"No. of cells",
8, true, 0, 1},

{"ncogseed",
"Seed Cognate Cells",
"No. of cells",
4, true, 0, 1},

{"nbnd",
"Bound Cognate Cells",
"No. of cells",
9, true, 0, 1},

{"act",
"Total DC antigen activity",
"",
2, true, 0, .01},

{"nDC",
"Antigen Presenting Cells",
"No. of cells",
1, true, 0, 1},

{"ncog_PER",
"Cognate T Cells in periphery",
"No. of cells",
6, true, 0, 1},

{"ncog_LN",
"Cognate T Cells in LN",
"No. of cells",
5, true, 0, 1},

{"ntot_LN",
"Total T Cell Population in LN",
"No. of cells",
3, true, 0, 1}

};
	nGraphs = sizeof(graphs)/sizeof(GRAPH_SET);
	graphList = new GRAPH_SET[nGraphs];
	for (int i=0; i<nGraphs; i++) {
		graphList[i] = graphs[i];
	}
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
	return graphList[k].title;
}

QString Graphs::get_yAxisTitle(int k)
{
	return graphList[k].yAxisTitle;
}

double Graphs::get_maxValue(int k) {
	return graphList[k].maxValue;
}

double Graphs::get_scaling(int k) {
	return graphList[k].scaling;
}

bool Graphs::isActive(int k)
{
	return graphList[k].active;
}

void Graphs::set_maxValue(int k, double v)
{
	graphList[k].maxValue = v;
}
