/****************************************************************************
 ABM_GUI
****************************************************************************/
//ABM

#include <QtGui>

#include "mainwindow.h"
#include "log.h"
#include "params.h"
#include "graphs.h"
#include "misc.h"
#include "plot.h"
#include "myvtk.h"
#include "transfer.h"

#ifdef _WIN32
#include "windows.h"
#define sleep(n) Sleep(1000 * n)
#endif

#ifdef __LINUX
#include <QTcpServer>
#else
#include <QTcpServer.h>
#endif

LOG_USE();

Params *parm;	// I don't believe this is the right way, but it works!
Graphs *grph;

bool showingVTK;
bool recordingVTK;
bool showingFACS;
bool recordingFACS;

int recordfrom, recordto;
int VTKbuffer[100];
//int TC_list[5*MAX_TC];
int *TC_list=NULL;
int nTC_list;
//int DC_list[5*MAX_DC];
int *DC_list=NULL;
int nDC_list;
//int bond_list[2*MAX_BOND];
int *bond_list=NULL;
int nbond_list;
double *FACS_data=NULL;
int nFACS_vars=6;   // CFSE, CD69, CD8, S1PR1, avidity, stimulation
int nFACS_cells=0;
int nFACS_dim=0;

QMutex mutex1, mutex2;

bool redimflag;

int summaryData[100];
int NX, NY, NZ, NBY;
int nt_vtk;
bool leftb;

double *profile_x[20];
double *profile_y[20];
int profile_n[20];

bool display_plots = true;


/*
QMyLabel::QMyLabel(QWidget *parent) : QLabel(parent)
{}
//--------------------------------------------------------------------------------------------------------
// Redefines mousePressEvent for QMyLabel, which extends QLabel.  This is used to display info about
// a model parameter.
//--------------------------------------------------------------------------------------------------------
void QMyLabel::mousePressEvent (QMouseEvent *event) {
	event->accept();
	QString sname = objectName().mid(6);
	QString text = "mousePressEvent";
	// Find which label_ sent the signal, and read its text
	for (int k=0; k<parm->nParams; k++) {
		PARAM_SET param = parm->get_param(k);
		if (sname.compare(param.tag) == 0)
			text = param.text;
	}
    emit labelClicked(text);
}
*/

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
MainWindow::MainWindow(QWidget *parent)
   : QMainWindow(parent)
{
	LOG_MSG("Started MainWindow");
	setupUi(this);
    showMaximized();
    // Some initializations
    nDistPts = 200;
	nTicks = 1000;
	tickVTK = 100;	// timer tick for VTK in milliseconds
	paramSaved = false;
	paused = false;
	posdata = false;
    DCmotion = false;
    done = false;
    first = true;
	started = false;
    firstVTK = true;
    recordingVTK = false;
    showingVTK = false;
    recordingFACS = false;
    showingFACS = false;
    nGraphCases = 0;
	for (int i=0; i<Plot::ncmax; i++) {
		graphResultSet[i] = 0;
	}
    for (int i=0; i<20; i++) {  // need to fix the hard-coded numbers !!!!!!!!
        profile_x[i] = (double *)malloc(1000*sizeof(double));
        profile_y[i] = (double *)malloc(1000*sizeof(double));
    }
	vtkfile = "basecase.pos";
	savepos_start = 0;
	ntimes = 0;
	hour = 0;

	param_to_sliderIndex = NULL;
	defaultInputFile = "basecase.inp";
	inputFile = defaultInputFile;

	parm = new Params();
	nParams = parm->nParams;
	grph = new Graphs();
    setupGraphSelector();
    setGraphsActive();

    for (int i=0; i<MAX_DATA; i++)
        pGraph[i] = NULL;
    LOG_QMSG("did Graphs");

	createLists();
    LOG_QMSG("did createLists");
    createActions();
    LOG_QMSG("did createActions");
    initDistPlots();
    LOG_QMSG("did initDistPlots");
    initFACSPlot();
    LOG_QMSG("did initFACSPlot");
    loadParams();
    LOG_QMSG("did loadParams");
    writeout();
    timer = new QTimer(this);
    QRect rect;
    rect = groupBox_run->geometry();
#ifdef __DISPLAY768
    rect.setHeight(480);
#else
    rect.setHeight(600);
#endif
    groupBox_run->setGeometry(rect);

    vtk = new MyVTK(mdiArea_VTK, test_page);
    vtk->init();

    rect.setX(50);
    rect.setY(30);
#ifdef __DISPLAY768
    rect.setHeight(642);
    rect.setWidth(642);
#else
    rect.setHeight(786);
    rect.setWidth(786);
#endif
    mdiArea_VTK->setGeometry(rect);
    /*
    rect.setX(10);
    rect.setY(0);
#ifdef __DISPLAY768
    rect.setHeight(500);
    rect.setWidth(1100);
#else
    rect.setHeight(700);
    rect.setWidth(1500);
#endif
    mdiArea->setGeometry(rect);
    */
    showmdiAreaSize();
    tabs->setCurrentIndex(0);

//Testing
    TC_list = (int *)malloc(5*MAX_TC*sizeof(int));
    DC_list = (int *)malloc(5*MAX_DC*sizeof(int));
    bond_list = (int *)malloc(2*MAX_BOND*sizeof(int));

    videoVTK = new QVideoOutput(this, VTK_SOURCE, vtk->renWin, NULL);
    videoFACS = new QVideoOutput(this, QWT_SOURCE, NULL, qpFACS);

//    sthread0 = NULL;
//    sthread1 = NULL;
//    exthread = NULL;
//    newR = NULL;
//    checkMemory("start");
    goToInputs();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::checkMemory(QString loc)
{
    QString system_info = loc + ": ";
    MEMORYSTATUSEX memory_status;
    ZeroMemory(&memory_status, sizeof(MEMORYSTATUSEX));
    memory_status.dwLength = sizeof(MEMORYSTATUSEX);
    if (GlobalMemoryStatusEx(&memory_status)) {
      system_info.append(
            QString("RAM: %1 MB")
            .arg(memory_status.ullTotalPhys / (1024 * 1024)));
    } else {
      system_info.append("Unknown RAM");
    }
    LOG_QMSG(system_info);
    QMessageBox msgBox;
    msgBox.setText(system_info);
    msgBox.exec();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::createActions()
{
	action_stop->setEnabled(false);
    action_pause->setEnabled(false);
    action_inputs->setEnabled(false);
    action_outputs->setEnabled(false);
	action_save_snapshot->setEnabled(false);
	text_more->setEnabled(false);
    connect(action_open_input, SIGNAL(triggered()), this, SLOT(readInputFile()));
    connect(action_load_results, SIGNAL(triggered()), this, SLOT(loadResultFile()));
    connect(action_saveAs, SIGNAL(triggered()), this, SLOT(saveAs()));
    connect(action_save, SIGNAL(triggered()), this, SLOT(save()));
    connect(action_inputs, SIGNAL(triggered()), SLOT(goToInputs()));
    connect(action_outputs, SIGNAL(triggered()), SLOT(goToOutputs()));
	connect(action_VTK, SIGNAL(triggered()), SLOT(goToVTK()));
    connect(action_FACS, SIGNAL(triggered()), SLOT(goToFACS()));
    connect(action_run, SIGNAL(triggered()), SLOT(runServer()));
    connect(action_pause, SIGNAL(triggered()), SLOT(pauseServer()));
    connect(action_stop, SIGNAL(triggered()), SLOT(stopServer()));
    connect(action_play_VTK, SIGNAL(triggered()), SLOT(playVTK()));
    connect(action_set_speed, SIGNAL(triggered()), SLOT(setVTKSpeed()));
    connect(actionShow_2D_gradient_field, SIGNAL(triggered()), this, SLOT(on_action_show_gradient2D_triggered()));
    connect(actionShow_3D_gradient_field, SIGNAL(triggered()), this, SLOT(on_action_show_gradient3D_triggered()));
    connect(buttonGroup_STAGED_CONTACT_RULE, SIGNAL(buttonClicked(QAbstractButton*)), this, SLOT(radioButtonChanged(QAbstractButton*)));
    connect(buttonGroup_ACTIVATION_MODE, SIGNAL(buttonClicked(QAbstractButton*)), this, SLOT(radioButtonChanged(QAbstractButton*)));
    connect(buttonGroup_EXIT_RULE, SIGNAL(buttonClicked(QAbstractButton*)), this, SLOT(radioButtonChanged(QAbstractButton*)));
    connect(buttonGroup_FACS_PLOT, SIGNAL(buttonClicked(QAbstractButton*)), this, SIGNAL(facs_update()));
    for (int i=0; i<nLabels; i++) {
		QLabel *label = label_list[i];
		QString label_str = label->objectName();
        if (label_str.startsWith("label_") && label->inherits("QMyLabel")) {
			connect((QObject *)label, SIGNAL(labelClicked(QString)), this, SLOT(showMore(QString)));
        }
	}
	// Graph menu
    connect(action_add_graph, SIGNAL(triggered()), this, SLOT(addGraph()));
    connect(action_remove_graph, SIGNAL(triggered()), this, SLOT(removeGraph()));
    connect(action_remove_all, SIGNAL(triggered()), this, SLOT(removeAllGraphs()));
    connect(action_save_snapshot, SIGNAL(triggered()), this, SLOT(saveSnapshot()));
    connect(actionStart_recording_VTK, SIGNAL(triggered()), this, SLOT(startRecorderVTK()));
    connect(actionStop_recording_VTK, SIGNAL(triggered()), this, SLOT(stopRecorderVTK()));
    connect(actionStart_recording_FACS, SIGNAL(triggered()), this, SLOT(startRecorderFACS()));
    connect(actionStop_recording_FACS, SIGNAL(triggered()), this, SLOT(stopRecorderFACS()));
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::createLists()
{        
	lineEdit_list = findChildren<QLineEdit *>();
	spin_list = findChildren<QSpinBox *>();
	combo_list = findChildren<QComboBox *>();
	checkbox_list = findChildren<QCheckBox *>();
	radiobutton_list = findChildren<QRadioButton *>();
	slider_list = findChildren<QSlider *>();
	label_list = findChildren<QLabel *>();

	for (int i=0; i<lineEdit_list.length(); i++) {
		widget_list.append(lineEdit_list[i]);
	}
	for (int i=0; i<spin_list.length(); i++) {
		widget_list.append(spin_list[i]);
	}
	for (int i=0; i<combo_list.length(); i++) {
		widget_list.append(combo_list[i]);
	}
	for (int i=0; i<checkbox_list.length(); i++) {
		widget_list.append(checkbox_list[i]);
	}
	for (int i=0; i<radiobutton_list.length(); i++) {
		widget_list.append(radiobutton_list[i]);
	}

	nWidgets = widget_list.length();
	nSliders = slider_list.length();
	nLabels = label_list.length();

	for (int i=0; i<nWidgets; i++) {
		QWidget *w = widget_list[i];
        QString wname = w->objectName();
		if (wname.startsWith("line_")) {
			connect(w, SIGNAL(textChanged(QString)), this, SLOT(changeParam()));
			connect(w, SIGNAL(textChanged(QString)), this, SLOT(redrawDistPlot()));
		}
		if (wname.startsWith("text_")) {
			connect(w, SIGNAL(textChanged(QString)), this, SLOT(changeParam()));
		}
		if (wname.startsWith("spin_")) {
			connect(w, SIGNAL(valueChanged(int)), this, SLOT(changeParam()));
		}
		if (wname.startsWith("comb_")) {
			connect(w, SIGNAL(activated(QString)), this, SLOT(changeParam()));
		}
		if (wname.startsWith("cbox_")) {
			connect(w, SIGNAL(toggled(bool)), this, SLOT(changeParam()));
		}
		if (wname.startsWith("rbut_")) {
			connect(w, SIGNAL(toggled(bool)), this, SLOT(changeParam()));
		}
	}

	QwtPlot *qp;

	qp = (QwtPlot *)qFindChild<QObject *>(this, "qwtPlot_TC_AVIDITY");
	distplot_list[0] = qp;
	qp = (QwtPlot *)qFindChild<QObject *>(this, "qwtPlot_DIVIDE1");
	distplot_list[1] = qp;
	qp = (QwtPlot *)qFindChild<QObject *>(this, "qwtPlot_DIVIDE2");
	distplot_list[2] = qp;
	qp = (QwtPlot *)qFindChild<QObject *>(this, "qwtPlot_DC_ANTIGEN");
	distplot_list[3] = qp;
	qp = (QwtPlot *)qFindChild<QObject *>(this, "qwtPlot_DC_LIFETIME");
	distplot_list[4] = qp;
    qp = (QwtPlot *)qFindChild<QObject *>(this, "qwtPlot_STIM_HILL");
    distplot_list[5] = qp;
    qp = (QwtPlot *)qFindChild<QObject *>(this, "qwtPlot_BINDTIME_HILL");
    distplot_list[6] = qp;
}



//-------------------------------------------------------------
// Loops through the workingParameterList and fills in the GUI.
//-------------------------------------------------------------
void MainWindow::loadParams()
{
	sprintf(msg,"nParams: %d nSliders: %d",nParams,nSliders);
    LOG_QMSG(msg);
	for (int k=0; k<nSliders; k++) {		// there must be a neater way to do this
		SliderPlus *splus = 0;
		sliderplus_list.append(splus);
		QWidget *w = 0;
		sliderParam.append(w);		
	}
	if (param_to_sliderIndex == NULL) {
		param_to_sliderIndex = new int[nParams];
		for (int i=0; i<nParams; i++)
			param_to_sliderIndex[i] = -1;
	}
	for (int i=0; i<nWidgets; i++) {
		QWidget *w = widget_list[i];							// w = widget_list[i] is the ith widget in the UI
        QString qsname = w->objectName();
		if (qsname.startsWith("line_") || qsname.startsWith("spin_")
			|| qsname.startsWith("comb_") || qsname.startsWith("cbox_")
			|| qsname.startsWith("rbut_") || qsname.startsWith("text_")) {
//            LOG_QMSG(qsname);
			QString wtag = qsname.mid(5);
			int rbutton_case = 0;
			if (qsname.startsWith("rbut_")) {
//				wtag = parse_rbutton(wtag,&rbutton_case);
                wtag = parse_rbutton(qsname,&rbutton_case);
            }
            // Find corresponding data in workingParameterList
            bool found = false;
			for (int k=0; k<nParams; k++) {
				PARAM_SET p = parm->get_param(k);
				QString ptag = p.tag;		// ptag is the tag of the kth parameter in the list
				if (wtag.compare(ptag) == 0) {
					double vmax = p.maxvalue;
					double vmin = p.minvalue;
                    // Update the widget (line_, spin_ or comb_) with data from the parameter list
                    // ---LineEdits
					if (qsname.startsWith("line_")) {
                        double val = p.value;
						QString val_str = QString::number(val);
						QLineEdit *w_l = (QLineEdit *)w;
                        w_l->setText(val_str);
						if (USE_RANGES) {
							// Set max and min values. If min=max=0, there're no restrictions.
							if (!(vmin == 0 && vmax == 0)) {
//	                            QValidator *aValidator = new QDoubleValidator(vmin, vmax, 10, w_l);
								QValidator *aValidator = new MyDoubleValidator(vmin, vmax, 8, w_l);
								w_l->setValidator(aValidator);
							}
						}						
					} else if (qsname.startsWith("spin_")) {
						double val = p.value;
						QSpinBox *w_s = (QSpinBox *)w;
                        w_s->setValue(val);
						if (!(vmin == 0 && vmax == 0)) {
                            w_s->setMinimum(vmin);
                            w_s->setMaximum(vmax);
						}
						if (qsname.contains("NCPU")) {
							ncpu = p.value;
						}
					} else if (qsname.startsWith("comb_")) {
                        int val = p.value - 1;	//0-based indexing
						QComboBox *w_c = (QComboBox *)w;
                        w_c->setCurrentIndex(val);
					} else if (qsname.startsWith("cbox_")) {
						QCheckBox *w_cb = (QCheckBox *)w;
						bool in_vitro = qsname.contains("IN_VITRO");
						if (p.value == 1) {
							w_cb->setChecked(true);
							if (in_vitro) 
								enableInVitro();
						} else {
							w_cb->setChecked(false);
							if (in_vitro) 
								disableInVitro();
						}
						bool dc_injection = qsname.contains("DC_INJECTION");
						if (p.value == 1) {
							w_cb->setChecked(true);
							if (dc_injection)
								enableDCInjection();
						} else {
							w_cb->setChecked(false);
							if (dc_injection)
								disableDCInjection();
						}
						bool use_traffic = qsname.contains("USE_TRAFFIC");
						if (p.value == 1) {
							w_cb->setChecked(true);
							if (use_traffic)
								enableUseTraffic();
						} else {
							w_cb->setChecked(false);
							if (use_traffic)
								disableUseTraffic();
						}
						bool use_exit_chemotaxis = qsname.contains("USE_EXIT_CHEMOTAXIS");
						if (p.value == 1) {
							w_cb->setChecked(true);
							if (use_exit_chemotaxis)
								enableUseExitChemotaxis();
						} else {
							w_cb->setChecked(false);
							if (use_exit_chemotaxis)
								disableUseExitChemotaxis();
						}
						bool use_DC_chemotaxis = qsname.contains("USE_DC_CHEMOTAXIS");
						if (p.value == 1) {
							w_cb->setChecked(true);
							if (use_DC_chemotaxis)
								enableUseDCChemotaxis();
						} else {
							w_cb->setChecked(false);
							if (use_DC_chemotaxis)
								disableUseDCChemotaxis();
						}
                        bool use_desens = qsname.contains("USE_DESENSITISATION");
                        if (p.value == 1) {
                            w_cb->setChecked(true);
                            if (use_desens)
                                enableUseDesensitisation();
                        } else {
                            w_cb->setChecked(false);
                            if (use_desens)
                                disableUseDesensitisation();
                        }
                        bool use_overstim = qsname.contains("USE_OVERSTIMULATION");
                        bool set_true = (p.value == 1);
                        w_cb->setChecked(set_true);
                        if (use_overstim) {
                            enableUseOverstimulation(set_true);
                        }
//                        if (p.value == 1) {
//                            w_cb->setChecked(true);
//                            if (use_overstim)
//                                enableUseOverstimulation();
//                        } else {
//                            w_cb->setChecked(false);
//                            if (use_overstim)
//                                disableUseOverstimulation();
//                        }

                        bool sim_periphery = qsname.contains("SIMULATE_PERIPHERY");
                        if (p.value == 1) {
                            w_cb->setChecked(true);
                            if (sim_periphery)
                                enableSimulatePeriphery();
                        } else {
                            w_cb->setChecked(false);
                            if (sim_periphery)
                                disableSimulatePeriphery();
                        }
                    } else if (qsname.startsWith("rbut_")) {
						QRadioButton *w_rb = (QRadioButton *)w;
						if (p.value == rbutton_case) {
							w_rb->setChecked(true);
						} else {
							w_rb->setChecked(false);
						}
					} else if (qsname.startsWith("text_")) {
						QLineEdit *w_l = (QLineEdit *)w;
						w_l->setText(p.label);
					}
					
					// Update Label text (except for "text_" variables)
                    // Get the corresponding label from the label list
                    QString labelString = "label_" + wtag;
					QLabel *label = NULL;
					bool foundLabel = false;
					for (int j=0; j<nLabels; j++) {
						label = label_list[j];
						if (!qsname.startsWith("text_") && labelString.compare(label->objectName()) == 0) {
							foundLabel = true;
							break;
						}
					}										// label is the pointer to the UI label for wtag and ptag
                    QString labelText = p.label;
                    
                    // Hardcode the distribution label names for now
                    if (wtag.compare("TC_AVIDITY_MEDIAN") == 0)
                        labelText = "Median";
					else if (wtag.compare("TC_AVIDITY_SHAPE") == 0)
                        labelText = "Shape";
                    else if (wtag.compare("DIVIDE1_MEDIAN") == 0)
                        labelText = "Median";
					else if (wtag.compare("DIVIDE1_SHAPE") == 0)
                        labelText = "Shape";
                    else if (wtag.compare("DIVIDE2_MEDIAN") == 0)
                        labelText = "Median";
					else if (wtag.compare("DIVIDE2_SHAPE") == 0)
                        labelText = "Shape";
                    else if (wtag.compare("DC_ANTIGEN_MEDIAN") == 0)
                        labelText = "Median";
					else if (wtag.compare("DC_ANTIGEN_SHAPE") == 0)
                        labelText = "Shape";
                    else if (wtag.compare("DC_LIFETIME_MEDIAN") == 0)
                        labelText = "Median";
					else if (wtag.compare("DC_LIFETIME_SHAPE") == 0)
                        labelText = "Shape";


					bool is_slider = false;
					int j;
					QSlider *s;
					QString sliderString;
					for (j=0; j<nSliders; j++) {
						sliderString = "slider_" + wtag;
						s = slider_list[j];
						if (sliderString.compare(s->objectName()) == 0) {
							is_slider = true;					// the jth slider in the list corresponds to wtag and ptag
							break;
						}
					}

					// Try this change to eliminate sliders except for distributions
					if (labelText.compare("Shape") != 0 && labelText.compare("Median") != 0) {
						is_slider = false;
					}

					if (is_slider) {
                        // If there is a slider corresponding to wtag, then just use the label.
                        if (foundLabel) {
	                        label->setText(labelText);
                        }
					} else {
						if (!(vmin == 0 && vmax == 0)) {
                            // If there is no slider, then add min and max values to the label text.
							QString min_str = QString::number(vmin);
							QString max_str = QString::number(vmax);
							if (foundLabel)
		                        label->setText(labelText + "  [ " + min_str + "-" + max_str + " ]");
						} else {
							if (foundLabel)
		                        label->setText(labelText);
						}
					}
						
                    // If there is a corresponding slider for this parameter, then apply settings.
					if (is_slider) {						
                        SliderPlus *splus = new SliderPlus(wtag,vmin,vmax,nTicks,k,i);
						sliderplus_list[j] = splus;
                        int ival = splus->val_to_int(p.value);
                        s->setMinimum(0);
                        s->setMaximum(splus->nTicks());
						s->setSliderPosition(ival);
						sliderParam[j] = w;
                        connect(s, SIGNAL(valueChanged(int)), this, SLOT(updateSliderBox())); //sliderReleased()   // valueChanged(int)
                        param_to_sliderIndex[k] = j;
					}                  
                    found = true;
                    break;

					if (!found) {
						sprintf(msg,"%s was not found in the parameter list",(wtag.toStdString()).data());
						LOG_MSG(msg);
					}
				}
			}
		}
	}
}

//--------------------------------------------------------------------------------------------------------
// This really should be changed.  Note that it requires that there is only one more "_" after "rbut_"
// We really want to split wtag into the parts before and after the last '_'
//--------------------------------------------------------------------------------------------------------
QString MainWindow::parse_rbutton(QString qsname, int *rbutton_case)
{
	// parse wtag into part before '_' and part after '_'
    QString wtag = qsname.mid(5);   // strips off "rbut_"
    int j = wtag.lastIndexOf('_');  // position of last '_'
	QString suffix = wtag.mid(j+1);
    // the prefix becomes wtag0, the suffix becomes rbutton_case, an integer 0,1,2,...
    QString wtag0 = wtag.mid(0,j);
	bool ok;
	*rbutton_case = suffix.toInt(&ok);
    return wtag0;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::reloadParams()
{
	for (int i=0; i<nWidgets; i++) {
		QWidget *w = widget_list[i];							// w = widget_list[i] is the ith widget in the UI
        QString qsname = w->objectName();
		if (qsname.startsWith("line_") || qsname.startsWith("spin_") 
			|| qsname.startsWith("comb_") || qsname.startsWith("cbox_")
            || qsname.startsWith("rbut_") || qsname.startsWith("text_")) {
			QString wtag = qsname.mid(5);
			int rbutton_case = 0;
			if (qsname.startsWith("rbut_")) {
//				wtag = parse_rbutton(wtag,&rbutton_case);
                wtag = parse_rbutton(qsname,&rbutton_case);
            }
            // Find corresponding data in workingParameterList
            bool found = false;
			for (int k=0; k<nParams; k++) {
				PARAM_SET p = parm->get_param(k);
				QString ptag = p.tag;		// ptag is the tag of the kth parameter in the list
				if (wtag.compare(ptag) == 0) {
					found = true;
                    // Update the widget (line_, spin_ or comb_) with data from the parameter list
					if (qsname.startsWith("line_")) {
                        double val = p.value;
						QString val_str = QString::number(val);
						QLineEdit *w_l = (QLineEdit *)w;
                        w_l->setText(val_str);
					} else if (qsname.startsWith("text_")) {
						QLineEdit *w_l = (QLineEdit *)w;
						w_l->setText(p.label);
					} else if (qsname.startsWith("spin_")) {
						double val = p.value;
						QSpinBox *w_s = (QSpinBox *)w;
                        w_s->setValue(val);
						if (qsname.contains("NCPU")) {
							ncpu = p.value;
						}
					} else if (qsname.startsWith("comb_")) {
                        int val = p.value - 1;	//0-based indexing
						QComboBox *w_c = (QComboBox *)w;
                        w_c->setCurrentIndex(val);
					} else if (qsname.startsWith("cbox_")) {
						QCheckBox *w_cb = (QCheckBox *)w;
						bool in_vitro = qsname.contains("IN_VITRO");
						if (p.value == 1) {
							w_cb->setChecked(true);
							if (in_vitro) 
								enableInVitro();
						} else {
							w_cb->setChecked(false);
							if (in_vitro) 
								disableInVitro();
						}
						bool dc_injection = qsname.contains("DC_INJECTION");
						if (p.value == 1) {
							w_cb->setChecked(true);
							if (dc_injection)
								enableDCInjection();
						} else {
							w_cb->setChecked(false);
							if (dc_injection)
								disableDCInjection();
						}
						bool use_traffic = qsname.contains("USE_TRAFFIC");
						if (p.value == 1) {
							w_cb->setChecked(true);
							if (use_traffic)
								enableUseTraffic();
						} else {
							w_cb->setChecked(false);
							if (use_traffic)
								disableUseTraffic();
						}
						bool use_exit_chemotaxis = qsname.contains("USE_EXIT_CHEMOTAXIS");
						if (p.value == 1) {
							w_cb->setChecked(true);
							if (use_exit_chemotaxis)
								enableUseExitChemotaxis();
						} else {
							w_cb->setChecked(false);
							if (use_exit_chemotaxis)
								disableUseExitChemotaxis();
						}
						bool use_DC_chemotaxis = qsname.contains("USE_DC_CHEMOTAXIS");
						if (p.value == 1) {
							w_cb->setChecked(true);
							if (use_DC_chemotaxis)
								enableUseDCChemotaxis();
						} else {
							w_cb->setChecked(false);
							if (use_DC_chemotaxis)
								disableUseDCChemotaxis();
						}
                        bool use_desens = qsname.contains("USE_DESENSITISATION");
                        if (p.value == 1) {
                            w_cb->setChecked(true);
                            if (use_desens)
                                enableUseDesensitisation();
                        } else {
                            w_cb->setChecked(false);
                            if (use_desens)
                                disableUseDesensitisation();
                        }
                        bool use_overstim = qsname.contains("USE_OVERSTIMULATION");
                        bool set_true = (p.value == 1);
                        w_cb->setChecked(set_true);
                        if (use_overstim) {
                            enableUseOverstimulation(set_true);
                        }
                        bool sim_periphery = qsname.contains("SIMULATE_PERIPHERY");
                        if (p.value == 1) {
                            w_cb->setChecked(true);
                            if (sim_periphery)
                                enableSimulatePeriphery();
                        } else {
                            w_cb->setChecked(false);
                            if (sim_periphery)
                                disableSimulatePeriphery();
                        }
                    } else if (qsname.startsWith("rbut_")) {
						QRadioButton *w_rb = (QRadioButton *)w;
						if (p.value == rbutton_case) {
							w_rb->setChecked(true);
						} else {
							w_rb->setChecked(false);
						}
					}
				}
			}
			if (!found) {
//				LOG_MSG("Widget tag not found:");
//				LOG_QMSG(qsname);
//				LOG_QMSG(wtag);
			}
		}
	}				
}

//--------------------------------------------------------------------------------------------------------
// I see no reason to hide the description when the same text is sent twice.
//--------------------------------------------------------------------------------------------------------
void MainWindow::showMore(QString moreText)
{
//	LOG_MSG("label clicked!");
//	LOG_QMSG(moreText);
	
//	if ((int)sender() != currentDescription) {
        text_more->setEnabled(true); // self.ui.text_description.setEnabled(1) #show()
        text_more->setText(moreText); // text_description
        currentDescription = (int)sender();
//    } else {
//        text_more->clear(); // text_description
//        text_more->setEnabled(false); // hide()#text_description
//        currentDescription = 0;
//	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::writeout()
{
	QString line;
    QFile file(inputFile);
    if (!file.open(QFile::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot write file %1:\n%2.")
                             .arg(inputFile)
                             .arg(file.errorString()));
		LOG_MSG("File open failed");
        return;
    }
    QTextStream out(&file);
	for (int k=0; k<parm->nParams; k++) {
		PARAM_SET p = parm->get_param(k);
		double val = p.value;
		if (p.tag.compare("INPUT_FILE") == 0)
			line = p.label;
		else if (p.tag.compare("SPECIAL_CASE_FILE") == 0)
			line = p.label;
		else if (p.tag.compare("DC_INJECTION_FILE") == 0)
			line = p.label;
		else if (val == int(val)) 	// whole number, write as integer
			line = QString::number(int(val));
		else
			line = QString::number(val);
		int nch = line.length();
		for (int i=0; i<max(12-nch,1); i++)
			line += " ";
		line += p.tag;
		line += "\n";
		out << line;
	}

    paramSaved = true;
	LOG_MSG("Input data saved");
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::readInputFile()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open ..."), ".", tr("Input Files (*.inp)"));
	if (fileName.compare("") == 0)
		return;
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot read file %1:\n%2.")
                             .arg(fileName)
                             .arg(file.errorString()));
        return;
    }

    QTextStream in(&file);
	QString line;
	for (int k=0; k<parm->nParams; k++) {
		line = in.readLine();
		QStringList data = line.split(" ",QString::SkipEmptyParts);
		PARAM_SET p = parm->get_param(k);
		QString ptag = p.tag;
		if (ptag.compare("INPUT_FILE") == 0) {
			parm->set_label(k,data[0]);
		} if (ptag.compare("SPECIAL_CASE_FILE") == 0) {
			parm->set_label(k,data[0]);
		} else if (ptag.compare("DC_INJECTION_FILE") == 0) {
				parm->set_label(k,data[0]);
		} else {
			parm->set_value(k,data[0].toDouble());
		}
	}
    reloadParams();
    paramSaved = true;
	inputFile = fileName;
	QLineEdit *ql = findChild<QLineEdit*>("line_SPECIAL_CASE");
	QLineEdit *qt = findChild<QLineEdit*>("text_SPECIAL_CASE_FILE");
	int ispecial_case = ql->text().toInt();
	sprintf(msg,"ispecial_case: %d",ispecial_case);
	LOG_MSG(msg);
	if (ispecial_case != 0)
		qt->setEnabled(true);
	else
		qt->setEnabled(false);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::loadResultFile()
{
	RESULT_SET *R;

	R = new RESULT_SET;
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open ..."), ".", tr("Result Files (*.res)"));
	if (fileName.compare("") == 0)
		return;
	R->casename = QFileInfo(fileName).baseName();
	for(int i=0; i<result_list.size(); i++) {
		if (R->casename.compare(result_list[i]->casename) == 0) {
			QMessageBox::warning(this, tr("Open results"),
                             tr("This result file is already loaded"));
			return;
		}
	}
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot read file %1:\n%2.")
                             .arg(fileName)
                             .arg(file.errorString()));
        return;
    }

    QTextStream in(&file);
	R->nsteps = 0;
	bool indata = false;
	QString line;
	do {
		line = in.readLine();
		if (line.length() > 0) {
			QStringList datalist = line.split(" ",QString::SkipEmptyParts);
			if (indata) {
				R->nsteps++;
			}
			if (datalist[0].contains("========")) {
				indata = true;
			}
		}
	} while (!line.isNull());
	in.seek(0);

	R->tnow = new double[R->nsteps];
    for (int i=0; i<nGraphs; i++) {
        if (!grph->isTimeseries(i)) continue;
        if (!grph->isActive(i)) continue;
        R->pData[i] = new double[R->nsteps];
    }

    step = -1;
	indata = false;
	do {
		line = in.readLine();
		if (line.length() > 0) {
			QStringList dataList = line.split(" ",QString::SkipEmptyParts);
			if (indata) {
				int ndata = dataList.length();
				double *data = new double[ndata];
				for (int k=0; k<ndata; k++)
					data[k] = dataList[k].toDouble();
				step++;
				if (step >= R->nsteps) {
					LOG_MSG("ERROR: loadResultFile: step >= nsteps_p");
					return;
				}
				R->tnow[step] = step;		//data[1];step

				for (int i=0; i<nGraphs; i++) {
                    if (!grph->isTimeseries(i)) continue;
                    if (!grph->isActive(i)) continue;
					int k = grph->get_dataIndex(i);
					R->pData[i][step] = data[k]*grph->get_scaling(i);
				}
			}
			if (dataList[0].contains("========")) {
				indata = true;
			}
		}
	} while (!line.isNull());

	// Compute the maxima
    for (int i=0; i<nGraphs; i++) {
        if (!grph->isActive(i)) continue;
        if (grph->isTimeseries(i)) {
            double maxval = getMaximum(R,R->pData[i]);
            grph->set_maxValue(i,maxval);
            R->maxValue[i] = maxval;
        }
    }
	// Now add the result set to the list
	result_list.append(R);
//	if (nGraphCases == 0) {
//		if (show_outputdata)
//			box_outputData = 0;
//		initializeGraphs(R);
//        drawGraphs();
//		goToOutputs();
//	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
bool MainWindow::save()
{
	writeout();
	return true;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
bool MainWindow::saveAs()
{
    // show the file dialog
	QString fileName = QFileDialog::getSaveFileName(this, tr("Select Input File"), ".", tr("Input Files (*.inp)"));    
	if (fileName.compare("") != 0) {
		LOG_MSG("Selected file:");
		LOG_QMSG(fileName);
		inputFile = fileName;
        writeout();
	}
    // Otherwise if user chooses cancel ...
	return true;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
double MainWindow::getMaximum(RESULT_SET *R, double *x)
{
	double maxx = 0;
	for (int i=0; i<R->nsteps; i++)
		maxx = max(maxx,x[i]);
	return maxx;
}

//-------------------------------------------------------------
// Switches to the input screen
// For when the outputs are being displayed
//-------------------------------------------------------------
void MainWindow::goToInputs()
{
    stackedWidget->setCurrentIndex(0);
    action_inputs->setEnabled(false);
    action_outputs->setEnabled(true);
    action_VTK->setEnabled(true);
    action_FACS->setEnabled(true);
    showingVTK = false;
    showingFACS = false;
}

//-------------------------------------------------------------
// Switches to the output screen
//-------------------------------------------------------------
void MainWindow::goToOutputs()
{
    stackedWidget->setCurrentIndex(1);    
    action_outputs->setEnabled(false);
    action_inputs->setEnabled(true);
    action_VTK->setEnabled(true);
    action_FACS->setEnabled(true);
    showingVTK = false;
    showingFACS = false;
}

//-------------------------------------------------------------
// Switches to the VTK screen
//-------------------------------------------------------------
void MainWindow::goToVTK()
{
//    if (started) {
		stackedWidget->setCurrentIndex(2);
		action_outputs->setEnabled(true);
		action_inputs->setEnabled(true);
		action_VTK->setEnabled(false);
        action_FACS->setEnabled(true);
        showingVTK = true;
        showingFACS = false;
//    }
}

//-------------------------------------------------------------
// Switches to the FACS screen
//-------------------------------------------------------------
void MainWindow::goToFACS()
{
    stackedWidget->setCurrentIndex(3);
    action_outputs->setEnabled(true);
    action_inputs->setEnabled(true);
    action_VTK->setEnabled(true);
    action_FACS->setEnabled(false);
    showingVTK = false;
    showingFACS = true;
}

//-------------------------------------------------------------
// Load and play stored cell position data
// Currently disabled
//-------------------------------------------------------------
void MainWindow::playVTK()
{
	LOG_MSG("playVTK");
	// Select a file to play
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open ..."), ".", tr("Cell Path Files (*.pos)"));
	if (fileName.compare("") == 0)
		return;
    QMessageBox::StandardButton reply;
    reply = QMessageBox::question(this, tr("Animation player"),
									tr("Save animation frames to image files (.jpg)?"),
                                    QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel);
	bool save_image;
	if (reply == QMessageBox::Yes)
		save_image = true;
	else
		save_image = false;
	started = true;
    goToVTK();
	if (!vtk->startPlayer(QFileInfo(fileName).absoluteFilePath(), timer, save_image)) {
		LOG_MSG("startPlayer failed");
		errorPopup("Open failure on this file");
		return;
//		exit(1);
	}
    connect(timer, SIGNAL(timeout()), this, SLOT(timer_update()));
    timer->start(tickVTK);
	action_stop->setEnabled(true);
    action_pause->setEnabled(true);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::timer_update()
{
	ntimes ++;
	if (!vtk->nextFrame()) {
		LOG_MSG("Player completed");
		action_stop->setEnabled(false);
		action_pause->setEnabled(false);
//		timer->stop();
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::setVTKSpeed()
{
    bool ok;
    int i = QInputDialog::getInt(this, tr("Set speed"),
		tr("Player timer tick (ms): "), tickVTK, 10, 10000, 1, &ok);
	if (ok) {
		tickVTK = i;
	}
	if (started) {
		timer->stop();
		timer->start(tickVTK);
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::setSavePosStart()
{
    bool ok;
    int i = QInputDialog::getInt(this, tr("Set savepos start"),
		tr("Start recording cell positions at (hours): "), savepos_start, 0, 1000, 1, &ok);
	if (ok) {
		savepos_start = i;
		sprintf(msg,"savepos_start: %d",savepos_start);
		LOG_MSG(msg);
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
bool MainWindow::getVideoFileInfo(int *nframes, QString *itemFormat, QString *itemCodec, QString *videoFileName)
{
    bool ok;

    int i = QInputDialog::getInteger(this, tr("Set nframes"),tr("Number of frames to capture: "), *nframes, 0, 10000, 1, &ok);
    if (ok) {
        *nframes = i;
    }
    if (!ok || nframes == 0) {
        return false;
    }

    QStringList formatItems;
    formatItems << tr("avi") << tr("mov") << tr("mpg");
    *itemFormat = QInputDialog::getItem(this, tr("QInputDialog::getItem()"),
                                         tr("Video file format:"), formatItems, 0, false, &ok);
    QStringList codecItems;
    codecItems << tr("h264") << tr("mpeg4") << tr("mpeg");
    *itemCodec = QInputDialog::getItem(this, tr("QInputDialog::getItem()"),
                                              tr("Codec:"), codecItems, 0, false, &ok);

    const char *prompt;
    if (itemFormat->contains("avi")) {
        prompt = "Videos (*.avi)";
    } else if (itemFormat->contains("mov")) {
        prompt = "Videos (*.mov)";
    } else if (itemFormat->contains("mpg")) {
        prompt = "Videos (*.mpg)";
    }
    *videoFileName = QFileDialog::getSaveFileName(this,
                                                    tr("Save File"),
                                                    QString(),
                                                    tr(prompt));
    return true;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow:: startRecorderVTK()
{
    bool ok;
    int nframes=0;
    QString itemFormat, itemCodec, videoFileName;

    ok = getVideoFileInfo(&nframes, &itemFormat, &itemCodec, &videoFileName);
    if (!ok) return;
    goToVTK();
    videoVTK->startRecorder(videoFileName,itemFormat,itemCodec,nframes);
    actionStart_recording_VTK->setEnabled(false);
    actionStop_recording_VTK->setEnabled(true);
    recordingVTK = true;
    started = true;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow:: stopRecorderVTK()
{
    videoVTK->stopRecorder();
    actionStart_recording_VTK->setEnabled(true);
    actionStop_recording_VTK->setEnabled(false);
    recordingVTK = false;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow:: startRecorderFACS()
{
    bool ok;
    int nframes=0;
    QString itemFormat, itemCodec, videoFileName;

    ok = getVideoFileInfo(&nframes, &itemFormat, &itemCodec, &videoFileName);
    if (!ok) return;
    goToFACS();
    videoFACS->startRecorder(videoFileName,itemFormat,itemCodec,nframes);
    actionStart_recording_FACS->setEnabled(false);
    actionStop_recording_FACS->setEnabled(true);
    recordingFACS = true;
    LOG_QMSG("startRecorderFACS");
    LOG_QMSG(videoFileName);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow:: stopRecorderFACS()
{
    videoFACS->stopRecorder();
    actionStart_recording_FACS->setEnabled(true);
    actionStop_recording_FACS->setEnabled(false);
    recordingFACS = false;
    LOG_QMSG("stopRecorderFACS");
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::saveSnapshot()
{
	LOG_MSG("saveSnapshot");
	QString fileName = QFileDialog::getSaveFileName(this, tr("Select image file"), ".", 
		tr("Image files (*.png *.jpg *.tif *.bmp)"));    
	if (fileName.compare("") == 0) {
		goToVTK();
		return;
	}
	QFileInfo fi(fileName);
	QString imgType = fi.suffix();
	goToVTK();
	vtk->saveSnapshot(fileName,imgType);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::runServer()
{
	if (paused) {
		if (vtk->playing) {
			vtk->playon();
		} else {
			exthread->unpause();
		}
        action_run->setEnabled(false);
        action_pause->setEnabled(true);
        action_stop->setEnabled(true);
		action_save_snapshot->setEnabled(false);
        paused = false;
        return;
	}
	
	if (!paramSaved) {
		int response = QMessageBox::critical(this, tr("ABM Model GUI"), \
					tr("The document has been modified.\nPlease save changes before continuing."), \
					QMessageBox::Save | QMessageBox::Cancel); // | Qt.QMessageBox.Discard
		if (response == QMessageBox::Save) {
            save();
		} else if (response == QMessageBox::Cancel) {
            return;
		}
	}
	
//    if (!first) {
//        int response = QMessageBox::question(this, tr("ABM Model GUI"),
//                        tr("Would you like to clear the graphs from the previous run?"),
//                        QMessageBox::Yes | QMessageBox::No);
//        if (response == QMessageBox::Yes)
//            mdiArea->closeAllSubWindows();
//        else if (response == QMessageBox::Cancel)
//            return;
//    }
    first = false;

    // Display the outputs screen
    if (showingVTK) {
        goToVTK();
    } else if(showingFACS) {
        goToFACS();
    } else {
        goToOutputs();
    }
    // Disable parts of the GUI
    action_run->setEnabled(false);
    action_pause->setEnabled(true);
    action_stop->setEnabled(true);
    action_inputs->setEnabled(true);
    action_VTK->setEnabled(true);
    action_FACS->setEnabled(true);
    action_save_snapshot->setEnabled(false);
    tab_T->setEnabled(false);
    tab_DC->setEnabled(false);
    tab_TCR->setEnabled(false);
    tab_run->setEnabled(false);

	if (show_outputdata)
	    box_outputData = new QTextBrowser();
	else
		box_outputData = 0;

//    if (sthread0) sthread0->deleteLater();
//    if (sthread1) sthread1->deleteLater();

	if (use_CPORT1) {

		// Port 5001
        LOG_MSG("sthread1");
        if (!sthread1)
            sthread1 = new SocketHandler(CPORT1);
		connect(sthread1, SIGNAL(sh_output(QString)), this, SLOT(outputData(QString)));
		sthread1->start();
	}

	// Port 5000
    sthread0 = new SocketHandler(CPORT0);
	connect(sthread0, SIGNAL(sh_output(QString)), box_outputLog, SLOT(append(QString))); //self.outputLog)
	connect(sthread0, SIGNAL(sh_connected()), this, SLOT(preConnection()));
	connect(sthread0, SIGNAL(sh_disconnected()), this, SLOT(postConnection()));
	sthread0->start();
	vtk->cleanup();
	Sleep(100);

	hours = 0;
	nt_vtk = 0;
	for (int k=0; k<parm->nParams; k++) {
		PARAM_SET p = parm->get_param(k);
		if (p.tag.compare("NDAYS") == 0) {
			hours = p.value*24;
		}
		if (p.tag.compare("NT_ANIMATION") == 0) {
			nt_vtk = p.value;
		}
	}
	started = true;
    if (cbox_record->isChecked()) {
        double h1 = lineEdit_record_hour1->text().toDouble();
        double h2 = lineEdit_record_hour2->text().toDouble();
        recordfrom = int(h1*60/DELTA_T);
        recordto = int(h2*60/DELTA_T);
        sprintf(msg,"Recording frames: from: %d to: %d",recordfrom,recordto);
        LOG_MSG(msg);
    } else {
        recordfrom = -1;
        recordto = -1;
    }
    LOG_MSG("exthread");
//    if (!exthread)
        exthread = new ExecThread(inputFile);
    connect(exthread, SIGNAL(display()), this, SLOT(displayScene()));
	connect(exthread, SIGNAL(summary()), this, SLOT(showSummary()));
    connect(exthread, SIGNAL(action_VTK()), this, SLOT(goToVTK()));
    connect(exthread, SIGNAL(facs_update()), this, SLOT(showFACS()));
    connect(this, SIGNAL(facs_update()), this, SLOT(showFACS()));
    connect(exthread, SIGNAL(redimension(int)), this, SLOT(redimensionCellArrays(int)));
    exthread->ncpu = ncpu;
	exthread->nsteps = int(hours*60/DELTA_T);
	exthread->paused = false;
	exthread->stopped = false;
    exthread->summary_interval = 4*lineEdit_summary_interval->text().toInt();
    exthread->start();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::preConnection()
{
	LOG_MSG("preConnection");

    double hours = 0;
	for (int k=0; k<parm->nParams; k++) {
		PARAM_SET p = parm->get_param(k);
		if (p.tag.compare("NDAYS") == 0) {
			hours = p.value*24;
			break;
		}
	}

    if (!display_plots) return;
	// We assume that the model output is at hourly intervals

//    if (newR) {
//        LOG_MSG("deleting newR");
//        for (int i=0; i<grph->nGraphs; i++) {
//            delete  newR->pData[i];
//            newR->pData[i] = NULL;
//            delete pGraph[i];
//            pGraph[i] = NULL;
//        }
//        delete newR->tnow;
//        delete newR;
//        newR = NULL;
//    }

    LOG_MSG("creating newR");
    newR = new RESULT_SET;
	QString casename = QFileInfo(inputFile).baseName();
	vtkfile = casename + ".pos";
	newR->casename = casename;
	int nsteps = int(hours+1.5);
	newR->nsteps = nsteps;
	newR->tnow = new double[nsteps];
    sprintf(msg,"allocate pData: nGraphs: %d nsteps: %d",grph->nGraphs,nsteps);
    LOG_MSG(msg);
    for (int i=0; i<grph->nGraphs; i++) {
        newR->pData[i] = new double[nsteps];
        newR->pData[i][0] = 0;
	}
    LOG_MSG("preconnection: Allocated result set arrays");

	newR->tnow[0] = 0;	// These are not the right initial values
	step = -1;

	// Initialize graphs
	initializeGraphs(newR);
    LOG_MSG("did initializeGraphs");
    posdata = false;
	if (cbox_savepos->isChecked()) {
//		savepos_start = 0;
//		setSavePosStart();
		if (QFile::exists(vtkfile))
			QFile::remove(vtkfile);
	}
//	action_stop->setChecked(false);
	LOG_MSG("preconnection: done");
}

void MainWindow::errorPopup(QString errmsg)
{
//	QMessageBox::warning(this, tr("DLL error"), tr((errmsg.toStdString()).data()));
//	QMessageBox::warning(this, tr("DLL error"), tr("Got an error message"));
	LOG_QMSG(errmsg);
	QMessageBox msgBox;
	msgBox.setText(errmsg);
	msgBox.exec();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::initializeGraphs(RESULT_SET *R)
{
    LOG_MSG("initializeGraphs");
    mdiArea->closeAllSubWindows();
    mdiArea->show();
    setGraphsActive();
    int non_ts = 0;
//    if (field->isConcPlot()) non_ts++;
//    if (field->isVolPlot()) non_ts++;
//    if (field->isOxyPlot()) non_ts++;
    grph->makeGraphList(non_ts);
    nGraphs = grph->nGraphs;
    if (nGraphCases > 0) {
        clearAllGraphs();
    }
    QString tag;
    QString title;
    QString yAxisTitle;
    for (int i=0; i<nGraphs; i++) {
        if (!grph->isTimeseries(i) && !grph->isProfile(i) && !grph->isDistribution(i)) continue;   // ???
        tag = grph->get_tag(i);
        title = grph->get_title(i);
        yAxisTitle = grph->get_yAxisTitle(i);
        if (pGraph[i] != NULL) {
            pGraph[i]->deleteLater();
            pGraph[i] = NULL;
        }
        if (pGraph[i] == NULL) {
            pGraph[i] = new Plot(tag,R->casename);
            pGraph[i]->setTitle(title);
            pGraph[i]->setAxisTitle(QwtPlot::yLeft, yAxisTitle);
//            LOG_QMSG(title);
//            LOG_QMSG(tag);
        }
    }

    nGraphCases = 1;
    graphResultSet[0] = R;

    for (int i=0; i<nGraphs; i++) {
        if (!grph->isTimeseries(i) && !grph->isProfile(i) && !grph->isDistribution(i)) continue;
        mdiArea->addSubWindow(pGraph[i]);
        pGraph[i]->show();
    }

    if (show_outputdata) {
        mdiArea->addSubWindow(box_outputData);	// Need another way of creating this window - should be floating
        box_outputData->show();
    }
/*
    if (field->isConcPlot())
        field->makeConcPlot(mdiArea);
    if (field->isVolPlot())
        field->makeVolPlot(mdiArea);
    if (field->isOxyPlot())
        field->makeOxyPlot(mdiArea);
*/
    mdiArea->tileSubWindows();

    for (int i=0; i<nGraphs; i++) {
        if (!grph->isTimeseries(i)) continue;
        pGraph[i]->setAxisScale(QwtPlot::xBottom, 0, R->nsteps, 0);
    }
//    Global::dist_nv = 20;
}

/*
//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::initializeGraphs(RESULT_SET *R)
{
	mdiArea->closeAllSubWindows();
	mdiArea->show();
    setGraphsActive();
    grph->makeGraphList();
    LOG_QMSG("did makeGraphList");
    nGraphs = grph->nGraphs;
    if (nGraphCases > 0) {
        clearAllGraphs();
        LOG_QMSG("did clearAllGraphs");
    }

    QString tag;
    QString title;
    QString yAxisTitle;
    int k = 0;
    for (int i=0; i<nGraphs; i++) {
        if (grph->isActive(i)) {
            tag = grph->get_tag(i);
            title = grph->get_title(i);
            yAxisTitle = grph->get_yAxisTitle(i);
            k++;
        } else {
            tag = "";
            title = "";
            yAxisTitle = "";
        }
        if (k > maxGraphs) break;
            if (pGraph[i] != NULL) {
                pGraph[i]->deleteLater();
                pGraph[i] = NULL;
            }
            pGraph[i] = new Plot(tag,R->casename);
            pGraph[i]->setTitle(title);
            pGraph[i]->setAxisTitle(QwtPlot::yLeft, yAxisTitle);
    }
    LOG_QMSG("did setTitles");

	for (int i=0; i<nGraphs; i++) {
//        if (grph->isTimeseries(i)) {
            LOG_QMSG("addSubWindow: " + grph->get_tag(i));
            mdiArea->addSubWindow(pGraph[i]);
            pGraph[i]->show();
//        }
	}

	if (show_outputdata) {
		mdiArea->addSubWindow(box_outputData);	// Need another way of creating this window - should be floating
		box_outputData->show();
	}

    mdiArea->tileSubWindows();

	for (int i=0; i<nGraphs; i++) {
        if (!grph->isActive(i)) continue;
//        sprintf(msg,"i: %d isActive",i);
//        LOG_MSG(msg);
        if (grph->isTimeseries(i)) {
            pGraph[i]->setAxisScale(QwtPlot::xBottom, 0, R->nsteps, 0);
            pGraph[i]->setAxisTitle(QwtPlot::xBottom, "Time (hours)");
        }
	}
//    showmdiAreaSize();
}
*/

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::showmdiAreaSize()
{
    QRect rect;
    rect = mdiArea->geometry();
    int h = rect.height();
    int w = rect.width();
    sprintf(msg,"mdiArea w,h: %d %d",w,h);
    LOG_MSG(msg);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::drawGraphs()
{
    RESULT_SET *R;
    for (int kres=0; kres<Plot::ncmax; kres++) {
        R = graphResultSet[kres];
        if (R != 0) {
            for (int i=0; i<nGraphs; i++) {
                if (!grph->isActive(i)) continue;
                if (grph->isTimeseries(i)) {
                    int k = grph->get_dataIndex(i);
                    QString tag = grph->get_tag(i);
                    double yscale = grph->get_yscale(i);
                    pGraph[i]->redraw(R->tnow, R->pData[i], R->nsteps, R->casename, tag, yscale, false);
                    if (k == 0) {
                        grph->set_maxValue(i,R->maxValue[i]);
                    } else {
                        double maxval = grph->get_maxValue(i);
                        double newmax = R->maxValue[i];
                        if (newmax > maxval) {
                            grph->set_maxValue(i,newmax);
                        }
                    }
                }
            }
        }
    }
    for (int i=0; i<nGraphs; i++) {
        if (!grph->isTimeseries(i)) continue;
        if (!grph->isActive(i)) continue;
        double maxval = grph->get_maxValue(i);
        pGraph[i]->setYScale(maxval);
        pGraph[i]->replot();
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::redimensionCellArrays(int nbond_size)
{
    LOG_MSG("redimensionCellArrays");
    nbond_list = nbond_size;
    if (bond_list) {
//        free(bond_list);
        LOG_MSG("freed bond_list");
    }
//    if (!bond_list) bond_list = (int *)malloc(10*(nbond_size+1)*sizeof(int));
    LOG_MSG("allocated bond_list");
    redimflag = true;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::displayScene()
{
//    LOG_MSG("displayScene");
	bool redo = false;	// need to understand this
	started = true;
//	mutex2.lock();
	bool fast = true;
	vtk->get_cell_positions(fast);
    vtk->renderCells(redo);
    if (videoVTK->record) {
        videoVTK->recorder();
    } else if (actionStop_recording_VTK->isEnabled()) {
        actionStart_recording_VTK->setEnabled(true);
        actionStop_recording_VTK->setEnabled(false);
    }
    // Old-style image saving
//    if (record) {
//        sprintf(msg,"record: framenum: %d",framenum);
//        LOG_MSG(msg);
//        QString basename = lineEdit_recordFileName->text();
//        vtk->record(basename, framenum);
//        framenum++;
//    }
    mutex2.unlock();
}

//--------------------------------------------------------------------------------------------------------
// Currently summaryData[] holds istep,ntot,nborn.  Hourly intervals, i.e. every 240 timesteps
//--------------------------------------------------------------------------------------------------------
void MainWindow::showSummary()
{
	char msg[128];
    double yscale;
//    LOG_MSG("showSummary");
    if (!display_plots) return;

	step++;
	if (step >= newR->nsteps) {
		LOG_MSG("ERROR: step >= nsteps");
		return;
	}

    mutex1.lock();

	hour = summaryData[1]*DELTA_T/60;
	progress = int(100.*hour/hours);
	progressBar->setValue(progress);
	QString hourstr = QString::number(int(hour));
	hour_display->setText(hourstr);
//    sprintf(msg,"showSummary: step: %d summaryData[7]: %d hour: %6.1f",step,summaryData[7],hour);
//    LOG_MSG(msg);

	QString casename = newR->casename;
	newR->tnow[step] = step;		//summaryData[0];

	for (int i=0; i<nGraphs; i++) {
        if (!grph->isActive(i)) continue;
        if (grph->isTimeseries(i)) {
            QString tag = grph->get_tag(i);
            int k = grph->get_dataIndex(i);
            double scale = grph->get_scaling(i);
            newR->pData[i][step] = summaryData[k]*scale;
//            LOG_QMSG(tag);
//            sprintf(msg,"scale: %f data: %d pt: %f",scale, summaryData[k], newR->pData[i][step]);
//            LOG_MSG(msg);
        }
	}
	for (int i=0; i<nGraphs; i++) {
        if (!grph->isActive(i)) continue;
        if (grph->isTimeseries(i)) {
            QString tag = grph->get_tag(i);
            yscale = grph->get_yscale(i);
            pGraph[i]->redraw(newR->tnow, newR->pData[i], step+1, casename, tag, yscale, false);
        }
	}
//    LOG_QMSG("did ts graphs");

    for (int i=0; i<nGraphs; i++) {
        if (!grph->isActive(i)) continue;
        if (grph->isProfile(i)) {
            double *x, *y;
            double xscale;
            int n;
            QString tag = grph->get_tag(i);
            int k = grph->get_dataIndex(i);
            x = profile_x[k];
            y = profile_y[k];
            n = profile_n[k];
            xscale = grph->get_xscale(x[n-1]);
            double maxval = 0;
            for (int j=0; j<n; j++) {
                if (y[j] > maxval) maxval = y[j];
            }
            yscale = pGraph[i]->calc_yscale_ts(maxval);
//            yscale = grph->get_yscale(i);
            pGraph[i]->setAxisScale(QwtPlot::xBottom, 0, xscale, 0);
            if (k == PROFILE_GENERATION_LN) {
                pGraph[i]->setAxisScale(QwtPlot::xBottom, 0, 20, 0);
            } else if (k == PROFILE_CFSE){
                pGraph[i]->setAxisScale(QwtPlot::xBottom, -20.0, 1.0, 0);
            }
            pGraph[i]->setAxisTitle(QwtPlot::xBottom, tag);
            pGraph[i]->setAxisTitle(QwtPlot::yLeft, grph->get_yAxisTitle(i));
            pGraph[i]->redraw(x, y, n, casename, tag, yscale, true);
        }
    }
//    LOG_QMSG("did profile graphs");

    mutex1.unlock();
}
//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::outputData(QString qdata)
{
//	if (qdata.compare("VTK") == 0) {
	if (qdata.startsWith("VTK")) {
		qdata.replace(0,3,"");
		bool savepos = cbox_savepos->isChecked();
		if (savepos) {
			if (step < savepos_start) {
				savepos = false;
			}
		}
		vtk->read_cell_positions(cellfile, vtkfile, savepos);
		started = true;
        if (showingVTK || firstVTK) {
			firstVTK = false;
			bool redo = false;
//			if (showingVTK == 1) {
//				redo = true;
//			}
            vtk->renderCells(redo);
		} 
	    posdata = true;
		if (qdata.length() == 0)
			return;
	}
//	if (qdata.contains("__EXIT__",Qt::CaseSensitive) || qdata.contains("Fortran") ) {
	if (quitMessage(qdata) || qdata.contains("Fortran") ) {
		return;
	}
	if (show_outputdata)
	    box_outputData->append(qdata);

    QStringList dataList = qdata.split(" ",QString::SkipEmptyParts);
	double data[11];
	for (int k=0; k<11; k++)
		data[k] = dataList[k].toDouble();
	step++;
	if (step >= newR->nsteps) {
		LOG_MSG("ERROR: step >= nsteps");
		return;
	}
	QString casename = newR->casename;
    newR->tnow[step] = step;		//data[1];
	for (int i=0; i<nGraphs; i++) {
        if (!grph->isTimeseries(i)) continue;
        if (!grph->isActive(i)) continue;
		int k = grph->get_dataIndex(i);
		newR->pData[i][step] = data[k]*grph->get_scaling(i);
	}

	for (int i=0; i<nGraphs; i++) {
        if (!grph->isActive(i)) continue;
        if (grph->isTimeseries(i)) {
            QString tag = grph->get_tag(i);
            double yscale = grph->get_yscale(i);
            pGraph[i]->redraw(newR->tnow, newR->pData[i], step+1, casename, tag, yscale, true);
        }
	}

}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::postConnection()
{
	LOG_MSG("postConnection");
	if (use_CPORT1) {
		sthread1->socket->close();
		sthread1->tcpServer->close();
		sthread1->quit();
		sthread1->wait(100);
		if (sthread1->isRunning()) {
			LOG_MSG("sthread1 did not terminate");
		}
	}

    if (actionStop_recording_VTK->isEnabled()) {
        stopRecorderVTK();
    }
    action_run->setEnabled(true);
    action_pause->setEnabled(false);
    action_stop->setEnabled(false);
	action_save_snapshot->setEnabled(true);
    tab_T->setEnabled(true);
    tab_DC->setEnabled(true);
    tab_TCR->setEnabled(true);
    tab_run->setEnabled(true);

    if (!display_plots) return;

	// Check if a result set of this name is already in the list, if so remove it
	for (int i=0; i<result_list.size(); i++) {
		if (newR->casename.compare(result_list[i]->casename) == 0) {
			result_list.removeAt(i);
		}
	}
	// Compute the maxima
	for (int i=0; i<nGraphs; i++) {
        if (!grph->isTimeseries(i)) continue;
        if (!grph->isActive(i)) continue;
		double maxval = getMaximum(newR,newR->pData[i]);
		newR->maxValue[i] = maxval;
	}

	// Add the new result set to the list
//	result_list.append(newR);
//	vtk->renderCells(true,true);		// for the case that the VTK page is viewed only after the execution is complete
    posdata = false;
	LOG_MSG("completed postConnection");
//    if (exthread->isFinished()) {
//        LOG_MSG("exthread isFinished: delete");
//        exthread->exit();
//        exthread->wait();
//        delete exthread;
//    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::close_sockets()
{
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::pauseServer()
{
	if (vtk->playing) {
		vtk->pause();
		LOG_MSG("Paused the player");
	} else {
		exthread->pause();
		LOG_MSG("Paused the ABM program.");
	}
	paused = true;
	action_run->setEnabled(true); 
	action_pause->setEnabled(false);
	action_stop->setEnabled(true);
	action_save_snapshot->setEnabled(true);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::stopServer()
{
    if (vtk->playing) {
		vtk->stop();
		LOG_MSG("Stopped the player");
	} else {
		LOG_MSG("stop ordered");
		if (paused) {
			LOG_MSG("was paused, runServer before stopping");
			runServer();
		}
        exthread->snapshot();
		exthread->stop();
        sleep(1);		// delay for Fortran to wrap up (does this help?)
		if (use_CPORT1) {
			sthread1->quit();
			sthread1->terminate();
		}
		sthread0->stop();
        if (display_plots)
            newR->nsteps = step+1;
	}
    action_run->setEnabled(true); 
    action_pause->setEnabled(false);
    action_stop->setEnabled(false);
	action_save_snapshot->setEnabled(true);
    if (actionStop_recording_VTK->isEnabled()) {
        stopRecorderVTK();
    }
    if (actionStop_recording_FACS->isEnabled()) {
        stopRecorderFACS();
    }
    delete exthread;
    exthread = NULL;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::clearAllGraphs()
{
	if (nGraphCases > 0) {
        if (nGraphs > 0) {
            for (int i=0; i<nGraphs; i++) {
                if (!grph->isTimeseries(i)) continue;
                if (!grph->isActive(i)) continue;
                LOG_QMSG(grph->get_tag(i));
//                pGraph[i]->removeAllCurves();
            }
        }
		nGraphCases = 0;
	}
	for (int i=0; i<Plot::ncmax; i++) {
		graphCaseName[i] = "";
		graphResultSet[i] = 0;
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
QString MainWindow::selectResultSet()
{
    QStringList items;
	for (int k=0; k<result_list.size(); k++) {
		RESULT_SET *R = result_list[k];
		if (R == 0) continue;
		bool inlist = false;
		for (int i=0; i<Plot::ncmax; i++) {
			if (graphResultSet[i] == 0) continue;
			if (R->casename.compare(graphResultSet[i]->casename) == 0) {
				inlist = true;
				break;
			}
		}
		if (!inlist)
			items << R->casename;
	}
	if (items.size() == 0) {
		QMessageBox::warning(this, tr("Select result case"),
			tr("No result sets available - use 'File > Load results'"));
		return QString("");
	}

    bool ok;
    QString item = QInputDialog::getItem(this, tr("Select result case"),
		tr("Case:"), items, 0, false, &ok);
    if (ok && !item.isEmpty())
		return item;
	else
		return QString("");
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::addGraph()
{
	// Need to select a result set from result_list then add the corresponding curves
        RESULT_SET *R = NULL;

	if (nGraphCases == Plot::ncmax) {
		QString mess = QString("The maximum number of cases is %1").arg(Plot::ncmax);
		QMessageBox::warning(this, tr("Add graph"),	tr((mess.toStdString()).data())); 
		return;
	}
	QString casename = selectResultSet();
	if (casename.compare("") == 0)
		return;

	for (int k=0; k<result_list.size(); k++) {
		if (casename.compare(result_list[k]->casename) == 0) {
			R = result_list[k];		// OK after doing a run or a load, followed by another load
			break;
		}
	}

	graphResultSet[nGraphCases] = R;
	nGraphCases++;
	// First add the curves
    for (int i=0; i<nGraphs; i++) {
        if (!grph->isActive(i)) continue;
        if (grph->isTimeseries(i)) {
            pGraph[i]->addCurve(R->casename);
            pGraph[i]->setAxisAutoScale(QwtPlot::xBottom);
        }
	}
	// Now redraw with the data
	drawGraphs();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
int MainWindow::selectGraphCase()
{
    QStringList items;
	for (int i=0; i<Plot::ncmax; i++) {
		if (graphResultSet[i] == 0) continue;
		items << graphResultSet[i]->casename;
	}
	if (items.size() == 0) {
		QMessageBox::warning(this, tr("Select graph case"),
			tr("No graph cases to remove"));
		return -1;
	}

    bool ok;
    QString item = QInputDialog::getItem(this, tr("Select graph case"),
		tr("Case:"), items, 0, false, &ok);
	if (ok && !item.isEmpty()) {
		for (int i=0; i<Plot::ncmax; i++) {
			if (graphResultSet[i] == 0) continue;
			if (item.compare(graphResultSet[i]->casename) == 0) {
				return i;
			}
		}
	}
	return -1;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::removeGraph()
{
	int i = selectGraphCase();
	if (i == -1) return;
	RESULT_SET *R = graphResultSet[i];
	// First remove the curves
	for (int i=0; i<nGraphs; i++) {
        if (!grph->isTimeseries(i)) continue;
        if (!grph->isActive(i)) continue;
		pGraph[i]->removeCurve(R->casename);
	}
	// Then remove the graph case
	graphResultSet[i] = 0;
	nGraphCases--;
	drawGraphs();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::removeAllGraphs()
{
	clearAllGraphs();
}

//---------------------------------------------------------------------
// Updates input parameter (QLineEdit) widgets according to the slider
//---------------------------------------------------------------------
void MainWindow::updateSliderBox()
{
    paramSaved = false;          // Keeps track of the fact that a param has been changed but not saved.
    // ---Get value from slider
    QString slider_str = sender()->objectName();
    QString stag = slider_str.mid(7);
    int ival = ((QSlider *)sender())->value();

    // --- Get param index from workingParameterList
	int k;
	for (k=0; k<nParams; k++) {
		PARAM_SET p = parm->get_param(k);
		if (stag.compare(p.tag) == 0)
			break;
	}
    int j = param_to_sliderIndex[k];
    SliderPlus *sp = sliderplus_list[j];
    double v = sp->int_to_val(ival);
    QString vstr = sp->val_to_str(v);
    if (((QSlider *)sender())->isSliderDown()) {    // This stops strange effects when changing the lineEdit field
        ((QLineEdit *)sliderParam[j])->setText(vstr);
    }
}

//------------------------------------------------------------------------------------------------------
// changeParam() is invoked in response to a signal sent when a value in a QLineEdit etc. widget
// is changed.  Note that when a QRadioButton widget is changed, signals are sent both the radiobuttons
// that change, but only one signal is used to change the parameter value.
//------------------------------------------------------------------------------------------------------
void MainWindow::changeParam()
{
	paramSaved = false;
    QObject *w = sender(); // Gets the pointer to the object that invoked the changeParam slot.
	if (w->isWidgetType()) {
		QString wname = w->objectName();
//		LOG_QMSG("changeParam:");
//		LOG_QMSG(wname);
		if (wname.contains("line_")) {
			QString wtag = wname.mid(5);
			QLineEdit *lineEdit = (QLineEdit *)w;
            QString text = lineEdit->displayText();
			// Determine if there is a slider associated with the sender widget
			for (int k=0; k<parm->nParams; k++) {
				PARAM_SET p = parm->get_param(k);
				if (wtag.compare(p.tag) == 0) {
					int j = param_to_sliderIndex[k];
					if (j >= 0) {
                        QSlider *slider = slider_list[j];
                        SliderPlus *sp = sliderplus_list[j];
                        double v = sp->str_to_val(text);
                        int ival = sp->val_to_int(v);
                        int ival_old = sp->val_to_int(p.value);
						if (ival != ival_old) {
                            slider->setSliderPosition(ival);
						}
					}
					parm->set_value(k,text.toDouble());
					break;
				}
			}
		} else if (wname.contains("text_")) {
			QString wtag = wname.mid(5);
			QLineEdit *lineEdit = (QLineEdit *)w;
			QString text = lineEdit->displayText();
			for (int k=0; k<parm->nParams; k++) {
				PARAM_SET p = parm->get_param(k);
				if (wtag.compare(p.tag) == 0) {
					parm->set_label(k,text);
					break;
				}
			}
		} else if (wname.contains("spin_")) {
			QSpinBox *spinBox = (QSpinBox *)w;
            int v = spinBox->value();
			QString wtag = wname.mid(5);
			for (int k=0; k<parm->nParams; k++) {
				PARAM_SET p = parm->get_param(k);
				if (wtag.compare(p.tag) == 0) {
					parm->set_value(k,v);
                    break;
				}
				if (wname.contains("NCPU")) {
					ncpu = v;
				}
			}
		} else if (wname.contains("cbox_")) {
			QCheckBox *checkBox = (QCheckBox *)w;
			bool in_vitro = wname.contains("IN_VITRO");
			int v;
			if (checkBox->isChecked()) {
				v = 1;
				if (in_vitro) 
					enableInVitro();
			} else {
				v = 0;
				if (in_vitro) 
					disableInVitro();
			}
			bool dc_injection = wname.contains("DC_INJECTION");
			if (checkBox->isChecked()) {
				v = 1;
				if (dc_injection)
					enableDCInjection();
			} else {
				v = 0;
				if (dc_injection)
					disableDCInjection();
			}
			bool use_traffic = wname.contains("USE_TRAFFIC");
			if (checkBox->isChecked()) {
				v = 1;
				if (use_traffic)
					enableUseTraffic();
			} else {
				v = 0;
				if (use_traffic)
					disableUseTraffic();
			}
			bool use_exit_chemotaxis = wname.contains("USE_EXIT_CHEMOTAXIS");
			if (checkBox->isChecked()) {
				v = 1;
				if (use_exit_chemotaxis)
                    enableUseExitChemotaxis();
			} else {
				v = 0;
				if (use_exit_chemotaxis)
					disableUseExitChemotaxis();
			}
			bool use_DC_chemotaxis = wname.contains("USE_DC_CHEMOTAXIS");
			if (checkBox->isChecked()) {
				v = 1;
				if (use_DC_chemotaxis)
					enableUseDCChemotaxis();
			} else {
				v = 0;
				if (use_DC_chemotaxis)
					disableUseDCChemotaxis();
			}
            bool use_desens = wname.contains("USE_DESENSITISATION");
            if (checkBox->isChecked()) {
                v = 1;
                if (use_desens)
                    enableUseDesensitisation();
            } else {
                v = 0;
                if (use_desens)
                    disableUseDesensitisation();
            }
            bool use_overstim = wname.contains("USE_OVERSTIMULATION");
            bool set_true = checkBox->isChecked();
            if (set_true)
                v = 1;
            else
                v = 0;
            if (use_overstim) {
                enableUseOverstimulation(set_true);
            }
            bool sim_periphery = wname.contains("SIMULATE_PERIPHERY");
            if (checkBox->isChecked()) {
                v = 1;
                if (sim_periphery)
                    enableSimulatePeriphery();
            } else {
                v = 0;
                if (sim_periphery)
                    disableSimulatePeriphery();
            }

			QString wtag = wname.mid(5);
			for (int k=0; k<parm->nParams; k++) {
				PARAM_SET p = parm->get_param(k);
				if (wtag.compare(p.tag) == 0) {
					parm->set_value(k,v);
                    break;
				}
			}
			if (wname.contains("savepos")) {
				if (checkBox->isChecked()) {
					setSavePosStart();
				}
			}
		} else if (wname.contains("comb_")) {
			QComboBox *comboBox = (QComboBox *)w;
            int v = comboBox->currentIndex();
			QString wtag = wname.mid(5);
//			sprintf(msg,"combo: %s  currentIndex: %d",wtag,v);
//			LOG_MSG(msg);
			for (int k=0; k<parm->nParams; k++) {
				PARAM_SET p = parm->get_param(k);
				if (wtag.compare(p.tag) == 0) {
					parm->set_value(k,v+1);
                    break;
				}
			}
		} else if (wname.contains("rbut_")) {
			QRadioButton *radioButton = (QRadioButton *)w;
            QString wtag = wname.mid(5);
            for (int k=0; k<parm->nParams; k++) {
                PARAM_SET p = parm->get_param(k);
                if (wtag.compare(p.tag) == 0) {
                    if (radioButton->isChecked()) {
//                        LOG_QMSG("is Checked");
                        parm->set_value(k,1);
//                        LOG_QMSG("set_value");
                    } else {
//                        LOG_QMSG("is not Checked");
                        parm->set_value(k,0);
//                        LOG_QMSG("set_value");
                    }
                    break;
                }
            }
        }
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::enableInVitro()
{
	for (int i=0; i<lineEdit_list.length(); i++) {
		QLineEdit *w = lineEdit_list[i];
		QString wname = w->objectName();
		if (wname.contains("line_IV_") || wname.contains("cbox_IV_SHOW_NONCOGNATE")) {
			w->setEnabled(true);
		}
	}
	for (int i=0; i<checkbox_list.length(); i++) {
		QCheckBox *w = checkbox_list[i];
		QString wname = w->objectName();
		if (wname.contains("cbox_IV_SHOW_NONCOGNATE")) {
			w->setEnabled(true);
		}
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::disableInVitro()
{
	for (int i=0; i<lineEdit_list.length(); i++) {
		QLineEdit *w = lineEdit_list[i];
		QString wname = w->objectName();
		if (wname.contains("line_IV_") || wname.contains("cbox_IV_SHOW_NONCOGNATE")) {
			w->setEnabled(false);
		}
	}
	for (int i=0; i<checkbox_list.length(); i++) {
		QCheckBox *w = checkbox_list[i];
		QString wname = w->objectName();
		if (wname.contains("cbox_IV_SHOW_NONCOGNATE")) {
			w->setEnabled(false);
		}
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::enableDCInjection()
{
	for (int i=0; i<lineEdit_list.length(); i++) {
		QLineEdit *w = lineEdit_list[i];
		QString wname = w->objectName();
		if (wname.contains("text_DC_INJECTION_FILE")) {
			w->setEnabled(true);
		}
		if (wname.contains("line_DCrate_100k")) {
			w->setEnabled(false);
		}
		if (wname.contains("line_T_DC1")) {
			w->setEnabled(false);
		}
		if (wname.contains("line_T_DC2")) {
			w->setEnabled(false);
		}
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::disableDCInjection()
{
	for (int i=0; i<lineEdit_list.length(); i++) {
		QLineEdit *w = lineEdit_list[i];
		QString wname = w->objectName();
		if (wname.contains("text_DC_INJECTION_FILE")) {
			w->setEnabled(false);
		}
		if (wname.contains("line_DCrate_100k")) {
			w->setEnabled(true);
		}
		if (wname.contains("line_T_DC1")) {
			w->setEnabled(true);
		}
		if (wname.contains("line_T_DC2")) {
			w->setEnabled(true);
		}	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::enableUseTraffic()
{
	for (int i=0; i<lineEdit_list.length(); i++) {
		QLineEdit *w = lineEdit_list[i];
		QString wname = w->objectName();
		if (wname.contains("line_RESIDENCE_TIME")) {	// || wname.contains("cbox_IV_SHOW_NONCOGNATE")) {
			w->setEnabled(true);
		}
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::disableUseTraffic()
{
	for (int i=0; i<lineEdit_list.length(); i++) {
		QLineEdit *w = lineEdit_list[i];
		QString wname = w->objectName();
		if (wname.contains("line_RESIDENCE_TIME")) {	//|| wname.contains("cbox_IV_SHOW_NONCOGNATE")) {
			w->setEnabled(false);
		}
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::enableUseExitChemotaxis()
{
	for (int i=0; i<lineEdit_list.length(); i++) {
		QLineEdit *w = lineEdit_list[i];
		QString wname = w->objectName();
		if (wname.contains("line_CHEMO_K_EXIT")) {	
			w->setEnabled(true);
		}
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::disableUseExitChemotaxis()
{
	for (int i=0; i<lineEdit_list.length(); i++) {
		QLineEdit *w = lineEdit_list[i];
		QString wname = w->objectName();
		if (wname.contains("line_CHEMO_K_EXIT")) {	
			w->setEnabled(false);
		}
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::enableUseDCChemotaxis()
{
	for (int i=0; i<lineEdit_list.length(); i++) {
		QLineEdit *w = lineEdit_list[i];
		QString wname = w->objectName();
		if (wname.contains("line_CHEMO_K_DC")) {	
			w->setEnabled(true);
		}
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::disableUseDCChemotaxis()
{
	for (int i=0; i<lineEdit_list.length(); i++) {
		QLineEdit *w = lineEdit_list[i];
		QString wname = w->objectName();
		if (wname.contains("line_CHEMO_K_DC")) {	
			w->setEnabled(false);
		}
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::enableUseDesensitisation()
{
    for (int i=0; i<lineEdit_list.length(); i++) {
        QLineEdit *w = lineEdit_list[i];
        QString wname = w->objectName();
        if (wname.contains("line_DESENS")) {
            w->setEnabled(true);
        }
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::disableUseDesensitisation()
{
    for (int i=0; i<lineEdit_list.length(); i++) {
        QLineEdit *w = lineEdit_list[i];
        QString wname = w->objectName();
        if (wname.contains("line_DESENS")) {
            w->setEnabled(false);
        }
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::enableUseOverstimulation(bool set)
{
    groupBox_overstimulation->setEnabled(set);
}

//--------------------------------------------------------------------------------------------------------
// No need for this at the moment
//--------------------------------------------------------------------------------------------------------
void MainWindow::enableSimulatePeriphery()
{
}

//--------------------------------------------------------------------------------------------------------
// No need for this at the moment
//--------------------------------------------------------------------------------------------------------
void MainWindow::disableSimulatePeriphery()
{
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_line_SPECIAL_CASE_textEdited(QString str)
{
	QLineEdit *ql = findChild<QLineEdit*>("text_SPECIAL_CASE_FILE");
	int num = str.toInt();
	if (num == 0) 
		ql->setEnabled(false);
	else
		ql->setEnabled(true);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::setupGraphSelector()
{
    QGridLayout *grid = new QGridLayout;
    int row[3];
    int col;
    row[0] = row[1] = row[2] = -1;

    cbox_ts = new QMyCheckBox*[grph->n_tsGraphs];
    for (int i=0; i<grph->n_tsGraphs; i++) {
        int itype = grph->tsGraphs[i].type;
        if (itype == 0) {
            if (i < 16)
                col = 0;
            else
                col = 1;
        } else {
            col = 2;
        }
        row[col]++;
        QString text = grph->tsGraphs[i].title;
        cbox_ts[i] = new QMyCheckBox;
        cbox_ts[i]->setText(text);
        cbox_ts[i]->setObjectName("cbox_"+grph->tsGraphs[i].tag);
        grid->addWidget(cbox_ts[i],row[col],col);
        connect((QObject *)cbox_ts[i], SIGNAL(checkBoxClicked(QString)), this, SLOT(showMore(QString)));
    }
    groupBox_graphselect->setLayout(grid);

    QRect rect = groupBox_graphselect->geometry();
#ifdef __DISPLAY768
    rect.setHeight(460);
#else
    rect.setHeight(500);
#endif
    groupBox_graphselect->setGeometry(rect);
}

/*
//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::setupGraphSelector-old()
{
//    QVBoxLayout *vbox = new QVBoxLayout;
    QGridLayout *grid = new QGridLayout;
    int row, col;
    int row0=-1, row1=-1;
    cbox_ts = new QCheckBox*[grph->n_tsGraphs];
    col = 0;
    for (int i=0; i<grph->n_tsGraphs; i++) {
//        if (grph->isTimeseries(i)) {
        if (grph->tsGraphs[i].ts) {
            row0++;

            row = row0;
        } else {
            if (row1 == -1) {
                col++;
            }
            row1++;
            row = row1;
        }
        QString text = grph->tsGraphs[i].title;
        cbox_ts[i] = new QCheckBox;
        cbox_ts[i]->setText(text);
        cbox_ts[i]->setObjectName("cbox_"+grph->tsGraphs[i].tag);
        cbox_ts[i]->setChecked(grph->tsGraphs[i].active);
        grid->addWidget(cbox_ts[i],row,col);
    }
//    groupBox_graphselect->setLayout(vbox);
    groupBox_graphselect->setLayout(grid);
    QRect rect;
    rect = groupBox_graphselect->geometry();
#ifdef __DISPLAY768
    rect.setHeight(480);
#else
    rect.setHeight(600);
#endif
    groupBox_graphselect->setGeometry(rect);
}
*/

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::setGraphsActive()
{
    for (int i=0; i<grph->n_tsGraphs; i++) {
        grph->tsGraphs[i].active = cbox_ts[i]->isChecked();
    }
}

//--------------------------------------------------------------------------------------------------------
void MainWindow:: initFACSPlot()
{
    qpFACS = (QwtPlot *)qFindChild<QObject *>(this, "qwtPlot_FACS");
    qpFACS->setTitle("FACS");
    QwtSymbol symbol = QwtSymbol( QwtSymbol::Rect, Qt::blue, Qt::NoPen, QSize( 2,2 ) );
    qpFACS->replot();
    connect((QObject *)groupBox_FACS,SIGNAL(groupBoxClicked(QString)),this,SLOT(processGroupBoxClick(QString)));
}

//--------------------------------------------------------------------------------------------------------
// Possible variables to plot against CFSE:
//   CD69, S1PR1, CD8, avidity, stimulation
// Note: nFACS_vars must be set to 1 + the number of variables (CFSE is 0)
//--------------------------------------------------------------------------------------------------------
void MainWindow::showFACS()
{
    double xmin, xmax, ymin, ymax, cfse, cvar, x, y;
    int i, k;
    int kvar;
    double scaling;
    QString ylabel;

    qpFACS = (QwtPlot *)qFindChild<QObject *>(this, "qwtPlot_FACS");
    qpFACS->clear();
    qpFACS->setTitle("FACS");
//    QwtSymbol symbol = QwtSymbol( QwtSymbol::Diamond, Qt::blue, Qt::NoPen, QSize( 2,4 ) );
    QwtSymbol symbol = QwtSymbol( QwtSymbol::Rect, Qt::blue, Qt::NoPen, QSize( 2,2 ) );
    xmin = 1.0e10;
    xmax = -1.0e10;
    ymin = 1.0e10;
    ymax = -1.0e10;
    if (radioButton_CD69->isChecked()) {
        kvar = 1;
        ylabel = "CD69";
        scaling = 1.0e4;
        ymin = 1;
        ymax = 1.0e4;
    } else if (radioButton_S1PR1->isChecked()) {
            kvar = 2;
            ylabel = "S1PR1";
            scaling = 1.0e4;
            ymin = 1;
            ymax = 1.0e4;
    } else if (radioButton_CD8->isChecked()) {
            kvar = 3;
            ylabel = "CD8";
            scaling = 0.1e4;
            ymin = 1;
            ymax = 1.0e4;
    } else if (radioButton_avidity->isChecked()) {
            kvar = 4;
            ylabel = "avidity";
            scaling = 1;
            ymin = 0.001;
            ymax = 1;
    } else if (radioButton_stimulation->isChecked()) {
            kvar = 5;
            ylabel = "stimulation";
            scaling = 1;
            ymin = 0.001;
            ymax = 1;
    }
    xmin = 1;
    xmax = 15000;
    for (i=0; i<nFACS_cells; i++) {
        cfse = FACS_data[nFACS_vars*i];
        cvar = scaling*FACS_data[nFACS_vars*i+kvar];
//        x = log(cfse)/log(2.);
//        y = cd69;
        x = 10000*cfse;
        y = max(cvar,1.01*ymin);
//        ymax = max(y,ymax);
        if (x >= xmin) {
            QwtPlotMarker* m = new QwtPlotMarker();
            m->setSymbol( symbol );
            m->setValue( QPointF( x,y ) );
            m->attach( qpFACS );
        }
    }
    qpFACS->setAxisScale(QwtPlot::yLeft, ymin, ymax, 0);
    qpFACS->setAxisTitle(QwtPlot::yLeft, ylabel);
    qpFACS->setAxisScaleEngine(QwtPlot::yLeft, new QwtLog10ScaleEngine);
    qpFACS->setAxisMaxMinor(QwtPlot::yLeft, 10);
    qpFACS->setAxisMaxMajor(QwtPlot::yLeft, 5);
    qpFACS->setAxisScale(QwtPlot::xBottom, xmin, xmax, 0);
    qpFACS->setAxisTitle(QwtPlot::xBottom, "CFSE");
    qpFACS->setAxisScaleEngine(QwtPlot::xBottom, new QwtLog10ScaleEngine);
    qpFACS->setAxisMaxMinor(QwtPlot::xBottom, 10);
    qpFACS->setAxisMaxMajor(QwtPlot::xBottom, 5);
    qpFACS->replot();
//    sprintf(msg,"x range: %f %f  y range: %f %f",xmin,xmax,ymin,ymax);
//    LOG_MSG(msg);

    if (videoFACS->record) {
        videoFACS->recorder();
    } else if (actionStop_recording_FACS->isEnabled()) {
        actionStart_recording_FACS->setEnabled(true);
        actionStop_recording_FACS->setEnabled(false);
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::processGroupBoxClick(QString text)
{
    LOG_QMSG("processGroupBoxClick: " + text);
    QwtPlot *plot;

//    if (text.compare("Histo") == 0) {
//        LOG_MSG("save Histo plot");
//        bool use_HistoBar = radioButton_histotype_1->isChecked();
//        if (use_HistoBar) {
//            plot = qpHistoBar;
//            qpHistoLine->hide();
//        } else {
//            plot = qpHistoLine;
//            qpHistoBar->hide();
//        }
    if (text.compare("FACS") == 0) {
        LOG_MSG("save FACS plot");
        plot = qpFACS;
    } else {
        return;
    }

    int w = plot->width();
    int h = plot->height();
    QPixmap pixmap(w, h);
    pixmap.fill(Qt::white); // Qt::transparent ?

    QwtPlotPrintFilter filter;
    int options = QwtPlotPrintFilter::PrintAll;
    options &= ~QwtPlotPrintFilter::PrintBackground;
    options |= QwtPlotPrintFilter::PrintFrameWithScales;
    filter.setOptions(options);

    plot->print(pixmap, filter);

//		QString fileName = getImageFile();
    QString fileName = QFileDialog::getSaveFileName(0,"Select image file", ".",
        "Image files (*.png *.jpg *.tif *.bmp)");
    if (fileName.isEmpty()) {
        return;
    }
    pixmap.save(fileName,0,-1);
}
//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow:: initDistPlots()
{
    QwtPlot *qp;
//    string name_str;

    for (int j=0; j<7; j++) {
        qp = distplot_list[j];
//        QString name = qp->objectName();
//        if (name.contains("TC_AVIDITY")) {
//            name_str = name.toStdString();
//			LOG_QMSG(name);
//        }
        if (j == 0) {
            qp->setTitle("TCR Avidity");
        } else if (j == 1) {
            qp->setTitle("First division (hrs)");
        } else if (j == 2) {
            qp->setTitle("Later division (hrs)");
        } else if (j == 3) {
            qp->setTitle("DC antigen density");
        } else if (j == 4) {
            qp->setTitle("DC lifetime (days)");
        } else if (j == 5) {
            qp->setTitle("Stimulation");
        } else if (j == 6) {
            qp->setTitle("Binding time");
        }

        QwtPlotCurve *curve = new QwtPlotCurve("title");
        curve->attach(qp);
        curve_list[j] = curve;
        qp->replot();
    }
}

//--------------------------------------------------------------------------------------------------------
// Redraws a probability distribution or Hill function plot when one of the parameters is changed.
// In the case of the distribution the change can be provided by a slider.
// If the name of the qwpplot object is qwtPlot_XXXX:
//   names of lineEdit objects with the median and shape must end with XXXX_median and XXXX_shape
//   names of lineEdit objects with the N and C must end with XXXX_N and XXXX_C
// The inital plotting is done when the lineEdit fields are initialized by loadParams()
//--------------------------------------------------------------------------------------------------------
void MainWindow::redrawDistPlot()
{
    QString sname = sender()->objectName();
//    LOG_QMSG(sname);
    for (int k=0; k<7; k++) {
		QwtPlot *qp = distplot_list[k];
        QString tag = qp->objectName().mid(8);
        QString tag_m = tag + "_MEDIAN";
        QString tag_s = tag + "_SHAPE";
        QString tag_thresh = tag + "_THRESHOLD";
        QString tag_N = tag + "_N";
        QString tag_C = tag + "_C";
        QString tag_min = tag + "_MIN";
        QString tag_max = tag + "_MAX";
        if (sname.endsWith(tag_m) || sname.endsWith(tag_s)) {
            int i_m = 0, i_s = 0;
            for (int i=0; i<nWidgets; i++) {
				QString wname = widget_list[i]->objectName();
                if (wname.endsWith(tag_m))
                    i_m = i;
                else if (wname.endsWith(tag_s))
                    i_s = i;
			}

            QString median_str = ((QLineEdit *)widget_list[i_m])->text();
            QString shape_str = ((QLineEdit *)widget_list[i_s])->text() ;
			if (median_str.compare("") == 0) return;
			if (shape_str.compare("") == 0) return;
            double median = median_str.toDouble();
            median = max(0.001, median);
            double shape = shape_str.toDouble();
            shape = max(1.0001, shape);
			
			double *x = new double[nDistPts];
			double *prob = new double[nDistPts];
            create_lognorm_dist(median,shape,nDistPts,x,prob);
            int n = dist_limit(prob,nDistPts);
            double xmax = x[n];
            curve_list[k]->setData(x, prob, n);
            qp->setAxisScale(QwtPlot::xBottom, 0.0, xmax, 0.0);
            qp->replot();
			delete [] x;
			delete [] prob;
		}
        if (sname.endsWith(tag_N) || sname.endsWith(tag_C) ||
                sname.endsWith(tag_thresh) || sname.endsWith(tag_min) || sname.endsWith(tag_max)) {
            int i_N = -1, i_C = -1, i_thresh = -1, i_min = -1, i_max = -1;
            for (int i=0; i<nWidgets; i++) {
                QString wname = widget_list[i]->objectName();
                if (wname.endsWith(tag_N))
                    i_N = i;
                else if (wname.endsWith(tag_C))
                    i_C = i;
                else if (wname.endsWith(tag_thresh))
                    i_thresh = i;
                else if (wname.endsWith(tag_min))
                    i_min = i;
                else if (wname.endsWith(tag_max))
                    i_max = i;
            }
            int hill_N;
            double hill_C, hill_thresh, hill_min, hill_max;
            if (i_N >= 0) {
                QString hill_N_str = ((QLineEdit *)widget_list[i_N])->text();
                hill_N = hill_N_str.toInt();
                hill_N = max(1, hill_N);
            } else {
                return;
            }
            if (i_C >= 0) {
                QString hill_C_str = ((QLineEdit *)widget_list[i_C])->text() ;
                hill_C = hill_C_str.toDouble();
                hill_C = max(0.01, hill_C);
                hill_C = min(1.0, hill_C);
            } else {
                return;
            }
            if (i_thresh >= 0) {
                QString hill_thresh_str = ((QLineEdit *)widget_list[i_thresh])->text() ;
                hill_thresh = hill_thresh_str.toDouble();
                hill_thresh = max(0.0, hill_thresh);
                hill_thresh = min(1.0, hill_thresh);
            } else {
                hill_thresh = 0;
            }
            if (i_min >= 0) {
                QString hill_min_str = ((QLineEdit *)widget_list[i_min])->text() ;
                hill_min = hill_min_str.toDouble();
                hill_min = max(0.0, hill_min);
                hill_min /= 60;     // convert mins to hrs
            } else {
                hill_min = 0.0;
            }
            if (i_max >= 0) {
                QString hill_max_str = ((QLineEdit *)widget_list[i_max])->text() ;
                hill_max = hill_max_str.toDouble();
                hill_max = max(1.0, hill_max);
            } else {
                hill_max = 1.0;
            }

            double *x = new double[nDistPts];
            double *hill = new double[nDistPts];
            create_hill_function(hill_N,hill_C,nDistPts,x,hill);
            if (i_thresh >= 0 || i_min >= 0 || i_max >= 0) {    // scale
                for (int i=0; i<nDistPts; i++) {
                    if (x[i] < hill_thresh) {
                        hill[i] = 0;
                    } else {
                        hill[i] = hill_min + (hill_max - hill_min)*hill[i];
                    }
                }
            }
            curve_list[k]->setData(x, hill, nDistPts);
            qp->setAxisScale(QwtPlot::xBottom, 0.0, 1.0, 0.0);
            qp->replot();
            delete [] x;
            delete [] hill;
        }
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
int MainWindow::dist_limit(double *p, int n)
{
	int i, imax;
    double pmax = 0;

	imax = n-1;
	for (i=0; i<n; i++) {
		if (p[i] > pmax) {
			pmax = p[i];
			imax = i;
		}
	}
    double plim = 0.01*pmax;
	for (i=n-1; i>0; i--) {
		if (p[i] > plim) {
            return min(n-1,max(i,2*imax));
		}
	}
	return 1;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
double MainWindow::erf(double z)
{
   double t = 1.0 / (1.0 + 0.5 * fabs(z));
   // use Horner's method
   double ans = 1 - t * exp( -z*z -  1.26551223 +
         t * ( 1.00002368 +
         t * ( 0.37409196 +
         t * ( 0.09678418 +
         t * (-0.18628806 +
         t * ( 0.27886807 +
         t * (-1.13520398 +
         t * ( 1.48851587 +
         t * (-0.82215223 +
         t * ( 0.17087277))))))))));
   if (z >= 0.0)
       return ans;
   else
       return -ans;
}

//-----------------------------------------------------------------------------------------
// When X is distributed N(mu,sig), this gives Prob{x1 < X <= x2}
//-----------------------------------------------------------------------------------------
double MainWindow::pnorm(double x1, double x2, double mu, double sig)
{
    double z1, z2, e1, e2;
		
	z1 = (x1-mu)/sig;
    e1 = erf(z1/sqrt(2.0))/2;
    z2 = (x2-mu)/sig;
    e2 = erf(z2/sqrt(2.0))/2;
    return e2 - e1;
}

//-----------------------------------------------------------------------------------------
// When log(X) is distributed N(mu,sig), this gives Prob{x1 < X <= x2}
//-----------------------------------------------------------------------------------------
double MainWindow::plognorm(double x1, double x2, double mu, double sig)
{
    double z1, z2, e1, e2;

    z1 = 0;
    z2 = 0;
    if (x1 == 0)
        e1 = -0.5;
	else {
        z1 = (log(x1)-mu)/sig;
        e1 = erf(z1/sqrt(2.0))/2;
	}
    if (x2 == 0)
        e2 = -0.5;
	else {
        z2 = (log(x2)-mu)/sig;
        e2 = erf(z2/sqrt(2.0))/2;
	}
    return e2 - e1;
}

//-----------------------------------------------------------------------------------------
// Create the lognormal distribution with median = p1, shape = p2
// at n points stored in x[], probability values stored in prob[].
// Note that x[0] = 0.
// The range of x is currently just less than 4*median.  This should be
// OK for values of shape < 2.
// Convert probability into probability density
//-----------------------------------------------------------------------------------------
void MainWindow::create_lognorm_dist(double p1, double p2,int n, double *x, double *prob)
{
	double xmax, dx, mu_l, sig_l, x1, x2;

    if (p1 >= 0.5)
        xmax = p1*4;
    else
        xmax = p1*8;
        
    dx = xmax/n;
    mu_l = log(p1);
    sig_l = log(p2);
	for (int ix=0; ix<n; ix++) {
        x1 = (ix - 0.5)*dx;
        x2 = x1 + dx;
		x1 = max(x1,0.0);
        x[ix] = (x1+x2)/2;
        prob[ix] = plognorm(x1,x2,mu_l,sig_l)/(x2-x1);
	}
}

void MainWindow::create_hill_function(int N, double C, int n, double *x, double *hill)
{
    double dx, cn;
    int i;

    dx = 1.0/(n-1);
    cn = pow(C,N);
    for (i=0; i<n; i++) {
        x[i] = i*dx;
        hill[i] = (1 + cn)*pow(x[i],N)/(pow(x[i],N) + cn);
    }
}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void MainWindow::on_rbut_ACTIVATION_MODE_0_toggled(bool checked)
{
    groupBox_STAGED->setEnabled(checked);
    groupBox_UNSTAGED->setEnabled(!checked);
}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void MainWindow::on_cbox_record_toggled(bool checked)
{
    lineEdit_recordFileName->setEnabled(checked);
    lineEdit_record_hour1->setEnabled(checked);
    lineEdit_record_hour2->setEnabled(checked);
    framenum = 0;
}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void MainWindow::on_checkBox_FACS_PLOT_toggled(bool checked)
{
    line_FACS_INTERVAL->setEnabled(checked);
    if (!checked) {
        line_FACS_INTERVAL->setText("0");
    }
}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void MainWindow::on_checkBox_EFFECTOR_FUNCTION_toggled(bool checked)
{
    line_CD8_EFFECTOR_PROB->setEnabled(checked);
    if (!checked) {
        line_CD8_EFFECTOR_PROB->setText("0");
    }
}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void MainWindow::exitRuleChanged()
{

}

//-------------------------------------------------------------
// Switches to the FACS screen
//-------------------------------------------------------------
//void MainWindow::on_action_FACS_triggered()
//{
//    stackedWidget->setCurrentIndex(3);
//    action_outputs->setEnabled(true);
//    action_inputs->setEnabled(true);
//    action_VTK->setEnabled(true);
//    action_FACS->setEnabled(false);
////    showingVTK = 1;
////    showingVTK += recording;
//}

//------------------------------------------------------------------------------------------------------
// This should be used for any radioButtonGroups
//------------------------------------------------------------------------------------------------------
void MainWindow::radioButtonChanged(QAbstractButton *b)
{
    QString wtag = b->objectName();
    int rbutton_case;
    if (b->isChecked()) {
        QString ptag = parse_rbutton(wtag,&rbutton_case);
        // Now need to reflect the change in the workingParameterList
        // Need to locate ptag
        for (int k=0; k<nParams; k++) {
            PARAM_SET p = parm->get_param(k);
            if (ptag.compare(p.tag) == 0) {
                parm->set_value(k,double(rbutton_case));
                break;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void MainWindow::on_action_show_gradient2D_triggered()
{
	LOG_QMSG("on_action_show_gradient2D_triggered");
	showGradient2D();
}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void MainWindow::on_action_show_gradient3D_triggered()
{
	LOG_QMSG("on_action_show_gradient3D_triggered");
	showGradient3D();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::showGradient2D()
{
    LOG_MSG("showGradient2D");
    SimpleView2D *mySimpleView2D = new SimpleView2D();
//    QSize size = mySimpleView2D->size();
//    sprintf(msg,"mySimpleView2D size: %d %d",size.height(),size.width());
//    LOG_MSG(msg);
    mySimpleView2D->chooseParameters();
    mySimpleView2D->create();
    mySimpleView2D->show();
    mySimpleView2D->aimCamera();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::showGradient3D()
{
    LOG_MSG("showGradient3D");
    SimpleView3D *mySimpleView3D = new SimpleView3D();
//    QSize size = mySimpleView3D->size();
//    sprintf(msg,"mySimpleView3D size: %d %d",size.height(),size.width());
//    LOG_MSG(msg);
    mySimpleView3D->chooseParameters();
    mySimpleView3D->create();
    mySimpleView3D->show();
    mySimpleView3D->aimCamera();
//    mySimpleView3D->GetRenderWindow()->SetSize(768,768);
}

//======================================================================================================
//------------------------------------------------------------------------------------------------------
// SliderPlus class member definitions
//------------------------------------------------------------------------------------------------------
SliderPlus::SliderPlus(QString aname, double valmin, double valmax, int nval, int iparam, int kwidget)
{
	int i;
    name = aname;
    pindex = iparam;
    windex = kwidget;
    dv = (valmax - valmin)/nval;
	for (i=10; i>-10; i--) {
		if (pow(10.0,i) < dv) {
            dv = pow(10.0,i);
            break;
		}
	}
    i = int(valmin/dv);
    vmin = dv*(i+1);
    int n1 = (int)((valmax - vmin)/dv + 0.5);	//round
	if (n1 > 5*nval) {
        dv = 5*dv;
        n1 = n1/5;
	}
	else if (n1 > 2*nval) {
        dv = 2*dv;
        n1 = n1/2;
	}
    i = int(valmin/dv);
    vmin = dv*i;
    n = n1;
    vmax = vmin + n*dv;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
SliderPlus::~SliderPlus()
{}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
int SliderPlus::val_to_int(double v) {
    int i = (int)((v-vmin)/dv + 0.5);	//round
    return i;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
double SliderPlus::int_to_val(int i) {
    double v = vmin + i*dv;
    return v;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
QString SliderPlus::val_to_str(double v) {
	int i = SliderPlus::val_to_int(v);
	QString vstr = QString::number(int_to_val(i));
    return vstr;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
double SliderPlus::str_to_val(QString vstr) {
    double v = vstr.toDouble();
    return v;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
int SliderPlus::pIndex() {
    return pindex;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
int SliderPlus::wIndex() {
    return windex;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
int SliderPlus::nTicks() {
    return n;
}



//==================================================================================================================
// Code below here is not used
//----------------------------
void MainWindow::closeEvent(QCloseEvent *event)
{
    event->accept();
	/*
    if (maybeSave()) {
        writeSettings();
        event->accept();
    } else {
        event->ignore();
    }
	*/
}

void MainWindow::newFile()
{
    if (maybeSave()) {
        textEdit->clear();
        setCurrentFile("");
    }
}

void MainWindow::open()
{
    if (maybeSave()) {
        QString fileName = QFileDialog::getOpenFileName(this);
        if (!fileName.isEmpty())
            loadFile(fileName);
    }
}

void MainWindow::about()
{
   QMessageBox::about(this, tr("About Application"),
            tr("The <b>Application</b> example demonstrates how to "
               "write modern GUI applications using Qt, with a menu bar, "
               "toolbars, and a status bar."));
}

void MainWindow::documentWasModified()
{
    setWindowModified(textEdit->document()->isModified());
}

void MainWindow::createMenus()
{
    fileMenu = menuBar()->addMenu(tr("&File"));
    fileMenu->addAction(newAct);
    fileMenu->addAction(openAct);
    fileMenu->addAction(saveAct);
    fileMenu->addAction(saveAsAct);
    fileMenu->addSeparator();
    fileMenu->addAction(exitAct);

    editMenu = menuBar()->addMenu(tr("&Edit"));
    editMenu->addAction(cutAct);
    editMenu->addAction(copyAct);
    editMenu->addAction(pasteAct);

    menuBar()->addSeparator();

    helpMenu = menuBar()->addMenu(tr("&Help"));
    helpMenu->addAction(aboutAct);
    helpMenu->addAction(aboutQtAct);
}

void MainWindow::createToolBars()
{
    fileToolBar = addToolBar(tr("File"));
    fileToolBar->addAction(newAct);
    fileToolBar->addAction(openAct);
    fileToolBar->addAction(saveAct);

    editToolBar = addToolBar(tr("Edit"));
    editToolBar->addAction(cutAct);
    editToolBar->addAction(copyAct);
    editToolBar->addAction(pasteAct);
}

void MainWindow::createStatusBar()
{
    statusBar()->showMessage(tr("Ready"));
}

void MainWindow::readSettings()
{
    QSettings settings("Trolltech", "Application Example");
    QPoint pos = settings.value("pos", QPoint(200, 200)).toPoint();
    QSize size = settings.value("size", QSize(400, 400)).toSize();
    resize(size);
    move(pos);
}

void MainWindow::writeSettings()
{
    QSettings settings("Trolltech", "Application Example");
    settings.setValue("pos", pos());
    settings.setValue("size", size());
}

bool MainWindow::maybeSave()
{
    if (textEdit->document()->isModified()) {
        QMessageBox::StandardButton ret;
        ret = QMessageBox::warning(this, tr("Application"),
                     tr("The document has been modified.\n"
                        "Do you want to save your changes?"),
                     QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel);
        if (ret == QMessageBox::Save)
            return save();
        else if (ret == QMessageBox::Cancel)
            return false;
    }
    return true;
}

void MainWindow::loadFile(const QString &fileName)
{
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot read file %1:\n%2.")
                             .arg(fileName)
                             .arg(file.errorString()));
        return;
    }

    QTextStream in(&file);
//#ifndef QT_NO_CURSOR
    QApplication::setOverrideCursor(Qt::WaitCursor);
//#endif
    textEdit->setPlainText(in.readAll());
//#ifndef QT_NO_CURSOR
    QApplication::restoreOverrideCursor();
//#endif

    setCurrentFile(fileName);
    statusBar()->showMessage(tr("File loaded"), 2000);
}
/*
bool MainWindow::saveFile(const QString &fileName)
{
    QFile file(fileName);
    if (!file.open(QFile::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot write file %1:\n%2.")
                             .arg(fileName)
                             .arg(file.errorString()));
        return false;
    }

    QTextStream out(&file);
//#ifndef QT_NO_CURSOR
    QApplication::setOverrideCursor(Qt::WaitCursor);
//#endif
	QString text = "This is a test string";
	text += "\n";
    out << text;
	text = "This is another test string";
	text += "\n";
    out << text;
//#ifndef QT_NO_CURSOR
    QApplication::restoreOverrideCursor();
//#endif

    return true;
}
*/
void MainWindow::setCurrentFile(const QString &fileName)
{
    curFile = fileName;
    textEdit->document()->setModified(false);
    setWindowModified(false);

    QString shownName = curFile;
    if (curFile.isEmpty())
        shownName = "untitled.txt";
    setWindowFilePath(shownName);
}

QString MainWindow::strippedName(const QString &fullFileName)
{
    return QFileInfo(fullFileName).fileName();
}
