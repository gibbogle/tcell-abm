/********************************************************************************
** Form generated from reading UI file 'ABM_GUI.ui'
**
** Created: Mon May 30 11:40:28 2011
**      by: Qt User Interface Compiler version 4.6.3
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_ABM_GUI_H
#define UI_ABM_GUI_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QComboBox>
#include <QtGui/QGridLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QMainWindow>
#include <QtGui/QMdiArea>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QProgressBar>
#include <QtGui/QRadioButton>
#include <QtGui/QSlider>
#include <QtGui/QSpacerItem>
#include <QtGui/QSpinBox>
#include <QtGui/QStackedWidget>
#include <QtGui/QStatusBar>
#include <QtGui/QTabWidget>
#include <QtGui/QTextBrowser>
#include <QtGui/QTextEdit>
#include <QtGui/QToolBar>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "qmylabel.h"
#include "qwt_plot.h"

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *action_saveAs;
    QAction *action_save;
    QAction *action_open_input;
    QAction *action_stop;
    QAction *action_run;
    QAction *action_inputs;
    QAction *action_outputs;
    QAction *action_VTK;
    QAction *action_pause;
    QAction *action_load_results;
    QAction *action_add_graph;
    QAction *action_remove_graph;
    QAction *action_remove_all;
    QAction *action_play_VTK;
    QAction *action_set_speed;
    QAction *action_save_snapshot;
    QWidget *centralwidget;
    QGridLayout *gridLayout_5;
    QStackedWidget *stackedWidget;
    QWidget *page_input;
    QVBoxLayout *verticalLayout;
    QTabWidget *tabs;
    QWidget *tab_T;
    QMyLabel *label_TC_AVIDITY_MEDIAN;
    QSlider *slider_TC_AVIDITY_MEDIAN;
    QLineEdit *line_TC_AVIDITY_MEDIAN;
    QMyLabel *label_TC_AVIDITY_SHAPE;
    QSlider *slider_TC_AVIDITY_SHAPE;
    QLineEdit *line_TC_AVIDITY_SHAPE;
    QWidget *layoutWidget;
    QGridLayout *gridLayout;
    QMyLabel *label_TC_COGNATE_FRACTION;
    QSpacerItem *horizontalSpacer_3;
    QSlider *slider_TC_COGNATE_FRACTION;
    QLineEdit *line_TC_COGNATE_FRACTION;
    QMyLabel *label_TC_STIM_RATE_CONSTANT;
    QSpacerItem *horizontalSpacer_4;
    QSlider *slider_TC_STIM_RATE_CONSTANT;
    QLineEdit *line_TC_STIM_RATE_CONSTANT;
    QMyLabel *label_TC_STIM_HALFLIFE;
    QSpacerItem *horizontalSpacer_5;
    QLineEdit *line_TC_STIM_HALFLIFE;
    QLabel *units_TC_STIM_HALFLIFE;
    QMyLabel *label_MOTILITY_BETA;
    QSpacerItem *horizontalSpacer_10;
    QLineEdit *line_MOTILITY_BETA;
    QMyLabel *label_MOTILITY_RHO;
    QSpacerItem *horizontalSpacer_11;
    QLineEdit *line_MOTILITY_RHO;
    QSlider *slider_TC_STIM_HALFLIFE;
    QSlider *slider_MOTILITY_BETA;
    QSlider *slider_MOTILITY_RHO;
    QwtPlot *qwtPlot_TC_AVIDITY;
    QwtPlot *qwtPlot_DIVIDE1;
    QMyLabel *label_DIVIDE1_MEDIAN;
    QSlider *slider_DIVIDE1_MEDIAN;
    QLineEdit *line_DIVIDE1_MEDIAN;
    QMyLabel *label_DIVIDE1_SHAPE;
    QSlider *slider_DIVIDE1_SHAPE;
    QLineEdit *line_DIVIDE1_SHAPE;
    QLabel *alabel_dist;
    QLineEdit *line_DIVIDE2_SHAPE;
    QSlider *slider_DIVIDE2_MEDIAN;
    QLineEdit *line_DIVIDE2_MEDIAN;
    QMyLabel *label_DIVIDE2_MEDIAN;
    QMyLabel *label_DIVIDE2_SHAPE;
    QSlider *slider_DIVIDE2_SHAPE;
    QwtPlot *qwtPlot_DIVIDE2;
    QWidget *tab_DC;
    QWidget *layoutWidget1;
    QGridLayout *gridLayout_2;
    QMyLabel *label_DC_BIND_DELAY;
    QLabel *units_DC_BIND_DELAY;
    QMyLabel *label_DC_DENS_HALFLIFE;
    QLineEdit *line_DC_DENS_HALFLIFE;
    QMyLabel *label_MAX_TC_BIND;
    QMyLabel *label_MAX_COG_BIND;
    QLabel *units_DC_DENS_HALFLIFE;
    QSpinBox *spin_MAX_TC_BIND;
    QSpinBox *spin_MAX_COG_BIND;
    QLineEdit *line_DC_BIND_DELAY;
    QMyLabel *label_DC_ANTIGEN_SHAPE;
    QwtPlot *qwtPlot_DC_ANTIGEN;
    QLineEdit *line_DC_ANTIGEN_SHAPE;
    QMyLabel *label_DC_ANTIGEN_MEDIAN;
    QLineEdit *line_DC_ANTIGEN_MEDIAN;
    QSlider *slider_DC_ANTIGEN_SHAPE;
    QSlider *slider_DC_ANTIGEN_MEDIAN;
    QMyLabel *label_DC_LIFETIME_SHAPE;
    QwtPlot *qwtPlot_DC_LIFETIME;
    QLineEdit *line_DC_LIFETIME_SHAPE;
    QMyLabel *label_DC_LIFETIME_MEDIAN;
    QLineEdit *line_DC_LIFETIME_MEDIAN;
    QSlider *slider_DC_LIFETIME_SHAPE;
    QSlider *slider_DC_LIFETIME_MEDIAN;
    QLabel *alabel_dist_2;
    QWidget *tab_TCR;
    QWidget *layoutWidget_2;
    QGridLayout *gridLayout_3;
    QMyLabel *label_IL2_THRESHOLD;
    QMyLabel *label_ACTIVATION_THRESHOLD;
    QMyLabel *label_FIRST_DIVISION_THRESHOLD;
    QMyLabel *label_DIVISION_THRESHOLD;
    QMyLabel *label_EXIT_THRESHOLD;
    QLineEdit *line_IL2_THRESHOLD;
    QLineEdit *line_ACTIVATION_THRESHOLD;
    QLineEdit *line_FIRST_DIVISION_THRESHOLD;
    QLineEdit *line_DIVISION_THRESHOLD;
    QLineEdit *line_EXIT_THRESHOLD;
    QMyLabel *label_STIMULATION_LIMIT;
    QLineEdit *line_STIMULATION_LIMIT;
    QLabel *label_40;
    QWidget *tab_run;
    QWidget *layoutWidget2;
    QGridLayout *gridLayout_4;
    QMyLabel *label_NX;
    QSpacerItem *horizontalSpacer_12;
    QSpinBox *spin_NX;
    QMyLabel *label_BLOB_RADIUS;
    QSpacerItem *horizontalSpacer_13;
    QLineEdit *line_BLOB_RADIUS;
    QMyLabel *label_TC_FRACTION;
    QSpacerItem *horizontalSpacer_14;
    QLineEdit *line_TC_FRACTION;
    QSpacerItem *horizontalSpacer_15;
    QLineEdit *line_FLUID_FRACTION;
    QMyLabel *label_DC_RADIUS;
    QSpacerItem *horizontalSpacer_16;
    QLineEdit *line_DC_RADIUS;
    QLabel *units_DC_RADIUS;
    QMyLabel *label_TC_TO_DC;
    QSpacerItem *horizontalSpacer_17;
    QMyLabel *label_DCrate_100k;
    QSpacerItem *horizontalSpacer_18;
    QLineEdit *line_DCrate_100k;
    QMyLabel *label_T_DC1;
    QSpacerItem *horizontalSpacer_19;
    QLineEdit *line_T_DC1;
    QLabel *units_T_DC1;
    QMyLabel *label_T_DC2;
    QSpacerItem *horizontalSpacer_20;
    QLineEdit *line_T_DC2;
    QLabel *units_T_DC2;
    QMyLabel *label_RESIDENCE_TIME;
    QSpacerItem *horizontalSpacer_21;
    QLineEdit *line_RESIDENCE_TIME;
    QLabel *units_RESIDENCE_TIME;
    QMyLabel *label_INFLAMM_DAYS1;
    QSpacerItem *horizontalSpacer_22;
    QLineEdit *line_INFLAMM_DAYS1;
    QMyLabel *label_INFLAMM_DAYS2;
    QSpacerItem *horizontalSpacer_23;
    QLineEdit *line_INFLAMM_DAYS2;
    QMyLabel *label_INFLAMM_LEVEL;
    QSpacerItem *horizontalSpacer_24;
    QLineEdit *line_INFLAMM_LEVEL;
    QMyLabel *label_EXIT_REGION;
    QSpacerItem *horizontalSpacer_26;
    QComboBox *comb_EXIT_REGION;
    QMyLabel *label_CHEMO_RADIUS;
    QSpacerItem *horizontalSpacer_27;
    QLineEdit *line_CHEMO_RADIUS;
    QMyLabel *label_CHEMO_K_EXIT;
    QSpacerItem *horizontalSpacer_28;
    QLineEdit *line_CHEMO_K_EXIT;
    QMyLabel *label_NDAYS;
    QSpacerItem *horizontalSpacer_29;
    QMyLabel *label_SEED1;
    QSpacerItem *horizontalSpacer_30;
    QSpinBox *spin_SEED1;
    QMyLabel *label_SEED2;
    QSpacerItem *horizontalSpacer_31;
    QSpinBox *spin_SEED2;
    QLineEdit *line_TC_TO_DC;
    QLineEdit *line_NDAYS;
    QMyLabel *label_FLUID_FRACTION;
    QMyLabel *label_NT_ANIMATION;
    QSpacerItem *horizontalSpacer_32;
    QSpinBox *spin_NT_ANIMATION;
    QMyLabel *label_CHEMO_K_DC;
    QLineEdit *line_CHEMO_K_DC;
    QSpacerItem *horizontalSpacer_33;
    QLabel *units_CHEMO_RADIUS;
    QLabel *units_NDAYS;
    QLabel *units_INFLAMM1;
    QLabel *units_INFLAMM2;
    QLabel *units_BLOB_RADIUS;
    QMyLabel *label_NCPU;
    QSpacerItem *horizontalSpacer_34;
    QSpinBox *spin_NCPU;
    QMyLabel *label_T_DC_INJECTION;
    QSpacerItem *horizontalSpacer_35;
    QLineEdit *line_T_DC_INJECTION;
    QLabel *units_T_DC_INJECTION;
    QCheckBox *cbox_savepos;
    QCheckBox *cbox_IN_VITRO;
    QWidget *layoutWidget3;
    QGridLayout *gridLayout_6;
    QMyLabel *label_IV_WELL_DIAMETER;
    QSpacerItem *horizontalSpacer;
    QLineEdit *line_IV_WELL_DIAMETER;
    QLabel *units_IV_WELL_DIAMETER;
    QMyLabel *label_IV_NTCELLS;
    QSpacerItem *horizontalSpacer_2;
    QLineEdit *line_IV_NTCELLS;
    QMyLabel *label_IV_COGNATE_FRACTION;
    QSpacerItem *horizontalSpacer_6;
    QLineEdit *line_IV_COGNATE_FRACTION;
    QCheckBox *cbox_IV_SHOW_NONCOGNATE;
    QCheckBox *cbox_DC_INJECTION;
    QCheckBox *cbox_USE_TRAFFIC;
    QRadioButton *rbut_SPECIES_1;
    QRadioButton *rbut_SPECIES_0;
    QCheckBox *cbox_USE_EXIT_CHEMOTAXIS;
    QCheckBox *cbox_USE_DC_CHEMOTAXIS;
    QCheckBox *cbox_COMPUTED_OUTFLOW;
    QMyLabel *label_INPUT_FILE;
    QLineEdit *text_INPUT_FILE;
    QLabel *label_input;
    QTextEdit *text_more;
    QWidget *page_output;
    QVBoxLayout *verticalLayout_2;
    QMdiArea *mdiArea;
    QTextBrowser *box_outputLog;
    QWidget *page_3D;
    QMdiArea *mdiArea_VTK;
    QProgressBar *progressBar;
    QLabel *label_hour;
    QLabel *hour_display;
    QMenuBar *menubar;
    QMenu *menuFile;
    QMenu *menuEdit;
    QMenu *menuABM;
    QMenu *menuGraphs;
    QMenu *menuPlayer;
    QMenu *menuSnapshot;
    QStatusBar *statusbar;
    QToolBar *toolBar;
    QButtonGroup *buttonGroup_SPECIES;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(1319, 977);
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/ABM.png"), QSize(), QIcon::Normal, QIcon::Off);
        MainWindow->setWindowIcon(icon);
        MainWindow->setStyleSheet(QString::fromUtf8("#tabs{\n"
"color: white;\n"
"background-color: QLinearGradient( x1: 0, y1: 0, x2: 0, y2: 1, stop: 0 #88d, stop: 0.1 #99e, stop: 0.49 #77c, stop: 0.5 #66b, stop: 1 #77c);\n"
"border-width: 1px;\n"
"border-color: #339;\n"
"border-style: solid;\n"
"border-radius: 15;\n"
"}\n"
"\n"
"#text_more{\n"
"padding: 8px\n"
"}"));
        action_saveAs = new QAction(MainWindow);
        action_saveAs->setObjectName(QString::fromUtf8("action_saveAs"));
        action_save = new QAction(MainWindow);
        action_save->setObjectName(QString::fromUtf8("action_save"));
        action_open_input = new QAction(MainWindow);
        action_open_input->setObjectName(QString::fromUtf8("action_open_input"));
        action_stop = new QAction(MainWindow);
        action_stop->setObjectName(QString::fromUtf8("action_stop"));
        action_stop->setCheckable(false);
        action_stop->setEnabled(true);
        QIcon icon1;
        icon1.addFile(QString::fromUtf8(":/icons/001_29.png"), QSize(), QIcon::Normal, QIcon::Off);
        action_stop->setIcon(icon1);
        action_run = new QAction(MainWindow);
        action_run->setObjectName(QString::fromUtf8("action_run"));
        QIcon icon2;
        icon2.addFile(QString::fromUtf8(":/icons/001_59.png"), QSize(), QIcon::Normal, QIcon::Off);
        action_run->setIcon(icon2);
        action_inputs = new QAction(MainWindow);
        action_inputs->setObjectName(QString::fromUtf8("action_inputs"));
        QIcon icon3;
        icon3.addFile(QString::fromUtf8(":/icons/001_45.png"), QSize(), QIcon::Normal, QIcon::Off);
        action_inputs->setIcon(icon3);
        action_outputs = new QAction(MainWindow);
        action_outputs->setObjectName(QString::fromUtf8("action_outputs"));
        QIcon icon4;
        icon4.addFile(QString::fromUtf8(":/icons/Display2.png"), QSize(), QIcon::Normal, QIcon::Off);
        action_outputs->setIcon(icon4);
        action_VTK = new QAction(MainWindow);
        action_VTK->setObjectName(QString::fromUtf8("action_VTK"));
        action_VTK->setEnabled(true);
        QIcon icon5;
        icon5.addFile(QString::fromUtf8(":/icons/cell2.png"), QSize(), QIcon::Normal, QIcon::Off);
        action_VTK->setIcon(icon5);
        action_pause = new QAction(MainWindow);
        action_pause->setObjectName(QString::fromUtf8("action_pause"));
        QIcon icon6;
        icon6.addFile(QString::fromUtf8(":/icons/001_07.png"), QSize(), QIcon::Normal, QIcon::Off);
        action_pause->setIcon(icon6);
        action_load_results = new QAction(MainWindow);
        action_load_results->setObjectName(QString::fromUtf8("action_load_results"));
        action_add_graph = new QAction(MainWindow);
        action_add_graph->setObjectName(QString::fromUtf8("action_add_graph"));
        action_remove_graph = new QAction(MainWindow);
        action_remove_graph->setObjectName(QString::fromUtf8("action_remove_graph"));
        action_remove_all = new QAction(MainWindow);
        action_remove_all->setObjectName(QString::fromUtf8("action_remove_all"));
        action_play_VTK = new QAction(MainWindow);
        action_play_VTK->setObjectName(QString::fromUtf8("action_play_VTK"));
        action_set_speed = new QAction(MainWindow);
        action_set_speed->setObjectName(QString::fromUtf8("action_set_speed"));
        action_save_snapshot = new QAction(MainWindow);
        action_save_snapshot->setObjectName(QString::fromUtf8("action_save_snapshot"));
        centralwidget = new QWidget(MainWindow);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        gridLayout_5 = new QGridLayout(centralwidget);
        gridLayout_5->setObjectName(QString::fromUtf8("gridLayout_5"));
        stackedWidget = new QStackedWidget(centralwidget);
        stackedWidget->setObjectName(QString::fromUtf8("stackedWidget"));
        page_input = new QWidget();
        page_input->setObjectName(QString::fromUtf8("page_input"));
        verticalLayout = new QVBoxLayout(page_input);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setContentsMargins(-1, 0, -1, -1);
        tabs = new QTabWidget(page_input);
        tabs->setObjectName(QString::fromUtf8("tabs"));
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Expanding);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(tabs->sizePolicy().hasHeightForWidth());
        tabs->setSizePolicy(sizePolicy);
        tabs->setMinimumSize(QSize(0, 0));
        tabs->setAutoFillBackground(false);
        tabs->setTabPosition(QTabWidget::West);
        tabs->setElideMode(Qt::ElideNone);
        tab_T = new QWidget();
        tab_T->setObjectName(QString::fromUtf8("tab_T"));
        label_TC_AVIDITY_MEDIAN = new QMyLabel(tab_T);
        label_TC_AVIDITY_MEDIAN->setObjectName(QString::fromUtf8("label_TC_AVIDITY_MEDIAN"));
        label_TC_AVIDITY_MEDIAN->setGeometry(QRect(760, 80, 61, 20));
        QSizePolicy sizePolicy1(QSizePolicy::Fixed, QSizePolicy::Preferred);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(label_TC_AVIDITY_MEDIAN->sizePolicy().hasHeightForWidth());
        label_TC_AVIDITY_MEDIAN->setSizePolicy(sizePolicy1);
        label_TC_AVIDITY_MEDIAN->setCursor(QCursor(Qt::ArrowCursor));
        label_TC_AVIDITY_MEDIAN->setMouseTracking(true);
        label_TC_AVIDITY_MEDIAN->setWordWrap(false);
        label_TC_AVIDITY_MEDIAN->setTextInteractionFlags(Qt::LinksAccessibleByMouse);
        slider_TC_AVIDITY_MEDIAN = new QSlider(tab_T);
        slider_TC_AVIDITY_MEDIAN->setObjectName(QString::fromUtf8("slider_TC_AVIDITY_MEDIAN"));
        slider_TC_AVIDITY_MEDIAN->setGeometry(QRect(823, 80, 61, 20));
        slider_TC_AVIDITY_MEDIAN->setOrientation(Qt::Horizontal);
        line_TC_AVIDITY_MEDIAN = new QLineEdit(tab_T);
        line_TC_AVIDITY_MEDIAN->setObjectName(QString::fromUtf8("line_TC_AVIDITY_MEDIAN"));
        line_TC_AVIDITY_MEDIAN->setGeometry(QRect(910, 80, 51, 20));
        line_TC_AVIDITY_MEDIAN->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_TC_AVIDITY_SHAPE = new QMyLabel(tab_T);
        label_TC_AVIDITY_SHAPE->setObjectName(QString::fromUtf8("label_TC_AVIDITY_SHAPE"));
        label_TC_AVIDITY_SHAPE->setGeometry(QRect(760, 120, 51, 20));
        sizePolicy1.setHeightForWidth(label_TC_AVIDITY_SHAPE->sizePolicy().hasHeightForWidth());
        label_TC_AVIDITY_SHAPE->setSizePolicy(sizePolicy1);
        label_TC_AVIDITY_SHAPE->setWordWrap(false);
        slider_TC_AVIDITY_SHAPE = new QSlider(tab_T);
        slider_TC_AVIDITY_SHAPE->setObjectName(QString::fromUtf8("slider_TC_AVIDITY_SHAPE"));
        slider_TC_AVIDITY_SHAPE->setGeometry(QRect(823, 120, 61, 20));
        slider_TC_AVIDITY_SHAPE->setOrientation(Qt::Horizontal);
        line_TC_AVIDITY_SHAPE = new QLineEdit(tab_T);
        line_TC_AVIDITY_SHAPE->setObjectName(QString::fromUtf8("line_TC_AVIDITY_SHAPE"));
        line_TC_AVIDITY_SHAPE->setGeometry(QRect(910, 120, 51, 20));
        line_TC_AVIDITY_SHAPE->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        layoutWidget = new QWidget(tab_T);
        layoutWidget->setObjectName(QString::fromUtf8("layoutWidget"));
        layoutWidget->setGeometry(QRect(0, 0, 391, 251));
        gridLayout = new QGridLayout(layoutWidget);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        gridLayout->setContentsMargins(0, 0, 0, 0);
        label_TC_COGNATE_FRACTION = new QMyLabel(layoutWidget);
        label_TC_COGNATE_FRACTION->setObjectName(QString::fromUtf8("label_TC_COGNATE_FRACTION"));
        sizePolicy1.setHeightForWidth(label_TC_COGNATE_FRACTION->sizePolicy().hasHeightForWidth());
        label_TC_COGNATE_FRACTION->setSizePolicy(sizePolicy1);
        label_TC_COGNATE_FRACTION->setWordWrap(false);

        gridLayout->addWidget(label_TC_COGNATE_FRACTION, 0, 0, 1, 1);

        horizontalSpacer_3 = new QSpacerItem(13, 17, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer_3, 0, 1, 1, 1);

        slider_TC_COGNATE_FRACTION = new QSlider(layoutWidget);
        slider_TC_COGNATE_FRACTION->setObjectName(QString::fromUtf8("slider_TC_COGNATE_FRACTION"));
        slider_TC_COGNATE_FRACTION->setOrientation(Qt::Horizontal);

        gridLayout->addWidget(slider_TC_COGNATE_FRACTION, 0, 2, 1, 1);

        line_TC_COGNATE_FRACTION = new QLineEdit(layoutWidget);
        line_TC_COGNATE_FRACTION->setObjectName(QString::fromUtf8("line_TC_COGNATE_FRACTION"));
        line_TC_COGNATE_FRACTION->setMaximumSize(QSize(90, 16777215));
        line_TC_COGNATE_FRACTION->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(line_TC_COGNATE_FRACTION, 0, 3, 1, 1);

        label_TC_STIM_RATE_CONSTANT = new QMyLabel(layoutWidget);
        label_TC_STIM_RATE_CONSTANT->setObjectName(QString::fromUtf8("label_TC_STIM_RATE_CONSTANT"));
        sizePolicy1.setHeightForWidth(label_TC_STIM_RATE_CONSTANT->sizePolicy().hasHeightForWidth());
        label_TC_STIM_RATE_CONSTANT->setSizePolicy(sizePolicy1);
        label_TC_STIM_RATE_CONSTANT->setMouseTracking(false);
        label_TC_STIM_RATE_CONSTANT->setWordWrap(false);

        gridLayout->addWidget(label_TC_STIM_RATE_CONSTANT, 1, 0, 1, 1);

        horizontalSpacer_4 = new QSpacerItem(13, 17, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer_4, 1, 1, 1, 1);

        slider_TC_STIM_RATE_CONSTANT = new QSlider(layoutWidget);
        slider_TC_STIM_RATE_CONSTANT->setObjectName(QString::fromUtf8("slider_TC_STIM_RATE_CONSTANT"));
        slider_TC_STIM_RATE_CONSTANT->setOrientation(Qt::Horizontal);

        gridLayout->addWidget(slider_TC_STIM_RATE_CONSTANT, 1, 2, 1, 1);

        line_TC_STIM_RATE_CONSTANT = new QLineEdit(layoutWidget);
        line_TC_STIM_RATE_CONSTANT->setObjectName(QString::fromUtf8("line_TC_STIM_RATE_CONSTANT"));
        line_TC_STIM_RATE_CONSTANT->setMaximumSize(QSize(90, 16777215));
        line_TC_STIM_RATE_CONSTANT->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(line_TC_STIM_RATE_CONSTANT, 1, 3, 1, 1);

        label_TC_STIM_HALFLIFE = new QMyLabel(layoutWidget);
        label_TC_STIM_HALFLIFE->setObjectName(QString::fromUtf8("label_TC_STIM_HALFLIFE"));
        sizePolicy1.setHeightForWidth(label_TC_STIM_HALFLIFE->sizePolicy().hasHeightForWidth());
        label_TC_STIM_HALFLIFE->setSizePolicy(sizePolicy1);
        label_TC_STIM_HALFLIFE->setMouseTracking(false);
        label_TC_STIM_HALFLIFE->setWordWrap(false);

        gridLayout->addWidget(label_TC_STIM_HALFLIFE, 2, 0, 1, 1);

        horizontalSpacer_5 = new QSpacerItem(13, 17, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer_5, 2, 1, 1, 1);

        line_TC_STIM_HALFLIFE = new QLineEdit(layoutWidget);
        line_TC_STIM_HALFLIFE->setObjectName(QString::fromUtf8("line_TC_STIM_HALFLIFE"));
        line_TC_STIM_HALFLIFE->setMaximumSize(QSize(90, 16777215));
        line_TC_STIM_HALFLIFE->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(line_TC_STIM_HALFLIFE, 2, 3, 1, 1);

        units_TC_STIM_HALFLIFE = new QLabel(layoutWidget);
        units_TC_STIM_HALFLIFE->setObjectName(QString::fromUtf8("units_TC_STIM_HALFLIFE"));

        gridLayout->addWidget(units_TC_STIM_HALFLIFE, 2, 4, 1, 1);

        label_MOTILITY_BETA = new QMyLabel(layoutWidget);
        label_MOTILITY_BETA->setObjectName(QString::fromUtf8("label_MOTILITY_BETA"));
        sizePolicy1.setHeightForWidth(label_MOTILITY_BETA->sizePolicy().hasHeightForWidth());
        label_MOTILITY_BETA->setSizePolicy(sizePolicy1);
        label_MOTILITY_BETA->setMouseTracking(false);
        label_MOTILITY_BETA->setWordWrap(false);

        gridLayout->addWidget(label_MOTILITY_BETA, 3, 0, 1, 1);

        horizontalSpacer_10 = new QSpacerItem(13, 17, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer_10, 3, 1, 1, 1);

        line_MOTILITY_BETA = new QLineEdit(layoutWidget);
        line_MOTILITY_BETA->setObjectName(QString::fromUtf8("line_MOTILITY_BETA"));
        line_MOTILITY_BETA->setMaximumSize(QSize(90, 16777215));
        line_MOTILITY_BETA->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(line_MOTILITY_BETA, 3, 3, 1, 1);

        label_MOTILITY_RHO = new QMyLabel(layoutWidget);
        label_MOTILITY_RHO->setObjectName(QString::fromUtf8("label_MOTILITY_RHO"));
        sizePolicy1.setHeightForWidth(label_MOTILITY_RHO->sizePolicy().hasHeightForWidth());
        label_MOTILITY_RHO->setSizePolicy(sizePolicy1);
        label_MOTILITY_RHO->setMouseTracking(false);
        label_MOTILITY_RHO->setWordWrap(false);

        gridLayout->addWidget(label_MOTILITY_RHO, 4, 0, 1, 1);

        horizontalSpacer_11 = new QSpacerItem(13, 17, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer_11, 4, 1, 1, 1);

        line_MOTILITY_RHO = new QLineEdit(layoutWidget);
        line_MOTILITY_RHO->setObjectName(QString::fromUtf8("line_MOTILITY_RHO"));
        line_MOTILITY_RHO->setMaximumSize(QSize(90, 16777215));
        line_MOTILITY_RHO->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(line_MOTILITY_RHO, 4, 3, 1, 1);

        slider_TC_STIM_HALFLIFE = new QSlider(layoutWidget);
        slider_TC_STIM_HALFLIFE->setObjectName(QString::fromUtf8("slider_TC_STIM_HALFLIFE"));
        slider_TC_STIM_HALFLIFE->setOrientation(Qt::Horizontal);

        gridLayout->addWidget(slider_TC_STIM_HALFLIFE, 2, 2, 1, 1);

        slider_MOTILITY_BETA = new QSlider(layoutWidget);
        slider_MOTILITY_BETA->setObjectName(QString::fromUtf8("slider_MOTILITY_BETA"));
        slider_MOTILITY_BETA->setOrientation(Qt::Horizontal);

        gridLayout->addWidget(slider_MOTILITY_BETA, 3, 2, 1, 1);

        slider_MOTILITY_RHO = new QSlider(layoutWidget);
        slider_MOTILITY_RHO->setObjectName(QString::fromUtf8("slider_MOTILITY_RHO"));
        slider_MOTILITY_RHO->setOrientation(Qt::Horizontal);

        gridLayout->addWidget(slider_MOTILITY_RHO, 4, 2, 1, 1);

        qwtPlot_TC_AVIDITY = new QwtPlot(tab_T);
        qwtPlot_TC_AVIDITY->setObjectName(QString::fromUtf8("qwtPlot_TC_AVIDITY"));
        qwtPlot_TC_AVIDITY->setGeometry(QRect(420, 50, 300, 170));
        qwtPlot_DIVIDE1 = new QwtPlot(tab_T);
        qwtPlot_DIVIDE1->setObjectName(QString::fromUtf8("qwtPlot_DIVIDE1"));
        qwtPlot_DIVIDE1->setGeometry(QRect(420, 250, 300, 170));
        label_DIVIDE1_MEDIAN = new QMyLabel(tab_T);
        label_DIVIDE1_MEDIAN->setObjectName(QString::fromUtf8("label_DIVIDE1_MEDIAN"));
        label_DIVIDE1_MEDIAN->setGeometry(QRect(760, 280, 61, 20));
        sizePolicy1.setHeightForWidth(label_DIVIDE1_MEDIAN->sizePolicy().hasHeightForWidth());
        label_DIVIDE1_MEDIAN->setSizePolicy(sizePolicy1);
        label_DIVIDE1_MEDIAN->setMouseTracking(false);
        label_DIVIDE1_MEDIAN->setWordWrap(false);
        slider_DIVIDE1_MEDIAN = new QSlider(tab_T);
        slider_DIVIDE1_MEDIAN->setObjectName(QString::fromUtf8("slider_DIVIDE1_MEDIAN"));
        slider_DIVIDE1_MEDIAN->setGeometry(QRect(820, 280, 61, 16));
        slider_DIVIDE1_MEDIAN->setOrientation(Qt::Horizontal);
        line_DIVIDE1_MEDIAN = new QLineEdit(tab_T);
        line_DIVIDE1_MEDIAN->setObjectName(QString::fromUtf8("line_DIVIDE1_MEDIAN"));
        line_DIVIDE1_MEDIAN->setGeometry(QRect(910, 280, 51, 20));
        line_DIVIDE1_MEDIAN->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_DIVIDE1_SHAPE = new QMyLabel(tab_T);
        label_DIVIDE1_SHAPE->setObjectName(QString::fromUtf8("label_DIVIDE1_SHAPE"));
        label_DIVIDE1_SHAPE->setGeometry(QRect(760, 320, 61, 20));
        sizePolicy1.setHeightForWidth(label_DIVIDE1_SHAPE->sizePolicy().hasHeightForWidth());
        label_DIVIDE1_SHAPE->setSizePolicy(sizePolicy1);
        label_DIVIDE1_SHAPE->setMouseTracking(false);
        label_DIVIDE1_SHAPE->setWordWrap(false);
        slider_DIVIDE1_SHAPE = new QSlider(tab_T);
        slider_DIVIDE1_SHAPE->setObjectName(QString::fromUtf8("slider_DIVIDE1_SHAPE"));
        slider_DIVIDE1_SHAPE->setGeometry(QRect(820, 320, 61, 16));
        slider_DIVIDE1_SHAPE->setOrientation(Qt::Horizontal);
        line_DIVIDE1_SHAPE = new QLineEdit(tab_T);
        line_DIVIDE1_SHAPE->setObjectName(QString::fromUtf8("line_DIVIDE1_SHAPE"));
        line_DIVIDE1_SHAPE->setGeometry(QRect(910, 320, 51, 20));
        line_DIVIDE1_SHAPE->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        alabel_dist = new QLabel(tab_T);
        alabel_dist->setObjectName(QString::fromUtf8("alabel_dist"));
        alabel_dist->setGeometry(QRect(550, 10, 291, 20));
        QFont font;
        font.setPointSize(14);
        font.setBold(true);
        font.setWeight(75);
        alabel_dist->setFont(font);
        line_DIVIDE2_SHAPE = new QLineEdit(tab_T);
        line_DIVIDE2_SHAPE->setObjectName(QString::fromUtf8("line_DIVIDE2_SHAPE"));
        line_DIVIDE2_SHAPE->setGeometry(QRect(910, 540, 51, 20));
        line_DIVIDE2_SHAPE->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        slider_DIVIDE2_MEDIAN = new QSlider(tab_T);
        slider_DIVIDE2_MEDIAN->setObjectName(QString::fromUtf8("slider_DIVIDE2_MEDIAN"));
        slider_DIVIDE2_MEDIAN->setGeometry(QRect(820, 500, 61, 16));
        slider_DIVIDE2_MEDIAN->setOrientation(Qt::Horizontal);
        line_DIVIDE2_MEDIAN = new QLineEdit(tab_T);
        line_DIVIDE2_MEDIAN->setObjectName(QString::fromUtf8("line_DIVIDE2_MEDIAN"));
        line_DIVIDE2_MEDIAN->setGeometry(QRect(910, 500, 51, 20));
        line_DIVIDE2_MEDIAN->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_DIVIDE2_MEDIAN = new QMyLabel(tab_T);
        label_DIVIDE2_MEDIAN->setObjectName(QString::fromUtf8("label_DIVIDE2_MEDIAN"));
        label_DIVIDE2_MEDIAN->setGeometry(QRect(760, 500, 61, 20));
        sizePolicy1.setHeightForWidth(label_DIVIDE2_MEDIAN->sizePolicy().hasHeightForWidth());
        label_DIVIDE2_MEDIAN->setSizePolicy(sizePolicy1);
        label_DIVIDE2_MEDIAN->setMouseTracking(false);
        label_DIVIDE2_MEDIAN->setWordWrap(false);
        label_DIVIDE2_SHAPE = new QMyLabel(tab_T);
        label_DIVIDE2_SHAPE->setObjectName(QString::fromUtf8("label_DIVIDE2_SHAPE"));
        label_DIVIDE2_SHAPE->setGeometry(QRect(760, 540, 61, 20));
        sizePolicy1.setHeightForWidth(label_DIVIDE2_SHAPE->sizePolicy().hasHeightForWidth());
        label_DIVIDE2_SHAPE->setSizePolicy(sizePolicy1);
        label_DIVIDE2_SHAPE->setMouseTracking(false);
        label_DIVIDE2_SHAPE->setWordWrap(false);
        slider_DIVIDE2_SHAPE = new QSlider(tab_T);
        slider_DIVIDE2_SHAPE->setObjectName(QString::fromUtf8("slider_DIVIDE2_SHAPE"));
        slider_DIVIDE2_SHAPE->setGeometry(QRect(820, 540, 61, 16));
        slider_DIVIDE2_SHAPE->setOrientation(Qt::Horizontal);
        qwtPlot_DIVIDE2 = new QwtPlot(tab_T);
        qwtPlot_DIVIDE2->setObjectName(QString::fromUtf8("qwtPlot_DIVIDE2"));
        qwtPlot_DIVIDE2->setGeometry(QRect(420, 450, 300, 170));
        QFont font1;
        font1.setFamily(QString::fromUtf8("Arial"));
        font1.setPointSize(4);
        qwtPlot_DIVIDE2->setFont(font1);
        tabs->addTab(tab_T, QString());
        tab_DC = new QWidget();
        tab_DC->setObjectName(QString::fromUtf8("tab_DC"));
        layoutWidget1 = new QWidget(tab_DC);
        layoutWidget1->setObjectName(QString::fromUtf8("layoutWidget1"));
        layoutWidget1->setGeometry(QRect(10, 20, 361, 181));
        gridLayout_2 = new QGridLayout(layoutWidget1);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        gridLayout_2->setContentsMargins(0, 0, 0, 0);
        label_DC_BIND_DELAY = new QMyLabel(layoutWidget1);
        label_DC_BIND_DELAY->setObjectName(QString::fromUtf8("label_DC_BIND_DELAY"));
        label_DC_BIND_DELAY->setMouseTracking(false);
        label_DC_BIND_DELAY->setWordWrap(false);

        gridLayout_2->addWidget(label_DC_BIND_DELAY, 0, 0, 1, 2);

        units_DC_BIND_DELAY = new QLabel(layoutWidget1);
        units_DC_BIND_DELAY->setObjectName(QString::fromUtf8("units_DC_BIND_DELAY"));

        gridLayout_2->addWidget(units_DC_BIND_DELAY, 0, 3, 1, 1);

        label_DC_DENS_HALFLIFE = new QMyLabel(layoutWidget1);
        label_DC_DENS_HALFLIFE->setObjectName(QString::fromUtf8("label_DC_DENS_HALFLIFE"));
        label_DC_DENS_HALFLIFE->setMouseTracking(false);
        label_DC_DENS_HALFLIFE->setWordWrap(false);

        gridLayout_2->addWidget(label_DC_DENS_HALFLIFE, 1, 0, 1, 2);

        line_DC_DENS_HALFLIFE = new QLineEdit(layoutWidget1);
        line_DC_DENS_HALFLIFE->setObjectName(QString::fromUtf8("line_DC_DENS_HALFLIFE"));
        line_DC_DENS_HALFLIFE->setMaximumSize(QSize(90, 16777215));
        line_DC_DENS_HALFLIFE->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(line_DC_DENS_HALFLIFE, 1, 2, 1, 1);

        label_MAX_TC_BIND = new QMyLabel(layoutWidget1);
        label_MAX_TC_BIND->setObjectName(QString::fromUtf8("label_MAX_TC_BIND"));
        label_MAX_TC_BIND->setMouseTracking(false);
        label_MAX_TC_BIND->setWordWrap(false);

        gridLayout_2->addWidget(label_MAX_TC_BIND, 2, 0, 1, 2);

        label_MAX_COG_BIND = new QMyLabel(layoutWidget1);
        label_MAX_COG_BIND->setObjectName(QString::fromUtf8("label_MAX_COG_BIND"));
        label_MAX_COG_BIND->setMouseTracking(false);
        label_MAX_COG_BIND->setWordWrap(false);

        gridLayout_2->addWidget(label_MAX_COG_BIND, 3, 0, 1, 2);

        units_DC_DENS_HALFLIFE = new QLabel(layoutWidget1);
        units_DC_DENS_HALFLIFE->setObjectName(QString::fromUtf8("units_DC_DENS_HALFLIFE"));

        gridLayout_2->addWidget(units_DC_DENS_HALFLIFE, 1, 3, 1, 1);

        spin_MAX_TC_BIND = new QSpinBox(layoutWidget1);
        spin_MAX_TC_BIND->setObjectName(QString::fromUtf8("spin_MAX_TC_BIND"));
        spin_MAX_TC_BIND->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(spin_MAX_TC_BIND, 2, 2, 1, 1);

        spin_MAX_COG_BIND = new QSpinBox(layoutWidget1);
        spin_MAX_COG_BIND->setObjectName(QString::fromUtf8("spin_MAX_COG_BIND"));
        spin_MAX_COG_BIND->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(spin_MAX_COG_BIND, 3, 2, 1, 1);

        line_DC_BIND_DELAY = new QLineEdit(layoutWidget1);
        line_DC_BIND_DELAY->setObjectName(QString::fromUtf8("line_DC_BIND_DELAY"));
        line_DC_BIND_DELAY->setMaximumSize(QSize(90, 16777215));
        line_DC_BIND_DELAY->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(line_DC_BIND_DELAY, 0, 2, 1, 1);

        label_DC_ANTIGEN_SHAPE = new QMyLabel(tab_DC);
        label_DC_ANTIGEN_SHAPE->setObjectName(QString::fromUtf8("label_DC_ANTIGEN_SHAPE"));
        label_DC_ANTIGEN_SHAPE->setGeometry(QRect(790, 150, 61, 20));
        sizePolicy1.setHeightForWidth(label_DC_ANTIGEN_SHAPE->sizePolicy().hasHeightForWidth());
        label_DC_ANTIGEN_SHAPE->setSizePolicy(sizePolicy1);
        label_DC_ANTIGEN_SHAPE->setMouseTracking(false);
        label_DC_ANTIGEN_SHAPE->setWordWrap(false);
        qwtPlot_DC_ANTIGEN = new QwtPlot(tab_DC);
        qwtPlot_DC_ANTIGEN->setObjectName(QString::fromUtf8("qwtPlot_DC_ANTIGEN"));
        qwtPlot_DC_ANTIGEN->setGeometry(QRect(430, 70, 300, 200));
        line_DC_ANTIGEN_SHAPE = new QLineEdit(tab_DC);
        line_DC_ANTIGEN_SHAPE->setObjectName(QString::fromUtf8("line_DC_ANTIGEN_SHAPE"));
        line_DC_ANTIGEN_SHAPE->setGeometry(QRect(940, 150, 51, 20));
        line_DC_ANTIGEN_SHAPE->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_DC_ANTIGEN_MEDIAN = new QMyLabel(tab_DC);
        label_DC_ANTIGEN_MEDIAN->setObjectName(QString::fromUtf8("label_DC_ANTIGEN_MEDIAN"));
        label_DC_ANTIGEN_MEDIAN->setGeometry(QRect(790, 110, 61, 20));
        sizePolicy1.setHeightForWidth(label_DC_ANTIGEN_MEDIAN->sizePolicy().hasHeightForWidth());
        label_DC_ANTIGEN_MEDIAN->setSizePolicy(sizePolicy1);
        label_DC_ANTIGEN_MEDIAN->setMouseTracking(false);
        label_DC_ANTIGEN_MEDIAN->setWordWrap(false);
        line_DC_ANTIGEN_MEDIAN = new QLineEdit(tab_DC);
        line_DC_ANTIGEN_MEDIAN->setObjectName(QString::fromUtf8("line_DC_ANTIGEN_MEDIAN"));
        line_DC_ANTIGEN_MEDIAN->setGeometry(QRect(940, 110, 51, 20));
        line_DC_ANTIGEN_MEDIAN->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        slider_DC_ANTIGEN_SHAPE = new QSlider(tab_DC);
        slider_DC_ANTIGEN_SHAPE->setObjectName(QString::fromUtf8("slider_DC_ANTIGEN_SHAPE"));
        slider_DC_ANTIGEN_SHAPE->setGeometry(QRect(850, 150, 61, 16));
        slider_DC_ANTIGEN_SHAPE->setOrientation(Qt::Horizontal);
        slider_DC_ANTIGEN_MEDIAN = new QSlider(tab_DC);
        slider_DC_ANTIGEN_MEDIAN->setObjectName(QString::fromUtf8("slider_DC_ANTIGEN_MEDIAN"));
        slider_DC_ANTIGEN_MEDIAN->setGeometry(QRect(850, 110, 61, 16));
        slider_DC_ANTIGEN_MEDIAN->setOrientation(Qt::Horizontal);
        label_DC_LIFETIME_SHAPE = new QMyLabel(tab_DC);
        label_DC_LIFETIME_SHAPE->setObjectName(QString::fromUtf8("label_DC_LIFETIME_SHAPE"));
        label_DC_LIFETIME_SHAPE->setGeometry(QRect(790, 430, 61, 20));
        sizePolicy1.setHeightForWidth(label_DC_LIFETIME_SHAPE->sizePolicy().hasHeightForWidth());
        label_DC_LIFETIME_SHAPE->setSizePolicy(sizePolicy1);
        label_DC_LIFETIME_SHAPE->setMouseTracking(false);
        label_DC_LIFETIME_SHAPE->setWordWrap(false);
        qwtPlot_DC_LIFETIME = new QwtPlot(tab_DC);
        qwtPlot_DC_LIFETIME->setObjectName(QString::fromUtf8("qwtPlot_DC_LIFETIME"));
        qwtPlot_DC_LIFETIME->setGeometry(QRect(430, 350, 300, 200));
        line_DC_LIFETIME_SHAPE = new QLineEdit(tab_DC);
        line_DC_LIFETIME_SHAPE->setObjectName(QString::fromUtf8("line_DC_LIFETIME_SHAPE"));
        line_DC_LIFETIME_SHAPE->setGeometry(QRect(940, 430, 51, 20));
        line_DC_LIFETIME_SHAPE->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_DC_LIFETIME_MEDIAN = new QMyLabel(tab_DC);
        label_DC_LIFETIME_MEDIAN->setObjectName(QString::fromUtf8("label_DC_LIFETIME_MEDIAN"));
        label_DC_LIFETIME_MEDIAN->setGeometry(QRect(790, 390, 61, 20));
        sizePolicy1.setHeightForWidth(label_DC_LIFETIME_MEDIAN->sizePolicy().hasHeightForWidth());
        label_DC_LIFETIME_MEDIAN->setSizePolicy(sizePolicy1);
        label_DC_LIFETIME_MEDIAN->setMouseTracking(false);
        label_DC_LIFETIME_MEDIAN->setWordWrap(false);
        line_DC_LIFETIME_MEDIAN = new QLineEdit(tab_DC);
        line_DC_LIFETIME_MEDIAN->setObjectName(QString::fromUtf8("line_DC_LIFETIME_MEDIAN"));
        line_DC_LIFETIME_MEDIAN->setGeometry(QRect(940, 390, 51, 20));
        line_DC_LIFETIME_MEDIAN->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        slider_DC_LIFETIME_SHAPE = new QSlider(tab_DC);
        slider_DC_LIFETIME_SHAPE->setObjectName(QString::fromUtf8("slider_DC_LIFETIME_SHAPE"));
        slider_DC_LIFETIME_SHAPE->setGeometry(QRect(850, 430, 61, 16));
        slider_DC_LIFETIME_SHAPE->setOrientation(Qt::Horizontal);
        slider_DC_LIFETIME_MEDIAN = new QSlider(tab_DC);
        slider_DC_LIFETIME_MEDIAN->setObjectName(QString::fromUtf8("slider_DC_LIFETIME_MEDIAN"));
        slider_DC_LIFETIME_MEDIAN->setGeometry(QRect(850, 390, 61, 16));
        slider_DC_LIFETIME_MEDIAN->setOrientation(Qt::Horizontal);
        alabel_dist_2 = new QLabel(tab_DC);
        alabel_dist_2->setObjectName(QString::fromUtf8("alabel_dist_2"));
        alabel_dist_2->setGeometry(QRect(580, 20, 291, 20));
        alabel_dist_2->setFont(font);
        tabs->addTab(tab_DC, QString());
        tab_TCR = new QWidget();
        tab_TCR->setObjectName(QString::fromUtf8("tab_TCR"));
        layoutWidget_2 = new QWidget(tab_TCR);
        layoutWidget_2->setObjectName(QString::fromUtf8("layoutWidget_2"));
        layoutWidget_2->setGeometry(QRect(20, 40, 421, 191));
        gridLayout_3 = new QGridLayout(layoutWidget_2);
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        gridLayout_3->setContentsMargins(0, 0, 0, 0);
        label_IL2_THRESHOLD = new QMyLabel(layoutWidget_2);
        label_IL2_THRESHOLD->setObjectName(QString::fromUtf8("label_IL2_THRESHOLD"));

        gridLayout_3->addWidget(label_IL2_THRESHOLD, 0, 0, 1, 1);

        label_ACTIVATION_THRESHOLD = new QMyLabel(layoutWidget_2);
        label_ACTIVATION_THRESHOLD->setObjectName(QString::fromUtf8("label_ACTIVATION_THRESHOLD"));

        gridLayout_3->addWidget(label_ACTIVATION_THRESHOLD, 1, 0, 1, 1);

        label_FIRST_DIVISION_THRESHOLD = new QMyLabel(layoutWidget_2);
        label_FIRST_DIVISION_THRESHOLD->setObjectName(QString::fromUtf8("label_FIRST_DIVISION_THRESHOLD"));

        gridLayout_3->addWidget(label_FIRST_DIVISION_THRESHOLD, 2, 0, 1, 1);

        label_DIVISION_THRESHOLD = new QMyLabel(layoutWidget_2);
        label_DIVISION_THRESHOLD->setObjectName(QString::fromUtf8("label_DIVISION_THRESHOLD"));
        QSizePolicy sizePolicy2(QSizePolicy::MinimumExpanding, QSizePolicy::Preferred);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(label_DIVISION_THRESHOLD->sizePolicy().hasHeightForWidth());
        label_DIVISION_THRESHOLD->setSizePolicy(sizePolicy2);
        label_DIVISION_THRESHOLD->setMouseTracking(false);
        label_DIVISION_THRESHOLD->setWordWrap(false);

        gridLayout_3->addWidget(label_DIVISION_THRESHOLD, 3, 0, 1, 1);

        label_EXIT_THRESHOLD = new QMyLabel(layoutWidget_2);
        label_EXIT_THRESHOLD->setObjectName(QString::fromUtf8("label_EXIT_THRESHOLD"));
        label_EXIT_THRESHOLD->setMouseTracking(false);
        label_EXIT_THRESHOLD->setWordWrap(false);

        gridLayout_3->addWidget(label_EXIT_THRESHOLD, 4, 0, 1, 1);

        line_IL2_THRESHOLD = new QLineEdit(layoutWidget_2);
        line_IL2_THRESHOLD->setObjectName(QString::fromUtf8("line_IL2_THRESHOLD"));
        QSizePolicy sizePolicy3(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy3.setHorizontalStretch(120);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(line_IL2_THRESHOLD->sizePolicy().hasHeightForWidth());
        line_IL2_THRESHOLD->setSizePolicy(sizePolicy3);
        line_IL2_THRESHOLD->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(line_IL2_THRESHOLD, 0, 1, 1, 1);

        line_ACTIVATION_THRESHOLD = new QLineEdit(layoutWidget_2);
        line_ACTIVATION_THRESHOLD->setObjectName(QString::fromUtf8("line_ACTIVATION_THRESHOLD"));
        sizePolicy3.setHeightForWidth(line_ACTIVATION_THRESHOLD->sizePolicy().hasHeightForWidth());
        line_ACTIVATION_THRESHOLD->setSizePolicy(sizePolicy3);
        line_ACTIVATION_THRESHOLD->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(line_ACTIVATION_THRESHOLD, 1, 1, 1, 1);

        line_FIRST_DIVISION_THRESHOLD = new QLineEdit(layoutWidget_2);
        line_FIRST_DIVISION_THRESHOLD->setObjectName(QString::fromUtf8("line_FIRST_DIVISION_THRESHOLD"));
        sizePolicy3.setHeightForWidth(line_FIRST_DIVISION_THRESHOLD->sizePolicy().hasHeightForWidth());
        line_FIRST_DIVISION_THRESHOLD->setSizePolicy(sizePolicy3);
        line_FIRST_DIVISION_THRESHOLD->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(line_FIRST_DIVISION_THRESHOLD, 2, 1, 1, 1);

        line_DIVISION_THRESHOLD = new QLineEdit(layoutWidget_2);
        line_DIVISION_THRESHOLD->setObjectName(QString::fromUtf8("line_DIVISION_THRESHOLD"));
        sizePolicy3.setHeightForWidth(line_DIVISION_THRESHOLD->sizePolicy().hasHeightForWidth());
        line_DIVISION_THRESHOLD->setSizePolicy(sizePolicy3);
        line_DIVISION_THRESHOLD->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(line_DIVISION_THRESHOLD, 3, 1, 1, 1);

        line_EXIT_THRESHOLD = new QLineEdit(layoutWidget_2);
        line_EXIT_THRESHOLD->setObjectName(QString::fromUtf8("line_EXIT_THRESHOLD"));
        sizePolicy3.setHeightForWidth(line_EXIT_THRESHOLD->sizePolicy().hasHeightForWidth());
        line_EXIT_THRESHOLD->setSizePolicy(sizePolicy3);
        line_EXIT_THRESHOLD->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(line_EXIT_THRESHOLD, 4, 1, 1, 1);

        label_STIMULATION_LIMIT = new QMyLabel(layoutWidget_2);
        label_STIMULATION_LIMIT->setObjectName(QString::fromUtf8("label_STIMULATION_LIMIT"));
        label_STIMULATION_LIMIT->setMouseTracking(false);
        label_STIMULATION_LIMIT->setWordWrap(false);

        gridLayout_3->addWidget(label_STIMULATION_LIMIT, 5, 0, 1, 1);

        line_STIMULATION_LIMIT = new QLineEdit(layoutWidget_2);
        line_STIMULATION_LIMIT->setObjectName(QString::fromUtf8("line_STIMULATION_LIMIT"));
        sizePolicy3.setHeightForWidth(line_STIMULATION_LIMIT->sizePolicy().hasHeightForWidth());
        line_STIMULATION_LIMIT->setSizePolicy(sizePolicy3);
        line_STIMULATION_LIMIT->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(line_STIMULATION_LIMIT, 5, 1, 1, 1);

        label_40 = new QLabel(tab_TCR);
        label_40->setObjectName(QString::fromUtf8("label_40"));
        label_40->setGeometry(QRect(10, 10, 211, 16));
        QFont font2;
        font2.setBold(true);
        font2.setWeight(75);
        label_40->setFont(font2);
        tabs->addTab(tab_TCR, QString());
        tab_run = new QWidget();
        tab_run->setObjectName(QString::fromUtf8("tab_run"));
        layoutWidget2 = new QWidget(tab_run);
        layoutWidget2->setObjectName(QString::fromUtf8("layoutWidget2"));
        layoutWidget2->setGeometry(QRect(10, 0, 431, 668));
        gridLayout_4 = new QGridLayout(layoutWidget2);
        gridLayout_4->setObjectName(QString::fromUtf8("gridLayout_4"));
        gridLayout_4->setContentsMargins(0, 0, 0, 0);
        label_NX = new QMyLabel(layoutWidget2);
        label_NX->setObjectName(QString::fromUtf8("label_NX"));
        label_NX->setWordWrap(false);

        gridLayout_4->addWidget(label_NX, 0, 0, 1, 1);

        horizontalSpacer_12 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_12, 0, 1, 1, 1);

        spin_NX = new QSpinBox(layoutWidget2);
        spin_NX->setObjectName(QString::fromUtf8("spin_NX"));
        spin_NX->setMaximumSize(QSize(120, 16777215));
        spin_NX->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        spin_NX->setMaximum(999999999);

        gridLayout_4->addWidget(spin_NX, 0, 2, 1, 1);

        label_BLOB_RADIUS = new QMyLabel(layoutWidget2);
        label_BLOB_RADIUS->setObjectName(QString::fromUtf8("label_BLOB_RADIUS"));
        label_BLOB_RADIUS->setWordWrap(false);

        gridLayout_4->addWidget(label_BLOB_RADIUS, 1, 0, 1, 1);

        horizontalSpacer_13 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_13, 1, 1, 1, 1);

        line_BLOB_RADIUS = new QLineEdit(layoutWidget2);
        line_BLOB_RADIUS->setObjectName(QString::fromUtf8("line_BLOB_RADIUS"));
        line_BLOB_RADIUS->setMaximumSize(QSize(120, 16777215));
        line_BLOB_RADIUS->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_BLOB_RADIUS, 1, 2, 1, 1);

        label_TC_FRACTION = new QMyLabel(layoutWidget2);
        label_TC_FRACTION->setObjectName(QString::fromUtf8("label_TC_FRACTION"));
        label_TC_FRACTION->setWordWrap(false);

        gridLayout_4->addWidget(label_TC_FRACTION, 2, 0, 1, 1);

        horizontalSpacer_14 = new QSpacerItem(38, 17, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_14, 2, 1, 1, 1);

        line_TC_FRACTION = new QLineEdit(layoutWidget2);
        line_TC_FRACTION->setObjectName(QString::fromUtf8("line_TC_FRACTION"));
        line_TC_FRACTION->setMaximumSize(QSize(120, 16777215));
        line_TC_FRACTION->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_TC_FRACTION, 2, 2, 1, 1);

        horizontalSpacer_15 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_15, 3, 1, 1, 1);

        line_FLUID_FRACTION = new QLineEdit(layoutWidget2);
        line_FLUID_FRACTION->setObjectName(QString::fromUtf8("line_FLUID_FRACTION"));
        line_FLUID_FRACTION->setMaximumSize(QSize(120, 16777215));
        line_FLUID_FRACTION->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_FLUID_FRACTION, 3, 2, 1, 1);

        label_DC_RADIUS = new QMyLabel(layoutWidget2);
        label_DC_RADIUS->setObjectName(QString::fromUtf8("label_DC_RADIUS"));
        label_DC_RADIUS->setMouseTracking(false);
        label_DC_RADIUS->setWordWrap(false);

        gridLayout_4->addWidget(label_DC_RADIUS, 4, 0, 1, 1);

        horizontalSpacer_16 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_16, 4, 1, 1, 1);

        line_DC_RADIUS = new QLineEdit(layoutWidget2);
        line_DC_RADIUS->setObjectName(QString::fromUtf8("line_DC_RADIUS"));
        line_DC_RADIUS->setMaximumSize(QSize(120, 16777215));
        line_DC_RADIUS->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_DC_RADIUS, 4, 2, 1, 1);

        units_DC_RADIUS = new QLabel(layoutWidget2);
        units_DC_RADIUS->setObjectName(QString::fromUtf8("units_DC_RADIUS"));

        gridLayout_4->addWidget(units_DC_RADIUS, 4, 3, 1, 1);

        label_TC_TO_DC = new QMyLabel(layoutWidget2);
        label_TC_TO_DC->setObjectName(QString::fromUtf8("label_TC_TO_DC"));
        label_TC_TO_DC->setMouseTracking(false);
        label_TC_TO_DC->setWordWrap(false);

        gridLayout_4->addWidget(label_TC_TO_DC, 5, 0, 1, 1);

        horizontalSpacer_17 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_17, 5, 1, 1, 1);

        label_DCrate_100k = new QMyLabel(layoutWidget2);
        label_DCrate_100k->setObjectName(QString::fromUtf8("label_DCrate_100k"));
        label_DCrate_100k->setMouseTracking(false);
        label_DCrate_100k->setWordWrap(false);

        gridLayout_4->addWidget(label_DCrate_100k, 6, 0, 1, 1);

        horizontalSpacer_18 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_18, 6, 1, 1, 1);

        line_DCrate_100k = new QLineEdit(layoutWidget2);
        line_DCrate_100k->setObjectName(QString::fromUtf8("line_DCrate_100k"));
        line_DCrate_100k->setMaximumSize(QSize(120, 16777215));
        line_DCrate_100k->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_DCrate_100k, 6, 2, 1, 1);

        label_T_DC1 = new QMyLabel(layoutWidget2);
        label_T_DC1->setObjectName(QString::fromUtf8("label_T_DC1"));
        label_T_DC1->setMouseTracking(false);
        label_T_DC1->setWordWrap(false);

        gridLayout_4->addWidget(label_T_DC1, 7, 0, 1, 1);

        horizontalSpacer_19 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_19, 7, 1, 1, 1);

        line_T_DC1 = new QLineEdit(layoutWidget2);
        line_T_DC1->setObjectName(QString::fromUtf8("line_T_DC1"));
        line_T_DC1->setMaximumSize(QSize(120, 16777215));
        line_T_DC1->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_T_DC1, 7, 2, 1, 1);

        units_T_DC1 = new QLabel(layoutWidget2);
        units_T_DC1->setObjectName(QString::fromUtf8("units_T_DC1"));

        gridLayout_4->addWidget(units_T_DC1, 7, 3, 1, 1);

        label_T_DC2 = new QMyLabel(layoutWidget2);
        label_T_DC2->setObjectName(QString::fromUtf8("label_T_DC2"));
        label_T_DC2->setMouseTracking(false);
        label_T_DC2->setWordWrap(false);

        gridLayout_4->addWidget(label_T_DC2, 8, 0, 1, 1);

        horizontalSpacer_20 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_20, 8, 1, 1, 1);

        line_T_DC2 = new QLineEdit(layoutWidget2);
        line_T_DC2->setObjectName(QString::fromUtf8("line_T_DC2"));
        line_T_DC2->setMaximumSize(QSize(120, 16777215));
        line_T_DC2->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_T_DC2, 8, 2, 1, 1);

        units_T_DC2 = new QLabel(layoutWidget2);
        units_T_DC2->setObjectName(QString::fromUtf8("units_T_DC2"));

        gridLayout_4->addWidget(units_T_DC2, 8, 3, 1, 1);

        label_RESIDENCE_TIME = new QMyLabel(layoutWidget2);
        label_RESIDENCE_TIME->setObjectName(QString::fromUtf8("label_RESIDENCE_TIME"));
        label_RESIDENCE_TIME->setMouseTracking(false);
        label_RESIDENCE_TIME->setWordWrap(false);

        gridLayout_4->addWidget(label_RESIDENCE_TIME, 10, 0, 1, 1);

        horizontalSpacer_21 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_21, 10, 1, 1, 1);

        line_RESIDENCE_TIME = new QLineEdit(layoutWidget2);
        line_RESIDENCE_TIME->setObjectName(QString::fromUtf8("line_RESIDENCE_TIME"));
        line_RESIDENCE_TIME->setMaximumSize(QSize(120, 16777215));
        line_RESIDENCE_TIME->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_RESIDENCE_TIME, 10, 2, 1, 1);

        units_RESIDENCE_TIME = new QLabel(layoutWidget2);
        units_RESIDENCE_TIME->setObjectName(QString::fromUtf8("units_RESIDENCE_TIME"));

        gridLayout_4->addWidget(units_RESIDENCE_TIME, 10, 3, 1, 1);

        label_INFLAMM_DAYS1 = new QMyLabel(layoutWidget2);
        label_INFLAMM_DAYS1->setObjectName(QString::fromUtf8("label_INFLAMM_DAYS1"));

        gridLayout_4->addWidget(label_INFLAMM_DAYS1, 11, 0, 1, 1);

        horizontalSpacer_22 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_22, 11, 1, 1, 1);

        line_INFLAMM_DAYS1 = new QLineEdit(layoutWidget2);
        line_INFLAMM_DAYS1->setObjectName(QString::fromUtf8("line_INFLAMM_DAYS1"));
        line_INFLAMM_DAYS1->setMaximumSize(QSize(120, 16777215));
        line_INFLAMM_DAYS1->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_INFLAMM_DAYS1, 11, 2, 1, 1);

        label_INFLAMM_DAYS2 = new QMyLabel(layoutWidget2);
        label_INFLAMM_DAYS2->setObjectName(QString::fromUtf8("label_INFLAMM_DAYS2"));

        gridLayout_4->addWidget(label_INFLAMM_DAYS2, 12, 0, 1, 1);

        horizontalSpacer_23 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_23, 12, 1, 1, 1);

        line_INFLAMM_DAYS2 = new QLineEdit(layoutWidget2);
        line_INFLAMM_DAYS2->setObjectName(QString::fromUtf8("line_INFLAMM_DAYS2"));
        line_INFLAMM_DAYS2->setMaximumSize(QSize(120, 16777215));
        line_INFLAMM_DAYS2->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_INFLAMM_DAYS2, 12, 2, 1, 1);

        label_INFLAMM_LEVEL = new QMyLabel(layoutWidget2);
        label_INFLAMM_LEVEL->setObjectName(QString::fromUtf8("label_INFLAMM_LEVEL"));

        gridLayout_4->addWidget(label_INFLAMM_LEVEL, 13, 0, 1, 1);

        horizontalSpacer_24 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_24, 13, 1, 1, 1);

        line_INFLAMM_LEVEL = new QLineEdit(layoutWidget2);
        line_INFLAMM_LEVEL->setObjectName(QString::fromUtf8("line_INFLAMM_LEVEL"));
        line_INFLAMM_LEVEL->setMaximumSize(QSize(120, 16777215));
        line_INFLAMM_LEVEL->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_INFLAMM_LEVEL, 13, 2, 1, 1);

        label_EXIT_REGION = new QMyLabel(layoutWidget2);
        label_EXIT_REGION->setObjectName(QString::fromUtf8("label_EXIT_REGION"));
        label_EXIT_REGION->setMouseTracking(false);
        label_EXIT_REGION->setWordWrap(false);

        gridLayout_4->addWidget(label_EXIT_REGION, 14, 0, 1, 1);

        horizontalSpacer_26 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_26, 14, 1, 1, 1);

        comb_EXIT_REGION = new QComboBox(layoutWidget2);
        comb_EXIT_REGION->setObjectName(QString::fromUtf8("comb_EXIT_REGION"));
        comb_EXIT_REGION->setEnabled(true);
        comb_EXIT_REGION->setMaximumSize(QSize(120, 16777215));

        gridLayout_4->addWidget(comb_EXIT_REGION, 14, 2, 1, 1);

        label_CHEMO_RADIUS = new QMyLabel(layoutWidget2);
        label_CHEMO_RADIUS->setObjectName(QString::fromUtf8("label_CHEMO_RADIUS"));

        gridLayout_4->addWidget(label_CHEMO_RADIUS, 15, 0, 1, 1);

        horizontalSpacer_27 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_27, 15, 1, 1, 1);

        line_CHEMO_RADIUS = new QLineEdit(layoutWidget2);
        line_CHEMO_RADIUS->setObjectName(QString::fromUtf8("line_CHEMO_RADIUS"));
        line_CHEMO_RADIUS->setMaximumSize(QSize(120, 16777215));
        line_CHEMO_RADIUS->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_CHEMO_RADIUS, 15, 2, 1, 1);

        label_CHEMO_K_EXIT = new QMyLabel(layoutWidget2);
        label_CHEMO_K_EXIT->setObjectName(QString::fromUtf8("label_CHEMO_K_EXIT"));
        label_CHEMO_K_EXIT->setMouseTracking(false);
        label_CHEMO_K_EXIT->setWordWrap(false);

        gridLayout_4->addWidget(label_CHEMO_K_EXIT, 16, 0, 1, 1);

        horizontalSpacer_28 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_28, 16, 1, 1, 1);

        line_CHEMO_K_EXIT = new QLineEdit(layoutWidget2);
        line_CHEMO_K_EXIT->setObjectName(QString::fromUtf8("line_CHEMO_K_EXIT"));
        line_CHEMO_K_EXIT->setMaximumSize(QSize(120, 16777215));
        line_CHEMO_K_EXIT->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_CHEMO_K_EXIT, 16, 2, 1, 1);

        label_NDAYS = new QMyLabel(layoutWidget2);
        label_NDAYS->setObjectName(QString::fromUtf8("label_NDAYS"));
        label_NDAYS->setMouseTracking(false);
        label_NDAYS->setWordWrap(false);

        gridLayout_4->addWidget(label_NDAYS, 18, 0, 1, 1);

        horizontalSpacer_29 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_29, 18, 1, 1, 1);

        label_SEED1 = new QMyLabel(layoutWidget2);
        label_SEED1->setObjectName(QString::fromUtf8("label_SEED1"));
        label_SEED1->setMouseTracking(false);
        label_SEED1->setWordWrap(false);

        gridLayout_4->addWidget(label_SEED1, 19, 0, 1, 1);

        horizontalSpacer_30 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_30, 19, 1, 1, 1);

        spin_SEED1 = new QSpinBox(layoutWidget2);
        spin_SEED1->setObjectName(QString::fromUtf8("spin_SEED1"));
        spin_SEED1->setMaximumSize(QSize(120, 16777215));
        spin_SEED1->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        spin_SEED1->setMaximum(999999999);

        gridLayout_4->addWidget(spin_SEED1, 19, 2, 1, 1);

        label_SEED2 = new QMyLabel(layoutWidget2);
        label_SEED2->setObjectName(QString::fromUtf8("label_SEED2"));
        label_SEED2->setMouseTracking(false);
        label_SEED2->setWordWrap(false);

        gridLayout_4->addWidget(label_SEED2, 20, 0, 1, 1);

        horizontalSpacer_31 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_31, 20, 1, 1, 1);

        spin_SEED2 = new QSpinBox(layoutWidget2);
        spin_SEED2->setObjectName(QString::fromUtf8("spin_SEED2"));
        spin_SEED2->setMaximumSize(QSize(120, 16777215));
        spin_SEED2->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        spin_SEED2->setMaximum(999999999);

        gridLayout_4->addWidget(spin_SEED2, 20, 2, 1, 1);

        line_TC_TO_DC = new QLineEdit(layoutWidget2);
        line_TC_TO_DC->setObjectName(QString::fromUtf8("line_TC_TO_DC"));
        line_TC_TO_DC->setMaximumSize(QSize(120, 16777215));
        line_TC_TO_DC->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_TC_TO_DC, 5, 2, 1, 1);

        line_NDAYS = new QLineEdit(layoutWidget2);
        line_NDAYS->setObjectName(QString::fromUtf8("line_NDAYS"));
        line_NDAYS->setMaximumSize(QSize(120, 16777215));
        line_NDAYS->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_NDAYS, 18, 2, 1, 1);

        label_FLUID_FRACTION = new QMyLabel(layoutWidget2);
        label_FLUID_FRACTION->setObjectName(QString::fromUtf8("label_FLUID_FRACTION"));
        label_FLUID_FRACTION->setMouseTracking(false);
        label_FLUID_FRACTION->setWordWrap(false);

        gridLayout_4->addWidget(label_FLUID_FRACTION, 3, 0, 1, 1);

        label_NT_ANIMATION = new QMyLabel(layoutWidget2);
        label_NT_ANIMATION->setObjectName(QString::fromUtf8("label_NT_ANIMATION"));

        gridLayout_4->addWidget(label_NT_ANIMATION, 22, 0, 1, 1);

        horizontalSpacer_32 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_32, 22, 1, 1, 1);

        spin_NT_ANIMATION = new QSpinBox(layoutWidget2);
        spin_NT_ANIMATION->setObjectName(QString::fromUtf8("spin_NT_ANIMATION"));
        spin_NT_ANIMATION->setMaximumSize(QSize(120, 16777215));
        spin_NT_ANIMATION->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        spin_NT_ANIMATION->setMaximum(999999999);

        gridLayout_4->addWidget(spin_NT_ANIMATION, 22, 2, 1, 1);

        label_CHEMO_K_DC = new QMyLabel(layoutWidget2);
        label_CHEMO_K_DC->setObjectName(QString::fromUtf8("label_CHEMO_K_DC"));

        gridLayout_4->addWidget(label_CHEMO_K_DC, 17, 0, 1, 1);

        line_CHEMO_K_DC = new QLineEdit(layoutWidget2);
        line_CHEMO_K_DC->setObjectName(QString::fromUtf8("line_CHEMO_K_DC"));
        line_CHEMO_K_DC->setMaximumSize(QSize(120, 16777215));
        line_CHEMO_K_DC->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_CHEMO_K_DC, 17, 2, 1, 1);

        horizontalSpacer_33 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_33, 17, 1, 1, 1);

        units_CHEMO_RADIUS = new QLabel(layoutWidget2);
        units_CHEMO_RADIUS->setObjectName(QString::fromUtf8("units_CHEMO_RADIUS"));

        gridLayout_4->addWidget(units_CHEMO_RADIUS, 15, 3, 1, 1);

        units_NDAYS = new QLabel(layoutWidget2);
        units_NDAYS->setObjectName(QString::fromUtf8("units_NDAYS"));

        gridLayout_4->addWidget(units_NDAYS, 18, 3, 1, 1);

        units_INFLAMM1 = new QLabel(layoutWidget2);
        units_INFLAMM1->setObjectName(QString::fromUtf8("units_INFLAMM1"));

        gridLayout_4->addWidget(units_INFLAMM1, 11, 3, 1, 1);

        units_INFLAMM2 = new QLabel(layoutWidget2);
        units_INFLAMM2->setObjectName(QString::fromUtf8("units_INFLAMM2"));

        gridLayout_4->addWidget(units_INFLAMM2, 12, 3, 1, 1);

        units_BLOB_RADIUS = new QLabel(layoutWidget2);
        units_BLOB_RADIUS->setObjectName(QString::fromUtf8("units_BLOB_RADIUS"));

        gridLayout_4->addWidget(units_BLOB_RADIUS, 1, 3, 1, 1);

        label_NCPU = new QMyLabel(layoutWidget2);
        label_NCPU->setObjectName(QString::fromUtf8("label_NCPU"));

        gridLayout_4->addWidget(label_NCPU, 21, 0, 1, 1);

        horizontalSpacer_34 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_34, 21, 1, 1, 1);

        spin_NCPU = new QSpinBox(layoutWidget2);
        spin_NCPU->setObjectName(QString::fromUtf8("spin_NCPU"));
        spin_NCPU->setMaximumSize(QSize(120, 16777215));
        spin_NCPU->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        spin_NCPU->setMaximum(999999999);

        gridLayout_4->addWidget(spin_NCPU, 21, 2, 1, 1);

        label_T_DC_INJECTION = new QMyLabel(layoutWidget2);
        label_T_DC_INJECTION->setObjectName(QString::fromUtf8("label_T_DC_INJECTION"));
        label_T_DC_INJECTION->setMouseTracking(false);
        label_T_DC_INJECTION->setWordWrap(false);

        gridLayout_4->addWidget(label_T_DC_INJECTION, 9, 0, 1, 1);

        horizontalSpacer_35 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_35, 9, 1, 1, 1);

        line_T_DC_INJECTION = new QLineEdit(layoutWidget2);
        line_T_DC_INJECTION->setObjectName(QString::fromUtf8("line_T_DC_INJECTION"));
        line_T_DC_INJECTION->setMaximumSize(QSize(120, 16777215));
        line_T_DC_INJECTION->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_T_DC_INJECTION, 9, 2, 1, 1);

        units_T_DC_INJECTION = new QLabel(layoutWidget2);
        units_T_DC_INJECTION->setObjectName(QString::fromUtf8("units_T_DC_INJECTION"));

        gridLayout_4->addWidget(units_T_DC_INJECTION, 9, 3, 1, 1);

        cbox_savepos = new QCheckBox(layoutWidget2);
        cbox_savepos->setObjectName(QString::fromUtf8("cbox_savepos"));

        gridLayout_4->addWidget(cbox_savepos, 23, 2, 1, 1);

        cbox_IN_VITRO = new QCheckBox(tab_run);
        cbox_IN_VITRO->setObjectName(QString::fromUtf8("cbox_IN_VITRO"));
        cbox_IN_VITRO->setGeometry(QRect(650, 10, 141, 18));
        layoutWidget3 = new QWidget(tab_run);
        layoutWidget3->setObjectName(QString::fromUtf8("layoutWidget3"));
        layoutWidget3->setGeometry(QRect(650, 30, 371, 152));
        gridLayout_6 = new QGridLayout(layoutWidget3);
        gridLayout_6->setObjectName(QString::fromUtf8("gridLayout_6"));
        gridLayout_6->setContentsMargins(0, 0, 0, 0);
        label_IV_WELL_DIAMETER = new QMyLabel(layoutWidget3);
        label_IV_WELL_DIAMETER->setObjectName(QString::fromUtf8("label_IV_WELL_DIAMETER"));

        gridLayout_6->addWidget(label_IV_WELL_DIAMETER, 0, 0, 1, 1);

        horizontalSpacer = new QSpacerItem(78, 33, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_6->addItem(horizontalSpacer, 0, 1, 1, 1);

        line_IV_WELL_DIAMETER = new QLineEdit(layoutWidget3);
        line_IV_WELL_DIAMETER->setObjectName(QString::fromUtf8("line_IV_WELL_DIAMETER"));

        gridLayout_6->addWidget(line_IV_WELL_DIAMETER, 0, 2, 1, 1);

        units_IV_WELL_DIAMETER = new QLabel(layoutWidget3);
        units_IV_WELL_DIAMETER->setObjectName(QString::fromUtf8("units_IV_WELL_DIAMETER"));

        gridLayout_6->addWidget(units_IV_WELL_DIAMETER, 0, 3, 1, 1);

        label_IV_NTCELLS = new QMyLabel(layoutWidget3);
        label_IV_NTCELLS->setObjectName(QString::fromUtf8("label_IV_NTCELLS"));

        gridLayout_6->addWidget(label_IV_NTCELLS, 1, 0, 1, 1);

        horizontalSpacer_2 = new QSpacerItem(78, 33, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_6->addItem(horizontalSpacer_2, 1, 1, 1, 1);

        line_IV_NTCELLS = new QLineEdit(layoutWidget3);
        line_IV_NTCELLS->setObjectName(QString::fromUtf8("line_IV_NTCELLS"));

        gridLayout_6->addWidget(line_IV_NTCELLS, 1, 2, 1, 1);

        label_IV_COGNATE_FRACTION = new QMyLabel(layoutWidget3);
        label_IV_COGNATE_FRACTION->setObjectName(QString::fromUtf8("label_IV_COGNATE_FRACTION"));

        gridLayout_6->addWidget(label_IV_COGNATE_FRACTION, 2, 0, 1, 1);

        horizontalSpacer_6 = new QSpacerItem(78, 33, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_6->addItem(horizontalSpacer_6, 2, 1, 1, 1);

        line_IV_COGNATE_FRACTION = new QLineEdit(layoutWidget3);
        line_IV_COGNATE_FRACTION->setObjectName(QString::fromUtf8("line_IV_COGNATE_FRACTION"));

        gridLayout_6->addWidget(line_IV_COGNATE_FRACTION, 2, 2, 1, 1);

        cbox_IV_SHOW_NONCOGNATE = new QCheckBox(tab_run);
        cbox_IV_SHOW_NONCOGNATE->setObjectName(QString::fromUtf8("cbox_IV_SHOW_NONCOGNATE"));
        cbox_IV_SHOW_NONCOGNATE->setGeometry(QRect(650, 190, 161, 18));
        cbox_DC_INJECTION = new QCheckBox(tab_run);
        cbox_DC_INJECTION->setObjectName(QString::fromUtf8("cbox_DC_INJECTION"));
        cbox_DC_INJECTION->setGeometry(QRect(450, 240, 120, 18));
        cbox_USE_TRAFFIC = new QCheckBox(tab_run);
        cbox_USE_TRAFFIC->setObjectName(QString::fromUtf8("cbox_USE_TRAFFIC"));
        cbox_USE_TRAFFIC->setGeometry(QRect(450, 270, 101, 18));
        rbut_SPECIES_1 = new QRadioButton(tab_run);
        buttonGroup_SPECIES = new QButtonGroup(MainWindow);
        buttonGroup_SPECIES->setObjectName(QString::fromUtf8("buttonGroup_SPECIES"));
        buttonGroup_SPECIES->addButton(rbut_SPECIES_1);
        rbut_SPECIES_1->setObjectName(QString::fromUtf8("rbut_SPECIES_1"));
        rbut_SPECIES_1->setGeometry(QRect(470, 30, 61, 18));
        rbut_SPECIES_0 = new QRadioButton(tab_run);
        buttonGroup_SPECIES->addButton(rbut_SPECIES_0);
        rbut_SPECIES_0->setObjectName(QString::fromUtf8("rbut_SPECIES_0"));
        rbut_SPECIES_0->setGeometry(QRect(470, 10, 71, 18));
        rbut_SPECIES_0->setChecked(true);
        cbox_USE_EXIT_CHEMOTAXIS = new QCheckBox(tab_run);
        cbox_USE_EXIT_CHEMOTAXIS->setObjectName(QString::fromUtf8("cbox_USE_EXIT_CHEMOTAXIS"));
        cbox_USE_EXIT_CHEMOTAXIS->setGeometry(QRect(450, 450, 131, 18));
        cbox_USE_DC_CHEMOTAXIS = new QCheckBox(tab_run);
        cbox_USE_DC_CHEMOTAXIS->setObjectName(QString::fromUtf8("cbox_USE_DC_CHEMOTAXIS"));
        cbox_USE_DC_CHEMOTAXIS->setGeometry(QRect(450, 480, 121, 18));
        cbox_COMPUTED_OUTFLOW = new QCheckBox(tab_run);
        cbox_COMPUTED_OUTFLOW->setObjectName(QString::fromUtf8("cbox_COMPUTED_OUTFLOW"));
        cbox_COMPUTED_OUTFLOW->setGeometry(QRect(450, 400, 141, 18));
        label_INPUT_FILE = new QMyLabel(tab_run);
        label_INPUT_FILE->setObjectName(QString::fromUtf8("label_INPUT_FILE"));
        label_INPUT_FILE->setGeometry(QRect(10, 680, 121, 16));
        text_INPUT_FILE = new QLineEdit(tab_run);
        text_INPUT_FILE->setObjectName(QString::fromUtf8("text_INPUT_FILE"));
        text_INPUT_FILE->setGeometry(QRect(130, 680, 151, 20));
        tabs->addTab(tab_run, QString());

        verticalLayout->addWidget(tabs);

        label_input = new QLabel(page_input);
        label_input->setObjectName(QString::fromUtf8("label_input"));
        label_input->setMaximumSize(QSize(16777215, 30));
        QFont font3;
        font3.setFamily(QString::fromUtf8("Segoe UI"));
        font3.setPointSize(16);
        font3.setBold(true);
        font3.setWeight(75);
        label_input->setFont(font3);
        label_input->setAlignment(Qt::AlignCenter);

        verticalLayout->addWidget(label_input);

        text_more = new QTextEdit(page_input);
        text_more->setObjectName(QString::fromUtf8("text_more"));
        QSizePolicy sizePolicy4(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy4.setHorizontalStretch(0);
        sizePolicy4.setVerticalStretch(0);
        sizePolicy4.setHeightForWidth(text_more->sizePolicy().hasHeightForWidth());
        text_more->setSizePolicy(sizePolicy4);
        text_more->setMinimumSize(QSize(480, 0));
        text_more->setMaximumSize(QSize(16777215, 100));
        text_more->setFrameShape(QFrame::StyledPanel);
        text_more->setFrameShadow(QFrame::Sunken);
        text_more->setReadOnly(true);

        verticalLayout->addWidget(text_more);

        stackedWidget->addWidget(page_input);
        page_output = new QWidget();
        page_output->setObjectName(QString::fromUtf8("page_output"));
        verticalLayout_2 = new QVBoxLayout(page_output);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        mdiArea = new QMdiArea(page_output);
        mdiArea->setObjectName(QString::fromUtf8("mdiArea"));
        QSizePolicy sizePolicy5(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy5.setHorizontalStretch(0);
        sizePolicy5.setVerticalStretch(5);
        sizePolicy5.setHeightForWidth(mdiArea->sizePolicy().hasHeightForWidth());
        mdiArea->setSizePolicy(sizePolicy5);
        mdiArea->setMinimumSize(QSize(451, 541));
        mdiArea->setSizeIncrement(QSize(1, 1));
        mdiArea->setBaseSize(QSize(1, 1));
        mdiArea->setFrameShape(QFrame::NoFrame);
        mdiArea->setFrameShadow(QFrame::Plain);
        mdiArea->setTabPosition(QTabWidget::North);

        verticalLayout_2->addWidget(mdiArea);

        box_outputLog = new QTextBrowser(page_output);
        box_outputLog->setObjectName(QString::fromUtf8("box_outputLog"));
        sizePolicy4.setHeightForWidth(box_outputLog->sizePolicy().hasHeightForWidth());
        box_outputLog->setSizePolicy(sizePolicy4);
        box_outputLog->setMinimumSize(QSize(451, 0));
        box_outputLog->setMaximumSize(QSize(16777215, 200));

        verticalLayout_2->addWidget(box_outputLog);

        stackedWidget->addWidget(page_output);
        page_3D = new QWidget();
        page_3D->setObjectName(QString::fromUtf8("page_3D"));
        mdiArea_VTK = new QMdiArea(page_3D);
        mdiArea_VTK->setObjectName(QString::fromUtf8("mdiArea_VTK"));
        mdiArea_VTK->setGeometry(QRect(50, 30, 900, 900));
        progressBar = new QProgressBar(page_3D);
        progressBar->setObjectName(QString::fromUtf8("progressBar"));
        progressBar->setGeometry(QRect(10, 30, 21, 391));
        progressBar->setValue(24);
        progressBar->setOrientation(Qt::Vertical);
        label_hour = new QLabel(page_3D);
        label_hour->setObjectName(QString::fromUtf8("label_hour"));
        label_hour->setGeometry(QRect(10, 0, 46, 20));
        QFont font4;
        font4.setPointSize(10);
        label_hour->setFont(font4);
        hour_display = new QLabel(page_3D);
        hour_display->setObjectName(QString::fromUtf8("hour_display"));
        hour_display->setGeometry(QRect(50, 0, 46, 20));
        hour_display->setFont(font4);
        stackedWidget->addWidget(page_3D);

        gridLayout_5->addWidget(stackedWidget, 0, 0, 1, 1);

        MainWindow->setCentralWidget(centralwidget);
        menubar = new QMenuBar(MainWindow);
        menubar->setObjectName(QString::fromUtf8("menubar"));
        menubar->setGeometry(QRect(0, 0, 1319, 18));
        menuFile = new QMenu(menubar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        menuEdit = new QMenu(menubar);
        menuEdit->setObjectName(QString::fromUtf8("menuEdit"));
        menuABM = new QMenu(menubar);
        menuABM->setObjectName(QString::fromUtf8("menuABM"));
        menuGraphs = new QMenu(menubar);
        menuGraphs->setObjectName(QString::fromUtf8("menuGraphs"));
        menuPlayer = new QMenu(menubar);
        menuPlayer->setObjectName(QString::fromUtf8("menuPlayer"));
        menuSnapshot = new QMenu(menubar);
        menuSnapshot->setObjectName(QString::fromUtf8("menuSnapshot"));
        MainWindow->setMenuBar(menubar);
        statusbar = new QStatusBar(MainWindow);
        statusbar->setObjectName(QString::fromUtf8("statusbar"));
        MainWindow->setStatusBar(statusbar);
        toolBar = new QToolBar(MainWindow);
        toolBar->setObjectName(QString::fromUtf8("toolBar"));
        toolBar->setEnabled(true);
        MainWindow->addToolBar(Qt::TopToolBarArea, toolBar);
        QWidget::setTabOrder(line_TC_COGNATE_FRACTION, line_TC_STIM_RATE_CONSTANT);
        QWidget::setTabOrder(line_TC_STIM_RATE_CONSTANT, line_TC_STIM_HALFLIFE);
        QWidget::setTabOrder(line_TC_STIM_HALFLIFE, line_MOTILITY_BETA);
        QWidget::setTabOrder(line_MOTILITY_BETA, line_MOTILITY_RHO);
        QWidget::setTabOrder(line_MOTILITY_RHO, line_DC_BIND_DELAY);
        QWidget::setTabOrder(line_DC_BIND_DELAY, line_DC_DENS_HALFLIFE);
        QWidget::setTabOrder(line_DC_DENS_HALFLIFE, spin_MAX_TC_BIND);
        QWidget::setTabOrder(spin_MAX_TC_BIND, spin_MAX_COG_BIND);
        QWidget::setTabOrder(spin_MAX_COG_BIND, spin_NX);
        QWidget::setTabOrder(spin_NX, line_BLOB_RADIUS);
        QWidget::setTabOrder(line_BLOB_RADIUS, line_TC_FRACTION);
        QWidget::setTabOrder(line_TC_FRACTION, line_FLUID_FRACTION);
        QWidget::setTabOrder(line_FLUID_FRACTION, line_DC_RADIUS);
        QWidget::setTabOrder(line_DC_RADIUS, line_TC_TO_DC);
        QWidget::setTabOrder(line_TC_TO_DC, line_DCrate_100k);
        QWidget::setTabOrder(line_DCrate_100k, line_T_DC1);
        QWidget::setTabOrder(line_T_DC1, line_T_DC2);
        QWidget::setTabOrder(line_T_DC2, line_RESIDENCE_TIME);
        QWidget::setTabOrder(line_RESIDENCE_TIME, line_INFLAMM_DAYS1);
        QWidget::setTabOrder(line_INFLAMM_DAYS1, line_INFLAMM_DAYS2);
        QWidget::setTabOrder(line_INFLAMM_DAYS2, line_INFLAMM_LEVEL);
        QWidget::setTabOrder(line_INFLAMM_LEVEL, comb_EXIT_REGION);
        QWidget::setTabOrder(comb_EXIT_REGION, line_CHEMO_RADIUS);
        QWidget::setTabOrder(line_CHEMO_RADIUS, line_CHEMO_K_EXIT);
        QWidget::setTabOrder(line_CHEMO_K_EXIT, line_NDAYS);
        QWidget::setTabOrder(line_NDAYS, spin_SEED1);
        QWidget::setTabOrder(spin_SEED1, spin_SEED2);

        menubar->addAction(menuFile->menuAction());
        menubar->addAction(menuEdit->menuAction());
        menubar->addAction(menuABM->menuAction());
        menubar->addAction(menuGraphs->menuAction());
        menubar->addAction(menuPlayer->menuAction());
        menubar->addAction(menuSnapshot->menuAction());
        menuFile->addAction(action_open_input);
        menuFile->addAction(action_load_results);
        menuFile->addSeparator();
        menuFile->addAction(action_save);
        menuFile->addAction(action_saveAs);
        menuABM->addAction(action_run);
        menuABM->addAction(action_pause);
        menuABM->addAction(action_stop);
        menuGraphs->addAction(action_add_graph);
        menuGraphs->addAction(action_remove_graph);
        menuGraphs->addAction(action_remove_all);
        menuPlayer->addAction(action_play_VTK);
        menuPlayer->addAction(action_set_speed);
        menuSnapshot->addAction(action_save_snapshot);
        toolBar->addAction(action_run);
        toolBar->addAction(action_pause);
        toolBar->addAction(action_stop);
        toolBar->addSeparator();
        toolBar->addAction(action_inputs);
        toolBar->addAction(action_outputs);
        toolBar->addAction(action_VTK);

        retranslateUi(MainWindow);

        stackedWidget->setCurrentIndex(0);
        tabs->setCurrentIndex(3);


        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "ABM", 0, QApplication::UnicodeUTF8));
        action_saveAs->setText(QApplication::translate("MainWindow", "Save &As", 0, QApplication::UnicodeUTF8));
        action_save->setText(QApplication::translate("MainWindow", "&Save", 0, QApplication::UnicodeUTF8));
        action_open_input->setText(QApplication::translate("MainWindow", "&Open input", 0, QApplication::UnicodeUTF8));
        action_open_input->setIconText(QApplication::translate("MainWindow", "Open input", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        action_open_input->setToolTip(QApplication::translate("MainWindow", "Open input file", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        action_stop->setText(QApplication::translate("MainWindow", "Stop", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        action_stop->setToolTip(QApplication::translate("MainWindow", "Stop ABM", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        action_stop->setShortcut(QApplication::translate("MainWindow", "F6", 0, QApplication::UnicodeUTF8));
        action_run->setText(QApplication::translate("MainWindow", "Run", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        action_run->setToolTip(QApplication::translate("MainWindow", "Run ABM", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        action_inputs->setText(QApplication::translate("MainWindow", "Inputs", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        action_inputs->setToolTip(QApplication::translate("MainWindow", "Switch to the Inputs", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        action_outputs->setText(QApplication::translate("MainWindow", "Outputs", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        action_outputs->setToolTip(QApplication::translate("MainWindow", "Switch to the Outputs", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        action_VTK->setText(QApplication::translate("MainWindow", "Animation", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        action_VTK->setToolTip(QApplication::translate("MainWindow", "Switch to animation", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        action_pause->setText(QApplication::translate("MainWindow", "Pause", 0, QApplication::UnicodeUTF8));
        action_load_results->setText(QApplication::translate("MainWindow", "&Load results", 0, QApplication::UnicodeUTF8));
        action_load_results->setIconText(QApplication::translate("MainWindow", "Load results", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        action_load_results->setToolTip(QApplication::translate("MainWindow", "Load result file", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        action_add_graph->setText(QApplication::translate("MainWindow", "Add graph", 0, QApplication::UnicodeUTF8));
        action_remove_graph->setText(QApplication::translate("MainWindow", "Remove graph", 0, QApplication::UnicodeUTF8));
        action_remove_all->setText(QApplication::translate("MainWindow", "Remove all", 0, QApplication::UnicodeUTF8));
        action_play_VTK->setText(QApplication::translate("MainWindow", "Play cell animation", 0, QApplication::UnicodeUTF8));
        action_set_speed->setText(QApplication::translate("MainWindow", "Set speed", 0, QApplication::UnicodeUTF8));
        action_save_snapshot->setText(QApplication::translate("MainWindow", "Save snapshot", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_TC_AVIDITY_MEDIAN->setToolTip(QApplication::translate("MainWindow", "TCR avidity has a lognormal distribution, described by the mean and shape parameters.", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_STATUSTIP
        label_TC_AVIDITY_MEDIAN->setStatusTip(QString());
#endif // QT_NO_STATUSTIP
#ifndef QT_NO_WHATSTHIS
        label_TC_AVIDITY_MEDIAN->setWhatsThis(QString());
#endif // QT_NO_WHATSTHIS
        label_TC_AVIDITY_MEDIAN->setText(QApplication::translate("MainWindow", "label_TC_AVIDITY_MEDIAN", 0, QApplication::UnicodeUTF8));
        label_TC_AVIDITY_SHAPE->setText(QApplication::translate("MainWindow", "label_TC_AVIDITY_SHAPE", 0, QApplication::UnicodeUTF8));
        label_TC_COGNATE_FRACTION->setText(QApplication::translate("MainWindow", "label_TC_COGNATE_FRACTION", 0, QApplication::UnicodeUTF8));
        label_TC_STIM_RATE_CONSTANT->setText(QApplication::translate("MainWindow", "label_TC_STIM_RATE_CONSTANT", 0, QApplication::UnicodeUTF8));
        label_TC_STIM_HALFLIFE->setText(QApplication::translate("MainWindow", "label_TC_STIM_HALFLIFE", 0, QApplication::UnicodeUTF8));
        units_TC_STIM_HALFLIFE->setText(QApplication::translate("MainWindow", "hours", 0, QApplication::UnicodeUTF8));
        label_MOTILITY_BETA->setText(QApplication::translate("MainWindow", "label_MOTILITY_BETA", 0, QApplication::UnicodeUTF8));
        label_MOTILITY_RHO->setText(QApplication::translate("MainWindow", "label_MOTILITY_RHO", 0, QApplication::UnicodeUTF8));
        label_DIVIDE1_MEDIAN->setText(QApplication::translate("MainWindow", "label_DIVIDE1_MEDIAN", 0, QApplication::UnicodeUTF8));
        label_DIVIDE1_SHAPE->setText(QApplication::translate("MainWindow", "label_DIVIDE1_SHAPE", 0, QApplication::UnicodeUTF8));
        alabel_dist->setText(QApplication::translate("MainWindow", "Probability distributions", 0, QApplication::UnicodeUTF8));
        label_DIVIDE2_MEDIAN->setText(QApplication::translate("MainWindow", "label_DIVIDE1_MEDIAN", 0, QApplication::UnicodeUTF8));
        label_DIVIDE2_SHAPE->setText(QApplication::translate("MainWindow", "label_DIVIDE1_SHAPE", 0, QApplication::UnicodeUTF8));
        tabs->setTabText(tabs->indexOf(tab_T), QApplication::translate("MainWindow", "T Cell", 0, QApplication::UnicodeUTF8));
        label_DC_BIND_DELAY->setText(QApplication::translate("MainWindow", "label_DC_BIND_DELAY", 0, QApplication::UnicodeUTF8));
        units_DC_BIND_DELAY->setText(QApplication::translate("MainWindow", "min", 0, QApplication::UnicodeUTF8));
        label_DC_DENS_HALFLIFE->setText(QApplication::translate("MainWindow", "label_DC_DENS_HALFLIFE", 0, QApplication::UnicodeUTF8));
        label_MAX_TC_BIND->setText(QApplication::translate("MainWindow", "label_MAX_TC_BIND", 0, QApplication::UnicodeUTF8));
        label_MAX_COG_BIND->setText(QApplication::translate("MainWindow", "label_MAX_COG_BIND", 0, QApplication::UnicodeUTF8));
        units_DC_DENS_HALFLIFE->setText(QApplication::translate("MainWindow", "hours", 0, QApplication::UnicodeUTF8));
        label_DC_ANTIGEN_SHAPE->setText(QApplication::translate("MainWindow", "label_DC_ANTIGEN_SHAPE", 0, QApplication::UnicodeUTF8));
        label_DC_ANTIGEN_MEDIAN->setText(QApplication::translate("MainWindow", "label_DC_ANTIGEN_MEDIAN", 0, QApplication::UnicodeUTF8));
        label_DC_LIFETIME_SHAPE->setText(QApplication::translate("MainWindow", "label_DC_LIFETIME_SHAPE", 0, QApplication::UnicodeUTF8));
        label_DC_LIFETIME_MEDIAN->setText(QApplication::translate("MainWindow", "label_DC_LIFETIME_MEDIAN", 0, QApplication::UnicodeUTF8));
        alabel_dist_2->setText(QApplication::translate("MainWindow", "Probability distributions", 0, QApplication::UnicodeUTF8));
        tabs->setTabText(tabs->indexOf(tab_DC), QApplication::translate("MainWindow", "DC", 0, QApplication::UnicodeUTF8));
        label_IL2_THRESHOLD->setText(QApplication::translate("MainWindow", "label_IL2_THRESHOLD", 0, QApplication::UnicodeUTF8));
        label_ACTIVATION_THRESHOLD->setText(QApplication::translate("MainWindow", "label_ACTIVATION_THRESHOLD", 0, QApplication::UnicodeUTF8));
        label_FIRST_DIVISION_THRESHOLD->setText(QApplication::translate("MainWindow", "label_FIRST_DIVISION_THRESHOLD", 0, QApplication::UnicodeUTF8));
        label_DIVISION_THRESHOLD->setText(QApplication::translate("MainWindow", "label_DIVISION_THRESHOLD", 0, QApplication::UnicodeUTF8));
        label_EXIT_THRESHOLD->setText(QApplication::translate("MainWindow", "label_EXIT_THRESHOLD", 0, QApplication::UnicodeUTF8));
        label_STIMULATION_LIMIT->setText(QApplication::translate("MainWindow", "label_STIMULATION_LIMIT", 0, QApplication::UnicodeUTF8));
        label_40->setText(QApplication::translate("MainWindow", "TCR activation thresholds", 0, QApplication::UnicodeUTF8));
        tabs->setTabText(tabs->indexOf(tab_TCR), QApplication::translate("MainWindow", "TCR Activation", 0, QApplication::UnicodeUTF8));
        label_NX->setText(QApplication::translate("MainWindow", "label_NX", 0, QApplication::UnicodeUTF8));
        label_BLOB_RADIUS->setText(QApplication::translate("MainWindow", "label_BLOB_RADIUS", 0, QApplication::UnicodeUTF8));
        label_TC_FRACTION->setText(QApplication::translate("MainWindow", "label_TC_FRACTION", 0, QApplication::UnicodeUTF8));
        label_DC_RADIUS->setText(QApplication::translate("MainWindow", "label_DC_RADIUS", 0, QApplication::UnicodeUTF8));
        units_DC_RADIUS->setText(QApplication::translate("MainWindow", "\302\265m", 0, QApplication::UnicodeUTF8));
        label_TC_TO_DC->setText(QApplication::translate("MainWindow", "label_TC_TO_DC", 0, QApplication::UnicodeUTF8));
        label_DCrate_100k->setText(QApplication::translate("MainWindow", "label_DCrate_100k", 0, QApplication::UnicodeUTF8));
        label_T_DC1->setText(QApplication::translate("MainWindow", "label_T_DC1", 0, QApplication::UnicodeUTF8));
        units_T_DC1->setText(QApplication::translate("MainWindow", "days", 0, QApplication::UnicodeUTF8));
        label_T_DC2->setText(QApplication::translate("MainWindow", "label_T_DC2", 0, QApplication::UnicodeUTF8));
        units_T_DC2->setText(QApplication::translate("MainWindow", "days", 0, QApplication::UnicodeUTF8));
        label_RESIDENCE_TIME->setText(QApplication::translate("MainWindow", "label_RESIDENCE_TIME", 0, QApplication::UnicodeUTF8));
        units_RESIDENCE_TIME->setText(QApplication::translate("MainWindow", "hours", 0, QApplication::UnicodeUTF8));
        label_INFLAMM_DAYS1->setText(QApplication::translate("MainWindow", "label_INFLAMM_DAYS1", 0, QApplication::UnicodeUTF8));
        label_INFLAMM_DAYS2->setText(QApplication::translate("MainWindow", "label_INFLAMM_DAYS2", 0, QApplication::UnicodeUTF8));
        label_INFLAMM_LEVEL->setText(QApplication::translate("MainWindow", "label_INFLAMM_LEVEL", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_EXIT_REGION->setToolTip(QApplication::translate("MainWindow", "determines blob region for cell exits", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_EXIT_REGION->setText(QApplication::translate("MainWindow", "label_EXIT_REGION", 0, QApplication::UnicodeUTF8));
        comb_EXIT_REGION->clear();
        comb_EXIT_REGION->insertItems(0, QStringList()
         << QApplication::translate("MainWindow", "everywhere", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "lower half of blob", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "blob portals", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "surface portals", 0, QApplication::UnicodeUTF8)
        );
        label_CHEMO_RADIUS->setText(QApplication::translate("MainWindow", "label_CHEMO_RADIUS", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_CHEMO_K_EXIT->setToolTip(QApplication::translate("MainWindow", "level of exit chemotactic influence", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_CHEMO_K_EXIT->setText(QApplication::translate("MainWindow", "label_CHEMO_K_EXIT", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_NDAYS->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        label_NDAYS->setText(QApplication::translate("MainWindow", "label_NDAYS", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_SEED1->setToolTip(QApplication::translate("MainWindow", "seed vector for RNGs", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_SEED1->setText(QApplication::translate("MainWindow", "label_SEED1", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_SEED2->setToolTip(QApplication::translate("MainWindow", "seed vector for RNGs", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_SEED2->setText(QApplication::translate("MainWindow", "label_SEED2", 0, QApplication::UnicodeUTF8));
        label_FLUID_FRACTION->setText(QApplication::translate("MainWindow", "label_FLUID_FRACTION", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_NT_ANIMATION->setToolTip(QApplication::translate("MainWindow", "animation update interval", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_NT_ANIMATION->setText(QApplication::translate("MainWindow", "label_NT_ANIMATION", 0, QApplication::UnicodeUTF8));
        label_CHEMO_K_DC->setText(QApplication::translate("MainWindow", "label_CHEMO_K_DC", 0, QApplication::UnicodeUTF8));
        units_CHEMO_RADIUS->setText(QApplication::translate("MainWindow", "\302\265m", 0, QApplication::UnicodeUTF8));
        units_NDAYS->setText(QApplication::translate("MainWindow", "days", 0, QApplication::UnicodeUTF8));
        units_INFLAMM1->setText(QApplication::translate("MainWindow", "days", 0, QApplication::UnicodeUTF8));
        units_INFLAMM2->setText(QApplication::translate("MainWindow", "days", 0, QApplication::UnicodeUTF8));
        units_BLOB_RADIUS->setText(QApplication::translate("MainWindow", "sites", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_NCPU->setToolTip(QApplication::translate("MainWindow", "animation update interval", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_NCPU->setText(QApplication::translate("MainWindow", "label_NCPU", 0, QApplication::UnicodeUTF8));
        label_T_DC_INJECTION->setText(QApplication::translate("MainWindow", "label_T_DC_INJECTION", 0, QApplication::UnicodeUTF8));
        units_T_DC_INJECTION->setText(QApplication::translate("MainWindow", "hours", 0, QApplication::UnicodeUTF8));
        cbox_savepos->setText(QApplication::translate("MainWindow", "Save cell paths", 0, QApplication::UnicodeUTF8));
        cbox_IN_VITRO->setText(QApplication::translate("MainWindow", "In vitro simulation", 0, QApplication::UnicodeUTF8));
        label_IV_WELL_DIAMETER->setText(QApplication::translate("MainWindow", "label_IV_WELL_DIAMETER", 0, QApplication::UnicodeUTF8));
        units_IV_WELL_DIAMETER->setText(QApplication::translate("MainWindow", "mm", 0, QApplication::UnicodeUTF8));
        label_IV_NTCELLS->setText(QApplication::translate("MainWindow", "label_IV_NTCELLS", 0, QApplication::UnicodeUTF8));
        label_IV_COGNATE_FRACTION->setText(QApplication::translate("MainWindow", "label_IV_COGNATE_FRACTION", 0, QApplication::UnicodeUTF8));
        cbox_IV_SHOW_NONCOGNATE->setText(QApplication::translate("MainWindow", "Display non-cognate T cells ", 0, QApplication::UnicodeUTF8));
        cbox_DC_INJECTION->setText(QApplication::translate("MainWindow", "DCs injected?", 0, QApplication::UnicodeUTF8));
        cbox_USE_TRAFFIC->setText(QApplication::translate("MainWindow", "T cell trafficking? ", 0, QApplication::UnicodeUTF8));
        rbut_SPECIES_1->setText(QApplication::translate("MainWindow", "Human", 0, QApplication::UnicodeUTF8));
        rbut_SPECIES_0->setText(QApplication::translate("MainWindow", "Mouse", 0, QApplication::UnicodeUTF8));
        cbox_USE_EXIT_CHEMOTAXIS->setText(QApplication::translate("MainWindow", "Use exit chemotaxis?", 0, QApplication::UnicodeUTF8));
        cbox_USE_DC_CHEMOTAXIS->setText(QApplication::translate("MainWindow", "Use DC chemotaxis?", 0, QApplication::UnicodeUTF8));
        cbox_COMPUTED_OUTFLOW->setText(QApplication::translate("MainWindow", "Compute T cell outflow?", 0, QApplication::UnicodeUTF8));
        label_INPUT_FILE->setText(QApplication::translate("MainWindow", "Auxiliary input data file", 0, QApplication::UnicodeUTF8));
        tabs->setTabText(tabs->indexOf(tab_run), QApplication::translate("MainWindow", "Run", 0, QApplication::UnicodeUTF8));
        label_input->setText(QApplication::translate("MainWindow", "Inputs", 0, QApplication::UnicodeUTF8));
        text_more->setDocumentTitle(QString());
        label_hour->setText(QApplication::translate("MainWindow", "Hour:", 0, QApplication::UnicodeUTF8));
        hour_display->setText(QString());
        menuFile->setTitle(QApplication::translate("MainWindow", "File", 0, QApplication::UnicodeUTF8));
        menuEdit->setTitle(QApplication::translate("MainWindow", "Edit", 0, QApplication::UnicodeUTF8));
        menuABM->setTitle(QApplication::translate("MainWindow", "ABM", 0, QApplication::UnicodeUTF8));
        menuGraphs->setTitle(QApplication::translate("MainWindow", "Graphs", 0, QApplication::UnicodeUTF8));
        menuPlayer->setTitle(QApplication::translate("MainWindow", "Player", 0, QApplication::UnicodeUTF8));
        menuSnapshot->setTitle(QApplication::translate("MainWindow", "Snapshot", 0, QApplication::UnicodeUTF8));
        toolBar->setWindowTitle(QApplication::translate("MainWindow", "toolBar", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_ABM_GUI_H
