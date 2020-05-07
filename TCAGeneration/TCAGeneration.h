#pragma once

#include <QtWidgets/QDialog>
#include "ui_TCAGeneration.h"
#include <QFileDialog>
#include "tca_network.h"
#include "tca_undirected_graph.h"
#include "tca_point_facility_catchment_areas.h"
#include "tca_mp_facility_catchment_areas.h"
using namespace tca;
//  
class TCAGeneration : public QDialog
{
	Q_OBJECT

public:
	TCAGeneration(QWidget *parent = Q_NULLPTR);

private:
	Ui::TCAGenerationClass ui;

private:

	//-------facility-----------------
	std::string ffacilities; //facility path
	std::string facilityID;   // 	
	bool isunifiedcutoff;
	bool ispointfacility;
	std::string facilitycutoff;
	std::string facilitylabel;
	QStringList facilityfieldnames;

	//-------network-----------------
	std::string froads;   //road path
	std::string egdeID;
	std::string source;
	std::string target;
	bool isundirectedroad;
	std::string direction;
	QStringList networkfieldnames;

	bool bFfromfacility;    // for directed road network, to specify if 
						   //  the direction "from_facility "(default) or "to_facility" direction

	double radius;         //radius for searching neartest points
	double unifiedcutoff;

	bool bAccedges;     // if output the accessible edges, the name of 
	bool bCatchment;    // if output the cacthment areas

	void GetParameters();  // get the parameters needed for generation	

private slots:

	//-------------------network--------------------------
	void on_pBtnOpenet_clicked();			//open road netwrok	
	// undirected and directed road network
	void on_rBtnNetundirected_clicked();
	void on_rBtnNetdirected_clicked();

	// from_facility or to_facility direction
	//void on_rBtnFromfacility_clicked();
	//void on_rBtnTofacility_clicked();

	//------------------facility------------------------------------- 
	void on_pBtnOpenfacility_clicked();    // open facility

	// individual and unified cut-off distances
	void on_rBtnCutoffidv_clicked();
	void on_rBtnCutoffunified_clicked();

	// point or multiple-points based facility
	void on_rBtnPoint_clicked();
	void on_rBtnMpoint_clicked();

	//-------------------output-----------------------------
	// output check box 
	void on_chkBoxAccedges_clicked();		// generate accessible edges
	void on_chkBoxCatchment_clicked();	    // generate catchment areas
	// final generation
	void on_pBtnGenerate_clicked();
};
