#include "TCAGeneration.h"

TCAGeneration::TCAGeneration(QWidget *parent)
	: QDialog(parent)
{
	ui.setupUi(this);

	facilitylabel = "label";   //defaulted label
	facilitycutoff = "delta"; //defaulted value
	unifiedcutoff = 1000;      //defaulted cutoff value
	direction = "directtion";
	bFfromfacility = true;	   // for directed road network only
	bAccedges = false;
	bCatchment = true;	

	ui.lEditNetpath->setReadOnly(true);
	ui.lEditFacilitypath->setReadOnly(true);

	ui.rBtnFromfacility->setChecked(true);
	ui.rBtnTofacility->setChecked(false);

	ui.rBtnNetundirected->setChecked(true);
	ui.rBtnPoint->setChecked(true);
	ui.rBtnCutoffunified->setChecked(true);
	ui.chkBoxCatchment->setChecked(true);
	ui.chkBoxAccedges->setChecked(false);

	QIntValidator* sradiusIntValidator = new QIntValidator;
	sradiusIntValidator->setRange(0, 5000);
	ui.lEditSearchingradius->setValidator(sradiusIntValidator); // nearest point should be a number within a range
	
	QIntValidator* cutoffIntValidator = new QIntValidator;
	cutoffIntValidator->setRange(0, 50000);
	ui.lEditUnifiedcutoff->setValidator(cutoffIntValidator);

	QIcon icon(":/TCAGeneration/Resources/TCAorange.ico");
	setWindowIcon(icon);

	//ispointfacility = true;
	//isundirectedroad = true;
}


//open facility 
void TCAGeneration::on_pBtnOpenfacility_clicked()
{
	QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"), "C://", "shp file(*.shp)");
	ui.lEditFacilitypath->setText(fileName);

	ffacilities = fileName.toStdString();

	GDALAllRegister();
	GDALDataset *poDS = static_cast<GDALDataset*>(GDALOpenEx(ffacilities.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL));
	if (poDS == NULL)
	{
		printf("Open failed.\n");
		exit(1);
	}
	OGRLayer *ogrlayer = poDS->GetLayer(0);

	// here we should restrict the type of file inside as "point"
	// Initialize network edges
	OGRFeatureDefn *ogrFDefn = ogrlayer->GetLayerDefn();
	OGRFeature *ogrFeature;
	int num_fileds = ogrFDefn->GetFieldCount();

	facilityfieldnames.clear();
	for (int i = 0; i < num_fileds; ++i)
	{
		std::string field = ogrFDefn->GetFieldDefn(i)->GetNameRef();
		QString qfield = QString::fromStdString(field);
		facilityfieldnames.append(qfield);
	}

	ui.cbBoxFacilityid->clear(); //  this is a must,in case users open a wrong file at the first operation
	ui.cbBoxFacilityid->addItems(facilityfieldnames);
	ui.cbBoxFacilityid->setSizeAdjustPolicy(QComboBox::AdjustToContents); // adjust the size depending on the 

	// only under the condition of individual 
	if (ui.rBtnCutoffidv->isChecked())
	{
		ui.cbBoxFacilitycutoff->clear();
		ui.cbBoxFacilitycutoff->addItems(facilityfieldnames);
		ui.cbBoxFacilitycutoff->setSizeAdjustPolicy(QComboBox::AdjustToContents);
		ui.lEditUnifiedcutoff->setReadOnly(false);
	}

	// only under the condition of multiple points
	if (ui.rBtnMpoint->isChecked())
	{
		ui.cbBoxFacilitylabel->clear();
		ui.cbBoxFacilitylabel->addItems(facilityfieldnames);
		ui.cbBoxFacilitylabel->setSizeAdjustPolicy(QComboBox::AdjustToContents);
	}
}

// open road netwrok
void TCAGeneration::on_pBtnOpenet_clicked()
{


	QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"), "C://", "shp file(*.shp)");
	ui.lEditNetpath->setText(fileName);

	froads = fileName.toStdString();

	GDALAllRegister();
	GDALDataset *poDS = static_cast<GDALDataset*>(GDALOpenEx(froads.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL));
	if (poDS == NULL)
	{
		printf("Open failed.\n");
		exit(1);
	}
	OGRLayer *ogrlayer = poDS->GetLayer(0);

	// here we should restrict the type of file inside as "point"
	// Initialize network edges
	OGRFeatureDefn *ogrFDefn = ogrlayer->GetLayerDefn();
	OGRFeature *ogrFeature;
	int num_fileds = ogrFDefn->GetFieldCount();

	facilityfieldnames.clear();
	for (int i = 0; i < num_fileds; ++i)
	{
		std::string field = ogrFDefn->GetFieldDefn(i)->GetNameRef();
		QString qfield = QString::fromStdString(field);
		facilityfieldnames.append(qfield);
	}

	ui.cbBoxNetid->clear(); //  this is a must,in case users open a wrong file at the first operation
	ui.cbBoxNetid->addItems(facilityfieldnames);
	ui.cbBoxNetid->setSizeAdjustPolicy(QComboBox::AdjustToContents); // adjust the size depending on the 

	ui.cbBoxNetsource->clear();
	ui.cbBoxNetsource->addItems(facilityfieldnames);
	ui.cbBoxNetsource->setSizeAdjustPolicy(QComboBox::AdjustToContents);

	ui.cbBoxNettarget->clear();
	ui.cbBoxNettarget->addItems(facilityfieldnames);

	if (ui.rBtnNetdirected->isChecked())
	{
		ui.cbBoxNetdirection->clear();
		ui.cbBoxNetdirection->addItems(facilityfieldnames);
		ui.cbBoxNetdirection->setSizeAdjustPolicy(QComboBox::AdjustToContents);
	}

}



void TCAGeneration::on_rBtnNetundirected_clicked()
{
	if (ui.rBtnNetundirected->isChecked())
	{
		ui.cbBoxNetdirection->clear();

		ui.rBtnFromfacility->setChecked(true);
		ui.rBtnTofacility->setCheckable(false);
	}
}

void TCAGeneration::on_rBtnNetdirected_clicked()
{
	if (ui.rBtnNetdirected->isChecked())
	{
		ui.cbBoxNetdirection->clear();
		ui.cbBoxNetdirection->addItems(facilityfieldnames);
		ui.cbBoxNetdirection->setSizeAdjustPolicy(QComboBox::AdjustToContents);
		ui.rBtnTofacility->setCheckable(true);
	}
}


void TCAGeneration::on_rBtnCutoffidv_clicked()
{
	if (ui.rBtnCutoffidv->isChecked())
	{
		ui.cbBoxFacilitycutoff->clear();
		ui.cbBoxFacilitycutoff->addItems(facilityfieldnames);
		ui.cbBoxFacilitycutoff->setSizeAdjustPolicy(QComboBox::AdjustToContents);
		ui.lEditUnifiedcutoff->setReadOnly(true);
	}

}

void TCAGeneration::on_rBtnCutoffunified_clicked()
{
	if (ui.rBtnCutoffunified->isChecked())
	{
		ui.cbBoxFacilitycutoff->clear();
		ui.lEditUnifiedcutoff->setReadOnly(false);
	}
}

// if the point is seleted 
void TCAGeneration::on_rBtnPoint_clicked()
{
	if (ui.rBtnPoint->isChecked())
	{
		ui.cbBoxFacilitylabel->clear();
	}
}

void TCAGeneration::on_rBtnMpoint_clicked()
{
	if (ui.rBtnMpoint->isChecked())
	{
		ui.cbBoxFacilitylabel->clear();
		ui.cbBoxFacilitylabel->addItems(facilityfieldnames);
		ui.cbBoxFacilitylabel->setSizeAdjustPolicy(QComboBox::AdjustToContents);
	}
}

// generate accessible edges
void TCAGeneration::on_chkBoxAccedges_clicked()
{
	if (ui.chkBoxAccedges->isChecked())
	{
		bAccedges = true;
	}
	else
	{
		bAccedges = false;
	}
}

// generate catchment areas
void TCAGeneration::on_chkBoxCatchment_clicked()
{
	if (ui.chkBoxAccedges->isChecked())
	{
		bCatchment = true;
	}
	else
	{
		bCatchment = false;
	}
}

void TCAGeneration::GetParameters()
{
	// facility
	facilityID = ui.cbBoxFacilityid->currentText().toStdString();
	isunifiedcutoff = ui.rBtnCutoffunified->isChecked();
	ispointfacility = ui.rBtnPoint->isChecked();
	bFfromfacility = ui.rBtnFromfacility->isChecked();

	if (isunifiedcutoff) // read unified cut-off from input
	{
		unifiedcutoff = ui.lEditUnifiedcutoff->text().toDouble();
	}
	else
	{
		facilitycutoff = ui.cbBoxFacilitycutoff->currentText().toStdString();
	}

	if (!ispointfacility) // if it is a non-point facility
	{
		facilitylabel = ui.cbBoxFacilitylabel->currentText().toStdString();
	}

	//network
	egdeID = ui.cbBoxNetid->currentText().toStdString();
	source = ui.cbBoxNetsource->currentText().toStdString();
	target = ui.cbBoxNettarget->currentText().toStdString();
	isundirectedroad = ui.rBtnNetundirected->isChecked();
	if (!isundirectedroad) // if itis a directed road graph
	{
		direction = ui.cbBoxNetdirection->currentText().toStdString();
	}

	//searchradius
	radius = ui.lEditSearchingradius->text().toDouble();

	//output
	bAccedges = ui.chkBoxAccedges->isChecked();
	bCatchment = ui.chkBoxCatchment->isChecked();

	//if (!bAccedges && !bCatchment)  
	//{
	//	QMessageBox::information(this, "Error", "You must select one output file");
	//}
}

// run the program and output 
void TCAGeneration::on_pBtnGenerate_clicked() 
{	
	//froads = "C://Users//dlint//source//repos//TCA_Nets//data_input//acticle_experiment//edges_undir.shp";
	//ffacilities = "C://Users//dlint//source//repos//TCA_Nets//data_input//acticle_experiment//line12_stations.shp";
	//MultipleFacilityCatchmentAreas multifacility_cas;
	//multifacility_cas.CalculateCatchmentAreas(froads, ffacilities);

	//---------get parameters based on user specification-----------------------
	GetParameters();
	if (ispointfacility)
	{
		PointFacilityCatchmentAreas multifacility_cas;
		
		if (isundirectedroad) // mode1 1: < point-based facility + undirected graph >
		{
			//----------for debuging convience------------------
			//froads = "D://Coding_projects//Munich_data//testnetwork1.shp";
			//egdeID = "ID";
			//source = "from";
			//target = "to";
			//ffacilities = "D://Coding_projects//Munich_data//teststation.shp";
			//facilityID = "id";
			//isunifiedcutoff = true;
			//unifiedcutoff = 800;
			//radius = 300;
			multifacility_cas.CalculateCatchmentAreas(froads, egdeID, source, target,
				ffacilities, facilityID, isunifiedcutoff, unifiedcutoff, facilitycutoff,
				radius, bAccedges, bCatchment);
		}
		else  // mode1 2: < point-based  facility + directed graph >          
		{
			multifacility_cas.CalculateCatchmentAreasDirected(bFfromfacility, froads, egdeID, source, target, direction,
				ffacilities, facilityID, isunifiedcutoff, unifiedcutoff, facilitycutoff,
				radius, bAccedges, bCatchment);
		}
		// the output are put at the same folder of facility
		std::string basicfpath = ffacilities.substr(0, ffacilities.length() - 4);
		if (bAccedges)
		{
			std::string accedgespath = basicfpath + "_accEdges.shp";
			multifacility_cas.write_accessible_edges(accedgespath);
		}
		if (bCatchment)
		{
			std::string catchmentpath = basicfpath + "_catchmentAreas.shp";
			multifacility_cas.write_contour_polygons(catchmentpath);
		}
	}
	else
	{
		MPFacilityCatchmentAreas mpfacilities;

		if (isundirectedroad) // mode1 3: <MP-point facility + undirected graph>
		{	
			mpfacilities.CalculateCatchmentAreas(froads, egdeID, source, target,
				ffacilities, facilityID, facilitylabel, isunifiedcutoff, unifiedcutoff, facilitycutoff,
				radius, bAccedges, bCatchment);
		}
		else      // mode1 4: <MP-point facility + directed graph>
		{
			mpfacilities.CalculateCatchmentAreasDirected(bFfromfacility, froads, egdeID, source, target, direction,
				ffacilities, facilityID, facilitylabel, isunifiedcutoff, unifiedcutoff, facilitycutoff,
				radius, bAccedges, bCatchment);
		}
		// the output are put at the same folder of facility
		std::string basicfpath = ffacilities.substr(0, ffacilities.length() - 4);
		if (bAccedges)
		{
			std::string accedgespath = basicfpath + "_accEdges.shp";
			mpfacilities.write_accessible_edges(accedgespath);
		}
		if (bCatchment)
		{
			std::string catchmentpath = basicfpath + "_catchmentAreas.shp";
			mpfacilities.write_contour_polygons(catchmentpath);
		}
	}
}
