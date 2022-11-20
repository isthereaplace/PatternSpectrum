#include "generalizedspectrumgui.h"

GeneralizedSpectrumGUI::GeneralizedSpectrumGUI(QWidget *parent)
	: QMainWindow(parent), spectrum(new GeneralizedSpectrum(0, 0.1, 0, 0))
{
	ui.setupUi(this);
	
	// отступ элементов
	ui.centralWidget->setContentsMargins(6, 6, 6, 6);

	// сопоставление раскладок
	// ui.centralWidget->setLayout(ui.horizontalLayout);

	connect(ui.openImageButton, SIGNAL(clicked()), this, SLOT(openImageButtonClicked()));

	connect(this, SIGNAL(newImageLoaded(QString)), ui.Painter, SLOT(imageChanged(QString)));

	connect(ui.circlesCB, SIGNAL(stateChanged(int)), this, SLOT(checkBoxesChanged(int)));
	connect(ui.bonesCB, SIGNAL(stateChanged(int)), this, SLOT(checkBoxesChanged(int)));
	connect(ui.contoursCB, SIGNAL(stateChanged(int)), this, SLOT(checkBoxesChanged(int)));
	connect(ui.imageCB, SIGNAL(stateChanged(int)), this, SLOT(checkBoxesChanged(int)));
	
	connect(ui.radius, SIGNAL(valueChanged(double)), this, SLOT(radiusChanged(double)));
	connect(ui.scale, SIGNAL(valueChanged(int)), this, SLOT(scaleChanged(int)));
	
	connect(ui.refreshSkeletonButton, SIGNAL(clicked()), this, SLOT(updateSkeleton()));

	connect(ui.saveImageButton, SIGNAL(clicked()), this, SLOT(saveImage()));

    connect(ui.savePolyButton, SIGNAL(clicked()), this, SLOT(savePolygons()));

	ui.CustomPlot->addGraph();
}

GeneralizedSpectrumGUI::~GeneralizedSpectrumGUI()
{
	if (spectrum != NULL)
		delete spectrum;
}

void GeneralizedSpectrumGUI::openImageButtonClicked()
{
	QFileDialog* fileDialog = new QFileDialog(this);
	
	fileDialog->setFileMode(QFileDialog::ExistingFile);
	fileDialog->setNameFilter("Изображения (*.png *.jpg *.bmp *.tif)");

	QStringList fileNames;

	if (fileDialog->exec()) {
		fileNames = fileDialog->selectedFiles();
	}

	if (fileNames.size() != 0) {
		
		filepath = QString(fileNames.at(0).toLocal8Bit().constData());
		QString filename = "Имя файла: ";
		int index = 0;
		
        for (int i = 0; i < filepath.size(); i++) {
			if (filepath[i] == '/') {
				index = i + 1;
			}
		}

        for (int i = index; i < filepath.size(); i++) {
			filename.push_back(filepath[i]);
		}

		ui.filenameLabel->setText(filename.toStdString().c_str());
		ui.radius->setValue(0.0);

        updateSkeleton();
	}

}

void GeneralizedSpectrumGUI::checkBoxesChanged(int)
{
	bool flag = false;
	if (ui.Painter->ready == true) {
		ui.Painter->ready = false;
		flag = true;
	}
	ui.Painter->drawBones = (ui.bonesCB->checkState() == Qt::Checked);
	ui.Painter->drawCircles = (ui.circlesCB->checkState() == Qt::Checked);
	ui.Painter->drawContours = (ui.contoursCB->checkState() == Qt::Checked);
	ui.Painter->drawImage = (ui.imageCB->checkState() == Qt::Checked);
	if (flag) {
		ui.Painter->ready = true;
	}
	ui.Painter->makeVisualization();
}

void GeneralizedSpectrumGUI::radiusChanged(double rad){
	ui.Painter->radius = rad;
	ui.Painter->makeVisualization();
}

void GeneralizedSpectrumGUI::scaleChanged(int sc){
	ui.Painter->scope = sc/100.0;
	ui.Painter->makeVisualization();
}

void GeneralizedSpectrumGUI::updateSkeleton()
{
	if (filepath.size() != 0) {

		spectrum->area = ui.area->value();
		spectrum->step = ui.step->value();
		spectrum->makeLacunas = ui.lacunasButton->isChecked();
		spectrum->SetDegree(ui.xDeg->value(), ui.yDeg->value());
		spectrum->ProcessImage(filepath, ui.invertCheckBox->isChecked());
		DrawPlot();

		// время выполнения
		spectrum->skeleton->RTab.Total =
			spectrum->skeleton->RTab.TimeTrace +
			spectrum->skeleton->RTab.TimeTree +
			spectrum->skeleton->RTab.TimeSkelet +
			spectrum->skeleton->RTab.TimePrun +
			spectrum->skeleton->RTab.TimeSpectrum;
		char buffer[500];
		sprintf(buffer,
			"TimeTrace: %d мс\n"
			"TimeTree: %d мс\n"
			"TimeSkelet: %d мс\n"
			"TimePrun: %d мс\n"
			"TimeSpectrum: %d мс\n"
			"Total wrap: %d мс\n"
			"Total: %d мс\n"
			"Connected comp.: %d \n"
			"Polygons: %d \n",
			spectrum->skeleton->RTab.TimeTrace,
			spectrum->skeleton->RTab.TimeTree,
			spectrum->skeleton->RTab.TimeSkelet,
			spectrum->skeleton->RTab.TimePrun,
			spectrum->skeleton->RTab.TimeSpectrum,
			spectrum->skeleton->RTab.Total - spectrum->skeleton->RTab.TimeSpectrum,
			spectrum->skeleton->RTab.Total,
			spectrum->skeleton->RTab.ConnectComp,
			spectrum->skeleton->RTab.Polygons
			);
		ui.timeLabel->setText(buffer);

		ui.Painter->ready = false;

		ui.Painter->skeleton = spectrum->skeleton;
		
		emit(newImageLoaded(spectrum->imagepath));

		ui.Painter->ready = true;
		ui.Painter->makeLacunas = spectrum->makeLacunas;
		ui.Painter->makeVisualization();
	}
}

void GeneralizedSpectrumGUI::DrawPlot(){

	ui.CustomPlot->graph(0)->setPen(QPen(Qt::red));
	ui.CustomPlot->graph(0)->setData(spectrum->radiuses, spectrum->values);
	ui.CustomPlot->graph(0)->rescaleAxes();

	ui.CustomPlot->xAxis->setLabel("Radius");
	ui.CustomPlot->yAxis->setLabel("Opening Area");
	ui.CustomPlot->replot();

}

void GeneralizedSpectrumGUI::saveImage()
{
	QPixmap qpimage = QPixmap::grabWidget(ui.Painter, 0, 0, spectrum->image.width(), spectrum->image.height());

	QFileDialog* fileDialog = new QFileDialog(this);
	
	fileDialog->setFileMode(QFileDialog::AnyFile);
    fileDialog->setNameFilter("Изображение (*.png)");
    fileDialog->setAcceptMode(QFileDialog::AcceptSave);

	QStringList fileNames;

	if (fileDialog->exec()) {
		fileNames = fileDialog->selectedFiles();
	}

	if (fileNames.size() != 0) {
        QString filen = fileNames.at(0);

        if (!filen.endsWith(".png")) {
            filen += ".png";
        }

		qpimage.save(filen);
	}
}

void GeneralizedSpectrumGUI::savePolygons()
{
	if (spectrum->skeleton != NULL) {
        QFileDialog* fileDialog = new QFileDialog(this);

        fileDialog->setFileMode(QFileDialog::AnyFile);
        fileDialog->setNameFilter("Полигоны (*.txt)");
        fileDialog->setAcceptMode(QFileDialog::AcceptSave);

        QStringList fileNames;

        if (fileDialog->exec()) {
            fileNames = fileDialog->selectedFiles();
        }

        if (fileNames.size() != 0) {
            QString filen = fileNames.at(0);

            if (!filen.endsWith(".txt")) {
                filen += ".txt";
            }

            // print polygons in file
            vector<string> edgesList;

			TContour* S = spectrum->skeleton->Boundary->first();
            while (S != NULL) {
                int cornersCount = S->ListPoints->cardinal();
                TPoint** points = new TPoint*[cornersCount];
                int i = 0;
                TPoint* Corn = S->ListPoints->first();
                while (Corn != NULL)
                {
                    points[i++] = Corn;
                    Corn = Corn->getNext();
                }
                for (int j = 0; j < cornersCount - 1; j++) {
                    stringstream ss;
                    ss << points[j]->X << " " << points[j]->Y << " " << points[j + 1]->X << " " << points[j + 1]->Y;
                    edgesList.push_back(ss.str());
                }
                stringstream ss;
                ss << points[cornersCount - 1]->X << " " << points[cornersCount - 1]->Y << " " << points[0]->X << " " << points[0]->Y;
                edgesList.push_back(ss.str());

                delete points;
                S = S->getNext();
            }

            ofstream f(filen.toStdString().c_str());
            f << "0" << std::endl;
            f << edgesList.size() << std::endl;
            for (ulong i = 0; i < edgesList.size(); ++i) {
                f << edgesList[i] << std::endl;
            }
            f.close();
        }
    }
}
