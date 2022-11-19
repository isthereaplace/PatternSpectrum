#include "mypainter.h"

#define _USE_MATH_DEFINES
#include <math.h>

MyPainter::MyPainter(QWidget *parent)
	: QGraphicsView(parent)
{
	ready = false;

	drawBones = true;
	drawCircles = false;
	drawContours = true;
	drawImage = true;
	radius = 0.0;
	scope = 1.0;
	scene = new QGraphicsScene();
	makeLacunas = false;
	skeleton = NULL;
	setScene(scene);
	setRenderHint(QPainter::Antialiasing);
}

MyPainter::~MyPainter()
{
	if (scene != NULL){
		delete scene;
	}
}

void MyPainter::makeVisualization()
{
	if (ready) {

		scene->clear();

		if (drawImage) {
			scene->addPixmap(QPixmap::fromImage(image));
		}
		if (skeleton == NULL) {
			resetMatrix();
			scale(scope, scope);
			show();
			return;
		}

		if (makeLacunas) {
			auto iNode = skeleton->AllNodes.begin();
			while (iNode != skeleton->AllNodes.end() && (*iNode)->r() < radius) {
				Lacuna* Lacn = (*iNode)->UpdateLacunas();
				if (Lacn != NULL)
					skeleton->Lacunas.push_back(Lacn);
				iNode = skeleton->AllNodes.erase(iNode);
			}
		}
		else {
			auto iNode = skeleton->AllNodes.begin();
			while (iNode != skeleton->AllNodes.end() && (*iNode)->r() < radius) {
				(*iNode)->RebuildNeighborhood();
				iNode = skeleton->AllNodes.erase(iNode);
			}
		}

		QPen pen(Qt::white, image.width() / 1000.0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);

		// Рисуем скелеты
		if (drawBones) {
			for (auto iBone : skeleton->AllBones) {
				QPointF Extr = iBone->GetExtremePoint(radius);
				if (iBone->dest->r() >= radius) {
					double xMin = Extr.x();
					double yMin = Extr.y();
					double xMax = iBone->dest->X();
					double yMax = iBone->dest->Y();
					scene->addLine(QLineF(xMin + 0.5, yMin + 0.5, xMax + 0.5, yMax + 0.5), pen);

					double angle = atan2(yMin - yMax, xMin - xMax);
					double len = image.width() / 100;
					double dx1 = len * cos(angle + M_PI / 6);
					double dy1 = len * sin(angle + M_PI / 6);
					double dx2 = len * cos(angle - M_PI / 6);
					double dy2 = len * sin(angle - M_PI / 6);
					scene->addLine(QLineF(xMax + 0.5, yMax + 0.5, xMax + dx1 + 0.5, yMax + dy1 + 0.5), pen);
					scene->addLine(QLineF(xMax + 0.5, yMax + 0.5, xMax + dx2 + 0.5, yMax + dy2 + 0.5), pen);
				}
			}

			for (auto iBone : skeleton->AllBones) {
				if (iBone->org->r() < radius) {
					QPointF Extr = iBone->GetExtremePoint(radius);
					pen.setColor(makeLacunas ? iBone->Lacuna->color : Qt::red);
					scene->addLine(QLineF(iBone->org->X() + 0.5, iBone->org->Y() + 0.5, Extr.x() + 0.5, Extr.y() + 0.5), pen);
				}
			}
		}

		if (drawCircles) {
			for (auto iBone : skeleton->AllBones) {
				if (iBone->org->r() < radius && iBone->dest->r() >= radius) {
					QPointF A = iBone->GetExtremePoint(radius);
					QPointF A1, A2;
					if (iBone->type == Linear) {
						A1 = ProjectPoint(A, (TEdge*)iBone->Sites[0]);
						A2 = ProjectPoint(A, (TEdge*)iBone->Sites[1]);
					}
					else if (iBone->type == Parabolic) {
						TVertex* V = (TVertex*)iBone->Sites[0];
						A1 = QPointF(V->p->X, V->p->Y);
						A2 = ProjectPoint(A, (TEdge*)iBone->Sites[1]);
					}
					else {
						TVertex* V1 = (TVertex*)iBone->Sites[0];
						TVertex* V2 = (TVertex*)iBone->Sites[1];
						A1 = QPointF(V1->p->X, V1->p->Y);
						A2 = QPointF(V2->p->X, V2->p->Y);
					}
					double alpha = atan2(-A1.y() + A.y(), A1.x() - A.x()) / M_PI * 180;
					double beta = atan2(-A2.y() + A.y(), A2.x() - A.x()) / M_PI * 180;
					if (abs(beta - alpha) > 180) {
						if (alpha < beta) {
							alpha += 360;
						}
						else {
							beta += 360;
						}
					}
					pen.setColor(makeLacunas ? iBone->Lacuna->color : Qt::red);
					QGraphicsEllipseItem* arc = new QGraphicsEllipseItem(A.x() - radius + 0.5, A.y() - radius + 0.5, 2 * radius, 2 * radius);
					arc->setStartAngle(16 * alpha);
					arc->setSpanAngle(16 * (beta - alpha));
					arc->setPen(pen);
					scene->addItem(arc);
				}
			}
		}
		
		pen.setWidth(2 * pen.width());
		if (makeLacunas) {
			// Рисуем триангуляции
			pen.setStyle(Qt::DotLine);
			for (auto iLacn : skeleton->Lacunas) {
				pen.setColor(iLacn->color);
				list<QPair<int, int>> Edges;
				vector<double> Vert;
				iLacn->TriangulateForViz(radius, Edges, Vert);
				for (auto edge : Edges) {
					scene->addLine(QLineF(Vert[2 * edge.first ] + 0.5, Vert[2 * edge.first  + 1] + 0.5,
							              Vert[2 * edge.second] + 0.5, Vert[2 * edge.second + 1] + 0.5), pen);
				}
			}
		}
		else{
			pen.setColor(Qt::green);
			for (auto iBone : skeleton->AllBones) {
				if (!iBone->Adjacent.empty()) {
					QPointF A = iBone->GetExtremePoint(radius);
					for (auto pair : iBone->Adjacent) {
						if (pair.Bone > iBone) {
							QPointF B = pair.Bone->GetExtremePoint(radius);
							scene->addLine(QLineF(A.x() + 0.5, A.y() + 0.5, B.x() + 0.5, B.y() + 0.5), pen);
						}
					}
				}
			}
		}

		resetMatrix();
		scale(scope, scope);
		show();
	}

}

void MyPainter::imageChanged(QString filepath_)
{
	filepath = filepath_;
	image = QImage(filepath);
}
