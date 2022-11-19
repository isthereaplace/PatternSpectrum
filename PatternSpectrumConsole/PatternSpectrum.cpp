#include "PatternSpectrum.h"


PatternSpectrum::PatternSpectrum(int a, double s)
	: skeleton(NULL), srcimg(NULL), area(a), step(s), rMax(-1), makeLacunas(false)
{}


PatternSpectrum::~PatternSpectrum()
{

}


void PatternSpectrum::ProcessImage(QString path, bool invert)
{
	imagepath = path;
	if (skeleton != NULL) {
		delete skeleton;
	}
	if (srcimg) {
		delete srcimg;
	}

	if (imagepath.right(4) == ".txt") {
		skeleton = new TPolFigure(imagepath);
		skeleton->MakeTriangDel();
	}
	else {
		image = QImage(imagepath);
		srcimg = new BitRaster(image.width(), image.height());
		for (int i = 0; i < image.height(); i++) {
			for (int j = 0; j < image.width(); j++) {
				double alpha = qAlpha(image.pixel(j, i)) / 255.0;
				int gray = qGray(image.pixel(j, i))*alpha + 255*(1-alpha);
				bool isBlack = gray < 128;
				srcimg->setBit(j, i, isBlack != invert);
			}
		}
		double pruning = 0;
		BondSkeletTrans(srcimg, pruning, area, skeleton);
	}

	uint TT = TimeInMilSecond();
	skeleton->MonotonicSubdivision();
	CalcSpectrum();
	skeleton->RTab.TimeSpectrum = TimeInMilSecond() - TT;

	int index = imagepath.lastIndexOf('.');
	QByteArray ba = (imagepath.left(index) + ".bin").toLocal8Bit();
	char* reportpath = ba.data();
	fid = fopen(reportpath, "wb");
	int len = values.size();
	fwrite(&len, sizeof(int), 1, fid);
	fwrite(radiuses.data(), sizeof(double), len, fid);
	fwrite(values.data(), sizeof(double), len, fid);
	if (imagepath.right(4) == ".txt") {
		ReportSkeleton();
		ReportContour();
	}
	fclose(fid);
}


bool CompareOrg(TBone* BoneA, TBone* BoneB){
	if (BoneA->dest == BoneB->org)
		return true;
	if (BoneA->org->r() == BoneB->org->r()){
		return BoneA->dest->r() < BoneB->dest->r();
	}
	return BoneA->org->r() < BoneB->org->r();
}


bool CompareDest(TBone* BoneA, TBone* BoneB){
	return BoneA->dest->r() < BoneB->dest->r();
}


bool CompareRadius(TNode* NodeA, TNode* NodeB) {
	return NodeA->r() < NodeB->r();
}

double PatternSpectrum::BoneValue(TBone* Bone, double rad){
	return Bone->Square(rad);
}

double PatternSpectrum::SectorValue(TBone* Bone, double rad){
	return Bone->SectorArea(rad);
}

double PatternSpectrum::LensValue(TBone* BoneA, TBone* BoneB, double rad){
	if (rad > BoneA->dest->r() || rad > BoneB->dest->r())
		return 0;
	QPointF A = BoneA->GetExtremePoint(rad);
	QPointF B = BoneB->GetExtremePoint(rad);
	return PairSquare(A, B, rad);
}

void PatternSpectrum::CalcSpectrum() {

	if (skeleton->Components->cardinal() == 0){
		radiuses.resize(1);
		values.resize(1);
		radiuses[0] = 0;
		values[0] = 0;
		return;
	}

	skeleton->AllNodes.sort(CompareRadius);
	skeleton->AllBones.sort(CompareOrg);

	rMax = (*max_element(skeleton->AllBones.begin(), skeleton->AllBones.end(), CompareDest))->dest->r();
	int nSteps = floor((rMax + 1e-6)/step) + 1;
	radiuses.resize(nSteps + 1);
	values.resize(nSteps + 1);

#ifdef SPECTRUMGUI
	list<TNode*> NodesBackup = skeleton->AllNodes;
	list<TBone*> BonesBackup = skeleton->AllBones;
#endif

	for (int iRad = 0; iRad < nSteps + 1; iRad++) {

		double radius = iRad*step;
		radiuses[iRad] = radius;
		values[iRad] = 0;

		if (makeLacunas) {
			auto iNode = skeleton->AllNodes.begin();
			while (iNode != skeleton->AllNodes.end() && (*iNode)->r() < radius) {
				Lacuna* Lacn = (*iNode)->UpdateLacunas();
				if (Lacn != NULL)
					skeleton->Lacunas.push_back(Lacn);
				iNode = skeleton->AllNodes.erase(iNode);
			}

			auto iLacn = skeleton->Lacunas.begin();
			while (iLacn != skeleton->Lacunas.end()) {
				values[iRad] -= (*iLacn)->Triangulate(radius);
				if ((*iLacn)->Bones.empty()) {
					delete *iLacn;
					iLacn = skeleton->Lacunas.erase(iLacn);
				}
				else iLacn++;
			}
		}
		else {
			auto iNode = skeleton->AllNodes.begin();
			while (iNode != skeleton->AllNodes.end() && (*iNode)->r() < radius) {
				(*iNode)->RebuildNeighborhood();
				iNode = skeleton->AllNodes.erase(iNode);
			}
		}

		auto iBone = skeleton->AllBones.begin();
		while (iBone != skeleton->AllBones.end()) {
			if ((*iBone)->dest->r() < radius) {
				iBone = skeleton->AllBones.erase(iBone);
			}
			else{
				values[iRad] += BoneValue(*iBone, radius);
				if ((*iBone)->org->r() < radius) {
					values[iRad] += SectorValue(*iBone, radius);

					if (!makeLacunas) {
						auto Pair = (*iBone)->Adjacent.begin();
						while (Pair != (*iBone)->Adjacent.end()) {
							double sub = LensValue(*iBone, Pair->Bone, radius);
							if (sub > 0) {
								values[iRad] -= sub / 2;
								Pair++;
							}
							else {
								Pair->Bone->Adjacent.erase(Pair->Link);
								Pair = (*iBone)->Adjacent.erase(Pair);
							}
						}
					}
				}
				iBone++;
			}
		}
	}

#ifdef SPECTRUMGUI
	skeleton->AllNodes = NodesBackup;
	skeleton->AllBones = BonesBackup;
#endif

	if (makeLacunas) {
		for (auto iLacn = skeleton->Lacunas.begin(); iLacn != skeleton->Lacunas.end(); iLacn++) {
			delete (*iLacn);
		}
		skeleton->Lacunas.clear();
		for (auto iBone = skeleton->AllBones.begin(); iBone != skeleton->AllBones.end(); iBone++)
			(*iBone)->Lacuna = NULL;
	}
}


void PatternSpectrum::ReportSkeleton(){

	int index = imagepath.lastIndexOf('.');
	QByteArray ba = (imagepath.left(index) + ".skel").toLocal8Bit();
	char* reportpath = ba.data();

	FILE* fid = fopen(reportpath, "wb");
	int n = skeleton->Components->cardinal();
	fwrite(&n, sizeof(int), 1, fid);
	for (TConnected* Com = skeleton->Components->first(); Com != NULL; Com = Com->getNext()){
		double params[14];
		int k = Com->Bones->cardinal();
		fwrite(&k, sizeof(int), 1, fid);
		for (TBone* Bone = Com->Bones->first(); Bone != NULL; Bone = Bone->getNext()){
			Bone->DetermineType();
			fwrite(&(Bone->type), sizeof(int), 1, fid);
			params[0] = Bone->org->X();
			params[1] = Bone->org->Y();
			params[2] = Bone->org->r();
			params[3] = Bone->dest->X();
			params[4] = Bone->dest->Y();
			params[5] = Bone->dest->r();
			int idxV, idxE;		// По каким индексам расположены сайт-точка и сайт-сегмент 

			switch (Bone->type){
			case Linear:
				params[6] = ((TEdge*)Bone->Sites[0])->org->X;
				params[7] = ((TEdge*)Bone->Sites[0])->org->Y;
				params[8] = ((TEdge*)Bone->Sites[0])->dest->X;
				params[9] = ((TEdge*)Bone->Sites[0])->dest->Y;
				params[10] = ((TEdge*)Bone->Sites[1])->org->X;
				params[11] = ((TEdge*)Bone->Sites[1])->org->Y;
				params[12] = ((TEdge*)Bone->Sites[1])->dest->X;
				params[13] = ((TEdge*)Bone->Sites[1])->dest->Y;
				fwrite(params, sizeof(double), 14, fid);
				break;
			case Parabolic:
				if (Bone->Sites[0]->IsVertex()){
					idxV = 0; idxE = 1;
				}
				else{
					idxV = 1; idxE = 0;
				}
				params[6] = ((TVertex*)Bone->Sites[idxV])->p->X;
				params[7] = ((TVertex*)Bone->Sites[idxV])->p->Y;
				params[8] = ((TEdge*)Bone->Sites[idxE])->org->X;
				params[9] = ((TEdge*)Bone->Sites[idxE])->org->Y;
				params[10] = ((TEdge*)Bone->Sites[idxE])->dest->X;
				params[11] = ((TEdge*)Bone->Sites[idxE])->dest->Y;
				fwrite(params, sizeof(double), 12, fid);
				break;
			case Hyperbolic:
				params[6] = ((TVertex*)Bone->Sites[0])->p->X;
				params[7] = ((TVertex*)Bone->Sites[0])->p->Y;
				params[8] = ((TVertex*)Bone->Sites[1])->p->X;
				params[9] = ((TVertex*)Bone->Sites[1])->p->Y;
				fwrite(params, sizeof(double), 10, fid);
				break;
			}
		}
	}
	fclose(fid);
}


void PatternSpectrum::ReportContour(){
	int index = imagepath.lastIndexOf('.');
	QByteArray ba = (imagepath.left(index) + ".cont").toLocal8Bit();
	char* reportpath = ba.data();

	FILE* fid = fopen(reportpath, "wb");
	int n = skeleton->Components->cardinal();
	fwrite(&n, sizeof(int), 1, fid);
	double coords[2];
	for (TConnected* Com = skeleton->Components->first(); Com != NULL; Com = Com->getNext()){
		int m = Com->HoleList.size() + 1;
		fwrite(&m, sizeof(int), 1, fid);
		int k = Com->Border->ListPoints->cardinal();
		fwrite(&k, sizeof(int), 1, fid);
		for (Point* Pt = Com->Border->ListPoints->first(); Pt != NULL; Pt = Pt->getNext()){
			coords[0] = Pt->X;
			coords[1] = Pt->Y;
			fwrite(coords, sizeof(double), 2, fid);
		}
		for (int i = 0; i < Com->HoleList.size(); i++){
			int k = Com->HoleList[i]->ListPoints->cardinal();
			fwrite(&k, sizeof(int), 1, fid);
			for (Point* Pt = Com->HoleList[i]->ListPoints->first(); Pt != NULL; Pt = Pt->getNext()){
				coords[0] = Pt->X;
				coords[1] = Pt->Y;
				fwrite(coords, sizeof(double), 2, fid);
			}
		}
	}
	fclose(fid);
}