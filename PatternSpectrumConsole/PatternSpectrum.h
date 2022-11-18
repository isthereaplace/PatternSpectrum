#include <QString>
#include <QVector>

#include "../SkeletonLib/BSTrans.h"

class PatternSpectrum
{
	friend class PatternSpectrumGUI;
	friend class GeneralizedSpectrumGUI;
public:
	PatternSpectrum(int area, double step);
	~PatternSpectrum();
	void ProcessImage(QString path, bool invert = false);

protected:
	virtual double BoneValue(TBone* Bone, double rad);
	virtual double SectorValue(TBone* Bone, double rad);
	virtual double LensValue(TBone* BoneA, TBone* BoneB, double rad);
	virtual void CalcSpectrum();

private:
	QString imagepath;
	FILE* fid;
	QImage image;
	BitRaster* srcimg;
	TPolFigure* skeleton;

	int area;
	double step;
	bool makeLacunas;
	double rMax;
	QVector<double> radiuses;
	QVector<double> values;

	void ReportSkeleton();
	void ReportContour();
};