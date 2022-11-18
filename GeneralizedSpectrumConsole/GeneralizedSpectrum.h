#include <QString>
#include <QVector>

#include "../SkeletonLib/BSTrans.h"
#include "../PatternSpectrumConsole/PatternSpectrum.h"

class GeneralizedSpectrum : public PatternSpectrum
{
public:
	GeneralizedSpectrum(int area, double step, int x, int y);
	~GeneralizedSpectrum();
	void SetDegree(int px, int py);

private:

	enum FuncType{
		Sin,
		Cos,
	};

	struct TrigSummand{
		double coef;		// Коэффициент при слагаемом
		FuncType func;		// Тип функции - синус или косинус
		int angle;			// Множитель при угле
		TrigSummand() : coef(1), func(Cos), angle(0) {}
		TrigSummand(double c, FuncType f, int a) : coef(c), func(f), angle(a) {}
	};

	int powx;
	int powy;
	vector<TrigSummand>** Code;

	void CollectSimilar(vector<TrigSummand> &Code);
	vector<TrigSummand> SinCosProd(int psin, int pcos);

	double PolygonMoment(vector<QPointF> Points);
	double SectorMoment(QPointF A, double r, double phi0, double phi1);

	double BoneValue(TBone* Bone, double rad);
	double SectorValue(TBone* Bone, double rad);
	double LensValue(TBone* BoneA, TBone* BoneB, double rad);
};