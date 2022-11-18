#include "GeneralizedSpectrum.h"

#define _USE_MATH_DEFINES
#include <math.h>


GeneralizedSpectrum::GeneralizedSpectrum(int a, double s, int x, int y)
: PatternSpectrum(a, s), powx(x), powy(y)
{
	Code = new vector<TrigSummand>*[powy + 1];
	for (int j = 0; j < powy + 1; j++){
		Code[j] = new vector<TrigSummand>[powx + 1];
		for (int i = 0; i < powx + 1; i++)
			Code[j][i] = SinCosProd(j, i);
	}
}


GeneralizedSpectrum::~GeneralizedSpectrum()
{
	for (int i = 0; i < powy + 1; i++){
		delete[] Code[i];
	}
	delete[] Code;
}


void GeneralizedSpectrum::SetDegree(int px, int py){
	if (px > powx || py >= powy){
		for (int i = 0; i < powy + 1; i++){
			delete[] Code[i];
		}
		delete[] Code;
		Code = new vector<TrigSummand>*[py + 1];
		for (int j = 0; j < py + 1; j++){
			Code[j] = new vector<TrigSummand>[px + 1];
			for (int i = 0; i < px + 1; i++)
				Code[j][i] = SinCosProd(j, i);
		}
	}
	powx = px;
	powy = py;
}


void GeneralizedSpectrum::CollectSimilar(vector<TrigSummand> &Code){
	for (auto it = Code.begin(); it != Code.end(); it++){
		if (it->func == Sin){
			if (it->angle < 0){
				it->angle = -it->angle;
				it->coef = -it->coef;
			}
		}
		else
			it->angle = abs(it->angle);
	}
	for (auto succ = Code.begin(); succ != Code.end(); succ++){
		bool exit = false;
		for (auto pred = Code.begin(); pred != succ && !exit; pred++){
			if (pred->func == succ->func && pred->angle == succ->angle){
				pred->coef += succ->coef;
				auto backup = succ - 1;
				Code.erase(succ);
				succ = backup;
				exit = true;
			}
		}
	}
}


vector<GeneralizedSpectrum::TrigSummand> GeneralizedSpectrum::SinCosProd(int psin, int pcos){

	// Разворачиваем синус
	vector<TrigSummand> CodeSin(1);
	int m = psin - 1;
	if (psin > 0)
		CodeSin[0] = TrigSummand(1, Sin, 1);
	else
		CodeSin[0] = TrigSummand(1, Cos, 0);	// Таким хитрым способом кодируем константу 1

	while (m > 0){
		vector<TrigSummand> CodeTemp(2 * CodeSin.size());
		for (int i = 0; i < CodeSin.size(); i++){
			if (CodeSin[i].func == Sin){
				// sin(a) * sin(b)
				CodeTemp[2 * i] = TrigSummand(-CodeSin[i].coef, Cos, 1 + CodeSin[i].angle);
				CodeTemp[2 * i + 1] = TrigSummand(CodeSin[i].coef, Cos, 1 - CodeSin[i].angle);
			}
			else{
				// sin(a) * cos(b)
				CodeTemp[2 * i] = TrigSummand(CodeSin[i].coef, Sin, 1 + CodeSin[i].angle);
				CodeTemp[2 * i + 1] = TrigSummand(CodeSin[i].coef, Sin, 1 - CodeSin[i].angle);
			}
		}
		CodeSin = CodeTemp;
		m--;
	}
	int denoms = CodeSin.size();
	CollectSimilar(CodeSin);

	// Разворачиваем косинус
	vector<TrigSummand> CodeCos(1);
	int n = pcos - 1;
	if (pcos > 0)
		CodeCos[0] = TrigSummand(1, Cos, 1);
	else
		CodeCos[0] = TrigSummand(1, Cos, 0);	// Таким хитрым способом кодируем константу 1

	while (n > 0){
		vector<TrigSummand> CodeTemp(2 * CodeCos.size());
		for (int i = 0; i < CodeCos.size(); i++){
			if (CodeCos[i].func == Sin){
				// cos(a) * sin(b)
				CodeTemp[2 * i] = TrigSummand(CodeCos[i].coef, Sin, 1 + CodeCos[i].angle);
				CodeTemp[2 * i + 1] = TrigSummand(-CodeCos[i].coef, Sin, 1 - CodeCos[i].angle);
			}
			else{
				// cos(a) * cos(b)
				CodeTemp[2 * i] = TrigSummand(CodeCos[i].coef, Cos, 1 + CodeCos[i].angle);
				CodeTemp[2 * i + 1] = TrigSummand(CodeCos[i].coef, Cos, 1 - CodeCos[i].angle);
			}
		}
		CodeCos = CodeTemp;
		n--;
	}
	int denomc = CodeCos.size();
	CollectSimilar(CodeCos);

	m = CodeSin.size();
	n = CodeCos.size();
	int idx = 0;
	vector<TrigSummand> CodeResult(2 * m*n);
	for (int i = 0; i < m; i++){
		for (int j = 0; j < n; j++){
			double a = CodeSin[i].coef * CodeCos[j].coef;
			int alpha = CodeSin[i].angle + CodeCos[j].angle;
			int beta = CodeSin[i].angle - CodeCos[j].angle;
			if (CodeSin[i].func == Sin){
				if (CodeCos[j].func == Sin){
					// sin(a) * sin(b)
					CodeResult[2 * idx] = TrigSummand(-a, Cos, alpha);
					CodeResult[2 * idx + 1] = TrigSummand(a, Cos, beta);
				}
				else{
					// sin(a) * cos(b)
					CodeResult[2 * idx] = TrigSummand(a, Sin, alpha);
					CodeResult[2 * idx + 1] = TrigSummand(a, Sin, beta);
				}
			}
			else{
				if (CodeCos[j].func == Sin){
					// cos(a) * sin(b)
					CodeResult[2 * idx] = TrigSummand(a, Sin, alpha);
					CodeResult[2 * idx + 1] = TrigSummand(-a, Sin, beta);
				}
				else{
					// cos(a) * cos(b)
					CodeResult[2 * idx] = TrigSummand(a, Cos, alpha);
					CodeResult[2 * idx + 1] = TrigSummand(a, Cos, beta);
				}
			}
			idx++;
		}
	}
	CollectSimilar(CodeResult);
	for (int i = 0; i < CodeResult.size(); i++)
		CodeResult[i].coef /= (2 * denoms*denomc);
	return CodeResult;

}


int ArrangementNumber(int n, int k)
{
	k = max(k, n - k);
	int res = 1;
	for (int i = k + 1; i <= n; i++)
		res *= i;
	for (int i = 2; i <= n - k; i++)
		res /= i;
	return res;
}


double GeneralizedSpectrum::PolygonMoment(vector<QPointF> Points){
	Points.push_back(Points[0]);
	double crv = 0;
	for (int i = 0; i < Points.size() - 1; i++){
		crv += (Points[i + 1].x()- Points[i].x()) * (Points[i + 1].y() + Points[i].y());
	}
	if (powx == 0 && powy == 0){
		return abs(crv) / 2;
	}
	if (abs(crv) < pow(10.0, 2*max(powx,powy)-6)){
		return 0;
	}

	double val = 0;
	double y0 = Points[0].y();

	for (int i = 0; i < Points.size() - 1; i++){
		
		double dx = Points[i + 1].x() - Points[i].x();
		double x1 = Points[i].x();
		double y1 = Points[i].y();
		double dy = Points[i + 1].y() - Points[i].y();
		double sum = 0;

		if (abs(dx) > pow(10.0, 2 * max(powx, powy) - 6)){
			double corr = 1.0 / (powy + 1);
			double a = dy / dx;
			double b = y1 - x1 * dy / dx;
			vector<double> coefs(powy + 2);
			for (int i = 0; i < powy + 2; i++){
				coefs[i] = ArrangementNumber(powy + 1, i) * pow(a, powy + 1 - i) * pow(b, i);;
			}
			coefs[powy + 1] -= pow(y0, powy + 1);

			for (int i = 0; i < powy + 2; i++){
				int powt = powx + powy + 2 - i;
				sum += coefs[i] * (pow(x1 + dx, powt) - pow(x1, powt)) / powt;
			}
			sum *= corr;

		}
		else{
			sum = dx * (y1 - y0 + dy / 2);
			sum *= dx * pow(x1 + dx / 2, powx) * pow((y0 + y1) / 2 + dy / 4, powy);
		}
		val += sum;
	}

	if (crv < 0)
		val *= -1;
	return val;
}


double GeneralizedSpectrum::BoneValue(TBone* Bone, double rad){
	if (rad > Bone->dest->r()){
		return 0;
	}
	if (rad > 0 && rad <= Bone->org->r())
		return Bone->s;

	vector<QPointF> Points;
	QPointF A = Bone->GetExtremePoint(rad);
	QPointF B = QPointF(Bone->dest->X(), Bone->dest->Y());
	Points.push_back(A);
	if (Bone->type == Linear){
		Points.push_back(ProjectPoint(A, (TEdge*)Bone->Sites[0]));
		Points.push_back(ProjectPoint(B, (TEdge*)Bone->Sites[0]));
		Points.push_back(B);
		Points.push_back(ProjectPoint(B, (TEdge*)Bone->Sites[1]));
		Points.push_back(ProjectPoint(A, (TEdge*)Bone->Sites[1]));
	}
	else if (Bone->type == Parabolic){
		Points.push_back(QPointF(((TVertex*)Bone->Sites[0])->p->X, ((TVertex*)Bone->Sites[0])->p->Y));
		Points.push_back(B);
		Points.push_back(ProjectPoint(B, (TEdge*)Bone->Sites[1]));
		Points.push_back(ProjectPoint(A, (TEdge*)Bone->Sites[1]));
	}
	else{
		Points.push_back(QPointF(((TVertex*)Bone->Sites[0])->p->X, ((TVertex*)Bone->Sites[0])->p->Y));
		Points.push_back(B);
		Points.push_back(QPointF(((TVertex*)Bone->Sites[1])->p->X, ((TVertex*)Bone->Sites[1])->p->Y));
	}
	double val = PolygonMoment(Points);
	if (rad == 0)
		Bone->s = val;

	return val;
}


double GeneralizedSpectrum::SectorMoment(QPointF A, double r, double phi0, double phi1){
	vector<double> CoefsX(powx + 1);
	for (int i = 0; i < powx + 1; i++){
		CoefsX[i] = ArrangementNumber(powx, i) * pow(A.x(), powx - i);
	}
	vector<double> CoefsY(powy + 1);
	for (int j = 0; j < powy + 1; j++){
		CoefsY[j] = ArrangementNumber(powy, j) * pow(A.y(), powy - j);
	}

	double val = 0;
	for (int i = 0; i < powx + 1; i++){
		for (int j = 0; j < powy + 1; j++){
			// Интегрируем r^(i+j+1) * (cos(a))^i * (sin(a))^j
			vector<GeneralizedSpectrum::TrigSummand> Temp = Code[j][i];
			double temp = 0;
			for (int k = 0; k < Temp.size(); k++){
				if (Temp[k].func == Sin){
					// Интегрируем синус
					if (Temp[k].angle > 0){
						temp -= Temp[k].coef / Temp[k].angle * (cos(Temp[k].angle * phi1) - cos(Temp[k].angle * phi0));
					}
				}
				else{
					// Интегрируем косинус
					if (Temp[k].angle > 0){
						temp += Temp[k].coef / Temp[k].angle * (sin(Temp[k].angle * phi1) - sin(Temp[k].angle * phi0));
					}
					else{
						temp += Temp[k].coef * (phi1 - phi0);
					}
				}
			}
			val += CoefsX[i] * CoefsY[j] / (i + j + 2) * pow(r, i + j + 2) * temp;
		}
	}
	return val;
}


double GeneralizedSpectrum::SectorValue(TBone* Bone, double rad){
	if (rad > Bone->dest->r())
		return 0;
	
	QPointF A = Bone->GetExtremePoint(rad);
	rad = max(rad, Bone->org->r());
	QPointF A1, A2;
	if (Bone->type == Linear){
		A1 = ProjectPoint(A, (TEdge*)Bone->Sites[0]);
		A2 = ProjectPoint(A, (TEdge*)Bone->Sites[1]);
	}
	else if (Bone->type == Parabolic){
		A1 = QPointF(((TVertex*)Bone->Sites[0])->p->X, ((TVertex*)Bone->Sites[0])->p->Y);
		A2 = ProjectPoint(A, (TEdge*)Bone->Sites[1]);
	}
	else{
		A1 = QPointF(((TVertex*)Bone->Sites[0])->p->X, ((TVertex*)Bone->Sites[0])->p->Y);
		A2 = QPointF(((TVertex*)Bone->Sites[1])->p->X, ((TVertex*)Bone->Sites[1])->p->Y);
	}
	double alpha1 = atan2(A1.y() - A.y(), A1.x() - A.x());
	double alpha2 = atan2(A2.y() - A.y(), A2.x() - A.x());
	if (alpha2 < alpha1){
		double temp = alpha2;
		alpha2 = alpha1;
		alpha1 = temp;
	}
	double val = (alpha2 - alpha1 > M_PI) ?
				 SectorMoment(A, rad, alpha2, alpha1 + 2 * M_PI) :
				 SectorMoment(A, rad, alpha1, alpha2);
	return val;
}


double GeneralizedSpectrum::LensValue(TBone* BoneA, TBone* BoneB, double rad){
	QPointF A = BoneA->GetExtremePoint(rad);
	QPointF B = BoneB->GetExtremePoint(rad);
	double len = sqrt(pow(B.y() - A.y(), 2) + pow(B.x() - A.x(), 2));
	if (len >= 2 * rad){
		return 0;
	}
	if (len == 0){
		return (SectorMoment(A, rad, 0, 2 * M_PI));
	}

	double h = sqrt(rad*rad - len*len / 4);
	double alpha = atan2(B.y() - A.y(), B.x() - A.x());
	double beta = asin(h / rad);

	vector<QPointF> Points;
	Points.push_back(A);
	Points.push_back(QPointF(A.x() + rad*cos(alpha + beta), A.y() + rad*sin(alpha + beta)));
	Points.push_back(B);
	Points.push_back(QPointF(A.x() + rad*cos(alpha - beta), A.y() + rad*sin(alpha - beta)));

	double val0 = SectorMoment(A, rad, alpha - beta, alpha + beta);
	double val1 = SectorMoment(B, rad, alpha + M_PI - beta, alpha + M_PI + beta);
	double val2 = PolygonMoment(Points);
	double val = val0 + val1 - val2;

	return val;
}