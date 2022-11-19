#include "GeneralizedSpectrum.h"

int main(int argc, char *argv[])
{
	int area = 0;
	double step = 0.1;
	double powx = 0;
	double powy = 0;

	int i = 1;
	bool ready = (argc <= 2);
	while (!ready){
		if (strcmp(argv[i], "-a") == 0){
			bool ok;
			int temp = QString(argv[i+1]).toInt(&ok);
			if (!ok){
				printf("Invalid ignorable area value\n");
				return -1;
			}
			else if (temp < 0){
				printf("Ignorable area value can't be negative\n");
				return -1;
			}
			else area = temp;
			i = i+2;
		}
		else if (strcmp(argv[i], "-s") == 0){
			bool ok;
			double temp = QString(argv[i+1]).toDouble(&ok);
			if (!ok){
				printf("Invalid step area value\n");
				return -1;
			}
			else if (temp <= 0){
				printf("Step value must be area positive\n");
				return -1;
			}
			else step = temp;
			i = i+2;
		}
		else if (strcmp(argv[i], "-x") == 0){
			bool ok;
			int temp = QString(argv[i + 1]).toInt(&ok);
			if (!ok){
				printf("Invalid moment power\n");
				return -1;
			}
			else if (temp < 0){
				printf("Moment power can't be negative\n");
				return -1;
			}
			else powx = temp;
			i = i + 2;
		}
		else if (strcmp(argv[i], "-y") == 0){
			bool ok;
			int temp = QString(argv[i + 1]).toInt(&ok);
			if (!ok){
				printf("Invalid moment power\n");
				return -1;
			}
			else if (temp < 0){
				printf("Moment power can't be negative\n");
				return -1;
			}
			else powy = temp;
			i = i + 2;
		}
		else ready = true;
		if (argc <= i+1)
			ready = true;
	}

	GeneralizedSpectrum Spectrum(area, step, powx, powy);

	ifstream file((argc > i) ? argv[i] : "files.txt");
	char* path = new char[1024];
	clock_t start = clock();
	while (file.getline(path, 1024)){
		if (path[0] != '\0')
			Spectrum.ProcessImage(QString(path));
	}
	clock_t finish = clock();
	printf("Time spent: %f\n", (finish - start) / double(CLOCKS_PER_SEC));

	return 0;
}
