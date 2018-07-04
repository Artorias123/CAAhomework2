#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include "stdlib.h"
#include <ctime>
using namespace std;
const extern double l[8][9] = {
	-1.0336594583151133048243995098609
	, 0.77334821015497883117383786684437
	, 0.10354329807491612738380850467365
	, 0.74470711700014343175142216908243
	, -0.75110975327374195273069575866623
	, -0.029089141516036408518894870307539
	, 0.79585772263400179410210097436473
	, -0.98896618563908585440636516561927
	, 0.38536819087993733606918578948877

	, -0.2688380322165280943113689951085
	, -0.59449085160088538876314135931884
	, 0.81676121363372113654145871100435
	, 0.30756527936720761587895264958188
	, -0.35137322062920934132278002891091
	, 0.023392092993567030485697368639718
	, 0.22276483229804155431248131539579
	, -0.23646652075988925760150865890276
	, 0.080685206913974744780208997619259

	, -0.1021167113291081406501145229248
	, 0.16600371560591905586959532233472
	, -1.4276528361300363985085513114456
	, 1.7716604371684574992014446890587
	, -0.30128708115648028162014569669466
	, -0.38332686956837030772559420615599
	, 0.52908115540069136345318445259839
	, -0.34051086514098804809117014951459
	, 0.08814905514991525807135142274386

	, -0.13206786354661002733709454493709
	, 0.44089854638509570495155846632431
	, -0.81323315762686436993460782207888
	, -0.40571411852638473275723616381851
	, 1.2547768197276136039904437659796
	, -0.56628503931084116545784212236297
	, 0.36486345445172571047826220600743
	, -0.17175644781106217311509826205209
	, 0.02851780625732744918161447693817

	, 1.03365945832, -0.773348210155, -0.103543298075, -0.744707117, 0.751109753274, 0.029089141516, -0.795857722634, 0.988966185639, -0.38536819088,
	0.268838032217, 0.594490851601, -0.816761213634, -0.307565279367, 0.351373220629, -0.0233920929936, -0.222764832298, 0.23646652076, -0.080685206914,
	0.102116711329, -0.166003715606, 1.42765283613, -1.77166043717, 0.301287081156, 0.383326869568, -0.529081155401, 0.340510865141, -0.0881490551499,
	0.132067863547, -0.440898546385, 0.813233157627, 0.405714118526, -1.25477681973, 0.566285039311, -0.364863454452, 0.171756447811, -0.0285178062573,
},
l0[9] = { 0.0084819701130100958845200048987371, -0.063181723581980368709508822732053, 0.25124052653115384932616857527417, -0.84686376276840697606189070194714, 0, 0.84686376276840697606189070194714, -0.25124052653115384932616857527417, 0.063181723581980368709508822732053, -0.0084819701130100958845200048987371 },
d[9] = { 0.0037600069514406009906645074097254, 0.011870988142434283917852572352166, -0.0011126188868847911213293909163933, -0.26187098814243428391785257235217, 0.49470522387088838026132976701334, -0.26187098814243428391785257235217, -0.0011126188868847911213293909163933, 0.011870988142434283917852572352166, 0.0037600069514406009906645074097254 },
b[4] = { 0.159177852630552, 0.333333333333333, 0.348310961405563, 0.159177852630552 },
a[4] = { 0, 0.5, 0.5, 1 },
M = 0.5,
v = -0.055,
sigma0 = 24,
belta=2;
int mod(int a, int b)
{
	return ((a + b) % b);
}
int sgn(int a)
{
	if (a<0) return -1;
	else return 1;
}
bool judge(int i,int p)
{
	if (i<p&&i>-p - 1) return 1;
	else return 0;
}
void cal_v(vector<double> &v, double x, double y)//����
{
	v[3] = exp(-log(2)*((x - 100)*(x - 100) + (y - 100)*(y - 100)) / 9);
	v[1] = 0.04*(y - 100)*exp(-log(2)*(((x - 100) - 67)*((x - 100) - 67) + (y - 100)*(y - 100)) / 25);
	v[2] = -0.04*((x - 100) - 67)*exp(-log(2)*(((x - 100) - 67)*((x - 100) - 67) + (y - 100)*(y - 100)) / 25);
	v[0] = exp(-log(2)*((x - 100)*(x - 100) + (y - 100)*(y - 100)) / 9) + 0.1*exp(-log(2)*(((x - 100) - 67)*((x - 100) - 67) + (y - 100)*(y - 100)) / 25);
	/*
	v[3] = exp(-log(2)*((x - 100)*(x - 100) + (y - 100)*(y - 100)) / 9);
	v[1] = 0;
	v[2] = 0;
	v[0] = exp(-log(2)*((x - 100)*(x - 100) + (y - 100)*(y - 100)) / 9);
	*/
}

class grid
{
public:
	vector<vector<vector<double> > > g;
	vector<vector<vector<double> > > q;
	double begin, end;
	grid(double xb, double xe, unsigned m);
	double space();
	void result(ofstream &fout, int m = 1);
	void prt(ofstream &fout, int m);
	void extent();
};
grid::grid(double xb, double xe, unsigned m)
{
	begin = xb;
	end = xe;
	double ds = (end - begin) / m;
	g.resize(m + 1 + 40 / ds, vector<vector<double> >(m + 1 + 40 / ds, vector<double>(4, 0)));
	q.resize(m + 1 + 40 / ds, vector<vector<double> >(m + 1 + 40 / ds, vector<double>(4, 0)));
}
double grid::space()
{
	return(end - begin + 40) / (g.size() - 1);
}
void grid::result(ofstream &fout, int m)
{
	int h = 1;
	fout << "VARIABLES = \"X\", \"Y\",\"��\"��\"u\",\"v\",\"p\"" << endl << "ZONE I=" << (g.size() - 9) / h + 1 << ", J=" << (g.size() - 9) / h + 1 << ", F=POINT" << endl;
	for (unsigned i = 20 / space(); i < g[0].size() - 20 / space(); i += h)
	{
		for (unsigned j = 20 / space(); j<g.size() - 20 / space(); j += h)
		{
			fout << double(i - 20 / space())*space() << '\t' << double(j - 20 / space())*space() << '\t';
			for (unsigned k = 0; k<4; k++) fout << g[i][j][k] << '\t';
			fout << endl;
		}
	}
}
void grid::prt(ofstream &fout, int m)
{
	const char *a = "ruvp";
	fout << a[m] << endl;
	int h = 1;
	for (unsigned j = 0; j < g[0].size(); j = j + h)
	{
		for (unsigned k = 0; k<g.size(); k = k + h) fout << g[k][j][m] << '\t';
		fout << endl;
	}
	fout << endl;
}
class calculator
{
public:
	grid *p;
	unsigned step, total_step;
	double d_t, t;
	calculator(grid &x, double dt, unsigned m, double it = 0);
	void cal_k(vector<vector<vector<double> > > &y, vector<vector<vector<double> > > &yq, vector<vector<vector<double> > > &k, vector<vector<vector<double> > > &kq, unsigned ki);
	void cal_k1(vector<vector<vector<double > > > &y, vector<vector<vector<double > > > &yq, vector<vector<vector<double> > > &k, vector<vector<vector<double> > > &kq, unsigned ki);
	void cal_k2(vector<vector<vector<double > > > &y, vector<vector<vector<double > > > &yq, vector<vector<vector<double> > > &k, vector<vector<vector<double> > > &kq, unsigned ki);
	void cal_k3(vector<vector<vector<double> > > &y, vector<vector<vector<double> > > &k, unsigned ki);
	void advance();
	void iteration();
};
calculator::calculator(grid &x, double dt, unsigned m, double it)
{
	p = &x;
	d_t = dt;
	t = it;
	total_step = m;
	step = 0;
}
void calculator::cal_k1(vector<vector<vector<double > > > &y, vector<vector<vector<double > > > &yq, vector<vector<vector<double> > > &k, vector<vector<vector<double> > > &kq, unsigned ki)
{
	int n = p->g.size();
	double ds = p->space(), dif[2][6], sigma[2] = { 0, 0 }, difq[2][6], dss = ds*ds;
	for (int i = -20 / ds; i<20 / ds; i++)
	{
		for (int j = -20 / ds; j<20 / ds; j++)
		{
			for (int m = 0; m < 6; m++) { dif[0][m] = 0; dif[1][m] = 0; difq[0][m]=0; difq[1][m]=0; }
			for (int m = 0; m<9; m++)
			{
				for (int s = 0; s<4; s++)
				{
					dif[0][s] += l0[m] * (p->g[mod(i+m-4,n)][mod(j,n)][s] + d_t*a[ki] * k[mod(i+m-4,n)][mod(j,n)][s]) / ds + v*d[m] * (p->g[mod(i+m-4,n)][mod(j,n)][s] + d_t*a[ki] * k[mod(i+m-4,n)][mod(j,n)][s]) / dss;
					dif[1][s] += l0[m] * (p->g[mod(i+m-4,n)][mod(j,n)][s] + d_t*a[ki] * k[mod(i+m-4,n)][mod(j,n)][s]) / ds - v*d[m] * (p->g[mod(i+m-4,n)][mod(j,n)][s] + d_t*a[ki] * k[mod(i+m-4,n)][mod(j,n)][s]) / dss;
				}
				for (int s = 0; s<2; s++)
				{
					dif[0][s + 4] += l0[m] * (p->g[mod(i,n)][mod(j+m-4,n)][s + 2] + d_t*a[ki] * k[mod(i,n)][mod(j+m-4,n)][2 + s]) / ds + v*d[m] * (p->g[mod(i,n)][mod(j+m-4,n)][2 + s] + d_t*a[ki] * k[mod(i,n)][mod(j+m-4,n)][2 + s]) / dss;
					dif[1][s + 4] += l0[m] * (p->g[mod(i,n)][mod(j+m-4,n)][s + 2] + d_t*a[ki] * k[mod(i,n)][mod(j+m-4,n)][2 + s]) / ds - v*d[m] * (p->g[mod(i,n)][mod(j+m-4,n)][2 + s] + d_t*a[ki] * k[mod(i,n)][mod(j+m-4,n)][2 + s]) / dss;
				}
			}
			y[mod(i,n)][mod(j,n)][0] = -(M*dif[0][0] + (M + 1)*dif[0][1] / 2 + (1 - M)*dif[0][3] / 2 + (1 - M)*dif[1][1] / 2 + (M - 1)*dif[1][3] / 2 + 0.5*dif[0][4] + 0.5*dif[0][5] + 0.5*dif[1][4] - 0.5*dif[1][5]);
			y[mod(i,n)][mod(j,n)][1] = -((M + 1)*dif[0][1] / 2 + (M + 1)*dif[0][3] / 2 + (M - 1)*dif[1][1] / 2 + (1 - M)*dif[1][3] / 2);
			y[mod(i,n)][mod(j,n)][2] = -(M*dif[0][2] + 0.5*dif[0][4] + 0.5*dif[0][5] - 0.5*dif[1][4] + 0.5*dif[1][5]);
			y[mod(i, n)][mod(j, n)][3] = -((1 + M)*dif[0][1] / 2 + (M + 1)*dif[0][3] / 2 + (1 - M)*dif[1][1] / 2 + (M - 1)*dif[1][3] / 2 + 0.5*dif[0][4] + 0.5*dif[0][5] + 0.5*dif[1][4] - 0.5*dif[1][5]);
		}
	}
}
void calculator::cal_k2(vector<vector<vector<double > > > &y, vector<vector<vector<double > > > &yq, vector<vector<vector<double> > > &k, vector<vector<vector<double> > > &kq, unsigned ki)
{
	int n = p->g.size();
	double ds = p->space(), dif[2][6], sigma[2] = {0,0}, difq[2][6], dss = ds*ds;
	for (int i = -20 / ds; i<20 / ds; i++)
	{
		for (int j = 20 / ds; j<n - 20 / ds; j++)
		{
			sigma[0] = sigma0*(1 - M*M)*pow((20 / ds - sgn(i)*i + (sgn(i) - 1) / 2)/(20/ds), belta);
			sigma[1] = 0;
			for (int m = 0; m < 6; m++) { dif[0][m] = 0; dif[1][m] = 0; difq[0][m]=0; difq[1][m]=0; }
			for (int m = 0; m<9; m++)
			{
				for (int s = 0; s<4; s++)
				{
					dif[0][s] += l0[m] * (p->g[mod(i+m-4,n)][j][s] + d_t*a[ki] * k[mod(i+m-4,n)][j][s]) / ds + v*d[m] * (p->g[mod(i+m-4,n)][j][s] + d_t*a[ki] * k[mod(i+m-4,n)][j][s]) / dss;
					dif[1][s] += l0[m] * (p->g[mod(i+m-4,n)][j][s] + d_t*a[ki] * k[mod(i+m-4,n)][j][s]) / ds - v*d[m] * (p->g[mod(i+m-4,n)][j][s] + d_t*a[ki] * k[mod(i+m-4,n)][j][s]) / dss;
					difq[0][s] += l0[m] * (p->q[mod(i + m - 4, n)][j][s] + d_t*a[ki] * kq[mod(i + m - 4, n)][j][s]) / ds + v*d[m] * (p->q[mod(i + m - 4, n)][j][s] + d_t*a[ki] * kq[mod(i + m - 4, n)][j][s]) / dss;
					difq[1][s] += l0[m] * (p->q[mod(i + m - 4, n)][j][s] + d_t*a[ki] * kq[mod(i + m - 4, n)][j][s]) / ds - v*d[m] * (p->q[mod(i + m - 4, n)][j][s] + d_t*a[ki] * kq[mod(i + m - 4, n)][j][s]) / dss;
				}
				for (int s = 0; s<2; s++)
				{
					dif[0][s + 4] += l0[m] * (p->g[mod(i,n)][j + m - 4][s + 2] + d_t*a[ki] * k[mod(i,n)][j + m - 4][2 + s]) / ds + v*d[m] * (p->g[mod(i,n)][j + m - 4][2 + s] + d_t*a[ki] * k[mod(i,n)][j + m - 4][2 + s]) / dss;
					dif[1][s + 4] += l0[m] * (p->g[mod(i,n)][j + m - 4][s + 2] + d_t*a[ki] * k[mod(i,n)][j + m - 4][2 + s]) / ds - v*d[m] * (p->g[mod(i,n)][j + m - 4][2 + s] + d_t*a[ki] * k[mod(i,n)][j + m - 4][2 + s]) / dss;
					difq[0][s + 4] += l0[m] * (p->q[mod(i, n)][j + m - 4][s + 2] + d_t*a[ki] * kq[mod(i, n)][j + m - 4][2 + s]) / ds + v*d[m] * (p->q[mod(i, n)][j + m - 4][2 + s] + d_t*a[ki] * kq[mod(i, n)][j + m - 4][2 + s]) / dss;
					difq[1][s + 4] += l0[m] * (p->q[mod(i, n)][j + m - 4][s + 2] + d_t*a[ki] * kq[mod(i, n)][j + m - 4][2 + s]) / ds - v*d[m] * (p->q[mod(i, n)][j + m - 4][2 + s] + d_t*a[ki] * kq[mod(i, n)][j + m - 4][2 + s]) / dss;
				}
			}
			y[mod(i,n)][j][0] = -(M*dif[0][0] + (M + 1)*dif[0][1] / 2 + (1 - M)*dif[0][3] / 2 + (1 - M)*dif[1][1] / 2 + (M - 1)*dif[1][3] / 2 + 0.5*dif[0][4] + 0.5*dif[0][5] + 0.5*dif[1][4] - 0.5*dif[1][5]);
			y[mod(i,n)][j][1] = -((M + 1)*dif[0][1] / 2 + (M + 1)*dif[0][3] / 2 + (M - 1)*dif[1][1] / 2 + (1 - M)*dif[1][3] / 2);
			y[mod(i,n)][j][2] = -(M*dif[0][2] + 0.5*dif[0][4] + 0.5*dif[0][5] - 0.5*dif[1][4] + 0.5*dif[1][5]);
			y[mod(i,n)][j][3] = -((1 + M)*dif[0][1] / 2 + (M + 1)*dif[0][3] / 2 + (1 - M)*dif[1][1] / 2 + (M - 1)*dif[1][3] / 2 + 0.5*dif[0][4] + 0.5*dif[0][5] + 0.5*dif[1][4] - 0.5*dif[1][5]);

			y[mod(i, n)][j][0] += -sigma[0] * (M*difq[0][0] + (M + 1)*difq[0][1] / 2 + (1 - M)*difq[0][3] / 2 + (1 - M)*difq[1][1] / 2 + (M - 1)*difq[1][3] / 2) - sigma[1] * (0.5*difq[0][4] + 0.5*difq[0][5] + 0.5*difq[1][4] - 0.5*difq[1][5]);
			y[mod(i, n)][j][1] += -sigma[0] * ((M + 1)*difq[0][1] / 2 + (M + 1)*difq[0][3] / 2 + (M - 1)*difq[1][1] / 2 + (1 - M)*difq[1][3] / 2);
			y[mod(i, n)][j][2] += -sigma[0] * (M*difq[0][2]) - sigma[1] * (0.5*difq[0][4] + 0.5*difq[0][5] - 0.5*difq[1][4] + 0.5*difq[1][5]);
			y[mod(i, n)][j][3] += -sigma[0] * ((1 + M)*difq[0][1] / 2 + (M + 1)*difq[0][3] / 2 + (1 - M)*difq[1][1] / 2 + (M - 1)*difq[1][3] / 2) - sigma[1] * (0.5*difq[0][4] + 0.5*difq[0][5] + 0.5*difq[1][4] - 0.5*difq[1][5]);

			y[mod(i, n)][j][0] += -((sigma[0] + sigma[1])*p->g[mod(i, n)][j][0] + (sigma[0] * sigma[1])*p->q[mod(i, n)][j][0]);
			y[mod(i, n)][j][1] += -((sigma[0] + sigma[1])*p->g[mod(i, n)][j][1] + (sigma[0] * sigma[1])*p->q[mod(i, n)][j][1]);
			y[mod(i, n)][j][2] += -((sigma[0] + sigma[1])*p->g[mod(i, n)][j][2] + (sigma[0] * sigma[1])*p->q[mod(i, n)][j][2]);
			y[mod(i, n)][j][3] += -((sigma[0] + sigma[1])*p->g[mod(i, n)][j][3] + (sigma[0] * sigma[1])*p->q[mod(i, n)][j][3]);

			yq[mod(i, n)][j][0] = p->g[mod(i, n)][j][0];
			yq[mod(i, n)][j][1] = p->g[mod(i, n)][j][1];
			yq[mod(i, n)][j][2] = p->g[mod(i, n)][j][2];
			yq[mod(i, n)][j][3] = p->g[mod(i, n)][j][3];
		}
	}
	for (int i = 20 / ds; i<n - 20 / ds; i++)
	{
		for (int j = -20 / ds; j<20 / ds; j++)
		{
			sigma[1] = sigma0*pow((20 / ds - sgn(j)*j + (sgn(j) - 1) / 2) / (20 / ds), belta);
			sigma[0] = 0;
			for (int m = 0; m < 6; m++) { dif[0][m] = 0; dif[1][m] = 0; difq[0][m]=0; difq[1][m]=0; }
			for (int m = 0; m<9; m++)
			{
				for (int s = 0; s<4; s++)
				{
					dif[0][s] += l0[m] * (p->g[i+m-4][mod(j,n)][s] + d_t*a[ki] * k[i+m-4][mod(j,n)][s]) / ds + v*d[m] * (p->g[i+m-4][mod(j,n)][s] + d_t*a[ki] * k[i+m-4][mod(j,n)][s]) / dss;
					dif[1][s] += l0[m] * (p->g[i+m-4][mod(j,n)][s] + d_t*a[ki] * k[i+m-4][mod(j,n)][s]) / ds - v*d[m] * (p->g[i+m-4][mod(j,n)][s] + d_t*a[ki] * k[i+m-4][mod(j,n)][s]) / dss;
					difq[0][s] += l0[m] * (p->q[i + m - 4][mod(j, n)][s] + d_t*a[ki] * kq[i + m - 4][mod(j, n)][s]) / ds + v*d[m] * (p->q[i + m - 4][mod(j, n)][s] + d_t*a[ki] * kq[i + m - 4][mod(j, n)][s]) / dss;
					difq[1][s] += l0[m] * (p->q[i + m - 4][mod(j, n)][s] + d_t*a[ki] * kq[i + m - 4][mod(j, n)][s]) / ds - v*d[m] * (p->q[i + m - 4][mod(j, n)][s] + d_t*a[ki] * kq[i + m - 4][mod(j, n)][s]) / dss;
				}
				for (int s = 0; s<2; s++)
				{
					dif[0][s + 4] += l0[m] * (p->g[i][mod(j+m-4,n)][s + 2] + d_t*a[ki] * k[i][mod(j+m-4,n)][2 + s]) / ds + v*d[m] * (p->g[i][mod(j+m-4,n)][2 + s] + d_t*a[ki] * k[i][mod(j+m-4,n)][2 + s]) / dss;
					dif[1][s + 4] += l0[m] * (p->g[i][mod(j+m-4,n)][s + 2] + d_t*a[ki] * k[i][mod(j+m-4,n)][2 + s]) / ds - v*d[m] * (p->g[i][mod(j+m-4,n)][2 + s] + d_t*a[ki] * k[i][mod(j+m-4,n)][2 + s]) / dss;
					difq[0][s + 4] += l0[m] * (p->q[i][mod(j + m - 4, n)][s + 2] + d_t*a[ki] * kq[i][mod(j + m - 4, n)][2 + s]) / ds + v*d[m] * (p->q[i][mod(j + m - 4, n)][2 + s] + d_t*a[ki] * kq[i][mod(j + m - 4, n)][2 + s]) / dss;
					difq[1][s + 4] += l0[m] * (p->q[i][mod(j + m - 4, n)][s + 2] + d_t*a[ki] * kq[i][mod(j + m - 4, n)][2 + s]) / ds - v*d[m] * (p->q[i][mod(j + m - 4, n)][2 + s] + d_t*a[ki] * kq[i][mod(j + m - 4, n)][2 + s]) / dss;
				}
			}
			y[i][mod(j,n)][0] = -(M*dif[0][0] + (M + 1)*dif[0][1] / 2 + (1 - M)*dif[0][3] / 2 + (1 - M)*dif[1][1] / 2 + (M - 1)*dif[1][3] / 2 + 0.5*dif[0][4] + 0.5*dif[0][5] + 0.5*dif[1][4] - 0.5*dif[1][5]);
			y[i][mod(j,n)][1] = -((M + 1)*dif[0][1] / 2 + (M + 1)*dif[0][3] / 2 + (M - 1)*dif[1][1] / 2 + (1 - M)*dif[1][3] / 2);
			y[i][mod(j,n)][2] = -(M*dif[0][2] + 0.5*dif[0][4] + 0.5*dif[0][5] - 0.5*dif[1][4] + 0.5*dif[1][5]);
			y[i][mod(j,n)][3] = -((1 + M)*dif[0][1] / 2 + (M + 1)*dif[0][3] / 2 + (1 - M)*dif[1][1] / 2 + (M - 1)*dif[1][3] / 2 + 0.5*dif[0][4] + 0.5*dif[0][5] + 0.5*dif[1][4] - 0.5*dif[1][5]);

			y[i][mod(j, n)][0] += -sigma[0] * (M*difq[0][0] + (M + 1)*difq[0][1] / 2 + (1 - M)*difq[0][3] / 2 + (1 - M)*difq[1][1] / 2 + (M - 1)*difq[1][3] / 2) - sigma[1] * (0.5*difq[0][4] + 0.5*difq[0][5] + 0.5*difq[1][4] - 0.5*difq[1][5]);
			y[i][mod(j,n)][1] += -sigma[0] * ((M + 1)*difq[0][1] / 2 + (M + 1)*difq[0][3] / 2 + (M - 1)*difq[1][1] / 2 + (1 - M)*difq[1][3] / 2);
			y[i][mod(j,n)][2] += -sigma[0] * (M*difq[0][2]) - sigma[1] * (0.5*difq[0][4] + 0.5*difq[0][5] - 0.5*difq[1][4] + 0.5*difq[1][5]);
			y[i][mod(j,n)][3] += -sigma[0] * ((1 + M)*difq[0][1] / 2 + (M + 1)*difq[0][3] / 2 + (1 - M)*difq[1][1] / 2 + (M - 1)*difq[1][3] / 2) - sigma[1] * (0.5*difq[0][4] + 0.5*difq[0][5] + 0.5*difq[1][4] - 0.5*difq[1][5]);

			y[i][mod(j,n)][0] += -((sigma[0] + sigma[1])*p->g[i][mod(j,n)][0] + (sigma[0] * sigma[1])*p->q[i][mod(j,n)][0]);
			y[i][mod(j,n)][1] += -((sigma[0] + sigma[1])*p->g[i][mod(j,n)][1] + (sigma[0] * sigma[1])*p->q[i][mod(j,n)][1]);
			y[i][mod(j,n)][2] += -((sigma[0] + sigma[1])*p->g[i][mod(j,n)][2] + (sigma[0] * sigma[1])*p->q[i][mod(j,n)][2]);
			y[i][mod(j,n)][3] += -((sigma[0] + sigma[1])*p->g[i][mod(j,n)][3] + (sigma[0] * sigma[1])*p->q[i][mod(j,n)][3]);

			yq[i][mod(j,n)][0] = p->g[i][mod(j,n)][0];
			yq[i][mod(j,n)][1] = p->g[i][mod(j,n)][1];
			yq[i][mod(j,n)][2] = p->g[i][mod(j,n)][2];
			yq[i][mod(j,n)][3] = p->g[i][mod(j,n)][3];
		}
	}
}
void calculator::cal_k3(vector<vector<vector<double > > > &y, vector<vector<vector<double> > > &k, unsigned ki)
{
	int n = p->g.size();
	double ds = p->space(), dif[2][6], dss = ds*ds;
	for (int i = 20 / ds; i<n - 20/ds; i++)
	{
		for (int j = 20/ds; j<n - 20/ds; j++)
		{
			for (int m = 0; m < 6; m++) { dif[0][m] = 0; dif[1][m] = 0; }
			for (int m = 0; m<9; m++)
			{
				for (int s = 0; s<4; s++)
				{
					dif[0][s] += l0[m] * (p->g[i + m - 4][j][s] + d_t*a[ki] * k[i + m - 4][j][s]) / ds + v*d[m] * (p->g[i + m - 4][j][s] + d_t*a[ki] * k[i + m - 4][j][s]) / dss;
					dif[1][s] += l0[m] * (p->g[i + m - 4][j][s] + d_t*a[ki] * k[i + m - 4][j][s]) / ds - v*d[m] * (p->g[i + m - 4][j][s] + d_t*a[ki] * k[i + m - 4][j][s]) / dss;
				}
				for (int s = 0; s<2; s++)
				{
					dif[0][s + 4] += l0[m] * (p->g[i][j + m - 4][s + 2] + d_t*a[ki] * k[i][j + m - 4][2 + s]) / ds + v*d[m] * (p->g[i][j + m - 4][2 + s] + d_t*a[ki] * k[i][j + m - 4][2 + s]) / dss;
					dif[1][s + 4] += l0[m] * (p->g[i][j + m - 4][s + 2] + d_t*a[ki] * k[i][j + m - 4][2 + s]) / ds - v*d[m] * (p->g[i][j + m - 4][2 + s] + d_t*a[ki] * k[i][j + m - 4][2 + s]) / dss;
				}
			}
			y[i][j][0] = -(M*dif[0][0] + (M + 1)*dif[0][1] / 2 + (1 - M)*dif[0][3] / 2 + (1 - M)*dif[1][1] / 2 + (M - 1)*dif[1][3] / 2 + 0.5*dif[0][4] + 0.5*dif[0][5] + 0.5*dif[1][4] - 0.5*dif[1][5]);
			y[i][j][1] = -((M + 1)*dif[0][1] / 2 + (M + 1)*dif[0][3] / 2 + (M - 1)*dif[1][1] / 2 + (1 - M)*dif[1][3] / 2);
			y[i][j][2] = -(M*dif[0][2] + 0.5*dif[0][4] + 0.5*dif[0][5] - 0.5*dif[1][4] + 0.5*dif[1][5]);
			y[i][j][3] = -((1 + M)*dif[0][1] / 2 + (M + 1)*dif[0][3] / 2 + (1 - M)*dif[1][1] / 2 + (M - 1)*dif[1][3] / 2 + 0.5*dif[0][4] + 0.5*dif[0][5] + 0.5*dif[1][4] - 0.5*dif[1][5]);
		}
	}
}
void calculator::cal_k(vector<vector<vector<double> > > &y, vector<vector<vector<double> > > &yq, vector<vector<vector<double> > > &k, vector<vector<vector<double> > > &kq, unsigned ki)
{
	unsigned n = p->g.size();
	double d = p->space();
	cal_k1(y, yq, k, kq, ki);
	cal_k2(y, yq, k, kq, ki);
	cal_k3(y, k, ki);
}
void calculator::advance()
{
	int n = p->g.size();
	vector<vector<vector<vector<double> > > > k(4, vector<vector<vector<double> > >(n, vector<vector<double> >(n, vector<double>(4,0))));
	vector<vector<vector<vector<double> > > > kq(4, vector<vector<vector<double> > >(n, vector<vector<double> >(n, vector<double>(4,0))));
	cal_k(k[0],kq[0],p->g,p->q, 0);
	for (int i = 1; i<4; i++)
	{
		cal_k(k[i],kq[i],k[i - 1],kq[i-1], i);
	}
	for (unsigned i = 4; i < n-4; i++)
	{
		for (unsigned j = 4; j<n-4; j++)
		{
			for (unsigned s = 0; s<4; s++) p->g[i][j][s] += d_t*(b[0] * k[0][i][j][s] + b[1] * k[1][i][j][s] + b[2] * k[2][i][j][s] + b[3] * k[3][i][j][s]);
			if ((i<int(20 / p->space()) && i>n-int(20 / p->space())) || (j<int(20 / p->space()) && j>n-int(20 / p->space()))) 
			for (unsigned s = 0; s<4; s++) p->q[i][j][s] += d_t*(b[0] * kq[0][i][j][s] + b[1] * kq[1][i][j][s] + b[2] * kq[2][i][j][s] + b[3] * kq[3][i][j][s]);
		}
	}
}
void calculator::iteration()
{
	int h = 2 / d_t, c = 0.1 / p->space(), k = 0;
	if (c <= 1)
	{
		c = 1;
	}
	ofstream fout("prt.dat", ios::out | ios::trunc);
	for (; step < total_step;)
	{
		if (abs(t - 10 * k)< d_t / 2)
		{
			/*string a = to_string(t);
			a = "result_tecplot\\t=" + a + ".dat";
			const char *name = a.data();
			ofstream fout(name, ios::out | ios::trunc);
			p->result(fout, c);
			*/
			p->prt(fout,0);
			k++;
		}
		advance();
		//x.result(c);
		step++;
		t = t + d_t;
	}
	fout.close();
}

int main()
{
	clock_t start_time = clock();
	grid g(0, 200, 400);
	calculator h(g, 0.05, 20 * 101);
	for (unsigned i = 4; i < g.g.size()-4; i++)
	{
		for (unsigned j = 4; j < g.g.size()-4; j++)
		{
			cal_v(g.g[i][j], (i - 20 / g.space())*g.space(), (j - 20 / g.space())*g.space());
			//for (unsigned k = 0; k < 4; k++)
			//g.g[i][j][k] = 0.00000001;
		}
	}
	h.iteration();
	clock_t end_time = clock();
	cout << "time=" << double(end_time - start_time) / CLOCKS_PER_SEC << endl;
	system("pause");
}
	