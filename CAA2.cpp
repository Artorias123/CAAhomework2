#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include "stdlib.h"
#include <ctime>
#include <sstream>
using namespace std;
const extern double l[4][9] = {
-2.3395528692935242494659180754011,5.6918127138916832780830908380831,- 7.8140885026933122049392966796901,9.165830709219167905429952721561,- 8.4792368080250403900208941339517,6.0085925546698365084198824702185,- 3.1179230191579330200435611157182,1.0561270335769628265427765635503,- 0.17156181218784065400603258865185 
,- 0.26307611979449402425238442133944,- 0.81144175298238348565007432899297,1.5546661479781045337371308273191,- 0.67932391760754101828126235980106,4942731959861391 / 18014398509481984,- 0.098200019222674648268879803961991,0.027994733031129730943101358437611,- 0.0055724442702692134181234282515191,0.00057663154869108046757149858546757
,0.033251572576462796976079620803714,- 0.37255324396035053810597639096766,- 0.73915860973343755186257984994498,1.716734042016982429558226598483,- 1.0345376058234758244688350156325,0.59792433878984890926005266212578,- 0.27393456457054981371544130042141,0.085388111491920751843608554593919,- 0.01311404078740115948513487903994
,- 0.02842881626632891647237213178739,0.20598437315146262448165237017793,- 0.83671060994077123477773893446026,0,0.91455656078673415120763361912676,- 0.37478672705161777701408269676911,0.15801498482251149257425030876093,- 0.044092453922085889862942933906439,0.0054626884200955498636003988575741
},
l0[9] = { 0.0084819701130100958845200048987371, -0.063181723581980368709508822732053, 0.25124052653115384932616857527417, -0.84686376276840697606189070194714, 0, 0.84686376276840697606189070194714, -0.25124052653115384932616857527417, 0.063181723581980368709508822732053, -0.0084819701130100958845200048987371 },
d[9] = { 0.0037600069514406009906645074097254, 0.011870988142434283917852572352166, -0.0011126188868847911213293909163933, -0.26187098814243428391785257235217, 0.49470522387088838026132976701334, -0.26187098814243428391785257235217, -0.0011126188868847911213293909163933, 0.011870988142434283917852572352166, 0.0037600069514406009906645074097254 },
b[4] = { 0.159177852630552, 0.333333333333333, 0.348310961405563, 0.159177852630552 },
a[4] = { 0, 0.5, 0.5, 1 },
M = 0.5,
v = -0.1;
int mod(int a, int b)
{
	return ((a + b) % b);
}
int sgn(int a)
{
	if (a<0) return -1;
	else return 1;
}
void cal_v(vector<double> &v, double x, double y)//����
{
	v[3] = exp(-log(2)*((x - 100)*(x - 100) + (y - 100)*(y - 100)) / 9);
	v[1] = 0.04*(y - 100)*exp(-log(2)*(((x - 100) - 67)*((x - 100) - 67) + (y - 100)*(y - 100)) / 25);
	v[2] = -0.04*((x - 100) - 67)*exp(-log(2)*(((x - 100) - 67)*((x - 100) - 67) + (y - 100)*(y - 100)) / 25);
	v[0] = exp(-log(2)*((x - 100)*(x - 100) + (y - 100)*(y - 100)) / 9) + 0.1*exp(-log(2)*(((x - 100) - 67)*((x - 100) - 67) + (y - 100)*(y - 100)) / 25);
}

class grid
{
public:
	vector<vector<vector<double> > > g;
	double begin, end;
	grid(double xb, double xe, unsigned m);
	double space();
	void result(ofstream &fout, int m = 1);
	void prt(ofstream &fout, int m);
};
grid::grid(double xb, double xe, unsigned m)
{
	begin = xb;
	end = xe;
	g.resize(m + 1, vector<vector<double> >(m + 1, vector<double>(4)));
}
double grid::space()
{
	return(end - begin) / (g.size() - 1);
}
void grid::result(ofstream &fout, int m)
{
	fout<<"VARIABLES = \"X\", \"Y\",\"ρ\"，\"u\",\"v\",\"p\""<<endl<<"ZONE I="<<g.size()<<", J="<<g.size()<<", F=POINT"<<endl;
	for (unsigned i = 0; i < g[0].size(); i++)
	{
		for (unsigned j = 0; j<g.size(); j++)
		{
			fout<<double(i)*space()<<'\t'<<double(j)*space()<<'\t';
			for(unsigned k=0;k<4;k++) fout<<g[i][j][k]<<'\t';
			fout<<endl;
		}
	}
}
void grid::prt(ofstream &fout, int m)
{
	const char *a = "ruvp";
	fout << a[m] << endl;
	for (unsigned j = 0; j < g[0].size(); j++)
	{
		for (unsigned k = 0; k<g.size(); k++) fout << g[k][j][m] << '\t';
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
	vector<vector<vector<double> > > cal_k(vector<vector<vector<double> > > &k, unsigned ki);
	void cal_k1(vector<vector<vector<double> > > &y, vector<vector<vector<double> > > &k, unsigned ki);
	void cal_k2(vector<vector<vector<double> > > &y, vector<vector<vector<double> > > &k, unsigned ki);
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
void calculator::cal_k1(vector<vector<vector<double > > > &y, vector<vector<vector<double> > > &k, unsigned ki)
{
	int n = p->g.size(), tj;
	double ds = p->space(), dif[8], sth, cth, r, vth;
	for (int i = -4; i < 4; i++)
	{
		for (int j = -4; j < 4; j++)
		{
			if (i == -1)
			{
				r = pow((100 - M*t)*(100 - M*t) + (mod(j, n) - 100)*ds*ds*(mod(j, n) - 100), 0.5);
				sth = (mod(j, n) - 100)*ds / r;
				cth = (100 - M*t) / r;
				vth = (M*cth + pow(1 - M*M*sth*sth, 0.5));
			}
			if (i == 0)
			{
				r = pow((100 + M*t)*(100 + M*t) + (mod(j, n) - 100)*ds*ds*(mod(j, n) - 100), 0.5);
				sth = (mod(j, n) - 100)*ds / r;
				cth = -(100 + M*t) / r;
				vth = (M*cth + pow(1 - M*M*sth*sth, 0.5));
			}
			if (j == -1)
			{
				r = pow((i*ds - 100 + M*t)*(i*ds - 100 + M*t) + 100 * ds*ds * 100, 0.5);
				sth = 100 * ds / r;
				cth = (i*ds - 100 + M*t) / r;
				vth = (M*cth + pow(1 - M*M*sth*sth, 0.5));
			}
			if (j == 0)
			{
				r = pow((i*ds - 100 + M*t)*(i*ds - 100 + M*t) + 100 * ds*ds * 100, 0.5);
				sth = -100 * ds / r;
				cth = (i*ds - 100 + M*t) / r;
				vth = (M*cth + pow(1 - M*M*sth*sth, 0.5));
			}
			for (int m = 0; m < 8; m++) dif[m] = 0;
			for (int m = 0; m < 9; m++)
			{
				for (unsigned s = 0; s < 4; s++)
				{
					dif[s] += sgn(i)*l[(abs(2 * i + 1) - 1) / 2][m] * (p->g[mod(sgn(i)*m + (sgn(i) - 1) / 2, n)][mod(j, n)][s] + d_t*a[ki] * k[mod(sgn(i)*m + (sgn(i) - 1) / 2, n)][mod(j, n)][s]) / ds;
					dif[s + 4] += sgn(j)*l[(abs(2 * j + 1) - 1) / 2][m] * (p->g[mod(i, n)][mod(sgn(j)*m + (sgn(j) - 1) / 2, n)][s] + d_t*a[ki] * k[mod(i, n)][mod(sgn(j)*m + (sgn(j) - 1) / 2, n)][s]) / ds;
				}
			}
			if (i != -1 && i != 0 && j != -1 && j != 0)
			{
				y[mod(i, n)][mod(j, n)][0] = -(M*dif[0] + dif[1] + dif[6]);
				y[mod(i, n)][mod(j, n)][1] = -(M*dif[1] + dif[3]);
				y[mod(i, n)][mod(j, n)][2] = -(M*dif[2] + dif[7]);
				y[mod(i, n)][mod(j, n)][3] = -(M*dif[3] + dif[1] + dif[6]);
			}
			else
			{
				if (i == -1)
				{
					y[mod(i, n)][mod(j, n)][1] = -M*dif[1] - dif[3];
					y[mod(i, n)][mod(j, n)][2] = -M*dif[2] - dif[7];
					y[mod(i, n)][mod(j, n)][3] = -vth*(cth*dif[3] + sth*dif[7] + p->g[mod(i, n)][mod(j, n)][3] / (2 * r));
					y[mod(i, n)][mod(j, n)][0] = -M*dif[0] + M*dif[3] + y[mod(i, n)][mod(j, n)][3];
				}
				else
				{
					y[mod(i, n)][mod(j, n)][0] = -vth*(cth*dif[0] + sth*dif[4] + p->g[mod(i, n)][mod(j, n)][0] / (2 * r));
					y[mod(i, n)][mod(j, n)][1] = -vth*(cth*dif[1] + sth*dif[5] + p->g[mod(i, n)][mod(j, n)][1] / (2 * r));
					y[mod(i, n)][mod(j, n)][2] = -vth*(cth*dif[2] + sth*dif[6] + p->g[mod(i, n)][mod(j, n)][2] / (2 * r));
					y[mod(i, n)][mod(j, n)][3] = -vth*(cth*dif[3] + sth*dif[7] + p->g[mod(i, n)][mod(j, n)][3] / (2 * r));
				}
			}
		}
	}
}
void calculator::cal_k2(vector<vector<vector<double > > > &y, vector<vector<vector<double > > > &k, unsigned ki)
{
	int n = p->g.size();
	double ds = p->space(), dif[8], dss = ds*ds, sth, cth, r, vth;
	for (int i = -4; i < 4; i++)
	{
		for (int j = 4; j < n - 4; j++)
		{
			if (i == -1)
			{
				r = pow((100 - M*t)*(100 - M*t) + (j - 100)*ds*ds*(j - 100), 0.5);
				sth = -(j - 100)*ds / r;
				cth = (100 - M*t) / r;
				vth = (M*cth + pow(1 - M*M*sth*sth, 0.5));
			}
			if (i == 0)
			{
				r = pow((100 + M*t)*(100 + M*t) + (j - 100)*ds*ds*(j - 100), 0.5);
				sth = (j - 100)*ds / r;
				cth = -(100 + M*t) / r;
				vth = (M*cth + pow(1 - M*M*sth*sth, 0.5));
			}
			for (int m = 0; m < 8; m++) dif[m] = 0;
			for (int m = 0; m < 9; m++)
			{
				for (unsigned s = 0; s < 4; s++)
				{
					dif[s] += sgn(i)*l[(abs(2 * i + 1) - 1) / 2][m] * (p->g[mod(sgn(i)*m + (sgn(i) - 1) / 2, n)][j][s] + d_t*a[ki] * k[mod(sgn(i)*m + (sgn(i) - 1) / 2, n)][j][s]) / ds;
					dif[s + 4] += l0[m] * (p->g[mod(i, n)][j + m - 4][s] + d_t*a[ki] * k[mod(i, n)][j + m - 4][s]) / ds + v*d[m]*(p->g[mod(i, n)][j + m - 4][s] + d_t*a[ki] * k[mod(i, n)][j + m - 4][s]) / dss;
				}
			}
			if (i != -1 && i != 0)
			{
				y[mod(i, n)][j][0] = -(M*dif[0] + dif[1] + dif[6]);
				y[mod(i, n)][j][1] = -(M*dif[1] + dif[3]);
				y[mod(i, n)][j][2] = -(M*dif[2] + dif[7]);
				y[mod(i, n)][j][3] = -(M*dif[3] + dif[1] + dif[6]);
			}
			else
			{
				if (i == -1)
				{
					y[mod(i, n)][j][1] = -M*dif[1] - dif[3];
					y[mod(i, n)][j][2] = -M*dif[2] - dif[7];
					y[mod(i, n)][j][3] = -vth*(cth*dif[3] + sth*dif[7] + p->g[mod(i, n)][j][3] / (2 * r));
					y[mod(i, n)][j][0] = -M*dif[0] + M*dif[3] + y[mod(i, n)][j][3];
				}
				else
				{
					y[mod(i, n)][j][1] = -vth*(cth*dif[1] + sth*dif[5] + p->g[mod(i, n)][j][1] / (2 * r));
					y[mod(i, n)][j][2] = -vth*(cth*dif[2] + sth*dif[6] + p->g[mod(i, n)][j][2] / (2 * r));
					y[mod(i, n)][j][3] = -vth*(cth*dif[3] + sth*dif[7] + p->g[mod(i, n)][j][3] / (2 * r));
					y[mod(i, n)][j][0] = -vth*(cth*dif[0] + sth*dif[4] + p->g[mod(i, n)][j][0] / (2 * r));
				}
			}
		}
	}
	for (int i = 4; i < n - 4; i++)
	{
		for (int j = -4; j < 4; j++)
		{
			for (int m = 0; m < 8; m++) dif[m] = 0;
			for (int m = 0; m < 9; m++)
			{
				if (j == -1)
				{
					r = pow((i*ds - 100 + M*t)*(i*ds - 100 + M*t) + 100 * ds*ds * 100, 0.5);
					sth = 100 * ds / r;
					cth = (i*ds - 100 + M*t) / r;
					vth = (M*cth + pow(1 - M*M*sth*sth, 0.5));
				}
				if (j == 0)
				{
					r = pow((i*ds - 100 + M*t)*(i*ds - 100 + M*t) + 100 * ds*ds * 100, 0.5);
					sth = -100 * ds / r;
					cth = (i*ds - 100 + M*t) / r;
					vth = (M*cth + pow(1 - M*M*sth*sth, 0.5));
				}
				for (unsigned s = 0; s < 4; s++)
				{
					dif[s + 4] += sgn(j)*l[(abs(2 * j + 1) - 1) / 2][m] * (p->g[i][mod(sgn(j)*m + (sgn(j) - 1) / 2, n)][s] + d_t*a[ki] * k[i][mod(sgn(j)*m + (sgn(j) - 1) / 2, n)][s]) / ds;
					dif[s] += l0[m] * (p->g[i + m - 4][mod(j, n)][s] + d_t*a[ki] * k[i + m - 4][mod(j, n)][s]) / ds + v*d[m] * (p->g[i + m - 4][mod(j, n)][s] + d_t*a[ki] * k[i + m - 4][mod(j, n)][s]) / dss;
				}
			}
			if (j != -1 && j != 0)
			{
				y[i][mod(j, n)][0] = -(M*dif[0] + dif[1] + dif[6]);
				y[i][mod(j, n)][1] = -(M*dif[1] + dif[3]);
				y[i][mod(j, n)][2] = -(M*dif[2] + dif[7]);
				y[i][mod(j, n)][3] = -(M*dif[3] + dif[1] + dif[6]);
			}
			else
			{
				y[i][mod(j, n)][0] = -vth*(cth*dif[0] + sth*dif[4] + p->g[i][mod(j, n)][0] / (2 * r));
				y[i][mod(j, n)][1] = -vth*(cth*dif[1] + sth*dif[5] + p->g[i][mod(j, n)][1] / (2 * r));
				y[i][mod(j, n)][2] = -vth*(cth*dif[2] + sth*dif[6] + p->g[i][mod(j, n)][2] / (2 * r));
				y[i][mod(j, n)][3] = -vth*(cth*dif[3] + sth*dif[7] + p->g[i][mod(j, n)][3] / (2 * r));
			}
		}
	}
}
void calculator::cal_k3(vector<vector<vector<double > > > &y, vector<vector<vector<double> > > &k, unsigned ki)
{
	int n = p->g.size();
	double ds = p->space(), dif[6], dss = ds*ds;
	for (int i = 4; i<n - 4; i++)
	{
		for (int j = 4; j<n - 4; j++)
		{
			for (int m = 0; m<6; m++) dif[m] = 0;
			for (int m = 0; m<9; m++)
			{
				for (int s = 0; s<4; s++)
				{
					dif[s] += l0[m] * (p->g[i + m - 4][j][s] + d_t*a[ki] * k[i + m - 4][j][s]) / ds + v*d[m] * (p->g[i + m - 4][j][s] + d_t*a[ki] * k[i + m - 4][j][s]) / dss;
				}
				for (int s = 0; s<2; s++)
				{
					dif[s + 4] += l0[m] * (p->g[i][j + m - 4][s + 2] + d_t*a[ki] * k[i][j + m - 4][2 + s]) / ds + v*d[m] * (p->g[i][j + m - 4][2 + s] + d_t*a[ki] * k[i][j + m - 4][2 + s]) / dss;
				}
			}
			y[i][j][0] = -(M*dif[0] + dif[1] + dif[4]);
			y[i][j][1] = -(M*dif[1] + dif[3]);
			y[i][j][2] = -(M*dif[2] + dif[5]);
			y[i][j][3] = -(M*dif[3] + dif[1] + dif[4]);
		}
	}
}
vector<vector<vector<double> > > calculator::cal_k(vector<vector<vector<double> > > &k, unsigned ki)
{
	unsigned n = p->g.size();
	double d = p->space(), dif[6] = { 0 };
	vector<vector<vector<double> > > y(n, vector<vector<double> >(n, vector<double>(4)));
	cal_k1(y, k, ki);
	cal_k2(y, k, ki);
	cal_k3(y, k, ki);
	return y;
}
void calculator::advance()
{
	int n = p->g.size();
	vector<vector<vector<vector<double> > > > k(4, vector<vector<vector<double> > >(n, vector<vector<double> >(n, vector<double>(4))));
	k[0] = cal_k(p->g, 0);
	for (int i = 1; i<4; i++)
	{
		k[i] = cal_k(k[i - 1], i);
	}
	for (unsigned i = 0; i < n; i++)
	{
		for (unsigned j = 0; j<n; j++)
		{
			for (unsigned s = 0; s<4; s++) p->g[i][j][s] += d_t*(b[0] * k[0][i][j][s] + b[1] * k[1][i][j][s] + b[2] * k[2][i][j][s] + b[3] * k[3][i][j][s]);
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
	for (; step < total_step;)
	{
		advance();
		//x.result(c);
		step++;
		t = t + d_t;
		if (step % h == h - 1)
		{
			string a=to_string(int(t));
			a="t="+a+".dat";
			const char *name=a.data();
			ofstream fout(name, ios::out | ios::trunc);
			p->result(fout, c);
			//p->prt(fouts,0);
			fout.close();
		}
	}
}

int main()
{
	clock_t start_time = clock();
	grid g(0, 200, 400);
	calculator h(g, 0.1, 1000);
	for (unsigned i = 0; i < g.g.size(); i++)
	{
		for (unsigned j = 0; j < g.g.size(); j++)
		{
			cal_v(g.g[i][j], i*g.space(), j*g.space());
		}
	}
	h.iteration();
	clock_t end_time = clock();
	cout << "time=" << static_cast<double>(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;
	return 0;
}