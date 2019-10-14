//一般的な蔵本モデルプログラム

#include<iostream>
#include<fstream>
#include<iomanip>
#include<string>
#include<vector>
#include<cmath>
#include<complex>

using namespace std;
#include"MT.h"

void RungeKutta(int ndim, double dt, double tt, vector<double>& fin, vector<double>& fout, void(*df)(int, double, vector<double>& , vector<double>&));
void function(int ndim, double tt, vector<double>& fin, vector<double>& df);
double sfunc(double theta);
double tfunc(double theta);

const double pi=M_PI;
double K;
vector <double> omega;

int main()
{
	int ndim;//次元
	//int nstep;//1周期における刻み幅(時間間隔)の数
	int nend;//最終時間間隔
	double dt;//刻み幅

	vector<double> fp;//初期値
	vector<double> fn;//解
	complex<double> R;//オーダーパラメータ
	complex<double> fsum;
	complex<double> i(0.0,1.0);//純虚数

	ifstream fin("input.dat");//入力ファイルを開く
	string str;
	ofstream fout("output1.dat");//出力ファイルを開く
	ofstream ffout("syokiti.dat");//出力ファイルを開く

	if (!fin)
	{
		cout<<"入力ファイルをオープンできません。"<<endl;
		return 1;
	}
	if(!fout)
	{
		cout<<"出力ファイルをオープンできません。"<<endl;
		return 1;
	}
	if(!ffout)
	{
		cout<<"出力ファイルをオープンできません。"<<endl;
		return 1;
	}

	while(getline(fin, str)){
		sscanf(str.data(), "%d, %d, %lf, %lf", &ndim, &nend, &dt, &K);
		cout<<"n(次元)="<<ndim<<endl;
		//cout<<"nstep(1周期における時間間隔の数)="<<nstep<<endl;
		cout<<"nend(最終時間間隔)="<<nend<<endl;
		cout<<"dt(刻み幅)"<<dt<<endl;
		//cout<<"omega0(位相周波数)="<<omega0<<endl;
		cout<<"K(振動子間の相互作用定数)="<<K<<endl;

		/*double tpr;
		tpr=2.0*pi/omega0;
		dt=tpr/(double)nstep;
		cout<<"dt(刻み幅)="<<dt<<endl;*/
	}

	fin.close();

	fp.resize(ndim);
	omega.resize(ndim);

	//乱数による振動数
	for(int n=0; n<ndim; n++){
		omega[n]=genrand_real1()*2-1;//[-1,1]
		cout<<"omega["<<n<<"](位相振動数)="<<omega[n]<<endl;
	}
	
	//乱数による初期位相
	int iseed=124;
	init_genrand(iseed);
	for(int n=0; n<ndim; n++){
		fp[n]=2*pi*genrand_real2()-pi;//[-pi, pi]
		cout<<"fp["<<n<<"](初期値)="<<fp[n]<<endl;
	}

	fout<<setw(5)<<"t";
	ffout<<setw(5)<<"t";
	for(int n=0; n<ndim; n++){
		fout<<setw(13)<<"f"<<n+1;
	}
	fout<<setw(14)<<"R";
	for(int n=0; n<ndim; n++){
		fout<<setw(25)<<"exp(i*fp["<<n+1<<"])";
		ffout<<setw(25)<<"exp(i*fp["<<n+1<<"])";
	}
	fout<<endl;
	ffout<<endl;

	double t0=0.0;//初期時刻
	fout<<setw(5)<<t0;
	ffout<<setw(5)<<t0;

	for(int n=0; n<ndim; n++){
		fout<<setw(14)<<fp[n];
	}
	//オーダーパラメータの定義
	fsum=0.0;
	for(int n=0; n<ndim; n++){
		fsum+=exp(i*fp[n]);
	}
	R=fsum/(double)ndim;
	fout<<setw(14)<<abs(R);
	for(int n=0; n<ndim; n++){
		fout<<setw(14)<<real(exp(i*fp[n]))<<setw(14)<<imag(exp(i*fp[n]));
		ffout<<setw(14)<<real(exp(i*fp[n]))<<setw(14)<<imag(exp(i*fp[n]));
	}
	fout<<endl;
	ffout<<endl;

	double tt;
	fn.resize(ndim);

	for(int nt=0; nt<nend; nt++){
		tt=t0+nt*dt;
		RungeKutta(ndim, dt, tt, fp, fn, &function);
		for(int n=0; n<ndim; n++){
			fn[n]=tfunc(fn[n]);
			fp[n]=fn[n];
		}

		fout<<setw(5)<<tt+dt;
		for(int n=0; n<ndim; n++){
			fout<<setw(14)<<fp[n];
		}
		//オーダーパラメータの定義
		fsum=0.0;
		for(int n=0; n<ndim; n++){
			fsum+=exp(i*fp[n]);
		}
		R=fsum/(double)ndim;
		fout<<setw(14)<<abs(R);

		for(int n=0; n<ndim; n++){
			fout<<setw(14)<<real(exp(i*fp[n]))<<setw(14)<<imag(exp(i*fp[n]));
		}
		fout<<endl;
	}

	fout.close();

	return 0;
}


void RungeKutta(int ndim, double dh, double hh, vector<double>& fin, vector<double>& fout, void(*df)(int, double, vector<double>& , vector<double>&))
{
	vector<double>df1;
	vector<double>df2;
	vector<double>df3;
	vector<double>df4;
	df1.resize(ndim);
	df2.resize(ndim);
	df3.resize(ndim);
	df4.resize(ndim);

	vector<double>ff;
	ff.resize(ndim);

	(*df)(ndim, hh, fin, df1);
	for(int i=0; i<ndim; i++){
		//k1[i]=dh*df1[i];
		ff[i]=fin[i]+0.5*dh*df1[i];
	}

	(*df)(ndim, hh+0.5*dh, ff, df2);
	for(int i=0; i<ndim; i++){
		//k2[i]=dh*df2[i];
		ff[i]=fin[i]+0.5*dh*df2[i];
	}

	(*df)(ndim, hh+0.5*dh, ff, df3);
	for(int i=0; i<ndim; i++){
		//k3[i]=dh*df3[i];
		ff[i]=fin[i]+dh*df3[i];
	}

	(*df)(ndim, hh+dh, ff, df4);
	for(int i=0; i<ndim; i++){
		//k4[i]=dh*df4[i];
		fout[i]=fin[i]+dh*(df1[i]/6.0+df2[i]/3.0+df3[i]/3.0+df4[i]/6.0);
	}
}

void function(int ndim, double tt, vector<double>& fin, vector<double>& df)
{
	/*vector<double>omega;
	omega.resize(ndim);

	for(int n=0; n<ndim; n++){
		omega[n]=omega0;
	}*/
	for(int n=0; n<ndim; n++){
		double couple=0.0;
		for(int m=0; m<ndim; m++){
			couple+=sfunc(fin[m]-fin[n]);
		}
		df[n]=omega[n]+K*couple/(double)ndim;
	}
}

double sfunc(double theta)
{
	double sfunc;

	sfunc=sin(theta);

	return sfunc;
}

double tfunc(double theta)
{
	double tfunc=theta;

	while(tfunc>pi){
		tfunc-=2.0*pi;
	}
	while(tfunc<-pi){
		tfunc+=2.0*pi;
	}
	
	return tfunc;
}

