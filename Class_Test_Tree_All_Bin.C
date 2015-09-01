#include "TROOT.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "Math/QuantFuncMathCore.h"
#include "TUUID.h"
#include "NDPdfTreeAllBin.h"
#include <iostream>
#include "TH2.h"


TRandom r;
TH2D *Hist;
TH2D *Hist_PDF;

Double_t X[50][50];
NDPdfTreeAllBin *A;





void Class_Test()
{
//_data=data;
Hist= new TH2D("Gaus_hist","Gaus_hist",50,80,120,50,80,120);
Hist_PDF= new TH2D("Gaus_hist_PDF","Gaus_hist_PDF",50,80,120,50,80,120);
std::vector<Double_t> data (20000);
std::vector<Double_t> points_x (10000);
std::vector<Double_t> points_y (10000);
std::vector<Double_t> weights (10000);

for (int i=0;i<10000;i++)
	{
		points_x[i]=r.Gaus(100,4);
		points_y[i]=r.Gaus(100,4);
		//std::cout<< points_x[i] <<endl;
		Hist->Fill(points_x[i],points_y[i]);
		data[2*i]=points_x[i];
		data[2*i+1]=points_y[i];
		weights[i]=r.Rndm();
		//weights[i]=1;
		
		
	}

int i_x,i_y;
Double_t eval;
TStopwatch w;
w.Start();
A= new NDPdfTreeAllBin(data,"a",1,3,false,10000,2,weights,1,10,5);
std::vector<double> X(2);
Double_t res=0;
Double_t res2=0;
w.Print();
w.Start();
for (int i=1; i<=50; i++)
{
for (int j=1;j<=50;j++)
{
  X[0] = Hist->GetXaxis()->GetBinCenter(i);
  X[1] = Hist->GetYaxis()->GetBinCenter(j);
  eval=A->evaluate(X.data());

//i_x=Hist->GetXaxis()->FindFixBin(X[0][i]);
//i_y=Hist->GetYaxis()->FindFixBin(X[1][i]);
//std::cout<<" i = " << i << " j = " << j << " eval = " << eval << " hist = " << Hist->GetBinContent(i,j)<< endl;
Hist_PDF->SetBinContent(i,j,eval);
Double_t g= TMath::Gaus(X[0],100,4, true)*TMath::Gaus(X[1],100,4,true)*10000*Hist->GetXaxis()->GetBinWidth(i)*Hist->GetYaxis()->GetBinWidth(j);
res+=eval-g;

if (g>0)
res2+=(eval-g)*(eval-g)/g;

}

}
w.Print();
cout <<"res=" << res << " res2= " <<res2<< "  " <<res2/2500 <<endl;
Hist->Draw("LEGO");
new TCanvas();
Hist_PDF->Draw("LEGO");


}
