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
TH2D *Hist_Interpol;
TH2D *DIST;

Double_t X[50][50];
NDPdfTreeAllBin *A;





void test()
{
    
    ROOT::Math::Random<ROOT::Math::GSLRngMT> r;
    TH1::SetDefaultSumw2();
//_data=data;
Hist= new TH2D("Gaus_hist","Gaus_hist",100,-10,10, 100,-10,10);
Hist_PDF= new TH2D("Gaus_hist_PDF","Gaus_hist_PDF",100,-10,10, 100,-10,10);
    Hist_Interpol= new TH2D("Gaus_hist_Interpol","Gaus_hist_Interpol",100,-10,10, 100,-10,10);
std::vector<Double_t> data (20000);
std::vector<Double_t> points_x (10000);
std::vector<Double_t> points_y (10000);
std::vector<Double_t> weights (10000);
    double sigmax = 2;
    double sigmay = 3;
    double rho = 0.5;
    double x,y;
for (int i=0;i<10000;i++)
	{
        r.Gaussian2D(sigmax,sigmay,0.5,x,y);
		//std::cout<< points_x[i] <<endl;
        data[2*i]=x;
        data[2*i+1]=y;
		
		//weights[i]=r.Rndm();
		weights[i]=1;
		Hist->Fill(x,y,weights[i]);
		
	}

int i_x,i_y;
    Double_t eval,interpol;
TStopwatch w;
    w.Start();
A= new NDPdfTreeAllBin(data,"a",1,3,true,10000,2,weights,3,10,1);
std::vector<double> X(2);
Double_t res=0;
Double_t res2=0;
w.Print();
    cout<<"Evauation&Interpolation time: "<<endl;
w.Start();
for (int i=1; i<=100; i++)
{
for (int j=1;j<=100;j++)
{
  X[0] = Hist->GetXaxis()->GetBinCenter(i);
  X[1] = Hist->GetYaxis()->GetBinCenter(j);
  eval=A->evaluate(X.data());
  interpol=A->Interpol(X.data());
//i_x=Hist->GetXaxis()->FindFixBin(X[0][i]);
//i_y=Hist->GetYaxis()->FindFixBin(X[1][i]);
//std::cout<<" i = " << i << " j = " << j << " eval = " << eval << " hist = " << Hist->GetBinContent(i,j)<< endl;
    Hist_PDF->SetBinContent(i,j,eval);
    Hist_Interpol->SetBinContent(i,j,interpol);
//Double_t g= TMath::Gaus(X[0],100,4, true)*TMath::Gaus(X[1],100,4,true)*10000*Hist->GetXaxis()->GetBinWidth(i)*Hist->GetYaxis()->GetBinWidth(j);
    Double_t g=Hist->GetBinContent(i,j);
res+=eval-g;

if (g>0)
res2+=(eval-g)*(eval-g)/g;
  //  cout<<"TEST"<<res2<<endl;
}

}
w.Print();
    
    Hist->Scale(1/Hist->Integral());
    Hist_PDF->Scale(1/Hist_PDF->Integral());
    Hist_Interpol->Scale(1/Hist_Interpol->Integral());
    Hist->SetLineColor(kRed);
    Hist_PDF->SetLineColor(kBlue);
    Hist_Interpol->SetLineColor(kMagenta);
//cout <<"res=" << res << " res2= " <<res2<< "  " <<res2/2500 <<endl;
Hist->Draw("LEGO");
new TCanvas();
Hist_PDF->Draw("LEGO");
new TCanvas();
Hist_Interpol->Draw("LEGO");
DIST= new TH2D("Gaus_DIST","Gaus_DIST",100,-10,10, 100,-10,10);
    DIST->Add(Hist,Hist_PDF,1,-1);
  
    for (int i=1; i<=100; i++)
        for (int j=1;j<=100;j++)
    {
        DIST->SetBinContent(i,j,(Hist->GetBinContent(i,j)-Hist_PDF->GetBinContent(i,j))/Hist->GetBinError(i,j));
    }
    
new TCanvas();
    DIST->Draw("COLZ");
}
