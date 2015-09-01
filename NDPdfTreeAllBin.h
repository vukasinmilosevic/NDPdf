#ifndef NDPdfTreeAllBin_H
#define NDPdfTreeAllBin_H
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include <map>
#include <vector>
#include <string>
#include "TKDTree.h"
#include "TStopwatch.h" 

class NDPdfTreeAllBin
{

	std::vector<Double_t> fData;
	std::vector<Double_t> fWgt;
    	mutable TString fOptions;
	mutable Int_t fNDim;
	mutable Int_t fNEvents;
	mutable Double_t fSqrt2pi;
	mutable Double_t fD;
 	mutable Double_t fN;
	mutable Double_t fNEventsW;
	mutable Int_t fMode;
	mutable Int_t fNBins;

	//mutable std::vector<Double_t> fX;
 	mutable std::vector<Double_t> fX0, fX1, fX2;
	mutable std::vector<Double_t> fMean, fSigma;
	mutable Double_t fWidthFactor;
	mutable Double_t fNSigma;
	mutable Bool_t fRotate;
	mutable std::vector<Double_t> fRho;

	mutable TMatrixDSym* fCovMat;
	mutable TMatrixDSym* fCorrMat; 
  	mutable TMatrixD* fRotMat; 
 	mutable TVectorD* fSigmaR; 
  	mutable TVectorD* fDx;
 	mutable Double_t fSigmaAvgR;

	mutable std::vector<std::vector<Double_t> > fDataPts;
	mutable std::vector<TVectorD> fDataPtsR;
	mutable std::vector<std::vector<Double_t> > fWeights0;
	mutable std::vector<std::vector<Double_t> > fWeights1;
	mutable std::vector<std::vector<Double_t> >* fWeights;
        mutable TKDTreeID *fTree; 
	mutable std::vector<Double_t> fAvrP;
	mutable std::vector<Double_t> fAvrW;
	mutable Double_t **fTData;
	mutable std::vector<std::vector<Double_t> > fAvrPoints;
	mutable std::vector<std::vector<Double_t> > fAvrWeights;// Average bandwith per bin
	mutable std::vector<Double_t> fAvrWgt;// Average weight of data per bin
	mutable int fTotalNodes;
	mutable int fNNodes;
	mutable int fSigmaRange;

	
	
 	
public:
  	NDPdfTreeAllBin(std::vector<Double_t>& data, TString options, Double_t rho, Double_t nSigma, Bool_t rotate,Int_t nEvents, Int_t nDim, Int_t mode, Int_t nbins );
	NDPdfTreeAllBin(std::vector<Double_t>& data, TString options, Double_t rho, Double_t nSigma, Bool_t rotate,Int_t nEvents, Int_t nDim, std::vector<Double_t>& wgt, Int_t mode, Int_t nbins );
  	NDPdfTreeAllBin(std::vector<Double_t>& data, TString options, Double_t rho, Double_t nSigma, Bool_t rotate,Int_t nEvents, Int_t nDim, Int_t mode, Int_t nbins, int sigmarange);
	NDPdfTreeAllBin(std::vector<Double_t>& data, TString options, Double_t rho, Double_t nSigma, Bool_t rotate,Int_t nEvents, Int_t nDim, std::vector<Double_t>& wgt, Int_t mode, Int_t nbins, int sigmarange );
	Double_t evaluate(double *x) const;
protected:

	void     CreatePdf(Bool_t firstCall=kTRUE) const;
	void     SetOptions() const;
	void     Initialize() const;
	void     LoadDataSet(Bool_t firstCall) const;
	void     CalculateBandWidth() const; 
	void     AverageW() const; 
	Double_t Gauss(double *x, std::vector<std::vector<Double_t> >& weights) const;
	Double_t GaussB(double *x, std::vector<std::vector<Double_t> >& weights) const;
	Double_t GaussAll(double *x, std::vector<std::vector<Double_t> >& weights) const;
	std::vector<Double_t> FindAveragePoint(double *x) const;
	std::vector<Double_t> FindAverageWeight(double *x) const;

};

#endif
