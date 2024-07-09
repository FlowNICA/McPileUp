#ifndef TOYMC_PILEUP_H
#define TOYMC_PILEUP_H

#include <iostream>

#include <TH1D.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <TRandom.h>

struct TMCParameters
{
  double f;
  int k;
  double mu;
  double p;
};

class ToyMc
{
private:
  TString fOutName;
  float fNpart, fNcoll, fB;
  std::unique_ptr<TTree> fInTree;
  TH1D hNbd, hMultAll, hMultPileUp, hMultSingle;
  TMCParameters fPars;
  bool isInputRead, isNbdInit;
  int fNev;
  std::function<int(double,double,double)> fNaFunc;
protected:
  double NBD(double n, double mu, double k);
  bool InitNbd();
  int GetNacestors(double f, double npart, double ncoll){ return fNaFunc(f,npart,ncoll); }
public:
  ToyMc();
  ToyMc(TMCParameters _pars);
  ToyMc(double f, int k, double mu, double p);
  virtual ~ToyMc();

  bool SetInput(std::unique_ptr<TTree> tree){ fInTree = std::move(tree); isInputRead = true; return true; }
  bool SetOutput(TString _outName){ fOutName = _outName;  return true; }
  bool SetParameters(TMCParameters _pars);
  bool SetParameters(double f, int k, double mu, double p);
  bool SetNevents(int _nev){ fNev = _nev; return true; }
  bool SetNancestors(std::function<int(double,double,double)> func){ fNaFunc = func; return true; } //wrapper for user-defined function

  bool Print();
  bool Run();
  bool Write();

ClassDef(ToyMc, 1);
};

#endif