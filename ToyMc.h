#ifndef TOYMC_PILEUP_H
#define TOYMC_PILEUP_H

#include <iostream>

#include <TH1D.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <TRandom.h>
#include <thread>
#include <mutex>
#include <random>
#include <boost/random.hpp>
#include <boost/math/distributions/negative_binomial.hpp>


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
  TH1D hNbd, hMultAll, hMultPileUp, hMultSingle, hTrigEff;
  TMCParameters fPars;
  bool isInputRead, isNbdInit, isTrEff;
  int fNev;
  std::function<int(double,double,double)> fNaFunc;
  std::mutex fMtx;
  unsigned int fNthreads{std::thread::hardware_concurrency()};
  std::vector<int> vNpart, vNcoll;
protected:
  double NBD(double n, double mu, double k);
  bool InitNbd();
  int GetNacestors(double f, double npart, double ncoll){ return fNaFunc(f,npart,ncoll); }
  bool BuildMultiplicity(double f, double mu, double k, double p, int i_start, int i_stop, int plp_start, int plp_stop);
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
  bool SetTriggerEfficiency(TH1D hist){ hTrigEff = hist; isTrEff = true; return true; }


  unsigned int GetNthreads(){ return fNthreads; }

  bool Print();
  bool Run();
  bool Write();

ClassDef(ToyMc, 1);
};

#endif