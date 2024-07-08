#include "ToyMc.h"

ClassImp(ToyMc);

ToyMc::ToyMc() : fInTree{nullptr},
                fOutName{},
                hNbd{"hNbd", "N_{part};N_{part};counts", 10000, 0., 10000.},
                hMultAll{"hMultAll", "N_{ch};N_{ch};counts", 10000, 0., 10000.},
                hMultPileUp{"hMultPileUp", "N_{ch};N_{ch};counts", 10000, 0., 10000.},
                hMultSingle{"hMultSingle", "N_{ch};N_{ch};counts", 10000, 0., 10000.},
                fPars{0.5, 10, 10.5, 0.},
                isInputRead{false},
                isNbdInit{false},
                isInit{false},
                fNev{-1}
{}

ToyMc::ToyMc(TMCParameters _pars) : fInTree{nullptr},
                fOutName{},
                hNbd{"hNbd", "N_{part};N_{part};counts", 10000, 0., 10000.},
                hMultAll{"hMultAll", "N_{ch};N_{ch};counts", 10000, 0., 10000.},
                hMultPileUp{"hMultPileUp", "N_{ch};N_{ch};counts", 10000, 0., 10000.},
                hMultSingle{"hMultSingle", "N_{ch};N_{ch};counts", 10000, 0., 10000.},
                fPars{},
                isInputRead{false},
                isNbdInit{false},
                isInit{false},
                fNev{-1}
{
  SetParameters(_pars);
}

ToyMc::ToyMc(double f, int k, double mu, double p) : fInTree{nullptr},
                fOutName{},
                hNbd{"hNbd", "N_{part};N_{part};counts", 10000, 0., 10000.},
                hMultAll{"hMultAll", "N_{ch};N_{ch};counts", 10000, 0., 10000.},
                hMultPileUp{"hMultPileUp", "N_{ch};N_{ch};counts", 10000, 0., 10000.},
                hMultSingle{"hMultSingle", "N_{ch};N_{ch};counts", 10000, 0., 10000.},
                fPars{},
                isInputRead{false},
                isNbdInit{false},
                isInit{false},
                fNev{-1}
{
  SetParameters(f, k, mu, p);
}

ToyMc::~ToyMc()
{
  fInTree.reset();
}

bool ToyMc::SetParameters(TMCParameters _pars)
{
  fPars = _pars;
  return true;
}

bool ToyMc::SetParameters(double f, int k, double mu, double p)
{
  fPars = {f, k, mu, p};
  return true;
}

double ToyMc::NBD(double n, double mu, double k)
{
  double func1, func2, func;

  if (n+k > 100.0) 
  {
    // log method for handling large numbers
    func1  = TMath::LnGamma(n + k)- TMath::LnGamma(n + 1.)- TMath::LnGamma(k);
    func2  = n * TMath::Log(mu/k) - (n + k) * TMath::Log(1.0 + mu/k);
    func = TMath::Exp(func1 + func2);
  } 
  else 
  {
    func1  = TMath::Gamma(n + k) / ( TMath::Gamma(n + 1.) * TMath::Gamma(k) );
    func2  = n * TMath::Log(mu/k) - (n + k) * TMath::Log(1.0 + mu/k);
    func = func1 * TMath::Exp(func2);
  }
  return func;
}

bool ToyMc::ReadInput()
{
  return true;
}

bool ToyMc::InitNbd()
{
  if (isNbdInit)
    return false;

  int nbins = ((fPars.mu+1.)*3 < 10) ? 10 : (int)(fPars.mu+1.)*3;
  for (int i=0; i<nbins; ++i){
    double val = NBD(i, fPars.mu, fPars.k);
    if (val>1e-10)
      hNbd.SetBinContent(i+1, val);
  }
  isNbdInit = true;

  return true;
}

bool ToyMc::Run()
{
  if (!fInTree)
    return false;

  float fB{-1.}, fNpart{-1}, fNcoll{-1};

  fInTree->SetBranchAddress("B",  &fB);
  fInTree->SetBranchAddress("Npart",  &fNpart);
  fInTree->SetBranchAddress("Ncoll",  &fNcoll);
  
  InitNbd();
  if (!isNbdInit){
    std::cerr << "ToyMc::Run: ERROR: NBD distribution was incorrect!" << std::endl;
    return false;
  }
  
  int nevtree = fInTree->GetEntriesFast();
  if (fNev < 0 || fNev > nevtree*(1.-fPars.p)) fNev = nevtree*(1.-fPars.p);

  int mult, pileup;
  for (int i=0; i<fNev; ++i){
    mult = 0;
    pileup = 0;

    if (fInTree->GetEntry(i) <= 0)
      continue;
    
    for (int j=0; j<GetNacestors(fNpart, fNcoll); ++j)
      mult += (int)hNbd.GetRandom();

    // Generating pile-up
    if (fPars.p > 0 && gRandom->Rndm() < fPars.p){
      if (fInTree->GetEntry(i+fNev) <= 0) 
        continue;
      for (int j=0; j<GetNacestors(fNpart, fNcoll); ++j)
        pileup += (int)hNbd.GetRandom();
      hMultPileUp.Fill(mult+pileup);
      fInTree->GetEntry(i);
    } else {
      hMultSingle.Fill(mult);
    }

    hMultAll.Fill(mult + pileup);

    std::cout << "ToyMc::Run: event [" << i << "/" << fNev << "]:" 
      << " Na = " << GetNacestors(fNpart, fNcoll)
      << ", Npart = " << fNpart
      << ", Ncoll = " << fNcoll
      << ", mult = " << mult 
      << ", pile-up = " << pileup 
      << ", total = " << mult + pileup 
      << ".\r" << std::flush;
      // << "." << std::endl;
  }

  std::cout << std::endl;

  return true;
}

bool ToyMc::Write()
{
  std::unique_ptr<TFile> fo{new TFile(fOutName, "recreate")};
  fo->cd();

  hMultAll.Write();
  hMultPileUp.Write();
  hMultSingle.Write();
  hNbd.Write();

  fo->Close();
  fo.reset();

  return true;
}

bool ToyMc::Print()
{
  std::cout << "ToyMc::Print:" << std::endl;
  std::cout << "\tInput parameters:" << std::endl;
  std::cout << "\t\tf  : " << fPars.f << std::endl;
  std::cout << "\t\tk  : " << fPars.k << std::endl;
  std::cout << "\t\tmu : " << fPars.mu << std::endl;
  std::cout << "\t\tp  : " << fPars.p << std::endl;
  std::cout << "\tOutput file: " << fOutName << std::endl;

  return true;
}