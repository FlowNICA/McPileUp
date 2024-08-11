#include "ToyMc.h"

ClassImp(ToyMc);

ToyMc::ToyMc() : fInTree{nullptr},
                 fOutName{},
                 fNpart{-1.},
                 fNcoll{-1.},
                 fB{-1.},
                 hNbd{"hNbd", "N_{part};N_{part};counts", 10000, 0., 10000.},
                 hMultAll{"hMultAll", "All;N_{ch};counts", 10000, 0., 10000.},
                 hMultPileUp{"hMultPileUp", "Pile-up;N_{ch};counts", 10000, 0., 10000.},
                 hMultSingle{"hMultSingle", "Single;N_{ch};counts", 10000, 0., 10000.},
                 hTrigEff{},
                 fPars{0.5, 10, 10.5, 0.},
                 isInputRead{false},
                 isNbdInit{false},
                 isTrEff{false},
                 fNev{-1}
{
}

ToyMc::ToyMc(TMCParameters _pars) : fInTree{nullptr},
                                    fOutName{},
                                    fNpart{-1.},
                                    fNcoll{-1.},
                                    fB{-1.},
                                    hNbd{"hNbd", "N_{part};N_{part};counts", 10000, 0., 10000.},
                                    hMultAll{"hMultAll", "All;N_{ch};counts", 10000, 0., 10000.},
                                    hMultPileUp{"hMultPileUp", "Pile-up;N_{ch};counts", 10000, 0., 10000.},
                                    hMultSingle{"hMultSingle", "Single;N_{ch};counts", 10000, 0., 10000.},
                                    hTrigEff{},
                                    fPars{},
                                    isInputRead{false},
                                    isNbdInit{false},
                                    isTrEff{false},
                                    fNev{-1}
{
  SetParameters(_pars);
}

ToyMc::ToyMc(double f, int k, double mu, double p) : fInTree{nullptr},
                                                     fOutName{},
                                                     fNpart{-1.},
                                                     fNcoll{-1.},
                                                     fB{-1.},
                                                     hNbd{"hNbd", "N_{part};N_{part};counts", 10000, 0., 10000.},
                                                     hMultAll{"hMultAll", "All;N_{ch};counts", 10000, 0., 10000.},
                                                     hMultPileUp{"hMultPileUp", "Pile-up;N_{ch};counts", 10000, 0., 10000.},
                                                     hMultSingle{"hMultSingle", "Single;N_{ch};counts", 10000, 0., 10000.},
                                                     hTrigEff{},
                                                     fPars{},
                                                     isInputRead{false},
                                                     isNbdInit{false},
                                                     isTrEff{false},
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

bool ToyMc::Run()
{
  if (!fInTree)
    return false;

  fInTree->SetBranchAddress("B", &fB);
  fInTree->SetBranchAddress("Npart", &fNpart);
  fInTree->SetBranchAddress("Ncoll", &fNcoll);

  float NpartMax = fInTree->GetMaximum("Npart");
  float NcollMax = fInTree->GetMaximum("Ncoll");

  int nevtree = fInTree->GetEntriesFast();
  if (fNev < 0 || fNev > nevtree)
    fNev = nevtree;

  int mult, pileup, plp_count = 0;
  for (int i = 0; i < fNev; ++i)
  {
    mult = 0;
    pileup = 0;

    if (fInTree->GetEntry(i) <= 0)
      continue;
    vNpart.push_back(fNpart);
    vNcoll.push_back(fNcoll);
  }

  std::vector<std::thread> v_thr;
  for (unsigned int i = 0; i < fNthreads; ++i)
  {
    int n_part = (int)(fNev / fNthreads);
    int i_start = i * n_part;
    int i_stop = (int)((i + 1) * n_part * (1. - fPars.p));
    int p_start = i_stop;
    int p_stop = (i + 1) * n_part;
    v_thr.emplace_back([&]
                       { ToyMc::BuildMultiplicity(fPars.f, fPars.mu, fPars.k, fPars.p, i_start, i_stop, p_start, p_stop); });
  }
  for (auto &thread : v_thr)
    thread.join();

  return true;
}

bool ToyMc::BuildMultiplicity(double f, double mu, double k, double p, int i_start, int i_stop, int plp_start, int plp_stop)
{

  boost::mt19937 rngnum;
  boost::random::negative_binomial_distribution<int> nbd(k, (double)(k / (k + mu)));
  std::uniform_real_distribution<double> unidist(0., 1.);
  int plp_counter = plp_start;
  for (int i = i_start; i < i_stop; i++)
  {
    const int Na = int(GetNacestors(f, vNpart.at(i), vNcoll.at(i)));
    int nHits{0}, nPlp{0};
    for (int j = 0; j < Na; j++)
      nHits += nbd(rngnum);
    if (p > 1e-10 && unidist(rngnum) <= p)
    {
      const int Na1 = int(GetNacestors(f, vNpart.at(plp_counter), vNcoll.at(plp_counter)));
      for (int jplp = 0; jplp < Na1; jplp++)
        nPlp += (int)nbd(rngnum);
      plp_counter++;
      nHits += nPlp;
      std::lock_guard<std::mutex> guard(fMtx);
      if (!isTrEff)
        hMultPileUp.Fill(nHits);
      if (isTrEff && unidist(rngnum) < hTrigEff.GetBinContent(hTrigEff.FindBin(nHits)))
        hMultPileUp.Fill(nHits);
    }
    else
    {
      std::lock_guard<std::mutex> guard(fMtx);
      if (!isTrEff)
        hMultSingle.Fill(nHits);
      if (isTrEff && unidist(rngnum) < hTrigEff.GetBinContent(hTrigEff.FindBin(nHits)))
        hMultSingle.Fill(nHits);
    }
    if (!isTrEff)
      hMultAll.Fill(nHits);
    if (isTrEff && unidist(rngnum) < hTrigEff.GetBinContent(hTrigEff.FindBin(nHits)))
      hMultAll.Fill(nHits);
  }
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
  if (isTrEff)
    hTrigEff.Write();

  fo->Close();
  fo.reset();

  return true;
}

bool ToyMc::Print()
{
  std::cout << "ToyMc::Print:" << std::endl;
  std::cout << "\tNumber of threads: " << fNthreads << std::endl;
  std::cout << "\tNumber of events:  " << fNev << std::endl;
  std::cout << "\tInput parameters:" << std::endl;
  std::cout << "\t\tf  : " << fPars.f << std::endl;
  std::cout << "\t\tk  : " << fPars.k << std::endl;
  std::cout << "\t\tmu : " << fPars.mu << std::endl;
  std::cout << "\t\tp  : " << fPars.p << std::endl;
  std::cout << "\tOutput file: " << fOutName << std::endl;

  return true;
}
