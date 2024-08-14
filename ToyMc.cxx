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
                 fNev{-1},
                 fUseNbd{true}
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
                                    fNev{-1},
                                    fUseNbd{true}
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
                                                     fNev{-1},
                                                     fUseNbd{true}
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

bool ToyMc::SetParameters(double f, double k, double mu, double p)
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
  std::deque<std::atomic<int>> v_progress;
  for (unsigned int i = 0; i < fNthreads; ++i)
  {
    int n_part = (int)(fNev / fNthreads);
    int i_start = i * n_part;
    int i_stop = (int)((i + 1) * n_part * (1. - fPars.p));
    int p_start = i_stop;
    int p_stop = (i + 1) * n_part;
    v_progress.emplace_back(0);
    v_thr.emplace_back([&]
                       { ToyMc::BuildMultiplicity(fPars.f, fPars.mu, fPars.k, fPars.p, i_start, i_stop, p_start, p_stop, std::ref(v_progress[i])); });
  }

  bool isOver = false;
  int tot_progress, tot_denum;
  while (not isOver)
  {
    isOver = true;
    tot_progress = 0;
    tot_denum = 0;
    std::cout << "\tToyMc::Run: Constructing multiplicity, progress: ";
    for (int j = 0; j < fNthreads; ++j)
    {
      if (v_progress[j].load() == 0)
        continue;
      tot_progress += v_progress[j].load();
      tot_denum++;
    }
    if (tot_denum > 0)
      tot_progress /= tot_denum;
    std::cout << tot_progress << "% \r" << std::flush;
    if (tot_progress < 100)
      isOver = false;
    std::chrono::milliseconds dura(200);
    std::this_thread::sleep_for(dura);
  }
  std::cout << "\t                                                                                                \r" << std::flush;

  for (auto &thread : v_thr)
    thread.join();

  return true;
}

bool ToyMc::BuildMultiplicity(double f, double mu, double k, double p, int i_start, int i_stop, int plp_start, int plp_stop, std::atomic<int> &_progress)
{
  std::random_device rd;
  std::mt19937 rngnum(rd());
  std::uniform_real_distribution<float> unidist(0., 1.);
  std::gamma_distribution<> gammadst((float)((mu * k) / (mu + k)), (float)((k + mu) / k));
  std::negative_binomial_distribution<> nbddst(k, (float)(k / (k + mu)));
  int plp_counter = plp_start;
  for (int i = i_start; i < i_stop; i++)
  {
    _progress.store((i - i_start)*100/(i_stop - i_start)+1);
    const int Na = int(GetNacestors(f, vNpart.at(i), vNcoll.at(i)));
    float nHits{0}, nPlp{0};
    for (int j = 0; j < Na; j++)
      if (fUseNbd) nHits += nbddst(rngnum);
      else nHits += gammadst(rngnum);
    if (p > 1e-10 && unidist(rngnum) <= p)
    {
      const int Na1 = int(GetNacestors(f, vNpart.at(plp_counter), vNcoll.at(plp_counter)));
      for (int jplp = 0; jplp < Na1; jplp++)
        if (fUseNbd) nPlp += nbddst(rngnum);
        else nPlp += gammadst(rngnum);
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
