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

  const int nevtree = fInTree->GetEntries();

  // Every main event fills the histograms once, so fNev main events are needed.
  // Pile-up partners are taken from a separate pool of ~p*fNev tree entries.
  if (fNev < 0 || fNev + (int)std::ceil(fPars.p * fNev) > nevtree)
  {
    const int nevmax = (int)(nevtree / (1. + fPars.p));
    if (fNev >= 0)
      std::cerr << "\tToyMc::Run: WARNING: requested " << fNev << " events, but the input tree (" << nevtree
                << " entries) allows only " << nevmax << " events with p = " << fPars.p << std::endl;
    fNev = nevmax;
  }
  const int nplp = (fPars.p > 1e-10) ? std::max(1, (int)std::ceil(fPars.p * fNev)) : 0;

  vNpart.clear();
  vNcoll.clear();
  for (int i = 0; i < nevtree && (int)vNpart.size() < fNev + nplp; ++i)
  {
    if (fInTree->GetEntry(i) <= 0)
      continue;
    vNpart.push_back(fNpart);
    vNcoll.push_back(fNcoll);
  }
  if ((int)vNpart.size() < fNev + nplp)
  {
    std::cerr << "\tToyMc::Run: ERROR: could only read " << vNpart.size() << " entries from the input tree" << std::endl;
    return false;
  }

  // Main events: [0, fNev), pile-up pool: [fNev, fNev + nplp)
  std::atomic<int> plp_counter{0};
  std::vector<std::thread> v_thr;
  std::deque<std::atomic<int>> v_progress;
  const int n_part = fNev / fNthreads;
  const int n_rest = fNev % fNthreads;
  for (unsigned int i = 0; i < fNthreads; ++i)
  {
    const int i_start = i * n_part + std::min<int>(i, n_rest);
    const int i_stop = i_start + n_part + ((int)i < n_rest ? 1 : 0);
    v_progress.emplace_back(0);
    std::atomic<int> &progress = v_progress.back();
    v_thr.emplace_back([this, i_start, i_stop, nplp, &plp_counter, &progress]
                       { BuildMultiplicity(fPars.f, fPars.mu, fPars.k, fPars.p, i_start, i_stop, fNev, fNev + nplp, plp_counter, progress); });
  }

  int tot_progress = 0;
  while (tot_progress < fNev)
  {
    tot_progress = 0;
    for (auto &prog : v_progress)
      tot_progress += prog.load();
    std::cout << "\tToyMc::Run: Constructing multiplicity, progress: " << (fNev > 0 ? (long)tot_progress * 100 / fNev : 100) << "% \r" << std::flush;
    if (tot_progress < fNev)
      std::this_thread::sleep_for(std::chrono::milliseconds(200));
  }
  std::cout << "\t                                                                                                \r" << std::flush;

  for (auto &thread : v_thr)
    thread.join();

  return true;
}

bool ToyMc::BuildMultiplicity(double f, double mu, double k, double p, int i_start, int i_stop, int plp_start, int plp_stop, std::atomic<int> &plp_counter, std::atomic<int> &_progress)
{
  std::random_device rd;
  std::mt19937 rngnum(rd());
  std::uniform_real_distribution<float> unidist(0., 1.);
  std::gamma_distribution<> gammadst((float)((mu * k) / (mu + k)), (float)((k + mu) / k));
  std::negative_binomial_distribution<> nbddst(k, (float)(k / (k + mu)));
  const int n_plp = plp_stop - plp_start;
  for (int i = i_start; i < i_stop; i++)
  {
    const int Na = int(GetNacestors(f, vNpart.at(i), vNcoll.at(i)));
    float nHits{0}, nPlp{0};
    for (int j = 0; j < Na; j++)
      if (fUseNbd) nHits += nbddst(rngnum);
      else nHits += gammadst(rngnum);
    const bool isPileUp = n_plp > 0 && p > 1e-10 && unidist(rngnum) <= p;
    if (isPileUp)
    {
      // The pool holds ~p*Nev entries; wrap around if fluctuations exceed it
      const int iplp = plp_start + plp_counter.fetch_add(1) % n_plp;
      const int Na1 = int(GetNacestors(f, vNpart.at(iplp), vNcoll.at(iplp)));
      for (int jplp = 0; jplp < Na1; jplp++)
        if (fUseNbd) nPlp += nbddst(rngnum);
        else nPlp += gammadst(rngnum);
      nHits += nPlp;
    }
    const bool isTriggered = !isTrEff || unidist(rngnum) < hTrigEff.GetBinContent(hTrigEff.FindBin(nHits));
    if (isTriggered)
    {
      std::lock_guard<std::mutex> guard(fMtx);
      if (isPileUp)
        hMultPileUp.Fill(nHits);
      else
        hMultSingle.Fill(nHits);
      hMultAll.Fill(nHits);
    }
    _progress.store(i - i_start + 1);
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
