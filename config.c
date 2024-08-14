void config()
{
  // Set up input file from MC-Glauber
  std::string inFileName = "../example_files/out_au3au3_30mb_500k.root";
  std::string inTreeName = "nt_Au3_Au3";

  // Set up outpu file
  std::string outFileName = "./out_test.root";

  // Set up file with trigger efficiency
  std::string inTrigFileName = "../example_files/trigger_eff_example.root";
  std::string inTrigHistName = "hEff";

  // Set up number of processed events
  int Nev = 500000;

  // Set up parameters for multiplicity sampling
  double f  = 0.9;
  double mu = 0.5;
  double k  = 1;

  // Set up ratio of pile-up events: 
  // (1. -> 100% of pile-up events, 0. -> 0% of pile-up events)
  double p  = 0.03;

  // Set Nancestors function
  // It should be set as func(double f, double npart, double ncoll)
  // Examples of the Na parametrisations:
  //    Default : f*npart + (1-f)*ncoll
  //    STAR    : f*ncoll + (1-f)*npart*0.5
  //    PSD     : f - npart
  //    Npart   : pow(npart, f)
  //    Ncoll   : pow(ncoll, f)
  // One can set a lambda function such as:
  //    auto funcNa{[](double f, double npart, double ncoll){ return (int)(f*npart + (1-f)*ncoll); }};
  // or a regular function before config():
  //    int funcNa(double f, double npart, double ncoll){ return (int)(f*npart + (1-f)*ncoll); }
  // And then set it using ToyMc::SetNancestors:
  //    mc.SetNancestors(funcNa);
  auto Na{[](double f, double npart, double ncoll){ return (int)(f*npart + (1-f)*ncoll); }};
  // auto Na{[](double f, double npart, double ncoll){ return (int)(npart); }};
  
  // Read input
  std::unique_ptr<TFile> fi{new TFile(inFileName.c_str(), "read")};
  if (!fi){
    std::cerr << "\nERROR: cannot read input file!" << std::endl;
    return;
  }
  std::unique_ptr<TTree> tree{(TTree*)fi->Get(inTreeName.c_str())};
  if (!tree){
    std::cerr << "\nERROR: cannot read input tree!" << std::endl;
    return;
  }
  std::unique_ptr<TFile> fiTrEff{new TFile(inTrigFileName.c_str(), "read")};
  if (!fi){
    std::cerr << "\nERROR: cannot read input file with trigger efficiency!" << std::endl;
    return;
  }
  std::unique_ptr<TH1D> hTrEff{(TH1D*)fiTrEff->Get(inTrigHistName.c_str())};
  if (!hTrEff){
    std::cerr << "\nERROR: cannot read input histogram with trigger efficiency!" << std::endl;
    return;
  }

  // Init main process
  ToyMc mc;
  mc.SetParameters(f, k, mu, p);
  mc.SetInput(std::move(tree));
  mc.SetNancestors(Na);
  mc.SetOutput(outFileName.c_str());
  mc.SetNevents(Nev);
  mc.SetTriggerEfficiency(*hTrEff);
  mc.UseGamma();

  mc.Print();
  mc.Run();

  mc.Write();
}
