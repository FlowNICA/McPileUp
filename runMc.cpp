#include <iostream>
#include<string>

#include <TROOT.h>
#include <TStopwatch.h>
#include <Math/SpecFuncMathMore.h>
#include <TRandom.h>
#include <TRandom3.h>

#include "ToyMc.h"

int main(int argc, char **argv)
{
  Int_t Nev, fK;
  TString inName, treeName, outName;
  Double_t fF, fMu, fP;
  int seed = -1;

  TStopwatch timer;
  if (argc < 17)
  {
    std::cerr << "./runGlauber -nev Nev -i INPUT -iname INPUT_TREE_NAME -f F_PARAMETER -k K_PARAMETER -mu MU_PARAMETER -p PILEUP_RATE -seed SEED -o OUTPUTFILE" << std::endl;
    return 10;
  }
  for (int i=1;i<argc;i++)
  {
    if (std::string(argv[i]) != "-nev" &&
        std::string(argv[i]) != "-o" &&
        std::string(argv[i]) != "-i" &&
        std::string(argv[i]) != "-iname" &&
        std::string(argv[i]) != "-f" &&
        std::string(argv[i]) != "-k" &&
        std::string(argv[i]) != "-mu" &&
        std::string(argv[i]) != "-seed" &&
        std::string(argv[i]) != "-p")
    {
      std::cerr << "\nUnknown parameter: " << argv[i] << std::endl;
      return 11;
    }
    else
    {
      if (std::string(argv[i]) == "-i" && i != argc-1)
      {
        inName = TString(argv[++i]);
      }
      if (std::string(argv[i]) == "-i" && i == argc-1)
      {
        std::cerr << "\nInput file is not defined!" << std::endl;
        return 12;
      }
      if (std::string(argv[i]) == "-o" && i != argc-1)
      {
        outName = TString(argv[++i]);
      }
      if (std::string(argv[i]) == "-o" && i == argc-1)
      {
        std::cerr << "\nOutput file is not defined!" << std::endl;
        return 13;
      }
      if (std::string(argv[i]) == "-f" && i != argc-1)
      {
        fF = std::stod(std::string(argv[++i]));
      }
      if (std::string(argv[i]) == "-f" && i == argc-1)
      {
        std::cerr << "\nParameter f is not defined!" << std::endl;
        return 14;
      }
      if (std::string(argv[i]) == "-p" && i != argc-1)
      {
        fP = std::stod(std::string(argv[++i]));
      }
      if (std::string(argv[i]) == "-p" && i == argc-1)
      {
        std::cerr << "\nParameter p is not defined!" << std::endl;
        return 15;
      }
      if (std::string(argv[i]) == "-mu" && i != argc-1)
      {
        fMu = std::stod(std::string(argv[++i]));
      }
      if (std::string(argv[i]) == "-mu" && i == argc-1)
      {
        std::cerr << "\nParameter mu is not defined!" << std::endl;
        return 16;
      }
      if (std::string(argv[i]) == "-k" && i != argc-1)
      {
        fK = std::stoi(std::string(argv[++i]));
      }
      if (std::string(argv[i]) == "-k" && i == argc-1)
      {
        std::cerr << "\nParameter k is not defined!" << std::endl;
        return 17;
      }
      if (std::string(argv[i]) == "-nev" && i != argc-1)
      {
        Nev = std::stoi(std::string(argv[++i]));
      }
      if (std::string(argv[i]) == "-nev" && i == argc-1)
      {
        std::cerr << "\nNumber of events is not defined!" << std::endl;
        return 18;
      }
      if (std::string(argv[i]) == "-seed" && i != argc-1)
      {
        seed = std::stoi(std::string(argv[++i]));
        if (seed <=0)
        {
          std::cerr << "\nWARNING: wrong value for seed! Set to default." << std::endl;
        }
      }
      if (std::string(argv[i]) == "-seed" && i == argc-1)
      {
        std::cerr << "\nSeed is not defined!" << std::endl;
        return 19;
      }
      if (std::string(argv[i]) == "-iname" && i != argc-1)
      {
        treeName = TString(argv[++i]);
      }
      if (std::string(argv[i]) == "-iname" && i == argc-1)
      {
        std::cerr << "\nInput file is not defined!" << std::endl;
        return 20;
      }
    }
  }
  if (outName == "" || inName == "" || treeName == "")
  {
    std::cerr << "\nOutput/Input/Input tree name has not been set properly!" << std::endl;
    return 100;
  }

  if (seed > 0)
    gRandom->SetSeed(seed);

  timer.Start();

  std::unique_ptr<TFile> fi{new TFile(inName, "read")};
  if (!fi){
    std::cerr << "\nERROR: cannot read input file!" << std::endl;
    return 101;
  }
  std::unique_ptr<TTree> tree{(TTree*)fi->Get(treeName)};
  if (!tree){
    std::cerr << "\nERROR: cannot read input tree!" << std::endl;
    return 102;
  }
  
  ToyMc mc;
  mc.SetParameters(fF, fK, fMu, fP);
  mc.SetInput(std::move(tree));
  mc.SetOutput(outName);
  mc.SetNevents(Nev);

  mc.Print();
  mc.Run();

  mc.Write();

  timer.Stop();
  timer.Print();

  return 0;
}
