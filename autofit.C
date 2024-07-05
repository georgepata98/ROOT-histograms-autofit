{
  #include <string>

  TString fileName;
  fstream ofile;
  Float_t particleMass, x1, x2, meanPt, meanPtErr;


  cout << "Enter the name of .root file: "; cin >> fileName;

  TFile *tfile = TFile::Open(fileName + ".root");
  if(!tfile || tfile->IsZombie())
  {
    return 1;
  }

  TCanvas *c1 = nullptr;
  TIter next(tfile->GetListOfKeys());
  TKey *key;

  //fitting functions
  TF1 *fitFunc[5][2];  //
  TF1 *meanFunc[5][2], *meanFunc2[5][2];

  fitFunc[0][0] = new TF1("fitFunc_11", "[0]/6.283/[1]/[1]*exp(-x/[1])", 0, 1.5);
  fitFunc[0][1] = new TF1("fitFunc_12", "[0]/6.283/[1]/[1]*exp(-x/[1])", 0, 2);
  fitFunc[1][0] = new TF1("fitFunc_21", "[0]*exp(-1*(sqrt(x*x+[1]*[1]))/[2])", 0, 1.5);
  fitFunc[1][1] = new TF1("fitFunc_22", "[0]*exp(-1*(sqrt(x*x+[1]*[1]))/[2])", 0, 2);
  fitFunc[2][0] = new TF1("fitFunc_31", "[0]*exp(-1*(sqrt(x*x+[1]*[1])-[1])/[2])", 0, 1.5);
  fitFunc[2][1] = new TF1("fitFunc_32", "[0]*exp(-1*(sqrt(x*x+[1]*[1])-[1])/[2])", 0, 2);
  fitFunc[3][0] = new TF1("fitFunc_41", "[0]*sqrt(x*x+[1]*[1])*exp(-sqrt(x*x+[1]*[1])/[2])", 0, 1.5);
  fitFunc[3][1] = new TF1("fitFunc_42", "[0]*sqrt(x*x+[1]*[1])*exp(-sqrt(x*x+[1]*[1])/[2])", 0, 2);
  fitFunc[4][0] = new TF1("fitFunc_51", "[0]*exp(-x*x/[2]/[2])+[3]*exp(-x*x/[4]/[4])", 0, 1.5);
  fitFunc[4][1] = new TF1("fitFunc_52", "[0]*exp(-x*x/[2]/[2])+[3]*exp(-x*x/[4]/[4])", 0, 2);



  vector<TH1F*> histArray;
  vector<TGraphErrors*> graphArray;

  TGraphErrors *graph;

  while((key = (TKey*)next()))
  {
    TObject *obj = key->ReadObj();
    TH1 *hist = dynamic_cast<TH1*>(obj);
    histArray.push_back(dynamic_cast<TH1F*>(hist));  //filling histArray vector
  }

  for(size_t i=0; i<histArray.size(); ++i)  //filling graphArray with histArray
  {
    graph = new TGraphErrors(histArray[i]->GetNbinsX());
    for(Int_t j=1; j<=histArray[i]->GetNbinsX(); j++)
    {
      Float_t x = histArray[i]->GetBinCenter(j);
      Float_t y = histArray[i]->GetBinContent(j);
      Float_t ex = histArray[i]->GetBinWidth(j)/2;
      Float_t ey = histArray[i]->GetBinError(j);

      graph->SetPoint(j-1, x, y);
      graph->SetPointError(j-1, ex, ey);
    }
    graphArray.push_back(graph);
    delete graph;
  }




  
  //Data writing
  ofile.open("output_data.txt", ios::out);  //write data to .txt file

  for(size_t i=0; i<histArray.size(); ++i)
  {
    cout << endl << "Enter the mass for '" << histArray[i]->GetName() << "' : "; cin >> particleMass;
    for(Int_t j=0; j<=4; j++)  //loop over all fit functions
    {
      if(j==0 || j==4) continue;
      for(Int_t k=0; k<=1; k++)  //loop over the two fitting intervals: (0,1.5) and (0,2)
      {
        x1=0;
        if(k==0) x2=1.5;
        if(k==1) x2=2.0;

        fitFunc[j][k]->FixParameter(1, particleMass);
        fitFunc[j][k]->SetParameter(0, 200);
        fitFunc[j][k]->SetParameter(2, 0.2);

        histArray[i]->Fit(fitFunc[j][k], "Q0", "", x1, x2);
        //graphArray[i]->Fit(fitFunc[j][k], "Q0", "", x1, x2);

        meanFunc[j][k] = new TF1("mean", ("6.283 * x * fitFunc_" + std::to_string(j+1) + std::to_string(k+1)).c_str(), 0, 10);
        meanFunc2[j][k] = new TF1("mean2", ("6.283 * x * x * fitFunc_" + std::to_string(j+1) + std::to_string(k+1)).c_str(), 0, 10);

        meanPt = meanFunc2[j][k]->Integral(0, 10) / meanFunc[j][k]->Integral(0, 10);
        meanPtErr = (meanFunc2[j][k]->Integral(0, 10) / meanFunc[j][k]->Integral(0, 10)) * sqrt(( meanFunc2[j][k]->IntegralError(0, 10) * meanFunc2[j][k]->IntegralError(0, 10)) / (meanFunc2[j][k]->Integral(0, 10) * meanFunc2[j][k]->Integral(0, 10)) + (meanFunc[j][k]->IntegralError(0, 10) * meanFunc[j][k]->IntegralError(0, 10)) / (meanFunc[j][k]->Integral(0, 10) * meanFunc[j][k]->Integral(0, 10)));


        ofile << histArray[i]->GetName() << "   " << ("fit_function" + std::to_string(j+1)).c_str() << "   " << "(" << x1 << "," << x2 << ")" << "   " << fitFunc[j][k]->GetParameter(0) << "+/-" << fitFunc[j][k]->GetParError(0) << "   " << fitFunc[j][k]->GetParameter(1) << "   " << fitFunc[j][k]->GetParameter(2) << "+/-" << fitFunc[j][k]->GetParError(2) << "   " << meanPt << "+/-" << meanPtErr << "   " << fitFunc[j][k]->GetChisquare() << "/" << fitFunc[j][k]->GetNDF() << endl;
      }
    }
  }

  ofile.close();
  tfile->Close();
}