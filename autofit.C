{
    #include <string>
    #include <cmath>
    #include <fstream>

    TCanvas *c1 = nullptr;
    TString fileName;
    fstream ofile;
    Float_t particleMass, x1, x2;
    Int_t n;  //hist or graph



    cout << endl << "Choose (1 or 2): " << endl << "1. TH1 (hist)" << endl << "2. TGraphErrors (graph)" << endl;
    cin >> n;
    cout << "Enter the name of .root file: "; 
    cin >> fileName;

    TFile *tfile = TFile::Open(fileName + ".root");
    if(!tfile || tfile->IsZombie())
    {
        return 1;
    }
    TIter next(tfile->GetListOfKeys());
    TKey *key;


    //----------FITTING FUNCTIONS----------
    TF1 *fitFunc[5][2];  //5 lines, 2 columns
    TF1 *meanFunc[5][2], *meanFunc2[5][2];
    
    fitFunc[0][0] = new TF1("fitFunc_11", "[0]/6.283/[1]/[1]*TMath::Exp(-x/[1])", 0, 1.5);
    fitFunc[0][1] = new TF1("fitFunc_12", "[0]/6.283/[1]/[1]*TMath::Exp(-x/[1])", 0, 2);
    fitFunc[1][0] = new TF1("fitFunc_21", "[0]*TMath::Exp(-(TMath::Sqrt(x*x+[1]*[1]))/[2])", 0, 1.5);
    fitFunc[1][1] = new TF1("fitFunc_22", "[0]*TMath::Exp(-(TMath::Sqrt(x*x+[1]*[1]))/[2])", 0, 2);
    fitFunc[2][0] = new TF1("fitFunc_31", "[0]*TMath::Exp(-(TMath::Sqrt(x*x+[1]*[1])-[1])/[2])", 0, 1.5);
    fitFunc[2][1] = new TF1("fitFunc_32", "[0]*TMath::Exp(-(TMath::Sqrt(x*x+[1]*[1])-[1])/[2])", 0, 2);
    fitFunc[3][0] = new TF1("fitFunc_41", "[0]*TMath::Sqrt(x*x+[1]*[1])*TMath::Exp(-TMath::Sqrt(x*x+[1]*[1])/[2])", 0, 1.5);
    fitFunc[3][1] = new TF1("fitFunc_42", "[0]*TMath::Sqrt(x*x+[1]*[1])*TMath::Exp(-TMath::Sqrt(x*x+[1]*[1])/[2])", 0, 2);
    fitFunc[4][0] = new TF1("fitFunc_51", "[0]*TMath::Exp(-x*x/[2]/[2])+[3]*TMath::Exp(-x*x/[4]/[4])", 0, 1.5);
    fitFunc[4][1] = new TF1("fitFunc_52", "[0]*TMath::Exp(-x*x/[2]/[2])+[3]*TMath::Exp(-x*x/[4]/[4])", 0, 2);



    vector<TH1F*> histArray;
    if(n==1)  //histograms
    {
        while((key = (TKey*)next()))
        {
            TObject *obj = key->ReadObj();
            TH1 *hist = dynamic_cast<TH1*>(obj);
            histArray.push_back(dynamic_cast<TH1F*>(hist));  //filling histArray vector
        }
    }


    vector<TGraphErrors*> graphArray;
    if(n==2)  //graphs with errors
    {
        while((key = (TKey*)next()))
        {
            if(key->IsFolder())
            {
                TDirectory *dir = (TDirectory*)key->ReadObj();
                dir->cd();

                TList *subKeyList = dir->GetListOfKeys();
                TIter subNext(subKeyList);
                TKey *subKey;

                while((subKey = (TKey*)subNext()))
                {
                    TGraphErrors *graph = (TGraphErrors*)subKey->ReadObj();
                    graphArray.push_back(graph);
                }
            }
        }
    }





  
    //----------DATA WRITING----------
    ofile.open("output_data.txt", ios::out);

    if(n==1)  //histograms
    {
        for(size_t i=0; i<histArray.size(); ++i)
        {
            cout << endl << "Enter the mass for '" << histArray[i]->GetName() << "' : "; cin >> particleMass;
            for(Int_t j=0; j<=4; j++)  //loop over all fit functions
            {
                if(j==0 || j==4) continue;  //fara prima si ultima functie de fit
                for(Int_t k=0; k<=1; k++)  //loop over the two fitting intervals (0,1.5) and (0,2)
                {
                    x1=0;
                    if(k==0) x2=1.5;
                    if(k==1) x2=2.0;

                    fitFunc[j][k]->FixParameter(1, particleMass);
                    fitFunc[j][k]->SetParameter(0, 200);
                    fitFunc[j][k]->SetParameter(2, 0.2);

                    histArray[i]->Fit(fitFunc[j][k], "Q0", "", x1, x2);

                    meanFunc[j][k] = new TF1("mean", ("6.283 * x * fitFunc_" + std::to_string(j+1) + std::to_string(k+1)).c_str(), 0, 10);
                    meanFunc2[j][k] = new TF1("mean2", ("6.283 * x * x * fitFunc_" + std::to_string(j+1) + std::to_string(k+1)).c_str(), 0, 10);

                    Float_t meanPt = meanFunc2[j][k]->Integral(0, 10) / meanFunc[j][k]->Integral(0, 10);
                    Float_t meanPtErr = (meanFunc2[j][k]->Integral(0, 10) / meanFunc[j][k]->Integral(0, 10)) * TMath::Sqrt((meanFunc2[j][k]->IntegralError(0, 10) * meanFunc2[j][k]->IntegralError(0, 10)) / (meanFunc2[j][k]->Integral(0, 10) * meanFunc2[j][k]->Integral(0, 10)) + (meanFunc[j][k]->IntegralError(0, 10) * meanFunc[j][k]->IntegralError(0, 10)) / (meanFunc[j][k]->Integral(0, 10) * meanFunc[j][k]->Integral(0, 10)));


                    ofile << histArray[i]->GetName() << "   " << ("fit_function" + std::to_string(j+1)).c_str() << "   " << "(" << x1 << "," << x2 << ")" << "   " << fitFunc[j][k]->GetParameter(0) << "+/-" << fitFunc[j][k]->GetParError(0) << "   " << fitFunc[j][k]->GetParameter(1) << "   " << fitFunc[j][k]->GetParameter(2) << "+/-" << fitFunc[j][k]->GetParError(2) << "   " << meanPt << "+/-" << meanPtErr << "   " << fitFunc[j][k]->GetChisquare() << "/" << fitFunc[j][k]->GetNDF() << endl;
                }
            }
        }
    }


    if(n==2)  //graphs with errors
    {
        for(size_t i=0; i<graphArray.size(); ++i)
        {
            cout << endl << "Enter the mass for '" << graphArray[i]->GetName() << "' : "; cin >> particleMass;
            for(Int_t j=0; j<=4; j++)  //loop over all fit functions
            {
                if(j==0 || j==4) continue;  //fara prima si ultima functie de fit
                for(Int_t k=0; k<=1; k++)  //loop over the two fitting intervals (0,1.5) and (0,2)
                {
                    x1=0;
                    if(k==0) x2=1.5;
                    if(k==1) x2=2.0;

                    fitFunc[j][k]->FixParameter(1, particleMass);
                    fitFunc[j][k]->SetParameter(0, 200);
                    fitFunc[j][k]->SetParameter(2, 0.2);

                    graphArray[i]->Fit(fitFunc[j][k], "Q0", "", x1, x2);

                    meanFunc[j][k] = new TF1("mean", ("6.283 * x * fitFunc_" + std::to_string(j+1) + std::to_string(k+1)).c_str(), 0, 10);
                    meanFunc2[j][k] = new TF1("mean2", ("6.283 * x * x * fitFunc_" + std::to_string(j+1) + std::to_string(k+1)).c_str(), 0, 10);

                    Float_t meanPt = meanFunc2[j][k]->Integral(0, 10) / meanFunc[j][k]->Integral(0, 10);
                    Float_t meanPtErr = (meanFunc2[j][k]->Integral(0, 10) / meanFunc[j][k]->Integral(0, 10)) * TMath::Sqrt((meanFunc2[j][k]->IntegralError(0, 10) * meanFunc2[j][k]->IntegralError(0, 10)) / (meanFunc2[j][k]->Integral(0, 10) * meanFunc2[j][k]->Integral(0, 10)) + (meanFunc[j][k]->IntegralError(0, 10) * meanFunc[j][k]->IntegralError(0, 10)) / (meanFunc[j][k]->Integral(0, 10) * meanFunc[j][k]->Integral(0, 10)));


                    ofile << graphArray[i]->GetName() << "   " << ("fit_function" + std::to_string(j+1)).c_str() << "   " << "(" << x1 << "," << x2 << ")" << "   " << fitFunc[j][k]->GetParameter(0) << "+/-" << fitFunc[j][k]->GetParError(0) << "   " << fitFunc[j][k]->GetParameter(1) << "   " << fitFunc[j][k]->GetParameter(2) << "+/-" << fitFunc[j][k]->GetParError(2) << "   " << meanPt << "+/-" << meanPtErr << "   " << fitFunc[j][k]->GetChisquare() << "/" << fitFunc[j][k]->GetNDF() << endl;
                }
            }
        }
    }


    ofile.close();
    tfile->Close();
}
