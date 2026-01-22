void analyze() {
    TFile *f = new TFile("output0.root");
    
    TTree *t = (TTree*)f->Get("DetectorHits");
    
    TH1F *h1 = new TH1F("h1", "Hits per Detector", 3, 0, 3); // 3 bins: 0, 1, 2
    
    t->Draw("DetectorID>>h1");
    
    int bottom = t->GetEntries("DetectorID==0");
    int middle = t->GetEntries("DetectorID==1");
    int top    = t->GetEntries("DetectorID==2");
    
    cout << "-----------------------------" << endl;
    cout << "Basement Hits: " << bottom << endl;
    cout << "Lab Hits: " << middle << endl;
    cout << "Top Hits:    " << top << endl;
    cout << "-----------------------------" << endl;
}
