void analyze2() {
    const int botA_ID = 10;
    const int botB_ID = 11;

    const int midA_ID = 20;
    const int midB_ID = 21;

    const int topA_ID = 30;
    const int topB_ID = 31;

    TFile *f = new TFile("output0.root");
    if (!f || f->IsZombie()) {
        cout << "Error: Could not open output0.root" << endl;
        return;
    }

    // 2. Get the Tree
    TTree *t = (TTree*)f->Get("DetectorHits");
    if (!t) {
        cout << "Error: Tree 'DetectorHits' not found!" << endl;
        return;
    }

    // 3. Set Branch Addresses to read data
    int eventID, detectorID;
    t->SetBranchAddress("EventID", &eventID);
    t->SetBranchAddress("DetectorID", &detectorID);

    // 4. Variables for Loop
    Long64_t nEntries = t->GetEntries();
    
    // Counters for valid "Coincidence" hits
    int countsBottom = 0;
    int countsMiddle = 0;
    int countsTop = 0;

    // Flags to track hits within the CURRENT event
    bool hitBotA = false; 
    bool hitBotB = false;
    bool hitMidA = false; 
    bool hitMidB = false;
    bool hitTopA = false; 
    bool hitTopB = false;

    int currentEventID = -1;

    // 5. Loop over all entries
    for (Long64_t i = 0; i < nEntries; i++) {
        t->GetEntry(i);

        // --- NEW EVENT DETECTED ---
        if (eventID != currentEventID) {
            // A. Analyze the PREVIOUS event (if it wasn't the startup -1)
            if (currentEventID != -1) {
                // Check Coincidence: Did we hit BOTH A and B?
                if (hitBotA && hitBotB) countsBottom++;
                if (hitMidA && hitMidB) countsMiddle++;
                if (hitTopA && hitTopB) countsTop++;
            }

            // B. Reset flags for the NEW event
            hitBotA = false; hitBotB = false;
            hitMidA = false; hitMidB = false;
            hitTopA = false; hitTopB = false;
            
            // C. Update Current Event ID
            currentEventID = eventID;
        }

        // --- RECORD HIT FLAGS FOR CURRENT EVENT ---
        if (detectorID == botA_ID) hitBotA = true;
        if (detectorID == botB_ID) hitBotB = true;

        if (detectorID == midA_ID) hitMidA = true;
        if (detectorID == midB_ID) hitMidB = true;

        if (detectorID == topA_ID) hitTopA = true;
        if (detectorID == topB_ID) hitTopB = true;
    }

    // 6. Don't forget to process the FINAL event after the loop finishes!
    if (hitBotA && hitBotB) countsBottom++;
    if (hitMidA && hitMidB) countsMiddle++;
    if (hitTopA && hitTopB) countsTop++;

    // 7. Print Results
    cout << "===========================================" << endl;
    cout << " ANALYSIS RESULT: COINCIDENCE HITS " << endl;
    cout << "===========================================" << endl;
    cout << "Bottom Telescope Hits: " << countsBottom << endl;
    cout << "Middle Telescope Hits: " << countsMiddle << endl;
    cout << "Top Telescope Hits:    " << countsTop << endl;
    cout << "===========================================" << endl;
}
