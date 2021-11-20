// Readers to access the data (delete the ones you do not need).
TTreeReader treereader;

TTreeReaderValue<UInt_t> runNumber = {treereader, "runNumber"};
TTreeReaderValue<ULong_t> eventNumber = {treereader, "eventNumber"};
TTreeReaderValue<Double_t> met = {treereader, "met"};
TTreeReaderValue<Double_t> met_phi = {treereader, "met_phi"};
TTreeReaderValue<Int_t> track_n = {treereader, "track_n"};
TTreeReaderValue<Int_t> emcluster_n = {treereader, "emcluster_n"};
TTreeReaderValue<Int_t> tile_n = {treereader, "tile_n"};
TTreeReaderValue<Int_t> truthvtx_n = {treereader, "truthvtx_n"};
TTreeReaderArray<unsigned int> track_fUniqueID = {treereader, "track.fUniqueID"};
TTreeReaderArray<unsigned int> track_fBits = {treereader, "track.fBits"};
TTreeReaderArray<TString> track_fName = {treereader, "track.fName"};
TTreeReaderArray<TString> track_fTitle = {treereader, "track.fTitle"};
TTreeReaderArray<Int_t> track_charge = {treereader, "track.charge"};
TTreeReaderArray<Double_t> track_pt = {treereader, "track.pt"};
TTreeReaderArray<Double_t> track_p = {treereader, "track.p"};
TTreeReaderArray<Double_t> track_phi = {treereader, "track.phi"};
TTreeReaderArray<Double_t> track_eta = {treereader, "track.eta"};
TTreeReaderArray<Double_t> track_d0 = {treereader, "track.d0"};
TTreeReaderArray<Double_t> track_calo_iso = {treereader, "track.calo_iso"};
TTreeReaderArray<Double_t> track_track_iso = {treereader, "track.track_iso"};
TTreeReaderArray<Int_t> track_truth_kind = {treereader, "track.truth_kind"};
TTreeReaderArray<Int_t> track_interest = {treereader, "track.interest"};
TTreeReaderArray<Double_t> track_x0 = {treereader, "track.x0"};
TTreeReaderArray<Double_t> track_y0 = {treereader, "track.y0"};
TTreeReaderArray<Double_t> track_z0 = {treereader, "track.z0"};
TTreeReaderArray<Double_t> track_x1 = {treereader, "track.x1"};
TTreeReaderArray<Double_t> track_y1 = {treereader, "track.y1"};
TTreeReaderArray<Double_t> track_z1 = {treereader, "track.z1"};
TTreeReaderArray<Int_t> track_hit_n = {treereader, "track.hit_n"};
TTreeReaderArray<vector<double>> track_hit_eta = {treereader, "track.hit_eta"};
TTreeReaderArray<vector<double>> track_hit_phi = {treereader, "track.hit_phi"};
TTreeReaderArray<vector<double>> track_hit_E = {treereader, "track.hit_E"};
TTreeReaderArray<unsigned int> emcluster_fUniqueID = {treereader, "emcluster.fUniqueID"};
TTreeReaderArray<unsigned int> emcluster_fBits = {treereader, "emcluster.fBits"};
TTreeReaderArray<TString> emcluster_fName = {treereader, "emcluster.fName"};
TTreeReaderArray<TString> emcluster_fTitle = {treereader, "emcluster.fTitle"};
TTreeReaderArray<Int_t> emcluster_charge = {treereader, "emcluster.charge"};
TTreeReaderArray<Double_t> emcluster_pt = {treereader, "emcluster.pt"};
TTreeReaderArray<Double_t> emcluster_p = {treereader, "emcluster.p"};
TTreeReaderArray<Double_t> emcluster_phi = {treereader, "emcluster.phi"};
TTreeReaderArray<Double_t> emcluster_eta = {treereader, "emcluster.eta"};
TTreeReaderArray<Double_t> emcluster_d0 = {treereader, "emcluster.d0"};
TTreeReaderArray<Double_t> emcluster_calo_iso = {treereader, "emcluster.calo_iso"};
TTreeReaderArray<Double_t> emcluster_track_iso = {treereader, "emcluster.track_iso"};
TTreeReaderArray<Int_t> emcluster_truth_kind = {treereader, "emcluster.truth_kind"};
TTreeReaderArray<Int_t> emcluster_interest = {treereader, "emcluster.interest"};
TTreeReaderArray<Double_t> emcluster_x0 = {treereader, "emcluster.x0"};
TTreeReaderArray<Double_t> emcluster_y0 = {treereader, "emcluster.y0"};
TTreeReaderArray<Double_t> emcluster_z0 = {treereader, "emcluster.z0"};
TTreeReaderArray<Double_t> emcluster_x1 = {treereader, "emcluster.x1"};
TTreeReaderArray<Double_t> emcluster_y1 = {treereader, "emcluster.y1"};
TTreeReaderArray<Double_t> emcluster_z1 = {treereader, "emcluster.z1"};
TTreeReaderArray<Int_t> emcluster_hit_n = {treereader, "emcluster.hit_n"};
TTreeReaderArray<vector<double>> emcluster_hit_eta = {treereader, "emcluster.hit_eta"};
TTreeReaderArray<vector<double>> emcluster_hit_phi = {treereader, "emcluster.hit_phi"};
TTreeReaderArray<vector<double>> emcluster_hit_E = {treereader, "emcluster.hit_E"};
TTreeReaderArray<unsigned int> tile_fUniqueID = {treereader, "tile.fUniqueID"};
TTreeReaderArray<unsigned int> tile_fBits = {treereader, "tile.fBits"};
TTreeReaderArray<Double_t> tile_E = {treereader, "tile.E"};
TTreeReaderArray<Double_t> tile_eta = {treereader, "tile.eta"};
TTreeReaderArray<Double_t> tile_phi = {treereader, "tile.phi"};
TTreeReaderArray<unsigned int> truthvtx_fUniqueID = {treereader, "truthvtx.fUniqueID"};
TTreeReaderArray<unsigned int> truthvtx_fBits = {treereader, "truthvtx.fBits"};
TTreeReaderArray<Double_t> truthvtx_x = {treereader, "truthvtx.x"};
TTreeReaderArray<Double_t> truthvtx_y = {treereader, "truthvtx.y"};
TTreeReaderArray<Double_t> truthvtx_z = {treereader, "truthvtx.z"};

//===========================================================
// here's the event processing...
void analyzeStage1()
{
   TFile* infile = TFile::Open("stage1.root");
   TTree* tree   = (TTree*)infile->Get("stage1");
   treereader.SetTree(tree);

   int ievent{0};
   while (treereader.Next()) 
   {
    if(*truthvtx_n==1)
    {
      ++ievent;
      cout<<"Event "<<ievent<<endl;
      cout<<"#Tracks: "<<*track_n<<endl;
      cout<<"#Truth Vertices: "<<*truthvtx_n<<endl;
      cout<<endl;
    }
   }

   infile->Close();
}