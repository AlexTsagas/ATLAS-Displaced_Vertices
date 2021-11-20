// Readers to access the data (delete the ones you do not need).
TTreeReader treereader;

TTreeReaderValue<ULong_t> eventNumber = {treereader, "eventNumber"};
TTreeReaderValue<Int_t> track_n = {treereader, "track_n"};
TTreeReaderValue<Int_t> truthvtx_n = {treereader, "truthvtx_n"};
TTreeReaderArray<Double_t> track_x0 = {treereader, "track.x0"};
TTreeReaderArray<Double_t> track_y0 = {treereader, "track.y0"};
TTreeReaderArray<Double_t> track_z0 = {treereader, "track.z0"};
TTreeReaderArray<Double_t> track_x1 = {treereader, "track.x1"};
TTreeReaderArray<Double_t> track_y1 = {treereader, "track.y1"};
TTreeReaderArray<Double_t> track_z1 = {treereader, "track.z1"};
// TTreeReaderArray<Double_t> emcluster_x0 = {treereader, "emcluster.x0"};
// TTreeReaderArray<Double_t> emcluster_y0 = {treereader, "emcluster.y0"};
// TTreeReaderArray<Double_t> emcluster_z0 = {treereader, "emcluster.z0"};
// TTreeReaderArray<Double_t> emcluster_x1 = {treereader, "emcluster.x1"};
// TTreeReaderArray<Double_t> emcluster_y1 = {treereader, "emcluster.y1"};
// TTreeReaderArray<Double_t> emcluster_z1 = {treereader, "emcluster.z1"};
TTreeReaderArray<Double_t> truthvtx_x = {treereader, "truthvtx.x"};
TTreeReaderArray<Double_t> truthvtx_y = {treereader, "truthvtx.y"};
TTreeReaderArray<Double_t> truthvtx_z = {treereader, "truthvtx.z"};

// // // // // // // // // // // // // // // // // // // // // // // // // // 

void myAnalyzeStage1()
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
      cout<<"#Truth Vertices: "<<*runNumber<<endl;
      cout<<endl;
    }
   }

   infile->Close();
}