// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

// Some functions to make my life easier

// Match TVector3 Class Object to Array in C++
double *makeArray(TVector3 tvector3)
{
    static double array[3];

    for(int i=0; i<3; i++)
    {
        array[i] = tvector3[i];
    }

    return array;
}

// Match Array in C++ to a TVector3 Class Object 
TVector3 makeTVector3(double *array)
{
    TVector3 tvector3;

    for(int i=0; i<3; i++)
    {
        tvector3[i] = array[i];
    }

    return tvector3;
}

// Print Components TVector3 Object
void printCoordinates(TVector3 vector)
{
    cout<<"("<<vector[0]<<", "<<vector[1]<<", "<<vector[2]<<")"<<endl;
}

// Define the norm of a vector a as \sqrt{\vec{a}\cdot\vec{a}}
double norm(TVector3 a)
{
    double norm = a.Dot(a);
    norm = pow(norm, 0.5);

    return norm;
}

// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// The main functions needed for the program to run

// Equation of Line in 3D
TVector3 line(TVector3 a, TVector3 b, double t)
{
    TVector3 N = b - a; 
    TVector3 n = N.Unit(); // Unit Vector in line's direction

    TVector3 x;

    for(int i=0; i<3; i++)
    {
        x[i] = a[i] + t*n[i];
    }

    return x;
}


// Function that returns values of which the minimun is the distance between line1 and line2
// Line1 contects a, b and line2 conects aa, bb vectors
// The idea is: whe find the distance from "every" point of line2 to line 1. While x[0] varies
// the point of which we compute the distance runs on line2.
double distancefunc(double *x, double *par)
{
    // These are the points of the two lines respectively
    TVector3 a(par[0], par[1], par[2]);
    TVector3 b(par[3], par[4], par[5]);
    TVector3 aa(par[6], par[7], par[8]);
    TVector3 bb(par[9], par[10], par[11]);

    // Calculate the unit vector in line1's direction
    TVector3 N = b - a;
    TVector3 n = N.Unit();

    // Equation of line2
    TVector3 p = line(aa, bb, x[0]);

    // The vector that its norm we need
    TVector3 d_vec;
    d_vec = (a - p).Cross(n);

    // Calculate the norm. The function d returns values of that represent the distance of
    // a point in line2 from line1.
    double d = d_vec.Dot(d_vec);
    d = pow(d , 0.5);

    return d;
}


// There we minimize the above function to compute the "distance" between line1 and line2
// The function puts the minimized value and the argument of that in the array x[2]
double *calculateDistance(TVector3 a, TVector3 b, TVector3 aa, TVector3 bb)
{
    // Parameters' Array
    double par[12];
    par[0] = a.X(); par[1] = a.Y(); par[2] = a.Z();
    par[3] = b.X(); par[4] = b.Y(); par[5] = b.Z();
    par[6] = aa.X(); par[7] = aa.Y(); par[8] = aa.Z();
    par[9] = bb.X(); par[10] = bb.Y(); par[11] = bb.Z();

    TF1 *func = new TF1("func", distancefunc, -1e2, 1e2, 12);

    for(int i=0; i<12; i++)
    {
        func->SetParameter(i, par[i]);
    }

    static double x[2];

    //The distance between line1 and line2
    double distance = func->GetMinimum(1e2, 1e2);
    x[0] = distance;

    // the parameter which gives the point of line2 from which the distance to line1 is minimum 
    double t_distance = func->GetMinimumX(1e2, 1e2); 
    x[1] = t_distance;

    return x;
}


// The function computes the point of the displaced vertex, where the long-lived particle decayed
TVector3 displacedVertex(TVector3 a, TVector3 b, TVector3 aa, TVector3 bb)
{
    // Compute the unit vector in the direction of line1
    TVector3 N = b - a;
    TVector3 n = N.Unit();

    // Compute the unit vector in the direction of line2
    TVector3 NN = bb - aa;
    TVector3 nn = NN.Unit();

    // Argument of point in line2
    double *x;
    x = calculateDistance(a, b, aa, bb);
    double t2 = x[1];

    // Argument of point in line1
    double t1 = t2 * nn.Dot(n) + (aa - a).Dot(n);

    // Point in line1 
    TVector3 r = line(a, b, t1);

    // Point in line2
    TVector3 rr = line(aa, bb, t2);

    // Displaced Vertex 
    TVector3 R = 0.5 * (r + rr);

    return R;
}


// Computes the minimum value of array distance[elementCount]
double minimumArrayValue(double *distance, int elementCount)
{
  double minimum = 0;

  for(int i = 0; i < elementCount; i++)
  {
    double element = distance[i];

    if(element < minimum)
    {
      minimum = element;
    }
  }

  return minimum;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

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


   while (treereader.Next()) 
   {
    if(*truthvtx_n==1)
    {
      TVector3 a(track_x0[0], track_y0[0], track_z0[0]);
      TVector3 b(track_x1[0], track_y1[0], track_z1[0]);

      TVector3 aa(track_x0[1], track_y0[1], track_z0[1]);
      TVector3 bb(track_x1[1], track_y1[1], track_z1[1]);

      double *distance;
      distance = calculateDistance(a, b, aa, bb);

      cout<<"The distance between line1 and line2 is: "<<distance[0]<<endl;

      TVector3 R = displacedVertex(a, b, aa, bb);

      cout<<"Displaced Vertex = "<<"("<<R.X()<<", "<<R.Y()<<", "<<R.Z()<<")"<<endl;
    }
   }

   infile->Close();
}