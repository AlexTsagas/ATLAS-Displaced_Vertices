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


// Computes the minimum value of array distance[elementCount][3]
double *minimumArrayValue(double distance[][3], int elementCount)
{
    static double minimum[3];
    minimum[0] = distance[0][0];
    minimum[1] = distance[0][1];
    minimum[2] = distance[0][2];

    double element;

    for(int i=1; i<elementCount; i++)
    {
        element = distance[i][0];

        if(minimum[0] > element)
        {
            minimum[0] = element;
            minimum[1] = distance[i][1];
            minimum[2] = distance[i][2];
        }
    }

    return minimum;
}


// Calculates the distance between Point1(x1, y1, z1) and (x2, y2, z2)
double Error(double x1, double y1, double z1, double x2, double y2, double z2)
{
    double error;

    error = sqrt(pow(x2-x1,2)+pow(y2-y1,2)+pow(z2-z1,2));

    return error;
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
TTreeReaderArray<Double_t> truthvtx_x = {treereader, "truthvtx.x"};
TTreeReaderArray<Double_t> truthvtx_y = {treereader, "truthvtx.y"};
TTreeReaderArray<Double_t> truthvtx_z = {treereader, "truthvtx.z"};

// // // // // // // // // // // // // // // // // // // // // // // // // // 

void myAnalyzeStage1()
{
    TH1D *h1 = new TH1D("h1", "Absolute Error;Error;Counts", 70, 0, 60);
    TH1D *h2 = new TH1D("hist", "Minimum Line Distance to Each Event;Distance;Counts", 40, 0, 0.6);

    TFile* infile = TFile::Open("stage1.root");
    TTree* tree   = (TTree*)infile->Get("stage1");
    treereader.SetTree(tree);

    // Line1 Points
    TVector3 a, b;

    // Line2 Points
    TVector3 aa, bb;

    // Array that gathers the distances between line i and j and stores them in the first column.
    //  The second and the third column store information about line_i and line_j, respectively
    double distance_ij[100][3];

    // Array that gathers the minimum values of distance_ij array for iterations to all lines_i and stores
    //  them in the first column. The second and the third columns store information about line_i and line_j,
    //  respectively
    double minDistanceArray[10][3];
    double *min_Distance_Array;

    // The minimum distance of all the possible pair of lines for every event (total events = 4299)
    // The first column contains the least distance for every event. The other two to which two lines
    // it is refering to.
    double leastDistance[4299][3];
    double *least_Distance;

    // Stores the coordinates of each displaced vertex of any event in rows
    double displacedVertexArray[4299][3];
    TVector3 displaced_Vertex;
    // Tvector3 needed fro the function displacedVertex
    TVector3 A, B, AA, BB;
    // Counters
    int I, J;

    // The distance of the truth displaced vertex from the calculated displaced vertex.
    // The m^th element is the calculated error of m^th event
    double absoluteError[4299];

    // Pointer to ascribe the distance between line1 and line2
    double *d;

    // Integers for for loops
    int i, j;

    // Counters for arrays' elemets
    int count_i, count_j;

    // Event Counter
    int event = 0;

    // Elements Counter
    int elementCount;

    while (treereader.Next()) 
    {
        if(*truthvtx_n==1) 
        {   
            count_i = 0;
            count_j = 0;

            for(i=0; i<*track_n; i++)
            {
                a.SetX(track_x0[i]); a.SetY(track_y0[i]); a.SetZ(track_z0[i]);
                b.SetX(track_x1[i]); b.SetY(track_y1[i]); b.SetZ(track_z1[i]);

                for(j=i+1; j<*track_n; j++)
                {   
                    aa.SetX(track_x0[j]); aa.SetY(track_y0[j]); aa.SetZ(track_z0[j]);
                    bb.SetX(track_x1[j]); bb.SetY(track_y1[j]); bb.SetZ(track_z1[j]);

                    d = calculateDistance(a, b, aa, bb);

                    distance_ij[count_j][0] = d[0];
                    distance_ij[count_j][1] = i;
                    distance_ij[count_j][2] = j;

                    // // 
                    // cout<<event<<". \t";
                    // cout<<count_j<<"("<<*track_n<<")"<<". ";
                    // for(int k=0; k<3; k++)
                    // {
                    //     cout<<distance_ij[count_j][k]<<" ";
                    // }
                    // cout<<endl;
                    // // 

                    // // 
                    // for(int k=1; k<3; k++)
                    // {
                    //     if(distance_ij[count_j][k]>=*track_n)
                    //     {
                    //         cout<<-1<<", ";
                    //     }
                    // }
                    // // 

                    count_j++;
                }

                // count_j -= 1; //so as count_j = elements in distance_ij rows

                // min_Distance_Array = minimumArrayValue(distance_ij, count_j);
                // minDistanceArray[count_i][0] = min_Distance_Array[0];
                // minDistanceArray[count_i][1] = min_Distance_Array[1];
                // minDistanceArray[count_i][2] = min_Distance_Array[2];

                count_i++;
            }

            elementCount = count_j;

            least_Distance = minimumArrayValue(distance_ij, elementCount);

            leastDistance[event][0] = least_Distance[0];
            leastDistance[event][1] = least_Distance[1];
            leastDistance[event][2] = least_Distance[2];

            h2->Fill(leastDistance[event][0]);

            // // 
            // cout<<event<<". ";
            // for(int k=0; k<3; k++)
            // {
            //     cout<<leastDistance[event][k]<<" ";
            // }
            // cout<<endl;
            // // 

            // //  3335e6897
            // for(int k=1; k<3; k++)
            // {
            //     if(leastDistance[event][k]>=*track_n)
            //     {
            //         //  leastDistance[event][k]=0;
            //         cout<<event<<". "<<-*track_n+leastDistance[event][k]<<endl;
            //     }
            // }
            // // 

            // count_i -= 1; //so as count_i = elements in minDistanceArray rows

            // least_Distance = minimumArrayValue(minDistanceArray, count_i); 
            // leastDistance[event][0] = least_Distance[0];
            // leastDistance[event][1] = least_Distance[1];
            // leastDistance[event][2] = least_Distance[2];

            I = leastDistance[event][1];
            J = leastDistance[event][2];

            // cout<<track_x0[0]<<", "<<track_y0[0]<<","<<track_z0[0]<<endl;
            // cout<<track_x0[J]<<", "<<track_y0[J]<<","<<track_z0[J]<<endl;

            A.SetX(track_x0[I]); A.SetY(track_y0[I]); A.SetZ(track_z0[I]);
            B.SetX(track_x1[I]); B.SetY(track_y1[I]); B.SetZ(track_z1[I]);
            AA.SetX(track_x0[J]); AA.SetY(track_y0[J]); AA.SetZ(track_z0[J]);
            BB.SetX(track_x1[J]); BB.SetY(track_y1[J]); BB.SetZ(track_z1[J]);

            displaced_Vertex = displacedVertex(A, B, AA, BB);
            displacedVertexArray[event][0] = displaced_Vertex.X();
            displacedVertexArray[event][1] = displaced_Vertex.Y();
            displacedVertexArray[event][2] = displaced_Vertex.Z();

            event++;

            absoluteError[event] = Error(displacedVertexArray[event][0], displacedVertexArray[event][1], displacedVertexArray[event][2], truthvtx_x[0], truthvtx_y[0], truthvtx_z[0]);

            h1->Fill(absoluteError[event]);

        }
    }

    h1->SetLineColor(kBlack);
    h1->SetLineWidth(2);
    // h1->Draw();

    h2->SetLineColor(kBlue);
    h2->SetLineWidth(2);
    // h2->Draw();

    cout<<endl<<"The program Works Fine!"<<endl<<endl;

    infile->Close();
    }
