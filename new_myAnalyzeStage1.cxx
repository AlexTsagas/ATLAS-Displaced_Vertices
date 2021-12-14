// Functions for Vector Operations //

// Print array[3] elements
void printCoordinates(double *array)
{
    cout<<"("<<array[0]<<", "<<array[1]<<", "<<array[2]<<")"<<endl;
}


// Input two vectors with three coordinates each and output their dot product
double dotProduct(double *a, double *b)
{
    double dotproduct;

    dotproduct = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];

    return dotproduct;
}


// Input two vectors with three coordinates each and output a pointer to the
// first element of the array with three elements of their cross product
double *crossProduct(double *a, double *b)
{
    static double crossproduct[3];

    crossproduct[0] = a[1]*b[2] - a[2]*b[1];
    crossproduct[1] = a[2]*b[0] - a[0]*b[2];
    crossproduct[2] = a[0]*b[1] - a[1]*b[0];

    return crossproduct;
}


// Triple product between vectors a, b, c with three coordinates each. Triple
// product is defined as a * (b x c).
double tripleProduct(double *a, double *b, double *c)
{
    double *crossp = crossProduct(b, c);
    double cross[3] = {crossp[0], crossp[1], crossp[2]};

    double triple = dotProduct(a, cross);

    return triple;
}


// Input an array with three elements and output the norm sqrt(a * a)
double norm(double *a)
{
    double norm;
    double dotp = dotProduct(a, a);

    norm = sqrt(dotp);

    return norm;
}


// Relative vector that points from a to b
double *relativeVector(double *a, double *b)
{
    static double relative[3];
    for(int i=0; i<3; i++)
    {
        relative[i] = b[i] - a[i];
    }

    return relative;
}


// Input a vector a with three elemets. Output the unit vector of a.
double *unitVector(double *a)
{
    static double unitvector[3];

    double Norm = norm(a);

    for(int i=0; i<3; i++)
    {
        unitvector[i] = a[i]/Norm;
    }

    return unitvector;
}


// Compute the coordinates of the unit vector pointing from a to b
double *relativeUnitVector(double *a, double *b)
{
    static double unitRelAB[3];

    double *rel_AB = relativeVector(a, b);
    double relAB[3] = {rel_AB[0], rel_AB[1], rel_AB[2]};

    double *unitRel_AB = unitVector(relAB);
    for(int i=0; i<3; i++)
    {
        unitRelAB[i] = unitRel_AB[i];
    }

    return unitRelAB;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

// Readers to access the data from .root file
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

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //


// Input two points of the line, a and b. Ouput a pointer to the first element of
// the array containing line points in relation with t
double *lineEquation(double *a, double *b, double t)
{   
    // Parallel to line unit vector pointing from a to b
    double *N = relativeVector(a, b);
    double *n = unitVector(N);
   
    // The line point
    static double lineVector[3];
    for(int i=0; i<3; i++)
    {
        lineVector[i] = a[i] + t*n[i];
    }

    return lineVector;
}


// Calculates the distance of line1 which passes through r1, r11 and line2
// which passes through r2, r22.
double distance(double *r1, double *rr1, double *r2, double *rr2)
{
    // Parallel to line1
    double *N1p = relativeVector(r1, rr1);
    double N1[3] = {N1p[0], N1p[1], N1p[2]};
    
    // Parallel to line2
    double *N2p = relativeVector(r2, rr2);
    double N2[3] = {N2p[0], N2p[1], N2p[2]};

    // The point where the plane passes: r0 = r2 - r1
    double *r0p = relativeVector(r1, r2);
    double r0[3] = {r0p[0], r0p[1], r0p[2]};

    // Perpendicular to the plane
    double *u = crossProduct(N2, N1);

    double d = abs(dotProduct(r0, u))/norm(u);

    return d;
}


// Calculates the coordinates of the displaced vertex as the midpoint of AB, where
// |AB| is the distance between the lines.
double *displacedVertex(double *r1, double *rr1, double *r2, double *rr2)
{
    // Parallel to line1
    double *N1p = relativeVector(r1, rr1);
    double N1[3] = {N1p[0], N1p[1], N1p[2]};

    // Parallel to line2
    double *N2p = relativeVector(r2, rr2);
    double N2[3] = {N2p[0], N2p[1], N2p[2]};

    // The point where the plane passes: r0 = r2 - r1
    double *r0p = relativeVector(r1, r2);
    double r0[3] = {r0p[0], r0p[1], r0p[2]};

    // Perpendicular to the plane
    double *up = crossProduct(N2, N1);
    double u[3] = {up[0], up[1], up[2]};

    // Argument to find A
    double t0 = tripleProduct(u, N2, r0)/(norm(u)*norm(u));

    // Line1 point A
    double OA[3];
    for(int i=0; i<3; i++)
    {
        OA[i] = r1[i] + t0*N1[i];
    }

    // Argument to find B
    double s0 = tripleProduct(u, N1, r0)/(norm(u)*norm(u));

    // Line2 point B
    double OB[3];
    for(int i=0; i<3; i++)
    {
        OB[i] = r2[i] + s0*N2[i];
    }

    static double displacedvertex[3];
    for(int i=0; i<3; i++)
    {
        displacedvertex[i] = 0.5*(OA[i] + OB[i]);
    }

    return displacedvertex;
}


// Computes the minimum value of array's first column distance[elementCount][6] and leaves other
// columns' elements intact.
double *minimumArrayValue(double distance[][6], int elementCount)
{
    // First element the minimum value. Second and third the indexes of disranceelementCount3]
    // in the same row.
    static double minimum[6];

    // Initialize with the first row of distance[elementCount][3]
    for(int i=0; i<6; i++)
    {
        minimum[i] = distance[0][i];
    }

    // Go throught the array
    double element;

    for(int i=1; i<elementCount; i++)
    {
        element = distance[i][0];

        if(minimum[0] >= element)
        {
            minimum[0] = element;

            for(int j=1; j<6; j++)
            {
                minimum[j] = distance[i][j];
            }
        }
    }

    return minimum;
}


// Calculates the distance between points: (x1, y1, z1) and (x2, y2, z2).
double Error(double x1, double y1, double z1, double x2, double y2, double z2)
{
    double error;
 
    error = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));

    return error;
}


// Input two points A, B, from where the line (ε) passes, and a point P to
// calculate the distance from (ε) to P.
double LinePointDistance(double *a, double *b, double *p)
{   
    // Vector perpendicular to the line
    double *n_p = relativeVector(a, b);
    double n[3] = {n_p[0], n_p[1], n_p[2]};

    // Vector from A to P
    double *AP_p = relativeVector(a, p);
    double AP[3] = {AP_p[0], AP_p[1], AP_p[2]};

    // Non Normalized Distance
    double *D_p = crossProduct(AP, n);
    double D[3] = {D_p[0], D_p[1], D_p[2]};

    // The actual distance (normalized)
    double distance = norm(D)/norm(n);

    return distance;
}


// Computes the minimum value from array[elementCount]'s elements
double minimumValuefromArrayElements(double *array, int elementCount)
{
    // initialize the minimum
    double minimum = array[0];
    
    double element;

    for(int i=1; i<elementCount; i++)
    {
        element = array[i];

        if(minimum >= element)
        {
            minimum = element;
        }
    }

    return minimum;
}


// Cases to choose if a point is a DV
double CaseDV(double *DV, double *a, double *b, double *aa, double *bb)
{
    // Relative unit vector from Dv to a
    double *unitRel_DVa = relativeUnitVector(DV, a);
    double unitRelDVa[3] = {unitRel_DVa[0], unitRel_DVa[1], unitRel_DVa[2]};

    // Relative unit vector from Dv to b
    double *unitRel_DVb = relativeUnitVector(DV, b);
    double unitRelDVb[3] = {unitRel_DVb[0], unitRel_DVb[1], unitRel_DVb[2]};

    // Relative unit vector from Dv to aa
    double *unitRel_DVaa = relativeUnitVector(DV, aa);
    double unitRelDVaa[3] = {unitRel_DVaa[0], unitRel_DVaa[1], unitRel_DVaa[2]};

    // Relative unit vector from Dv to bb
    double *unitRel_DVbb = relativeUnitVector(DV, bb);
    double unitRelDVbb[3] = {unitRel_DVbb[0], unitRel_DVbb[1], unitRel_DVbb[2]};

    // cos(θ_i) where θ_i is the angle between DVa and DVb vectors
    double dotProd_i = dotProduct(unitRelDVa, unitRelDVb);
    // cos(θ_j) where θ_j is the angle between DVaa and DVbb vectors
    double dotProd_j = dotProduct(unitRelDVaa, unitRelDVbb);

    double prodSum = dotProd_i + dotProd_j;

    return prodSum;
}


// Calculates the distance of the closest third trajectory to the DV
double ThirdTrajectoryDistance(double *DV, int trackNumber, int j, int k)
{
    int elementNumber = 0;
    double a[3], b[3], aa[3], bb[3];
    double DvTrajectoryDistance[20];

    for(int i=0; i<trackNumber; i++)
    {
        if(i!=j && i!=k)
        {
            a[0] = track_x0[i]; a[1] = track_y0[i]; a[2] = track_z0[i];
            b[0] = track_x1[i]; b[1] = track_y1[i]; b[2] = track_z1[i];

            DvTrajectoryDistance[elementNumber] =  LinePointDistance(a, b, DV);

            elementNumber++;
        }
    }

    double DvTrajectory = minimumValuefromArrayElements(DvTrajectoryDistance, elementNumber);

    return DvTrajectory;
}


// Returns false if there is an element of array equal to i and true otherwise
bool TrajectoryUsed(int *array, int i)
{
    int used = 1;

    // 30 elements in array
    for(int k=0; k<30; k++)
    {
        if(i == array[k])
        {
            used = 0;
            break;
        }
    }

    if(used==0)
    {
        return false;
    }
    else
    {
        return true;
    }
}


// ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ //

void new_myAnalyzeStage1()
{
    // For time capture
    clock_t tStart = clock();

    // Histograms
    // 3D //
    TH2 *H = new TH2D("H", "Absolute Error in Relation to (Minimum) Distance of Trajectories;Distance;Absolute Error;Count", 20, 0, 0.15, 30, 0, 6.5);
    TH1 *HH = new TH2D("HH", "Distance R (of DV) from Detector's Center with Respect to Z Coordinate;R;Z;Counts", 40, 0, 40, 40, -40, 40);
    // 2D //
    TH1 *h1 = new TH1D("h1", "Absolute Error;Error;Counts", 50, 0, 6.5);
    TH1 *h2 = new TH1D("h2", "Minimum Trajectory Distance to Each Event;Distance;Counts", 40, 0, 0.15);
    TH1 *h3 = new TH1D("h3", "Distance R of DV From Detector's Center;R;Counts", 40, 0, 40);
    TH1 *h4 = new TH1D("h4", "Z Coordinate of DV;Z;Counts", 40, -40, 40);
    TH1 *h5 = new TH1D("h5", "Smallest Distance D of Third Trajectory from Dv;D;Counts", 50, 0, 0.6);

    TFile* infile = TFile::Open("stage1.root");
    TTree* tree   = (TTree*)infile->Get("stage1");
    treereader.SetTree(tree);


    // Line_i Points
    double a[3], b[3];
    // Line_j Points
    double aa[3], bb[3];

    // Array that gathers the distances between line i and j and stores them in the first column.
    // The second, third and forth column store the coordinates of the displaced vertex resulting from
    // line_i and line_j. The fifth and sixth store the i-th and j-th lines' indexes, respactively.
    double distance_ij[100][6];
    // Elements Counter for distance_ij
    int elementCount;
    // Counters for distance_ij array elemets
    int count_j;


    // The minimum distance of all the possible pair of lines for every event (total events = 4299)
    // The first column contains the least distance for every event. The other three store the coordinates
    // of the displaced vertex resulting from the lines that produce the minimum distance. The fifth and 
    // sixth store the i-th and j-th lines' indexes, respactively.
    double leastDistance[6];
    double *least_Distance;


    // Stores the coordinates of each displaced vertex of any event in rows
    double displacedVertexArray[3];
    double *displaced_Vertex;
    // Coordinates for displaced vertex produced by line_i and line_j
    double DV[3];

    // To apple the case for DVs
    double prodSum;
    // Limits
    double theta = 90;
    double epsilon1 = cos(M_PI * theta/180);
    double epsilon2 = cos(0);


    // The distance between truth DV and computed DV.
    // The m^th element is the calculated error of m^th event
    double absoluteError;

    // The minimun distance of the third trajectory to the DV
    double DvTrajectory;
    
    // Event Counter
    int event = 0;

    // Integers for for loops
    int i, j;

    // For Histograms
    double distance_xyz; // Distance of DV from begining of axis
    double DV_Z; // z coordinate of DV


    // Stores indexes of lines that have been used to calculate a DV
    int usedLineIndex[30];
    // Initialize all values to -1 so as not to coincide with other events
    for(int k=0; k<30; k++)
    {
        usedLineIndex[k] = -1;
    }
    // Counter for usedLineIndex[30] elements
    int countLine;

    // ------------------------------------------------------------------------------------------------------------------------------------------------------ //

    // Event Loop
    while (treereader.Next()) 
    {
        // Loop in events with 1 DV
        if(*truthvtx_n==1)
        {   
            countLine = 0;
            // Initialize all values to -1 so as not to coincide with other events
            for(int k=0; k<30; k++)
            {
                usedLineIndex[k] = -1;
            }

            // Loop to find each DV. In total we have *truthvtx_n DVs
            for(int DvNumber=0; DvNumber<*truthvtx_n; DvNumber++)
            {
                count_j = 0;

                // Find the DV
                for(i=0; i<*track_n; i++)
                {
                    if(TrajectoryUsed(usedLineIndex, i))
                    {
                        a[0] = track_x0[i]; a[1] = track_y0[i]; a[2] = track_z0[i];
                        b[0] = track_x1[i]; b[1] = track_y1[i]; b[2] = track_z1[i];

                        for(j=i+1; j<*track_n; j++)
                        {  
                            if(TrajectoryUsed(usedLineIndex, j))
                            {
                                aa[0] = track_x0[j]; aa[1] = track_y0[j]; aa[2] = track_z0[j];
                                bb[0] = track_x1[j]; bb[1] = track_y1[j]; bb[2] = track_z1[j];

                                // Displaced Vertex Coordinates
                                displaced_Vertex = displacedVertex(a, b, aa, bb);
                                DV[0] = displaced_Vertex[0];
                                DV[1] = displaced_Vertex[1];
                                DV[2] = displaced_Vertex[2];

                                // Apply case for DVs
                                prodSum = CaseDV(DV, a, b, aa, bb);
                                
                                // epsilon1 < cosΘ_i <= epsilon2, epsilon1 < cosΘ_j <= epsilon2, 
                                if(prodSum>=2*epsilon1 && prodSum<2*epsilon2)
                                {
                                    distance_ij[count_j][0] = distance(a, b, aa, bb);

                                    distance_ij[count_j][1] = DV[0];
                                    distance_ij[count_j][2] = DV[1];
                                    distance_ij[count_j][3] = DV[2];

                                    distance_ij[count_j][4] = i;
                                    distance_ij[count_j][5] = j;
                                    
                                    count_j++;      
                                }
                            } 
                        }
                    }
                }   

                // The number of elements in distance_ij columns
                elementCount = count_j;

                // leastDistance[6] = (Distacne, DV_x, DV_y, DV_z, indexes_i, index_j)
                least_Distance = minimumArrayValue(distance_ij, elementCount);
                for(int k=0; k<6; k++)
                {
                    leastDistance[k] = least_Distance[k]; 
                }

                h2->Fill(leastDistance[0]);

                distance_xyz = sqrt(leastDistance[1]*leastDistance[1] + leastDistance[2]*leastDistance[2] + leastDistance[3]*leastDistance[3]);
                // distance_xyz = sqrt(truthvtx_x[0]*truthvtx_x[0] + truthvtx_y[0]*truthvtx_y[0] + truthvtx_z[0]*truthvtx_z[0]);

                h3->Fill(distance_xyz);

                DV_Z = leastDistance[3];

                h4->Fill(DV_Z);
                HH->Fill(distance_xyz, DV_Z);
                
                // Assigning coordinates to DV Array
                displacedVertexArray[0] = leastDistance[1]; // DV_x
                displacedVertexArray[1] = leastDistance[2]; // DV_y
                displacedVertexArray[2] = leastDistance[3]; // DV_z

                // TODO: Make a loop to find all possible errors and store the minimum for every displayed vertex excluding one DV_truth every time

                absoluteError = Error(displacedVertexArray[0], displacedVertexArray[1], displacedVertexArray[2], truthvtx_x[0], truthvtx_y[0], truthvtx_z[0]);

                h1->Fill(absoluteError); // Distance of calculated DV from truth DV (Error)
                H->Fill(leastDistance[0], absoluteError);

                // TODO: Check if this is viable when you work with more than 1 DVs

                // Calculating the distance from DV of the third trajectory (if there is one)
                if(*track_n>2)
                {
                    DvTrajectory = ThirdTrajectoryDistance(displacedVertexArray, *track_n, leastDistance[4], leastDistance[5]);
                }

                h5->Fill(DvTrajectory);

                // Store line indexes that have been used to calculate a DV
                for(int k=0; k<2; k++)
                {
                    usedLineIndex[countLine] = least_Distance[4+k];
                    countLine++;
                }
            }
            event++;
        }
    }

    // Canvas 1
    TCanvas *c1 = new TCanvas("c1", "DV Errors - Distance Between Trajectories - Distance of Dv from Third Trajectory - 3D Histogram", 900, 700);
    c1->Divide(2,2);

    gStyle->SetOptStat(1111111);

    c1->cd(1);
    h1->SetFillColor(kBlue-2);
    h1->SetMinimum(0);
    h1->Draw();

    c1->cd(3);
    h2->SetFillColor(kGreen-7);
    h2->SetMinimum(0);
    h2->Draw();

    c1->cd(2);
    H->SetFillColor(kRed);
    H->SetMinimum(0);
    H->Draw("LEGO1");

    c1->cd(4);
    h5->SetFillColor(kViolet);
    h5->SetMinimum(0);
    h5->Draw();

    c1->Print();

    // Canvas 2
    TCanvas *c2 = new TCanvas("c2", "Distance of Trajectories and DV from Detector's Center - Z Coordinate - 3D Histogram", 1300, 400);
    c2->Divide(3,1);

    c2->cd(1);
    h3->SetFillColor(kYellow);
    h3->SetMinimum(0);
    h3->Draw();

    c2->cd(2);
    h4->SetFillColor(kOrange+7);
    h4->SetMinimum(0);
    h4->Draw();

    c2->cd(3);
    HH->SetFillColor(kBlue-9);
    HH->SetMinimum(0);
    HH->Draw("LEGO1");

    c2->Print();

    // Print time needed for the program to complete
    printf("\nTime taken: %.2fs\n\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

    infile->Close();
    }
