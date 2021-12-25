// Functions for Vector Operations //

// Print array[3] elements - Used for coordinates of vectors in 3D space
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
TTreeReaderValue<Int_t> track_n = {treereader, "track_n"};              // Number of tracks
TTreeReaderValue<Int_t> truthvtx_n = {treereader, "truthvtx_n"};        // Number of DVs
// First point of trajectory
TTreeReaderArray<Double_t> track_x0 = {treereader, "track.x0"};         // x coordinate
TTreeReaderArray<Double_t> track_y0 = {treereader, "track.y0"};         // y coordinate
TTreeReaderArray<Double_t> track_z0 = {treereader, "track.z0"};         // z coordinate
// Second point of trajectory
TTreeReaderArray<Double_t> track_x1 = {treereader, "track.x1"};         // x coordinate
TTreeReaderArray<Double_t> track_y1 = {treereader, "track.y1"};         // y coordinate
TTreeReaderArray<Double_t> track_z1 = {treereader, "track.z1"};         // z coordinate
// The DV_truth
TTreeReaderArray<Double_t> truthvtx_x = {treereader, "truthvtx.x"};     // x coordinate
TTreeReaderArray<Double_t> truthvtx_y = {treereader, "truthvtx.y"};     // y coordinate
TTreeReaderArray<Double_t> truthvtx_z = {treereader, "truthvtx.z"};     // z coordinate

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


// Computes the minimum value of array's first column distance[elementCount][6] and
// leaves other columns' elements intact.
double *minimumArrayValueSix(double distance[][6], int elementCount)
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


// Cases to choose if a point is a DV or not
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


// Returns true if at least one element in array equals to i
// and false if all elements in array are different to i.
bool IndexUsed(int *array, int i)
{
    bool used = false;

    // 30 elements in array
    for(int k=0; k<30; k++)
    {
        if(i == array[k])
        {
            used = true;
            break;
        }
    }

    return used;
}


// Computes the minimum value of array's first column error[elementCount][2] and 
// leaves other columns' elements intact. //! The min must be different than -1
double *minimumArrayValueTwo(double error[][2], int elementCount)
{
    // First element the minimum value. Second the index of DV_truth used
    static double minimum[2];

    // Initialize with the first row of error[elementCount][2]
    minimum[0] = error[0][0];
    minimum[1] = error[0][1];

    // Go throught the error array
    double element;

    for(int i=1; i<elementCount; i++)
    {
        element = error[i][0];

        if(element != -1 && minimum[0] >= element)
        {
            minimum[0] = element;
            minimum[1] = error[i][1];
        }
    }

    return minimum;
}


// ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ //

void new_myAnalyzeStage1()
{
    // For time capture
    clock_t tStart = clock();

    // Histograms
    TH1 *error_XYZ = new TH1D("error_XYZ", "Error in 3D Space;Error;Counts", 100, 0, 35);
    TH1 *error_XY = new TH1D("error_XY", "Error in xy Plane;Error;Counts", 50, 0, 14);

    TFile* infile = TFile::Open("stage1.root");
    TTree* tree   = (TTree*)infile->Get("stage1");
    treereader.SetTree(tree);


    //! Search for DVs !//
    // Condition to decide if a trajectory belongs to a DV without constructing it 
    double Dcut1 = 0.2;
    // Condition to decide if two trajectories form a DV
    double Dcut2 = 0.35;

    // Line_i Points
    double a[3], b[3];
    // Line_j Points
    double aa[3], bb[3];

    // Array that gathers the distances between line i and j and stores them in the first column.
    // The second, third and forth columns store the coordinates of the displaced vertex resulting from
    // line_i and line_j. The fifth and sixth store the i-th and j-th lines' indexes, respactively.
    double distance_ij[100][6];
    // Greater than zero element counter for distance_ij array
    int elementCount;
    // Counter for distance_ij array elemets
    int count_j;

    // The minimum distance of all the possible pair of trajectories for every event
    // The first column contains the least distance for every event. The other three store the coordinates
    // of the displaced vertex resulting from the lines that produce the minimum distance. The fifth and 
    // sixth store the i-th and j-th lines' indexes, respactively.
    double leastDistance[6];
    double *least_Distance;
    // The minimum of all possible pair of trajectories
    double leastdistance;

    // Stores indexes of lines that have been used to calculate a DV
    int usedLineIndex[30];
    // Counter for usedLineIndex[30] elements
    int countLine;

    // Stores the coordinates of each displaced vertex of any event
    double displacedVertexArray[3];
    double *displaced_Vertex;
    // Coordinates for displaced vertex produced by line_i and line_j
    double DV[3];

    // To apply the case for DVs
    double prodSum;
    // Limits
    double theta = 90;
    double epsilon1 = cos(M_PI * theta/180);
    double epsilon2 = cos(0);


    //! Errors in xy Plane !//
    // First Column: Errors of DV_truth[i] with DV_reco (in i^th row),
    // Second Column: Index of DV_truth used in rows, respectively.
    double errorXY[10][2];
    // First Column: The minimum error for each DV_reco,
    // Second Column: Index of DV_truth used.
    double minErrorXY[2];
    double *min_ErrorXY;

    //! Errors in 3D Space !//
    // First Column: Errors of DV_truth[i] with DV_reco (in i^th row),
    // Second Column: Index of DV_truth used in rows, respectively.
    double errorXYZ[10][2];
    // Counter for errorXYZ array elements
    int errorCounter;
    // First Column: The minimum error for each DV_reco,
    // Second Column: Index of DV_truth used.
    double minErrorXYZ[2];
    double *min_ErrorXYZ;

    // Stores the indexes of the used DV_thuth
    int usedErrorIndex[30];
    // Counter for usedErrorIndex array elements
    int indexCounter;


    //!  !//
    // The minimun distance of the third trajectory to the DV
    double DvTrajectory;
    //!  !//


    //! Other Parameters !//
    // Event Counter
    int event = 0;

    // Integers for for loops
    int i, j;


    // ------------------------------------------------------------------------------------------------------------------------------------------------------ //

    // Event Loop
    while (treereader.Next()) 
    {
        // Loop in events with 1 DV
        if(*truthvtx_n==1)
        {   
            // Renew for every event
            countLine = 0;
            indexCounter = 0;
            // Initialize all values to -1 so as not to coincide with other events
            for(int k=0; k<30; k++)
            {
                usedLineIndex[k] = -1;
                usedErrorIndex[k] = -1;
            }

            // Loop to find multiple DVs. If leastdistance <= Dcut2 
            do
            {
                count_j = 0;
                leastdistance = 1; // TODO: Check if this is needed
                errorCounter = 0;
                for(int k=0; k<10; k++)
                {
                    errorXYZ[k][0] = -1;
                    errorXYZ[k][1] = -1;
                    // 
                    errorXY[k][0] = -1;
                    errorXY[k][1] = -1;
                }

                //! Find the DV
                for(i=0; i<*track_n; i++)
                {
                    if(!IndexUsed(usedLineIndex, i))
                    {
                        a[0] = track_x0[i]; a[1] = track_y0[i]; a[2] = track_z0[i];
                        b[0] = track_x1[i]; b[1] = track_y1[i]; b[2] = track_z1[i];

                        for(j=i+1; j<*track_n; j++)
                        {  
                            if(!IndexUsed(usedLineIndex, j))
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
                least_Distance = minimumArrayValueSix(distance_ij, elementCount);
                for(int k=0; k<6; k++)
                {
                    leastDistance[k] = least_Distance[k]; 
                }

                // Condition to decide if there is a DV
                leastdistance = leastDistance[0];
                if(leastdistance > Dcut2)
                {
                    break;
                }
                
                // Store line indexes that have been used to calculate a DV
                for(int k=4; k<6; k++)
                {
                    usedLineIndex[countLine] = least_Distance[k];
                    countLine++;
                }
                
                // Assigning coordinates to DV Array
                displacedVertexArray[0] = leastDistance[1]; // DV_x
                displacedVertexArray[1] = leastDistance[2]; // DV_y
                displacedVertexArray[2] = leastDistance[3]; // DV_z

                //! Compute Errors
                for(int k=0; k<*truthvtx_n; k++)
                {
                    if(!IndexUsed(usedErrorIndex, k))
                    {
                        // Error of DV_truth[k] and DV_reco
                        errorXYZ[errorCounter][0] = Error(displacedVertexArray[0], displacedVertexArray[1], displacedVertexArray[2], truthvtx_x[k], truthvtx_y[k], truthvtx_z[k]); 
                        errorXY[errorCounter][0] = Error(displacedVertexArray[0], displacedVertexArray[1], 0, truthvtx_x[k], truthvtx_y[k], 0); 
                        // Index of DV_truth used
                        errorXYZ[errorCounter][1] = k;
                        errorXY[errorCounter][1] = k;

                        errorCounter++;
                    }
                }

                min_ErrorXYZ = minimumArrayValueTwo(errorXYZ, 10);
                minErrorXYZ[0] = min_ErrorXYZ[0]; // Minimum Error 
                minErrorXYZ[1] = min_ErrorXYZ[1]; // Index of DV_truth used

                min_ErrorXY = minimumArrayValueTwo(errorXY, 10);
                minErrorXY[0] = min_ErrorXY[0]; // Minimum Error 
                minErrorXY[1] = min_ErrorXY[1]; // Index of DV_truth used

                usedErrorIndex[indexCounter] = minErrorXYZ[1];
                indexCounter++;

                error_XYZ->Fill(minErrorXYZ[0]); // Distance of calculated DV from truth DV (Error)
                error_XY->Fill(minErrorXY[0]); // Distance of calculated DV from truth DV (Error)

                //! Calculating the distance from DV of the third trajectory (if there is one)
                if(*track_n>2)
                {
                    DvTrajectory = ThirdTrajectoryDistance(displacedVertexArray, *track_n, leastDistance[4], leastDistance[5]);
                }

            } while(leastdistance > Dcut2);

            event++;
        }
    }

    // Canvas 1
    TCanvas *c1 = new TCanvas("c1", "DV Errors in XYZ Space and xy Plane", 1400, 400);
    c1->Divide(2,1);

    gStyle->SetOptStat(1111111);

    c1->cd(1);
    error_XYZ->SetFillColor(kBlue-2);
    error_XYZ->SetMinimum(0);
    error_XYZ->Draw();

    c1->cd(2);
    error_XY->SetFillColor(kRed);
    error_XY->SetMinimum(0);
    error_XY->Draw();

    c1->Print();

    // Print time needed for the program to complete
    printf("\nTime taken: %.2fs\n\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

    infile->Close();
    }
