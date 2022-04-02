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


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

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

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //


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

        if(minimum[0] >= element && element != -1)
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

        if(minimum >= element && element != -1)
        {
            minimum = element;
        }
    }

    return minimum;
}


// Cases to choose if a point is a DV or not with respect to relative angles between the lines which
// are defined by DV-A, DV-AA and DV-B, DV-BB. Returns the sum of angles in degrees.
double *CaseDV_RelativeAngles(double *DV, double *a, double *b, double *aa, double *bb)
{
    // Relative Angles
    static double theta[2];

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

    // Θ_i and Θ_j in degrees
    theta[0] = acos(dotProd_i) * 180.0 / M_PI;
    theta[1] = acos(dotProd_j) * 180.0 / M_PI;

    return theta;
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


// Computes the minimum value of array's first column array[elementCount][2] and 
// leaves other columns' elements intact. //! The min must be different than -1
double *minimumArrayValueTwo(double array[][2], int elementCount)
{
    // First element the minimum value. Second the index of DV_truth used
    static double minimum[2];

    // Initialize with the first row of array[elementCount][2]
    minimum[0] = array[0][0];
    minimum[1] = array[0][1];

    // Go throught the array array
    double element;

    for(int i=1; i<elementCount; i++)
    {
        element = array[i][0];

        if(element != -1 && minimum[0] >= element)
        {
            minimum[0] = element;
            minimum[1] = array[i][1];
        }
    }

    return minimum;
}


// Sorts the array elements in ascending order
void ArraySorting(double array[][2], int elementCount)
{
    double temp0;
    double temp1;

    for(int i=0; i<elementCount; i++)
    {
        for(int j=i; j<elementCount; j++)
        {
            temp0 = array[i][0];

            if(array[j][0] < temp0)
            {
                temp0 = array[j][0]; 
                temp1 = array[j][1];

                array[j][0] = array[i][0];
                array[j][1] = array[i][1];
                
                array[i][0] = temp0;
                array[i][1] = temp1;
            }
        }
    }
}


// Maximum Angle Between vector ODV and A_iAA_i, B_jBB_j, where A_i, AA_i are points of line_i and
// B_j, BB_j are points of line_j
double AngleDVTracks(double *DV, double *a, double *aa, double *b, double *bb)
{
    // Angles Between ODV and A_iAA_i, B_iBB_i
    double theta[2];
    // The maximun of angles
    double theta_max;

    // Center of Axis - Interaction Point
    double O[3] = {0,0,0};

    // Unit vector from O to DV
    double *unit_ODV = relativeUnitVector(O, DV);
    double unitODV[3] = {unit_ODV[0], unit_ODV[1], unit_ODV[2]};

    // Relative unit vector from a to aa
    double *unitRel_aaa = relativeUnitVector(a, aa);
    double unitRelaaa[3] = {unitRel_aaa[0], unitRel_aaa[1], unitRel_aaa[2]};

    // Relative unit vector from b to bb
    double *unitRel_bbb = relativeUnitVector(b, bb);
    double unitRelbbb[3] = {unitRel_bbb[0], unitRel_bbb[1], unitRel_bbb[2]};

    // cos(θ_a) where θ_a is the angle between ODV and aaa vectors
    double dotProd_a = dotProduct(unitODV, unitRelaaa);
    // cos(θ_b) where θ_b is the angle between ODV and bbb vectors
    double dotProd_b = dotProduct(unitODV, unitRel_bbb);

    // Θ_a and Θ_b in degrees
    theta[0] = acos(dotProd_a) * 180.0 / M_PI;
    theta[1] = acos(dotProd_b) * 180.0 / M_PI;

    if(theta[0]>=theta[1]) 
    {
        theta_max = theta[0];
    }
    else
    {
        theta_max = theta[1];
    }

    return theta_max;
}

// ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ // ~~~~ //

void myAnalyzeStage1()
{
    // For time capture
    clock_t tStart = clock();

    //! Histograms for Canvas 1 - Errors in XYZ Space and xy Plane !//
    TH1 *error_XYZ = new TH1D("error_XYZ", "Error in 3D Space;Error;Counts", 100, 0, 35);
    TH1 *error_XYZ_OneDV = new TH1D("error_XYZ_OneDV", "Error in 3D Space for One DV;Error;Counts", 100, 0, 35);
    TH1 *error_XYZ_TwoDVs = new TH1D("error_XYZ_TwoDVs", "Error in 3D Space for Two DVs;Error;Counts", 100, 0, 35);
    TH1 *error_XY = new TH1D("error_XY", "Error in xy Plane;Error;Counts", 50, 0, 14);
    TH1 *error_XY_OneDV = new TH1D("error_XY_OneDV", "Error in xy Plane for One DV;Error;Counts", 50, 0, 14);
    TH1 *error_XY_TwoDVs = new TH1D("error_XY_TwoDVs", "Error in xy Plane for Two DVs;Error;Counts", 50, 0, 14);

    //! Histograms for Canvas 2 - Number of DV_total and DV_close !//
    TH1 *DvNumber = new TH1D("DvNumber", "Total Number of DV_true;DV_true;Counts", 41, 0, 4);
    TH1 *clarity = new TH1D("clarity", "DV_reco Independent of Distance from DV_true;DV_reco;Counts", 61, 0, 6);
    TH1 *performance = new TH1D("performance", "DV_reco that are Close to DV_true;DV_reco;Counts", 61, 0, 6);
    TH1 *OffErrorPerformance = new TH1D("OffErrorPerformance", "DV_reco Out of Limit Territory;DV_reco;Counts", 61, 0, 6);
    // One DV
    TH1 *DvNumber_OneDV = new TH1D("DvNumber_OneDV", "Number of DV_true = 1;DV_true;Counts", 41, 0, 4);
    TH1 *clarity_OneDV = new TH1D("clarity_OneDV", "DV_reco Independent of Distance from DV_true - One DV;DV_reco;Counts", 61, 0, 6);
    TH1 *performance_OneDV = new TH1D("performance_OneDV", "DV_reco that are Close to DV_true - One DV;DV_reco;Counts", 61, 0, 6);
    TH1 *OffErrorPerformance_OneDV = new TH1D("OffErrorPerformance_OneDv", "DV_reco Out of Limit Territory - One DV;DV_reco;Counts", 61, 0, 6);
    // Two DVs
    TH1 *DvNumber_TwoDVs = new TH1D("DvNumber_TwoDVs", "Number of DV_true = 2;DV_true;Counts", 41, 0, 4);
    TH1 *clarity_TwoDVs = new TH1D("clarity_TwoDVs", "DV_reco Independent of Distance from DV_true - Two DVs;DV_reco;Counts", 61, 0, 6);
    TH1 *performance_TwoDVs = new TH1D("performance_TwoDVs", "DV_reco that are Close to DV_true - Two DVs;DV_reco;Counts", 61, 0, 6);
    TH1 *OffErrorPerformance_TwoDVs = new TH1D("OffErrorPerformance_TwoDvs", "DV_reco Out of Limit Territory - Two DVs;DV_reco;Counts", 61, 0, 6);

    //! Histograms for Canvas 3 - Relative Number of DV_reco with respect to DV_true !//
    TH1 *RelativeNumber = new TH1D("RelativeNumber", "Relative Number of DV_reco and DV_true;DV_true-DV_reco;Counts", 101, -5, 5);
    TH1 *RelativeNumber_OneDV = new TH1D("RelativeNumber_OneDV", "Relative Number of DV_reco and DV_true - One DV;DV_true-DV_reco;Counts", 101, -5, 5);
    TH1 *RelativeNumber_TwoDVs = new TH1D("RelativeNumber_TwoDVs", "Relative Number of DV_reco and DV_true - Two DVs;DV_true-DV_reco;Counts", 101, -5, 5);

    //! Histograms for Canvas 4 - Distances Between DVtrue in Events with Two DVs !//
    TH1 *DistanceBetweenTwoDVs = new TH1D("DistanceBetweenTwoDVs", "Distance Between DVs in Events with Two DVtrue;Distance;Counts", 100, 0, 100);

    //! Histograms for Canvas 5 - Minimum Distance Between DVreco and Beginning of Lines Used to Reconstruct It !//
    TH1 *DVreco_RecoLinesMinDistance_Matched = new TH1D("DVreco_RecoLinesMinDistance_Matched", "Minimum Distance Between DVreco and Reconstructing Lines (Matched);Distance;Counts", 200, 0, 100);
    TH1 *DVreco_RecoLinesMinDistance_NotMatched = new TH1D("DVreco_RecoLinesMinDistance_NotMatched", "Minimum Distance Between DVreco and Reconstructing Lines (Not Matched);Distance;Counts", 200, 0, 100);

    //! Histograms for Canvas 6 -  Maximum Angle Between ODV and A_iAA_i, B_jBB_j vectors !//
    TH1 *Angle_ODV_AB_max_Matched = new TH1D("Angle_ODV_AB_max_Matched", "Maximum Angle Between ODV and A_iAA_i, B_jBB_j vectors (Matched);Angle;Counts", 89, 0, 180);
    TH1 *Angle_ODV_AB_max_NotMatched = new TH1D("Angle_ODV_AB_max_NotMatched", "Maximum Angle Between ODV and A_iAA_i, B_jBB_j vectors (Not Matched);Angle;Counts", 89, 0, 180);

    //! Histogram for Canvas 7 - Distance Between Tracks Used to Reconstruct the DVreco !//
    TH1 *Track_Distance_Matched = new TH1D("Track_Distance_Matched", "Distance Between Tracks Used to Reconstruct the DVreco (Mathced);Track Distance;Counts", 200, 0, 1);
    TH1 *Track_Distance_NotMatched = new TH1D("Track_Distance_NotMatched", "Distance Between Tracks Used to Reconstruct the DVreco (Not Mathced);Track Distance;Counts", 200, 0, 1);

    //! Histogram for Canvas 8 - Distance of Closest to DVreco Additional Track (If It Exists) !//
    TH1 *ClosestTrackDV_Matched = new TH1D("ClosestTrackDV_Matched", "Distance of Closest to DVreco Additional Track (Matched);Track to DV Distance;Counts", 200, 0, 100);
    TH1 *ClosestTrackDV_NotMatched = new TH1D("ClosestTrackDV_NotMatched", "Distance of Closest to DVreco Additional Track (Not Mathced);Track to DV Distance;Counts", 200, 0, 100);

    //! Histogram for Canvas 9 - Difference of Distances of DVreco (R_DVreco) and Point of Reconstructing Tracks (Rmin) from IT !//
    TH1 *RDVreco_Rmin_Matched = new TH1D("RDVreco_Rmin_Matched", "Difference of Distances of DVreco and Point of Reconstructing Tracks from IT (Matched);R_DVreco - Rmin;Counts", 200, -20, 20);
    TH1 *RDVreco_Rmin_NotMatched = new TH1D("RDVreco_Rmin_NotMatched", "Difference of Distances of DVreco and Point of Reconstructing Tracks from IT (Not Matched);R_DVreco - Rmin;Counts", 200, -40, 40);

    //! Histogram for Canvas 10 - Maximum Relative Angle Between DVreco and Given Points of Reconstructed Tracks !//
    TH1 *Relative_Angle_Matched = new TH1D("Relative_Angle_Matched", "Maximum Relative Angle Between DVreco and Given Points of Reconstructed Tracks (Matched);Relative Angle;Counts", 100, 0, 180);
    TH1 *Relative_Angle_NotMatched = new TH1D("Relative_Angle_NotMatched", "Maximum Relative Angle Between DVreco and Given Points of Reconstructed Tracks (Not Matched);Relative Angle;Counts", 100, 0, 180);

    TFile* infile = TFile::Open("stage1.root");
    TTree* tree   = (TTree*)infile->Get("stage1");
    treereader.SetTree(tree);

    //! Search for DVs !//
    // Line_i Points
    double a[3], b[3];
    // Line_j Points
    double aa[3], bb[3];

    // First Column: Distance between line i and j.
    // Second - Fourth Column: Coordinates of the displaced vertex resulting from line_i and line_j. 
    // Fifth and Sixth Column: The i-th and j-th lines' indexes, respactively.
    double distance_ij[100][6];
    // Counter for distance_ij array elemets
    int count_j;

    // First Column: the least distance for every event.
    // Second - Fourth Column: Coordinates of the DV resulting from the lines that produce the minimum distance. 
    // Fifth and Sixth Column: Lines indexes that produce the least distance.
    double leastDistance[6];
    double *least_Distance;

    // Stores indexes of lines that have been used to calculate a DV
    int usedLineIndex[30];
    // Counter for usedLineIndex[30] elements
    int countLine;

    // Stores the coordinates of each displaced vertex of any event
    double displacedVertexArray[3];
    double *displaced_Vertex;
    // Coordinates for displaced vertex produced by line_i and line_j
    // Used for application of the condition for the DVs
    double DV[3];

    //! Application of Restrictions  !//
    // Condition to decide if two trajectories form a DV //? DVcut = 0.08
    double DVcut = 1;
    // Condition to decide if a trajectory belongs to a DV without constructing it 
    double TrajectoryCut = DVcut/2;

    // Relative Angles
    double *Angles_Rel;
    double AnglesRel[2];
    // Limits
    double thetaRel_max = 90;
    double thetaRel_min = 0;

    //! Errors in xy Plane !//
    // First Column: Errors of DV_truth[i] with DV_reco (in i^th row),
    // Second Column: Index of DV_truth used in rows, respectively.
    double errorXY[10][2];
    // First Column: The minimum error for each DV_reco,
    // Second Column: Index of DV_truth used.
    double minErrorXY[2];
    double *min_ErrorXY;

    //! DV Counters !//
    // Counter for DV_truth // 

    // The total number of DV_truth in all events
    int DV_Total = 0;
    // The total number of DV_truth in event with one DV
    int DVnumber_OneDV_Total = 0;
    // The total number of DV_truth in event with two DVs
    int DVnumber_TwoDVs_Total = 0;
    // The total number of DV_truth in event with three DVs
    int DVnumber_ThreeDVs_Total = 0;
    // The total number of DV_truth in event with four DVs
    int DVnumber_FourDVs_Total = 0;
    // Total number of DV_true matched by DV_reco that respect limit in xy plane //
    int DV_true_with_match_XY = 0;
    int DV_true_with_match_XY_OneDV = 0;
    int DV_true_with_match_XY_TwoDVs = 0;
    // Total number of DV_true matched by DV_reco that respect limit in sz space //
    int DV_true_with_match_SZ = 0;
    int DV_true_with_match_SZ_OneDV = 0;
    int DV_true_with_match_SZ_TwoDVs = 0; 

    // Counter for DV_reco //

    // The total number of DV_reco in all events
    int DV_reco_total = 0;
    // The total number of DV_reco in events with one DV
    int DV_reco_OneDV_Total = 0;
    // The total number of DV_reco in events with two DVs
    int DV_reco_TwoDVs_Total = 0;
    // Number of DV_reco that respect xy limit
    int DV_reco_xy = 0;
    // Number of DV_reco that respect xy limit- Events with one DV
    int DV_reco_xy_OneDV = 0;
    // Number of DV_reco that respect xy limit - Events with one DV
    int DV_reco_xy_TwoDVs = 0;
    // Number of DV_reco that respect sz limit
    int DV_reco_sz = 0;
    // Number of DV_reco that respect sz limit - Events with one DV
    int DV_reco_sz_OneDV = 0;
    // Number of DV_reco that respect sz limit - Events with two DVs
    int DV_reco_sz_TwoDVs = 0;

    // Counter for DVs in single events //

    // The total number of DV_reco in an event
    int DVnumber_Total;
    // The number of DV_reco that respect both limits (xy and sz) in an event
    int DVnumber_Close;
    // The total numbers of DV_reco that do not respect any limit (neither xy nor sz)
    int DVnumber_Far;


    //! Efficiency: DV_reco_matched/DV_truth_total and Purity: DV_reco_matched/DV_reco_total !//
    // DV_reco that respect limits in xy plane
    double Eff_xy;
    double Eff_xy_OneDV;
    double Eff_xy_TwoDVs;
    // DV_reco that respect limits in sz space
    double Eff_sz;
    double Eff_sz_OneDV;
    double Eff_sz_TwoDVs;
    // DV_reco that respect limits in xy plane
    double Pu_xy;
    double Pu_xy_OneDV;
    double Pu_xy_TwoDVs;
    // DV_reco that respect limits in sz space
    double Pu_sz;
    double Pu_sz_OneDV;
    double Pu_sz_TwoDVs;

    //! Accuracy = Errors_Respect_Limits/Errors_Computed !//
    double accuracy_total;
    double accuracy_OneDV;
    double accuracy_TwoDVs;

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

    // Stores the indexes of the used DV_truth
    int usedErrorIndex[30];
    // Counter for usedErrorIndex array elements
    int errorIndexCounter;

    //! Distance from a DV to a different Trajectory from the two that constructed it !//
    // Points that define the trajectory
    double A[3], B[3];
    // First Column: Stores the distances between trajectory_i and DV
    // Second Column: Stores the index of the trajectory used, respectively
    double DvTrajectoryDistance[20][2];
    // DvTrajectoryDistance array element number
    int DvTrajectoryCounter;
    // Distance of "third" trajectory to DV
    double thirdTrajectoryDistance;

    //! DVreco Distance from Lines that Reconstructed It !//
    double DVrecoLine[4];
    double DVrecoLine_min;

    //! Angle Between ODV and A_iAA_i, B_iBB_i vectors !//
    double A_i[3], AA_i[3];
    double B_j[3], BB_j[3];
    double theta_ODVAB_max;

    //! Distance Between Tracks Used to Reconstruct the DVreco !//
    double Track_Distance_DVrecoMatch;
    double Track_Distance_DVrecoNoMatch;

    //! Distance of Closest to DVreco Additional Track (If It Exists) !//
    double ClosestTrack_DV;

    //! Difference of Distances of DVreco (R_DVreco) and Point of Reconstructing Tracks (Rmin) from IT !//
    double R_DVreco;
    double R_track[4];
    double Rmin;

    //! Maximum Relative Angle Between DVreco and Given Points of Reconstructed Tracks !//
    double *Relative_Angle;
    double RelativeAngle[2];

    //! Miscellaneous !//
    // Conditions for DV_reco that are "close"
    int limitXYZ = 35;
    int limitXY = 14;
    double distanceBetweenTwoDVs;
    int ErrorsComputed = 0;
    int ErrorsComputed_OneDV = 0;
    int ErrorsComputed_TwoDVs = 0;

    // Event Counter
    int event = 1;

    // Integers for for loops
    int line_i, line_j, line_k, er;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

    // Event Loop
    while (treereader.Next()) 
    {
        // Loop in events with multiple DVs
        if(*truthvtx_n >= 1)
        {   
            //! Renew for every event !//
            // Total number of DV_reco
            DVnumber_Total = 0;  
            // Number of DV_reco that respect the limits in errors
            DVnumber_Close = 0; 
            // Number of DV_reco that are out of error limits
            DVnumber_Far = 0; 
            // element counter for usedLineIndex array
            countLine = 0; 
            // element counter for usedErrorIndex array
            errorIndexCounter = 0; 

            // Initialize all values to -1 so as not to coincide with other events
            for(int k=0; k<30; k++)
            {
                usedLineIndex[k] = -1;
                usedErrorIndex[k] = -1;
            }

            //! Loop to find multiple DVs !//
            do
            {
                //! Renew for every DV !//
                // DvTrajectoryDistance element counter
                DvTrajectoryCounter = 0;
                // Array distance_ij element counter
                count_j = 0;
                // Arrays errorXYZ and errorXY element counter
                errorCounter = 0;
                // Initialize to -1 to avoid coincidence with other DVs
                ClosestTrack_DV = -1;
                // Initialize to -1 to avoid coincidence with other DVs
                RelativeAngle[0] = -1000; RelativeAngle[1] = -1000;

                // Initialize fot Multiple Trajectories
                for(int k=0; k<20; k++)
                {
                    DvTrajectoryDistance[k][0] = -1;
                    DvTrajectoryDistance[k][1] = -1;
                }
                // Initialize all to -1 so as not to coincide with other DVs
                for(int k=0; k<10; k++)
                {
                    errorXYZ[k][0] = -1;
                    errorXYZ[k][1] = -1;
                    errorXY[k][0] = -1;
                    errorXY[k][1] = -1;
                }
                // Initianlize trajectories distances
                for(int k=0; k<100; k++)
                {
                    for(int l=0; l<6; l++)
                    {
                        distance_ij[k][l] = -1;
                    }
                }

                //! As the number of DVreco increases should be exponentially more difficult to find other DVs !//
                // DVcut = 1;
                // for(int k=1; k<7; k++)
                // {
                //     if(DVnumber_Total == k)
                //     {
                //         DVcut = DVcut * exp(-2*k);
                //         TrajectoryCut = DVcut/2;
                //     }
                // }

                //! Find the DV !//
                for(line_i=0; line_i<*track_n; line_i++)
                {
                    if(!IndexUsed(usedLineIndex, line_i))
                    {
                        // line_i
                        a[0] = track_x0[line_i]; a[1] = track_y0[line_i]; a[2] = track_z0[line_i];
                        b[0] = track_x1[line_i]; b[1] = track_y1[line_i]; b[2] = track_z1[line_i];

                        for(line_j=line_i+1; line_j<*track_n; line_j++)
                        {  
                            if(!IndexUsed(usedLineIndex, line_j))
                            {
                                // line_j
                                aa[0] = track_x0[line_j]; aa[1] = track_y0[line_j]; aa[2] = track_z0[line_j];
                                bb[0] = track_x1[line_j]; bb[1] = track_y1[line_j]; bb[2] = track_z1[line_j];

                                // DV Coordinates formed by line_i and line_j
                                displaced_Vertex = displacedVertex(a, b, aa, bb);
                                DV[0] = displaced_Vertex[0];
                                DV[1] = displaced_Vertex[1];
                                DV[2] = displaced_Vertex[2];

                                //! Application of Restrictions !//
                                // Relative Angle 
                                Angles_Rel = CaseDV_RelativeAngles(DV, a, b, aa, bb);
                                AnglesRel[0] = Angles_Rel[0]; AnglesRel[1] = Angles_Rel[1]; 
                            
                                if(AnglesRel[0]>=thetaRel_min && AnglesRel[0]<=thetaRel_max && AnglesRel[1]>=thetaRel_min && AnglesRel[1]<=thetaRel_max)
                                {
                                    // Distance between line_i and line_j
                                    distance_ij[count_j][0] = distance(a, b, aa, bb);

                                    // DV coordinates
                                    distance_ij[count_j][1] = DV[0];
                                    distance_ij[count_j][2] = DV[1];
                                    distance_ij[count_j][3] = DV[2];

                                    // Used lines indexes
                                    distance_ij[count_j][4] = line_i;
                                    distance_ij[count_j][5] = line_j;
                                    
                                    count_j++;      
                                }
                            } 
                        }
                    }
                }  

                // leastDistance[6] = (Distacne, DV_x, DV_y, DV_z, indexes_i, index_j)
                least_Distance = minimumArrayValueSix(distance_ij, count_j);
                for(int k=0; k<6; k++)  leastDistance[k] = least_Distance[k]; 

                //! Condition to decide if there is a DVreco !//
                if(leastDistance[0] < DVcut && leastDistance[0] >= 0)
                {
                    // If it passes the previous condition it means that we have a DV
                    DVnumber_Total++;

                    // Store line indexes that have been used to calculate the DV
                    for(int k=4; k<6; k++)
                    {
                        usedLineIndex[countLine] = least_Distance[k];
                        countLine++;
                    }
                    
                    // Assigning coordinates to DV Array
                    displacedVertexArray[0] = leastDistance[1]; // DV_x
                    displacedVertexArray[1] = leastDistance[2]; // DV_y
                    displacedVertexArray[2] = leastDistance[3]; // DV_z

                    //! Minimum Distance Between DVreco and Beginning of Lines Used to Reconstruct It !//
                    // Compute distance between DVreco and lines that reconstructed it
                    // Line_i
                    DVrecoLine[0] = Error(leastDistance[1], leastDistance[2], leastDistance[3], track_x0[leastDistance[4]], track_y0[leastDistance[4]], track_z0[leastDistance[4]]);
                    DVrecoLine[1] = Error(leastDistance[1], leastDistance[2], leastDistance[3], track_x1[leastDistance[4]], track_y1[leastDistance[4]], track_z1[leastDistance[4]]);
                    // Line_j
                    DVrecoLine[2] = Error(leastDistance[1], leastDistance[2], leastDistance[3], track_x0[leastDistance[5]], track_y0[leastDistance[5]], track_z0[leastDistance[5]]);
                    DVrecoLine[3] = Error(leastDistance[1], leastDistance[2], leastDistance[3], track_x1[leastDistance[5]], track_y1[leastDistance[5]], track_z1[leastDistance[5]]);
                    // Minimum Distance
                    DVrecoLine_min = minimumValuefromArrayElements(DVrecoLine, 4);

                    //! Maximum Angle Between ODV and A_iAA_i, B_iBB_i vectors !//   
                    // Line_i
                    A_i[0] = track_x0[leastDistance[4]]; A_i[1] = track_y0[leastDistance[4]]; A_i[2] = track_z0[leastDistance[4]];
                    AA_i[0] = track_x1[leastDistance[4]]; AA_i[1] = track_y1[leastDistance[4]]; AA_i[2] = track_z1[leastDistance[4]];
                    // Line_j
                    B_j[0] = track_x0[leastDistance[5]]; B_j[1] = track_y0[leastDistance[5]]; B_j[2] = track_z0[leastDistance[5]];
                    BB_j[0] = track_x1[leastDistance[5]]; BB_j[1] = track_y1[leastDistance[5]]; BB_j[2] = track_z1[leastDistance[5]];
                    // Maximum Angle
                    theta_ODVAB_max = AngleDVTracks(displacedVertexArray, A_i, AA_i, B_j, BB_j);

                    //! Difference of Distances of DVreco (R_DVreco) and Point of Reconstructing Tracks (Rmin) from IT !//
                    // Distance of DVreco to IT
                    R_DVreco = Error(leastDistance[1],leastDistance[2],leastDistance[3],0,0,0);
                    // Distance of Line_i Points to IT
                    R_track[0] = Error(track_x0[leastDistance[4]],track_y0[leastDistance[4]],track_z0[leastDistance[4]],0,0,0);
                    R_track[1] = Error(track_x1[leastDistance[4]],track_y1[leastDistance[4]],track_z1[leastDistance[4]],0,0,0);
                    // Distance of Line_j Points to IT
                    R_track[2] = Error(track_x0[leastDistance[5]],track_y0[leastDistance[5]],track_z0[leastDistance[5]],0,0,0);
                    R_track[3] = Error(track_x1[leastDistance[5]],track_y1[leastDistance[5]],track_z1[leastDistance[5]],0,0,0);
                    // Minimun of Lines Points Distance from IT
                    Rmin = minimumValuefromArrayElements(R_track, 4);

                    //! Maximum Relative Angle Between DVreco and Given Points of Reconstructed Tracks !//
                    Relative_Angle = CaseDV_RelativeAngles(displacedVertexArray, A_i, AA_i, B_j, BB_j);
                    RelativeAngle[0] = Relative_Angle[0]; RelativeAngle[1] = Relative_Angle[1];

                    //! Condition to take into consideration multiple trajectories that might belong to the same DV !//
                    if(*track_n - countLine >= 1)
                    {
                        for(line_k=0; line_k<*track_n; line_k++)
                        {
                            if(!IndexUsed(usedLineIndex, line_k))
                            {
                                // line_i
                                A[0] = track_x0[line_k]; A[1] = track_y0[line_k]; A[2] = track_z0[line_k];
                                B[0] = track_x1[line_k]; B[1] = track_y1[line_k]; B[2] = track_z1[line_k];

                                // Distance from the DV to line_i
                                DvTrajectoryDistance[DvTrajectoryCounter][0] =  LinePointDistance(A, B, DV);
                                // Index of line_i
                                DvTrajectoryDistance[DvTrajectoryCounter][1] = line_k;

                                DvTrajectoryCounter++;
                            }
                        }

                        // Sort the DvTrajectoryDistance elements in ascending order with respect to distances 
                        ArraySorting(DvTrajectoryDistance, DvTrajectoryCounter);

                        // The closest distance of an additional track to DVreco
                        ClosestTrack_DV = DvTrajectoryDistance[0][0];

                        for(int k=0; k<DvTrajectoryCounter; k++)
                        {
                            // The distance of third trajectory
                            thirdTrajectoryDistance = DvTrajectoryDistance[k][0];

                            if(thirdTrajectoryDistance <= TrajectoryCut)
                            {
                                // The index of trajectory used
                                usedLineIndex[countLine] = DvTrajectoryDistance[k][1];
                                countLine++;
                            }
                            else
                            {
                                break;
                            }
                        }
                    }

                    //! Compute Errors !//
                    if(DVnumber_Total <= *truthvtx_n)
                    {
                        ErrorsComputed++;
                        if(*truthvtx_n == 1) ErrorsComputed_OneDV++;
                        if(*truthvtx_n == 2) ErrorsComputed_TwoDVs++;

                        for(er=0; er<*truthvtx_n; er++)
                        {
                            if(!IndexUsed(usedErrorIndex, er))
                            {
                                // Error of DV_truth[k] and DV_reco
                                errorXYZ[errorCounter][0] = Error(displacedVertexArray[0], displacedVertexArray[1], displacedVertexArray[2], truthvtx_x[er], truthvtx_y[er], truthvtx_z[er]); 
                                errorXY[errorCounter][0] = Error(displacedVertexArray[0], displacedVertexArray[1], 0, truthvtx_x[er], truthvtx_y[er], 0); 
                                // Index of DV_truth used
                                errorXYZ[errorCounter][1] = er;
                                errorXY[errorCounter][1] = er;

                                errorCounter++;
                            }
                        }
                        min_ErrorXYZ = minimumArrayValueTwo(errorXYZ, errorCounter);
                        min_ErrorXY = minimumArrayValueTwo(errorXY, errorCounter);
                        // Minimum Error
                        minErrorXYZ[0] = min_ErrorXYZ[0];  
                        minErrorXY[0] = min_ErrorXY[0]; 
                        // Index of DV_truth used
                        minErrorXYZ[1] = min_ErrorXYZ[1]; 
                        minErrorXY[1] = min_ErrorXY[1];

                        usedErrorIndex[errorIndexCounter] = minErrorXYZ[1];
                        errorIndexCounter++;

                        error_XYZ->Fill(minErrorXYZ[0]); // Distance of calculated DV from truth DV (Error)
                        error_XY->Fill(minErrorXY[0]); // Distance of calculated DV from truth DV (Error)

                        if(*truthvtx_n == 1)
                        {
                            error_XYZ_OneDV->Fill(minErrorXYZ[0]); // Distance of calculated DV from truth DV (Error)
                            error_XY_OneDV->Fill(minErrorXY[0]); // Distance of calculated DV from truth DV (Error)
                        }
                        else if(*truthvtx_n == 2)
                        {
                            error_XYZ_TwoDVs->Fill(minErrorXYZ[0]); // Distance of calculated DV from truth DV (Error)
                            error_XY_TwoDVs->Fill(minErrorXY[0]); // Distance of calculated DV from truth DV (Error)
                        }

                        //! Calculate the number of DV_truth with match that respect the limits !//
                        // DV_true_matched that respect limit on sz space
                        if(minErrorXYZ[0] <= limitXYZ)
                        {
                            DV_true_with_match_SZ++;
                            if(*truthvtx_n==1) DV_true_with_match_SZ_OneDV++;
                            if(*truthvtx_n==2) DV_true_with_match_SZ_TwoDVs++;
                        } 

                        // DV_true_matched that respect limit on xy plane
                        if(minErrorXY[0] <= limitXY)
                        {
                            DV_true_with_match_XY++;
                            if(*truthvtx_n==1) DV_true_with_match_XY_OneDV++;
                            if(*truthvtx_n==2) DV_true_with_match_XY_TwoDVs++;
                        } 

                        //! Condition to calculate the DV_reco that respects or not the limits !//
                        // DV_reco that respect limit in sz space
                        if(minErrorXYZ[0] <= limitXYZ)
                        {
                            DV_reco_sz++;
                            if(*truthvtx_n==1) DV_reco_sz_OneDV++;
                            if(*truthvtx_n==2) DV_reco_sz_TwoDVs++;
                        }

                        // DV_reco that respect limit in xy plane
                        if(minErrorXY[0] <= limitXY)
                        {
                            DV_reco_xy++;
                            if(*truthvtx_n==1) DV_reco_xy_OneDV++;
                            if(*truthvtx_n==2) DV_reco_xy_TwoDVs++;
                        }

                        //! Matched Dvreco  - Repsect Both Limits!//
                        if(minErrorXY[0] <= limitXY && minErrorXYZ[0] <= limitXYZ)
                        {
                            // Canvas 5 - Minimum Distance Between DVreco and Beginning of Lines Used to Reconstruct It
                            DVreco_RecoLinesMinDistance_Matched->Fill(DVrecoLine_min); 

                            // Canvas 6 - Maximum Angle Between ODV and A_iAA_i, B_jBB_j vectors
                            Angle_ODV_AB_max_Matched->Fill(theta_ODVAB_max);

                            // Canvas 7 - Distance Between Tracks Used to Reconstruct the DVreco
                            Track_Distance_DVrecoMatch = leastDistance[0];
                            Track_Distance_Matched->Fill(Track_Distance_DVrecoMatch);

                            // Canvas 8 - Distance of Closest to DVreco Additional Track (If It Exists)
                            if(ClosestTrack_DV>0) ClosestTrackDV_Matched->Fill(ClosestTrack_DV);

                            // Canvas 9 - Difference of Distances of DVreco (R_DVreco) and Point of Reconstructing Tracks (Rmin) from IT
                            RDVreco_Rmin_Matched->Fill(R_DVreco-Rmin);

                            // Canvas 10 -Maximum Relative Angle Between DVreco and Given Points of Reconstructed Tracks
                            if(RelativeAngle[0]>-1000 && RelativeAngle[1]>-1000)
                            {
                                Relative_Angle_Matched->Fill(RelativeAngle[0]);
                                Relative_Angle_Matched->Fill(RelativeAngle[1]);   
                            }
                        }

                        //! No Match Dvreco  - Does not Repsect Both Limits!//
                        if(minErrorXY[0] > limitXY || minErrorXYZ[0] > limitXYZ)
                        {
                            // Canvas 5 - Minimum Distance Between DVreco and Beginning of Lines Used to Reconstruct It
                            DVreco_RecoLinesMinDistance_NotMatched->Fill(DVrecoLine_min); 

                            // Canvas 6 - Maximum Angle Between ODV and A_iAA_i, B_jBB_j vectors
                            Angle_ODV_AB_max_NotMatched->Fill(theta_ODVAB_max);

                            // Canvas 7 - Distance Between Tracks Used to Reconstruct the DVreco
                            Track_Distance_DVrecoNoMatch = leastDistance[0];
                            if(Track_Distance_DVrecoNoMatch>0) Track_Distance_NotMatched->Fill(Track_Distance_DVrecoNoMatch);

                            // Canvas 8 - Distance of Closest to DVreco Additional Track (If It Exists)
                            if(ClosestTrack_DV>0) ClosestTrackDV_NotMatched->Fill(ClosestTrack_DV);

                            // Canvas 9 - Difference of Distances of DVreco (R_DVreco) and Point of Reconstructing Tracks (Rmin) from IT
                            RDVreco_Rmin_NotMatched->Fill(R_DVreco-Rmin);

                            // Canvas 10 - Maximum Relative Angle Between DVreco and Given Points of Reconstructed Tracks
                            if(Relative_Angle[0]>-1000 && Relative_Angle[1]>-1000)
                            {
                                Relative_Angle_NotMatched->Fill(Relative_Angle[0]);
                                Relative_Angle_NotMatched->Fill(Relative_Angle[1]);   
                            }
                        }
                    }

                    // DV_reco that respect both limits
                    if(minErrorXYZ[0] <= limitXYZ && minErrorXY[0] <= limitXY) DVnumber_Close++;
                        
                    // DV_reco that do not respect at least one limit
                    if((minErrorXYZ[0] > limitXYZ || minErrorXY[0] > limitXY)) DVnumber_Far++;

                    // Distances Between DVs in Events with Two DVs //
                    if(*truthvtx_n == 2)
                    {
                        DistanceBetweenTwoDVs->Fill(Error(truthvtx_x[0], truthvtx_y[0], truthvtx_z[0], truthvtx_x[1], truthvtx_y[1], truthvtx_z[1]));
                    }
                }else break;

            } while(countLine <= *track_n-2);

            //! Histogram 2 !//
            // Number of DVs 
            DvNumber->Fill(*truthvtx_n);
            // One DV
            if(*truthvtx_n == 1) DvNumber_OneDV->Fill(*truthvtx_n);
            // Two DV
            if(*truthvtx_n == 2) DvNumber_TwoDVs->Fill(*truthvtx_n);

            if(DVnumber_Total!=0)
            {
                // Takes into consideration all the DVs
                clarity->Fill(DVnumber_Total);
                // One DV
                if(*truthvtx_n == 1) clarity_OneDV->Fill(DVnumber_Total);
                // Two DV
                if(*truthvtx_n == 2) clarity_TwoDVs->Fill(DVnumber_Total);
            }

            if(DVnumber_Close!=0)
            {
                // Takes into consideration the DVs that respect limits
                performance->Fill(DVnumber_Close);
                // One DV
                if(*truthvtx_n == 1) performance_OneDV->Fill(DVnumber_Close);
                // Two DVs
                if(*truthvtx_n == 2) performance_TwoDVs->Fill(DVnumber_Close);
            }

            if(DVnumber_Far!=0)
            {
                // Takes into consideration the DVs that do not respect the limits
                OffErrorPerformance->Fill(DVnumber_Far);
                // One DV
                if(*truthvtx_n == 1) OffErrorPerformance_OneDV->Fill(DVnumber_Far);
                // Two DV
                if(*truthvtx_n == 2) OffErrorPerformance_TwoDVs->Fill(DVnumber_Far);
            }

            // ! Histogram 3 !//
            RelativeNumber->Fill(*truthvtx_n - DVnumber_Total);
            if(*truthvtx_n == 1) RelativeNumber_OneDV->Fill(*truthvtx_n - DVnumber_Total);
            if(*truthvtx_n == 2) RelativeNumber_TwoDVs->Fill(*truthvtx_n - DVnumber_Total);

            // Event Counter
            event++;

            //! DV_truth counters !//
            // Total number of Dvs in all events
            DV_Total += *truthvtx_n;
            // Total number of DVs in events with one DV
            if(*truthvtx_n == 1) DVnumber_OneDV_Total += *truthvtx_n;
            // Total number of DVs in events with two DVs
            if(*truthvtx_n == 2) DVnumber_TwoDVs_Total += *truthvtx_n;
            // Total number of DVs in events with three DVs
            if(*truthvtx_n == 3) DVnumber_ThreeDVs_Total += *truthvtx_n;
            // Total number of DVs in events with four DVs
            if(*truthvtx_n == 4) DVnumber_FourDVs_Total += *truthvtx_n;

            //! DV_reco counters !//
            DV_reco_total += DVnumber_Total;
            if(*truthvtx_n == 1) DV_reco_OneDV_Total += DVnumber_Total;
            if(*truthvtx_n == 2) DV_reco_TwoDVs_Total += DVnumber_Total;
        }
    }

    //! Canvas 1 !//
    TCanvas *c1 = new TCanvas("c1", "DV Errors in XYZ Space and xy Plane", 800, 900);
    c1->Divide(2,3);

    gStyle->SetOptStat(1111111);

    c1->cd(1);
    error_XYZ->SetFillColor(kAzure+1);
    error_XYZ->SetMinimum(0);
    error_XYZ->Draw();
    gPad->SetLogx();
    gPad->Update();

    c1->cd(2);
    error_XY->SetFillColor(kRed);
    error_XY->SetMinimum(0);
    error_XY->Draw();
    gPad->SetLogx();
    gPad->Update();

    c1->cd(3);
    error_XYZ_OneDV->SetFillColor(kOrange+7);
    error_XYZ_OneDV->SetMinimum(0);
    error_XYZ_OneDV->Draw();
    gPad->SetLogx();
    gPad->Update();

    c1->cd(4);
    error_XY_OneDV->SetFillColor(kGreen);
    error_XY_OneDV->SetMinimum(0);
    error_XY_OneDV->Draw();
    gPad->SetLogx();
    gPad->Update();

    c1->cd(5);
    error_XYZ_TwoDVs->SetFillColor(kMagenta);
    error_XYZ_TwoDVs->SetMinimum(0);
    error_XYZ_TwoDVs->Draw();
    gPad->SetLogx();
    gPad->Update();

    c1->cd(6);
    error_XY_TwoDVs->SetFillColor(kYellow);
    error_XY_TwoDVs->SetMinimum(0);
    error_XY_TwoDVs->Draw();
    gPad->SetLogx();
    gPad->Update();

    c1->Print();

    //! Canvas 2 !//
    TCanvas *c2 = new TCanvas("c2", "Performance and Clarity", 1400, 900);
    c2->Divide(4,3);

    gStyle->SetOptStat(1111111);

    c2->cd(1);
    DvNumber->SetFillColor(kMagenta);
    DvNumber->SetMinimum(0);
    DvNumber->SetMaximum(3500);
    DvNumber->Draw();

    c2->cd(2);
    clarity->SetFillColor(kRed);
    clarity->SetMinimum(0);
    clarity->SetMaximum(3500);
    clarity->Draw();

    c2->cd(3);
    performance->SetFillColor(kAzure+1);
    performance->SetMinimum(0);
    performance->SetMaximum(3500);
    performance->Draw();

    c2->cd(4);
    OffErrorPerformance->SetFillColor(kGreen);
    OffErrorPerformance->SetMinimum(0);
    OffErrorPerformance->SetMaximum(3500);
    OffErrorPerformance->Draw();

    c2->cd(5);
    DvNumber_OneDV->SetFillColor(kMagenta);
    DvNumber_OneDV->SetMinimum(0);
    DvNumber_OneDV->SetMaximum(3500);
    DvNumber_OneDV->Draw();

    c2->cd(6);
    clarity_OneDV->SetFillColor(kRed);
    clarity_OneDV->SetMinimum(0);
    clarity_OneDV->SetMaximum(3500);
    clarity_OneDV->Draw();

    c2->cd(7);
    performance_OneDV->SetFillColor(kAzure+1);
    performance_OneDV->SetMinimum(0);
    performance_OneDV->SetMaximum(3500);
    performance_OneDV->Draw();

    c2->cd(8);
    OffErrorPerformance_OneDV->SetFillColor(kGreen);
    OffErrorPerformance_OneDV->SetMinimum(0);
    OffErrorPerformance_OneDV->SetMaximum(3500);
    OffErrorPerformance_OneDV->Draw();

    c2->cd(9);
    DvNumber_TwoDVs->SetFillColor(kMagenta);
    DvNumber_TwoDVs->SetMinimum(0);
    DvNumber_TwoDVs->SetMaximum(3500);
    DvNumber_TwoDVs->Draw();

    c2->cd(10);
    clarity_TwoDVs->SetFillColor(kRed);
    clarity_TwoDVs->SetMinimum(0);
    clarity_TwoDVs->SetMaximum(3500);
    clarity_TwoDVs->Draw();

    c2->cd(11);
    performance_TwoDVs->SetFillColor(kAzure+1);
    performance_TwoDVs->SetMinimum(0);
    performance_TwoDVs->SetMaximum(3500);
    performance_TwoDVs->Draw();

    c2->cd(12);
    OffErrorPerformance_TwoDVs->SetFillColor(kGreen);
    OffErrorPerformance_TwoDVs->SetMinimum(0);
    OffErrorPerformance_TwoDVs->SetMaximum(3500);
    OffErrorPerformance_TwoDVs->Draw();

    c2->Print();

    //! Canvas 3 !//
    TCanvas *c3 = new TCanvas("c3", "Relative Number of DV_reco with Respect to DV_truth", 1200, 300);
    c3->Divide(3,1);

    gStyle->SetOptStat(1111111);

    c3->cd(1);
    RelativeNumber->SetFillColor(kRed);
    RelativeNumber->SetMinimum(0);
    RelativeNumber->SetMaximum(3000);
    RelativeNumber->Draw();

    c3->cd(2);
    RelativeNumber_OneDV->SetFillColor(kAzure+1);
    RelativeNumber_OneDV->SetMinimum(0);
    RelativeNumber->SetMaximum(3000);
    RelativeNumber_OneDV->Draw();

    c3->cd(3);
    RelativeNumber_TwoDVs->SetFillColor(kGreen);
    RelativeNumber_TwoDVs->SetMinimum(0);
    RelativeNumber->SetMaximum(3000);
    RelativeNumber_TwoDVs->Draw();
    
    c3->Print();


    //! Canvas 4 !//
    TCanvas *c4 = new TCanvas("c4", "Distances Between DVs in Events with Two DVs", 400, 300);
    // c4->Divide(3,1);

    c4->cd(1);
    DistanceBetweenTwoDVs->SetFillColor(kGreen);
    DistanceBetweenTwoDVs->SetMinimum(0);
    DistanceBetweenTwoDVs->Draw();
    gPad->SetLogx();
    gPad->Update();

    c4->Print();

    //! Canvas 5 !//
    TCanvas *c5 = new TCanvas("c5", "Minimum Distance Between DVreco and Reconstructing Lines", 800, 300);
    c5->Divide(2,1);

    c5->cd(1);
    DVreco_RecoLinesMinDistance_Matched->SetFillColor(kAzure+1);
    DVreco_RecoLinesMinDistance_Matched->SetMinimum(0);
    DVreco_RecoLinesMinDistance_Matched->Draw();
    gPad->SetLogx();
    gPad->Update();

    c5->cd(2);
    DVreco_RecoLinesMinDistance_NotMatched->SetFillColor(kRed);
    DVreco_RecoLinesMinDistance_NotMatched->SetMinimum(0);
    DVreco_RecoLinesMinDistance_NotMatched->Draw();
    gPad->SetLogx();
    gPad->Update();

    c5->Print();

    //! Canvas 6 !//
    TCanvas *c6 = new TCanvas("c6", "Maximum Angle Between ODV and A_iAA_i, B_jBB_j vectors", 800, 300);
    c6->Divide(2,1);

    c6->cd(1);
    Angle_ODV_AB_max_Matched->SetFillColor(kAzure+1);
    Angle_ODV_AB_max_Matched->SetMinimum(0);
    Angle_ODV_AB_max_Matched->Draw();

    c6->cd(2);
    Angle_ODV_AB_max_NotMatched->SetFillColor(kRed);
    Angle_ODV_AB_max_NotMatched->SetMinimum(0);
    Angle_ODV_AB_max_NotMatched->Draw();

    c6->Print();

    //! Canvas 7 !//
    TCanvas *c7 = new TCanvas("c7", "Track Distance Between Tracks Used to Reconstruct the DVreco", 800, 300);
    c7->Divide(2,1);

    gStyle->SetOptStat(1111111);

    c7->cd(1);
    Track_Distance_Matched->SetFillColor(kAzure+1);
    Track_Distance_Matched->SetMinimum(0);
    Track_Distance_Matched->Draw();
    gPad->SetLogx();
    gPad->Update();

    c7->cd(2);
    Track_Distance_NotMatched->SetFillColor(kRed);
    Track_Distance_NotMatched->SetMinimum(0);
    Track_Distance_NotMatched->Draw();
    gPad->SetLogx();
    gPad->Update();

    c7->Print();

    //! Canvas 8 !//
    TCanvas *c8 = new TCanvas("c8", "Distance of Closest to DVreco Additional Track (If It Exists)", 800, 300);
    c8->Divide(2,1);

    gStyle->SetOptStat(1111111);

    c8->cd(1);
    ClosestTrackDV_Matched->SetFillColor(kAzure+1);
    ClosestTrackDV_Matched->SetMinimum(0);
    ClosestTrackDV_Matched->Draw();
    gPad->SetLogx();
    gPad->Update();

    c8->cd(2);
    ClosestTrackDV_NotMatched->SetFillColor(kRed);
    ClosestTrackDV_NotMatched->SetMinimum(0);
    ClosestTrackDV_NotMatched->Draw();
    gPad->SetLogx();
    gPad->Update();

    c8->Print();

    //! Canvas 9 !//
    TCanvas *c9 = new TCanvas("c9", "Difference of Distances of DVreco (R_DVreco) and Point of Reconstructing Tracks (Rmin) from IT", 800, 300);
    c9->Divide(2,1);

    gStyle->SetOptStat(1111111);

    c9->cd(1);
    RDVreco_Rmin_Matched->SetFillColor(kAzure+1);
    RDVreco_Rmin_Matched->SetMinimum(0);
    RDVreco_Rmin_Matched->Draw();

    c9->cd(2);
    RDVreco_Rmin_NotMatched->SetFillColor(kRed);
    RDVreco_Rmin_NotMatched->SetMinimum(0);
    RDVreco_Rmin_NotMatched->Draw();

    c9->Print();

    //! Canvas 10 !//
    TCanvas *c10 = new TCanvas("c10", "Maximum Relative Angle Between DVreco and Given Points of Reconstructed Tracks", 800, 300);
    c10->Divide(2,1);

    gStyle->SetOptStat(1111111);

    c10->cd(1);
    Relative_Angle_Matched->SetFillColor(kAzure+1);
    Relative_Angle_Matched->SetMinimum(0);
    Relative_Angle_Matched->Draw();
    gPad->SetLogx();
    gPad->Update();

    c10->cd(2);
    Relative_Angle_NotMatched->SetFillColor(kRed);
    Relative_Angle_NotMatched->SetMinimum(0);
    Relative_Angle_NotMatched->Draw();
    gPad->SetLogx();
    gPad->Update();

    c10->Print();

    //! Efficiency !//
    Eff_xy = 1.*DV_true_with_match_XY/DV_Total;
    Eff_xy_OneDV = 1.*DV_true_with_match_XY_OneDV/DVnumber_OneDV_Total;
    Eff_xy_TwoDVs = 1.*DV_true_with_match_XY_TwoDVs/DVnumber_TwoDVs_Total;
    Eff_sz = 1.*DV_true_with_match_SZ/DV_Total;
    Eff_sz_OneDV = 1.*DV_true_with_match_SZ_OneDV/DVnumber_OneDV_Total;
    Eff_sz_TwoDVs = 1.*DV_true_with_match_SZ_TwoDVs/DVnumber_TwoDVs_Total;

    //! Purity !//
    Pu_xy = 1.*DV_reco_xy/DV_reco_total;
    Pu_xy_OneDV = 1.*DV_reco_xy_OneDV/DV_reco_OneDV_Total;
    Pu_xy_TwoDVs = 1.*DV_reco_xy_TwoDVs/DV_reco_TwoDVs_Total;
    Pu_sz = 1.*DV_reco_sz/DV_reco_total;
    Pu_sz_OneDV = 1.*DV_reco_sz_OneDV/DV_reco_OneDV_Total;
    Pu_sz_TwoDVs = 1.*DV_reco_sz_TwoDVs/DV_reco_TwoDVs_Total;

    //! Accuracy !//
    accuracy_total = 1.*DV_reco_sz/ErrorsComputed;
    accuracy_OneDV = 1.*DV_reco_sz_OneDV/ErrorsComputed_OneDV;
    accuracy_TwoDVs = 1.*DV_reco_sz_TwoDVs/ErrorsComputed_TwoDVs;

    cout<<endl<<"Errors Computated: "<<ErrorsComputed<<endl;

    cout<<endl<<"Events: "<<event<<endl;
    cout<<"Total Number of Dvs: "<<DV_Total<<endl;
    cout<<"Total Number of Dvs in events with one DV: "<<DVnumber_OneDV_Total<<endl;
    cout<<"Total Number of Dvs in events with two DVs: "<<DVnumber_TwoDVs_Total<<endl;
    cout<<"Total Number of Dvs in events with three DVs: "<<DVnumber_ThreeDVs_Total<<endl;
    cout<<"Total Number of Dvs in events with four DVs: "<<DVnumber_FourDVs_Total<<endl;

    cout<<endl;
    for(int k = 0; k<50; k++) cout<<"~";
    cout<<endl;

    cout<<endl<<"Efficiency:"<<endl<<endl;
    cout<<"A_xy: "<<Eff_xy<<endl;
    cout<<"One DV - A_xy: "<<Eff_xy_OneDV<<endl;
    cout<<"Two DVs - A_xy: "<<Eff_xy_TwoDVs<<endl;
    cout<<endl<<"A_ρz: "<<Eff_sz<<endl;
    cout<<"One DV - A_ρz: "<<Eff_sz_OneDV<<endl;
    cout<<"Two DVs - A_ρz: "<<Eff_sz_TwoDVs<<endl;

    cout<<endl;
    for(int k = 0; k<50; k++) cout<<"~";
    cout<<endl;

    cout<<endl<<"Purity:"<<endl<<endl;
    cout<<"Pu_xy: "<<Pu_xy<<endl;
    cout<<"One DV - Pu_xy: "<<Pu_xy_OneDV<<endl;
    cout<<"Two DVs - Pu_xy: "<<Pu_xy_TwoDVs<<endl;
    cout<<endl<<"Pu_ρz: "<<Pu_sz<<endl;
    cout<<"One DV - Pu_ρz: "<<Pu_sz_OneDV<<endl;
    cout<<"Two DVs - Pu_ρz: "<<Pu_sz_TwoDVs<<endl;

    cout<<endl;
    for(int k = 0; k<50; k++) cout<<"~";
    cout<<endl;

    cout<<endl<<"Accuraty:"<<endl<<endl;
    cout<<"Total: "<<accuracy_total<<endl;
    cout<<"Event with One DV: "<<accuracy_OneDV<<endl;
    cout<<"Event with Two DVs: "<<accuracy_TwoDVs<<endl;

    cout<<endl;
    for(int k = 0; k<50; k++) cout<<"~";
    cout<<endl;

    // Print time needed for the program to complete
    printf("\nTime taken: %.2fs\n\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

    infile->Close();
    }