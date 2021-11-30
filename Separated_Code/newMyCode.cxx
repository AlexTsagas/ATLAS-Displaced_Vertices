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
// product is defined as (a, b x c).
double tripleProduct(double *a, double *b, double *c)
{
    double *crossp = crossProduct(b, c);
    double cross[3] = {crossp[0], crossp[1], crossp[2]};

    double triple = dotProduct(a, cross);

    return triple;
}


// Input an array with three elements and output the norm sqrt{\vb{a}\cdot\vb{a}}
double norm(double *a)
{
    double norm;
    double dotp = dotProduct(a, a);

    norm = sqrt(dotp);

    return norm;
}


// Relative Vector that points from a to b
double *relativeVector(double *a, double *b)
{
    static double relative[3];
    for(int i=0; i<3; i++)
    {
        relative[i] = b[i] - a[i];
    }

    return relative;
}

// Input two vectors a, b with three elemets each. Output the unit vector that points
// from a to b.
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


//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //


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


//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //


void NewMyCode()
{
    // line1
    double a[3] = {0, 0, 0};
    double b[3] = {1, 0, 0};

    // line2
    double aa[3] = {0, 0, 1};
    double bb[3] = {0, 1, 1};

    displacedVertex(a, b, aa, bb);

    // Distance
    double d = distance(a, b, aa, bb);
    cout<<"Distance of two lines is: "<<d<<endl;

    // Displaced Vertex 
    double *DV = displacedVertex(a, b, aa, bb);
    cout<<"The coordinates of displaced vertex are: ";
    printCoordinates(DV);

    //  ~~~~~~~~ // ~~~~~~~~ // ~~~~~~~~ // ~~~~~~~~ // ~~~~~~~~ // ~~~~~~~~ // ~~~~~~~~ //

    // Graph to Visualize Trajectories and the Displayced Vertex

    // Make Points For Plot

    // number of points
    int n = 150; 
    // t in (begin, end) 
    int begin = -30, end = 30;
    // Step
    double dt = (end - begin)/(1.*n);
    // Argument value
    double t;

    // line1 points' coordinates
    double line1_x[n], line1_y[n], line1_z[n];
    // line2 points' coordinates
    double line2_x[n], line2_y[n], line2_z[n]; 
    // dv's points' coordinates
    double dv_x[1], dv_y[1], dv_z[1];

    // Assign values to coordinates of lines
    for(int i=0; i<n; i++)
    {
        t = begin + dt*i;

        double *line1 = lineEquation(a, b, t);
        line1_x[i] = line1[0];
        line1_y[i] = line1[1];
        line1_z[i] = line1[2];

        double *line2 = lineEquation(aa, bb, t);
        line2_x[i] = line2[0];
        line2_y[i] = line2[1];
        line2_z[i] = line2[2];
    }

    // Assining values to coordinates of dv
    double *dv = displacedVertex(a, b, aa, bb);
    dv_x[0] = dv[0];
    dv_y[0] = dv[1];
    dv_z[0] = dv[2];

    // Points for the Plot
    double x[2*n+1], y[2*n+1], z[2*n+1];

    for (int i = 0; i < n; i++)
    {
        x[i] = line1_x[i];
        x[i+n] = line2_x[i];

        y[i] = line1_y[i];
        y[i+n] = line2_y[i];

        z[i] = line1_z[i];
        z[i+n] = line2_z[i];
    }
    
    x[2*n] = dv_x[0];
    y[2*n] = dv_y[0];
    z[2*n] = dv_z[0];

    // Make the canvas and specify title and width, height
    TCanvas *c = new TCanvas("c", "Displaced Vertex;x;y;z");

    TGraph2D *gr = new TGraph2D(2*n+1, x, y, z);
    gr->SetTitle("Displaced Vertex and Particles' Trajectories;x;y;z");

    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.5);
    gr->SetMarkerColor(kBlack);

    gr->Draw("P");
}   
