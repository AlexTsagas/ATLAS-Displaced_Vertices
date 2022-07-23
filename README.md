# In Search for Displaced Vertices

## How to Run the Algorithm

For the production of the algorithm root framework has been used. More information can be found in the root's website:
[root.cern](https://root.cern).

## Goal of Project

The purpose of the project is to find long lived particles/Displaced Vertices (DVs) which are depicted as points where
trajectories of its products (stable particles) converge. The data used for the trajectories points are simulated open
data from ATLAS Experiment. The DVs are calculated as follows:

## Event Processing

Specifically, for every trajectory two points are given, the beginning and the end point of the trajectory. Considering
that the trajectories are linear, two points are enough to specify the whole trajectory.

### Events with One DV

  * The trajectories that lead us to it need to be "close together". So, we choose the closest trajectories of the event
    to calculate the DV. Specifically, the distance vector from line_i to line_j is calculated and only values less than
    DVcut are allowed for its length. 

  * Then, we define the DV as the midpoint of the distance vector of the closest trajectories. Also, in order to narrow the
    false positives an "if statement" is used to exclude some points. The requirement for a point P to be calculated from line_i 
    and line_j is the vectors: PA_i, PB_i, and PA_j, PB_j to create an angle less or equal than π/2, where A and B are the begin and
    end points of each line. Also, the DV must be closer to the detector's center than any point from the used trajectories.
    The distance between the real DV and the reconstructed is presented as error in histogram.

  * Furthermore, to process events with more than one DVs the closest distance of DV to the third trajectory of every event (if it exists)
    is calculated. The histogram produced by the results offers information about the cases that must be considered. To take them into 
    account an "if statement" should be used that would count trajectories which belong to the same DV but are not used to construct it.
  
### Events with Multiple DVs

  * The first step is to specify the number of DVs that exists in every event. Thus, the goal is to find the DV_reco(nstructed) with the 
    most accurate approximation possible. When the total number of DV_reco is calculated, the code produces histograms that compare it
    with the total number of real DVs (DV_true).

  * For every DV_reco calculated, two trajectories are used:
      * those that are closest together and construct the DV_reco, and
      * those which respect the "if statement" of relative angles and distance from detector's center mentioned in the previous paragraph. 
    Then, an if statement is used to exclude more trajectories that might belong to the same DV_reco, using an upper bound (TrajectoryCut)
    for the distance between them and the DV_reco. If another DV_reco exists in the same event, we exclude the previous trajectories used 
    and we search the other with the remaining ones. To continue searching for DVs the algorithm checks:
      * if there is a sufficient number of trajectories left to construct a DV (more than 2), and
      * the minimum distance of the closest trajectories remaining unused (check for minimum distances lower than the upper bound DVcut).

  * The error between the DV_true and DV_reco is defined as the distance between them. For the first calculated DV_reco of an event 
    every possible error with DV_trues is computed and, then, the minimum is picked to represent the value. Then, the DV_true used to 
    calculate the error of the first DV_reco is excluded. The same process is repeated to compute the other DV_reco's errors (if they exist).
