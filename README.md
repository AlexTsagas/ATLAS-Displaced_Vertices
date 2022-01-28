# Project - In Search for Displaced Vertices

## Goal of Project

The purpose of the project is to find long lived particles/Displaved Vertices (DVs) which are depicted as points where
trajectories of its products (stable particles) converge. To specify the trajectories we are given two points that belong to them. The DVs 
are calculated as follows:

## Even Proccessing

### Events with One DV

  • The trajectories that lead us to it need to be "close together". So we choose the closest trajectories of the event
    to calculate the DV.

  • Then, we define the DV as the midpoint of the distance vector of the closest trajectories. Also, to narrow the false 
    positives an if statement is used to exclude some points. The requirement for a point P calculated from line_i and line_j
    is vectors PA_i, PB_i, and PA_j, PB_j to created angle less or equal than π/2, where A and B are the begin and end
    points of each line. The distance between the real DV and the recommended is presented as error in histogram.

  • Also, the closest distance to the third trajectory of every event is calculated to help find an if statement to apply
    in order to procces events with more than one DV. That if statement would count trajectories that belong to the same DV
    but are not used to construct it.
  
### Events with Multiple DVs

  • The first step is to specify the number of DVs that exists in every event. So the goal is to find the DV_reco(mmended) with the 
    most accurate approximation possible. When the total number of DV_reco is calculated the code produces histograms that compare it
    with the total number of real DVs (DV_true).

  • For every DV_reco calculated two trajectories are used, the one that are closest together and that respect if statement of relative
    angles mentioned in the previous paragraph. Then, an if statement is used to exclude more trajectories that might belong to the same
    DV_reco, using an upper bound for the distance between them and the DV_reco. If another DV_reco exists in the same event we exclude the 
    previous trajectories used, and we search the other with the remaining ones. To continue searching for DVs the algorith checks if there
    are sufficient number of trajectories left to construct a DV (more than 2) and the minimum distance of the closest trajectories remaining
    unused (checks for minimum distances lower than an upper bound).

  • The error between the DV_truth and DV_reco is defined as the distance between them. For the first calculated DV_reco of an event 
    every possible error with DV_truths is computed and, then, the minimum is picked to represent the value. Then, the DV_truth used to 
    calculated the error of the first DV_reco is excluded. The same proccess is repeated to compute the other DV_reco's erros (if they exist).
