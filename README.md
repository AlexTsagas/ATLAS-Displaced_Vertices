# Project - In Search for Displaced Vertices

The purpose of the project is to find long lived particles/Displaved Vertices (DVs) which are depicted as points where
trajectories of its products converge. To specify the trajectories we are given two points that belong to them. The DVs 
are calculated as follows:

First, we proccess events with one DV.

• The trajectories that leads us to it need to be "close together". So we choose the closest trajectories of the event
to calculate the DV. The code produces histogram of the closest distances of every event.
  
• Then, we define the DV as the midpoint of the distance vector of the closest trajectories. Also, to narrow the false 
positives an if statement is used to exclude some points. The requirement for a point P calculated from line_i and line_j
is that vectors PA_i, PB_i, and PA_j, PB_j to created angle less or equal than π/2, where A and B are the begin and end
points of each line. The distance between the real DV and the recommended is presented as absolute error in a histogram.
  
• Also, the closest distance to the third trajectory of every event is calculated and presented in a histogram to help
  find an if statement to apply in order to procces events with more than one DV.
  


