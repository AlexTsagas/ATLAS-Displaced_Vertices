# Project - In Search for Displaced Vertices

The purpose of the project is to find long lived particles/Displaved Vertices (DVs) which are depicted as points where
trajectories of its products converge. To specify the trajectories we are given two points that belong to them. The DVs 
are calculated as follows:

• The trajectories that leads us to them need to be "close together". So we choose the closest trajectories of an event
  to calculate the DV. The code produces histogramm of the closest distances of every event.
  
• Then, we use define the DV as the midpoint of the distance vector of the closest trajectories of the event. The distance
  between the real DV and the recommended is presented as absolute error in a histogram.
  
• Also, there is calculated the closest distance to the third trajectory of every event and presented in a histogram to help
  find an if statement to apply in order to procces events with more than one DV.
  


