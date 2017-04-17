# SelfDrivingCar2-UKF
Sigma-Point filtering for Moving Object Tracking using Lidar and Radar processed neasurements

# AN EKF filter for position tracking using LIDAR and RADAR processed measurements


## Compiling and Executing
Compling and running should be straightforward (I have included the test files in the build directory):
```
cd build
cmake ..
make
./UnscentedKF sample-laser-radar-measurement-data-1.txt output1.txt
./UnscentedKF sample-laser-radar-measurement-data-2.txt output2.txt
```
## How the filter works
The filter processes the LIDAR and RADAR measurements interchangeably with the same update function. The prediction constructs a block-diagonal matrix using the prior covariance and the covariance of the measurement model as diagonal blocks in the order given. This yields a joint covariance matrix from which we can sample the sigma-points along the directions of the Cholesky factor columns (strangely, this works better than the orthogonal direction obtaind with the SVD). Thus, the state posterior is  the marginal of the state variable in the joint, which can be obtained by the standard UKF weighted sums for the first 5 components of the transformed sigma-points. 

Now, the Lidar measurement update is linear and can be done analytically, but since this is an exercise, I decided to use the usual sigma-point summation in this as well. Similarly, the Radar update is done using the typical UKF formulas. I attempted to use _iterative sigma-point_ updates a Levenberg-Marquardt adaptive regularization factor (which is similar to the method [here](http://robotics.usc.edu/publications/media/uploads/pubs/500.pdf) but NOT the same, as I had to modify the formulas to account for the LM regularization factor), in a sense, is a way to adapt the locations of the sigma-points along the columns of the Cholesky factor of the joint distribution, but, the iteration proved somewhat unstable (i.e., the descent was erratic). Although I can't exclude the possibility of a bug, it seems that iterative approaches suffer from the same weaknesses that standard UKFs do, which is initial tuning.

## Tuning

The sigma-point method is a nice trick, but has certain disadavantages and they all pertain to **how to pick the sigma points along the axes of the covariance ellipsoid**. In this code, the Cholesky factor columns are used, which are clearly not orthogonal directions as implied/stated in the sigma-point papers. Surprisingly, the scaled, non-orthogonal directions in Cholesky factor gave much better results than the SVD-based sigma points. Of course, the latter may require different tuning of the original parameters using α - squared, β, and κ instead of th λ suggested in the lectures. In any case, it is a dark path to follow with many-many trial-and-error actions. 

## A Few Implementation Comments

The update functions are static, so they contain duplicate constant designations of  constants that are already declared in the UKF class. They were made static because I wanted to run them in various tests separately from the UKF class. Thus, the weights have to be computed inside these  functions. There is a single multiplication constant-time additional burden inside a loop that needs to run anyway in order to transform the sigma-points and it is basically negligible. As explained earlier, I initially intended to explore the posibility of configuring the UKF with its standard parameters α, β, κ, but I quickly realized that λ=3-dim_x is a very good option; since the updates, being static functions, receive α, β, κ in their parameters, I simply modified the code inside the updates to work-out α in order to achive λ = dim_x. All in all, the code was adjusted for plenty of experimentation, but it doesn't replicate trails of computations.   
