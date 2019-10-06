# Functional Thin Films on Surfaces

The attached code provides a partial implementation of the method for integrating the thin films equation on curved sufraces as described in
https://dl.acm.org/citation.cfm?id=2786793

![Alt text](images/ftf_code_img_hr.png?raw=true "Teaser")

# Contents

folders structure:

code/  tbsp; tbsp; tbsp; tbsp;               contains the code for our method\
code/experiments/ 		the saved simulations are dropped here\
data/ 					      geometry files\
external/ 				    code of other methods that we used\

files:

gen_* 					      scripts generating specific simulations\
LINFTF 					      this class implements the simulation part (integration in time)\
MESH 					        this class is responsible for the geometry (integration in space)\
MESH_VIS 				      helper class for visualizing functions and vector fields in MATLAB\
MESH_READER 			    helper class to read geometry files\
quadratic_* 			    scripts for the initial condition computations
