# Functional Thin Films on Surfaces

The attached code provides a partial implementation of the method for integrating the thin films equation on curved sufraces as described in
https://dl.acm.org/citation.cfm?id=2786793

![Alt text](images/ftf_code_img_hr.png?raw=true "Teaser")

# Contents

folders structure:

code/  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;               contains the code for our method\
code/experiments/ &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 		the saved simulations are dropped here\
data/ 					&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;       geometry files\
external/ 			&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 	    code of other methods that we used\

files:

gen_* 			&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 		      scripts generating specific simulations\
LINFTF 			&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 		      this class implements the simulation part (integration in time)\
MESH 				&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 	        this class is responsible for the geometry (integration in space)\
MESH_VIS 		&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 		      helper class for visualizing functions and vector fields in MATLAB\
MESH_READER 	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 		    helper class to read geometry files\
quadratic_* 		&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 	    scripts for the initial condition computations
