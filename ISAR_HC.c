#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>


/*-----------------------------------------------------------------------------
Generates perfect samples of Hard-core Gibbs process using grid based importance sampling proposed in Moka et al (2017) (arXiv version).
Where the Hard-core process is defined on [0,1]^2 with radius r/lambda^eta and absolutely continuous with respect to a lambda-hamogeneous PPP
-------------------------------------------------------------------------------*/
/*----------------  Steps to run this program -----------------------

$ gcc ISAR_HC_Euclidean_v3.c -o isar -Wall -lgsl -lgslcblas -lm
$ ./isar
---------------------------------------------------------------------*/


/*----------------------------------------------------------------------
 Possible improvements:
	   1. Don't normalize the distribution of M. Instead multiply the uniform random variable with the total mass
	   2. Use linked list to store blocking cells instead of grid array which is taking a lot of memory
	   3. Use binary tree for checking overlap.
	   4. Truncate the distribution of M using a decent upper bound on the close packing density
------------------------------------------------------------------------*/


void currentTime(void)
{
    time_t current_time;
    char* c_time_string;

    /* Obtain current time. */
    current_time = time(NULL);

    if (current_time == ((time_t)-1))
    {
        (void) printf("Failure to obtain the current time.\n");
        exit(EXIT_FAILURE);
    }

    /* Convert to local time format. */
    c_time_string = ctime(&current_time);

    if (c_time_string == NULL)
    {
        (void) fprintf(stderr, "Failure to convert the current time.\n");
        exit(EXIT_FAILURE);
    }

    /* Print to stdout. ctime() has already added a terminating newline character. */
    (void) printf("Current time is %s", c_time_string);
}
// Returns Euclidean distance between the points x and y
double euclideanDistance(double *x, double *y) {
	return sqrt((x[0] - y[0])*(x[0] - y[0]) + (x[1] - y[1])*(x[1] - y[1]));
}

// Returns distance between x and y under the assumption the underlying cube is a torus
double torusDistance(double *x, double *y) {

        double min_dist, t[2];
        min_dist = euclideanDistance(x, y);

        t[0] = x[0] - 1;
        t[1] = x[1] + 1;
        min_dist = fmin(min_dist, euclideanDistance(t, y));

        t[0] = x[0];
        t[1] = x[1] + 1;
        min_dist = fmin(min_dist, euclideanDistance(t, y));

        t[0] = x[0] + 1;
        t[1] = x[1] + 1;
        min_dist = fmin(min_dist, euclideanDistance(t, y));

        t[0] = x[0] - 1;
        t[1] = x[1];
        min_dist = fmin(min_dist, euclideanDistance(t, y));

        t[0] = x[0] + 1;
        t[1] = x[1];
        min_dist = fmin(min_dist, euclideanDistance(t, y));

        t[0] = x[0] - 1;
        t[1] = x[1] - 1;
        min_dist = fmin(min_dist, euclideanDistance(t, y));

        t[0] = x[0];
        t[1] = x[1] - 1;
        min_dist = fmin(min_dist, euclideanDistance(t, y));

        t[0] = x[0] + 1;
        t[1] = x[1] - 1;
        min_dist = fmin(min_dist, euclideanDistance(t, y));

        return(min_dist);
}


// Definition and functions related to list of points
typedef struct node node_t;

struct node {
        double pnt[2];
        struct node *next;
};

// Function for creating a node
node_t * createNode(){
        node_t *temp;
        temp = malloc(sizeof(node_t));
        temp->next = NULL;
        return temp;
}

// Function for freeing the memory
void freeMem(node_t *head) {
        node_t *temp;
        while (head != NULL) {
                temp = head->next;
                free(head);
                head = temp;
        }
}

// Function for adding a point to the list of existing points
node_t *addNode(node_t *head, double *a) {
	node_t *temp;
	temp = createNode();
	temp->pnt[0] = a[0];
	temp->pnt[1] = a[1];
	temp->next = head;
	head = temp;
	//printf("Point check in addNode fn (%lf, %lf)\n", head->pnt[0]*57, head->pnt[1]*57);
	return(head);
}

// Returns true if a overlaps with any of the existing spheres
bool overlap(node_t *head, double *a, double radius) {
	node_t *temp;
	temp = head;
	while (temp != NULL) {
		if (torusDistance(temp->pnt, a) <= 2*radius) {
			//printf("There is a overlap\n");
			return(true);
		}
		temp = temp->next;
	}
	return(false);
}


//Function for checking non-overlap property of the generated sample
bool OverlapOfSample(node_t *ptr, double radius) {
	if (ptr == NULL) {
		return(false);
	}
	node_t *temp;
	while (ptr != NULL) {
		temp = ptr->next;
		while (temp != NULL) {
			if (torusDistance(ptr->pnt, temp->pnt) < 2*radius) {
				return(true);
			}
			temp = temp->next;
		}
		ptr = ptr->next;
	}
	return(false);
}


// Function for printing the list of points in the perfect sample
void printList(node_t *head) {
	if (head == NULL) {
		printf("The list is empty\n");
		return;
	}

	node_t *temp;
	int i;

	temp = head;
	i = 0;
	while (temp != NULL) {
		printf("(%lf, %lf)\n", temp->pnt[0], temp->pnt[1]);
		i++;
		temp = temp->next;
	}
	printf("Total number of points: %d\n", i);
	return;
}

// Function for counting the number points in the list
int countPoints(node_t *ptr) {
	int i = 0;
	while (ptr != NULL) {
		ptr = ptr->next;
		i++;
	}
	return(i);
}




// The main program to generate a perfect sample of the target Hard-core process
int  main() {

	// Parameters of the HC process
  	double lambda = 50; // Intensity of the PPP
	double eta = 0.25; // Parameter associated with regime
	double r = 1; // r as in the radius r/lambda^eta


	// Variables
	double rad, vol, sqrt_2; // rad is actual radius r/lambda^eta and vol is the volume of a sphere with radius rad. See the grid method in the paper.
  	int  grid_size, sigma_len, Itot; // Number of non-zero sigma values
	bool stop;


	// One time computing of sqrt of 2.
	sqrt_2 = sqrt(2);

	// Variable assignments
	Itot = 1000;
	rad = r/((double) pow(lambda, eta));
	vol = M_PI*rad*rad;
	sigma_len =  fmin(floor( 1/((double) vol) + 2 ), ceil(M_PI/((double) 2*sqrt(3)*vol)));
	grid_size = ceil(2*sqrt_2/((double)rad)); // This assignment guarantees that every generated sphere covered by blocking cells
	printf("  Intensity = %lf \n  eta = %lf \n  r = %lf\n  Radius = %lf \n  Sigma_len: %d\n  Grid_size: %d\n  Itot: %d\n", lambda, eta, r, rad, sigma_len, grid_size, Itot);


	// Setting up random number gererator
	const gsl_rng_type *rng_type;
	gsl_rng *rng;
	gsl_rng_env_setup();
	rng_type = gsl_rng_default;
	rng = gsl_rng_alloc (rng_type);
	gsl_rng_set(rng, 0);

	// Exit loop to handle a sphere overlapping with itself
	if (2*rad > 1) {
		printf("Sphere radius r/lambda^eta should be less than 1/2 to avoid any sphere overlapping with itself\n");
		printf("Exiting the program ....\n");
		exit(0);
		}

	// Defining other variables required for IS method
	double sigma[sigma_len], dist_M[sigma_len], cum_dist_M[sigma_len],  epsilon, block_area, temp_var, npoints; // sigma stores the corresponding sigma values.
	int i, j, ii, jj, k, n, M, iter, nonblock, x_nonblock[grid_size];
	bool temp_bool, grid[grid_size][grid_size]; // This array takes a lot of memory. We should find a better way to write this.


	// Initializing the values of sigma to 1
	for (n = 0; n < sigma_len; n++) {
		sigma[n] = 1.0;
		dist_M[n] = (double) gsl_ran_poisson_pdf(n, lambda);
	}

	// Compute both sigma and dist_M
	temp_var = dist_M[0] + dist_M[1];
	for (n = 2; n < sigma_len; n++) {

		sigma[n] = sigma[n-1]*(1 - (n-1)*M_PI*pow(rad,2));

/*		if (sigma[n] < 0) { // This loop can be removed after debugging
			printf("Sigma is negative. Exiting the program ....\n");
			exit(0);
		} */

		dist_M[n] = dist_M[n]*sigma[n];
		temp_var = temp_var + dist_M[n];

		//printf("\n--------n : %u -------------- \n dist_M: %0.100lf", dist_M[n]);
	}

	// Normalization of dist_M comuting cummulative dist_M
	dist_M[0] = dist_M[0]/((double) temp_var);
	cum_dist_M[0] = dist_M[0];

  	for (n = 1; n < sigma_len; n++) {
		dist_M[n] = dist_M[n]/((double) temp_var);
		cum_dist_M[n] = cum_dist_M[n-1] + dist_M[n];
	}


/*	for (n = 0; n < sigma_len; n++) {
		printf("cum_dist_M[%u] = %lf \n", n, cum_dist_M[n]);
	}
*/

	// Initializing pointer to the list of points of HC Model and define some other variables
  	node_t *head;
	double temp_pt[2], Exp_pts_generated;
	// After generating a cicle with the center in a cell, say (x,y),
	// we check only cells ranging from i = x - range_check, ..., x +  range_check + 1, and j = y - range_check, ... y + range_check
	int range_check, corner_x, corner_y;
    int count_pts_generated;


	range_check = ceil(2*rad*grid_size);
	epsilon = 1/((double) grid_size); // Eadge length of each cell

	// Loop for generating a perfect sample starts here
	currentTime();
    clock_t start_time = clock();
	npoints = 0;
    Exp_pts_generated = 0;
	for (iter = 0; iter < Itot; iter++) {
        count_pts_generated = 0;
		stop = false;
		while (!stop) { // While loop 1

			// Generate M from dist_M
			temp_var = (double) gsl_rng_uniform_pos(rng);
			M = 0;
			while (1) {
				if (temp_var <= cum_dist_M[M]) {
					break;
				}
				M++;
			}

			//printf("M = %d\n", M);
			// Empty configuration is immediately accepted
		  	if (M == 0) {
	      		head = NULL;
				break;
	  		}

			// First point is generated
		    head = createNode();
		    head->pnt[0] = (double) gsl_rng_uniform_pos(rng);
		    head->pnt[1] = (double) gsl_rng_uniform_pos(rng);
		    head->next = NULL;
            count_pts_generated++;

			// Terminate the loop if only one point is needed
			if (M == 1) {
				break;
			}

			// Clearing the grid
			nonblock = grid_size*grid_size; // No of non-blocking cubes available (all of them at the beginning)
			for (i = 0; i < grid_size; i++) {
				x_nonblock[i] = grid_size; // No of non-blocking cubes avaiable in the i^th row (all of them at the beginning)
				for (j = 0; j < grid_size; j++) {
					grid[i][j] = false;

				}
			}

			// Generate the remaining points
			block_area = 0.0;

			for (n = 2; n <= M; n++) { // for loop 1 . Recursively generates M points
				corner_x = floor(((double) grid_size)*head->pnt[0]); // Cell values corresponds to the last generated circle
				corner_y = floor(((double) grid_size)*head->pnt[1]);
				//printf("Corners : (%d, %d) and range_check : %d\n", corner_x, corner_y, range_check);
				//printf("corner_x - range_check: %d\n", corner_x - range_check);
				//printf("corner_x + range_check: %d\n", corner_x + range_check);
				temp_var = 0.0;
				for (ii = corner_x - range_check; ii <= corner_x + range_check + 1; ii++) {
					//printf("column %d\n", ii);
					for (jj = corner_y - range_check; jj <= corner_y + range_check + 1; jj++) {
						//printf("row %d\n", jj);

						// On a torus, cells should be wrapped around
						if (ii < 0) {
							i = grid_size + ii;
						} else if (ii >= grid_size) {
							i = ii - grid_size;
						} else {
							i = ii;
						}

						if (jj < 0) {
							j = grid_size + jj;
						} else if (jj >= grid_size) {
							j = jj - grid_size;
						} else {
							j = jj;
						}
						//printf("Cell ckecking (%d, %d)\n", i,j);
						if (grid[i][j] == false) { // Check only when the cell is not blocked.

							temp_pt[0] = (double) (i + 0.5)*epsilon;
							temp_pt[1] = (double) (j + 0.5)*epsilon;
							temp_bool = false; // default the cell is assumed to be within 2rad from the center

							// If mid point of a cell is at least a distance of 2*rad - epsilon/2 away, the cell should be marked non-blocked.
							// If mid point of a cell is within a distance of 2*rad - sqrt(2)*epsilon/2, the entire cell should be within the sphere of radius 2r.
							//printf("torusDistance(temp_pt, head->pnt) %lf\n", torusDistance(temp_pt, head->pnt));
							if (torusDistance(temp_pt, head->pnt) > 2*rad - 0.5*epsilon) {
								//printf("Condition 1 holds\n");
								temp_bool = true;
							} else if (torusDistance(temp_pt, head->pnt) > 2*rad - 0.5*sqrt_2*epsilon ) {
							  // For the other case, check all the four corners
								temp_pt[0] = ((double) i)*epsilon;
								temp_pt[1] = ((double) j)*epsilon;
								if (torusDistance(temp_pt, head->pnt) > 2*rad ) {
									temp_bool = true;
								}

								temp_pt[1] = ((double) (j + 1))*epsilon;
								if (!temp_bool && torusDistance(temp_pt, head->pnt) > 2*rad ) {
									temp_bool = true;
								}

								temp_pt[0] = ((double) (i + 1))*epsilon;
								if (!temp_bool && torusDistance(temp_pt, head->pnt) > 2*rad ) {
									temp_bool = true;
								}

								temp_pt[1] = ((double) j)*epsilon;
								if (!temp_bool && torusDistance(temp_pt, head->pnt) > 2*rad ) {
									temp_bool = true;
								}
							}

							// Updating values related to the grid
							if (!temp_bool) {
								//printf("Blocked cell (%d, %d)\n", i, j);
								temp_var = temp_var + epsilon*epsilon;
								grid[i][j] = true;
						       	nonblock--;
		                		x_nonblock[i]--;
								block_area = block_area + epsilon*epsilon;

				                // Sanity checks
				              	/*if (nonblock < 0 || x_nonblock[i] < 0 ) {
					                printf("Volume of each cell: %lf\n ", epsilon*epsilon);
									printf("Blocked area: %0.100lf, Nonblock: %d, x_nonblock: %d, row: %d\n", block_area, nonblock, x_nonblock[i],i);
			                      	exit(0);
								}*/
							}
						}
					}
				}
				//printf("Area blocked by %d th circle is %0.10lf\n", n, temp_var);
				// Generate new circle on the non-blocking cells
				temp_var = ((double) gsl_rng_uniform_pos(rng))*nonblock;

				i = 0;
				k = x_nonblock[0];
				while (1) {
					if (temp_var <= k) {
						break;
					}
					i++;
					k = k + x_nonblock[i];
				}
				temp_var = ((double) gsl_rng_uniform_pos(rng))*x_nonblock[i];

				j = 0;
				k = 1 - grid[i][0];
				while (1) {
					if (temp_var <= k) {
						break;
					}
					j++;
					k = k + 1 - grid[i][j];
				}
				//printf("(%d, %d) th cell is %d\n", i, j, grid[i][j]);
				temp_pt[0] = ((double) i + (double) gsl_rng_uniform_pos(rng))*epsilon;
				temp_pt[1] = ((double) j + (double) gsl_rng_uniform_pos(rng))*epsilon;
                count_pts_generated++;
				//printf("New point generated (%lf, %lf)\n", temp_pt[0]*grid_size, temp_pt[1]*grid_size);

				if (((double) gsl_rng_uniform_pos(rng))*((double) (1 - (n-1)*M_PI*pow(rad,2))) > (1 - block_area) || overlap(head, temp_pt, rad)) {
					//printf("Not accepted\n");
					freeMem(head);
					break;
				} else {
					//printf("Accepted\n");
					head = addNode(head, temp_pt);
					//printf("New point cross check (%lf, %lf)\n", head->pnt[0]*grid_size,  head->pnt[1]*grid_size);
					if (n == M) {
						stop = true;
					}
				}
			/*	// Sanity check
				if ( (1 - block_area) > ((double) (1 - (n-1)*M_PI*pow(rad,2)))) {
					printf("Non-blocking area : %lf and the sigma: %lf\n", (1 - block_area), (double) (1 - (n-1)*M_PI*pow(rad,2)));
					printf("Error: when n = %u, the non-blocking area is more than sigma value....\n", n);
					exit(0);
				}*/
			} // End of the for loop 1

		} // End of the while loop 1


		// Verifying the generated points
		/*if (OverlapOfSample(head, rad)){
			printf("Not a valid Hard Core configuration\n");
		}*/
		// Printing the list of point in a perfect sample of the HC model
		// printList(head);

		npoints = (iter*npoints)/((double) (iter + 1)) + countPoints(head)/((double) (iter + 1));
        Exp_pts_generated = (iter*Exp_pts_generated)/((double) (iter + 1)) + count_pts_generated/((double) (iter + 1));
        //printf("\n------------- Iteration  %d -------------\n", iter + 1);
		//printf("Expected number of points of the HS model : %lf\n", npoints);
        //printf("Expected number of points generated : %lf\n", Exp_pts_generated);

		// Freeing the memory
		freeMem(head);
	}
    clock_t end_time = clock();
    printf("  Intensity = %lf \n  eta = %lf \n  r = %lf\n  Radius = %lf \n  Sigma_len: %d\n  Grid_size: %d\n", lambda, eta, r, rad, sigma_len, grid_size);
    printf("  Expected number of points of the HS model : %lf\n", npoints);
    printf("  Expected number of points generated : %lf\n", Exp_pts_generated);
    printf("  Running time for %d sample : %lf\n", Itot, (double) (end_time - start_time)/((double) CLOCKS_PER_SEC));
	currentTime();

	gsl_rng_free (rng);
} // End of the main function
