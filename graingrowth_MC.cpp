// graingrowth.hpp
// Algorithms for 2D and 3D isotropic Monte Carlo grain growth
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef GRAINGROWTH_UPDATE
#define GRAINGROWTH_UPDATE
#include <iomanip>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <pthread.h>
#include "rdtsc.h"
#include"graingrowth_MC.hpp"
#include"MMSP.hpp"
#include"tessellate.hpp"
#include"output.cpp"

double lambda = 1.0e-3/256; //This is fixed from Monte Carlo simulation, so do not change it.  here 10^-3mm is the domain size, 10^-3mm = 1 um, so each pixel is 1 nm, all length unit should be with mm.
double L_initial = 30*1.0e-3/256; // initially 30 nm diameter
double L0 = 3*1.0e-3/256; 
double K1 = 0.94514608;
double n1 = 0.48921977;
double Q = 1.2552e5; //fitted from Mcbee, William C., and John A. McComb. "Grain growth in thin aluminum-4% copper alloy films." Thin Solid Films 30.1 (1975): 137-143.
double n = 1/0.18; 
double K_ = 5.3137e-17; // length in mm, time in second  fitted from Mcbee, William C., and John A. McComb. "Grain growth in thin aluminum-4% copper alloy films." Thin Solid Films 30.1 (1975): 137-143.
double R = 8.314;

void print_progress(const int step, const int steps, const int iterations);

namespace MMSP
{
template <int dim>
unsigned long generate(MMSP::grid<dim,int >*& grid, int seeds, int nthreads)
{
	#if (defined CCNI) && (!defined MPI_VERSION)
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	#ifdef MPI_VERSION
	int np = MPI::COMM_WORLD.Get_size();
	#endif

	unsigned long timer=0;
	if (dim == 2) {
		const int edge = 256;
		int number_of_fields(seeds);
		if (number_of_fields==0) number_of_fields = static_cast<int>(float(edge*edge)/(M_PI*1.5*1.5)); /* average grain is a disk of radius XXX
, XXX cannot be smaller than 0.1, or BGQ will abort.*/
		#ifdef MPI_VERSION
		while (number_of_fields % np) --number_of_fields; 
		#endif
		grid = new MMSP::grid<dim,int>(0, 0, edge, 0, edge);

		#ifdef MPI_VERSION
		number_of_fields /= np;
		#endif

		#if (!defined MPI_VERSION) && ((defined CCNI) || (defined BGQ))
		std::cerr<<"Error: CCNI requires MPI."<<std::endl;
		std::exit(1);
		#endif
		timer = tessellate<dim,int>(*grid, number_of_fields, nthreads);
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
	} else if (dim == 3) {
		const int edge = 256;
		int number_of_fields(seeds);
		if (number_of_fields==0) number_of_fields = static_cast<int>(float(edge*edge*edge)/(4./3*M_PI*10.*10.*10.)); // Average grain is a sphere of radius 10 voxels
		#ifdef MPI_VERSION
		while (number_of_fields % np) --number_of_fields;
		#endif
		grid = new MMSP::grid<dim,int>(0, 0, edge, 0, edge, 0, edge);

		#ifdef MPI_VERSION
		number_of_fields /= np;
		#endif

		timer = tessellate<dim,int >(*grid, number_of_fields, nthreads);
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
	}
/*------------------Initial tmc----------------------*/
  double tmc_initial = pow(L_initial/K1/lambda,1.0/n1);
  vector<int> coords (dim,0);
  if(dim==2){
    for(int codx=x0(*grid, 0); codx <= x1(*grid, 0); codx++) 
      for(int cody=x0(*grid, 1); cody <= x1(*grid, 1); cody++){
        coords[0] = codx;
        coords[1] = cody;
        (*grid).AccessToTmc(coords) = tmc_initial;
          (*grid).AccessToTmp(coords) = 1.0e6; //set the initial temp to be large enough such that initial physical time = 0.
      }
  }
  else if(dim==3){
    for(int codx=x0(*grid, 0); codx <= x1(*grid, 0); codx++) 
      for(int cody=x0(*grid, 1); cody <= x1(*grid, 1); cody++) 
        for(int codz=x0(*grid, 2); codz <= x1(*grid, 2); codz++){
          coords[0] = codx;
          coords[1] = cody;
          coords[2] = codz;
          (*grid).AccessToTmc(coords) = tmc_initial;
          (*grid).AccessToTmp(coords) = 1.0e6;
        }
  }
  
/*---------------------------------------------------*/
	return timer;
}

unsigned long generate(int dim, char* filename, int seeds, int nthreads)
{
	#if (defined CCNI) && (!defined MPI_VERSION)
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	int rank=0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

	unsigned long timer = 0;
	if (dim == 2) {
		MMSP::grid<2,int>* grid2=NULL;
		timer = generate<2>(grid2,seeds,nthreads);
		assert(grid2!=NULL);
		#ifdef BGQ
		output_bgq(*grid2, filename);
		#else
		output(*grid2, filename);
		#endif
		#ifndef SILENT
		if (rank==0) std::cout<<"Wrote initial file to "<<filename<<"."<<std::endl;
		#endif
	}

	if (dim == 3) {
		MMSP::grid<3,int>* grid3=NULL;
		timer = generate<3>(grid3,seeds,nthreads);
		assert(grid3!=NULL);
		#ifdef BGQ
		output_bgq(*grid3, filename);
		#else
		output(*grid3, filename);
		#endif
		#ifndef SILENT
		if (rank==0) std::cout<<"Wrote initial file to "<<filename<<"."<<std::endl;
		#endif
	}
	return timer;
}

template <int dim>
unsigned long generate(MMSP::grid<dim,int >*& grid, const char* filename)
{
	#if (defined CCNI) && (!defined MPI_VERSION)
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	#ifdef MPI_VERSION
	int np = MPI::COMM_WORLD.Get_size();
	#endif

	unsigned long timer=0;
	if (dim == 2) {
		const int edge = 256;
		grid = new MMSP::grid<dim,int>(0, 0, edge, 0, edge);
		(*grid).input(filename, 1, false);
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
	} else if (dim == 3) {
		const int edge = 256;
		grid = new MMSP::grid<dim,int>(0, 0, edge, 0, edge, 0, edge);
		(*grid).input(filename, 1, false);
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
	}
	return timer;
}

template <int dim>
unsigned long growthexperiment(MMSP::grid<dim,int >*& grid, const char* filename)
{
	#if (defined CCNI) && (!defined MPI_VERSION)
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	#ifdef MPI_VERSION
	int np = MPI::COMM_WORLD.Get_size();
	#endif
  int nx,ny,nz;
	unsigned long timer=0;
	if (dim == 2) {
	  nx = 256;
		ny = 256;
		grid = new MMSP::grid<dim,int>(0, 0, nx, 0, ny);
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
	} else if (dim == 3) {
		std::cerr << " Initial experiment grain structure input only available for 2D simulation" << std::endl;
	  exit(-1);
	}

/*---------------read grain ID------------------*/
	std::ifstream idtxt(filename, std::ios::in);
	if (!idtxt) {
		std::cerr << "File idtxt error: could not open ";
		std::cerr << filename << "." << std::endl;
	  exit(-1);
  }

  int* ids = new int[nx*ny];
  int grainID;
  int count = 0;
  int fields = 1;
  while(idtxt>>grainID){
    ids[count] = grainID;
    ++count; 
    fields = (grainID>fields)?grainID:fields;
  }
	MPI::COMM_WORLD.Barrier();
  idtxt.close();
/*----------------------------------------------*/

/*------------------Initial tmc----------------------*/
  double tmc_initial = pow(L_initial/K1/lambda,1.0/n1);
  vector<int> coords (dim,0);
  if(dim==2){
    for(int codx=x0(*grid, 0); codx <= x1(*grid, 0); codx++) 
      for(int cody=x0(*grid, 1); cody <= x1(*grid, 1); cody++){
        coords[0] = codx;
        coords[1] = cody;
        (*grid).AccessToTmc(coords) = tmc_initial;
        (*grid).AccessToTmp(coords) = 1.0e6; //set the initial temp to be large enough such that initial physical time = 0.
        (*grid)(coords) = ids[ny*codx+cody];
      }
  }
  else if(dim==3){
    std::cerr<<"readin experimental initial grains only support 2D"<<std::endl;
    exit(1);
  }
/*--------------------------------------------------*/
MPI::COMM_WORLD.Barrier();
  delete [] ids;
  ids = NULL;

	return timer;
}

int LargeNearestInteger(int a, int b){
  if(a%b==0) return a/b;
  else return a/b+1;
}

template <int dim> struct flip_index {
	MMSP::grid<dim, int>* grid;
  int num_of_cells_in_thread;
	int sublattice;
  int num_of_points_to_flip;
  int cell_coord[dim];
  int lattice_cells_each_dimension[dim];
  double Pdenominator;
};

template <int dim> void* flip_index_helper_uniformly( void* s )
{
  srand(time(NULL)); /* seed random number generator */
	flip_index<dim>* ss = static_cast<flip_index<dim>*>(s);
  double kT=0.0;
	vector<int> x (dim,0);
  int first_cell_start_coordinates[dim];
  for(int kk=0; kk<dim; kk++) first_cell_start_coordinates[kk] = x0(*(ss->grid), kk);
  for(int i=0; i<dim; i++){
    if(x0(*(ss->grid), i)%2!=0) first_cell_start_coordinates[i]--;
  }
  int cell_coords_selected[dim];
	for (int hh=0; hh<ss->num_of_points_to_flip; hh++) {
	  // choose a random cell to flip
    int cell_numbering_in_thread = rand()%(ss->num_of_cells_in_thread); //choose a cell to flip, from 0 to num_of_cells_in_thread-1
    if(dim==2){
      cell_coords_selected[dim-1]=((ss->cell_coord)[dim-1]+cell_numbering_in_thread)%(ss->lattice_cells_each_dimension)[dim-1];//1-indexed
      cell_coords_selected[0]=(ss->cell_coord)[0]+(((ss->cell_coord)[dim-1]+cell_numbering_in_thread)/(ss->lattice_cells_each_dimension)[dim-1]);
    }else if(dim==3){
      cell_coords_selected[dim-1]=((ss->cell_coord)[dim-1]+cell_numbering_in_thread)%(ss->lattice_cells_each_dimension)[dim-1];//1-indexed
      cell_coords_selected[1]=(  (ss->cell_coord)[1]+ ((ss->cell_coord)[dim-1]+cell_numbering_in_thread)/(ss->lattice_cells_each_dimension)[dim-1]  )%(ss->lattice_cells_each_dimension)[1];
      cell_coords_selected[0]=(ss->cell_coord)[0]+ ( (ss->cell_coord)[1] + ((ss->cell_coord)[dim-1]+cell_numbering_in_thread)/(ss->lattice_cells_each_dimension)[dim-1] ) /(ss->lattice_cells_each_dimension)[1];
    }
    for(int i=0; i<dim; i++){
      x[i]=first_cell_start_coordinates[i]+2*cell_coords_selected[i];
    }
    if(dim==2){
      switch(ss->sublattice){
        case 0:break;// 0,0
        case 1:x[1]++; break; //0,1
        case 2:x[0]++; break; //1,0
        case 3:x[0]++; x[1]++; break; //1,1
      }
    }else if(dim==3){
      switch(ss->sublattice){
        case 0:break;// 0,0,0
        case 1:x[2]++; break; //0,0,1
        case 2:x[1]++; break; //0,1,0
        case 3:x[2]++; x[1]++; break; //0,1,1
        case 4:x[0]++; break; //1,0,0
        case 5:x[2]++; x[0]++; break; //1,0,1
        case 6:x[1]++; x[0]++; break; //1,1,0
        case 7:x[2]++; x[1]++; x[0]++; //1,1,1
      }
    }
    
    bool site_out_of_domain = false;
    for(int i=0; i<dim; i++){
      if(x[i]<x0(*(ss->grid), i) || x[i]>x1(*(ss->grid), i)){
//      if(x[i]<x0(*(ss->grid), i) || x[i]>x1(*(ss->grid), i)-1){
        site_out_of_domain = true;
        break;//break from the for int i loop
      }
    }
    if(site_out_of_domain == true){
      hh--;
      continue; //continue the int hh loop
    }

    int rank = MPI::COMM_WORLD.Get_rank();
    double temperature=273.0;

		int spin1 = (*(ss->grid))(x);
		// determine neighboring spins
    vector<int> r(dim,0);
    std::vector<int> neighbors;
    neighbors.clear();
    int number_of_same_neighours = 0;
    if(dim==2){
      for (int i=-1; i<=1; i++) {
        for (int j=-1; j<=1; j++) {
          if(!(i==0 && j==0)){
            r[0] = x[0] + i;
            r[1] = x[1] + j;
            int spin = (*(ss->grid))(r);

            neighbors.push_back(spin);
            if(spin==spin1) 
              number_of_same_neighours++;
          }
        }
      }
    }else if(dim==3){
		  for (int i=-1; i<=1; i++){
			  for (int j=-1; j<=1; j++){
			    for (int k=-1; k<=1; k++){
            if(!(i==0 && j==0 && k==0)){
				      r[0] = x[0] + i;
				      r[1] = x[1] + j;
				      r[2] = x[2] + k;
				      int spin = (*(ss->grid))(r);
              neighbors.push_back(spin);
              if(spin==spin1) 
                number_of_same_neighours++;
            }
			    }
        }
		  }
    }
    
    //check if inside a grain
    if(number_of_same_neighours==neighbors.size()){//inside a grain
      continue;//continue for
    }
    //choose a random neighbor spin
    int spin2 = neighbors[rand()%neighbors.size()];

		if (spin1!=spin2){
			// compute energy change
			double dE = 0.0;
      if(dim==2){
			  for (int i=-1; i<=1; i++){
				  for (int j=-1; j<=1; j++){
            if(!(i==0 && j==0)){
  					  r[0] = x[0] + i;
	  				  r[1] = x[1] + j;
    				  int spin = (*(ss->grid))(r);
				      dE += 1.0/2*((spin!=spin2)-(spin!=spin1));
            }// if(!(i==0 && j==0))
				  }
        }
      }
      if(dim==3){
			  for (int i=-1; i<=1; i++){ 
				  for (int j=-1; j<=1; j++){ 
    	      for (int k=-1; k<=1; k++){ 
              if(!(i==0 && j==0 && k==0)){
					      r[0] = x[0] + i;
					      r[1] = x[1] + j;
					      r[2] = x[2] + k;
					      int spin = (*(ss->grid))(r);
					      dE += (spin!=spin2)-(spin!=spin1);
              }
				    }
          }
        }
      }
			// attempt a spin flip
			double r = double(rand())/double(RAND_MAX);
      kT = 1.3806488e-23*temperature;
			if (dE <= 0.0) {
        (*(ss->grid))(x) = spin2;
      }
  	  else if (r<exp(-dE/kT)) (*(ss->grid))(x) = spin2;
		}
	}
	pthread_exit(0);
	return NULL;
}

template <int dim> void* flip_index_helper( void* s )
{
  double cos_tolerance=0.999; //cos^-1(0.999)=2.563 deg   cos^-1(0.9999)=0.8103 deg
  srand(time(NULL)); /* seed random number generator */
	flip_index<dim>* ss = static_cast<flip_index<dim>*>(s);
  double kT=0.0;
	vector<int> x (dim,0);
  int first_cell_start_coordinates[dim];
  for(int kk=0; kk<dim; kk++) first_cell_start_coordinates[kk] = x0(*(ss->grid), kk);
  for(int i=0; i<dim; i++){
    if(x0(*(ss->grid), i)%2!=0) first_cell_start_coordinates[i]--;
  }
  int cell_coords_selected[dim];
	for (int hh=0; hh<ss->num_of_points_to_flip; hh++) {
	  // choose a random cell to flip
    int cell_numbering_in_thread = rand()%(ss->num_of_cells_in_thread); //choose a cell to flip, from 0 to num_of_cells_in_thread-1
    if(dim==2){
      cell_coords_selected[dim-1]=((ss->cell_coord)[dim-1]+cell_numbering_in_thread)%(ss->lattice_cells_each_dimension)[dim-1];//1-indexed
      cell_coords_selected[0]=(ss->cell_coord)[0]+(((ss->cell_coord)[dim-1]+cell_numbering_in_thread)/(ss->lattice_cells_each_dimension)[dim-1]);
    }else if(dim==3){
      cell_coords_selected[dim-1]=((ss->cell_coord)[dim-1]+cell_numbering_in_thread)%(ss->lattice_cells_each_dimension)[dim-1];//1-indexed
      cell_coords_selected[1]=(  (ss->cell_coord)[1]+ ((ss->cell_coord)[dim-1]+cell_numbering_in_thread)/(ss->lattice_cells_each_dimension)[dim-1]  )%(ss->lattice_cells_each_dimension)[1];
      cell_coords_selected[0]=(ss->cell_coord)[0]+ ( (ss->cell_coord)[1] + ((ss->cell_coord)[dim-1]+cell_numbering_in_thread)/(ss->lattice_cells_each_dimension)[dim-1] ) /(ss->lattice_cells_each_dimension)[1];
    }
    for(int i=0; i<dim; i++){
      x[i]=first_cell_start_coordinates[i]+2*cell_coords_selected[i];
    }
    if(dim==2){
      switch(ss->sublattice){
        case 0:break;// 0,0
        case 1:x[1]++; break; //0,1
        case 2:x[0]++; break; //1,0
        case 3:x[0]++; x[1]++; break; //1,1
      }
    }else if(dim==3){
      switch(ss->sublattice){
        case 0:break;// 0,0,0
        case 1:x[2]++; break; //0,0,1
        case 2:x[1]++; break; //0,1,0
        case 3:x[2]++; x[1]++; break; //0,1,1
        case 4:x[0]++; break; //1,0,0
        case 5:x[2]++; x[0]++; break; //1,0,1
        case 6:x[1]++; x[0]++; break; //1,1,0
        case 7:x[2]++; x[1]++; x[0]++; //1,1,1
      }
    }

    bool site_out_of_domain = false;
    for(int i=0; i<dim; i++){
      if(x[i]<x0(*(ss->grid), i) || x[i]>=x1(*(ss->grid), i)){
//      if(x[i]<x0(*(ss->grid), i) || x[i]>x1(*(ss->grid), i)-1){
        site_out_of_domain = true;
        break;//break from the for int i loop
      }
    }
    if(site_out_of_domain == true){
      hh--;
      continue; //continue the int hh loop
    }

    int rank = MPI::COMM_WORLD.Get_rank();
//if(rank==0) std::cout<<"num_of_points_to_flip is "<<ss->num_of_points_to_flip<<std::endl;
//getchar();
    double temperature=(*(ss->grid)).AccessToTmp(x);
//    double initial_physical_time=1.0/K_/exp(-Q/R/temperature)*(pow(L_initial,n)-pow(L0,n));
//    double t_mcs_initial = pow(L_initial/K1/lambda,1.0/n1);
    double t_mcs = (*(ss->grid)).AccessToTmc(x);
    double Pnumerator = exp(-Q/R/temperature)/pow(t_mcs,(n*n1-1));
    double site_selection_probability = Pnumerator/ss->Pdenominator;
//if(rank==0) std::cout<<"site_selection_probability is "<<site_selection_probability<<std::endl;
	  double rd = double(rand())/double(RAND_MAX);
    if(rd>site_selection_probability){
//      hh--;// no need to guarantee that N times selection is performed in a configurational MC step, so hh is ony need to be the same as in uniform temp case for the highest temp site
      continue;//this site wont be selected
    }

		int spin1 = (*(ss->grid))(x);
		// determine neighboring spins
    vector<int> r(dim,0);
    std::vector<int> neighbors;
    neighbors.clear();
    int number_of_same_neighours = 0;
    if(dim==2){
      for (int i=-1; i<=1; i++) {
        for (int j=-1; j<=1; j++) {
          if(!(i==0 && j==0)){
            r[0] = x[0] + i;
            r[1] = x[1] + j;
            if(r[0]<g0(*(ss->grid), 0) || r[0]>=g1(*(ss->grid), 0) || r[1]<g0(*(ss->grid), 1) || r[1]>=g1(*(ss->grid), 1) )
              continue;// neighbour outside the global boundary, skip it. 

            int spin = (*(ss->grid))(r);
            neighbors.push_back(spin);
            if(spin==spin1) 
              number_of_same_neighours++;
          }
        }
      }
    }else if(dim==3){
		  for (int i=-1; i<=1; i++){
			  for (int j=-1; j<=1; j++){
			    for (int k=-1; k<=1; k++) {
            if(!(i==0 && j==0 && k==0)){
				      r[0] = x[0] + i;
				      r[1] = x[1] + j;
				      r[2] = x[2] + k;

              if(r[0]<g0(*(ss->grid), 0) || r[0]>g1(*(ss->grid), 0) || r[1]<g0(*(ss->grid), 1) || r[1]>g1(*(ss->grid), 1) ||
                 r[2]<g0(*(ss->grid), 2) || r[2]>g1(*(ss->grid), 2))
                continue;// neighbour outside the global boundary, skip it. 
				      int spin = (*(ss->grid))(r);
              neighbors.push_back(spin);
              if(spin==spin1) 
                number_of_same_neighours++;
            }
			    }
        }
		  }
    }
    
    //check if inside a grain
    if(number_of_same_neighours==neighbors.size()){//inside a grain
      continue;//continue for
    }
    //     choose a random neighbor spin
    int spin2 = neighbors[rand()%neighbors.size()];
		// choose a random spin from Q states
 //       int spin2 = rand()%200;
		if (spin1!=spin2){
			// compute energy change
			double dE = 0.0;
      if(dim==2){
/*
        double film_thickness = 1.0e-6;
        double h2=sin(((ss->grain_orientations)[spin2]).psi)*sin(((ss->grain_orientations)[spin2]).phi_two); 
        double k2=sin(((ss->grain_orientations)[spin2]).psi)*cos(((ss->grain_orientations)[spin2]).phi_two);  
        double l2=cos(((ss->grain_orientations)[spin2]).psi);
        double h1=sin(((ss->grain_orientations)[spin1]).psi)*sin(((ss->grain_orientations)[spin1]).phi_two);   
        double k1=sin(((ss->grain_orientations)[spin1]).psi)*cos(((ss->grain_orientations)[spin1]).phi_two);  
        double l1=cos(((ss->grain_orientations)[spin1]).psi);
*/
	//      		  	  dE += 2*( SurfaceEnergy(h2, k2, l2, temperature)-SurfaceEnergy(h1, k1, l1, temperature)); //surface and interface energy
 //           + film_thickness*(StrainEnergyDenstiy(h2, k2, l2, temperature)-StrainEnergyDenstiy(h1, k1, l1, temperature));//elastic strain energy
			  for (int i=-1; i<=1; i++){
				  for (int j=-1; j<=1; j++){
            if(!(i==0 && j==0)){
  					  r[0] = x[0] + i;
	  				  r[1] = x[1] + j;

	            if(r[0]<g0(*(ss->grid), 0) || r[0]>=g1(*(ss->grid), 0) || r[1]<g0(*(ss->grid), 1) || r[1]>=g1(*(ss->grid), 1) )
              continue;// neighbour outside the global boundary, skip it. 
    				  int spin = (*(ss->grid))(r);
				  dE += 1.0/2*((spin!=spin2)-(spin!=spin1));
				    //				                dE += 1.0/2*((spin!=spin2)-(spin!=spin1))*film_thickness*LargeAngleGrainBoundaryEnergy(temperature); // grain boundary energy  
/*
              bool incoherent_twin_boundary = false, coherent_twin_boundary = false;
              double max,medium,min;
              //check with index before
              double h=sin(((ss->grain_orientations)[spin]).psi)*sin(((ss->grain_orientations)[spin]).phi_two);  
              double k=sin(((ss->grain_orientations)[spin]).psi)*cos(((ss->grain_orientations)[spin]).phi_two);  
              double l=cos(((ss->grain_orientations)[spin]).psi);
              if(CheckSixtyDegreeRotation(h,k,l,h1,k1,l1)){// if h k l is 60 degree away from h1 k1 l1 w.r.t. rotating around <111>
                incoherent_twin_boundary = true;
                if(i==0){//i=0 j!=0   x[0] along rolling direction
                  double transverse_direction[3];
                  //check self
                  transverse_direction[0]=sin(((ss->grain_orientations)[spin1]).phi_one)*cos(((ss->grain_orientations)[spin1]).phi_two)
                                         +cos(((ss->grain_orientations)[spin1]).phi_one)*sin(((ss->grain_orientations)[spin1]).phi_two)
                                          *cos(((ss->grain_orientations)[spin1]).psi);  
                  transverse_direction[1]=-sin(((ss->grain_orientations)[spin1]).phi_one)*sin(((ss->grain_orientations)[spin1]).phi_two)
                                         +cos(((ss->grain_orientations)[spin1]).phi_one)*cos(((ss->grain_orientations)[spin1]).phi_two)
                                          *cos(((ss->grain_orientations)[spin1]).psi);  
                  transverse_direction[2]=-cos(((ss->grain_orientations)[spin1]).phi_one)*sin(((ss->grain_orientations)[spin1]).psi); 
                  transverse_direction[0]=fabs(transverse_direction[0]);
                  transverse_direction[1]=fabs(transverse_direction[1]);
                  transverse_direction[2]=fabs(transverse_direction[2]);
                  max=MaxOfThreeNumber(transverse_direction[0],transverse_direction[1],transverse_direction[2]);
                  medium=MediumOfThreeNumber(transverse_direction[0],transverse_direction[1],transverse_direction[2]);
                  min=MinOfThreeNumber(transverse_direction[0],transverse_direction[1],transverse_direction[2]);
                  double cosine_td_trione=fabs((-max+medium+min))/sqrt(3.0)/sqrt(max*max+medium*medium+min*min);
                  //if coheret twin
                  //check neighbour
                  if(cosine_td_trione>cos_tolerance){
                    transverse_direction[0]=sin(((ss->grain_orientations)[spin]).phi_one)*cos(((ss->grain_orientations)[spin]).phi_two)
                                           +cos(((ss->grain_orientations)[spin]).phi_one)*sin(((ss->grain_orientations)[spin]).phi_two)
                                            *cos(((ss->grain_orientations)[spin]).psi);  
                    transverse_direction[1]=-sin(((ss->grain_orientations)[spin]).phi_one)*sin(((ss->grain_orientations)[spin]).phi_two)
                                           +cos(((ss->grain_orientations)[spin]).phi_one)*cos(((ss->grain_orientations)[spin]).phi_two)
                                            *cos(((ss->grain_orientations)[spin]).psi);  
                    transverse_direction[2]=-cos(((ss->grain_orientations)[spin]).phi_one)*sin(((ss->grain_orientations)[spin]).psi); 
                    transverse_direction[0]=fabs(transverse_direction[0]);
                    transverse_direction[1]=fabs(transverse_direction[1]);
                    transverse_direction[2]=fabs(transverse_direction[2]);
                    max=MaxOfThreeNumber(transverse_direction[0],transverse_direction[1],transverse_direction[2]);
                    medium=MediumOfThreeNumber(transverse_direction[0],transverse_direction[1],transverse_direction[2]);
                    min=MinOfThreeNumber(transverse_direction[0],transverse_direction[1],transverse_direction[2]);
                    cosine_td_trione=fabs((-max+medium+min))/sqrt(3.0)/sqrt(max*max+medium*medium+min*min);
                    if(cosine_td_trione>cos_tolerance){
                      incoherent_twin_boundary = false;
                      coherent_twin_boundary = true;
                      dE += -0.5*film_thickness*CoherentTwinBoundaryEnergy(temperature);
                    }
                  }
                }else{// i!=0 j=0  x[1] along rolling direction
                  double rolling_direction[3];
                  double cosine_rd_trione;
                  //check self
                  rolling_direction[0]=cos(((ss->grain_orientations)[spin1]).phi_one)*cos(((ss->grain_orientations)[spin1]).phi_two)
                                      -sin(((ss->grain_orientations)[spin1]).phi_one)*sin(((ss->grain_orientations)[spin1]).phi_two)
                                       *cos(((ss->grain_orientations)[spin1]).psi);  
                  rolling_direction[1]=-cos(((ss->grain_orientations)[spin1]).phi_one)*sin(((ss->grain_orientations)[spin1]).phi_two)
                                      -sin(((ss->grain_orientations)[spin1]).phi_one)*cos(((ss->grain_orientations)[spin1]).phi_two)
                                       *cos(((ss->grain_orientations)[spin1]).psi); 
                  rolling_direction[2]=sin(((ss->grain_orientations)[spin1]).phi_one)*sin(((ss->grain_orientations)[spin1]).psi); 
                  rolling_direction[0]=fabs(rolling_direction[0]);
                  rolling_direction[1]=fabs(rolling_direction[1]);
                  rolling_direction[2]=fabs(rolling_direction[2]);
                  double max=MaxOfThreeNumber(rolling_direction[0],rolling_direction[1],rolling_direction[2]);
                  double medium=MediumOfThreeNumber(rolling_direction[0],rolling_direction[1],rolling_direction[2]);
                  double min=MinOfThreeNumber(rolling_direction[0],rolling_direction[1],rolling_direction[2]);
                  cosine_rd_trione=fabs((-max+medium+min))/sqrt(3.0)/sqrt(max*max+medium*medium+min*min);
                  //check neighbour
                  if(cosine_rd_trione>cos_tolerance){
                    rolling_direction[0]=cos(((ss->grain_orientations)[spin1]).phi_one)*cos(((ss->grain_orientations)[spin1]).phi_two)
                                        -sin(((ss->grain_orientations)[spin1]).phi_one)*sin(((ss->grain_orientations)[spin1]).phi_two)
                                         *cos(((ss->grain_orientations)[spin1]).psi);  
                    rolling_direction[1]=-cos(((ss->grain_orientations)[spin1]).phi_one)*sin(((ss->grain_orientations)[spin1]).phi_two)
                                        -sin(((ss->grain_orientations)[spin1]).phi_one)*cos(((ss->grain_orientations)[spin1]).phi_two)
                                         *cos(((ss->grain_orientations)[spin1]).psi); 
                    rolling_direction[2]=sin(((ss->grain_orientations)[spin1]).phi_one)*sin(((ss->grain_orientations)[spin1]).psi); 
                    rolling_direction[0]=fabs(rolling_direction[0]);
                    rolling_direction[1]=fabs(rolling_direction[1]);
                    rolling_direction[2]=fabs(rolling_direction[2]);
                    max=MaxOfThreeNumber(rolling_direction[0],rolling_direction[1],rolling_direction[2]);
                    medium=MediumOfThreeNumber(rolling_direction[0],rolling_direction[1],rolling_direction[2]);
                    min=MinOfThreeNumber(rolling_direction[0],rolling_direction[1],rolling_direction[2]);
                    cosine_rd_trione=fabs((-max+medium+min))/sqrt(3.0)/sqrt(max*max+medium*medium+min*min);
                    cosine_rd_trione=fabs((-max+medium+min))/sqrt(3.0)/sqrt(max*max+medium*medium+min*min);
                    if(cosine_rd_trione>cos_tolerance){
                      incoherent_twin_boundary = false;
                      coherent_twin_boundary = true;
                      dE += -0.5*film_thickness*CoherentTwinBoundaryEnergy(temperature);
                    }
                  }
                }
              }
              if(incoherent_twin_boundary == true){
                dE += -0.5*film_thickness*IncoherentTwinBoundaryEnergy(temperature);
              }
              else if(incoherent_twin_boundary == false && coherent_twin_boundary == false){
                dE += -0.5*(spin!=spin1)*film_thickness*LargeAngleGrainBoundaryEnergy(temperature); // grain boundary energy  
              } 
              //check with index after flip
              incoherent_twin_boundary = false;
              coherent_twin_boundary = false;
              if(CheckSixtyDegreeRotation(h,k,l,h2,k2,l2)){//60 deg rotation from ND to ND
                incoherent_twin_boundary = true;
                if(i==0){//i=0 j!=0   x[0] along rolling direction
                  double transverse_direction[3];
                  //check self
                  transverse_direction[0]=sin(((ss->grain_orientations)[spin2]).phi_one)*cos(((ss->grain_orientations)[spin2]).phi_two)
                                         +cos(((ss->grain_orientations)[spin2]).phi_one)*sin(((ss->grain_orientations)[spin2]).phi_two)
                                          *cos(((ss->grain_orientations)[spin2]).psi);  
                  transverse_direction[1]=-sin(((ss->grain_orientations)[spin2]).phi_one)*sin(((ss->grain_orientations)[spin2]).phi_two)
                                         +cos(((ss->grain_orientations)[spin2]).phi_one)*cos(((ss->grain_orientations)[spin2]).phi_two)
                                          *cos(((ss->grain_orientations)[spin2]).psi);  
                  transverse_direction[2]=-cos(((ss->grain_orientations)[spin2]).phi_one)*sin(((ss->grain_orientations)[spin2]).psi); 
                  transverse_direction[0]=fabs(transverse_direction[0]);
                  transverse_direction[1]=fabs(transverse_direction[1]);
                  transverse_direction[2]=fabs(transverse_direction[2]);
                  max=MaxOfThreeNumber(transverse_direction[0],transverse_direction[1],transverse_direction[2]);
                  medium=MediumOfThreeNumber(transverse_direction[0],transverse_direction[1],transverse_direction[2]);
                  min=MinOfThreeNumber(transverse_direction[0],transverse_direction[1],transverse_direction[2]);
                  double cosine_td_trione=fabs((-max+medium+min))/sqrt(3.0)/sqrt(max*max+medium*medium+min*min);
                  //check neighbour
                  if(cosine_td_trione>cos_tolerance){
                    transverse_direction[0]=sin(((ss->grain_orientations)[spin]).phi_one)*cos(((ss->grain_orientations)[spin]).phi_two)
                                           +cos(((ss->grain_orientations)[spin]).phi_one)*sin(((ss->grain_orientations)[spin]).phi_two)
                                            *cos(((ss->grain_orientations)[spin]).psi);  
                    transverse_direction[1]=-sin(((ss->grain_orientations)[spin]).phi_one)*sin(((ss->grain_orientations)[spin]).phi_two)
                                           +cos(((ss->grain_orientations)[spin]).phi_one)*cos(((ss->grain_orientations)[spin]).phi_two)
                                            *cos(((ss->grain_orientations)[spin]).psi);  
                    transverse_direction[2]=-cos(((ss->grain_orientations)[spin]).phi_one)*sin(((ss->grain_orientations)[spin]).psi); 
                    transverse_direction[0]=fabs(transverse_direction[0]);
                    transverse_direction[1]=fabs(transverse_direction[1]);
                    transverse_direction[2]=fabs(transverse_direction[2]);
                    max=MaxOfThreeNumber(transverse_direction[0],transverse_direction[1],transverse_direction[2]);
                    medium=MediumOfThreeNumber(transverse_direction[0],transverse_direction[1],transverse_direction[2]);
                    min=MinOfThreeNumber(transverse_direction[0],transverse_direction[1],transverse_direction[2]);
                    cosine_td_trione=fabs((-max+medium+min))/sqrt(3.0)/sqrt(max*max+medium*medium+min*min);
                    if(cosine_td_trione>cos_tolerance){
                      incoherent_twin_boundary = false;
                      coherent_twin_boundary = true;
                      dE += 0.5*film_thickness*CoherentTwinBoundaryEnergy(temperature);
                    }
                  }
                }else{// i!=0 j=0  x[1] along rolling direction
                  double rolling_direction[3];
                  double cosine_rd_trione;
                  //check self
                  rolling_direction[0]=cos(((ss->grain_orientations)[spin2]).phi_one)*cos(((ss->grain_orientations)[spin2]).phi_two)
                                      -sin(((ss->grain_orientations)[spin2]).phi_one)*sin(((ss->grain_orientations)[spin2]).phi_two)
                                       *cos(((ss->grain_orientations)[spin2]).psi);  
                  rolling_direction[1]=-cos(((ss->grain_orientations)[spin2]).phi_one)*sin(((ss->grain_orientations)[spin2]).phi_two)
                                      -sin(((ss->grain_orientations)[spin2]).phi_one)*cos(((ss->grain_orientations)[spin2]).phi_two)
                                       *cos(((ss->grain_orientations)[spin2]).psi); 
                  rolling_direction[2]=sin(((ss->grain_orientations)[spin2]).phi_one)*sin(((ss->grain_orientations)[spin2]).psi); 
                  rolling_direction[0]=fabs(rolling_direction[0]);
                  rolling_direction[1]=fabs(rolling_direction[1]);
                  rolling_direction[2]=fabs(rolling_direction[2]);
                  double max=MaxOfThreeNumber(rolling_direction[0],rolling_direction[1],rolling_direction[2]);
                  double medium=MediumOfThreeNumber(rolling_direction[0],rolling_direction[1],rolling_direction[2]);
                  double min=MinOfThreeNumber(rolling_direction[0],rolling_direction[1],rolling_direction[2]);
                  cosine_rd_trione=fabs((-max+medium+min))/sqrt(3.0)/sqrt(max*max+medium*medium+min*min);
                  //check neighbour
                  if(cosine_rd_trione>cos_tolerance){
                    rolling_direction[0]=cos(((ss->grain_orientations)[spin2]).phi_one)*cos(((ss->grain_orientations)[spin2]).phi_two)
                                        -sin(((ss->grain_orientations)[spin2]).phi_one)*sin(((ss->grain_orientations)[spin2]).phi_two)
                                         *cos(((ss->grain_orientations)[spin2]).psi);  
                    rolling_direction[1]=-cos(((ss->grain_orientations)[spin2]).phi_one)*sin(((ss->grain_orientations)[spin2]).phi_two)
                                        -sin(((ss->grain_orientations)[spin2]).phi_one)*cos(((ss->grain_orientations)[spin2]).phi_two)
                                         *cos(((ss->grain_orientations)[spin2]).psi); 
                    rolling_direction[2]=sin(((ss->grain_orientations)[spin2]).phi_one)*sin(((ss->grain_orientations)[spin2]).psi); 
                    rolling_direction[0]=fabs(rolling_direction[0]);
                    rolling_direction[1]=fabs(rolling_direction[1]);
                    rolling_direction[2]=fabs(rolling_direction[2]);
                    max=MaxOfThreeNumber(rolling_direction[0],rolling_direction[1],rolling_direction[2]);
                    medium=MediumOfThreeNumber(rolling_direction[0],rolling_direction[1],rolling_direction[2]);
                    min=MinOfThreeNumber(rolling_direction[0],rolling_direction[1],rolling_direction[2]);
                    cosine_rd_trione=fabs((-max+medium+min))/sqrt(3.0)/sqrt(max*max+medium*medium+min*min);
                    cosine_rd_trione=fabs((-max+medium+min))/sqrt(3.0)/sqrt(max*max+medium*medium+min*min);
                    if(cosine_rd_trione>cos_tolerance){
                      incoherent_twin_boundary = false;
                      coherent_twin_boundary = true;
                      dE += 0.5*film_thickness*CoherentTwinBoundaryEnergy(temperature);
                    }
                  }
                }
              }
              if(incoherent_twin_boundary == true){
                dE += 0.5*film_thickness*IncoherentTwinBoundaryEnergy(temperature);
              }
              else if(incoherent_twin_boundary == false && coherent_twin_boundary == false){
                dE += 0.5*(spin!=spin2)*film_thickness*LargeAngleGrainBoundaryEnergy(temperature); // grain boundary energy  
		          }*/
            }// if(!(i==0 && j==0))
				  }
        }
      }
      if(dim==3){
			  for (int i=-1; i<=1; i++){ 
				  for (int j=-1; j<=1; j++){ 
    	      for (int k=-1; k<=1; k++){ 
              if(!(i==0 && j==0 && k==0)){
					      r[0] = x[0] + i;
					      r[1] = x[1] + j;
					      r[2] = x[2] + k;

                if(r[0]<g0(*(ss->grid), 0) || r[0]>g1(*(ss->grid), 0) || r[1]<g0(*(ss->grid), 1) || r[1]>g1(*(ss->grid), 1) ||
                   r[2]<g0(*(ss->grid), 2) || r[2]>g1(*(ss->grid), 2))
                  continue;// neighbour outside the global boundary, skip it.
					      int spin = (*(ss->grid))(r);
					      dE += (spin!=spin2)-(spin!=spin1);
              }
				    }
          }
        }
      }
			// attempt a spin flip
			double r = double(rand())/double(RAND_MAX);
      kT = 1.3806488e-23*temperature;
//	      int rank = MPI::COMM_WORLD.Get_rank();
			if (dE <= 0.0) {
        (*(ss->grid))(x) = spin2;
//        if(rank==0 ) {std::cout<<x[0]<<" "<<x[1]<<"  dE is "<<dE<<std::endl;}
      }
  	  else if (r<exp(-dE/kT)) (*(ss->grid))(x) = spin2;
		}
	}
	pthread_exit(0);
	return NULL;
}

template <int dim> bool OutsideDomainCheck(MMSP::grid<dim, int>& grid, vector<int>* x){
  bool outside_domain=false;
  for(int i=0; i<dim; i++){
    if((*x)[i]<x0(grid, i) || (*x)[i]>x1(grid, i)){
//    if((*x)[i]<x0(grid, i) || (*x)[i]>x1(grid, i)-1){
      outside_domain=true;
      break;
    }
  }
  return outside_domain;
}

template <int dim> unsigned long update_uniformly(MMSP::grid<dim, int>& grid, int steps, int nthreads)
{
	#if (!defined MPI_VERSION) && ((defined CCNI) || (defined BGQ))
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	int rank=0;
	unsigned int np=0;
//	#ifdef MPI_VERSION
	rank=MPI::COMM_WORLD.Get_rank();
	np=MPI::COMM_WORLD.Get_size();
	MPI::COMM_WORLD.Barrier();
//	#endif

	unsigned long update_timer = 0;

	pthread_t* p_threads = new pthread_t[nthreads];
	flip_index<dim>* mat_para = new flip_index<dim> [nthreads];
	pthread_attr_t attr;
	pthread_attr_init (&attr);

	#ifndef SILENT
	static int iterations = 1;
	if (rank==0) print_progress(0, steps, iterations);
	#endif

/*
  int edge = 256;
  int number_of_fields = static_cast<int>(float(edge*edge)/(M_PI*10.*10.)); // average grain is a disk of radius 10
  #ifdef MPI_VERSION
	while (number_of_fields % np) --number_of_fields;
  #endif
*/
/*
  EulerAngles *grain_orientations = new EulerAngles[200];
  if(rank==0){
    ReadData(grain_orientations);
  }
	MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Bcast(grain_orientations, 200*3, MPI_DOUBLE, 0);
*/
/*-----------------------------------------------*/
/*---------------generate cells------------------*/
/*-----------------------------------------------*/
  int dimension_length=0, number_of_lattice_cells=1;
  int lattice_cells_each_dimension[dim];
  for(int i=0; i<dim; i++){
    dimension_length = x1(grid, i)-x0(grid, i);
//    dimension_length = x1(grid, i)-1-x0(grid, i);
    if(x0(grid, 0)%2==0)
      lattice_cells_each_dimension[i] = dimension_length/2+1;
    else
      lattice_cells_each_dimension[i] = 1+LargeNearestInteger(dimension_length,2);
    number_of_lattice_cells *= lattice_cells_each_dimension[i];
  }
  #ifndef SILENT
  if(rank==0){ 
    std::cout<<"lattice_cells_each_dimension are:  ";
    for(int i=0; i<dim; i++)
      std::cout<<lattice_cells_each_dimension[i]<<"  ";
    std::cout<<"\n";
  }
  #endif
//----------assign cells for each pthreads
  int num_of_cells_in_thread = number_of_lattice_cells/nthreads;
  //check if num of the pthread is too large, if so, reduce it.                                                                                                                                         
  if (num_of_cells_in_thread<1) {
    std::cerr<<"ERROR: number of pthread is too large, please reduce it to a value <= "<<number_of_lattice_cells<<std::endl;
    exit(0);
  }

	vector<int> x (dim,0);
	vector<int> x_prim (dim,0);
  int coordinates_of_cell[dim];
  int initial_coordinates[dim];
  int **cell_coord = new int*[nthreads];//record the start coordinates of each pthread domain.
  for(int i=0; i<nthreads; i++){
    cell_coord[i] = new int[dim];
    for(int j=0; j<dim; j++){
      cell_coord[i][j]=0;
    }
  }

  int **num_of_grids_to_flip = new int*[nthreads];
  for(int i=0; i<nthreads; i++){
    num_of_grids_to_flip[i] = new int[( static_cast<int>(pow(2,dim)) )];
    for(int j=0; j<pow(2,dim); j++){
      num_of_grids_to_flip[i][j]=0;
    }
  }

  for(int k=0; k<dim; k++) 
    initial_coordinates[k] = x0(grid, k);
  for(int i=0; i<dim; i++){
    if(x0(grid, i)%2!=0) 
      initial_coordinates[i]--;
  }

  for (int i=0; i<nthreads; i++) {
        int cell_numbering = num_of_cells_in_thread*i; //0-indexed, celling_numbering is the start cell numbering
        if(dim==2){
          cell_coord[i][dim-1]=cell_numbering%lattice_cells_each_dimension[dim-1];//0-indexed
          cell_coord[i][0]=(cell_numbering/lattice_cells_each_dimension[dim-1]);
	        if(cell_coord[i][0]>=lattice_cells_each_dimension[0]){
	          std::cerr<<"the cell coordinates is wrong!"<<std::endl;
	          exit(1);
	        }
        }else if(dim==3){
          cell_coord[i][dim-1]=cell_numbering%lattice_cells_each_dimension[dim-1];//0-indexed
          cell_coord[i][1]=(cell_numbering/lattice_cells_each_dimension[dim-1])%lattice_cells_each_dimension[1];
          cell_coord[i][0]=(cell_numbering/lattice_cells_each_dimension[dim-1])/lattice_cells_each_dimension[1];
	        if(cell_coord[i][0]>=lattice_cells_each_dimension[0]){
	          std::cerr<<"the cell coordinates is wrong!"<<std::endl;
	          exit(1);
	        }
        }

    mat_para[i].grid = &grid;
    if(i==(nthreads-1)) 
      mat_para[i].num_of_cells_in_thread = number_of_lattice_cells - num_of_cells_in_thread*(nthreads-1);
    else 
      mat_para[i].num_of_cells_in_thread = num_of_cells_in_thread;

    #ifndef SILENT
    if(rank==0) std::cout<<"num_of_cells_in_thread is "<<mat_para[i].num_of_cells_in_thread<<" in thread "<<i<<"\n";
    #endif

    for(int k=0; k<dim; k++) 
      mat_para[i].lattice_cells_each_dimension[k]=lattice_cells_each_dimension[k];

//    mat_para[i].temperature_along_x = temperature_along_x;
//    mat_para[i].grain_orientations = grain_orientations;

    for(int j=0; j<mat_para[i].num_of_cells_in_thread; j++){
      int start_cell_numbering = num_of_cells_in_thread*i;
      if(dim==2){
        coordinates_of_cell[dim-1]=(start_cell_numbering+j)%lattice_cells_each_dimension[dim-1];//0-indexed
        coordinates_of_cell[0]=(start_cell_numbering+j)/lattice_cells_each_dimension[dim-1];
      }else if(dim==3){
        coordinates_of_cell[dim-1]=(start_cell_numbering+j)%lattice_cells_each_dimension[dim-1];//0-indexed
        coordinates_of_cell[1]=((start_cell_numbering+j)/lattice_cells_each_dimension[dim-1])%lattice_cells_each_dimension[1];
        coordinates_of_cell[0]=((start_cell_numbering+j)/lattice_cells_each_dimension[dim-1])/lattice_cells_each_dimension[1];
      }
      for(int ii=0; ii<dim; ii++){
        x[ii]=initial_coordinates[ii]+2*coordinates_of_cell[ii];
      }

      if(dim==2){
        x_prim = x;
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][0]+=1;

        x_prim = x;
        x_prim[1]=x[1]+1; //0,1 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][1]+=1;

        x_prim = x;
        x_prim[0]=x[0]+1; //1,0
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][2]+=1;

        x_prim = x;
        x_prim[0]=x[0]+1;
        x_prim[1]=x[1]+1; //1,1 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][3]+=1;
      }else if(dim==3){
        x_prim = x;//0,0,0
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][0]+=1;

        x_prim = x;
        x_prim[2]=x[2]+1; //0,0,1 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][1]+=1;

        x_prim = x;
        x_prim[1]=x[1]+1; //0,1,0
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][2]+=1;

        x_prim = x;
        x_prim[2]=x[2]+1;
        x_prim[1]=x[1]+1; //0,1,1 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][3]+=1;

        x_prim = x;
        x_prim[0]=x[0]+1; //1,0,0 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][4]+=1;

        x_prim = x;
        x_prim[2]=x[2]+1;
        x_prim[0]=x[0]+1; //1,0,1 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][5]+=1;

        x_prim = x;
        x_prim[1]=x[1]+1;
        x_prim[0]=x[0]+1; //1,1,0 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][6]+=1;

        x_prim = x;
        x_prim[2]=x[2]+1;
        x_prim[1]=x[1]+1;
        x_prim[0]=x[0]+1; //1,1,1 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][7]+=1;
      }
    }// for int j 
  }//for int i

	for (int step=0; step<steps; step++){
		unsigned long start = rdtsc();
    int num_of_sublattices=0;
    if(dim==2) num_of_sublattices = 4; 
    else if(dim==3) num_of_sublattices = 8;
		for (int sublattice=0; sublattice < num_of_sublattices; sublattice++) {
			for (int i=0; i!= nthreads ; i++) {
				mat_para[i].sublattice=sublattice;
				mat_para[i].num_of_points_to_flip=num_of_grids_to_flip[i][sublattice];
        for(int k=0; k<dim; k++) mat_para[i].cell_coord[k]=cell_coord[i][k];
				pthread_create(&p_threads[i], &attr, flip_index_helper_uniformly<dim>, (void*) &mat_para[i] );
			}//loop over threads

			for (int ii=0; ii!= nthreads ; ii++)
				pthread_join(p_threads[ii], NULL);

			#ifdef MPI_VERSION
			MPI::COMM_WORLD.Barrier();
			#endif

			ghostswap(grid, sublattice); // once looped over a "color", ghostswap.
//			ghostswap(grid); // once looped over a "color", ghostswap.
      #ifdef MPI_VERSION
			MPI::COMM_WORLD.Barrier();
                #endif

		}//loop over color
		#ifndef SILENT
		if (rank==0) print_progress(step+1, steps, iterations);
		#endif
		update_timer += rdtsc()-start;
	}//loop over step
	#ifndef SILENT
	++iterations;
	#endif

  for(int i=0; i<nthreads; i++){
    delete [] num_of_grids_to_flip[i];
    num_of_grids_to_flip[i]=NULL;
    delete [] cell_coord[i];
    cell_coord[i]=NULL;
  }
  delete num_of_grids_to_flip; 
  num_of_grids_to_flip=NULL; 
  delete cell_coord;
  cell_coord=NULL;
//	delete [] grain_orientations;
//	grain_orientations=NULL;
//	delete [] temperature_along_x;
//	temperature_along_x=NULL;
	delete [] p_threads;
	p_threads=NULL;
	delete [] mat_para;
	mat_para=NULL;
	unsigned long total_update_time=update_timer;
	#ifdef MPI_VERSION
	MPI::COMM_WORLD.Allreduce(&update_timer, &total_update_time, 1, MPI_UNSIGNED_LONG, MPI_SUM);
	#endif
	return total_update_time/np; // average update time
}

template <int dim> double PdenominatorMax(MMSP::grid<dim, int>& grid, double& tmc_at_PdenominatorMax){
   double Pdenominator_max = 0.0;
   vector<int> coords (dim,0);
   if(dim==2){
       for(int codx=x0(grid, 0); codx <= x1(grid, 0); codx++) 
         for(int cody=x0(grid, 1); cody <= x1(grid, 1); cody++){
           coords[0] = codx;
           coords[1] = cody;
           double temperature = grid.AccessToTmp(coords);
           double Pdenominator = exp(-Q/R/temperature)/pow(grid.AccessToTmc(coords), (n*n1-1));
           if(Pdenominator > Pdenominator_max){
             Pdenominator_max = Pdenominator;
             tmc_at_PdenominatorMax = grid.AccessToTmc(coords);
           }
         }
   }
   else if(dim==3){
     for(int codx=x0(grid, 0); codx <= x1(grid, 0); codx++) 
       for(int cody=x0(grid, 1); cody <= x1(grid, 1); cody++) 
         for(int codz=x0(grid, 2); codz <= x1(grid, 2); codz++){
           coords[0] = codx;
           coords[1] = cody;
           coords[2] = codz;
           double temperature = grid.AccessToTmp(coords);
           double Pdenominator = exp(-Q/R/temperature)/pow(grid.AccessToTmc(coords), (n*n1-1));
           if(Pdenominator > Pdenominator_max){
             Pdenominator_max = Pdenominator;
             tmc_at_PdenominatorMax = grid.AccessToTmc(coords);
           }
        }
   }
   return Pdenominator_max;
}

template <int dim> void UpdateLocalTmc(MMSP::grid<dim, int>& grid, double t_inc){
   vector<int> coords (dim,0);
   if(dim==2){
       for(int codx=x0(grid, 0); codx <= x1(grid, 0); codx++) 
         for(int cody=x0(grid, 1); cody <= x1(grid, 1); cody++){
           coords[0] = codx;
           coords[1] = cody;
           long double exp_eqn_rhs = pow(K1*lambda*pow(grid.AccessToTmc(coords), n1), n) - pow(L0, n);
           long double temperature = grid.AccessToTmp(coords);
           exp_eqn_rhs += K_*t_inc*exp(-Q/R/temperature);
           grid.AccessToTmc(coords) = pow(1.0/K1/lambda*pow(exp_eqn_rhs+pow(L0, n),1.0/n), 1.0/n1);  
         }
   }
   else if(dim==3){
     for(int codx=x0(grid, 0); codx <= x1(grid, 0); codx++)  
         for(int cody=x0(grid, 1); cody <= x1(grid, 1); cody++) 
           for(int codz=x0(grid, 2); codz <= x1(grid, 2); codz++){
             coords[0] = codx;
             coords[1] = cody;
             coords[2] = codz;
             long double exp_eqn_rhs = pow(K1*lambda*pow(grid.AccessToTmc(coords), n1), n) - pow(L0, n);
             long double temperature = grid.AccessToTmp(coords);
             exp_eqn_rhs += K_*t_inc*exp(-Q/R/temperature);
             grid.AccessToTmc(coords) = pow(1.0/K1/lambda*pow(exp_eqn_rhs+pow(L0, n),1.0/n), 1.0/n1);  
           }
   }
}

template <int dim> void UpdateLocalTmp(MMSP::grid<dim, int>& grid, long double physical_time, double* temp){
   vector<int> coords (dim,0);
   if(dim==2){
       for(int codx=x0(grid, 0); codx <= x1(grid, 0); codx++) 
         for(int cody=x0(grid, 1); cody <= x1(grid, 1); cody++){
           coords[0] = codx;
           coords[1] = cody;
           grid.AccessToTmp(coords) = temp[1]+(temp[0]-temp[1])/256*codx;
         }
/*-----------------------
if(codx<=0.25*256)
  grid.AccessToTmp(coords) = 100; 
else if(0.25*256<codx<=0.5*256)
  grid.AccessToTmp(coords) = 300; 
else if(0.5*256<codx<=0.75*256)
  grid.AccessToTmp(coords) = 256; 
else if(0.75*256<codx<=256)
  grid.AccessToTmp(coords) = 700; 
-----------------------*/
   }
   else if(dim==3){
     for(int codx=x0(grid, 0); codx <= x1(grid, 0); codx++)  
         for(int cody=x0(grid, 1); cody <= x1(grid, 1); cody++) 
           for(int codz=x0(grid, 2); codz <= x1(grid, 2); codz++){
             coords[0] = codx;
             coords[1] = cody;
             coords[2] = codz;
             grid.AccessToTmp(coords) = 473 + 273.0*sin(3.14/256*codx + 3.14/(1.0e5)*physical_time);
           }
   }
}

template <int dim> unsigned long update(MMSP::grid<dim, int>& grid, int steps, int steps_finished, int nthreads, int step_to_nonuniform, long double &physical_time, double* temp)
{
	#if (!defined MPI_VERSION) && ((defined CCNI) || (defined BGQ))
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	int rank=0;
	unsigned int np=0;
//	#ifdef MPI_VERSION
	rank=MPI::COMM_WORLD.Get_rank();
	np=MPI::COMM_WORLD.Get_size();
	MPI::COMM_WORLD.Barrier();
//	#endif

	unsigned long update_timer = 0;

	pthread_t* p_threads = new pthread_t[nthreads];
	flip_index<dim>* mat_para = new flip_index<dim> [nthreads];
	pthread_attr_t attr;
	pthread_attr_init (&attr);

	#ifndef SILENT
	static int iterations = 1;
	if (rank==0) print_progress(0, steps, iterations);
	#endif
//	ghostswap(grid); 
/*
  int edge = 256;
  int number_of_fields = static_cast<int>(float(edge*edge)/(M_PI*10.*10.)); // average grain is a disk of radius 10
  #ifdef MPI_VERSION
	while (number_of_fields % np) --number_of_fields;
  #endif
*/
/*
  EulerAngles *grain_orientations = new EulerAngles[200];
  if(rank==0){
    ReadData(grain_orientations);
  }
	MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Bcast(grain_orientations, 200*3, MPI_DOUBLE, 0);
*/
/*-----------------------------------------------*/
/*---------------generate cells------------------*/
/*-----------------------------------------------*/
  int dimension_length=0, number_of_lattice_cells=1;
  int lattice_cells_each_dimension[dim];
  for(int i=0; i<dim; i++){
    dimension_length = x1(grid, i)-x0(grid, i);
//    dimension_length = x1(grid, i)-1-x0(grid, i);
    if(x0(grid, 0)%2==0)
      lattice_cells_each_dimension[i] = dimension_length/2+1;
    else
      lattice_cells_each_dimension[i] = 1+LargeNearestInteger(dimension_length,2);
    number_of_lattice_cells *= lattice_cells_each_dimension[i];
  }
  #ifndef SILENT
  if(rank==0){ 
    std::cout<<"lattice_cells_each_dimension are:  ";
    for(int i=0; i<dim; i++)
      std::cout<<lattice_cells_each_dimension[i]<<"  ";
    std::cout<<"\n";
  }
  #endif
//----------assign cells for each pthreads
  int num_of_cells_in_thread = number_of_lattice_cells/nthreads;
  //check if num of the pthread is too large, if so, reduce it.                                                                                                                                         
  if (num_of_cells_in_thread<1) {
    std::cerr<<"ERROR: number of pthread is too large, please reduce it to a value <= "<<number_of_lattice_cells<<std::endl;
    exit(0);
  }
  int model_dimension=(g1(grid, 0)-g0(grid, 0)+1);
//  double *temperature_along_x = new double[model_dimension];

//  if(rank==0){
//    ReadTemperature(temperature_along_x, model_dimension);
//  }
//	MPI::COMM_WORLD.Barrier();
//  MPI::COMM_WORLD.Bcast(temperature_along_x, model_dimension, MPI_DOUBLE, 0);

	vector<int> x (dim,0);
	vector<int> x_prim (dim,0);
  int coordinates_of_cell[dim];
  int initial_coordinates[dim];
  int **cell_coord = new int*[nthreads];//record the start coordinates of each pthread domain.
  for(int i=0; i<nthreads; i++){
    cell_coord[i] = new int[dim];
    for(int j=0; j<dim; j++){
      cell_coord[i][j]=0;
    }
  }

  int **num_of_grids_to_flip = new int*[nthreads];
  for(int i=0; i<nthreads; i++){
    num_of_grids_to_flip[i] = new int[( static_cast<int>(pow(2,dim)) )];
    for(int j=0; j<pow(2,dim); j++){
      num_of_grids_to_flip[i][j]=0;
    }
  }

  for(int k=0; k<dim; k++) 
    initial_coordinates[k] = x0(grid, k);
  for(int i=0; i<dim; i++){
    if(x0(grid, i)%2!=0) 
      initial_coordinates[i]--;
  }

  for (int i=0; i<nthreads; i++) {
        int cell_numbering = num_of_cells_in_thread*i; //0-indexed, celling_numbering is the start cell numbering
        if(dim==2){
          cell_coord[i][dim-1]=cell_numbering%lattice_cells_each_dimension[dim-1];//0-indexed
          cell_coord[i][0]=(cell_numbering/lattice_cells_each_dimension[dim-1]);
	        if(cell_coord[i][0]>=lattice_cells_each_dimension[0]){
	          std::cerr<<"the cell coordinates is wrong!"<<std::endl;
	          exit(1);
	        }
        }else if(dim==3){
          cell_coord[i][dim-1]=cell_numbering%lattice_cells_each_dimension[dim-1];//0-indexed
          cell_coord[i][1]=(cell_numbering/lattice_cells_each_dimension[dim-1])%lattice_cells_each_dimension[1];
          cell_coord[i][0]=(cell_numbering/lattice_cells_each_dimension[dim-1])/lattice_cells_each_dimension[1];
	        if(cell_coord[i][0]>=lattice_cells_each_dimension[0]){
	          std::cerr<<"the cell coordinates is wrong!"<<std::endl;
	          exit(1);
	        }
        }

    mat_para[i].grid = &grid;
    if(i==(nthreads-1)) 
      mat_para[i].num_of_cells_in_thread = number_of_lattice_cells - num_of_cells_in_thread*(nthreads-1);
    else 
      mat_para[i].num_of_cells_in_thread = num_of_cells_in_thread;

    #ifndef SILENT
    if(rank==0) std::cout<<"num_of_cells_in_thread is "<<mat_para[i].num_of_cells_in_thread<<" in thread "<<i<<"\n";
    #endif

    for(int k=0; k<dim; k++) 
      mat_para[i].lattice_cells_each_dimension[k]=lattice_cells_each_dimension[k];

//    mat_para[i].temperature_along_x = temperature_along_x;
//    mat_para[i].grain_orientations = grain_orientations;

    for(int j=0; j<mat_para[i].num_of_cells_in_thread; j++){
      int start_cell_numbering = num_of_cells_in_thread*i;
      if(dim==2){
        coordinates_of_cell[dim-1]=(start_cell_numbering+j)%lattice_cells_each_dimension[dim-1];//0-indexed
        coordinates_of_cell[0]=(start_cell_numbering+j)/lattice_cells_each_dimension[dim-1];
      }else if(dim==3){
        coordinates_of_cell[dim-1]=(start_cell_numbering+j)%lattice_cells_each_dimension[dim-1];//0-indexed
        coordinates_of_cell[1]=((start_cell_numbering+j)/lattice_cells_each_dimension[dim-1])%lattice_cells_each_dimension[1];
        coordinates_of_cell[0]=((start_cell_numbering+j)/lattice_cells_each_dimension[dim-1])/lattice_cells_each_dimension[1];
      }
      for(int ii=0; ii<dim; ii++){
        x[ii]=initial_coordinates[ii]+2*coordinates_of_cell[ii];
      }

      if(dim==2){
        x_prim = x;
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][0]+=1;

        x_prim = x;
        x_prim[1]=x[1]+1; //0,1 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][1]+=1;

        x_prim = x;
        x_prim[0]=x[0]+1; //1,0
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][2]+=1;

        x_prim = x;
        x_prim[0]=x[0]+1;
        x_prim[1]=x[1]+1; //1,1 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][3]+=1;
      }else if(dim==3){
        x_prim = x;//0,0,0
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][0]+=1;

        x_prim = x;
        x_prim[2]=x[2]+1; //0,0,1 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][1]+=1;

        x_prim = x;
        x_prim[1]=x[1]+1; //0,1,0
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][2]+=1;

        x_prim = x;
        x_prim[2]=x[2]+1;
        x_prim[1]=x[1]+1; //0,1,1 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][3]+=1;

        x_prim = x;
        x_prim[0]=x[0]+1; //1,0,0 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][4]+=1;

        x_prim = x;
        x_prim[2]=x[2]+1;
        x_prim[0]=x[0]+1; //1,0,1 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][5]+=1;

        x_prim = x;
        x_prim[1]=x[1]+1;
        x_prim[0]=x[0]+1; //1,1,0 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][6]+=1;

        x_prim = x;
        x_prim[2]=x[2]+1;
        x_prim[1]=x[1]+1;
        x_prim[0]=x[0]+1; //1,1,1 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][7]+=1;
      }
    }// for int j 
  }//for int i

	for (int step=0; step<steps; step++){

    double tmc_at_PdenominatorMax = 0.0;
    double Pdenominator_max_partition = PdenominatorMax(grid,tmc_at_PdenominatorMax); 
    MPI::COMM_WORLD.Barrier();
    double Pdenominator_max_global = 0.0;

    MPI::COMM_WORLD.Allreduce(&Pdenominator_max_partition, &Pdenominator_max_global, 1, MPI_DOUBLE, MPI_MAX);
    MPI::COMM_WORLD.Barrier();

    double tmc_at_PdenominatorMax_global = 0.0;
    if(Pdenominator_max_partition != Pdenominator_max_global){
      tmc_at_PdenominatorMax = 0.0;
    }
    MPI::COMM_WORLD.Allreduce(&tmc_at_PdenominatorMax, &tmc_at_PdenominatorMax_global, 1, MPI_DOUBLE, MPI_MAX);
    MPI::COMM_WORLD.Barrier();
std::cout<<"tmc_at_PdenominatorMax_global is "<<tmc_at_PdenominatorMax_global<<std::endl;
//    double t_inc = n*n1*pow(K1*lambda,n)/K_/Pdenominator_max_global;

    vector<int> coords (dim,0);
    coords[0] = x0(grid, 0);
    coords[1] = x0(grid, 1);
    coords[dim] = x0(grid, dim);
    double uniform_temperature = grid.AccessToTmp(coords);
//std::cout<<"uniform_temperature is "<<uniform_temperature<<std::endl;
    double t_inc = ( pow(K1*lambda*pow(tmc_at_PdenominatorMax_global+1,n1), n) - 
                     pow(K1*lambda*pow(tmc_at_PdenominatorMax_global,n1), n) )/K_/exp(-Q/R/uniform_temperature);
if(rank==0)
std::cout<<"~~ "<<pow(K1*lambda*pow(tmc_at_PdenominatorMax_global+1,n1), n)<<"  ~~  "<<pow(K1*lambda*pow(tmc_at_PdenominatorMax_global,n1), n)<<std::endl;

		unsigned long start = rdtsc();
    int num_of_sublattices=0;
    if(dim==2) num_of_sublattices = 4; 
    else if(dim==3) num_of_sublattices = 8;
		for (int sublattice=0; sublattice < num_of_sublattices; sublattice++) {
			for (int i=0; i!= nthreads ; i++) {
				mat_para[i].sublattice=sublattice;
				mat_para[i].num_of_points_to_flip=num_of_grids_to_flip[i][sublattice];
        mat_para[i].Pdenominator = Pdenominator_max_global;
        for(int k=0; k<dim; k++) mat_para[i].cell_coord[k]=cell_coord[i][k];
				pthread_create(&p_threads[i], &attr, flip_index_helper<dim>, (void*) &mat_para[i] );
			}//loop over threads

			for (int ii=0; ii!= nthreads ; ii++)
				pthread_join(p_threads[ii], NULL);

			#ifdef MPI_VERSION
			MPI::COMM_WORLD.Barrier();
			#endif

			ghostswap(grid, sublattice); // once looped over a "color", ghostswap.
//			ghostswap(grid); // once looped over a "color", ghostswap.
            #ifdef MPI_VERSION
			MPI::COMM_WORLD.Barrier();
                #endif
		}//loop over color

	  MPI::COMM_WORLD.Barrier();
    physical_time += t_inc;
std::cout<<"physical_time is "<<physical_time<<std::endl;
	  MPI::COMM_WORLD.Barrier();
    UpdateLocalTmc(grid, t_inc);
	  MPI::COMM_WORLD.Barrier();
    UpdateLocalTmp(grid, physical_time, temp);
	  MPI::COMM_WORLD.Barrier();

		#ifndef SILENT
		if (rank==0) print_progress(step+1, steps, iterations);
		#endif
		update_timer += rdtsc()-start;
	}//loop over step
	#ifndef SILENT
	++iterations;
	#endif

  for(int i=0; i<nthreads; i++){
    delete [] num_of_grids_to_flip[i];
    num_of_grids_to_flip[i]=NULL;
    delete [] cell_coord[i];
    cell_coord[i]=NULL;
  }
  delete num_of_grids_to_flip; 
  num_of_grids_to_flip=NULL; 
  delete cell_coord;
  cell_coord=NULL;
//	delete [] grain_orientations;
//	grain_orientations=NULL;
//	delete [] temperature_along_x;
//	temperature_along_x=NULL;
	delete [] p_threads;
	p_threads=NULL;
	delete [] mat_para;
	mat_para=NULL;

	unsigned long total_update_time=update_timer;
	#ifdef MPI_VERSION
	MPI::COMM_WORLD.Allreduce(&update_timer, &total_update_time, 1, MPI_UNSIGNED_LONG, MPI_SUM);
  MPI::COMM_WORLD.Barrier();
	#endif
	return total_update_time/np; // average update time
}

}


#ifndef SILENT
void print_progress(const int step, const int steps, const int iterations)
{
	char* timestring;
	static unsigned long tstart;
	struct tm* timeinfo;

	if (step==0) {
		tstart = time(NULL);
		std::time_t rawtime;
		std::time( &rawtime );
		timeinfo = std::localtime( &rawtime );
		timestring = std::asctime(timeinfo);
		timestring[std::strlen(timestring)-1] = '\0';
		std::cout<<"Pass "<<std::setw(3)<<std::right<<iterations<<": "<<timestring<<" ["<<std::flush;
	} else if (step==steps) {
		unsigned long deltat = time(NULL)-tstart;
		std::cout << "•] "
							<<std::setw(2)<<std::right<<deltat/3600<<"h:"
							<<std::setw(2)<<std::right<<(deltat%3600)/60<<"m:"
							<<std::setw(2)<<std::right<<deltat%60<<"s"
							<<" (File "<<std::setw(5)<<std::right<<iterations*steps<<")."<<std::endl;
	} else if ((20 * step) % steps == 0) std::cout<<"• "<<std::flush;
}
#endif

#endif

#include"MMSP.main.hpp"

// Formatted using astyle:
//  astyle --style=linux --indent-col1-comments --indent=tab --indent-preprocessor --pad-header --align-pointer=type --keep-one-line-blocks --suffix=none
