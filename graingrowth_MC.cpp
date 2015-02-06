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

double lambda = 1e-4/1024;
double L_initial = 1.2e-6; 
double L0 = 1.1e-6;
double K1 = 0.12336665;
double n1 = 0.52822367;
double Q = 1.413e5; //fitted from Gangulee, A. ”Structure of electroplated and vapordeposited copper films. III. Recrystallization and grain growth.” Journal of Applied Physics 45.9 (1974): 3749-3756.
double n = 2; 
double K_ = 3.554e-5; //fitted from Gangulee, A. ”Structure of electroplated and vapordeposited copper films. III. Recrystallization and grain growth.” Journal of Applied Physics 45.9 (1974): 3749-3756.
double R = 8.314;

struct EulerAngles{
  double phi_one;
  double psi;
  double phi_two;  
}; 

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
		const int edge = 1024;
		int number_of_fields(seeds);
		if (number_of_fields==0) number_of_fields = static_cast<int>(float(edge*edge)/(M_PI*5*5)); /* average grain is a disk of radius XXX
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
		const int edge = 1024;
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
        (*grid).AccessToTmp(coords) = 273.0;
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
          (*grid).AccessToTmp(coords) = 273.0;
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
//  double* temperature_along_x; 
//  EulerAngles* grain_orientations;
  double t_mcs_max;
//  double t_s;
};

double MaxOfThreeNumber(const double h, const double k, const double l){
  double max_of_hkl = fabs(h)>fabs(k)?fabs(h):fabs(k);
  max_of_hkl = fabs(l)>max_of_hkl?fabs(l):max_of_hkl;
  return max_of_hkl;
}

double MediumOfThreeNumber(const double h, const double k, const double l){
  double medium = (fabs(h)+fabs(k)+fabs(l))/3;
  double medium_of_hkl = fabs(fabs(h)-medium)<fabs(fabs(k)-medium)?fabs(h):fabs(k);
  medium_of_hkl = fabs(fabs(l)-medium)<fabs(medium_of_hkl-medium)?fabs(l):medium_of_hkl;
  return medium_of_hkl;  
}

double MinOfThreeNumber(const double h, const double k, const double l){
  double min_of_hkl = fabs(h)<fabs(k)?fabs(h):fabs(k);
  min_of_hkl = fabs(l)<min_of_hkl?fabs(l):min_of_hkl;
  return min_of_hkl;
}

double SurfaceEnergy(const double h, const double k, const double l, const double temperature){//copper (100)|ND thin film texture
  // A.J.W Moore, Metal surfaces: structure, energetics and kinetic, charpter 5, thermal faceting, Page 166-169. American society for metals, Metals Park, Ohio, 1963. 
  // surface energy changing rate rate w.r.t temperature is set to -5.0e-4 for all hkl planes (From D.P. Field's paper (2005))
  double surface_energy = 0.0;
  double cosine_one = MaxOfThreeNumber(h, k, l)/sqrt(h*h+k*k+l*l); //  {100} w.r.t ND
  if(cosine_one>0.907){// choose {100} as the low index plane.  note in D.P. Field's paper (2005), only this low index plane is considered.
    surface_energy = 2.8465/cosine_one + temperature*(-5.0e-4); //2.610-473.0*(-5.0e-4) = 2.8465
  }else{
    double cosine_two = (fabs(h)+fabs(k)+fabs(l))/sqrt(h*h+k*k+l*l)/sqrt(3.0); // {111} w.r.t ND
    if(cosine_two>0.835){// choose {111} as the low index plane. 
      surface_energy = 2.7085/cosine_two + temperature*(-5.0e-4); //2.472-473.0*(-5.0e-4) = 2.7085
    }else{
      double cosine_three = (MaxOfThreeNumber(h,k,l) + MediumOfThreeNumber(h,k,l))/sqrt(h*h+k*k+l*l)/sqrt(2.0);//{110} w.r.t ND
      surface_energy = 3.0885/cosine_three + temperature*(-5.0e-4); //2.852-473.0*(-5.0e-4) = 3.0885
    } 
  }
  return surface_energy;
}

void ReadTemperature(double* temperature_along_x, const int model_dimension){
  std::ifstream ifs("temperature.txt", std::ios::in); // sample_temperature.txt is a sample temperature file. temperature is assumed to be a function (only) of x coordinate. 1 st column means x coordinate and 2nd column means temperature. The 1st row of the file (current time is XXXX) is the time record from heater transfer simulation.
  std::vector<std::pair<double, double> > coords_and_temperatures;
  double x_coordinate, temperature;
  ifs.ignore(200, '\n');
  while(ifs>>x_coordinate>>temperature){//start from the second line
    coords_and_temperatures.push_back(std::make_pair(x_coordinate*(model_dimension-1), temperature));
  }
  ifs.close();
/*  for(unsigned int i=0; i<coords_and_temperatures.size(); i++){
    std::cout<<coords_and_temperatures[i].first<<"   "<<coords_and_temperatures[i].second<<std::endl;
  }*/
  int loop_temperatures_count=0;
  unsigned int i_last_step=0;
  unsigned int i=0;
  while(loop_temperatures_count<model_dimension){
    i=i_last_step;
    while(i<coords_and_temperatures.size()){
      if(i!=coords_and_temperatures.size()-1 && loop_temperatures_count>(coords_and_temperatures[i].first) && (coords_and_temperatures[i+1].first)>loop_temperatures_count){
	        temperature_along_x[loop_temperatures_count]=coords_and_temperatures[i].second+(coords_and_temperatures[i+1].second-coords_and_temperatures[i].second)/(coords_and_temperatures[i+1].first-coords_and_temperatures[i].first)*((1.0*loop_temperatures_count)-coords_and_temperatures[i].first);
        i_last_step = i;
        i++;
        break;
      }else if((1.0*loop_temperatures_count)==coords_and_temperatures[i].first){
	      temperature_along_x[loop_temperatures_count]=coords_and_temperatures[i].second;
        i_last_step = i;
        i++;
        break;
      }
      i++;
    }// while i
    loop_temperatures_count++;
  }
/*
  for(int i=0; i<loop_temperatures_count;i++){
    std::cout<<"x="<<i<<"  temp="<<temperature_along_x[i]<<std::endl;
  }  
*/
}

double LargeAngleGrainBoundaryEnergy(const double temperature){//copper (100)|ND thin film texture
  return 0.625+(temperature-1198.0)*(-1.0e-4);
}

double StrainEnergyDenstiy(const double h, const double k, const double l, const double temperature){
  double c11=170.2e9+(temperature-298)*(-0.0353e9);//Pa^-1  N.J. Park, D.P. Field, Predicting thickness dependent twin boundary formation in sputtered Cu films, Scripta Materialia, 54.6, 2006, pp. 999-1003;
  double c12=123.2e9+(temperature-298)*(-0.0153e9);//Pa^-1
  double c44=75.4e9+(temperature-298)*(-0.0277e9);//Pa^-1
  double K=(2*c44-c11+c12)*(h*h*k*k+k*k*l*l+l*l*h*h);
  double M=c11+c12+K-2*(c12-K)*(c12-K)/(c11-2*K);
  double deposite_temperature=298.0;
  double elastic_strain=(3e-6-17e-6)*(temperature-deposite_temperature);
  double strain_energy_density=elastic_strain*elastic_strain*M;
  return strain_energy_density;
}

double CoherentTwinBoundaryEnergy(const double temperature){
  double coherent_twin_boundary_energy = 0.498+(-1.0e-4)*(temperature-1223);
  std::cout<<"CoherentTwinBoundaryEnergy fuction has been called"<<std::endl;
  return coherent_twin_boundary_energy;
}

double IncoherentTwinBoundaryEnergy(const double temperature){
  double incoherent_twin_boundary_energy = 0.024+(-2.0e-5)*(temperature-1073);
//  std::cout<<"IncoherentTwinBoundaryEnergy fuction has been called"<<std::endl;
  return incoherent_twin_boundary_energy;
}
 
bool CheckSixtyDegreeRotation(const double h, const double k, const double l, const double hh, const double kk, const double ll){//hkl is the neighbour, hh kk ll is the site to be/after fliped
  bool SixtyDegreeRotationAboutTriones = false;
  double cosine_hkl_abc, cosine_hhkkll_abc, hn1, kn1, ln1, hn2, kn2, ln2, hkl_magnitude, hhkkll_magnitude, abc_magnitude, cosine_two_prjectors;
  const double rotation_tolerance=1.0e-2;  // cos^-1(0.51)=59.34 deg   cos^-1(0.49)=60.66 deg
  for(double a=-1.0; a<=1.0; a+=2.0){
    for(double b=-1.0; b<=1.0; b+=2.0){
      for(double c=-1.0; c<=1.0; c+=2.0){
        abc_magnitude = sqrt(a*a+b*b+c*c);
        hkl_magnitude = sqrt(h*h+k*k+l*l);
        hhkkll_magnitude = sqrt(hh*hh+kk*kk+ll*ll);
        cosine_hkl_abc = (h*a+k*b+l*c)/hkl_magnitude/abc_magnitude;
        cosine_hhkkll_abc = (hh*a+kk*b+ll*c)/hhkkll_magnitude/abc_magnitude;
        hn1 = h - hkl_magnitude*cosine_hkl_abc*a/abc_magnitude; // h, k, l preject to the equator plan as hn1 kn1 ln1
        kn1 = k - hkl_magnitude*cosine_hkl_abc*b/abc_magnitude;
        ln1 = l - hkl_magnitude*cosine_hkl_abc*c/abc_magnitude;
        hn2 = hh - hhkkll_magnitude*cosine_hhkkll_abc*a/abc_magnitude; // hh, kk, ll preject to the equator plan as hn1 kn1 ln1
        kn2 = kk - hhkkll_magnitude*cosine_hhkkll_abc*b/abc_magnitude;
        ln2 = ll - hhkkll_magnitude*cosine_hhkkll_abc*c/abc_magnitude;
        cosine_two_prjectors = (hn1*hn2+kn1*kn2+ln1*ln2)/sqrt(hn1*hn1+kn1*kn1+ln1*ln1)/sqrt(hn2*hn2+kn2*kn2+ln2*ln2);
        if(fabs(cosine_two_prjectors-0.5)<rotation_tolerance){
          SixtyDegreeRotationAboutTriones = true;
          return SixtyDegreeRotationAboutTriones;
        }
      }
    }
  }
  return SixtyDegreeRotationAboutTriones;
}

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
			    for (int k=-1; k<=1; k++) {
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
//if(rank==0) std::cout<<"num_of_points_to_flip is "<<ss->num_of_points_to_flip<<std::endl;
//getchar();
    double temperature=(*(ss->grid)).AccessToTmp(x);
//    double initial_physical_time=1.0/K_/exp(-Q/R/temperature)*(pow(L_initial,n)-pow(L0,n));
//    double t_mcs_initial = pow(L_initial/K1/lambda,1.0/n1);
    double t_mcs = (*(ss->grid)).AccessToTmc(x);
    double t_mcs_max = ss->t_mcs_max;
    double site_selection_probability = t_mcs/t_mcs_max;
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

void ReadData(EulerAngles *grain_orientations){
  EulerAngles euler_angles;
  std::ifstream ifs("RandomNew.td", std::ios::in);  
  int i=0;
  while(ifs>>euler_angles.phi_one>>euler_angles.psi>>euler_angles.phi_two){
    grain_orientations[i]=euler_angles; 
    i++;
  }
  ifs.close();
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
  int edge = 1024;
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

//			ghostswap(grid, sublattice); // once looped over a "color", ghostswap.
			ghostswap(grid); // once looped over a "color", ghostswap.
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

template <int dim> double TmcMax(MMSP::grid<dim, int>& grid, double *temp_at_max_tmc){
   double tmc_max = 0.0;
   vector<int> coords (dim,0);
   if(dim==2){
       for(int codx=x0(grid, 0); codx <= x1(grid, 0); codx++) 
         for(int cody=x0(grid, 1); cody <= x1(grid, 1); cody++){
           coords[0] = codx;
           coords[1] = cody;
/*           if(grid.AccessToTmc(coords) > 1000){
             std::cout<<coords[0]<<"  "<<coords[1]<<";   "<< grid.AccessToTmc(coords)  <<std::endl;
             getchar();
           }*/
           if(grid.AccessToTmc(coords) > tmc_max)
             tmc_max = grid.AccessToTmc(coords);
         }
   }
   else if(dim==3){
     for(int codx=x0(grid, 0); codx <= x1(grid, 0); codx++) 
       for(int cody=x0(grid, 1); cody <= x1(grid, 1); cody++) 
         for(int codz=x0(grid, 2); codz <= x1(grid, 2); codz++){
           coords[0] = codx;
           coords[1] = cody;
           coords[2] = codz;
           if(grid.AccessToTmc(coords) > tmc_max)
             tmc_max = grid.AccessToTmc(coords);
         }
   }
   (*temp_at_max_tmc) = grid.AccessToTmp(coords);
   return tmc_max;
}

template <int dim> void UpdateLocalTmc(MMSP::grid<dim, int>& grid, long double t_inc){
   vector<int> coords (dim,0);
   if(dim==2){
       for(int codx=x0(grid, 0); codx <= x1(grid, 0); codx++) 
         for(int cody=x0(grid, 1); cody <= x1(grid, 1); cody++){
           coords[0] = codx;
           coords[1] = cody;
           long double exp_eqn_rhs = pow(K1*lambda*pow(grid.AccessToTmc(coords), n1), n) - pow(L0, n);
           long double temperature = grid.AccessToTmp(coords);
           exp_eqn_rhs += K_*t_inc*exp(-Q/R/temperature);
           int rank=MPI::COMM_WORLD.Get_rank();
//if(rank==0) std::cout<<"t_inc is "<<t_inc<<"  temperature is "<<temperature<<"  before grid.AccessToTmc(coords)  is "<<grid.AccessToTmc(coords)<<std::endl;
           grid.AccessToTmc(coords) = pow(1.0/K1/lambda*pow(exp_eqn_rhs+pow(L0, n),1.0/n), 1.0/n1);  
//             if(rank==0) std::cout<<"after grid.AccessToTmc(coords)  is "<<grid.AccessToTmc(coords)<<std::endl;
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

template <int dim> void UpdateLocalTmp(MMSP::grid<dim, int>& grid, long double physical_time){
   vector<int> coords (dim,0);
   if(dim==2){
       for(int codx=x0(grid, 0); codx <= x1(grid, 0); codx++) 
         for(int cody=x0(grid, 1); cody <= x1(grid, 1); cody++){
           coords[0] = codx;
           coords[1] = cody;
/*-----------------------*/
if(codx<=0.25*1024)
  grid.AccessToTmp(coords) = 100; 
else if(0.25*1024<codx<=0.5*1024)
  grid.AccessToTmp(coords) = 300; 
else if(0.5*1024<codx<=0.75*1024)
  grid.AccessToTmp(coords) = 500; 
else if(0.75*1024<codx<=1024)
  grid.AccessToTmp(coords) = 700; 
/*-----------------------*/
//           grid.AccessToTmp(coords) = 473 + 273.0*sin(3.14/1024*codx + 3.14/(1.0e5)*physical_time);
         }
   }
   else if(dim==3){
     for(int codx=x0(grid, 0); codx <= x1(grid, 0); codx++)  
         for(int cody=x0(grid, 1); cody <= x1(grid, 1); cody++) 
           for(int codz=x0(grid, 2); codz <= x1(grid, 2); codz++){
             coords[0] = codx;
             coords[1] = cody;
             coords[2] = codz;
             grid.AccessToTmp(coords) = 473 + 273.0*sin(3.14/1024*codx + 3.14/(1.0e5)*physical_time);
           }
   }
}

template <int dim> unsigned long update(MMSP::grid<dim, int>& grid, int steps, int steps_finished, int nthreads, int step_to_nonuniform, long double &physical_time)
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
  int edge = 1024;
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

/* calculate tmc_max
  double tmc_max;
	#ifdef MPI_VERSION
	MPI::COMM_WORLD.Allreduce(&update_timer, &total_update_time, 1, MPI_UNSIGNED_LONG, MPI_SUM);
	#endif
 calculate tmc_max*/

    double temp_at_max_tmc = 0.0;
    double tmc_max_partition = TmcMax(grid, &temp_at_max_tmc); // tmc includes tmc_initial
    MPI::COMM_WORLD.Barrier();
    double tmc_max_global = 0.0;
//if(rank==0) std::cout<<"temp_at_max_tmc is "<<temp_at_max_tmc<<std::endl;
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Allreduce(&tmc_max_partition, &tmc_max_global, 1, MPI_DOUBLE, MPI_MAX);
    MPI::COMM_WORLD.Barrier();

//if(rank==0) std::cout<<"steps_finished is "<<steps_finished<<"step_to_nonuniform is "<<step_to_nonuniform<<std::endl;
//    std::cout<<"at rank "<<rank << "  tmc_max_partition is "<<tmc_max_partition<<"tmc_max_global is "<<tmc_max_global<<std::endl;

    if(steps_finished != step_to_nonuniform && tmc_max_partition == tmc_max_global){
      MPI::COMM_WORLD.Bcast(&temp_at_max_tmc, 1, MPI_DOUBLE, rank);
    }
    MPI::COMM_WORLD.Barrier();
//    if(rank==2) std::cout<<"temp_at_max_tmc is "<<temp_at_max_tmc<<std::endl;
    double t_inc = ( pow(K1*lambda*pow(tmc_max_global+1,n1), n) - pow(K1*lambda*pow(tmc_max_global,n1), n) )/K_/exp(-Q/R/temp_at_max_tmc);

		unsigned long start = rdtsc();
    int num_of_sublattices=0;
    if(dim==2) num_of_sublattices = 4; 
    else if(dim==3) num_of_sublattices = 8;
		for (int sublattice=0; sublattice < num_of_sublattices; sublattice++) {
			for (int i=0; i!= nthreads ; i++) {
				mat_para[i].sublattice=sublattice;
				mat_para[i].num_of_points_to_flip=num_of_grids_to_flip[i][sublattice];
        mat_para[i].t_mcs_max = tmc_max_global; //tmc_max is the MC time steps counted from the beginning of simulation
        for(int k=0; k<dim; k++) mat_para[i].cell_coord[k]=cell_coord[i][k];
				pthread_create(&p_threads[i], &attr, flip_index_helper<dim>, (void*) &mat_para[i] );
			}//loop over threads

			for (int ii=0; ii!= nthreads ; ii++)
				pthread_join(p_threads[ii], NULL);

			#ifdef MPI_VERSION
			MPI::COMM_WORLD.Barrier();
			#endif

//			ghostswap(grid, sublattice); // once looped over a "color", ghostswap.
			ghostswap(grid); // once looped over a "color", ghostswap.
            #ifdef MPI_VERSION
			MPI::COMM_WORLD.Barrier();
                #endif
		}//loop over color
    //after 1 global tmc, update all the local tmc
	  MPI::COMM_WORLD.Barrier();
    UpdateLocalTmc(grid, t_inc);
    physical_time += t_inc;
	  MPI::COMM_WORLD.Barrier();
//    if(rank==0) std::cout<<"physical_time is "<< physical_time << ";  t_inc is" << t_inc << std::endl; 
    UpdateLocalTmp(grid, physical_time);
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
