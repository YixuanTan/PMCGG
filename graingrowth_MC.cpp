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
		const int edge = 64;
		int number_of_fields(seeds);
		if (number_of_fields==0) number_of_fields = static_cast<int>(float(edge*edge)/(M_PI*2.*2.)); // average grain is a disk of radius 10
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
		const int edge = 512;
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
  double* temperature_along_x; 
  EulerAngles* grain_orientations;
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

void ReadTemperature(double* temperature_along_x, int size){
  std::ifstream ifs("temperature.txt", std::ios::in); // sample_temperature.txt is a sample temperature file. temperature is assumed to be a function (only) of x coordinate. 1 st column means x coordinate and 2nd column means temperature. The 1st row of the file (current time is XXXX) is the time record from heater transfer simulation.
  std::vector<std::pair<double, double> > coords_and_temperatures;
  double x_coordinate, temperature;
  ifs.ignore(200, '\n');
  while(ifs>>x_coordinate>>temperature){//start from the second line
    coords_and_temperatures.push_back(std::make_pair(x_coordinate*(size-1), temperature));
  }
  ifs.close();
/*  for(unsigned int i=0; i<coords_and_temperatures.size(); i++){
    std::cout<<coords_and_temperatures[i].first<<"   "<<coords_and_temperatures[i].second<<std::endl;
  }*/
  int loop_temperatures_count=0;
  unsigned int i_last_step=0;
  unsigned int i=0;
  while(loop_temperatures_count<size){
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
}

double LargeAngleGrainBoundaryEnergy(const double temperature){//copper (100)|ND thin film texture
  return 0.625+(temperature-1198.0)*(-1.0e-4);
}

double StrainEnergyDenstiy(double h, double k, double l, double temperature){
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

template <int dim> void* flip_index_helper( void* s )
{
  double rotation_tolerance=0.5e-1;  // cos^-1(0.51)=59.34 deg   cos^-1(0.49)=60.66 deg
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
        site_out_of_domain = true;
        break;//break from the for int i loop
      }
    }
    if(site_out_of_domain == true){
      hh--;
      continue; //continue the int hh loop
    }
		int spin1 = (*(ss->grid))(x)%200;
		// determine neighboring spins
    vector<int> r(dim,0);
		sparse<bool> neighbors;
    if(dim==2){
		  for (int i=-1; i<=1; i++) {
			  for (int j=-1; j<=1; j++) {
				  r[0] = x[0] + i;
				  r[1] = x[1] + j;
				  int spin = (*(ss->grid))(r)%200;
				  set(neighbors,spin) = true;
        }
			}
    }else if(dim==3){
		  for (int i=-1; i<=1; i++){
			  for (int j=-1; j<=1; j++){
			    for (int k=-1; k<=1; k++) {
				    r[0] = x[0] + i;
				    r[1] = x[1] + j;
				    r[2] = x[2] + k;
				    int spin = (*(ss->grid))(r);
				    set(neighbors,spin) = true;
			    }
        }
		  }
    }

    //check if inside a grain
    int number_of_same_neighours = 0;
    for(int j=0; j<length(neighbors)-1; j++){
		  if(index(neighbors,j)!=index(neighbors,j+1)) break; //not inside a grain, break from for int j
      number_of_same_neighours++;
    }
    if(number_of_same_neighours==length(neighbors)-1){//inside a grain
      continue;//continue for 
    }

		// choose a random neighbor spin
		int spin2 = index(neighbors,rand()%length(neighbors))%200;
int rank = MPI::COMM_WORLD.Get_rank();
		// choose a random spin from Q states
//    int spin2 = rand()%200;
		if (spin1!=spin2){
			// compute energy change
			double dE = 0.0;
      if(dim==2){
        double film_thickness = 1.0e-6;
        double temperature=((ss->temperature_along_x))[x[0]];

        double h2=sin(((ss->grain_orientations)[spin2]).psi)*sin(((ss->grain_orientations)[spin2]).phi_two); 
        double k2=sin(((ss->grain_orientations)[spin2]).psi)*cos(((ss->grain_orientations)[spin2]).phi_two);  
        double l2=cos(((ss->grain_orientations)[spin2]).psi);
        double h1=sin(((ss->grain_orientations)[spin1]).psi)*sin(((ss->grain_orientations)[spin1]).phi_two);   
        double k1=sin(((ss->grain_orientations)[spin1]).psi)*cos(((ss->grain_orientations)[spin1]).phi_two);  
        double l1=cos(((ss->grain_orientations)[spin1]).psi);
//	  	  dE += 2*( SurfaceEnergy(h2, k2, l2, temperature)-SurfaceEnergy(h1, k1, l1, temperature)) //surface and interface energy
//            + film_thickness*(StrainEnergyDenstiy(h2, k2, l2, temperature)-StrainEnergyDenstiy(h1, k1, l1, temperature));//elastic strain energy
			  for (int i=-1; i<=1; i++){
				  for (int j=-1; j<=1; j++){
					  r[0] = x[0] + i;
					  r[1] = x[1] + j;
					  int spin = (*(ss->grid))(r)%200;
                dE += 0.5*((spin!=spin2)-(spin!=spin1)); // grain boundary energy  
 //             std::cout<<"dE is "<<dE<<" dE is "<<dE/unit_grain_boundary_area<<"\n";
				  }
        }
//        if(rank==0 && dE<0) std::cout<<"dE is "<<dE<<std::endl;
      }
      if(dim==3){
			  for (int i=-1; i<=1; i++){ 
				  for (int j=-1; j<=1; j++){ 
    	      for (int k=-1; k<=1; k++){ 
					    r[0] = x[0] + i;
					    r[1] = x[1] + j;
					    r[2] = x[2] + k;
					    int spin = (*(ss->grid))(r);
					    dE += (spin!=spin2)-(spin!=spin1);
				    }
          }
        }
      }
			// attempt a spin flip
			double r = double(rand())/double(RAND_MAX);
      kT = 1.3806488e-23*((ss->temperature_along_x))[x[0]];

			if (dE <= 0.0) {
        (*(ss->grid))(x) = spin2;
//	      int rank = MPI::COMM_WORLD.Get_rank();
//        if(rank==0 && dE<0) std::cout<<"dE is "<<dE<<std::endl;
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

template <int dim> unsigned long update(MMSP::grid<dim, int>& grid, int steps, int nthreads)
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

  EulerAngles *grain_orientations = new EulerAngles[200];
  if(rank==0){
    ReadData(grain_orientations);
/*    double phi, phi_two;
    for(int i=0; i<200; i++){
      psi=grain_orientations[i].k;
      phi_two=grain_orientations[i].l;

      h=sin(grain_orientations[i].psi)*sin(grain_orientations[i].phi_two);
      k=sin(grain_orientations[i].psi)*cos(grain_orientations[i].phi_two);
      l=cos(grain_orientations[i].psi);
//    std::cout<<grain_orientations[i].h<<" "<<grain_orientations[i].k<<" "<<grain_orientations[i].l<<std::endl;
    }*/
  }
	MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Bcast(grain_orientations, 200*3, MPI_DOUBLE, 0);

/*-----------------------------------------------*/
/*---------------generate cells------------------*/
/*-----------------------------------------------*/
  int dimension_length=0, number_of_lattice_cells=1;
  int lattice_cells_each_dimension[dim];
  for(int i=0; i<dim; i++){
    dimension_length = x1(grid, i)-x0(grid, i);
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
  int size=(g1(grid, 0)-g0(grid, 0)+1);
  double *temperature_along_x = new double[size];
  if(rank==0){
    ReadTemperature(temperature_along_x, size);
  }
	MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Bcast(temperature_along_x, size, MPI_DOUBLE, 0);
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

    mat_para[i].temperature_along_x = temperature_along_x;
    mat_para[i].grain_orientations = grain_orientations;

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
				pthread_create(&p_threads[i], &attr, flip_index_helper<dim>, (void*) &mat_para[i] );

			}//loop over threads

			for (int ii=0; ii!= nthreads ; ii++)
				pthread_join(p_threads[ii], NULL);

			#ifdef MPI_VERSION
			MPI::COMM_WORLD.Barrier();
			#endif

//			ghostswap(grid, sublattice); // once looped over a "color", ghostswap.
			ghostswap(grid); // once looped over a "color", ghostswap.
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
	delete [] grain_orientations;
	grain_orientations=NULL;
	delete [] temperature_along_x;
	temperature_along_x=NULL;
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
