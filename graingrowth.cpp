// graingrowth.cpp
// Coarsening algorithms for 2D and 3D sparse phase field (sparsePF) methods
// Questions/comments to kellet@rpi.edu (Trevor Keller)

#ifndef GRAINGROWTH_UPDATE
#define GRAINGROWTH_UPDATE

#include <iomanip>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include"graingrowth.hpp"
#include"MMSP.hpp"
#include"tessellate.hpp"
#include"output.cpp"

#ifndef SILENT
void print_progress(const int step, const int steps, const int iterations);
#endif

namespace MMSP {

template <int dim>
unsigned long generate(MMSP::grid<dim,MMSP::sparse<float> >*& grid, int seeds, int nthreads)
{
	#if (defined CCNI) && (!defined MPI_VERSION)
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	int rank=0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	int np = MPI::COMM_WORLD.Get_size();
	#endif
	unsigned long timer=0;
	if (dim == 2) {
		const int edge = 128;
		int number_of_fields(seeds);
		if (number_of_fields==0) number_of_fields = static_cast<int>(float(edge*edge)/(M_PI*10.*10.)); // average grain is a disk of radius 10
		#ifdef MPI_VERSION
		while (number_of_fields % np) --number_of_fields;
		#endif
		grid = new MMSP::grid<dim,MMSP::sparse<float> >(0, 0, edge, 0, edge);
		#ifndef SILENT
		if (rank==0) std::cout<<"Grid origin: ("<<g0(*grid,0)<<','<<g0(*grid,1)<<"),"
												<<" dimensions: "<<g1(*grid,0)-g0(*grid,0)<<" × "<<g1(*grid,1)-g0(*grid,1)
												<<" with "<<number_of_fields<<" grains."<<std::endl;
		#endif
		#ifdef MPI_VERSION
		number_of_fields /= np;
		#endif

		#if (!defined MPI_VERSION) && ((defined CCNI) || (defined BGQ))
		std::cerr<<"Error: CCNI requires MPI."<<std::endl;
		std::exit(1);
		#endif
		timer = tessellate<dim,float>(*grid, number_of_fields, nthreads);
		#ifndef SILENT
		if (rank==0) std::cout<<"Tessellation complete."<<std::endl;
		#endif
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
	} else if (dim == 3) {
		const int edge = 64;
		int number_of_fields(seeds);
		if (number_of_fields==0) number_of_fields = static_cast<int>(float(edge*edge*edge)/(4./3*M_PI*10.*10.*10.)); // Average grain is a sphere of radius 10 voxels
		#ifdef MPI_VERSION
		while (number_of_fields % np) --number_of_fields;
		#endif
		grid = new MMSP::grid<dim,MMSP::sparse<float> >(0,0,edge,0,edge,0,edge);
		#ifndef SILENT
		if (rank==0) std::cout<<"Grid origin: ("<<g0(*grid,0)<<','<<g0(*grid,1)<<','<<g0(*grid,2)<<"),"
												<<" dimensions: "<<g1(*grid,0)-g0(*grid,0)<<" × "<<g1(*grid,1)-g0(*grid,1)<<" × "<<g1(*grid,2)-g0(*grid,2)
												<<" with "<<number_of_fields<<" grains."<<std::endl;
		#endif
		#ifdef MPI_VERSION
		number_of_fields /= np;
		#endif

		#if (!defined MPI_VERSION) && ((defined CCNI) || (defined BGQ))
		std::cerr<<"Error: CCNI requires MPI."<<std::endl;
		std::exit(1);
		#endif
		timer = tessellate<dim,float>(*grid, number_of_fields, nthreads);
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
		#ifndef SILENT
		if (rank==0) std::cout<<"Tessellation complete."<<std::endl;
		#endif
	}
	return timer;
}


unsigned long generate(int dim, char* filename, int seeds, int nthreads) {
	#if (defined CCNI) && (!defined MPI_VERSION)
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	unsigned long timer = 0;
	int rank=0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif
	if (dim == 2) {
		MMSP::grid<2,MMSP::sparse<float> >* grid2=NULL;
		timer=generate<2>(grid2,seeds,nthreads);
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
		MMSP::grid<3,MMSP::sparse<float> >* grid3=NULL;
		timer=generate<3>(grid3,seeds,nthreads);
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
struct update_thread_para {
	MMSP::grid<dim,sparse<float> >* grid;
	MMSP::grid<dim,sparse<float> >* update;
	unsigned long nstart;
	unsigned long nend;
};

template <int dim>
void* update_threads_helper( void * s )
{
	update_thread_para<dim>* ss = ( update_thread_para<dim>* ) s ;

	const float dt = 0.01;
	const float width = 10.0;
	const float gamma = 1.0;
	const float eps = 4.0 / acos(-1.0) * sqrt(0.5 * gamma * width);
	const float w = 4.0 * gamma / width;
	const float mu = 1.0;
	const float epsilon = 1.0e-8;

	for (unsigned int i = ss->nstart; i < ss->nend; i++) {
		vector<int> x = position((*ss->grid), i);

		// determine nonzero fields within
		// the neighborhood of this node
		// (2 adjacent voxels along each cardinal direction)
		sparse<int> s;
		for (int j = 0; j < dim; j++)
			for (int k = -1; k <= 1; k++) {
				x[j] += k;
				for (int h = 0; h < length((*ss->grid)(x)); h++) {
					int index = MMSP::index((*ss->grid)(x), h);
					set(s, index) = 1;
				}
				x[j] -= k;
			}
		float S = float(length(s));

		// if only one field is nonzero,
		// then copy this node to update
		if (S < 2.0) (*ss->update)(i) = (*ss->grid)(i);
		else {
			// compute laplacian of each field
			sparse<float> lap = laplacian((*ss->grid), i);

			// compute variational derivatives
			sparse<float> dFdp;
			for (int h = 0; h < length(s); h++) {
				int hindex = MMSP::index(s, h);
				for (int j = h + 1; j < length(s); j++) {
					int jindex = MMSP::index(s, j);
					// Update dFdp_h and dFdp_j, so the inner loop can be over j>h instead of j≠h
					set(dFdp, hindex) += 0.5 * eps * eps * lap[jindex] + w * (*ss->grid)(i)[jindex];
					set(dFdp, jindex) += 0.5 * eps * eps * lap[hindex] + w * (*ss->grid)(i)[hindex];
				}
			}

			// compute time derivatives
			sparse<float> dpdt;
			for (int h = 0; h < length(s); h++) {
				int hindex = MMSP::index(s, h);
				for (int j = h + 1; j < length(s); j++) {
					int jindex = MMSP::index(s, j);
					set(dpdt, hindex) -= mu * (dFdp[hindex] - dFdp[jindex]);
					set(dpdt, jindex) -= mu * (dFdp[jindex] - dFdp[hindex]);
				}
			}

			// compute update values
			float sum = 0.0;
			for (int h = 0; h < length(s); h++) {
				int index = MMSP::index(s, h);
				float value = (*ss->grid)(i)[index] + dt * (2.0 / S) * dpdt[index]; // Extraneous factor of 2?
				if (value > 1.0) value = 1.0;
				if (value < 0.0) value = 0.0;
				if (value > epsilon) set((*ss->update)(i), index) = value;
				sum += (*ss->update)(i)[index];
			}

			// project onto Gibbs simplex (enforce Σφ=1)
			float rsum = 0.0;
			if (fabs(sum) > 0.0) rsum = 1.0 / sum;
			for (int h = 0; h < length((*ss->update)(i)); h++) {
				int index = MMSP::index((*ss->update)(i), h);
				set((*ss->update)(i), index) *= rsum;
			}
		}
	} // Loop over nodes(grid)


	pthread_exit(0);
	return NULL;
}


template <int dim>
unsigned long update(MMSP::grid<dim, sparse<float> >& grid, int steps, int nthreads)
{

    #if (!defined MPI_VERSION) && ((defined CCNI) || (defined BGQ))
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	int rank=0;
	unsigned int np=0;
	#ifdef MPI_VERSION
 	rank=MPI::COMM_WORLD.Get_rank();
 	np=MPI::COMM_WORLD.Get_size();
	#endif

	#ifndef SILENT
	static int iterations = 1;
 	if (rank==0) print_progress(0, steps, iterations);
 	#endif

	unsigned long timer = 0;
	for (int step = 0; step < steps; step++) {
		unsigned long comptime=rdtsc();
		// update grid must be overwritten each time
		MMSP::grid<dim, sparse<float> > update(grid);
		ghostswap(grid);

		pthread_t * p_threads = new pthread_t[ nthreads];
        update_thread_para<dim>* update_para = new update_thread_para<dim>[nthreads];
        pthread_attr_t attr;
        pthread_attr_init (&attr);
        unsigned long nincr = nodes(grid)/nthreads;
        unsigned long ns = 0;

        for(int i=0; i<nthreads; i++) {
    	    update_para[i].nstart=ns;
            ns+=nincr;
            update_para[i].nend=ns;

            update_para[i].grid= &grid;
            update_para[i].update= &update;

            pthread_create(&p_threads[i], &attr, update_threads_helper<dim>, (void *) &update_para[i] );
        }

        for(int i=0; i!= nthreads ; i++) {
            pthread_join(p_threads[i], NULL);
        }

        delete [] p_threads ;
		delete [] update_para ;

		#ifndef SILENT
		if (rank==0) print_progress(step+1, steps, iterations);
		#endif
		swap(grid, update);
		timer += rdtsc() - comptime;
	} // Loop over steps
	ghostswap(grid);
	#ifndef SILENT
	++iterations;
	#endif

	unsigned long total_update_time=timer;
	#ifdef MPI_VERSION
	MPI::COMM_WORLD.Allreduce(&timer, &total_update_time, 1, MPI_UNSIGNED_LONG, MPI_SUM);
	#endif
	return total_update_time/np; // average update time
}

template <int dim>
unsigned long update_old(MMSP::grid<dim, sparse<float> >& grid, int steps, int nthreads)
{
	#if (!defined MPI_VERSION) && ((defined CCNI) || (defined BGQ))
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	unsigned long timer=0;
	int rank=0;
	#ifdef MPI_VERSION
 	rank=MPI::COMM_WORLD.Get_rank();
	#endif
	const float dt = 0.01;
	const float width = 10.0;
	const float gamma = 1.0;
	const float eps = 4.0 / acos(-1.0) * sqrt(0.5 * gamma * width);
	const float w = 4.0 * gamma / width;
	const float mu = 1.0;
	const float epsilon = 1.0e-8;

	#ifndef SILENT
	static int iterations = 1;
	#ifdef DEBUG
	if (iterations==1 && rank==0)
		printf("CFL condition Co=%2.2f.\n", mu*eps*eps*dt/(dx(grid, 0)*dx(grid,0)));
	#endif
	#endif

	#ifndef SILENT
 	if (rank==0) print_progress(0, steps, iterations);
 	#endif

	for (int step = 0; step < steps; step++) {
		// update grid must be overwritten each time
		MMSP::grid<dim, sparse<float> > update(grid);
		ghostswap(grid);

		unsigned long comp_time=rdtsc();
		for (int i = 0; i < nodes(grid); i++) {
			vector<int> x = position(grid, i);

			// determine nonzero fields within
			// the neighborhood of this node
			// (2 adjacent voxels along each cardinal direction)
			sparse<int> s;
			for (int j = 0; j < dim; j++)
				for (int k = -1; k <= 1; k++) {
					x[j] += k;
					for (int h = 0; h < length(grid(x)); h++) {
						int index = MMSP::index(grid(x), h);
						set(s, index) = 1;
					}
					x[j] -= k;
				}
			float S = float(length(s));

			// if only one field is nonzero,
			// then copy this node to update
			if (S < 2.0) update(i) = grid(i);
			else {
				// compute laplacian of each field
				sparse<float> lap = laplacian(grid, i);

				// compute variational derivatives
				sparse<float> dFdp;
				for (int h = 0; h < length(s); h++) {
					int hindex = MMSP::index(s, h);
					for (int j = h + 1; j < length(s); j++) {
						int jindex = MMSP::index(s, j);
						// Update dFdp_h and dFdp_j, so the inner loop can be over j>h instead of j≠h
						set(dFdp, hindex) += 0.5 * eps * eps * lap[jindex] + w * grid(i)[jindex];
						set(dFdp, jindex) += 0.5 * eps * eps * lap[hindex] + w * grid(i)[hindex];
					}
				}

				// compute time derivatives
				sparse<float> dpdt;
				for (int h = 0; h < length(s); h++) {
					int hindex = MMSP::index(s, h);
					for (int j = h + 1; j < length(s); j++) {
						int jindex = MMSP::index(s, j);
						set(dpdt, hindex) -= mu * (dFdp[hindex] - dFdp[jindex]);
						set(dpdt, jindex) -= mu * (dFdp[jindex] - dFdp[hindex]);
					}
				}

				// compute update values
				float sum = 0.0;
				for (int h = 0; h < length(s); h++) {
					int index = MMSP::index(s, h);
					float value = grid(i)[index] + dt * (2.0 / S) * dpdt[index]; // Extraneous factor of 2?
					if (value > 1.0) value = 1.0;
					if (value < 0.0) value = 0.0;
					if (value > epsilon) set(update(i), index) = value;
					sum += update(i)[index];
				}

				// project onto Gibbs simplex (enforce Σφ=1)
				float rsum = 0.0;
				if (fabs(sum) > 0.0) rsum = 1.0 / sum;
				for (int h = 0; h < length(update(i)); h++) {
					int index = MMSP::index(update(i), h);
					set(update(i), index) *= rsum;
				}
			}
		} // Loop over nodes(grid)
		swap(grid, update);
		comp_time=rdtsc()-comp_time;
		timer+=comp_time;
		#ifndef SILENT
		if (rank==0) print_progress(step+1, steps, iterations);
		#endif
	} // Loop over steps
	ghostswap(grid);
	#ifndef SILENT
	++iterations;
	#endif
	return timer;
}

} // namespace MMSP

#ifndef SILENT
void print_progress(const int step, const int steps, const int iterations) {
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
//	astyle --style=linux --indent-col1-comments --indent=tab --indent-preprocessor --pad-header --align-pointer=type --keep-one-line-blocks --suffix=none
