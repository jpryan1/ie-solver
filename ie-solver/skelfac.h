#ifndef _SKELFAC_H_
#define _SKELFAC_H_
#include <mutex>
#include <omp.h>
#include "common.h"
#include "quadtree.h"
#include <iostream>
#include <cstdarg>
#include <stdlib.h>
#include <stdio.h>

namespace ie_solver{

class Skelfac {

	public:
		Skelfac(){}
		Skelfac(double id_tol, std::vector<double> points, 
			std::vector<double> normals, std::vector<double> weights, bool s) {
			this->id_tol = id_tol;
			this->points = points;
			this->normals = normals;
			this->weights = weights;
			is_stokes = s;
		}
		~Skelfac(){}

		// The actual work functions
		void GetXMatrices(ie_Mat&, ie_Mat&, ie_Mat&, std::vector<unsigned int>&, 
			std::vector<unsigned int>&, std::vector<unsigned int>&);
		void SchurUpdate(ie_Mat&, ie_Mat&, ie_Mat&, ie_Mat&, 
			QuadTreeNode* node);
		int InterpolativeDecomposition(ie_Mat&, ie_Mat& Z, QuadTreeNode*);
		
		void get_descendents_updates(ie_Mat&, std::vector<unsigned int>&, 
			QuadTreeNode*);
		void get_all_schur_updates(  ie_Mat&, std::vector<unsigned int>&, 
			QuadTreeNode*);
		void get_update( ie_Mat&, std::vector<unsigned int>&, QuadTreeNode*);
			
		void ApplySweepMatrix(ie_Mat&, ie_Mat&, std::vector<unsigned int>, 
			std::vector<unsigned int>, bool);
	    void ApplyDiagMatrix(ie_Mat& , ie_Mat&, std::vector<unsigned int>);
		void ApplyDiagInvMatrix(ie_Mat& , ie_Mat&, std::vector<unsigned int>);
		void SparseMatVec(ie_Mat&, QuadTree&, ie_Mat&, ie_Mat&);
	   
		void Solve(ie_Mat& F, QuadTree& tree, ie_Mat& x, ie_Mat& b);
	    void Skeletonize(ie_Mat& F, QuadTree& tree);
	    void Conjugate_Gradient(ie_Mat& F, QuadTree& tree, ie_Mat& phi, ie_Mat& f);
		void remove_inactive_dofs(QuadTreeNode*);
		void make_proxy_mat(ie_Mat& pxy, double cntr_x, double cntr_y, double r, 
			std::vector<unsigned int>& range);
		void make_stokes_proxy_mat(ie_Mat& pxy, double cntr_x, double cntr_y, 
			double r, std::vector<unsigned int>& range);

		Clock set_get;

    	/**
		 * The parameters to be used in computation
		 */
		double id_tol = 0.001;
		
		/**
		 * How much should be outputted
		 */
		int verbosity = 0;

		std::vector<double> points, normals, weights; //necessary for proxy creation

        std::mutex lock;
	private:
		void set_rs_ranges(Box&, std::vector<unsigned int>&, unsigned int, unsigned int);
		void set_skelnear_range(Box& box);
		bool is_stokes;
		/**
		 * Quick and easy logging function based on verbosity
		 */
		double gemm_time = 0;
};

} // namespace ie_solver

#endif
