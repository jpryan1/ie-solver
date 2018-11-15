#ifndef _IE_SOLVER_TOOLS_H_
#define _IE_SOLVER_TOOLS_H_
#include <mutex>
#include <omp.h>
#include "common.h"
#include "quadtree.h"
#include <iostream>
#include <cstdarg>
#include <stdlib.h>
#include <stdio.h>

#define MIN_DOFS_TO_COMPRESS 16

namespace ie_solver{

class IeSolverTools {

	public:

		Clock set_get;
		double id_tol;
		//necessary for proxy creation
		std::vector<double> points, normals, weights; 
		std::mutex lock;
		bool is_stokes;

		IeSolverTools(){}
		IeSolverTools(double id_tol, const std::vector<double> points, 
			const std::vector<double> normals, 
			const std::vector<double> weights, bool is_stokes);
		~IeSolverTools(){}

		// The actual work functions
		void GetXMatrices(ie_Mat& K, const ie_Mat& Z, ie_Mat&,
			const std::vector<unsigned int>& r, 
			const std::vector<unsigned int>& s,
			const std::vector<unsigned int>& n);
		
		void SchurUpdate(const ie_Mat& K, const ie_Mat& Z, ie_Mat& L, ie_Mat& U, 
			QuadTreeNode* node);
		int InterpolativeDecomposition(const ie_Mat& K, ie_Mat& Z, 
			QuadTreeNode* node);
		
		void get_all_schur_updates( ie_Mat& updates, 
			const std::vector<unsigned int>& BN, const QuadTreeNode*);
		void get_descendents_updates(ie_Mat& updates, 
			const std::vector<unsigned int>& BN, const QuadTreeNode* node);
		void get_update( ie_Mat& updates, const std::vector<unsigned int>& BN, 
			const QuadTreeNode* node);
			
		void ApplySweepMatrix(const ie_Mat& mat, ie_Mat& vec, 
			const std::vector<unsigned int>& a, 
			const std::vector<unsigned int>& b, bool transpose);
	    void ApplyDiagMatrix(const ie_Mat& mat, ie_Mat& vec, 
			const std::vector<unsigned int>& range);
		void ApplyDiagInvMatrix(const ie_Mat& mat, ie_Mat& vec, 
			const std::vector<unsigned int>& range);

		void Skeletonize(const ie_Mat& K, QuadTree& tree);

		void SparseMatVec(const ie_Mat& K, QuadTree& tree, const ie_Mat& x, 
			ie_Mat& b);
		void Solve(const ie_Mat& K, QuadTree& tree, ie_Mat& x, const ie_Mat& b);


		void remove_inactive_dofs(QuadTreeNode*);
		// TODO when generalizing, let someone pass these functions in initing
		// IeSolverTools
		void make_proxy_mat(ie_Mat& pxy, double cntr_x, double cntr_y, double r, 
			const std::vector<unsigned int>& range);
		void make_stokes_proxy_mat(ie_Mat& pxy, double cntr_x, double cntr_y, 
			double r, const std::vector<unsigned int>& range);

		void set_rs_ranges(InteractionLists& interaction_lists, const std::vector<unsigned int>& prm, 
			unsigned int sk, unsigned int rd);
		void set_skelnear_range(InteractionLists& interaction_lists);
		// void Conjugate_Gradient(ie_Mat& F, QuadTree& tree, ie_Mat& phi, 
	    // 	ie_Mat& f);

};

} // namespace ie_solver

#endif
