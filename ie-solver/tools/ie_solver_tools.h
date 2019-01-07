#ifndef _IE_SOLVER_TOOLS_H_
#define _IE_SOLVER_TOOLS_H_

#include "ie-solver/ie_mat.h"
#include "ie-solver/quadtree.h"
#include "ie-solver/kernel.h"

#define MIN_DOFS_TO_COMPRESS 16

namespace ie_solver{

class IeSolverTools {

	public:

		double id_tol;
		bool is_stokes, strong_admissibility;

		IeSolverTools(){}
		IeSolverTools(double id_tol, bool strong_admissibility,
			bool is_stokes);
		~IeSolverTools(){}

		// The actual work functions
		void get_x_matrices(ie_Mat& K, const ie_Mat& Z, ie_Mat&,
			const std::vector<unsigned int>& r, 
			const std::vector<unsigned int>& s,
			const std::vector<unsigned int>& n);
		
		void schur_update(const Kernel& K, const ie_Mat& Z, ie_Mat& L, ie_Mat& U, 
			QuadTreeNode* node);
		int interpolative_decomposition(const Kernel& K, QuadTree& tree,
			QuadTreeNode* node);
		void make_id_mat(const Kernel& K, ie_Mat& pxy, QuadTree& tree, 
			QuadTreeNode* node);

		
		void get_all_schur_updates( ie_Mat& updates, 
			const std::vector<unsigned int>& BN, const QuadTreeNode*);
		void get_descendents_updates(ie_Mat& updates, 
			const std::vector<unsigned int>& BN, const QuadTreeNode* node);
		void get_update( ie_Mat& updates, const std::vector<unsigned int>& BN, 
			const QuadTreeNode* node);
			
		void apply_sweep_matrix(const ie_Mat& mat, ie_Mat& vec, 
			const std::vector<unsigned int>& a, 
			const std::vector<unsigned int>& b, bool transpose);
	    void apply_diag_matrix(const ie_Mat& mat, ie_Mat& vec, 
			const std::vector<unsigned int>& range);
		void apply_diag_inv_matrix(const ie_Mat& mat, ie_Mat& vec, 
			const std::vector<unsigned int>& range);

		void skeletonize(const Kernel& K, QuadTree& tree);

		void sparse_matvec(const Kernel& K, QuadTree& tree, const ie_Mat& x, 
			ie_Mat& b);
		void solve(const Kernel& K, QuadTree& tree, ie_Mat& x, const ie_Mat& b);


		void populate_active_box(QuadTreeNode*);
		// TODO when generalizing, let someone pass these functions in initing
		// IeSolverTools
		void make_proxy_mat(ie_Mat& pxy, double cntr_x, double cntr_y, double r, 
			QuadTree& tree, const std::vector<unsigned int>& box_indices);
		void make_stokes_proxy_mat(ie_Mat& pxy, double cntr_x, double cntr_y, 
			double r, const std::vector<unsigned int>& box_indices);

		void set_rs_ranges(InteractionLists& interaction_lists, const std::vector<unsigned int>& prm, 
			unsigned int sk, unsigned int rd);
		void set_skelnear_range(InteractionLists& interaction_lists);
		// void Conjugate_Gradient(ie_Mat& F, QuadTree& tree, ie_Mat& phi, 
	    // 	ie_Mat& f);

};

} // namespace ie_solver

#endif
