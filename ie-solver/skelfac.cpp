#include "ie_solver_tools.h"

namespace ie_solver{

IeSolverTools::IeSolverTools(double id_tol, const std::vector<double> points, 
	const std::vector<double> normals, const std::vector<double> weights, 
	bool is_stokes_) {

	this->id_tol = id_tol;
	this->points = points;
	this->normals = normals;
	this->weights = weights;
	is_stokes = is_stokes_;
}

void IeSolverTools::GetXMatrices(ie_Mat& K, const ie_Mat& Z, ie_Mat& Xrr, 
	const std::vector<unsigned int>& r, const std::vector<unsigned int>& s, 
	const std::vector<unsigned int>& n) {
    
	unsigned int r_size = r.size();
	unsigned int s_size = s.size();
	unsigned int n_size = n.size();
    ie_Mat  Xrs(r_size, s_size), Xsr(s_size, r_size) Xrn(r_size, n_size), 
    Xnr(n_size, r_size);
	
	//this is just for readability
	Xrr = K(r,r);
	Xrs = K(r,s);
	Xsr = K(s,r);
	if(n_size > 0){
		Xrn = K(r,n);
		Xnr = K(n,r);
	}

	Matmul::ie_gemm(TRANSPOSE, NORMAL, -1., Z, 	  K(s,r), 1., Xrr);
	Matmul::ie_gemm(TRANSPOSE, NORMAL, -1., Z, 	  K(s,s), 1., Xrs);
	Matmul::ie_gemm(NORMAL, 	NORMAL, -1., Xrs, 	  Z, 	  1., Xrr);
	Matmul::ie_gemm(NORMAL, 	NORMAL, -1., K(s,s), Z, 	  1., Xsr);
	if(n_size > 0){
		Matmul::ie_gemm(TRANSPOSE, NORMAL, -1., Z, 	  K(s,n), 1., Xrn);
		Matmul::ie_gemm(NORMAL, 	NORMAL, -1., K(n,s), Z, 	  1., Xnr);
	}

	set_get.tic();
	K.set_submatrix(r, s, Xrs);
	//K.set_submatrix(r, r, Xrr);
	K.set_submatrix(s, r, Xsr);
	if(n_size > 0){
		K.set_submatrix(r, n, Xrn);
		K.set_submatrix(n, r, Xnr);
	}
	set_get.toc();
	// TODO this seems like an awful lot of stores, can we avoid this? 
}


// TODO These involved multiplying into a buffer, then copying the buffer into K
// Can we do this in place? Or use just one buffer and resize each time?
// Does it matter?
// Also, we might not need to do ALL of these matmuls, since some of the 
// blocks will be eliminated anyways. This should be worked out on a whiteboard
void IeSolverTools::SchurUpdate(const ie_Mat& K, const ie_Mat& Z, ie_Mat& L, 
	ie_Mat& U, QuadTreeNode* node) {
    //height of Z is number of skeleton columns
	unsigned int num_redundant = Z.width();
	unsigned int num_skel 	  = Z.height();
	unsigned int num_near 	  = node->interaction_lists.near.size();

	//GENERATE K_BN,BN
	std::vector<unsigned int> BN;
	for(unsigned int idx : node->interaction_lists.box)  BN.push_back(idx);
	for(unsigned int idx : node->interaction_lists.near) BN.push_back(idx);

	// Note that BN has all currently deactivated DoFs removed. 
	ie_Mat K_BN = K(BN, BN); //que bien!
	ie_Mat update(BN.size(), BN.size());
	get_all_schur_updates(update, BN, node);
	K_BN -= update;
	
	// Generate various index ranges within BN
	std::vector<unsigned int> s, r, n, sn;
	for(unsigned int i=0; i<num_skel; i++){
		s.push_back(node->interaction_lists.permutation[i]);
		sn.push_back(node->interaction_lists.permutation[i]);
	}
	
	for(unsigned int i=0; i<num_redundant; i++){
		r.push_back(node->interaction_lists.permutation[i+num_skel]);
	}
	// num_near may be zero. In this case, GetXMatrices must act 
	// accordingly
	for(unsigned int i=0; i<num_near; i++){
		n.push_back(i+num_redundant+num_skel);
		sn.push_back(i+num_redundant+num_skel);
	}

	ie_Mat  Xrr(num_redundant, num_redundant);
	GetXMatrices(K_BN, Z, Xrr, r,s,n);
	node->D_r = Xrr;

	// Generate left and right schur complement matrices
	
	L = ie_Mat(num_skel + num_near, num_redundant);
	U = ie_Mat(num_redundant, num_skel + num_near);
	
	Xrr.right_multiply_inverse(K_BN(sn, r), L);
	Xrr.left_multiply_inverse( K_BN(r, sn), U);

	ie_Mat schur(node->interaction_lists.skelnear.size(), 
		node->interaction_lists.skelnear.size());
	Matmul::ie_gemm(NORMAL, NORMAL, 1.0, L, 
		K_BN(r,sn), 0., schur);
	
	//set schur update
	node->schur_update = schur;
	node->schur_updated = true;
}


int IeSolverTools::InterpolativeDecomposition(const ie_Mat& K, ie_Mat& Z, 
	QuadTreeNode* node) {

	double cntr_x = node->corners[0] + node->side_length/2.0;
	double cntr_y = node->corners[1] + node->side_length/2.0;
	ie_Mat pxy;

	if(is_stokes){
		pxy = ie_Mat(200, node->interaction_lists.box.size());
		make_stokes_proxy_mat(pxy, cntr_x, cntr_y, node->side_length*2, 
			node->interaction_lists.box );
	}
	else{
		pxy = ie_Mat(100, node->interaction_lists.box.size());
		make_proxy_mat(pxy, cntr_x, cntr_y, node->side_length*2, 
			node->interaction_lists.box);
	}
	std::vector<unsigned int> p;
	unsigned int numskel = pxy.id(p, Z, id_tol);
	if(numskel==0) return 0;
	set_rs_ranges(node->interaction_lists, p, Z.height(), Z.width());
	
	set_skelnear_range(node->interaction_lists);
	
	return Z.width();
}


void IeSolverTools::get_all_schur_updates(ie_Mat& updates, 
	const std::vector<unsigned int>& BN, const QuadTreeNode* node){
	if(!node->is_leaf) get_descendents_updates(updates, BN, node);
	for(QuadTreeNode* neighbor : node->neighbors){
		if(neighbor->schur_updated) get_update(updates, BN, neighbor);
		if(!neighbor->is_leaf) get_descendents_updates(updates, BN, neighbor);
	}
	// TODO make sure the following is actually supposed to be commented out
	// for(QuadTreeNode* neighbor : node->far_neighbors){
	// 	if(neighbor->schur_updated) get_update(updates, BN, neighbor);
	// 	if(!neighbor->is_leaf) get_descendents_updates(updates, BN, neighbor);
	// }
}
	

void IeSolverTools::get_descendents_updates(ie_Mat& updates, 
	const std::vector<unsigned int>& BN, const QuadTreeNode* node){
	//by assumption, node is not a leaf	
	for(QuadTreeNode* child : node->children){
		if(child->schur_updated) get_update(updates, BN, child);
		if(!child->is_leaf) get_descendents_updates(updates, BN, child);
	}
}


void IeSolverTools::get_update(ie_Mat& update, const std::vector<unsigned int>& BN, 
	const QuadTreeNode* node){
	// node needs to check all its dofs against BN, enter interactions into 
	// corresponding locations
	// node only updated its own BN dofs, and the redundant ones are no longer 
	// relevant, so we only care about child's SN dofs
	std::vector<unsigned int> sn = node->interaction_lists.skelnear;
	// First create a list of Dofs that are also in node's skelnear, 
	// and with each one give the index in skelnear and the index in BN
	std::vector<unsigned int> BN_indices;
	std::vector<unsigned int> sn_indices;
	for(unsigned int i=0; i<sn.size(); i++){
		for(unsigned int j=0; j<BN.size(); j++){
			if(BN[j] == sn[i]){
				sn_indices.push_back(i);
				BN_indices.push_back(j);
			}
		}
	}	
	for(unsigned int i=0; i<BN_indices.size(); i++){
		for(unsigned int j=0; j<BN_indices.size(); j++){
			update.addset(BN_indices[i], BN_indices[j],
				node->schur_update.get(sn_indices[i], sn_indices[j]));
		}
	}
}


// Sets vec(b) = vec(b) + mat*vec(a)
void IeSolverTools::ApplySweepMatrix(const ie_Mat& mat, ie_Mat& vec, 
		const std::vector<unsigned int>& a, const std::vector<unsigned int>& b, 
		bool transpose = false) {
	if(a.size()*b.size()==0) return;

	//This vector is just used for indexing an Vector	
	std::vector<unsigned int> ZERO_VECTOR;
	ZERO_VECTOR.push_back(0);
	ie_Mat temp = vec(a, ZERO_VECTOR);
	ie_Mat product(b.size(), 1);

	if(transpose){
		Matmul::ie_gemv(TRANSPOSE, 1., mat, temp, 0., product);
	}else{
		Matmul::ie_gemv(NORMAL, 1., mat, temp, 0., product);
	}
	product += vec(b, ZERO_VECTOR);
	vec.set_submatrix( b, ZERO_VECTOR, product);
}


// Sets vec(range) = mat * vec(range)
void IeSolverTools::ApplyDiagMatrix(const ie_Mat& mat, ie_Mat& vec, 
	const std::vector<unsigned int>& range){
	if(range.size()==0) return;

	std::vector<unsigned int> ZERO_VECTOR;
	ZERO_VECTOR.push_back(0);
	ie_Mat temp = vec(range, ZERO_VECTOR);
	ie_Mat product(range.size(), 1);
	Matmul::ie_gemv(NORMAL, 1., mat, temp, 0., product);
		
	vec.set_submatrix( range, ZERO_VECTOR, product);
}


void IeSolverTools::ApplyDiagInvMatrix(const ie_Mat& mat, ie_Mat& vec, 
	const std::vector<unsigned int>& range){
	if(range.size()==0) return;

	std::vector<unsigned int> ZERO_VECTOR;
	ZERO_VECTOR.push_back(0);
	ie_Mat temp = vec(range, ZERO_VECTOR);
	ie_Mat product(range.size(), 1);
	mat.left_multiply_inverse(temp, product);

	vec.set_submatrix( range, ZERO_VECTOR, product);
}


void IeSolverTools::Skeletonize(const ie_Mat& K, QuadTree& tree){
	
	unsigned int lvls = tree.levels.size();
	for(unsigned int level = lvls-1; level>0; level--){	
		QuadTreeLevel* current_level = tree.levels[level];
		for(unsigned int n = 0; n<current_level->nodes.size(); n++){
			QuadTreeNode* current_node = current_level->nodes[n];
			
			//First, get rid of inactive dofs
			remove_inactive_dofs(current_node);

			if(current_node->interaction_lists.box.size() 
				< MIN_DOFS_TO_COMPRESS){
				continue;
			}

			int redundants = InterpolativeDecomposition(K, current_node->T, 
				current_node);
			if(redundants == 0) continue;

			SchurUpdate(K, current_node->T, current_node->L, current_node->U, 
				current_node);
		}
	}	
	remove_inactive_dofs(tree.root);
}


void IeSolverTools::SparseMatVec(const ie_Mat& K, QuadTree& tree, const ie_Mat& x, 
	ie_Mat& b){
	
	b = ie_Mat(x.height(), 1);
	x.copy(b);
	int lvls = tree.levels.size();
	for(int level = lvls-1; level>=0; level--){
		QuadTreeLevel* current_level = tree.levels[level];
		for (QuadTreeNode* current_node : current_level->nodes) {
		
			// First we need to apply L_T inverse
			// L_T inverse changes the skel elements - it makes them equal to
			// T times the redundant elements + the skeleton elements. 
			ApplySweepMatrix(current_node->T, b, 
				current_node->interaction_lists.redundant, current_node->interaction_lists.skel,
				 false);
			// Next we need to apply U inverse
			// U inverse changes the redundant elements - it makes them equal to 
			// L transpose times the skelnear elements + the redundant elements
			ApplySweepMatrix(current_node->U, b, 
				current_node->interaction_lists.skelnear, 
				current_node->interaction_lists.redundant, false);	
		}
	}
	
	//This can go through the tree in any order, is parallelizable 
	for(int level = lvls-1; level>=0; level--){
		QuadTreeLevel* current_level = tree.levels[level];
		for (QuadTreeNode* current_node : current_level->nodes) {
		
			ApplyDiagMatrix(current_node->D_r, b,
			 current_node->interaction_lists.redundant);
		}
	}
	
	// We need all of the skeleton indices. This is just the negation of
	// [0,b.size()] and the redundant DoFs
	// with that in mind...
	std::vector<unsigned int> allskel = tree.root->interaction_lists.box;
	
	ie_Mat allskel_mat(allskel.size(), allskel.size());
	get_all_schur_updates(allskel_mat, allskel, tree.root);
	allskel_mat*=-1;
	allskel_mat += K(allskel, allskel);
	ApplyDiagMatrix(allskel_mat, b, allskel);
	for(int level = 0; level<lvls; level++){
		QuadTreeLevel* current_level = tree.levels[level];
		// TODO record in notes and explore the following observation:
		// changing the order of this for loop affects the accuracy of the 
		// sparsematvec on a random vector, BUT NOT THE SOLUTION ERROR
		for(int n = current_level->nodes.size() - 1; n >= 0; n--){
		// for(unsigned int n = 0; n < current_level->nodes.size(); n++){
		
			QuadTreeNode* current_node = current_level->nodes[n];
		
			// Next we need to apply L inverse
			// L inverse changes the skelnear elements - it makes them equal to 
			// L times the redundant elements + the skelnear elements
			ApplySweepMatrix(current_node->L, b, 
				current_node->interaction_lists.redundant, 
				current_node->interaction_lists.skelnear, false);
			// Finally we need to apply U_T inverse
			// U_T inverse changes the redundant elements - it makes them equal
			// to T transpose times the skeleton elements + the redundant
			// elements
	 		ApplySweepMatrix(current_node->T, b, current_node->interaction_lists.skel,
	 		 current_node->interaction_lists.redundant, true);
	 	}
	}
}


void IeSolverTools::Solve(const ie_Mat& K, QuadTree& tree, ie_Mat& x, 
	const ie_Mat& b){

	int lvls = tree.levels.size();
	x = ie_Mat(b.height(), 1);
	b.copy(x);
	
	for (int level = lvls-1; level>=0; level--){//level>=0; level--){
		QuadTreeLevel* current_level = tree.levels[level];
		for (QuadTreeNode* current_node : current_level->nodes) {
			// Next we need to apply L inverse
			// L inverse changes the skelnear elements - it makes them equal to
			// L times the redundant elements + the skelnear elements
			current_node->L *= -1;
			current_node->T *= -1;
			// Finally we need to apply U_T inverse
			// U_T inverse changes the redundant elements - it makes them equal
			// to T transpose times the skeleton elements + the redundant
			// elements
	 		ApplySweepMatrix(current_node->T, x, current_node->interaction_lists.skel,
	 		 current_node->interaction_lists.redundant, true);
	 		ApplySweepMatrix(current_node->L, x, 
	 			current_node->interaction_lists.redundant, 
	 			current_node->interaction_lists.skelnear, false);
	 		
	 		current_node->L *= -1;
			current_node->T *= -1;
		}
	}
	//This can go through the tree in any order, is parallelizable 
	
	for(int level = lvls-1; level>=0; level--){//level>=0; level--){
		QuadTreeLevel* current_level = tree.levels[level];
		for (QuadTreeNode* current_node : current_level->nodes) {
			if(current_node->interaction_lists.redundant.size()==0) continue;
			ApplyDiagInvMatrix(current_node->D_r, x, 
				current_node->interaction_lists.redundant);
		}
	}

	// We need all of the skeleton indices. This is just the negation of 
	// [0,b.size()] and the redundant DoFs
	// with that in mind...
	std::vector<unsigned int> allskel = tree.root->interaction_lists.box;
	ie_Mat allskel_mat(allskel.size(), allskel.size());
	get_all_schur_updates(allskel_mat, allskel, tree.root);

	allskel_mat *= -1;
	allskel_mat += K(allskel, allskel);

	ApplyDiagInvMatrix(allskel_mat, x, allskel);
	for(int level = 0; level < lvls; level++){
		QuadTreeLevel* current_level = tree.levels[level];
		for(int n = current_level->nodes.size() - 1; n >= 0; n--){
		// for(unsigned int n = 0; n < current_level->nodes.size(); n++){
		
			QuadTreeNode* current_node = current_level->nodes[n];
		
			current_node->U *= -1;
			current_node->T *= -1;
			//First we need to apply L_T inverse
			// L_T inverse changes the skel elements - it makes them equal to T
			// times the redundant elements + the skeleton elements. 
			// Next we need to apply U inverse
			// U inverse changes the redundant elements - it makes them equal to 
			// L transpose times the skelnear elements + the redundant elements
			ApplySweepMatrix(current_node->U, x, 
				current_node->interaction_lists.skelnear, 
				current_node->interaction_lists.redundant, false);	
			ApplySweepMatrix(current_node->T, x, 
				current_node->interaction_lists.redundant, current_node->interaction_lists.skel,
				 false);
		
			current_node->U *= -1;
			current_node->T *= -1;
	 	}
	}	
}


void IeSolverTools::remove_inactive_dofs(QuadTreeNode* node){

	// this function removes from the box any DoFs which have already been made 
	// redundant. It involves a bunch of annoying C++ functions and probably 
	// would look nicer in matlab.

	std::vector<unsigned int> active_box, active_near;

	// populate active_box
	if(!node->is_leaf){
		for(QuadTreeNode* child: node->children){
			if(child->schur_updated){
				for(unsigned int i : child->interaction_lists.skel){
					active_box.push_back(i);
				}
			}else{
				for(unsigned int i : child->interaction_lists.box){
					active_box.push_back(i);
				}
			}
		}
		node->interaction_lists.box = active_box;
	}

	// populate active_near
	// note: this matters because the near dofs are used in the schur updating

	for(QuadTreeNode* neighbor: node->neighbors){
		if(neighbor->schur_updated){
			for(unsigned int i : neighbor->interaction_lists.skel){
				active_near.push_back(i);
			}
		}
		else if(!neighbor->is_leaf){
			for(QuadTreeNode* child: neighbor->children){
				if(child->schur_updated){
					for(unsigned int i : child->interaction_lists.skel){
						active_near.push_back(i);
					}
				}else{
					for(unsigned int i : child->interaction_lists.box){
						active_near.push_back(i);
					}
				}
			}
		}else{
			for(unsigned int i : neighbor->interaction_lists.box){
				active_near.push_back(i);
			}
		}
	}
	node->interaction_lists.near = active_near;	
}


void IeSolverTools::make_proxy_mat(ie_Mat& pxy, double cntr_x, double cntr_y, 
	double r, const std::vector<unsigned int>& range){
	
	//each row is a pxy point, cols are box dofs

	double scale = 1.0 / (2*M_PI);

	for(int i=0; i<100; i++){
		
		double ang = 2*M_PI*i*0.01;
		Vec2 p(cntr_x+r*cos(ang), cntr_y+r*sin(ang));
		for(unsigned int j_=0; j_<range.size(); j_++){
			unsigned int j = 2*range[j_];
			Vec2 q(points[j], points[j+1]);

			Vec2 r = p-q;
			Vec2 n(normals[j], normals[j+1]);

			double potential = -weights[j/2]*scale*(r.dot(n))/(r.dot(r));
			pxy.set(i, j_, potential);
		}
	}
}


void IeSolverTools::make_stokes_proxy_mat(ie_Mat& pxy, double cntr_x, double cntr_y, 
	double r, const std::vector<unsigned int>& range){

	double scale = 1.0 / (M_PI);

	for(int i=0; i<200; i+=2){
		
		double ang = 2*M_PI*i*0.01;
		Vec2 p(cntr_x + r*cos(ang), cntr_y + r*sin(ang));

		for(unsigned int j_=0; j_<range.size(); j_++){
			unsigned int j = range[j_];
			Vec2 q(points[j], points[j+1]);

			Vec2 r = p-q;
			Vec2 n(normals[j], normals[j+1]);
			double r0 = r.a[0];
			double r1 = r.a[1];
			
			double potential = weights[j/2]*scale*(r.dot(n))/(pow(r.dot(r),2));
		
			if( j % 2 == 0 ){
				pxy.set(i  , j_, potential*r0*r0);
				pxy.set(i+1, j_, potential*r0*r1);


			}else{
				pxy.set(i  , j_, potential*r0*r1);
				pxy.set(i+1, j_, potential*r1*r1);
			}
		}
	}
}

void IeSolverTools::set_rs_ranges(InteractionLists& interaction_lists, const std::vector<unsigned int>& prm, 
	unsigned int sk, unsigned int rd) {

	assert(prm.size() == sk+rd);

	for(unsigned int i = 0; i<sk; i++){
    	interaction_lists.skel.push_back(interaction_lists.box[prm[i]]);
    	interaction_lists.permutation.push_back(prm[i]);	
    }
    for(unsigned int i = sk; i<sk+rd; i++){
        interaction_lists.redundant.push_back(interaction_lists.box[prm[i]]);
        interaction_lists.permutation.push_back(prm[i]);
    }
}

void IeSolverTools::set_skelnear_range(InteractionLists& interaction_lists){
	for(unsigned int i = 0; i<interaction_lists.skel.size(); i++){
        interaction_lists.skelnear.push_back(interaction_lists.skel[i]);
    }
    for(unsigned int i=0; i<interaction_lists.near.size(); i++){
        interaction_lists.skelnear.push_back(interaction_lists.near[i]);
    }
}



// void IeSolverTools::Conjugate_Gradient(ie_Mat& K, QuadTree& tree, ie_Mat& phi, 
// 	ie_Mat& f){
// 	int iteration = 0;
// 	Zero(phi);
// 	ie_Mat r(f.height(), f.width());
// 	Copy(f, r);
// 	ie_Mat p(f.height(), f.width());
// 	Copy(f,p);

// 	ie_Mat Ap(p.height(), 1);


// 	double r_dot = Dot(r, r);
		
// 	for(; iteration<100; iteration++){
// 		//r^T r
		
// 		SparseMatVec(F, tree, p, Ap);
// 		double p_dot = Dot(p, Ap);

// 		double alpha = r_dot / p_dot;

// 		//fix this stupidity
// 		p *= alpha;
// 		phi += p;
// 		p*= 1.0/alpha;

// 		Ap *= alpha;
// 		r -= Ap;
// 		Ap *= 1.0/alpha;



// 		double new_rdot = Dot(r, r);

// 		if(new_rdot < 0.0000001){
// 			printf("Done in %d iterations\n", iteration);
// 			break;
// 		}
// 		double beta = new_rdot / r_dot;

// 		r_dot = new_rdot;

// 		p *= beta;
// 		p += r;


// 	}
	
// }

} // namespace ie_solver