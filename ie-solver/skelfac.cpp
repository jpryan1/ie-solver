#include "skelfac.h"

namespace ie_solver{

void Skelfac::GetXMatrices(ie_Mat& K, ie_Mat& Z, ie_Mat& Xrr, 
	std::vector<unsigned int>& r, std::vector<unsigned int>& s, 
	std::vector<unsigned int>& n) {
    
	unsigned int rs = r.size();
	unsigned int ss = s.size();
	unsigned int ns = n.size();
    ie_Mat  Xrs(rs, ss), Xsr(ss, rs),
	Xrn(rs, ns), Xnr(ns, rs);
	
	//this is just for readability
	Xrr = K(r,r);
	Xrs = K(r,s);
	Xrn = K(r,n);
	Xsr = K(s,r);
	Xnr = K(n,r);

	Matmul::ie_gemm(TRANSPOSE_, NORMAL_, -1., Z, 	  K(s,r), 1., Xrr);
	Matmul::ie_gemm(TRANSPOSE_, NORMAL_, -1., Z, 	  K(s,s), 1., Xrs);
	Matmul::ie_gemm(TRANSPOSE_, NORMAL_, -1., Z, 	  K(s,n), 1., Xrn);
	Matmul::ie_gemm(NORMAL_, 	NORMAL_, -1., Xrs, 	  Z, 	  1., Xrr);
	Matmul::ie_gemm(NORMAL_, 	NORMAL_, -1., K(s,s), Z, 	  1., Xsr);
	Matmul::ie_gemm(NORMAL_, 	NORMAL_, -1., K(n,s), Z, 	  1., Xnr);

	set_get.tic();
	K.set_submatrix(r, s, Xrs);
	K.set_submatrix(r, n, Xrn);
	//K.set_submatrix(r, r, Xrr);
	K.set_submatrix(s, r, Xsr);
	K.set_submatrix(n, r, Xnr);
	set_get.toc();
	// TODO this seems like an awful lot of stores, can we avoid this? 
}


//TODO These involved multiplying into a buffer, then copying the buffer into K
	//Can we do this in place? Or use just one buffer and resize each time?
	//Does it matter?
	//Also, we might not need to do ALL of these matmuls, since some of the 
	//blocks will be eliminated anyways. This should be worked out on a whiteboard
void Skelfac::SchurUpdate(ie_Mat& K, ie_Mat& Z, ie_Mat& L, ie_Mat& U, QuadTreeNode* node) {
    //height of Z is number of skeleton columns
	unsigned int num_box 	  = node->box.box_range.size();
	unsigned int num_redundant = Z.width();
	unsigned int num_skel 	  = Z.height();
	unsigned int num_near 	  = node->box.near_range.size();

	//GENERATE K_BN,BN
	std::vector<unsigned int> BN;
	for(unsigned int i=0; i<num_box ; i++) BN.push_back(node->box.box_range[i]);
	for(unsigned int i=0; i<num_near; i++) BN.push_back(node->box.near_range[i]);

	ie_Mat K_BN = K(BN, BN); //que bien!
	ie_Mat update(BN.size(), BN.size());
	get_all_schur_updates(update, BN, node);
	K_BN -= update;
	//printf("Update norm is %f\n", testing.norm2());

	// Generate various index ranges within BN
	std::vector<unsigned int> s, r, n, sn;
	for(unsigned int i=0; i<num_skel; i++){
		s.push_back(node->box.p[i]);
		sn.push_back(node->box.p[i]);
	}
	
	for(unsigned int i=0; i<num_redundant; i++){
		r.push_back(node->box.p[i+num_skel]);
	}

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
	
	Xrr.right_multiply_inverse(K_BN(sn, r), L);//F.K(node->box.skelnear_range, node->box.redundant_range), L);
	Xrr.left_multiply_inverse( K_BN(r, sn), U);//F.K(node->box.redundant_range, node->box.skelnear_range), U);

	ie_Mat schur(node->box.skelnear_range.size(), node->box.skelnear_range.size());
	Matmul::ie_gemm(NORMAL_, NORMAL_, 1.0, L, 
		K_BN(r,sn), 0., schur);
	
	//set schur update
	node->schur_update = schur;
	node->schur_updated = true;
}


int Skelfac::InterpolativeDecomposition( ie_Mat& K, ie_Mat& Z, 
	QuadTreeNode* node) {

	double cntr_x = node->corners[0] + node->side/2.0;
	double cntr_y = node->corners[1] + node->side/2.0;
	ie_Mat pxy;

	if(is_stokes){
		pxy = ie_Mat(200, node->box.box_range.size());
		make_stokes_proxy_mat(pxy, cntr_x, cntr_y, node->side*2, node->box.box_range );
	}
	else{
		pxy = ie_Mat(100, node->box.box_range.size());
		make_proxy_mat(pxy, cntr_x, cntr_y, node->side*2, node->box.box_range );
	}
	std::vector<unsigned int> p;
	unsigned int numskel = pxy.id(p, Z, id_tol);
	if(numskel==0) return 0;
	set_rs_ranges(node->box, p, Z.height(), Z.width());
	
	set_skelnear_range(node->box);
	
	return Z.width();
}


void Skelfac::Skeletonize(ie_Mat& K, QuadTree& tree){
	
	unsigned int lvls = tree.levels.size();
	for(unsigned int level = lvls-1; level>0; level--){	
		QuadTreeLevel* current_level = tree.levels[level];
		for(unsigned int n = 0; n<current_level->nodes.size(); n++){
			QuadTreeNode* current_node = current_level->nodes[n];

			if(current_node->box.box_range.size()==0){
				continue;
			}
			
			ie_Mat Z,L, U;
			//First, get rid of inactive dofs
			remove_inactive_dofs(current_node);

			
			if(current_node->box.box_range.size()<=1){
				continue;
			}

			int redundants = InterpolativeDecomposition(K, Z, current_node);
			if(redundants == 0) continue;

			current_node->T = Z;
			SchurUpdate(K, Z, L, U, current_node);

			set_get.tic();
			current_node->L = L;
			current_node->U = U;
			set_get.toc();
		}
	}	
	remove_inactive_dofs(tree.root);
}


//This function sets vec(b) = vec(b) + mat*vec(a)
void Skelfac::ApplySweepMatrix(ie_Mat& mat, ie_Mat& vec, 
	std::vector<unsigned int> a, std::vector<unsigned int> b, 
	bool transpose = false) {
	if(a.size()*b.size()==0) return;

	//This vector is just used for indexing an Vector	
	std::vector<unsigned int> ZERO_VECTOR;
	ZERO_VECTOR.push_back(0);

	ie_Mat temp = vec(a, ZERO_VECTOR);
	ie_Mat product(b.size(), 1);
	if(transpose){
		Matmul::ie_gemv(TRANSPOSE_, 1., mat, temp, 0., product);
	}else{
		Matmul::ie_gemv(NORMAL_, 1., mat, temp, 0., product);
	}
	product += vec(b, ZERO_VECTOR);
	vec.set_submatrix( b, ZERO_VECTOR, product);
}


//This function sets vec(range) = mat * vec(range)
void Skelfac::ApplyDiagMatrix(ie_Mat& K, ie_Mat& vec, std::vector<unsigned int> range){
	if(range.size()==0) return;

	std::vector<unsigned int> ZERO_VECTOR;
	ZERO_VECTOR.push_back(0);
	ie_Mat temp = vec(range, ZERO_VECTOR);
	ie_Mat product(range.size(), 1);
	Matmul::ie_gemv(NORMAL_, 1., K, temp, 0., product);
		
	vec.set_submatrix( range, ZERO_VECTOR, product);
}


void Skelfac::ApplyDiagInvMatrix(ie_Mat& K, ie_Mat& vec, std::vector<unsigned int> range){
	if(range.size()==0) return;

	std::vector<unsigned int> ZERO_VECTOR;
	ZERO_VECTOR.push_back(0);
	ie_Mat temp = vec(range, ZERO_VECTOR);
	ie_Mat product(range.size(), 1);
	K.left_multiply_inverse(temp, product);
	//Matmul::ie_gemv(NORMAL_, 1., K, temp, 0., product);
		
	vec.set_submatrix( range, ZERO_VECTOR, product);
}


void Skelfac::SparseMatVec(ie_Mat& K, QuadTree& tree, ie_Mat& x, ie_Mat& b){

	LOG::INFO("Beginning sparse matrix vector multiply...");
	
	b = ie_Mat(x.height(), 1);
	x.copy( b);
	int lvls = tree.levels.size();
	LOG::INFO("Begin sweep up...");
	for(int level = lvls-1; level>=0; level--){
		LOG::INFO("Level " + std::to_string(level));
		QuadTreeLevel* current_level = tree.levels[level];
		for (QuadTreeNode* current_node : current_level->nodes) {
		
			// First we need to apply L_T inverse
			// L_T inverse changes the skel elements - it makes them equal to T times the redundant elements + the skeleton elements. 
			ApplySweepMatrix(current_node->T, b, current_node->box.redundant_range, current_node->box.skel_range, false);
			// Next we need to apply U inverse
			// U inverse changes the redundant elements - it makes them equal to L transpose times the skelnear elements + the redundant elements
			ApplySweepMatrix(current_node->U, b, current_node->box.skelnear_range, current_node->box.redundant_range, false);	
		}
	}
	LOG::INFO("End sweep up.");
	
	//This can go through the tree in any order, is parallelizable 
	LOG::INFO("Begin block diagonal multiply...");
	for(int level = lvls-1; level>=0; level--){
		LOG::INFO("Level " + std::to_string(level));
		QuadTreeLevel* current_level = tree.levels[level];
		for (QuadTreeNode* current_node : current_level->nodes) {
		
			ApplyDiagMatrix(current_node->D_r, b, current_node->box.redundant_range);
		}
	}
	LOG::INFO("End block diagonal multiply.");
	
	// This is just for initial testing - once the tree structure is working properly, this will be
	// implemented in a smarter fashion

	// We need all of the skeleton indices. This is just the negation of [0,b.size()] and the redundant DoFs
	// with that in mind...
	std::vector<unsigned int> allskel = tree.root->box.box_range;
	
	ie_Mat allskel_mat(allskel.size(), allskel.size());
	get_all_schur_updates(allskel_mat, allskel, tree.root);
	allskel_mat*=-1;
	allskel_mat += K(allskel, allskel);
	ApplyDiagMatrix(allskel_mat, b, allskel);
	LOG::INFO("Begin sweep down...");
	for(int level = 0; level<lvls; level++){
		LOG::INFO("Level " + std::to_string(level));
		QuadTreeLevel* current_level = tree.levels[level];
		for (QuadTreeNode* current_node : current_level->nodes) {
		
			// Next we need to apply L inverse
			// L inverse changes the skelnear elements - it makes them equal to L times the redundant elements + the skelnear elements
			ApplySweepMatrix(current_node->L, b, current_node->box.redundant_range, current_node->box.skelnear_range,false);
			// Finally we need to apply U_T inverse
			// U_T inverse changes the redundant elements - it makes them equal to T transpose times the skeleton elements + the redundant elements
	 		ApplySweepMatrix(current_node->T, b, current_node->box.skel_range, current_node->box.redundant_range, true);
	 	}
	}
	LOG::INFO("End sweep down.");
	LOG::INFO("End sparse matrix vector multiply.");
	
}


void Skelfac::Solve( ie_Mat& K, QuadTree& tree, ie_Mat& x, ie_Mat& b){

	int lvls = tree.levels.size();

	x = ie_Mat(b.height(), 1);

	b.copy(x);

	for (int level = lvls-1; level>=0; level--){//level>=0; level--){
		QuadTreeLevel* current_level = tree.levels[level];
		for (QuadTreeNode* current_node : current_level->nodes) {
		//Next we need to apply L inverse
		//L inverse changes the skelnear elements - it makes them equal to L times the redundant elements + the skelnear elements
			current_node->L*=-1;
			current_node->T*=-1;
			
		//Finally we need to apply U_T inverse
		//U_T inverse changes the redundant elements - it makes them equal to T transpose times the skeleton elements + the redundant elements
	 		ApplySweepMatrix(current_node->T, x, current_node->box.skel_range, current_node->box.redundant_range, true);
	 		ApplySweepMatrix(current_node->L, x, current_node->box.redundant_range, current_node->box.skelnear_range,false);
	 		current_node->L*=-1;
			current_node->T*=-1;
		}
	}
	//This can go through the tree in any order, is parallelizable 
	
	for(int level = lvls-1; level>=0; level--){//level>=0; level--){
		QuadTreeLevel* current_level = tree.levels[level];
		for (QuadTreeNode* current_node : current_level->nodes) {
			if(current_node->box.redundant_range.size()==0) continue;

			ie_Mat red_mat = current_node->D_r;
			//red_mat.inverse();
	
			ApplyDiagInvMatrix(red_mat, x, current_node->box.redundant_range);
			

		}
	}

	// This is just for initial testing - once the tree structure is working properly, this will be
	// implemented in a smarter fashion

	// We need all of the skeleton indices. This is just the negation of [0,b.size()] and the redundant DoFs
	// with that in mind...
	std::vector<unsigned int> allskel = tree.root->box.box_range;
	
	ie_Mat allskel_mat(allskel.size(), allskel.size());
	get_all_schur_updates(allskel_mat, allskel, tree.root);

	allskel_mat *=-1;
	allskel_mat += K(allskel, allskel);

	// if(is_stokes){
	// 	int buddy_counter = 0;

	// 	for(int i=0; i<allskel.size()-1; i++){
	// 		if(i%2==1) continue;
	// 		if(allskel[i+1] == allskel[i] + 1) buddy_counter++;
	// 	}
	// 	//printf("Of the %lu skeleton columns that remain, there are %d DoF pairs, hence %.0f%% are solo\n", 
	// 	//		allskel.size(), buddy_counter, 100*(allskel.size()-2*buddy_counter+0.0)/allskel.size());
	// }
	//allskel_mat.inverse();
	ApplyDiagInvMatrix(allskel_mat, x, allskel);
	
	for(int level = 0; level < lvls; level++){
		QuadTreeLevel* current_level = tree.levels[level];
		for (QuadTreeNode* current_node : current_level->nodes) {
		
			current_node->U*=-1;
			current_node->T*=-1;
	//First we need to apply L_T inverse
		// //L_T inverse changes the skel elements - it makes them equal to T times the redundant elements + the skeleton elements. 
			// //Next we need to apply U inverse
		// //U inverse changes the redundant elements - it makes them equal to L transpose times the skelnear elements + the redundant elements
			ApplySweepMatrix(current_node->U, x, current_node->box.skelnear_range, current_node->box.redundant_range, false);	
			ApplySweepMatrix(current_node->T, x, current_node->box.redundant_range, current_node->box.skel_range, false);
		
			current_node->U*=-1;
			current_node->T*=-1;
	 	}
	}		
}


void Skelfac::remove_inactive_dofs(QuadTreeNode* node){

	//this function removes from the box any DoFs which have already been made redundant. 
	//It involves a bunch of annoying C++ functions and probably would look nicer in matlab

	if(!node->is_leaf){
		std::vector<unsigned int> active_box;
		for(QuadTreeNode* child: node->children){
			if(child->schur_updated){
				for(unsigned int i : child->box.skel_range) active_box.push_back(i);
			}else{
				for(unsigned int i : child->box.box_range) active_box.push_back(i);
			}
		}
		node->box.box_range = active_box;
	}

	// std::sort(box.box_range.begin(), box.box_range.end());
	// std::set_difference(node->box.box_range.begin(), node->box.box_range.end(),
	// 	F.redundant_dofs.begin(), F.redundant_dofs.end(), std::back_inserter(active_box));

	std::vector<unsigned int> active_near;	
	for(QuadTreeNode* neighbor: node->neighbors){
		if(neighbor->schur_updated){
			for(unsigned int i : neighbor->box.skel_range) active_near.push_back(i);
		}
		else if(!neighbor->is_leaf){
			for(QuadTreeNode* child: neighbor->children){
				if(child->schur_updated){
					for(unsigned int i : child->box.skel_range) active_near.push_back(i);
				}else{
					for(unsigned int i : child->box.box_range) active_near.push_back(i);
				}
			}
		}else{
			for(unsigned int i : neighbor->box.box_range) active_near.push_back(i);
		}
	}
	node->box.near_range = active_near;	
}


void Skelfac::make_proxy_mat(ie_Mat& pxy, double cntr_x, double cntr_y, double r, std::vector<unsigned int>& range){
	
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


void Skelfac::make_stokes_proxy_mat(ie_Mat& pxy, double cntr_x, double cntr_y, 
	double r, std::vector<unsigned int>& range){

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

void Skelfac::set_rs_ranges(Box& box, std::vector<unsigned int>& prm, unsigned int sk, unsigned int rd) {

	assert(prm.size() == sk+rd);

	for(unsigned int i = 0; i<sk; i++){
       box.skel_range.push_back(box.box_range[prm[i]]);
    	box.p.push_back(prm[i]);	
    }
    for(unsigned int i = sk; i<sk+rd; i++){
        box.redundant_range.push_back(box.box_range[prm[i]]);
        box.p.push_back(prm[i]);
    }
}

void Skelfac::set_skelnear_range(Box& box){
	for(unsigned int i = 0; i<box.skel_range.size(); i++){
        box.skelnear_range.push_back(box.skel_range[i]);
    }
    for(unsigned int i=0; i<box.near_range.size(); i++){
        box.skelnear_range.push_back(box.near_range[i]);
    }
}


void Skelfac::get_all_schur_updates(ie_Mat& updates, std::vector<unsigned int>& BN, 
	QuadTreeNode* node){
	//	printf("Getting all neighbors and descendents on %d...\n",node->id);
	if(!node->is_leaf) get_descendents_updates(updates, BN, node);
	//Check each neighbor's neighbor, if updated, grab its schur update, then recurse on everyone's children
	for(QuadTreeNode* neighbor : node->neighbors){
		if(neighbor->schur_updated) get_update(updates, BN, neighbor);
		if(!neighbor->is_leaf) get_descendents_updates(updates, BN, neighbor);
	}

	for(QuadTreeNode* neighbor : node->far_neighbors){
		if(neighbor->schur_updated) get_update(updates, BN, neighbor);
		if(!neighbor->is_leaf) get_descendents_updates(updates, BN, neighbor);
	}
	//printf("Done with neighbors descendents on %d!\n", node->id);
}
	

void Skelfac::get_descendents_updates(ie_Mat& updates, std::vector<unsigned int>& BN, QuadTreeNode* node){
	//by assumption, node is not a leaf
	
	for(QuadTreeNode* child : node->children){
		if(child->schur_updated) get_update(updates, BN, child);
		if(!child->is_leaf) get_descendents_updates(updates, BN, child);
	}
}


void Skelfac::get_update(ie_Mat& update, std::vector<unsigned int>& BN, QuadTreeNode* node){
	
	//node needs to check all its dofs against BN, enter interactions into corresponding locations
	
	//node only updated its own BN dofs, and the redundant ones are no longer relevant, so we only care 
	//about child's SN dofs

	std::vector<unsigned int> sn = node->box.skelnear_range;

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
			update.addset(BN_indices[i], BN_indices[j],node->schur_update.get(sn_indices[i], sn_indices[j]));
		}
	}
}


void Skelfac::Conjugate_Gradient(ie_Mat& K, QuadTree& tree, ie_Mat& phi, ie_Mat& f){
	// int iteration = 0;
	// Zero(phi);
	// ie_Mat r(f.height(), f.width());
	// Copy(f, r);
	// ie_Mat p(f.height(), f.width());
	// Copy(f,p);

	// ie_Mat Ap(p.height(), 1);


	// double r_dot = Dot(r, r);
		
	// for(; iteration<100; iteration++){
	// 	//r^T r
		
	// 	SparseMatVec(F, tree, p, Ap);
	// 	double p_dot = Dot(p, Ap);

	// 	double alpha = r_dot / p_dot;

	// 	//fix this stupidity
	// 	p *= alpha;
	// 	phi += p;
	// 	p*= 1.0/alpha;

	// 	Ap *= alpha;
	// 	r -= Ap;
	// 	Ap *= 1.0/alpha;



	// 	double new_rdot = Dot(r, r);

	// 	if(new_rdot < 0.0000001){
	// 		printf("Done in %d iterations\n", iteration);
	// 		break;
	// 	}
	// 	double beta = new_rdot / r_dot;

	// 	r_dot = new_rdot;

	// 	p *= beta;
	// 	p += r;


	// }
	
}

} // namespace ie_solver