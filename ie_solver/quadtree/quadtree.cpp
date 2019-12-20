// Copyright 2019 John Paul Ryan
#include <omp.h>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <iostream>
#include <utility>
#include <boost/functional/hash.hpp>
#include "ie_solver/log.h"
#include "ie_solver/quadtree/quadtree.h"
#include "ie_solver/boundaries/circle.h"
#include "ie_solver/boundaries/annulus.h"
#include "ie_solver/boundaries/rounded_square.h"
#include "ie_solver/boundaries/cubic_spline.h"
#include "ie_solver/boundaries/donut.h"
#include "ie_solver/boundaries/ex1boundary.h"
#include "ie_solver/boundaries/ex2boundary.h"
#include "ie_solver/boundaries/ex3boundary.h"

namespace ie_solver {

typedef std::pair<double, double> pair;
int QuadTreeNode::id_count = 0;

QuadTree::~QuadTree() {
  for (QuadTreeLevel* level : levels) {
    for (QuadTreeNode* node : level->nodes) {
      delete node;
    }
  }
  for (QuadTreeLevel* level : levels) {
    if (level) {
      delete level;
    }
  }
  levels.clear();
}


void QuadTree::initialize_tree(Boundary* boundary_,
                               const std::vector<double>& domain_points_,
                               int solution_dimension_,
                               int domain_dimension_) {
  // double tree_start_time = omp_get_wtime();
  assert(boundary_->points.size() > 0
         && "number of boundary->points to init tree cannot be 0.");
  // todo(john) later we can assert that the boundary->points are a multiple of
  // dimension
  this->boundary = boundary_;
  this->domain_points = domain_points_;
  this->solution_dimension = solution_dimension_;
  this->domain_dimension = domain_dimension_;
  min = boundary->points[0];
  max = boundary->points[0];

  for (double point : boundary->points) {
    if (point < min) min = point;
    if (point > max) max = point;
  }
  double tree_min = min - 0.01;
  double tree_max = max + 0.01;

  root = new QuadTreeNode();
  root->level = 0;
  root->parent = nullptr;
  // bl
  root->corners[0] = tree_min;
  root->corners[1] = tree_min;
  // tl
  root->corners[2] = tree_min;
  root->corners[3] = tree_max;
  // tr
  root->corners[4] = tree_max;
  root->corners[5] = tree_max;
  // br
  root->corners[6] = tree_max;
  root->corners[7] = tree_min;
  root->side_length = tree_max - tree_min;
  QuadTreeLevel* level_one = new QuadTreeLevel();
  level_one->nodes.push_back(root);
  levels.push_back(level_one);

  for (int i = 0; i < boundary->points.size(); i += domain_dimension) {
    recursive_add(this->root, boundary->points[i], boundary->points[i + 1],
                  i / domain_dimension, true);
  }
  for (int i = 0; i < domain_points.size(); i += domain_dimension) {
    recursive_add(this->root, domain_points[i], domain_points[i + 1],
                  i / domain_dimension, false);
  }

  // make neighbor lists in a stupid way
  for (int level = 0; level < levels.size(); level++) {
    QuadTreeLevel* current_level = levels[level];
    for (int k = 0; k < current_level->nodes.size(); k++) {
      QuadTreeNode* node_a = current_level->nodes[k];

      // Each node is neighbors with all its siblings
      if (node_a->parent != nullptr) {
        for (QuadTreeNode* sibling : node_a->parent->children) {
          if (sibling->id != node_a->id) node_a->neighbors.push_back(sibling);
        }

        // Now check all parents' neighbors' children
        for (QuadTreeNode* parents_neighbor : node_a->parent->neighbors) {
          for (QuadTreeNode* cousin : parents_neighbor->children) {
            if (cousin == nullptr) continue;
            if (cousin->level != node_a->level) continue;
            double dist = sqrt(pow(node_a->corners[0] - cousin->corners[0], 2)
                               + pow(node_a->corners[1] - cousin->corners[1], 2));
            // just need to check if the distance of the bl corners
            // is <=s*sqrt(2)
            if (dist < node_a->side_length * sqrt(2) + 1e-5) {
              node_a->neighbors.push_back(cousin);
            }
          }
        }
      }

      // now if it is a leaf, check against nodes in all subsequent levels
      if (node_a->is_leaf) {
        for (int n = 0; n < node_a->neighbors.size(); n++) {
          QuadTreeNode* neighbor =  node_a->neighbors[n];
          // make sure this isn't a neighbor from a higher level
          if (neighbor->level != node_a->level) {
            continue;
          }
          for (QuadTreeNode* child : neighbor->children) {
            if (child != nullptr) {
              get_descendent_neighbors(node_a, child);
            }
          }
        }
      }
    }
  }

  // double tree_end_time = omp_get_wtime();
  // std::cout << "timing: tree_init " << (tree_end_time - tree_start_time) <<
  //           std::endl;
}


// adds neighbors to leaf which are on lower levels, by recursing and checking
// if the corners are along the wall.
void QuadTree::get_descendent_neighbors(QuadTreeNode* big,
                                        QuadTreeNode* small) {
  assert(big->level < small->level);
  double top = big->corners[3];
  double bottom = big->corners[1];
  double left = big->corners[0];
  double right = big->corners[4];
  for (int i = 0; i < 8; i += 2) {
    double x = small->corners[i];
    double y = small->corners[i + 1];
    // note: in theory, the equalities should be exact (i think), but to be
    // overcareful we will allow for machine error
    if (x > left && x < right) {
      if (fabs(y - top) < 1e-14) {
        big->neighbors.push_back(small);
        small->neighbors.push_back(big);
        break;
      }
      if (fabs(y - bottom) < 1e-14) {
        big->neighbors.push_back(small);
        small->neighbors.push_back(big);
        break;
      }
    }
    if (y < top && y < bottom) {
      if (fabs(x - left) < 1e-14) {
        big->neighbors.push_back(small);
        small->neighbors.push_back(big);
        break;
      }
      if (fabs(x - right) < 1e-14) {
        big->neighbors.push_back(small);
        small->neighbors.push_back(big);
        break;
      }
    }
  }
  for (QuadTreeNode* child : small->children) {
    if (child != nullptr) {
      get_descendent_neighbors(big, child);
    }
  }
}


void QuadTree::recursive_add(QuadTreeNode* node, double x, double y,
                             int point_ind, bool is_boundary) {
  assert(node != nullptr && "recursive_add fails on null node.");
  for (int i = 0; i < solution_dimension; i++) {
    if (is_boundary) {
      node->src_dof_lists.original_box.push_back(
        solution_dimension * point_ind + i);
    } else {
      node->tgt_dof_lists.original_box.push_back(
        solution_dimension * point_ind + i);
    }
  }
  // figure out which child
  double midx = ((node->corners[6] - node->corners[0]) / 2.0)
                + node->corners[0];
  double midy = ((node->corners[3] - node->corners[1]) / 2.0)
                + node->corners[1];
  QuadTreeNode* child;
  if (x < midx && y < midy) {
    child = node->bl;
  } else if (x < midx && y >= midy) {
    child = node->tl;
  } else if (x >= midx && y < midy) {
    child = node->br;
  } else {
    child = node->tr;
  }
  // does that child exist?
  if (child) {
    recursive_add(child, x, y, point_ind, is_boundary);
  } else {
    // do we need one?
    // if this node is exploding and needs children
    if (node->is_leaf
        && node->src_dof_lists.original_box.size() +
        node->tgt_dof_lists.original_box.size() > MAX_LEAF_DOFS) {
      node_subdivide(node);
    }
  }
}


// 1) gives node its four children
// 2) puts these children in their proper level
// 3) gives these children their corners
void QuadTree::node_subdivide(QuadTreeNode* node) {
  assert(node != nullptr && "node_subdivide fails on null node.");
  node->is_leaf = false;
  QuadTreeNode *bl, *br, *tl, *tr;
  double midx = ((node->corners[6] - node->corners[0]) / 2.0)
                + node->corners[0];
  double midy = ((node->corners[3] - node->corners[1]) / 2.0)
                + node->corners[1];

  bl = new QuadTreeNode();
  bl->corners[0] = node->corners[0];
  bl->corners[1] = node->corners[1];
  bl->corners[2] = node->corners[2];
  bl->corners[3] = midy;
  bl->corners[4] = midx;
  bl->corners[5] = midy;
  bl->corners[6] = midx;
  bl->corners[7] = node->corners[7];

  tl = new QuadTreeNode();
  tl->corners[0] = node->corners[0];
  tl->corners[1] = midy;
  tl->corners[2] = node->corners[2];
  tl->corners[3] = node->corners[3];
  tl->corners[4] = midx;
  tl->corners[5] = node->corners[5];
  tl->corners[6] = midx;
  tl->corners[7] = midy;

  tr = new QuadTreeNode();
  tr->corners[0] = midx;
  tr->corners[1] = midy;
  tr->corners[2] = midx;
  tr->corners[3] = node->corners[3];
  tr->corners[4] = node->corners[4];
  tr->corners[5] = node->corners[5];
  tr->corners[6] = node->corners[6];
  tr->corners[7] = midy;

  br = new QuadTreeNode();
  br->corners[0] = midx;
  br->corners[1] = node->corners[1];
  br->corners[2] = midx;
  br->corners[3] = midy;
  br->corners[4] = node->corners[4];
  br->corners[5] = midy;
  br->corners[6] = node->corners[6];
  br->corners[7] = node->corners[7];

  bl->level = node->level + 1;
  br->level = node->level + 1;
  tl->level = node->level + 1;
  tr->level = node->level + 1;

  bl->side_length = node->side_length / 2.0;
  br->side_length = node->side_length / 2.0;
  tl->side_length = node->side_length / 2.0;
  tr->side_length = node->side_length / 2.0;

  bl->parent = node;
  br->parent = node;
  tl->parent = node;
  tr->parent = node;
  if (levels.size() < node->level + 2) {
    QuadTreeLevel* new_level = new QuadTreeLevel();
    levels.push_back(new_level);
    new_level->nodes.push_back(bl);
    new_level->nodes.push_back(br);
    new_level->nodes.push_back(tr);
    new_level->nodes.push_back(tl);
  } else {
    levels[node->level + 1]->nodes.push_back(bl);
    levels[node->level + 1]->nodes.push_back(br);
    levels[node->level + 1]->nodes.push_back(tr);
    levels[node->level + 1]->nodes.push_back(tl);
  }
  node->bl = bl;
  node->tl = tl;
  node->tr = tr;
  node->br = br;

  node->children[0] = bl;
  node->children[1] = tl;
  node->children[2] = tr;
  node->children[3] = br;

  // Now we bring the indices from the parent's box down into its childrens
  // boxes
  for (int index = 0; index < node->src_dof_lists.original_box.size();
       index += solution_dimension) {
    int matrix_index = node->src_dof_lists.original_box[index];
    int points_vec_index = (matrix_index / solution_dimension) *
                           domain_dimension;
    // So we are trying to be general, but x,y implies 2D. To generalize later
    double x = boundary->points[points_vec_index];
    double y = boundary->points[points_vec_index + 1];
    if (x < midx && y < midy) {
      for (int i = 0; i < solution_dimension; i++) {
        bl->src_dof_lists.original_box.push_back(matrix_index + i);
      }
    } else if (x < midx && y >= midy) {
      for (int i = 0; i < solution_dimension; i++) {
        tl->src_dof_lists.original_box.push_back(matrix_index + i);
      }
    } else if (x >= midx && y < midy) {
      for (int i = 0; i < solution_dimension; i++) {
        br->src_dof_lists.original_box.push_back(matrix_index + i);
      }
    } else {
      for (int i = 0; i < solution_dimension; i++) {
        tr->src_dof_lists.original_box.push_back(matrix_index + i);
      }
    }
  }
  // Bring the tgt indices down too
  for (int index = 0; index < node->tgt_dof_lists.original_box.size();
       index += solution_dimension) {
    int matrix_index = node->tgt_dof_lists.original_box[index];
    int points_vec_index = (matrix_index / solution_dimension) *
                           domain_dimension;
    // So we are trying to be general, but x,y implies 2D. To generalize later
    double x = domain_points[points_vec_index];
    double y = domain_points[points_vec_index + 1];
    if (x < midx && y < midy) {
      for (int i = 0; i < solution_dimension; i++) {
        bl->tgt_dof_lists.original_box.push_back(matrix_index + i);
      }
    } else if (x < midx && y >= midy) {
      for (int i = 0; i < solution_dimension; i++) {
        tl->tgt_dof_lists.original_box.push_back(matrix_index + i);
      }
    } else if (x >= midx && y < midy) {
      for (int i = 0; i < solution_dimension; i++) {
        br->tgt_dof_lists.original_box.push_back(matrix_index + i);
      }
    } else {
      for (int i = 0; i < solution_dimension; i++) {
        tr->tgt_dof_lists.original_box.push_back(matrix_index + i);
      }
    }
  }
  if (bl->src_dof_lists.original_box.size() +
      bl->tgt_dof_lists.original_box.size() > MAX_LEAF_DOFS) {
    node_subdivide(bl);
  }
  if (tl->src_dof_lists.original_box.size() +
      tl->tgt_dof_lists.original_box.size() > MAX_LEAF_DOFS) {
    node_subdivide(tl);
  }
  if (tr->src_dof_lists.original_box.size() +
      tr->tgt_dof_lists.original_box.size() > MAX_LEAF_DOFS) {
    node_subdivide(tr);
  }
  if (br->src_dof_lists.original_box.size() +
      br->tgt_dof_lists.original_box.size() > MAX_LEAF_DOFS) {
    node_subdivide(br);
  }
}


void QuadTree::mark_neighbors_and_parents(QuadTreeNode * node) {
  if (node == nullptr) return;
  //TODO(John) make all this a fn of nodes.
  node->compressed = false;
  node->X_rr_is_LU_factored = false;
  node->src_dof_lists.active_box.clear();
  node->src_dof_lists.skel.clear();
  node->src_dof_lists.skelnear.clear();
  node->src_dof_lists.redundant.clear();
  node->src_dof_lists.permutation.clear();
  node->T = ie_Mat(0, 0);
  node->U = ie_Mat(0, 0);
  node->L = ie_Mat(0, 0);
  for (QuadTreeNode* neighbor : node->neighbors) {
    neighbor->compressed = false;
    neighbor->X_rr_is_LU_factored = false;
    neighbor->src_dof_lists.active_box.clear();
    neighbor->src_dof_lists.skel.clear();
    neighbor->src_dof_lists.skelnear.clear();
    neighbor->src_dof_lists.redundant.clear();
    neighbor->src_dof_lists.permutation.clear();
    neighbor->T = ie_Mat(0, 0);
    neighbor->U = ie_Mat(0, 0);
    neighbor->L = ie_Mat(0, 0);
  }
  mark_neighbors_and_parents(node->parent);
}


void QuadTree::consolidate_node(QuadTreeNode* node) {

  // Need to
  //  Move leaf child dofs into my original box
  //  erase all descendents from levels
  //  delete immediate descentdents

  node->src_dof_lists.original_box.clear();

  // This can be parallelized
  std::vector<QuadTreeNode*> remove_from_lvl;
  std::vector<QuadTreeNode*> queue;
  queue.push_back(node);
  for (int i = 0; i < queue.size(); i++) {
    QuadTreeNode* current = queue[i];
    if (current->is_leaf) {
      node->src_dof_lists.original_box.insert(
        node->src_dof_lists.original_box.end(),
        current->src_dof_lists.original_box.begin(),
        current->src_dof_lists.original_box.end());
    } else {
      for (QuadTreeNode* child : current->children) {
        queue.push_back(child);
      }
    }
    if (current != node) {
      remove_from_lvl.push_back(current);
    }
  }

  for (QuadTreeNode* erase : remove_from_lvl) {
    QuadTreeLevel* erase_level = levels[erase->level];
    for (int i = 0; i < erase_level->nodes.size(); i++) {
      if (erase_level->nodes[i]->id == erase->id) {
        erase_level->nodes.erase(erase_level->nodes.begin() + i);
        break;
      }
    }
  }

  node->tl = nullptr;
  node->tr = nullptr;
  node->bl = nullptr;
  node->br = nullptr;
  node->children[0] = nullptr;
  node->children[1] = nullptr;
  node->children[2] = nullptr;
  node->children[3] = nullptr;
  node->is_leaf = true;
}



// BUG : this might confuse interior holes if they get too close to each other
void QuadTree::perturb(const Boundary & perturbed_boundary) {
  // 1) create mapping, storing vectors of additions/deletions
  // 2) go to every node, marking those with additions and deletions
  // these are vectors of point indices (p_0, p_1, etc)
  std::vector<double> additions;
  std::vector<double> deletions;
  // first we do it with vectors, optimize later
  std::vector<double> old_points = boundary->points;
  std::vector<double> new_points = perturbed_boundary.points;
  // now create mapping of new_points to their point index in the new vec
  std::unordered_map<pair, int, boost::hash<pair>> point_to_new_index;

  // double a = omp_get_wtime();
  for (int i = 0; i < new_points.size(); i += 2) {
    pair new_point(new_points[i], new_points[i + 1]);
    point_to_new_index[new_point] = i / 2;
  }

  std::vector<bool> found_in_old(new_points.size() / 2);
  for (int i = 0; i < found_in_old.size(); i++) {
    found_in_old[i] = false;
  }
  // Mapping from point index in old points vec to point index in new points vec
  std::unordered_map<int, int> old_index_to_new_index;
  for (int i = 0; i < old_points.size(); i += 2) {
    pair old_point(old_points[i], old_points[i + 1]);
    // Is this point also in the new points vec?
    std::unordered_map<pair, int, boost::hash<pair>>::const_iterator element =
          point_to_new_index.find(old_point);
    if (element != point_to_new_index.end()) {
      old_index_to_new_index[i / 2] = element->second;
      found_in_old[element->second] = true;
    } else {
      // If it's in the old vec but not the new vec, it was deleted
      deletions.push_back(i / 2);
    }
  }
  for (int i = 0; i < found_in_old.size(); i++) {
    if (!found_in_old[i]) {
      additions.push_back(i);
    }
  }

  int num_compressed = 0;
  int num_total = 0;
  for (QuadTreeLevel* level : levels) {
    for (QuadTreeNode* node : level->nodes) {
      num_total++;
      if (node->compressed) {
        num_compressed++;
      }
    }
  }

  // std::cout << "Before perturb, " << num_compressed << " of " << num_total <<
  //           " are compressed." << std::endl;
  // go through all leaf original box vectors and apply mapping.
  // (if there is a deletion it will be processed later)
  // each node will be one of three things
  //   1) unmarked, in which case the below is a perfectly good mapping
  //   2) marked non-leaf, the below is irrelevant, everything will be dumped
  //   3) marked leaf, only the leaf portion of the below is relevant.
  for (QuadTreeLevel* level : levels) {
    for (QuadTreeNode* node : level->nodes) {
      std::vector<int> ob, ab, s, r, sn, n;
      if (node->is_leaf) {
        for (int idx : node->src_dof_lists.original_box) {
          int point_index = idx / solution_dimension;
          std::unordered_map<int, int>::const_iterator element =
            old_index_to_new_index.find(point_index);
          if (element != old_index_to_new_index.end()) {
            ob.push_back(solution_dimension * element->second
                         + idx % solution_dimension);
          }
        }
        node->src_dof_lists.original_box = ob;
      }
      for (int idx : node->src_dof_lists.active_box) {
        int point_index = idx / solution_dimension;
        std::unordered_map<int, int>::const_iterator element =
          old_index_to_new_index.find(point_index);
        if (element != old_index_to_new_index.end()) {
          ab.push_back(solution_dimension * element->second
                       + idx % solution_dimension);
        }
      }
      for (int idx : node->src_dof_lists.skel) {
        int point_index = idx / solution_dimension;
        std::unordered_map<int, int>::const_iterator element =
          old_index_to_new_index.find(point_index);
        if (element != old_index_to_new_index.end()) {
          s.push_back(solution_dimension * element->second
                      + idx % solution_dimension);
        }
      }
      for (int idx : node->src_dof_lists.redundant) {
        int point_index = idx / solution_dimension;
        std::unordered_map<int, int>::const_iterator element =
          old_index_to_new_index.find(point_index);
        if (element != old_index_to_new_index.end()) {
          r.push_back(solution_dimension * element->second
                      + idx % solution_dimension);
        }
      }
      for (int idx : node->src_dof_lists.skelnear) {
        int point_index = idx / solution_dimension;
        std::unordered_map<int, int>::const_iterator element =
          old_index_to_new_index.find(point_index);
        if (element != old_index_to_new_index.end()) {
          sn.push_back(solution_dimension * element->second
                       + idx % solution_dimension);
        }
      }
      node->src_dof_lists.active_box = ab;
      node->src_dof_lists.skel = s;
      node->src_dof_lists.skelnear = sn;
      node->src_dof_lists.redundant = r;
    }
  }

  // go through all additions, find their leaves, make addition and call mark
  // function
  std::vector<QuadTreeNode*> maybe_bursting;
  for (int i = 0; i < additions.size(); i++) {
    double newx = new_points[2 * additions[i]];
    double newy = new_points[2 * additions[i] + 1];
    QuadTreeNode* current = root;
    while (!current->is_leaf) {
      double midx = ((current->corners[6] - current->corners[0]) / 2.0)
                    + current->corners[0];
      double midy = ((current->corners[3] - current->corners[1]) / 2.0)
                    + current->corners[1];
      if (newx < midx && newy < midy) {
        current = current->bl;
      } else if (newx < midx && newy >= midy) {
        current = current->tl;
      } else if (newx >= midx && newy < midy) {
        current = current->br;
      } else {
        current = current->tr;
      }
    }
    for (int j = 0; j < solution_dimension; j++) {
      current->src_dof_lists.original_box.push_back(solution_dimension
          * additions[i] + j);
    }
    maybe_bursting.push_back(current);
    mark_neighbors_and_parents(current);
  }

  for (QuadTreeLevel* level : levels) {
    for (QuadTreeNode* node : level->nodes) {
      if (node->is_leaf) {
        node->dofs_below = node->src_dof_lists.original_box.size();
      } else {
        node->dofs_below = 0;
      }
    }
  }
  for (int l = levels.size() - 1; l >= 1; l--) {
    QuadTreeLevel* level = levels[l];
    for (QuadTreeNode* node : level->nodes) {
      node->parent->dofs_below += node->dofs_below;
    }
  }

  // go through all deletions, find their leaves, make deletion and call mark
  // function
  std::unordered_map<QuadTreeNode*, bool> sparse;

  for (int i = 0; i < deletions.size(); i++) {
    double oldx = old_points[2 * deletions[i]];
    double oldy = old_points[2 * deletions[i] + 1];
    QuadTreeNode* current = root;
    bool path_marked = false;
    while (!current->is_leaf) {
      if (current->dofs_below < MAX_LEAF_DOFS && !path_marked) {
        path_marked = true;
        sparse[current] = true;
      }
      double midx = ((current->corners[6] - current->corners[0]) / 2.0)
                    + current->corners[0];
      double midy = ((current->corners[3] - current->corners[1]) / 2.0)
                    + current->corners[1];
      if (oldx < midx && oldy < midy) {
        current = current->bl;
      } else if (oldx < midx && oldy >= midy) {
        current = current->tl;
      } else if (oldx >= midx && oldy < midy) {
        current = current->br;
      } else {
        current = current->tr;
      }
    }

    mark_neighbors_and_parents(current);
  }

  // TODO(John) this should be a boundary copy routine
  boundary->points = perturbed_boundary.points;
  boundary->normals = perturbed_boundary.normals;
  boundary->weights = perturbed_boundary.weights;
  boundary->curvatures = perturbed_boundary.curvatures;
  boundary->boundary_values = perturbed_boundary.boundary_values;
  boundary->perturbation_parameters = perturbed_boundary.perturbation_parameters;
  boundary->holes = perturbed_boundary.holes;

  // If any nodes are bursting now, subdivide them.
  for (QuadTreeNode* node : maybe_bursting) {
    if (node->is_leaf
        && node->src_dof_lists.original_box.size() +
        node->tgt_dof_lists.original_box.size() > MAX_LEAF_DOFS) {
      node_subdivide(node);
    }
  }

  // If we can consolidate nodes into their parent, do that.

  for (auto it = sparse.begin(); it != sparse.end(); ++it) {
    consolidate_node(it->first);
  }

  num_compressed = 0;
  num_total = 0;
  for (QuadTreeLevel* level : levels) {
    for (QuadTreeNode* node : level->nodes) {
      node->neighbors.clear();
      num_total++;
      if (node->compressed) {
        num_compressed++;
      }
    }
  }
  // For now, just recalculate neighbors the same as above.
  for (int level = 0; level < levels.size(); level++) {
    QuadTreeLevel* current_level = levels[level];
    for (int k = 0; k < current_level->nodes.size(); k++) {
      QuadTreeNode* node_a = current_level->nodes[k];

      // Each node is neighbors with all its siblings
      if (node_a->parent != nullptr) {
        for (QuadTreeNode* sibling : node_a->parent->children) {
          if (sibling->id != node_a->id) node_a->neighbors.push_back(sibling);
        }

        // Now check all parents' neighbors' children
        for (QuadTreeNode* parents_neighbor : node_a->parent->neighbors) {
          for (QuadTreeNode* cousin : parents_neighbor->children) {
            if (cousin == nullptr) continue;
            if (cousin->level != node_a->level) continue;
            double dist = sqrt(pow(node_a->corners[0] - cousin->corners[0], 2)
                               + pow(node_a->corners[1] - cousin->corners[1], 2));
            // just need to check if the distance of the bl corners
            // is <=s*sqrt(2)
            if (dist < node_a->side_length * sqrt(2) + 1e-5) {
              node_a->neighbors.push_back(cousin);
            }
          }
        }
      }

      // now if it is a leaf, check against nodes in all subsequent levels
      if (node_a->is_leaf) {
        for (int n = 0; n < node_a->neighbors.size(); n++) {
          QuadTreeNode* neighbor =  node_a->neighbors[n];
          // make sure this isn't a neighbor from a higher level
          if (neighbor->level != node_a->level) {
            continue;
          }
          for (QuadTreeNode* child : neighbor->children) {
            if (child != nullptr) {
              get_descendent_neighbors(node_a, child);
            }
          }
        }
      }
    }
  }

  // std::cout << "After perturb, " << num_compressed << " of " << num_total <<
  //           " are compressed." << std::endl;
}


// All that is not copied is id, tltrblbr, parent, children, and neighbors
// TODO(John) see if copying over id would be fine, I feel like it would be
void copy_info(QuadTreeNode* old_node, QuadTreeNode* new_node) {
  new_node->level = old_node->level;
  new_node->dofs_below = old_node->dofs_below;

  new_node->is_leaf = old_node->is_leaf;
  // DEBUGGING
  new_node->id = old_node->id;
  new_node->X_rr_is_LU_factored = old_node->X_rr_is_LU_factored;
  new_node->compressed = old_node->compressed;

  new_node->side_length = old_node->side_length;
  new_node->compression_ratio = old_node->compression_ratio;
  new_node->compress_time = old_node->compress_time;

  new_node->src_dof_lists = old_node->src_dof_lists;
  new_node->tgt_dof_lists = old_node->tgt_dof_lists;

  new_node->T = old_node->T;
  new_node->L = old_node->L;
  new_node->U = old_node->U;
  new_node->X_rr = old_node->X_rr;
  new_node->schur_update = old_node->schur_update;
  new_node->X_rr_lu = old_node->X_rr_lu;
  new_node->X_rr_piv = old_node->X_rr_piv;

  new_node->src_T = old_node->src_T;
  new_node->tgt_T = old_node->tgt_T;
  new_node->X_rs = old_node->X_rs;
  new_node->X_sr = old_node->X_sr;
  for (int i = 0; i < 8; i++) {
    new_node->corners[i] = old_node->corners[i];
  }
}


void QuadTree::copy_into(QuadTree* new_tree) const {
  // The strategy here is going to be to create a new node for every old node,
  // then keep a mapping from new to old. With that, we'll copy all the data
  // over, including connections, levels, and matrices.
  *new_tree = QuadTree();
  std::vector < QuadTreeNode*> new_nodes;
  std::unordered_map<QuadTreeNode*, QuadTreeNode*> old_to_new;
  std::unordered_map<QuadTreeNode*, QuadTreeNode*> new_to_old;

  for (int lvl = 0; lvl < levels.size(); lvl++) {
    QuadTreeLevel* level = levels[lvl];
    for (int n = 0; n < level->nodes.size(); n++) {
      QuadTreeNode* old_node = level->nodes[n];
      QuadTreeNode* new_node = new QuadTreeNode();
      new_nodes.push_back(new_node);
      copy_info(old_node, new_node);
      old_to_new[old_node] = new_node;
      new_to_old[new_node] = old_node;
    }
  }

  for (int n = 0; n < new_nodes.size(); n++) {
    QuadTreeNode* new_node = new_nodes[n];
    new_node->parent = old_to_new[new_to_old[new_node]->parent];
    if (!new_to_old[new_node]->is_leaf) {
      new_node->tl = old_to_new[new_to_old[new_node]->tl];
      new_node->tr = old_to_new[new_to_old[new_node]->tr];
      new_node->bl = old_to_new[new_to_old[new_node]->bl];
      new_node->br = old_to_new[new_to_old[new_node]->br];
      for (int c = 0; c < 4; c++) {
        new_node->children[c] =  old_to_new[new_to_old[new_node]->children[c]];
      }
    }
    for (int nbr = 0; nbr < new_to_old[new_node]->neighbors.size(); nbr++) {
      QuadTreeNode* neighbor = new_to_old[new_node]->neighbors[nbr];
      new_node->neighbors.push_back(old_to_new[neighbor]);
    }

  }

  new_tree->root = old_to_new[root];

  new_tree->solution_dimension = solution_dimension;
  new_tree->domain_dimension = domain_dimension;

  new_tree->no_proxy_level = no_proxy_level;

  new_tree->min = min;
  new_tree->max = max;

  new_tree->domain_points = domain_points;


  new_tree->allskel_mat = allskel_mat;
  new_tree->allskel_mat_lu = allskel_mat_lu;
  new_tree->U = U;
  new_tree->Psi = Psi;
  new_tree->S_LU = S_LU;
  new_tree->allskel_mat_piv = allskel_mat_piv;
  new_tree->S_piv = S_piv;
  // ie_Mat allskel_mat, allskel_mat_lu, U, Psi, S_LU;
  // std::vector<lapack_int> allskel_mat_piv, S_piv;


  // Boundary* boundary;
  // std::vector<QuadTreeLevel*> levels;

// TODO(John) again, should be replaced with good copy routine
  // For now, a hack - we'll initialize it, then completely redefine it.
  // smart pointer, make a circle or something (or switch like elsewhere)

  if (new_tree->boundary) {
    delete new_tree->boundary;
  }

  switch (boundary->boundary_shape) {
    case Boundary::BoundaryShape::CIRCLE:
      new_tree->boundary = new Circle();
      break;
    case Boundary::BoundaryShape::ROUNDED_SQUARE:
      new_tree->boundary = new RoundedSquare();
      break;
    case Boundary::BoundaryShape::ANNULUS:
      new_tree->boundary = new Annulus();
      break;
    case Boundary::BoundaryShape::DONUT:
      new_tree->boundary = new Donut();
      break;
    case Boundary::BoundaryShape::CUBIC_SPLINE:
      new_tree->boundary = new CubicSpline();
      break;
    case Boundary::BoundaryShape::EX1:
      new_tree->boundary = new Ex1Boundary();
      break;
    case Boundary::BoundaryShape::EX2:
      new_tree->boundary = new Ex2Boundary();
      break;
    case Boundary::BoundaryShape::EX3:
      new_tree->boundary = new Ex3Boundary();
      break;
  }

  new_tree->boundary->points = boundary->points;
  new_tree->boundary->normals = boundary->normals;
  new_tree->boundary->weights = boundary->weights;
  new_tree->boundary->curvatures = boundary->curvatures;
  new_tree->boundary->boundary_values = boundary->boundary_values;
  new_tree->boundary->boundary_shape = boundary->boundary_shape;
  new_tree->boundary->perturbation_parameters =
    boundary->perturbation_parameters;
  new_tree->boundary->holes = boundary->holes;
  new_tree->boundary->num_outer_nodes = boundary->num_outer_nodes;

  for (int lvl = 0; lvl < levels.size(); lvl++) {
    QuadTreeLevel* old_level = levels[lvl];
    QuadTreeLevel* new_level = new QuadTreeLevel();
    for (int n = 0; n < old_level->nodes.size(); n++) {
      new_level->nodes.push_back(old_to_new[old_level->nodes[n]]);
    }
    new_tree->levels.push_back(new_level);
  }
}


void QuadTree::reset() {
  if (root) {
    delete root;
  }
  for (QuadTreeLevel* level : levels) {
    if (level) {
      delete level;
    }
  }
  levels.clear();
  initialize_tree(boundary, domain_points, solution_dimension,
                  domain_dimension);
}


void QuadTree::reset(Boundary * boundary_) {
  if (root) {
    delete root;
  }
  for (QuadTreeLevel* level : levels) {
    if (level) {
      delete level;
    }
  }
  levels.clear();
  initialize_tree(boundary_, domain_points, solution_dimension,
                  domain_dimension);
}


void QuadTree::remove_inactive_dofs_at_all_boxes() {
  int lvls = levels.size();
  for (int level = lvls - 1; level >= 0; level--) {
    remove_inactive_dofs_at_level(level);
  }
}


void QuadTree::remove_inactive_dofs_at_level(int level) {
  QuadTreeLevel* current_level = levels[level];
  // First, get all active dofs from children
  for (QuadTreeNode * node : current_level->nodes) {
    if (node->compressed) continue;
    remove_inactive_dofs_at_box(node);
  }
  // Next, get all active near dofs from neighbors
  for (QuadTreeNode* node_a : current_level->nodes) {
    // if (node_a->compressed) continue;
    node_a->src_dof_lists.near.clear();
    for (QuadTreeNode* neighbor : node_a->neighbors) {
      // Some neighbors are smaller boxes from higher levels, we don't
      // care about those, their parents have the updated information.
      if (neighbor->level > node_a->level) {
        continue;
      }
      if (neighbor->is_leaf) {
        for (int idx : neighbor->src_dof_lists.original_box) {
          node_a->src_dof_lists.near.push_back(idx);
        }
      } else {
        for (int idx : neighbor->src_dof_lists.active_box) {
          node_a->src_dof_lists.near.push_back(idx);
        }
      }
    }
  }
}


void QuadTree::remove_inactive_dofs_at_box(QuadTreeNode* node) {
  // this function removes from the box any DoFs which have already been made
  // redundant. It involves a bunch of annoying C++ functions and probably
  // would look nicer in matlab.

  // populate active_box
  node->src_dof_lists.skel.clear();
  node->src_dof_lists.skelnear.clear();
  node->src_dof_lists.redundant.clear();
  node->src_dof_lists.active_box.clear();

  if (!node->is_leaf) {
    for (QuadTreeNode* child : node->children) {
      if (child->compressed) {
        for (int i : child->src_dof_lists.skel) {
          node->src_dof_lists.active_box.push_back(i);
        }
      } else {
        for (int i : child->src_dof_lists.active_box) {
          node->src_dof_lists.active_box.push_back(i);
        }
      }
    }
  } else {
    node->src_dof_lists.active_box = node->src_dof_lists.original_box;
  }

  node->tgt_dof_lists.skel.clear();
  node->tgt_dof_lists.skelnear.clear();
  node->tgt_dof_lists.redundant.clear();
  node->tgt_dof_lists.active_box.clear();
  if (!node->is_leaf) {
    for (QuadTreeNode* child : node->children) {
      if (child->compressed) {
        for (int i : child->tgt_dof_lists.skel) {
          node->tgt_dof_lists.active_box.push_back(i);
        }
      } else {
        for (int i : child->tgt_dof_lists.active_box) {
          node->tgt_dof_lists.active_box.push_back(i);
        }
      }
    }
  } else {
    node->tgt_dof_lists.active_box = node->tgt_dof_lists.original_box;
  }
}

}  // namespace ie_solver
