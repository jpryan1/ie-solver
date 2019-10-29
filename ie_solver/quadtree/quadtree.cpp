// Copyright 2019 John Paul Ryan
#include <omp.h>
#include <fstream>
#include <cassert>
#include <cmath>
#include <unordered_map>
#include <iostream>
#include <utility>
#include <boost/functional/hash.hpp>
#include "ie_solver/log.h"
#include "ie_solver/quadtree/quadtree.h"

namespace ie_solver {

typedef std::pair<double, double> pair;
unsigned int QuadTreeNode::id_count = 0;

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
  double tree_start_time = omp_get_wtime();
  assert(boundary_->points.size() > 0
         && "number of boundary->points to init tree cannot be 0.");
  // todo(john) later we can assert that the boundary->points are a multiple of
  // dimension
  this->boundary = boundary_;
  this->domain_points = domain_points_;
  this->solution_dimension = solution_dimension_;
  this->domain_dimension = domain_dimension_;
  QuadTreeNode::id_count = 0;
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
  // all near and far ranges should be taken care of by this guy
  for (unsigned int i = 0; i < boundary->points.size(); i += domain_dimension) {
    recursive_add(this->root, boundary->points[i], boundary->points[i + 1],
                  i / domain_dimension, true);
  }
  for (unsigned int i = 0; i < domain_points.size(); i += domain_dimension) {
    recursive_add(this->root, domain_points[i], domain_points[i + 1],
                  i / domain_dimension, false);
  }
  // make neighbor lists in a stupid way
  for (unsigned int level = 0; level < levels.size(); level++) {
    QuadTreeLevel* current_level = levels[level];
    for (unsigned int k = 0; k < current_level->nodes.size(); k++) {
      QuadTreeNode* node_a = current_level->nodes[k];
      // if (level > no_proxy_level && node_a->src_dof_lists.original_box.size()
      //     > 0.25 * solution_dimension *
      //     (boundary->points.size() / domain_dimension)) {
      //   no_proxy_level = level;
      //   std::cout << "No proxy level "
      //             << no_proxy_level << " threshold " <<
      //             0.25 * solution_dimension *
      //             (boundary->points.size() /
      //               domain_dimension)  << " total " <<
      //             node_a->src_dof_lists.original_box.size()
      //             << std::endl;
      // } else {
      //   std::cout << level << " " << no_proxy_level << std::endl;
      //   std::cout << node_a->src_dof_lists.original_box.size()  <<
      //  " > " << 0.25 *
      //             solution_dimension * (boundary->points.size() /
      //                                   domain_dimension) << std::endl;
      // }
      // first check against all nodes on this level
      for (unsigned int l = k + 1; l < current_level->nodes.size(); l++) {
        QuadTreeNode* node_b = current_level->nodes[l];
        double dist = sqrt(pow(node_a->corners[0] - node_b->corners[0], 2)
                           + pow(node_a->corners[1] - node_b->corners[1], 2));
        // just need to check if the distance of the bl corners
        // is <=s*sqrt(2)
        if (dist < node_a->side_length * sqrt(2) + 1e-5) {
          node_a->neighbors.push_back(node_b);
          node_b->neighbors.push_back(node_a);
        }
      }
      // now if it is a leaf, check against nodes in all subsequent levels
      if (node_a->is_leaf) {
        for (unsigned int n = 0; n < node_a->neighbors.size(); n++) {
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
  double tree_end_time = omp_get_wtime();
  std::cout << "timing: tree_init " << (tree_end_time - tree_start_time) <<
            std::endl;
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
                             unsigned int point_ind, bool is_boundary) {
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
  // in the stokes case, the original_box contains pairs of consecutive
  // indices into the kernel matrix. therefore, its size had better be even
  // Now we bring the indices from the parent's box down into its childrens
  // boxes
  for (unsigned int index = 0; index < node->src_dof_lists.original_box.size();
       index += solution_dimension) {
    unsigned int matrix_index = node->src_dof_lists.original_box[index];
    unsigned int points_vec_index = (matrix_index / solution_dimension) *
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
  for (unsigned int index = 0; index < node->tgt_dof_lists.original_box.size();
       index += solution_dimension) {
    unsigned int matrix_index = node->tgt_dof_lists.original_box[index];
    unsigned int points_vec_index = (matrix_index / solution_dimension) *
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
  node->compressed = false;
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
  for (unsigned int i = 0; i < new_points.size(); i += 2) {
    pair new_point(new_points[i], new_points[i + 1]);
    point_to_new_index[new_point] = i / 2;
  }
  std::vector<bool> found_in_old(new_points.size() / 2);
  for (unsigned int i = 0; i < found_in_old.size(); i++) {
    found_in_old[i] = false;
  }
  // Mapping from point index in old points vec to point index in new points vec
  std::unordered_map<int, int> old_index_to_new_index;
  for (unsigned int i = 0; i < old_points.size(); i += 2) {
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
  for (unsigned int i = 0; i < found_in_old.size(); i++) {
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
  std::cout << "Before perturb, " << num_compressed << " of " << num_total <<
            " are compressed." << std::endl;
  // TODO(John) the below needs to be changed for stokes
  // go through all leaf original box vectors and apply mapping.
  // (if there is a deletion it will be processed later)
  // each node will be one of three things
  //   1) unmarked, in which case the below is a perfectly good mapping
  //   2) marked non-leaf, the below is irrelevant, everything will be dumped
  //   3) marked leaf, only the leaf portion of the below is relevant.
  for (QuadTreeLevel* level : levels) {
    for (QuadTreeNode* node : level->nodes) {
      std::vector<unsigned int> ob, ab, s, r, sn, n;
      if (node->is_leaf) {
        for (unsigned int idx : node->src_dof_lists.original_box) {
          unsigned int point_index = idx / solution_dimension;
          std::unordered_map<int, int>::const_iterator element =
            old_index_to_new_index.find(point_index);
          if (element != old_index_to_new_index.end()) {
            ob.push_back(solution_dimension * element->second
                         + idx % solution_dimension);
          }
        }

        node->src_dof_lists.original_box = ob;
      }
      for (unsigned int idx : node->src_dof_lists.active_box) {
        unsigned int point_index = idx / solution_dimension;
        std::unordered_map<int, int>::const_iterator element =
          old_index_to_new_index.find(point_index);
        if (element != old_index_to_new_index.end()) {
          ab.push_back(solution_dimension * element->second
                       + idx % solution_dimension);
        }
      }
      for (unsigned int idx : node->src_dof_lists.skel) {
        unsigned int point_index = idx / solution_dimension;
        std::unordered_map<int, int>::const_iterator element =
          old_index_to_new_index.find(point_index);
        if (element != old_index_to_new_index.end()) {
          s.push_back(solution_dimension * element->second
                      + idx % solution_dimension);
        }
      }
      for (unsigned int idx : node->src_dof_lists.redundant) {
        unsigned int point_index = idx / solution_dimension;
        std::unordered_map<int, int>::const_iterator element =
          old_index_to_new_index.find(point_index);
        if (element != old_index_to_new_index.end()) {
          r.push_back(solution_dimension * element->second
                      + idx % solution_dimension);
        }
      }
      for (unsigned int idx : node->src_dof_lists.skelnear) {
        unsigned int point_index = idx / solution_dimension;
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
  std::vector<bool> found(additions.size());
  for (int i = 0; i < additions.size(); i++) {
    found[i] = false;
  }
  for (QuadTreeLevel* level : levels) {
    for (QuadTreeNode* node : level->nodes) {
      if (!node->is_leaf) {
        continue;
      }
      for (unsigned int i = 0; i < additions.size(); i++) {
        double difx = new_points[2 * additions[i]] - node->corners[0];
        double dify = new_points[2 * additions[i] + 1] - node->corners[1];
        if (!found[i] && difx < node->side_length && dify < node->side_length
            && difx >= 0 && dify >= 0) {
          found[i] = true;
          for (int j = 0; j < solution_dimension; j++) {
            node->src_dof_lists.original_box.push_back(solution_dimension
                * additions[i] + j);
          }
          mark_neighbors_and_parents(node);
        }
      }
    }
  }
// go through all deletions, find their leaves, make deletion and call mark
// function
  for (QuadTreeLevel* level : levels) {
    for (QuadTreeNode* node : level->nodes) {
      if (!node->is_leaf) {
        continue;
      }
      for (unsigned int i = 0; i < deletions.size(); i++) {
        double difx = old_points[2 * deletions[i]] - node->corners[0];
        double dify = old_points[2 * deletions[i] + 1] - node->corners[1];
        if (difx < node->side_length && dify < node->side_length && difx >= 0
            && dify >= 0) {
          mark_neighbors_and_parents(node);
        }
      }
    }
  }
  boundary->points = perturbed_boundary.points;
  boundary->normals = perturbed_boundary.normals;
  boundary->weights = perturbed_boundary.weights;
  boundary->curvatures = perturbed_boundary.curvatures;
  boundary->boundary_values = perturbed_boundary.boundary_values;
  boundary->perturbation_parameters[0] =
    perturbed_boundary.perturbation_parameters[0];
  boundary->holes = perturbed_boundary.holes;

  // If any nodes are bursting now, subdivide them.
  for (int l = levels.size() - 1; l >= 0; l--) {
    QuadTreeLevel* level = levels[l];
    for (QuadTreeNode* node : level->nodes) {
      if (node->is_leaf
          && node->src_dof_lists.original_box.size() +
          node->tgt_dof_lists.original_box.size() > MAX_LEAF_DOFS) {
        for (int hh = 0; hh < 4; hh++) {
          if (node->children[hh] != nullptr) {
            std::cout << "leaf with children?" << std::endl;
            exit(0);
          }
        }
        node_subdivide(node);
      }
    }
  }
  // If we can consolidate nodes into their parent, do that.
  for (int l = levels.size() - 1; l >= 0; l--) {
    QuadTreeLevel* level = levels[l];
    for (QuadTreeNode* node : level->nodes) {
      if (!node->is_leaf) {
        bool all_children_leaves = true;

        for (QuadTreeNode* child : node->children) {
          if (!child->is_leaf) {
            all_children_leaves = false;
          }
        }
        if (all_children_leaves) {
          int num_child_dofs = 0;
          for (QuadTreeNode* child : node->children) {
            num_child_dofs += child->src_dof_lists.original_box.size();
          }
          if (num_child_dofs < MAX_LEAF_DOFS) {
            assert(!node->compressed);
            node->src_dof_lists.original_box.clear();
            for (QuadTreeNode* child : node->children) {
              for (unsigned int idx : child->src_dof_lists.original_box) {
                node->src_dof_lists.original_box.push_back(idx);
              }
            }

            for (QuadTreeNode* child : node->children) {
              QuadTreeLevel* child_level = levels[child->level];
              for (int i = 0; i < child_level->nodes.size(); i++) {
                if (child_level->nodes[i]->id == child->id) {
                  child_level->nodes.erase(child_level->nodes.begin() + i);
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
        }
      }
    }
  }

  for (QuadTreeLevel* level : levels) {
    for (QuadTreeNode* node : level->nodes) {
      if (!node->is_leaf) {
        int num_child_dofs = 0;
        for (QuadTreeNode* child : node->children) {
          num_child_dofs += child->src_dof_lists.original_box.size();
        }
      }
    }
  }

  num_compressed = 0;
  num_total = 0;
  for (QuadTreeLevel* level : levels) {
    for (QuadTreeNode* node : level->nodes) {
      num_total++;
      if (node->compressed) {
        num_compressed++;
      }
    }
  }
  std::cout << "After perturb, " << num_compressed << " of " << num_total <<
            " are compressed." << std::endl;
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
  QuadTreeNode::id_count = 0;
  initialize_tree(boundary, domain_points, solution_dimension,
                  domain_dimension);
}


void QuadTree::reset(Boundary * boundary_) {
  if (root) {delete root; } for (QuadTreeLevel* level : levels) {
    if (level) {
      delete level;
    }
  }
  levels.clear();
  QuadTreeNode::id_count = 0;
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
    if (node_a->compressed) continue;
    node_a->src_dof_lists.near.clear();
    for (QuadTreeNode* neighbor : node_a->neighbors) {
      // Some neighbors are smaller boxes from higher levels, we don't
      // care about those, their parents have the updated information.
      if (neighbor->level > node_a->level) {
        continue;
      }
      for (unsigned int idx : neighbor->src_dof_lists.active_box) {
        node_a->src_dof_lists.near.push_back(idx);
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
        for (unsigned int i : child->src_dof_lists.skel) {
          node->src_dof_lists.active_box.push_back(i);
        }
      } else {
        for (unsigned int i : child->src_dof_lists.active_box) {
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
        for (unsigned int i : child->tgt_dof_lists.skel) {
          node->tgt_dof_lists.active_box.push_back(i);
        }
      } else {
        for (unsigned int i : child->tgt_dof_lists.active_box) {
          node->tgt_dof_lists.active_box.push_back(i);
        }
      }
    }
  } else {
    node->tgt_dof_lists.active_box = node->tgt_dof_lists.original_box;
  }
}

}  // namespace ie_solver
