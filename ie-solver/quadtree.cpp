// Copyright 2019 John Paul Ryan
#include <omp.h>
#include <fstream>
#include <cassert>
#include <cmath>
#include <unordered_map>
#include <iostream>
#include <utility>
#include <boost/functional/hash.hpp>
#include "ie-solver/log.h"
#include "ie-solver/quadtree.h"

namespace ie_solver {

typedef std::pair<double, double> pair;

unsigned int QuadTreeNode::id_count = 0;

void QuadTree::initialize_tree(Boundary* boundary_, int solution_dimension_) {
  assert(boundary_->points.size() > 0
         && "number of boundary->points to init tree cannot be 0.");
  // todo(john) later we can assert that the boundary->points are a multiple of
  // dimension
  this->boundary = boundary_;
  this->solution_dimension = solution_dimension_;
  QuadTreeNode::id_count = 0;

  min = boundary->points[0];
  max = boundary->points[0];

  double tree_min = 0;  // 1 - tree_max_;//0.05;
  double tree_max = 1.0;  // m_pi;

  // double tree_min = -m_pi;
  // double tree_max = m_pi;

  for (double point : boundary->points) {
    if (point < min) min = point;
    if (point > max) max = point;
  }
  assert(min > tree_min && max < tree_max);

  root = new QuadTreeNode();
  root->level = 0;
  root->parent = nullptr;
  // root is the box [0,1]x[0,1]
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
  for (unsigned int i = 0; i < boundary->points.size(); i += 2) {
    recursive_add(this->root, boundary->points[i], boundary->points[i + 1],
                  i / 2);
  }
  // make neighbor lists in a stupid way

  for (unsigned int level = 0; level < levels.size(); level++) {
    QuadTreeLevel* current_level = levels[level];
    for (unsigned int k = 0; k < current_level->nodes.size(); k++) {
      QuadTreeNode* node_a = current_level->nodes[k];

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

  for (unsigned int j = 0; j < levels.size(); j++) {
    QuadTreeLevel* current_level = levels[j];
    for (QuadTreeNode* node_a : current_level->nodes) {
      for (QuadTreeNode* neighbor : node_a->neighbors) {
        neighbor->interaction_lists.near.insert(
          neighbor->interaction_lists.near.end(),
          node_a->interaction_lists.original_box.begin(),
          node_a->interaction_lists.original_box.end());
      }
    }
  }
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
                             unsigned int point_ind) {
  assert(node != nullptr && "recursive_add fails on null node.");
  for (int i = 0; i < solution_dimension; i++) {
    node->interaction_lists.original_box.push_back(
      solution_dimension * point_ind + i);
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
    recursive_add(child, x, y, point_ind);
  } else {
    // do we need one?
    // if this node is exploding and needs children
    if (node->is_leaf
        && node->interaction_lists.original_box.size() > MAX_LEAF_DOFS) {
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
  if (solution_dimension == 2) { // For the love of god come back to this
    assert(node->interaction_lists.original_box.size() % 2 == 0);

    for (unsigned int index = 0;
         index < node->interaction_lists.original_box.size(); index += 2) {
      unsigned int stokes_index = node->interaction_lists.original_box[index];
      double x = boundary->points[stokes_index];
      double y = boundary->points[stokes_index + 1];
      if (x < midx && y < midy) {
        bl->interaction_lists.original_box.push_back(stokes_index);
        bl->interaction_lists.original_box.push_back(stokes_index + 1);
      } else if (x < midx && y >= midy) {
        tl->interaction_lists.original_box.push_back(stokes_index);
        tl->interaction_lists.original_box.push_back(stokes_index + 1);
      } else if (x >= midx && y < midy) {
        br->interaction_lists.original_box.push_back(stokes_index);
        br->interaction_lists.original_box.push_back(stokes_index + 1);
      } else {
        tr->interaction_lists.original_box.push_back(stokes_index);
        tr->interaction_lists.original_box.push_back(stokes_index + 1);
      }
    }
  } else {
    // Laplace case
    for (unsigned int index : node->interaction_lists.original_box) {
      double x = boundary->points[2 * index];
      double y = boundary->points[2 * index + 1];
      if (x < midx && y < midy) {
        bl->interaction_lists.original_box.push_back(index);
      } else if (x < midx && y >= midy) {
        tl->interaction_lists.original_box.push_back(index);
      } else if (x >= midx && y < midy) {
        br->interaction_lists.original_box.push_back(index);
      } else {
        tr->interaction_lists.original_box.push_back(index);
      }
    }
  }

  if (bl->interaction_lists.original_box.size() > MAX_LEAF_DOFS) {
    node_subdivide(bl);
  }
  if (tl->interaction_lists.original_box.size() > MAX_LEAF_DOFS) {
    node_subdivide(tl);
  }
  if (tr->interaction_lists.original_box.size() > MAX_LEAF_DOFS) {
    node_subdivide(tr);
  }
  if (br->interaction_lists.original_box.size() > MAX_LEAF_DOFS) {
    node_subdivide(br);
  }
}


void QuadTree::write_quadtree_to_file() {
  std::ofstream output;
  output.open("output/data/ie_solver_tree.txt");
  if (output.is_open()) {
    for (QuadTreeLevel* level : levels) {
      for (QuadTreeNode* node : level->nodes) {
        output << node->corners[0] << "," << node->corners[1]
               << "," << node->side_length << "," << node->id << "," <<
               node->interaction_lists.original_box.size() <<
               std::endl;
        //}
      }
    }
    output.close();
  } else {
    LOG::ERROR("failed to open output file!");
  }
}

void QuadTree::mark_neighbors_and_parents(QuadTreeNode * node) {
  if (node == nullptr) return;
  node->schur_updated = false;
  node->interaction_lists.active_box.clear();
  node->interaction_lists.skel.clear();
  node->interaction_lists.skelnear.clear();
  node->interaction_lists.redundant.clear();
  node->interaction_lists.permutation.clear();

  for (QuadTreeNode* neighbor : node->neighbors) {
    neighbor->schur_updated = false;
    neighbor->interaction_lists.active_box.clear();
    neighbor->interaction_lists.skel.clear();
    neighbor->interaction_lists.skelnear.clear();
    neighbor->interaction_lists.redundant.clear();
    neighbor->interaction_lists.permutation.clear();
  }
  mark_neighbors_and_parents(node->parent);
}


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

  //TODO(John) the below needs to be changed for stokes

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
        for (unsigned int idx : node->interaction_lists.original_box) {
          std::unordered_map<int, int>::const_iterator element =
            old_index_to_new_index.find(idx);
          if (element != old_index_to_new_index.end()) {
            ob.push_back(element->second);
          }
        }
        node->interaction_lists.original_box = ob;
      }

      for (unsigned int idx : node->interaction_lists.active_box) {
        std::unordered_map<int, int>::const_iterator element =
          old_index_to_new_index.find(idx);
        if (element != old_index_to_new_index.end()) {
          ab.push_back(element->second);
        }
      }

      for (unsigned int idx : node->interaction_lists.skel) {
        std::unordered_map<int, int>::const_iterator element =
          old_index_to_new_index.find(idx);
        if (element != old_index_to_new_index.end()) {
          s.push_back(element->second);
        }
      }

      for (unsigned int idx : node->interaction_lists.redundant) {
        std::unordered_map<int, int>::const_iterator element =
          old_index_to_new_index.find(idx);
        if (element != old_index_to_new_index.end()) {
          r.push_back(element->second);
        }
      }

      for (unsigned int idx : node->interaction_lists.skelnear) {
        std::unordered_map<int, int>::const_iterator element =
          old_index_to_new_index.find(idx);
        if (element != old_index_to_new_index.end()) {
          sn.push_back(element->second);
        }
      }

      node->interaction_lists.active_box = ab;
      node->interaction_lists.skel = s;
      node->interaction_lists.skelnear = sn;
      node->interaction_lists.redundant = r;
    }
  }

  // go through all additions, find their leaves, make addition and call mark
  // function
  for (QuadTreeLevel* level : levels) {
    for (QuadTreeNode* node : level->nodes) {
      if (!node->is_leaf) {
        continue;
      }
      for (int i = 0; i < additions.size(); i++) {
        double difx = new_points[2 * additions[i]] - node->corners[0];
        double dify = new_points[2 * additions[i] + 1]
                      - node->corners[1];
        if (difx < node->side_length && dify < node->side_length && difx > 0
            && dify > 0) {
          node->interaction_lists.original_box.push_back(additions[i]);
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
      for (int i = 0; i < deletions.size(); i++) {
        double difx = old_points[2 * deletions[i]] - node->corners[0];
        double dify = old_points[2 * deletions[i] + 1]
                      - node->corners[1];
        if (difx < node->side_length && dify < node->side_length && difx > 0
            && dify > 0) {
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
  boundary->perturbation_size = perturbed_boundary.perturbation_size;
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
  initialize_tree(boundary, solution_dimension);
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
  QuadTreeNode::id_count = 0;
  initialize_tree(boundary_, solution_dimension);
}

}  // namespace ie_solver
