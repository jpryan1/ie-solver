// Copyright 2019 John Paul Ryan
#include "ie_solver/io/io.h"

#include <cassert>
#include <cstring>
#include <iostream>
#include "ie_solver/boundaries/boundary.h"
#include "ie_solver/log.h"

namespace ie_solver {

void io::write_boundary_to_file(const std::string& filename,
                                const std::vector<double>& points) {
  assert(points.size() % 2 == 0);
  std::ofstream output;
  output.open(filename);
  if (output.is_open()) {
    for (unsigned int i = 0; i < points.size(); i += 2) {
      output << points[i] << "," << points[i + 1] << std::endl;
    }
    output.close();
  } else {
    LOG::ERROR("Failed to open boundary output file!");
  }
}


void io::write_potential_to_file() {
  std::cout << "unimplemented" << std::endl;
}


void io::write_times_to_files(int* scale_n,
                              const std::vector<double>& n_times,
                              double* scale_eps,
                              const std::vector<double>& eps_times) {
  std::ofstream n_output, e_output;
  n_output.open("output/data/ie_solver_n_scaling.txt");
  e_output.open("output/data/ie_solver_e_scaling.txt");

  if (n_output.is_open()) {
    for (unsigned int i = 0; i < n_times.size(); i++) {
      n_output << scale_n[i] << "," << n_times[i] << std::endl;
    }
    n_output.close();
  } else {
    printf("Failed to open n output file!\n");
  }

  if (e_output.is_open()) {
    for (unsigned int i = 0; i < eps_times.size(); i++) {
      e_output << scale_eps[i] << "," << eps_times[i] << std::endl;
    }
    e_output.close();
  } else {
    printf("Failed to open e output file!\n");
  }
}


void io::write_solution_to_file(const std::string& filename,
                                const ie_Mat& domain,
                                const std::vector<double>& domain_points,
                                int solution_dimension) {
  assert(domain.height() > 0 && domain.width() == 1);
  std::ofstream output;
  output.open(filename);

  int points_index = 0;

  if (output.is_open()) {
    for (unsigned int i = 0; i < domain.height(); i += solution_dimension) {
      output << domain_points[points_index] << "," <<
             domain_points[points_index + 1] << ",";
      points_index += 2;  // depends on domain dimension
      for (int dim = 0; dim < solution_dimension; dim++) {
        output << domain.get(i + dim, 0);
        if (dim < solution_dimension - 1) {
          output << ",";
        }
      }
      output << std::endl;
    }
    output.close();
  } else {
    printf("Failed to open solution output file!\n");
  }
}


void io::write_quadtree_to_file(const std::string& filename,
                                const QuadTree& quadtree) {
  std::ofstream output;
  output.open(filename);
  if (output.is_open()) {
    for (int lvl = 0; lvl < quadtree.levels.size(); lvl++) {
      QuadTreeLevel* level = quadtree.levels[lvl];
      output << "level " << lvl << std::endl;
      for (QuadTreeNode* node : level->nodes) {
        output << node->corners[0] << "," << node->corners[1] << ","
               << node->side_length << "," << node->compression_ratio << ", " <<
               node->compress_time << std::endl;
      }
    }
    output.close();
  } else {
    LOG::ERROR("failed to open output file!");
  }
}


void io::write_ex2_gradients_to_file(const std::string& filename,
                                     std::vector<double> angs_and_grads) {
  std::ofstream output;
  output.open(filename);
  if (output.is_open()) {
    for (int i = 0; i < angs_and_grads.size() / 3; i++) {
      output << angs_and_grads[3 * i] << " " << angs_and_grads[3 * i + 1] << " " <<
             angs_and_grads[3 * i + 2] << std::endl;
    }
    output.close();
  } else {
    LOG::ERROR("failed to open output file!");
  }
}

void io::write_ex3_flows_to_file(const std::string& filename,
                                 std::vector<double> ang_and_flow) {
  std::ofstream output;
  output.open(filename);
  if (output.is_open()) {
    for (int i = 0; i < ang_and_flow.size() / 2; i++) {
      output << ang_and_flow[2 * i] << " " << ang_and_flow[2 * i + 1] <<
             std::endl;
    }
    output.close();
  } else {
    LOG::ERROR("failed to open output file!");
  }
}

// TODO(John) this needs to be updated following removal of STOKES
int io::parse_input_into_config(int argc, char** argv,
                                ie_solver_config* config) {

  std::string usage = "\n\tusage: ./ie_solver "
                      "-pde {LAPLACE|LAPLACE_NEUMANN|STOKES} "
                      "-boundary {CIRCLE|ANNULUS|"
                      "CUBIC_SPLINE|EX1|EX2|EX3} "
                      "-boundary_condition {SINGLE_ELECTRON|ALL_ONES|"
                      "TANGENT_VEC|NORMAL_VEC|NO_SLIP} "
                      "-N {number of nodes} "
                      "-D {side length of square domain grid} "
                      "-e {ID error tolerance} "
                      "{-scaling} {-strong} {-animation}"
                      "\nOmitting an arg triggers a default value.";
  std::string boundary_name = "CIRCLE";
  std::string pde_name = "LAPLACE";
  std::string boundary_condition_name = "SINGLE_ELECTRON";
  // The default boundary is Circle, set it here.

  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-pde")) {
      if (i < argc - 1) {
        if (!strcmp(argv[i + 1], "STOKES")) {
          config->pde = ie_solver_config::STOKES;
        } else if (!strcmp(argv[i + 1], "LAPLACE")) {
          config->pde = ie_solver_config::LAPLACE;
        } else if (!strcmp(argv[i + 1], "LAPLACE_NEUMANN")) {
          config->pde = ie_solver_config::LAPLACE_NEUMANN;
        } else {
          LOG::ERROR("Unrecognized pde: " + std::string(argv[i + 1])
                     + usage);
          return -1;
        }
        pde_name = argv[i + 1];
      }
      i++;
    } else if (!strcmp(argv[i], "-N")) {
      if (i < argc - 1) {
        config->num_boundary_points = std::stoi(argv[i + 1]);
      }
      i++;
    } else if (!strcmp(argv[i], "-D")) {
      if (i < argc - 1) {
        config->domain_size = std::stoi(argv[i + 1]);
      }
      i++;
    } else if (!strcmp(argv[i], "-strong")) {
      config->is_strong_admissibility = true;
    } else if (!strcmp(argv[i], "-scaling")) {
      config->scaling = true;
    } else if (!strcmp(argv[i], "-animation")) {
      config->animation = true;
    } else if (!strcmp(argv[i], "-h")) {
      LOG::log_level_ = LOG::LOG_LEVEL::INFO_;
      LOG::INFO(usage);
      return -1;
    } else if (!strcmp(argv[i], "-e")) {
      if (i < argc - 1) {
        config->id_tol = std::stof(argv[i + 1]);
      }
      i++;
    } else if (!strcmp(argv[i], "-boundary")) {
      if (i < argc - 1) {
        if (!strcmp(argv[i + 1], "CIRCLE")) {
          config->boundary_shape = Boundary::BoundaryShape::CIRCLE;
        } else if (!strcmp(argv[i + 1], "ANNULUS")) {
          config->boundary_shape = Boundary::BoundaryShape::ANNULUS;
        } else if (!strcmp(argv[i + 1], "CUBIC_SPLINE")) {
          config->boundary_shape = Boundary::BoundaryShape::CUBIC_SPLINE;
        }else if (!strcmp(argv[i + 1], "EX1")) {
          config->boundary_shape = Boundary::BoundaryShape::EX1;
        }else if (!strcmp(argv[i + 1], "EX2")) {
          config->boundary_shape = Boundary::BoundaryShape::EX2;
        }else if (!strcmp(argv[i + 1], "EX3")) {
          config->boundary_shape = Boundary::BoundaryShape::EX3;
        } else {
          LOG::ERROR("Unrecognized boundary: " + std::string(argv[i + 1])
                     + usage);
          return -1;
        }
        boundary_name = argv[i + 1];
      }
      i++;
    } else if (!strcmp(argv[i], "-boundary_condition")) {
      if (i < argc - 1) {
        if (!strcmp(argv[i + 1], "SINGLE_ELECTRON")) {
          config->boundary_condition =
            BoundaryCondition::SINGLE_ELECTRON;
        } else if (!strcmp(argv[i + 1], "ALL_ONES")) {
          config->boundary_condition = BoundaryCondition::ALL_ONES;
        } else if (!strcmp(argv[i + 1], "TANGENT_VEC")) {
          config->boundary_condition = BoundaryCondition::TANGENT_VEC;
        } else if (!strcmp(argv[i + 1], "NORMAL_VEC")) {
          config->boundary_condition = BoundaryCondition::NORMAL_VEC;
        } else if (!strcmp(argv[i + 1], "NO_SLIP")) {
          config->boundary_condition = BoundaryCondition::NO_SLIP;
        } else {
          LOG::ERROR("Unrecognized boundary_condition: " +
                     std::string(argv[i + 1]) + usage);
          return -1;
        }
        boundary_condition_name = argv[i + 1];
      }
      i++;
    } else {
      LOG::ERROR("Unrecognized argument: " + std::string(argv[i]) + usage);
      return -1;
    }
  }

  if (config->pde == ie_solver::ie_solver_config::STOKES) {
    config->solution_dimension = 2;
  }

  if(config->solution_dimension == 2){
    if(config->boundary_condition ==
            BoundaryCondition::SINGLE_ELECTRON ||config->boundary_condition ==
            BoundaryCondition::ALL_ONES){
      LOG::ERROR("Boundary condition dimension does not match PDE solution dimention.");
      return -1;
    }
  }else if(config->solution_dimension == 1){
    if(config->boundary_condition !=
            BoundaryCondition::SINGLE_ELECTRON && config->boundary_condition !=
            BoundaryCondition::ALL_ONES){
      LOG::ERROR("Boundary condition dimension does not match PDE solution dimention.");
      return -1;
    }
  }

  LOG::INFO("PDE: " + pde_name);
  LOG::INFO("Boundary: " + boundary_name);
  LOG::INFO("Boundary condition: " + boundary_condition_name);
  LOG::INFO("Number of nodes: " + std::to_string(config->num_boundary_points));
  LOG::INFO("Side length of square domain grid: " + std::to_string(
              config->domain_size));
  LOG::INFO("ID error tolerance: " + std::to_string(config->id_tol));
  if (config->scaling) {
    LOG::INFO("Scaling run");
  }
  if (config->animation) {
    LOG::INFO("Animation run");
  }
  if (config->is_strong_admissibility) {
    LOG::INFO("Strong admissibility");
  } else {
    LOG::INFO("Weak admissibility");
  }

  return 0;
}

}  // namespace ie_solver
