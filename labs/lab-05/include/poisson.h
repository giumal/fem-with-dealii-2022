/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2020 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 *
 * Authors: Wolfgang Bangerth, 1999,
 *          Guido Kanschat, 2011
 *          Luca Heltai, 2021
 */

// Make sure we don't redefine things
#ifndef poisson_include_file
#define poisson_include_file

#include <deal.II/base/function.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/parsed_convergence_table.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/grid/grid_refinement.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

// Forward declare the tester class
template <typename Integral>
class PoissonTester;

using namespace dealii;

template <int dim>
class Poisson : ParameterAcceptor
{
public:
  Poisson();
  void
  run();

  void
  initialize(const std::string &filename);

  void
  parse_string(const std::string &par);

protected:
  void
  make_grid();
  void
  refine_grid();
  void
  setup_system();
  void
  assemble_system();
  void
  solve();
  void
  output_results(const unsigned cycle) const;

  Triangulation<dim>         triangulation;
  std::unique_ptr<FE_Q<dim>> fe;
  DoFHandler<dim>            dof_handler;
  AffineConstraints<double>  constraints;
  SparsityPattern            sparsity_pattern;
  SparseMatrix<double>       system_matrix;
  Vector<double>             solution;
  Vector<double>             system_rhs;

  FunctionParser<dim> exact_term;
  FunctionParser<dim> forcing_term;
  FunctionParser<dim> dirichlet_boundary_condition;
  FunctionParser<dim> neumann_boundary_condition;


  unsigned int fe_degree           = 1;
  unsigned int n_refinements       = 4;
  unsigned int n_refinement_cycles = 1;
  std::string  output_filename     = "poisson";

  std::set<types::boundary_id> dirichlet_ids = {0};
  std::set<types::boundary_id> neumann_ids;

  std::string                   exact_expression                  = "1";
  std::string                   forcing_term_expression                  = "1";
  std::string                   dirichlet_boundary_conditions_expression = "0";
  std::string                   neumann_boundary_conditions_expression   = "0";
  std::map<std::string, double> constants;

  std::string grid_generator_function  = "hyper_cube";
  std::string grid_generator_arguments = "0: 1: false";

  ParsedConvergenceTable error_table;

  template <typename Integral>
  friend class PoissonTester;
};

#endif