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

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/base/parsed_convergence_table.h>


#include <fstream>
#include <iostream>

// Forward declare the tester class
// class PoissonTester;

using namespace dealii;

template <int dim>
class Poisson : public ParameterAcceptor
{
public:
  Poisson();
  void
  run();
  void 
  initialize(const std::string &);
  void
  parse_string(const std::string &);

protected:
  void
  make_grid();
  void
  setup_system();
  void
  assemble_system();
  void
  solve();
  void
  output_results(const unsigned ) const;


  Triangulation<dim>   triangulation;
  FE_Q<dim>            fe;
  DoFHandler<dim>      dof_handler;
  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;
  Vector<double>       solution;
  Vector<double>       system_rhs;

  //parameter
  unsigned int n_refinement_cycles=1;
  std::string exact_expression= "1";
  std::string rhs_expression= "1";

  std::map<std::string,double> constants;

  std::unique_ptr<FunctionParser<dim>> exact;
  std::unique_ptr<FunctionParser<dim>> rhs ;



  ParsedConvergenceTable table;

  /*
    std::unique_ptr<FunctionParse<2>> rhs;
    std::unique_ptr<FunctionParse<2>> exact;
    in file .cc 
    rhs=std::make_unique<FunctionParse<2>>:(rhs_expression);
  */

  template<typename Integral>
  friend class PoissonTester;
};

#endif