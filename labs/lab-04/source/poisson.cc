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
#include "poisson.h"

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/base/parameter_handler.h>



#include <fstream>
#include <iostream>

using namespace dealii;

template <int dim>
Poisson<dim>::Poisson()
  : ParameterAcceptor("/")
  , fe(1)
  , dof_handler(triangulation)
{
  // ParameterHandler prm;
  /*
  prm.declare_entry("Number of global refinements ","1",Patterns::Integer(0));
  prm.declare_entry("Rhs expression","0",Patterns::Anything());
  prm.parse_input("poisson-prm");

  global_refinements=prm.get_integer("Number of global refinements");
  rhs_expression=prm.get("Rhs expression");
  */

  // Option2 
  this-> add_parameter("Number rifiniment cycles",n_refinement_cycles);

  this-> add_parameter("Exact expression",exact_expression);

  this-> add_parameter("Rhs expression",rhs_expression);

  this-> add_parameter("Problem constants",constants);


  this->prm.enter_subsection("Table");
  table.add_parameters(this->prm);
  this->prm.leave_subsection();
  // this-> add_parameter("Exact solution",exact_solution_expression);

  //if you want to read a list of point 
  // this->add_parameter("Some points",some_points);
  // this->add_parameter("Some constants",constants);

}

template <int dim>
void
Poisson<dim>::initialize(const std::string &filename)
{
  ParameterAcceptor::initialize(filename);
}


template <int dim>
void
Poisson<dim>::parse_string(const std::string &parameters)
{
  ParameterAcceptor::prm.parse_input_from_string(parameters);
  ParameterAcceptor::parse_all_parameters();
}


template <int dim>
void
Poisson<dim>::make_grid()
{
/**
   * @brief Hyper cube generation
   */

  GridGenerator::hyper_cube(triangulation, -1, 1);
  
  /**
   * @brief First way to apply the Bc conditions 
   */
  triangulation.begin_active()->face(0)->set_boundary_id(1);

  /**
   * @brief Second way to apply the Bc conditions 
   */
  // for (auto &face : triangulation.active_face_iterators())
  //   if (std::fabs(face->center()(1) - (-1.0)) < 1e-12 ||
  //       std::fabs(face->center()(1) - (1.0)) < 1e-12)
  //     face->set_boundary_id(-0.5);

  triangulation.refine_global(5);

  // /**
  //  * @brief Hyper L generation 
  //  * 
  //  */

  // const Point<2> corner(0, 0);

  // GridGenerator::hyper_L(triangulation, -1, 1, false);

  // for (auto &face : triangulation.active_face_iterators())
  //   {
  //     for(const auto v :face->vertex_indices()){
  //         if(face->vertex(v)==corner && face->at_boundary()==true)
  //         {
  //           face->set_boundary_id(1);
  //           break;
  //         }
  //     }
  //   }

  // for (unsigned int step = 0; step < 5; ++step)
  // {

  //   for (auto &cell : triangulation.active_cell_iterators())
  //     {
  //       for (const auto v : cell->vertex_indices())
  //         {
  //           const double distance_from_center =
  //             corner.distance(cell->vertex(v));

  //           if (std::fabs(distance_from_center) <=2/3)
  //             {
  //               cell->set_refine_flag();
  //               break;
  //             }
  //         }
  //     }

  //   triangulation.execute_coarsening_and_refinement();
  // }
  // std::cout << "Number of active cells: " << triangulation.n_active_cells()
  //           << std::endl;
}


template <int dim>
void
Poisson<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);
  std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;
  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);
  system_matrix.reinit(sparsity_pattern);
  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());

  exact=std::make_unique<FunctionParser<dim>>();
  rhs=std::make_unique<FunctionParser<dim>>();
  
  std::cout<<"eccomi "<<std::endl;
  exact-> initialize(dim == 1 ? "x" :
                    dim == 2 ? "x,y" :
                                "x,y,z",
                    exact_expression,
                    constants);
  std::cout<<"eccomi2 "<<std::endl;

  rhs-> initialize(dim == 1 ? "x" :
                  dim == 2 ? "x,y" :
                              "x,y,z",
                  rhs_expression,
                  constants);


}


template <int dim>
void
Poisson<dim>::assemble_system()
{
  QGauss<dim>          quadrature_formula(fe.degree + 1);
  FEValues<dim>        fe_values(fe,
                        quadrature_formula,
                        update_values | update_quadrature_points | update_gradients | update_JxW_values);
  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      cell_matrix = 0;
      cell_rhs    = 0;
      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
          for (const unsigned int i : fe_values.dof_indices())
            for (const unsigned int j : fe_values.dof_indices())
              cell_matrix(i, j) +=
                (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                 fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                 fe_values.JxW(q_index));           // dx
          for (const unsigned int i : fe_values.dof_indices())
            cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                            rhs->value(fe_values.quadrature_point(q_index)) *  // f(x_q)
                            fe_values.JxW(q_index));            // dx
        }
      cell->get_dof_indices(local_dof_indices);
      for (const unsigned int i : fe_values.dof_indices())
        for (const unsigned int j : fe_values.dof_indices())
          system_matrix.add(local_dof_indices[i],
                            local_dof_indices[j],
                            cell_matrix(i, j));
      for (const unsigned int i : fe_values.dof_indices())
        system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }
  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           *exact,
                                           boundary_values);
  MatrixTools::apply_boundary_values(boundary_values,
                                     system_matrix,
                                     solution,
                                     system_rhs);
}


template <int dim>
void
Poisson<dim>::solve()
{
  SolverControl            solver_control(1000, 1e-12);
  SolverCG<Vector<double>> solver(solver_control);
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
}


template <int dim>
void
Poisson<dim>::output_results(const unsigned cycle) const
{
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  data_out.build_patches();
  std::string   fname = "poisson_2d_" + std::to_string(cycle) + ".vtu";
  std::ofstream output(fname);
  data_out.write_vtu(output);
}


template <int dim>
void
Poisson<dim>::run()
{
  
  make_grid();
  for(unsigned int cycle=0; cycle <n_refinement_cycles; ++cycle){
    setup_system();
    assemble_system();
    solve();
    table.error_from_exact(dof_handler, solution, *exact);
    output_results(cycle);

    if (cycle < n_refinement_cycles - 1)
      triangulation.refine_global(1);
  }

  table.output_table(std::cout);
}


template class Poisson<1>;
template class Poisson<2>;
template class Poisson<3>;
