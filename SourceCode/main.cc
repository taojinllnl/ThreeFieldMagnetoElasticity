/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2010 - 2022 by the deal.II authors
 *                          and Tao Jin
 *
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
*/

/* Source code for paper draft "A three-field based finite element analysis 
 * for a class of magnetoelastic materials"
 *
 * Author: Tao Jin, PhD (https://www.taojinlab.com/)
 *         Department of Mechanical Engineering
 *         University of Ottawa
 *         Ottawa, Ontario K1N 6N5, Canada
 * Starting date: Oct. 2022
 * Release date: Feb. 2023
 */

/* Reference:
 *
 * This work is based on the constitutive formulation proposed in the paper:
 * Zhao, R., Kim, Y., Chester, S.A., Sharma, P., Zhao, X., 2019. Mechanics of
 * hard-magnetic soft materials. Journal of the Mechanics and Physics of Solids
 * 124, 244–263. doi:https://doi.org/10.1016/j.jmps.2018.10.008.
 *
 * The three-field finite element implementation is developed by using Step-44
 * of the deal.II library as the basis:
 * Pelteret, Jean-Paul, & McBride, Andrew. (2012). The deal.II tutorial step-44:
 * Three-field formulation for non-linear solid mechanics. Zenodo.
 * https://doi.org/10.5281/zenodo.439772
 */


/* Main features of the code:
 * Finite deformation
 * Three-field formulation
 * Modified neo-Hookean model:
 *   total_energy = volumetric_energy + deviatoric_energy
 *   volumetric_energy =
 *                type 1: m_K/2 * pow((m_det_F - 1.0), 2);
 *                type 2: m_K/4 * [J^2 - 1.0 - 2*ln(J)];
 *   I1_bar = pow(m_det_F, -2.0/3) * trace(m_C) = pow(m_det_F, -2.0/3) * trace(m_b);
 *   deviatoric_energy = m_G/2 * (I1_bar - 3) = c1 (I1_bar - 3);
 * Adaptive mesh refinement for the first time step
 * Multi-thread for system assembly (TBB)
 */

/**********************************************************
 * WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * Due to the existence of magnetic field, the total Cauchy
 * stress and the tangent modulus might not be symmetric
 * anymore. Be careful about the type of tensors used to
 * represent the above two quantities.
 **********************************************************/

/**********************************************************
 * operator* overloaded for different type
 * of tensors in deal.ii:
 *
 * 1. for two symmetric 2nd order tensors A and B, A * B means
 *    double contraction between A and B, that is,
 *    A : B  = A_{ij}B_{ij}, which generates a scalar;
 *
 * 2. for a symmetric 2nd order tensor A and a regular 2nd
 *    order tensor B, A * B or B * A means tensor product,
 *    that is, A * B = A_{ik} B_{kj} or B * A = B_{ik} A_{kj}
 *    which generates another 2nd order tensor;
 *
 * 3. for two regular 2nd order tensors A and B, A * B or
 *    B * A means tensor product,
 *    that is, A * B = A_{ik} B_{kj} or B * A = B_{ik} A_{kj}
 *    which generates another 2nd order tensor;
 *
 * 4. for two regular 2nd order tensors A and B, in order to
 *    perform the double contraction operation, call
 *    c = scalar_product(A, B), which c is a scalar. That is,
 *    c = A_{ij} B_{ij};
 *
 * 5. for two regular 2nd order tensors A and B, another way
 *    to perform double contraction A_{ij} B_{ij} is
 *    c = double_contract<0,0,1,1>(A, B), where c is a scalar.
 *    For double contraction A_{ij} B_{ji}, we call
 *    c = double_contract<0,1,1,0>(A, B);
 *
 * 6. for a 2nd order symmetric tensor A and a 4th order
 *    symmetric tensor C, A * C = A_{ij} C_{ijkl} generates
 *    another 2nd order symmetric tensor.
 *    C * A = C_{ijkl} A_{kl} also generates a 2nd order
 *    symmetric tensor. If either A or C is a regular tensor,
 *    then "*" means tensor product, that is,
 *    A * C = A_{im} C_{mjkl},
 *    C * A = C_{ijkm} A_{ml} generates another 4th order
 *    tensor;
 *
 * 7. for a 2nd order general tensor A and a 4th order general
 *    tensor C, in order to perform the double contraction,
 *    A : C = A_{ij} C_{ijkl}, we call
 *    double_contract<0,0,1,1>(A, C) to generate a second order
 *    tensor. Similarly, for C : A = C_{ijkl} A_{kl}, we call
 *    double_contract<2,0,3,1>(C, A) to generate a second order
 *    tensor;
 *
 * 8. Whenever call double_contract<>(A,B), both A and B need to
 *    be regular tensor and CANNOT be symmetric tensor;
 *
 * 9. Whenever call scalar_product(A,B), both A and B can be either
 *    a regular or symmetric 2nd order tensor. Since the function
 *    interfaces exist for all cases;
 *
 **********************************************************/

#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/base/quadrature_point_data.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition_selector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_selector.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/packaged_operation.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/physics/elasticity/kinematics.h>
#include <deal.II/physics/elasticity/standard_tensors.h>

#include <deal.II/numerics/error_estimator.h>
#include <deal.II/grid/grid_refinement.h>

#include <iostream>
#include <fstream>


namespace FiniteDeformationMagnetoElasticityThreeField
{
  using namespace dealii;

  template <int dim>
  std::vector<types::global_dof_index> get_vertex_dofs(
    const typename Triangulation<dim>::active_vertex_iterator &vertex,
    const DoFHandler<dim> &dof_handler)
  {
    DoFAccessor<0, dim, dim, false> vertex_dofs(
        &(dof_handler.get_triangulation()),
        vertex->level(),
        vertex->index(),
        &dof_handler);
    const unsigned int n_dofs = dof_handler.get_fe().dofs_per_vertex;
    std::vector<types::global_dof_index> dofs(n_dofs);
    for (unsigned int i = 0; i < n_dofs; ++i)
    {
      dofs[i] = vertex_dofs.vertex_dof_index(0, i);
    }
    return dofs;
  }

  // body force
  template <int dim>
  void right_hand_side(const std::vector<Point<dim>> &points,
				   std::vector<Tensor<1, dim>> &  values,
				   const double fx,
				   const double fy,
				   const double fz)
  {
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));
    Assert(dim >= 2, ExcNotImplemented());

    for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
      {
	if (dim == 2)
	  {
	    values[point_n][0] = fx;
	    values[point_n][1] = fy;
	  }
	else
	  {
	    values[point_n][0] = fx;
	    values[point_n][1] = fy;
	    values[point_n][2] = fz;
	  }
      }
  }

  namespace Parameters
  {
    struct Scenario
    {
      unsigned int m_scenario;
      unsigned int m_global_refine;
      unsigned int m_adaptive_refine;
      unsigned int m_total_material_regions;
      std::string m_material_file_name;
      unsigned int m_volumetric_type;

      static void declare_parameters(ParameterHandler &prm);
      void parse_parameters(ParameterHandler &prm);
    };

    void Scenario::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Scenario");
      {
        prm.declare_entry("Scenario number",
                          "1",
                          Patterns::Integer(0),
                          "Geometry, loading and boundary conditions scenario");
        prm.declare_entry("Global refinement",
                          "2",
                          Patterns::Integer(0),
                          "Number of global refinement");
        prm.declare_entry("Adaptive refinement",
                          "0",
                          Patterns::Integer(0),
                          "Number of adaptive refinement");
        prm.declare_entry("Material regions",
                          "1",
                          Patterns::Integer(0),
                          "Number of material regions");
        prm.declare_entry("Material data file",
                          "1",
                          Patterns::FileName(Patterns::FileName::input),
                          "Material data file");
        prm.declare_entry("Material volumetric type",
                          "1",
                          Patterns::Integer(0),
                          "Type of volumetric strain energy function (1 or 2)");
      }
      prm.leave_subsection();
    }

    void Scenario::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Scenario");
      {
        m_scenario = prm.get_integer("Scenario number");
        std::cout << "scenario number = " << m_scenario << std::endl;
        m_global_refine = prm.get_integer("Global refinement");
        std::cout << "global refinement number = " << m_global_refine << std::endl;
        m_adaptive_refine = prm.get_integer("Adaptive refinement");
        std::cout << "adaptive refinement number = " << m_adaptive_refine << std::endl;
        m_total_material_regions = prm.get_integer("Material regions");
        std::cout << "total number of material types = " << m_total_material_regions << std::endl;
        m_material_file_name = prm.get("Material data file");
        std::cout << "material data file name = " << m_material_file_name << std::endl;
        m_volumetric_type = prm.get_integer("Material volumetric type");
        std::cout << "type of material volumetric function = " << m_volumetric_type << std::endl;
      }
      prm.leave_subsection();
    }

    struct FESystem
    {
      unsigned int m_poly_degree;
      unsigned int m_quad_order;

      static void declare_parameters(ParameterHandler &prm);

      void parse_parameters(ParameterHandler &prm);
    };


    void FESystem::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Finite element system");
      {
        prm.declare_entry("Polynomial degree",
                          "2",
                          Patterns::Integer(0),
                          "Displacement system polynomial order");

        prm.declare_entry("Quadrature order",
                          "3",
                          Patterns::Integer(0),
                          "Gauss quadrature order");
      }
      prm.leave_subsection();
    }

    void FESystem::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Finite element system");
      {
        m_poly_degree = prm.get_integer("Polynomial degree");
        m_quad_order  = prm.get_integer("Quadrature order");
      }
      prm.leave_subsection();
    }


    struct MaterialsMagnetism
    {
      double m_permeability;
      // B_applied in the deformed configuration
      Tensor<1, 3> m_B_applied;
      unsigned int m_magnetic_type;

      static void declare_parameters(ParameterHandler &prm);
      void parse_parameters(ParameterHandler &prm);
    };

    void MaterialsMagnetism::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Material magnetic properties");
      {
        prm.declare_entry("Material magnetic type",
                          "1",
                          Patterns::Integer(0),
                          "Type of magnetic strain energy function (1-F or 2-Fbar)");

        prm.declare_entry("Material permeability",
                          "1.2566370614e-6",
                          Patterns::Double(),
                          "Permeability mu_0 (default value 1.2566370614e-6 N/A^2)");

        prm.declare_entry("Applied magnetic field x component",
			  "1.0",
			  Patterns::Double(),
			  "B_applied x-component in the deformed configuration");

        prm.declare_entry("Applied magnetic field y component",
			  "1.0",
			  Patterns::Double(),
			  "B_applied y-component in the deformed configuration");

        prm.declare_entry("Applied magnetic field z component",
			  "1.0",
			  Patterns::Double(),
			  "B_applied z-component in the deformed configuration");
      }
      prm.leave_subsection();
    }

    void MaterialsMagnetism::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Material magnetic properties");
      {
        m_magnetic_type = prm.get_integer("Material magnetic type");
        std::cout << "type of material magnetic function (1-F; 2-Fbar) = " << m_magnetic_type << std::endl;

        m_permeability = prm.get_double("Material permeability");
        std::cout << "permeability = " << m_permeability << std::endl;

        double B_x = prm.get_double("Applied magnetic field x component");
        double B_y = prm.get_double("Applied magnetic field y component");
        double B_z = prm.get_double("Applied magnetic field z component");

        m_B_applied[0] = B_x;
        m_B_applied[1] = B_y;
        m_B_applied[2] = B_z;

        std::cout << "Applied magnetic field in deformed configuration = ("
                  << m_B_applied[0] << ", "
		  << m_B_applied[1] << ", "
		  << m_B_applied[2] << ")"
		  << std::endl;
      }
      prm.leave_subsection();
    }

    struct LinearSolver
    {
      std::string m_type_lin;
      double      m_tol_lin;
      double      m_max_iterations_lin;
      bool        m_use_static_condensation;
      std::string m_preconditioner_type;
      double      m_preconditioner_relaxation;

      static void declare_parameters(ParameterHandler &prm);

      void parse_parameters(ParameterHandler &prm);
    };

    void LinearSolver::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Linear solver");
      {
        prm.declare_entry("Solver type",
                          "CG",
                          Patterns::Selection("CG|Direct"),
                          "Type of solver used to solve the linear system");

        prm.declare_entry("Residual",
                          "1e-6",
                          Patterns::Double(0.0),
                          "Linear solver residual (scaled by residual norm)");

        prm.declare_entry(
          "Max iteration multiplier",
          "1",
          Patterns::Double(0.0),
          "Linear solver iterations (multiples of the system matrix size)");

        prm.declare_entry("Use static condensation",
                          "true",
                          Patterns::Bool(),
                          "Solve the full block system or a reduced problem");

        prm.declare_entry("Preconditioner type",
                          "ssor",
                          Patterns::Selection("jacobi|ssor"),
                          "Type of preconditioner");

        prm.declare_entry("Preconditioner relaxation",
                          "0.65",
                          Patterns::Double(0.0),
                          "Preconditioner relaxation value");
      }
      prm.leave_subsection();
    }

    void LinearSolver::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Linear solver");
      {
        m_type_lin                  = prm.get("Solver type");
        m_tol_lin                   = prm.get_double("Residual");
        m_max_iterations_lin        = prm.get_double("Max iteration multiplier");
        m_use_static_condensation   = prm.get_bool("Use static condensation");
        m_preconditioner_type       = prm.get("Preconditioner type");
        m_preconditioner_relaxation = prm.get_double("Preconditioner relaxation");
      }
      prm.leave_subsection();
    }


    struct NonlinearSolver
    {
      unsigned int m_max_iterations_NR;
      double       m_tol_f;
      double       m_tol_u;

      static void declare_parameters(ParameterHandler &prm);

      void parse_parameters(ParameterHandler &prm);
    };

    void NonlinearSolver::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Nonlinear solver");
      {
        prm.declare_entry("Max iterations Newton-Raphson",
                          "10",
                          Patterns::Integer(0),
                          "Number of Newton-Raphson iterations allowed");

        prm.declare_entry("Tolerance force",
                          "1.0e-9",
                          Patterns::Double(0.0),
                          "Force residual tolerance");

        prm.declare_entry("Tolerance displacement",
                          "1.0e-6",
                          Patterns::Double(0.0),
                          "Displacement error tolerance");
      }
      prm.leave_subsection();
    }

    void NonlinearSolver::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Nonlinear solver");
      {
        m_max_iterations_NR = prm.get_integer("Max iterations Newton-Raphson");
        m_tol_f             = prm.get_double("Tolerance force");
        m_tol_u             = prm.get_double("Tolerance displacement");
      }
      prm.leave_subsection();
    }

    struct TimeInfo
    {
      double m_end_time;
      std::string m_time_file_name;

      static void declare_parameters(ParameterHandler &prm);

      void parse_parameters(ParameterHandler &prm);
    };

    void TimeInfo::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Time");
      {
        prm.declare_entry("End time", "1", Patterns::Double(), "End time");

        prm.declare_entry("Time data file",
                          "1",
                          Patterns::FileName(Patterns::FileName::input),
                          "Time data file");
      }
      prm.leave_subsection();
    }

    void TimeInfo::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Time");
      {
        m_end_time = prm.get_double("End time");
        std::cout << "End time = " << m_end_time << std::endl;

        m_time_file_name = prm.get("Time data file");
        std::cout << "Time data file name = " << m_time_file_name << std::endl;
      }
      prm.leave_subsection();
    }


    // body force in the reference configuration (N/m^3)
    struct BodyForce
    {
      double m_x_component;
      double m_y_component;
      double m_z_component;

      static void declare_parameters(ParameterHandler &prm);

      void parse_parameters(ParameterHandler &prm);
    };

    void BodyForce::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Body force in the reference configuration");
      {
        prm.declare_entry("Body force x component",
			  "0.0",
			  Patterns::Double(),
			  "Body force x-component in the reference configuration (N/m^3)");

        prm.declare_entry("Body force y component",
			  "0.0",
			  Patterns::Double(),
			  "Body force y-component in the reference configuration (N/m^3)");

        prm.declare_entry("Body force z component",
			  "0.0",
			  Patterns::Double(),
			  "Body force z-component in the reference configuration (N/m^3)");
      }
      prm.leave_subsection();
    }

    void BodyForce::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Body force in the reference configuration");
      {
        m_x_component = prm.get_double("Body force x component");
        m_y_component = prm.get_double("Body force y component");
        m_z_component = prm.get_double("Body force z component");

        std::cout << "Body force = (" << m_x_component << ", "
                                      << m_y_component << ", "
	                              << m_z_component << ") (N/m^3)"
				      << std::endl;
      }
      prm.leave_subsection();
    }


    struct AllParameters : public Scenario,
	                   public FESystem,
			   public MaterialsMagnetism,
                           public LinearSolver,
                           public NonlinearSolver,
                           public TimeInfo,
			   public BodyForce

    {
      AllParameters(const std::string &input_file);

      static void declare_parameters(ParameterHandler &prm);

      void parse_parameters(ParameterHandler &prm);
    };

    AllParameters::AllParameters(const std::string &input_file)
    {
      ParameterHandler prm;
      declare_parameters(prm);
      prm.parse_input(input_file);
      parse_parameters(prm);
    }

    void AllParameters::declare_parameters(ParameterHandler &prm)
    {
      Scenario::declare_parameters(prm);
      FESystem::declare_parameters(prm);
      MaterialsMagnetism::declare_parameters(prm);
      LinearSolver::declare_parameters(prm);
      NonlinearSolver::declare_parameters(prm);
      TimeInfo::declare_parameters(prm);
      BodyForce::declare_parameters(prm);
    }

    void AllParameters::parse_parameters(ParameterHandler &prm)
    {
      Scenario::parse_parameters(prm);
      FESystem::parse_parameters(prm);
      MaterialsMagnetism::parse_parameters(prm);
      LinearSolver::parse_parameters(prm);
      NonlinearSolver::parse_parameters(prm);
      TimeInfo::parse_parameters(prm);
      BodyForce::parse_parameters(prm);
    }
  } // namespace Parameters


  class Time
  {
  public:
    Time(const double time_end)
      : m_timestep(0)
      , m_time_current(0.0)
      , m_time_end(time_end)
      , m_delta_t(0.0)
    {}

    virtual ~Time() = default;

    double current() const
    {
      return m_time_current;
    }
    double end() const
    {
      return m_time_end;
    }
    double get_delta_t() const
    {
      return m_delta_t;
    }
    unsigned int get_timestep() const
    {
      return m_timestep;
    }
    void increment(std::vector<std::array<double, 3>> time_table)
    {
      double t_1, t_delta;
      for (auto & time_group : time_table)
        {
	  t_1 = time_group[1];
	  t_delta = time_group[2];

	  if (m_time_current < t_1 - 1.0e-6*t_delta)
	    {
	      m_delta_t = t_delta;
	      break;
	    }
        }

      m_time_current += m_delta_t;
      ++m_timestep;
    }

  private:
    unsigned int m_timestep;
    double       m_time_current;
    const double m_time_end;
    double m_delta_t;
  };


  template <int dim>
  class MaterialModifiedNeoHookeanMagnetoThreeField
  {
  public:
    MaterialModifiedNeoHookeanMagnetoThreeField(const double shear_modulus,
					        const double bulk_modulus,
						const double vacuum_permeability,
						Tensor<1, 3> Br_tilde,
						Tensor<1, 3> B_applied,
						unsigned int volumetric_type,
						unsigned int magnetic_type)
      : m_kappa(bulk_modulus)
      , m_c_1(shear_modulus / 2.0)
      , m_mu0(vacuum_permeability)
      , m_Br_tilde(Br_tilde)
      , m_B_applied(B_applied)
      , m_volumetric_type(volumetric_type)
      , m_magnetic_type(magnetic_type)
      , m_det_F(1.0)
      , m_p_tilde(0.0)
      , m_J_tilde(1.0)
      , m_b_bar(Physics::Elasticity::StandardTensors<dim>::I)
      , m_F(Physics::Elasticity::StandardTensors<dim>::I)
    {
      //const double poisson_ratio =
      //  (3*bulk_modulus-2*shear_modulus)/(2*(3*bulk_modulus+shear_modulus));
      Assert( ((3*bulk_modulus-2*shear_modulus)/(2*(3*bulk_modulus+shear_modulus)) <= 0.5)
	     &((3*bulk_modulus-2*shear_modulus)/(2*(3*bulk_modulus+shear_modulus)) >=-1.0),
	       ExcInternalError());
    }

    void update_material_data(const Tensor<2, dim> &F,
                              const double          p_tilde_in,
                              const double          J_tilde_in)
    {
      m_det_F                      = determinant(F);
      m_F                          = F;
      const Tensor<2, dim> F_bar = Physics::Elasticity::Kinematics::F_iso(F);
      m_b_bar                      = Physics::Elasticity::Kinematics::b(F_bar);
      m_p_tilde                    = p_tilde_in;
      m_J_tilde                    = J_tilde_in;

      Assert(m_det_F > 0, ExcInternalError());
    }

    double get_volumetric_strain_energy() const
    {
      double volumetric_energy = 0.0;
      if (m_volumetric_type == 1)
	{
	  volumetric_energy = m_kappa/2.0 * std::pow(m_J_tilde - 1.0, 2);
	}
      else if (m_volumetric_type == 2)
	{
	  volumetric_energy = m_kappa/4.0 * (m_J_tilde*m_J_tilde - 1.0 - 2.0*std::log(m_J_tilde));
	}
      else
	{
	  Assert(false, ExcMessage("Volumetric function type is not implemented!"));
	}
      return volumetric_energy;
    }

    double get_isotropic_strain_energy() const
    {
      const double I1_bar = trace(m_b_bar);

      double isotropic_energy;
      isotropic_energy = m_c_1 * (I1_bar - 3);

      return isotropic_energy;
    }

    double get_mechanical_strain_energy() const
    {
      return get_volumetric_strain_energy() + get_isotropic_strain_energy();
    }


    double get_magnetic_strain_energy() const
    {
      double magnetic_strain_energy = 0.0;
      if (m_magnetic_type == 1)
	{
	  magnetic_strain_energy = -1.0/m_mu0 * (m_F * m_Br_tilde) * m_B_applied;
	}
      else if (m_magnetic_type == 2)
	{
	  Tensor<2, dim> Fbar = std::pow(m_det_F, -1.0/3) * m_F;
	  magnetic_strain_energy = -1.0/m_mu0 * (Fbar * m_Br_tilde) * m_B_applied;
	}
      else
	{
	  Assert(false, ExcMessage("Magnetic function type is not implemented!"));
	}
      return magnetic_strain_energy;
    }

    double get_total_strain_energy() const
    {
      return get_mechanical_strain_energy() + get_magnetic_strain_energy();
    }

    SymmetricTensor<2, dim> get_tau()
    {
      return get_tau_iso() + get_tau_vol();
    }

    // get the Kirchoff stress tau = J * Cauchy_stress from the magnetic
    // strain energy
    Tensor<2, dim> get_tau_magnetic()
    {
      Tensor<2, dim> tau_magnetic;
      if (m_magnetic_type == 1)
	{
	  tau_magnetic = -1.0/m_mu0 * outer_product(m_B_applied, m_F * m_Br_tilde);
	}
      else if (m_magnetic_type == 2)
	{
	  Tensor<2, dim> Fbar = std::pow(m_det_F, -1.0/3) * m_F;
	  double coeff = 1.0/3/m_mu0 * (Fbar * m_Br_tilde) * m_B_applied;
	  Tensor<2, dim> temp1 = Physics::Elasticity::StandardTensors<dim>::I;
	  tau_magnetic = coeff * temp1 - 1.0/m_mu0 * outer_product(m_B_applied, Fbar * m_Br_tilde);
	}
      else
	{
	  Assert(false, ExcMessage("Magnetic function type is not implemented!"));
	}

      return tau_magnetic;
    }

    Tensor<4, dim> get_Jc_magnetic()
    {
      Tensor<4, dim> Jc_magnetic;

      if (m_magnetic_type == 1)
	{
	  Tensor<4, dim> A;
          A.clear(); // may not be necessary
          Jc_magnetic = A;
	}
      else if (m_magnetic_type == 2)
	{
	  Tensor<2, dim> Fbar = std::pow(m_det_F, -1.0/3) * m_F;
          Tensor<2, dim> temp1 = 1.0/3/m_mu0 * outer_product(m_B_applied, Fbar*m_Br_tilde);
          Tensor<2, dim> temp2 = Physics::Elasticity::StandardTensors<dim>::I;

          Tensor<4, dim> term2 = outer_product(temp1, temp2);
          Tensor<4, dim> term1_2 = outer_product(temp2, temp1);

          double coeff = -1.0/9/m_mu0 * (Fbar * m_Br_tilde) * m_B_applied;
          Tensor<4, dim> term1_1 = coeff * outer_product(temp2, temp2);

          double coeff1 = -1.0/3/m_mu0 * (Fbar * m_Br_tilde) * m_B_applied;
          // Second order identity tensor, \delta
          const SymmetricTensor<2,dim> id = Physics::Elasticity::StandardTensors<dim>::I;
          Tensor<4,dim> ITranspose; // \mathcal{\bar{I}}
          for (unsigned int a=0; a<dim; ++a)
            for (unsigned int b=0; b<dim; ++b)
              for (unsigned int c=0; c<dim; ++c)
                for (unsigned int d=0; d<dim; ++d)
                  ITranspose[a][b][c][d] = id[a][d] * id[b][c];

          Tensor<4, dim> term1_3 = coeff1 * ITranspose;

          Jc_magnetic = term2 + term1_1 + term1_2 + term1_3;
	}
      else
	{
	  Assert(false, ExcMessage("Magnetic function type is not implemented!"));
	}

      return Jc_magnetic;
    }

    SymmetricTensor<4, dim> get_Jc() const
    {
      return get_Jc_vol() + get_Jc_iso();
    }

    double get_dPsi_vol_dJ() const
    {
      if (m_volumetric_type == 1)
	{
          return m_kappa * (m_J_tilde - 1.0);
	}
      else if (m_volumetric_type == 2)
	{
          return m_kappa / 2.0 * (m_J_tilde - 1.0 / m_J_tilde);
	}
      else
	{
	  Assert(false, ExcMessage("Volumetric function type is not implemented!"));
	  return 0.0;
	}
    }

    double get_d2Psi_vol_dJ2() const
    {
      if (m_volumetric_type == 1)
	{
          return m_kappa * 1.0;
	}
      else if (m_volumetric_type == 2)
	{
          return m_kappa / 2.0 * ( 1.0 + 1.0/(m_J_tilde * m_J_tilde) );
	}
      else
	{
	  Assert(false, ExcMessage("Volumetric function type is not implemented!"));
	  return 0.0;
	}
    }

    double get_det_F() const
    {
      return m_det_F;
    }

    double get_p_tilde() const
    {
      return m_p_tilde;
    }

    double get_J_tilde() const
    {
      return m_J_tilde;
    }

  protected:
    const double m_kappa;
    const double m_c_1;
    const double m_mu0; //vacuum permeability (default value 1.2566370614e-6 N/A^2 Newton per square Ampere)
    Tensor<1, dim> m_Br_tilde;
    Tensor<1, dim> m_B_applied;
    unsigned int m_volumetric_type;
    unsigned int m_magnetic_type;

    double                  m_det_F;
    double                  m_p_tilde;
    double                  m_J_tilde;
    SymmetricTensor<2, dim> m_b_bar;
    Tensor<2, dim>          m_F;

    SymmetricTensor<2, dim> get_tau_vol() const
    {
      return m_p_tilde * m_det_F * Physics::Elasticity::StandardTensors<dim>::I;
    }

    SymmetricTensor<2, dim> get_tau_iso() const
    {
      return Physics::Elasticity::StandardTensors<dim>::dev_P * get_tau_bar();
    }

    SymmetricTensor<2, dim> get_tau_bar() const
    {
      return 2.0 * m_c_1 * m_b_bar;
    }

    SymmetricTensor<4, dim> get_Jc_vol() const
    {
      return m_p_tilde * m_det_F *
             (Physics::Elasticity::StandardTensors<dim>::IxI -
              (2.0 * Physics::Elasticity::StandardTensors<dim>::S));
    }

    SymmetricTensor<4, dim> get_Jc_iso() const
    {
      const SymmetricTensor<2, dim> tau_bar = get_tau_bar();
      const SymmetricTensor<2, dim> tau_iso = get_tau_iso();
      const SymmetricTensor<4, dim> tau_iso_x_I =
        outer_product(tau_iso, Physics::Elasticity::StandardTensors<dim>::I);
      const SymmetricTensor<4, dim> I_x_tau_iso =
        outer_product(Physics::Elasticity::StandardTensors<dim>::I, tau_iso);
      const SymmetricTensor<4, dim> c_bar = get_c_bar();

      return (2.0 / dim) * trace(tau_bar) *
               Physics::Elasticity::StandardTensors<dim>::dev_P -
             (2.0 / dim) * (tau_iso_x_I + I_x_tau_iso) +
             Physics::Elasticity::StandardTensors<dim>::dev_P * c_bar *
               Physics::Elasticity::StandardTensors<dim>::dev_P;
    }

    SymmetricTensor<4, dim> get_c_bar() const
    {
      return SymmetricTensor<4, dim>();
    }
  };


  template <int dim>
  class PointHistory
  {
  public:
    PointHistory()
      : m_F_inv(Physics::Elasticity::StandardTensors<dim>::I)
      , m_tau(SymmetricTensor<2, dim>())
      , m_tau_magnetic(Tensor<2, dim>())
      , m_d2Psi_vol_dJ2(0.0)
      , m_dPsi_vol_dJ(0.0)
      , m_Jc(SymmetricTensor<4, dim>())
      , m_Jc_magnetic(Tensor<4, dim>())
    {}

    virtual ~PointHistory() = default;

    void setup_lqp(const Parameters::AllParameters &parameters,
		   const double shear_modulus,
		   const double bulk_modulus,
		   Tensor<1, 3> Br_tilde)
    {
      m_material =
        std::make_shared<MaterialModifiedNeoHookeanMagnetoThreeField<dim>>(shear_modulus,
            	                                                           bulk_modulus,
            								   parameters.m_permeability,
            								   Br_tilde,
            								   parameters.m_B_applied,
									   parameters.m_volumetric_type,
									   parameters.m_magnetic_type);
      update_values(Tensor<2, dim>(), 0.0, 1.0);
    }

    void update_values(const Tensor<2, dim> &Grad_u_n,
                       const double          p_tilde,
                       const double          J_tilde)
    {
      const Tensor<2, dim> F = Physics::Elasticity::Kinematics::F(Grad_u_n);
      m_material->update_material_data(F, p_tilde, J_tilde);

      m_F_inv         = invert(F);
      m_tau           = m_material->get_tau();
      m_tau_magnetic  = m_material->get_tau_magnetic();
      m_Jc            = m_material->get_Jc();
      m_dPsi_vol_dJ   = m_material->get_dPsi_vol_dJ();
      m_d2Psi_vol_dJ2 = m_material->get_d2Psi_vol_dJ2();
      m_Jc_magnetic   = m_material->get_Jc_magnetic();
    }

    double get_J_tilde() const
    {
      return m_material->get_J_tilde();
    }

    double get_det_F() const
    {
      return m_material->get_det_F();
    }

    const Tensor<2, dim> &get_F_inv() const
    {
      return m_F_inv;
    }

    double get_p_tilde() const
    {
      return m_material->get_p_tilde();
    }

    const SymmetricTensor<2, dim> &get_tau() const
    {
      return m_tau;
    }

    const Tensor<2, dim> &get_tau_magnetic() const
    {
      return m_tau_magnetic;
    }

    double get_dPsi_vol_dJ() const
    {
      return m_dPsi_vol_dJ;
    }

    double get_d2Psi_vol_dJ2() const
    {
      return m_d2Psi_vol_dJ2;
    }

    const SymmetricTensor<4, dim> &get_Jc() const
    {
      return m_Jc;
    }

    const Tensor<4, dim> &get_Jc_magnetic() const
    {
      return m_Jc_magnetic;
    }

    double get_volumetric_strain_energy() const
    {
      return m_material->get_volumetric_strain_energy();
    }

    double get_isotropic_strain_energy() const
    {
      return m_material->get_isotropic_strain_energy();
    }

    double get_mechanical_strain_energy() const
    {
      return m_material->get_mechanical_strain_energy();
    }

    double get_magnetic_strain_energy() const
    {
      return m_material->get_magnetic_strain_energy();
    }

    double get_total_strain_energy() const
    {
      return m_material->get_total_strain_energy();
    }

  private:
    std::shared_ptr<MaterialModifiedNeoHookeanMagnetoThreeField<dim>> m_material;

    Tensor<2, dim> m_F_inv;

    SymmetricTensor<2, dim> m_tau;
    Tensor<2, dim>          m_tau_magnetic;
    double                  m_d2Psi_vol_dJ2;
    double                  m_dPsi_vol_dJ;

    SymmetricTensor<4, dim> m_Jc;
    Tensor<4, dim> m_Jc_magnetic;
  };



  template <int dim>
  class Solid
  {
  public:
    Solid(const std::string &input_file);

    void run();

  private:
    struct PerTaskData_ASM;
    struct ScratchData_ASM;

    struct PerTaskData_SC;
    struct ScratchData_SC;

    struct PerTaskData_UQPH;
    struct ScratchData_UQPH;

    void make_grid();
    void make_grid_case_1();
    void make_grid_case_2();
    void make_grid_case_3();
    void make_grid_case_4();
    void make_grid_case_5();
    void make_grid_case_6();
    void make_grid_case_7();
    void create_2d_grid(const std::vector<Point<2>> & vertices,
			const std::vector<std::array<unsigned int,
			                  GeometryInfo<2>::vertices_per_cell>> & vertex_indices,
			Triangulation<2> & coarse_grid,
			const unsigned int material_id);

    void refine_grid();

    void system_setup();

    void determine_component_extractors();

    void make_constraints(const unsigned int it_nr);

    void assemble_system();

    void assemble_system_one_cell(
      const typename DoFHandler<dim>::active_cell_iterator &cell,
      ScratchData_ASM &                                     scratch,
      PerTaskData_ASM &                                     data) const;

    void assemble_sc();

    void assemble_sc_one_cell(
      const typename DoFHandler<dim>::active_cell_iterator &cell,
      ScratchData_SC &                                      scratch,
      PerTaskData_SC &                                      data);

    void copy_local_to_global_sc(const PerTaskData_SC &data);

    void setup_qph();

    void update_qph_incremental(const BlockVector<double> &solution_delta);

    void update_qph_incremental_one_cell(
      const typename DoFHandler<dim>::active_cell_iterator &cell,
      ScratchData_UQPH &                                    scratch,
      PerTaskData_UQPH &                                    data);

    void copy_local_to_global_UQPH(const PerTaskData_UQPH & /*data*/)
    {}

    void solve_nonlinear_timestep(BlockVector<double> &solution_delta);

    std::pair<unsigned int, double>
    solve_linear_system(BlockVector<double> &newton_update);

    BlockVector<double>
    get_total_solution(const BlockVector<double> &solution_delta) const;

    void output_results() const;

    void read_time_data(const std::string &data_file,
    		        std::vector<std::array<double, 3>> & time_table) const;

    // Should not make this function const
    void read_material_data(const std::string &data_file,
			    const unsigned int total_material_regions);

    Parameters::AllParameters m_parameters;

    double m_vol_reference;

    Triangulation<dim> m_triangulation;

    Time                m_time;
    mutable TimerOutput m_timer;

    CellDataStorage<typename Triangulation<dim>::cell_iterator,
                    PointHistory<dim>>
      m_quadrature_point_history;

    const unsigned int               m_degree;
    const FESystem<dim>              m_fe;
    DoFHandler<dim>                  m_dof_handler;
    const unsigned int               m_dofs_per_cell;
    const FEValuesExtractors::Vector m_u_fe;
    const FEValuesExtractors::Scalar m_p_fe;
    const FEValuesExtractors::Scalar m_J_fe;

    static const unsigned int m_n_blocks          = 3;
    static const unsigned int m_n_components      = dim + 2;
    static const unsigned int m_first_u_component = 0;
    static const unsigned int m_p_component       = dim;
    static const unsigned int m_J_component       = dim + 1;

    enum
    {
      m_u_dof = 0,
      m_p_dof = 1,
      m_J_dof = 2
    };

    std::vector<types::global_dof_index> m_dofs_per_block;
    std::vector<types::global_dof_index> m_element_indices_u;
    std::vector<types::global_dof_index> m_element_indices_p;
    std::vector<types::global_dof_index> m_element_indices_J;

    const QGauss<dim>     m_qf_cell;
    const QGauss<dim - 1> m_qf_face;
    const unsigned int    m_n_q_points;
    const unsigned int    m_n_q_points_f;

    AffineConstraints<double> m_constraints;
    BlockSparsityPattern      m_sparsity_pattern;
    BlockSparseMatrix<double> m_tangent_matrix;
    BlockVector<double>       m_system_rhs;
    BlockVector<double>       m_solution_n;
    std::map<unsigned int, std::vector<double>> m_material_data;

    struct Errors
    {
      Errors()
        : m_norm(1.0)
        , m_u(1.0)
        , m_p(1.0)
        , m_J(1.0)
      {}

      void reset()
      {
        m_norm = 1.0;
        m_u    = 1.0;
        m_p    = 1.0;
        m_J    = 1.0;
      }
      void normalize(const Errors &rhs)
      {
        if (rhs.m_norm != 0.0)
          m_norm /= rhs.m_norm;
        if (rhs.m_u != 0.0)
          m_u /= rhs.m_u;
        if (rhs.m_p != 0.0)
          m_p /= rhs.m_p;
        if (rhs.m_J != 0.0)
          m_J /= rhs.m_J;
      }

      double m_norm, m_u, m_p, m_J;
    };

    Errors m_error_residual, m_error_residual_0, m_error_residual_norm, m_error_update,
      m_error_update_0, m_error_update_norm;

    void get_error_residual(Errors &error_residual);

    void get_error_update(const BlockVector<double> &newton_update,
                          Errors &                   error_update);

    std::pair<double, double> get_error_dilation() const;

    double compute_vol_current() const;

    struct Energy
    {
      Energy()
      : m_total_strain_energy(0.0)
      , m_volumetric_strain_energy(0.0)
      , m_isotropic_strain_energy(0.0)
      , m_mechanical_strain_energy(0.0)
      , m_magnetic_strain_energy(0.0)
      {}

      double m_total_strain_energy, m_volumetric_strain_energy, m_isotropic_strain_energy,
      m_mechanical_strain_energy, m_magnetic_strain_energy;
    };

    Energy m_energy;

    void compute_energy_current(Energy &energy) const;

    static void print_conv_header();

    void print_conv_footer();
  };

  template <int dim>
  Solid<dim>::Solid(const std::string &input_file)
    : m_parameters(input_file)
    , m_vol_reference(0.)
    , m_triangulation(Triangulation<dim>::maximum_smoothing)
    , m_time(m_parameters.m_end_time)
    , m_timer(std::cout, TimerOutput::summary, TimerOutput::wall_times)
    , m_degree(m_parameters.m_poly_degree)
    , m_fe(FE_Q<dim>(m_parameters.m_poly_degree),
           dim, // displacement
           FE_DGP<dim>(m_parameters.m_poly_degree - 1),
           1, // pressure
           FE_DGP<dim>(m_parameters.m_poly_degree - 1),
           1) // dilatation
    , m_dof_handler(m_triangulation)
    , m_dofs_per_cell(m_fe.n_dofs_per_cell())
    , m_u_fe(m_first_u_component)
    , m_p_fe(m_p_component)
    , m_J_fe(m_J_component)
    , m_dofs_per_block(m_n_blocks)
    , m_qf_cell(m_parameters.m_quad_order)
    , m_qf_face(m_parameters.m_quad_order)
    , m_n_q_points(m_qf_cell.size())
    , m_n_q_points_f(m_qf_face.size())
  {
    Assert(dim == 2 || dim == 3,
           ExcMessage("This problem only works in 2 or 3 space dimensions."));
    determine_component_extractors();
  }


  template <int dim>
  void Solid<dim>::run()
  {
    if (m_parameters.m_use_static_condensation && (m_parameters.m_adaptive_refine > 0))
      {
	std::cout << "If static condensation is used, then adaptive refinement cannot be used!" << std::endl;
	std::exit(0);
      }

    read_material_data(m_parameters.m_material_file_name,
		       m_parameters.m_total_material_regions);

    std::vector<std::array<double, 3>> time_table;

    read_time_data(m_parameters.m_time_file_name, time_table);

    make_grid();
    system_setup();
    {
      AffineConstraints<double> constraints;
      constraints.close();

      const ComponentSelectFunction<dim> J_mask(m_J_component, m_n_components);

      VectorTools::project(
        m_dof_handler, constraints, QGauss<dim>(m_degree + 2), J_mask, m_solution_n);
    }
    output_results();
    m_time.increment(time_table);


    // The first time step, we adaptively refine the mesh
    if (m_time.get_timestep() == 1)
      {
	const unsigned int n_cycle = m_parameters.m_adaptive_refine;
	for (unsigned int cycle = 0; cycle <= n_cycle; ++cycle)
	  {
	    if (cycle > 0)
	      {
		refine_grid();
		system_setup();
		{
		  AffineConstraints<double> constraints;
		  constraints.close();

		  const ComponentSelectFunction<dim> J_mask(m_J_component, m_n_components);

		  VectorTools::project(
		    m_dof_handler, constraints, QGauss<dim>(m_degree + 2), J_mask, m_solution_n);
		}
	      }
	    BlockVector<double> solution_delta(m_dofs_per_block);
	    solution_delta = 0.0;
	    solve_nonlinear_timestep(solution_delta);
	    m_solution_n += solution_delta;
	  }
	output_results();
	m_time.increment(time_table);
      }

    BlockVector<double> solution_delta(m_dofs_per_block);
    while(m_time.current() < m_time.end() + m_time.get_delta_t()*1.0e-6)
      {
	solution_delta = 0.0;
	solve_nonlinear_timestep(solution_delta);
	m_solution_n += solution_delta;
	output_results();
	m_time.increment(time_table);
      }
  }

  template <int dim>
  struct Solid<dim>::PerTaskData_ASM
  {
    FullMatrix<double>                   m_cell_matrix;
    Vector<double>                       m_cell_rhs;
    std::vector<types::global_dof_index> m_local_dof_indices;

    PerTaskData_ASM(const unsigned int dofs_per_cell)
      : m_cell_matrix(dofs_per_cell, dofs_per_cell)
      , m_cell_rhs(dofs_per_cell)
      , m_local_dof_indices(dofs_per_cell)
    {}

    void reset()
    {
      m_cell_matrix = 0.0;
      m_cell_rhs    = 0.0;
    }
  };


  template <int dim>
  struct Solid<dim>::ScratchData_ASM
  {
    FEValues<dim>     m_fe_values;
    FEFaceValues<dim> m_fe_face_values;

    std::vector<std::vector<double>>                  m_Nx;  // shape function values for J_tilde and p_tilde
    std::vector<std::vector<Tensor<1, dim>>>          m_Nx_disp; // shape function values for displacement field
    std::vector<std::vector<Tensor<2, dim>>>          m_grad_Nx;
    std::vector<std::vector<SymmetricTensor<2, dim>>> m_symm_grad_Nx;

    ScratchData_ASM(const FiniteElement<dim> &fe_cell,
                    const QGauss<dim> &       qf_cell,
                    const UpdateFlags         uf_cell,
                    const QGauss<dim - 1> &   qf_face,
                    const UpdateFlags         uf_face)
      : m_fe_values(fe_cell, qf_cell, uf_cell)
      , m_fe_face_values(fe_cell, qf_face, uf_face)
      , m_Nx(qf_cell.size(), std::vector<double>(fe_cell.n_dofs_per_cell()))
      , m_Nx_disp(qf_cell.size(), std::vector<Tensor<1, dim>>(fe_cell.n_dofs_per_cell()))
      , m_grad_Nx(qf_cell.size(),
                std::vector<Tensor<2, dim>>(fe_cell.n_dofs_per_cell()))
      , m_symm_grad_Nx(qf_cell.size(),
                     std::vector<SymmetricTensor<2, dim>>(
                       fe_cell.n_dofs_per_cell()))
    {}

    ScratchData_ASM(const ScratchData_ASM &rhs)
      : m_fe_values(rhs.m_fe_values.get_fe(),
                  rhs.m_fe_values.get_quadrature(),
                  rhs.m_fe_values.get_update_flags())
      , m_fe_face_values(rhs.m_fe_face_values.get_fe(),
                       rhs.m_fe_face_values.get_quadrature(),
                       rhs.m_fe_face_values.get_update_flags())
      , m_Nx(rhs.m_Nx)
      , m_Nx_disp(rhs.m_Nx_disp)
      , m_grad_Nx(rhs.m_grad_Nx)
      , m_symm_grad_Nx(rhs.m_symm_grad_Nx)
    {}

    void reset()
    {
      const unsigned int n_q_points      = m_Nx.size();
      const unsigned int n_dofs_per_cell = m_Nx[0].size();
      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        {
          Assert(m_Nx[q_point].size() == n_dofs_per_cell, ExcInternalError());
          Assert(m_grad_Nx[q_point].size() == n_dofs_per_cell,
                 ExcInternalError());
          Assert(m_symm_grad_Nx[q_point].size() == n_dofs_per_cell,
                 ExcInternalError());
          for (unsigned int k = 0; k < n_dofs_per_cell; ++k)
            {
              m_Nx[q_point][k]           = 0.0;
	      m_Nx_disp[q_point][k]      = 0.0;
              m_grad_Nx[q_point][k]      = 0.0;
              m_symm_grad_Nx[q_point][k] = 0.0;
            }
        }
    }
  };


  template <int dim>
  struct Solid<dim>::PerTaskData_SC
  {
    FullMatrix<double>                   m_cell_matrix;
    std::vector<types::global_dof_index> m_local_dof_indices;

    FullMatrix<double> m_k_orig;
    FullMatrix<double> m_k_pu;
    FullMatrix<double> m_k_pJ;
    FullMatrix<double> m_k_JJ;
    FullMatrix<double> m_k_pJ_inv;
    FullMatrix<double> m_k_bbar;
    FullMatrix<double> m_A;
    FullMatrix<double> m_B;
    FullMatrix<double> m_C;

    PerTaskData_SC(const unsigned int dofs_per_cell,
                   const unsigned int n_u,
                   const unsigned int n_p,
                   const unsigned int n_J)
      : m_cell_matrix(dofs_per_cell, dofs_per_cell)
      , m_local_dof_indices(dofs_per_cell)
      , m_k_orig(dofs_per_cell, dofs_per_cell)
      , m_k_pu(n_p, n_u)
      , m_k_pJ(n_p, n_J)
      , m_k_JJ(n_J, n_J)
      , m_k_pJ_inv(n_p, n_J)
      , m_k_bbar(n_u, n_u)
      , m_A(n_J, n_u)
      , m_B(n_J, n_u)
      , m_C(n_p, n_u)
    {}

    void reset()
    {}
  };


  template <int dim>
  struct Solid<dim>::ScratchData_SC
  {
    void reset()
    {}
  };


  template <int dim>
  struct Solid<dim>::PerTaskData_UQPH
  {
    void reset()
    {}
  };


  template <int dim>
  struct Solid<dim>::ScratchData_UQPH
  {
    const BlockVector<double> & m_solution_total;

    std::vector<Tensor<2, dim>> m_solution_grads_u_total;
    std::vector<double>         m_solution_values_p_total;
    std::vector<double>         m_solution_values_J_total;

    FEValues<dim> m_fe_values;

    ScratchData_UQPH(const FiniteElement<dim> & fe_cell,
                     const QGauss<dim> &        qf_cell,
                     const UpdateFlags          uf_cell,
                     const BlockVector<double> &solution_total)
      : m_solution_total(solution_total)
      , m_solution_grads_u_total(qf_cell.size())
      , m_solution_values_p_total(qf_cell.size())
      , m_solution_values_J_total(qf_cell.size())
      , m_fe_values(fe_cell, qf_cell, uf_cell)
    {}

    ScratchData_UQPH(const ScratchData_UQPH &rhs)
      : m_solution_total(rhs.m_solution_total)
      , m_solution_grads_u_total(rhs.m_solution_grads_u_total)
      , m_solution_values_p_total(rhs.m_solution_values_p_total)
      , m_solution_values_J_total(rhs.m_solution_values_J_total)
      , m_fe_values(rhs.m_fe_values.get_fe(),
                  rhs.m_fe_values.get_quadrature(),
                  rhs.m_fe_values.get_update_flags())
    {}

    void reset()
    {
      const unsigned int n_q_points = m_solution_grads_u_total.size();
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          m_solution_grads_u_total[q]  = 0.0;
          m_solution_values_p_total[q] = 0.0;
          m_solution_values_J_total[q] = 0.0;
        }
    }
  };

  template <int dim>
  void Solid<dim>::refine_grid()
  {
    Vector<float> estimated_error_per_cell(m_triangulation.n_active_cells());

    FEValuesExtractors::Vector displacements(0);
    ComponentMask displacements_mask = m_fe.component_mask(displacements);

    KellyErrorEstimator<dim>::estimate(m_dof_handler,
                                       QGauss<dim - 1>(m_fe.degree + 1),
                                       {},
                                       m_solution_n,
                                       estimated_error_per_cell,
				       displacements_mask);

    GridRefinement::refine_and_coarsen_fixed_number(m_triangulation,
                                                    estimated_error_per_cell,
                                                    0.3,
                                                    0.03);

    m_triangulation.execute_coarsening_and_refinement();
  }

  template <int dim>
  void Solid<dim>::make_grid()
  {
    if (dim == 3)
      {
	if (m_parameters.m_scenario == 1)
	  {
	    make_grid_case_1();
	  }
	else if (m_parameters.m_scenario == 2)
	  {
	    make_grid_case_2();
	  }
	else if (m_parameters.m_scenario == 3)
	  {
	    make_grid_case_3();
	  }
	else if (m_parameters.m_scenario == 4)
	  {
	    make_grid_case_4();
	  }
	else if (m_parameters.m_scenario == 5)
	  {
	    make_grid_case_5();
	  }
	else if (m_parameters.m_scenario == 6)
	  {
	    make_grid_case_6();
	  }
	else if (m_parameters.m_scenario == 7)
	  {
	    make_grid_case_7();
	  }
	else
	  Assert(false, ExcMessage("The scenario has not been implemented!"));
      }
    else
      Assert(false, ExcMessage("Dimension has to be 3!"));

    std::ofstream out("original_mesh.vtu");
    GridOut       grid_out;
    grid_out.write_vtu(m_triangulation, out);

    m_triangulation.refine_global(m_parameters.m_global_refine);

    m_vol_reference = GridTools::volume(m_triangulation);
    std::cout << "Grid:\n\t Reference volume: " << m_vol_reference << std::endl;
  }


  template <int dim>
  void Solid<dim>::make_grid_case_1()
  {
    for (unsigned int i = 0; i < 80; ++i)
      std::cout << "*";
    std::cout << std::endl;
    std::cout << "Unit cube undergoes homogeneous deformation in the x-direction." << std::endl;
    for (unsigned int i = 0; i < 80; ++i)
      std::cout << "*";
    std::cout << std::endl;

    GridGenerator::hyper_rectangle(m_triangulation, Point<dim>(-0.5, -0.5, 0.0)
						  , Point<dim>( 0.5,  0.5, 1.0));
    /*
    for (const auto &cell : m_triangulation.active_cell_iterators())
      for (const auto &face : cell->face_iterators())
	{
	  if (face->at_boundary() == true)
	    {
	      if (std::fabs(face->center()[2] - 0.0) < 1.0e-9)
		face->set_boundary_id(6);
	      else if (std::fabs(face->center()[0] + 0.5) < 1.0e-9)
		face->set_boundary_id(7);
	      else if (std::fabs(face->center()[1] + 0.5) < 1.0e-9)
		face->set_boundary_id(8);
	      else if (std::fabs(face->center()[2] - 1.0) < 1.0e-9)
		face->set_boundary_id(9);
	      else
		face->set_boundary_id(0);
	    }
	}
    */
  }

  template <int dim>
  void Solid<dim>::make_grid_case_2()
  {
    for (unsigned int i = 0; i < 80; ++i)
      std::cout << "*";
    std::cout << std::endl;
    std::cout << "Beam undergoes small bending in the y-direction." << std::endl;
    for (unsigned int i = 0; i < 80; ++i)
      std::cout << "*";
    std::cout << std::endl;

    std::vector<unsigned int> repetitions(dim, 1);
    repetitions[0] = 50;
    repetitions[1] = 32;
    repetitions[2] = 8;

    GridGenerator::subdivided_hyper_rectangle(m_triangulation,
					      repetitions,
					      Point<dim>( 0.0,      0.0,     0.0    ),
					      Point<dim>( 17.5e-3,  1.75e-3, 5.0e-3 ) );
    for (const auto &cell : m_triangulation.active_cell_iterators())
      for (const auto &face : cell->face_iterators())
	{
	  if (face->at_boundary() == true)
	    {
	      if (std::fabs(face->center()[0] - 0.0) < 1.0e-9)
		face->set_boundary_id(6);
	      else
		face->set_boundary_id(0);
	    }
	}
  }

  template <int dim>
  void Solid<dim>::make_grid_case_3()
  {
    for (unsigned int i = 0; i < 80; ++i)
      std::cout << "*";
    std::cout << std::endl;
    std::cout << "A quarter of hollow cross" << std::endl;
    for (unsigned int i = 0; i < 80; ++i)
      std::cout << "*";
    std::cout << std::endl;

    const std::vector<Point<2>> vertices_1{
                                           {  0.0,  0.0},
                                           { 11.0,  0.0},
                                           {  0.0,  6.0},
                                           {  5.0,  6.0}
                                          };
    const std::vector<std::array<unsigned int,
                      GeometryInfo<2>::vertices_per_cell>>
        cell_vertices_1 = {
                           {{0, 1, 2, 3}},
                          };
    Triangulation<2> triangulation_2d_1;
    const unsigned int material_id_1 = 1;
    create_2d_grid(vertices_1, cell_vertices_1, triangulation_2d_1, material_id_1);
    Triangulation<3> triangulation_3d_1;
    GridGenerator::extrude_triangulation(triangulation_2d_1,
                                         2,
                                         0.41,
                                         triangulation_3d_1);

    const std::vector<Point<2>> vertices_2{
                                           {  5.0,  6.0},
                                           { 11.0,  0.0},
                                           {  5.0, 11.0},
                                           { 11.0, 11.0}
                                          };
    const std::vector<std::array<unsigned int,
                      GeometryInfo<2>::vertices_per_cell>>
        cell_vertices_2 = {
                           {{0, 1, 2, 3}},
                          };
    Triangulation<2> triangulation_2d_2;
    const unsigned int material_id_2 = 2;
    create_2d_grid(vertices_2, cell_vertices_2, triangulation_2d_2, material_id_2);
    Triangulation<3> triangulation_3d_2;
    GridGenerator::extrude_triangulation(triangulation_2d_2,
                                         2,
                                         0.41,
                                         triangulation_3d_2);
    GridGenerator::merge_triangulations(triangulation_3d_1,
                                        triangulation_3d_2,
                                        m_triangulation);

    for (const auto &cell : m_triangulation.active_cell_iterators())
      for (const auto &face : cell->face_iterators())
	{
	  if (face->at_boundary() == true)
	    {
	      if (std::fabs(face->center()[0] - 0.0) < 1.0e-9)
		face->set_boundary_id(6);
	      else
		face->set_boundary_id(0);
	    }
	}
  }

  template <int dim>
  void Solid<dim>::make_grid_case_4()
  {
    for (unsigned int i = 0; i < 80; ++i)
      std::cout << "*";
    std::cout << std::endl;
    std::cout << "Hollow cross" << std::endl;
    for (unsigned int i = 0; i < 80; ++i)
      std::cout << "*";
    std::cout << std::endl;

    Triangulation<2> triangulation_2d;

    GridIn<2> gridin;
    gridin.attach_triangulation(triangulation_2d);
    std::ifstream f("HollowCross.msh");
    gridin.read_msh(f);

    const double thickness = 0.41;
    const unsigned int n_layer = 4;
    GridGenerator::extrude_triangulation(triangulation_2d, n_layer, thickness, m_triangulation);
    /*
    for (const auto &cell : m_triangulation.active_cell_iterators())
      for (const auto &face : cell->face_iterators())
	{
	  if (face->at_boundary() == true)
	    {
	      if (    (std::fabs(face->center()[1] - 0.0) < 1.0e-9)
		   && (std::fabs(face->center()[0] - 0.0) < 0.25) )
		{
		  face->set_boundary_id(6);
		  std::cout << "1" << std::endl;
		}
	      else
		face->set_boundary_id(0);
	    }
	}
    */
  }

  template <int dim>
  void Solid<dim>::make_grid_case_5()
  {
    for (unsigned int i = 0; i < 80; ++i)
      std::cout << "*";
    std::cout << std::endl;
    std::cout << "Miura-ori fold" << std::endl;
    for (unsigned int i = 0; i < 80; ++i)
      std::cout << "*";
    std::cout << std::endl;

    Triangulation<2> triangulation_2d;

    GridIn<2> gridin;
    gridin.attach_triangulation(triangulation_2d);
    std::ifstream f("MiuraOri.msh");
    gridin.read_msh(f);

    const double thickness = 0.41;
    const unsigned int n_layer = 4;
    GridGenerator::extrude_triangulation(triangulation_2d, n_layer, thickness, m_triangulation);
  }

  template <int dim>
  void Solid<dim>::make_grid_case_6()
  {
    for (unsigned int i = 0; i < 80; ++i)
      std::cout << "*";
    std::cout << std::endl;
    std::cout << "Annulus" << std::endl;
    for (unsigned int i = 0; i < 80; ++i)
      std::cout << "*";
    std::cout << std::endl;

    Triangulation<2> triangulation_2d;

    GridIn<2> gridin;
    gridin.attach_triangulation(triangulation_2d);
    std::ifstream f("Annulus.msh");
    gridin.read_msh(f);

    const double thickness = 0.41;
    const unsigned int n_layer = 4;
    GridGenerator::extrude_triangulation(triangulation_2d, n_layer, thickness, m_triangulation);

    const Tensor<1, 3>           axis({0.0, 0.0, thickness});
    const Point<3>               axial_point(0.0, 0.0, 0.0);
    const CylindricalManifold<3> outer_cylinder(axis, axial_point);
    const types::manifold_id     outer_cylinder_id = 1;
    m_triangulation.set_manifold(outer_cylinder_id, outer_cylinder);
    const double outer_radius = 20.0;

    const CylindricalManifold<3> inner_cylinder(axis, axial_point);
    const types::manifold_id     inner_cylinder_id = 2;
    m_triangulation.set_manifold(inner_cylinder_id, inner_cylinder);
    const double inner_radius = 12.5;

    for (const auto &cell : m_triangulation.active_cell_iterators())
      for (const auto &face : cell->face_iterators() )
        {
	  const Point<3> face_center = face->center();
	  if (std::abs(face_center[0]) < outer_radius + 1.0e-9)
	    face->set_all_manifold_ids(outer_cylinder_id);
	  if (std::abs(face_center[0]) < inner_radius + 1.0e-9)
	    face->set_all_manifold_ids(inner_cylinder_id);
        }
  }

  template <int dim>
  void Solid<dim>::make_grid_case_7()
  {
    for (unsigned int i = 0; i < 80; ++i)
      std::cout << "*";
    std::cout << std::endl;
    std::cout << "Quadrupedal" << std::endl;
    for (unsigned int i = 0; i < 80; ++i)
      std::cout << "*";
    std::cout << std::endl;

    Triangulation<2> triangulation_2d;

    GridIn<2> gridin;
    gridin.attach_triangulation(triangulation_2d);
    std::ifstream f("Quadrupedal.msh");
    gridin.read_msh(f);

    const double thickness = 0.41;
    const unsigned int n_layer = 4;
    GridGenerator::extrude_triangulation(triangulation_2d, n_layer, thickness, m_triangulation);
  }


  template <int dim>
  void Solid<dim>::create_2d_grid(const std::vector<Point<2>> & vertices,
				  const std::vector<std::array<unsigned int,
				                    GeometryInfo<2>::vertices_per_cell>> & vertex_indices,
				  Triangulation<2> & coarse_grid,
				  const unsigned int material_id)
  {
    std::vector<CellData<2>> cells(vertex_indices.size());
    for (unsigned int i = 0; i < cells.size(); ++i)
      {
        for (unsigned int j = 0; j < vertex_indices[i].size(); ++j)
  	{
          cells[i].vertices[j] = vertex_indices[i][j];
          cells[i].material_id = material_id;
  	}
      }

    coarse_grid.create_triangulation(vertices, cells, SubCellData());
  }

  template <int dim>
  void Solid<dim>::system_setup()
  {
    m_timer.enter_subsection("Setup system");

    std::vector<unsigned int> block_component(m_n_components,
                                              m_u_dof); // Displacement
    block_component[m_p_component] = m_p_dof;             // Pressure
    block_component[m_J_component] = m_J_dof;             // Dilatation

    m_dof_handler.distribute_dofs(m_fe);
    DoFRenumbering::Cuthill_McKee(m_dof_handler);
    DoFRenumbering::component_wise(m_dof_handler, block_component);

    m_constraints.clear();
    DoFTools::make_hanging_node_constraints(m_dof_handler, m_constraints);
    m_constraints.close();

    m_dofs_per_block =
      DoFTools::count_dofs_per_fe_block(m_dof_handler, block_component);

    std::cout << "Triangulation:"
              << "\n\t Number of active cells: "
              << m_triangulation.n_active_cells()
              << "\n\t Number of degrees of freedom: " << m_dof_handler.n_dofs()
              << std::endl;

    m_tangent_matrix.clear();
    {
      BlockDynamicSparsityPattern dsp(m_dofs_per_block, m_dofs_per_block);

      Table<2, DoFTools::Coupling> coupling(m_n_components, m_n_components);
      for (unsigned int ii = 0; ii < m_n_components; ++ii)
        for (unsigned int jj = 0; jj < m_n_components; ++jj)
          if (((ii < m_p_component) && (jj == m_J_component)) ||
              ((ii == m_J_component) && (jj < m_p_component)) ||
              ((ii == m_p_component) && (jj == m_p_component)))
            coupling[ii][jj] = DoFTools::none;
          else
            coupling[ii][jj] = DoFTools::always;
      DoFTools::make_sparsity_pattern(
        m_dof_handler, coupling, dsp, m_constraints, false);
      m_sparsity_pattern.copy_from(dsp);
    }

    m_tangent_matrix.reinit(m_sparsity_pattern);

    m_system_rhs.reinit(m_dofs_per_block);
    m_solution_n.reinit(m_dofs_per_block);

    setup_qph();

    m_timer.leave_subsection();
  }


  template <int dim>
  void Solid<dim>::determine_component_extractors()
  {
    m_element_indices_u.clear();
    m_element_indices_p.clear();
    m_element_indices_J.clear();

    for (unsigned int k = 0; k < m_fe.n_dofs_per_cell(); ++k)
      {
        const unsigned int k_group = m_fe.system_to_base_index(k).first.first;
        if (k_group == m_u_dof)
          m_element_indices_u.push_back(k);
        else if (k_group == m_p_dof)
          m_element_indices_p.push_back(k);
        else if (k_group == m_J_dof)
          m_element_indices_J.push_back(k);
        else
          {
            Assert(k_group <= m_J_dof, ExcInternalError());
          }
      }
  }

  template <int dim>
  void Solid<dim>::setup_qph()
  {
    std::cout << "     Setting up quadrature point data ("
	      << m_n_q_points
	      << " points per cell)" << std::endl;

    m_quadrature_point_history.clear();
    for (auto const & cell : m_triangulation.active_cell_iterators())
      {
	m_quadrature_point_history.initialize(cell, m_n_q_points);
      }

    unsigned int material_id;
    double shear_modulus = 0.0;
    double bulk_modulus = 0.0;
    Tensor<1, 3> Br_tilde;

    for (const auto &cell : m_triangulation.active_cell_iterators())
      {
        material_id = cell->material_id();
        if (m_material_data.find(material_id) != m_material_data.end())
          {
	    shear_modulus = m_material_data[material_id][0];
	    bulk_modulus  = m_material_data[material_id][1];
	    Br_tilde[0]   = m_material_data[material_id][2];
	    Br_tilde[1]   = m_material_data[material_id][3];
	    Br_tilde[2]   = m_material_data[material_id][4];
	  }
        else
          {
            std::cout << "Could not find material data for material id: " << material_id << std::endl;
            Assert(false, ExcMessage("Could not find material data for material id."));
          }

        const std::vector<std::shared_ptr<PointHistory<dim>>> lqph =
          m_quadrature_point_history.get_data(cell);
        Assert(lqph.size() == m_n_q_points, ExcInternalError());

        for (unsigned int q_point = 0; q_point < m_n_q_points; ++q_point)
          lqph[q_point]->setup_lqp(m_parameters, shear_modulus, bulk_modulus, Br_tilde);
      }
  }

  template <int dim>
  void
  Solid<dim>::update_qph_incremental(const BlockVector<double> &solution_delta)
  {
    m_timer.enter_subsection("Update QPH data");
    std::cout << " UQPH " << std::flush;

    const BlockVector<double> solution_total(
      get_total_solution(solution_delta));

    const UpdateFlags uf_UQPH(update_values | update_gradients);
    PerTaskData_UQPH  per_task_data_UQPH;
    ScratchData_UQPH  scratch_data_UQPH(m_fe, m_qf_cell, uf_UQPH, solution_total);

    WorkStream::run(m_dof_handler.active_cell_iterators(),
                    *this,
                    &Solid::update_qph_incremental_one_cell,
                    &Solid::copy_local_to_global_UQPH,
                    scratch_data_UQPH,
                    per_task_data_UQPH);

    m_timer.leave_subsection();
  }


  template <int dim>
  void Solid<dim>::update_qph_incremental_one_cell(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    ScratchData_UQPH &                                    scratch,
    PerTaskData_UQPH & /*data*/)
  {
    const std::vector<std::shared_ptr<PointHistory<dim>>> lqph =
      m_quadrature_point_history.get_data(cell);
    Assert(lqph.size() == m_n_q_points, ExcInternalError());

    Assert(scratch.m_solution_grads_u_total.size() == m_n_q_points,
           ExcInternalError());
    Assert(scratch.m_solution_values_p_total.size() == m_n_q_points,
           ExcInternalError());
    Assert(scratch.m_solution_values_J_total.size() == m_n_q_points,
           ExcInternalError());

    scratch.reset();

    scratch.m_fe_values.reinit(cell);
    scratch.m_fe_values[m_u_fe].get_function_gradients(
      scratch.m_solution_total, scratch.m_solution_grads_u_total);
    scratch.m_fe_values[m_p_fe].get_function_values(
      scratch.m_solution_total, scratch.m_solution_values_p_total);
    scratch.m_fe_values[m_J_fe].get_function_values(
      scratch.m_solution_total, scratch.m_solution_values_J_total);

    for (const unsigned int q_point :
         scratch.m_fe_values.quadrature_point_indices())
      lqph[q_point]->update_values(scratch.m_solution_grads_u_total[q_point],
                                   scratch.m_solution_values_p_total[q_point],
                                   scratch.m_solution_values_J_total[q_point]);
  }



  template <int dim>
  void Solid<dim>::solve_nonlinear_timestep(BlockVector<double> &solution_delta)
  {
    std::cout << std::endl
              << "Timestep " << m_time.get_timestep() << " @ " << m_time.current()
              << 's' << std::endl;

    BlockVector<double> newton_update(m_dofs_per_block);

    m_error_residual.reset();
    m_error_residual_0.reset();
    m_error_residual_norm.reset();
    m_error_update.reset();
    m_error_update_0.reset();
    m_error_update_norm.reset();

    print_conv_header();

    unsigned int newton_iteration = 0;
    for (; newton_iteration < m_parameters.m_max_iterations_NR; ++newton_iteration)
      {
        std::cout << ' ' << std::setw(2) << newton_iteration << ' '
                  << std::flush;

        make_constraints(newton_iteration);
        assemble_system();

        get_error_residual(m_error_residual);
        if (newton_iteration == 0)
          m_error_residual_0 = m_error_residual;

        m_error_residual_norm = m_error_residual;
        m_error_residual_norm.normalize(m_error_residual_0);

        if (newton_iteration > 0 && m_error_update_norm.m_u <= m_parameters.m_tol_u &&
            m_error_residual_norm.m_u <= m_parameters.m_tol_f)
          {
            std::cout << " CONVERGED! " << std::endl;
            print_conv_footer();

            break;
          }

        const std::pair<unsigned int, double> lin_solver_output =
          solve_linear_system(newton_update);

        get_error_update(newton_update, m_error_update);
        if (newton_iteration == 0)
          m_error_update_0 = m_error_update;

        m_error_update_norm = m_error_update;
        m_error_update_norm.normalize(m_error_update_0);

        solution_delta += newton_update;
        update_qph_incremental(solution_delta);

        std::cout << " | " << std::fixed << std::setprecision(3) << std::setw(7)
                  << std::scientific << lin_solver_output.first << "  "
                  << lin_solver_output.second << "  "
                  << m_error_residual_norm.m_norm << "  " << m_error_residual_norm.m_u
                  << "  " << m_error_residual_norm.m_p << "  "
                  << m_error_residual_norm.m_J << "  " << m_error_update_norm.m_norm
                  << "  " << m_error_update_norm.m_u << "  " << m_error_update_norm.m_p
                  << "  " << m_error_update_norm.m_J << "  " << std::endl;
      }

    AssertThrow(newton_iteration < m_parameters.m_max_iterations_NR,
                ExcMessage("No convergence in nonlinear solver!"));
  }



  template <int dim>
  void Solid<dim>::print_conv_header()
  {
    static const unsigned int l_width = 150;

    for (unsigned int i = 0; i < l_width; ++i)
      std::cout << '_';
    std::cout << std::endl;

    std::cout << "               SOLVER STEP               "
              << " |  LIN_IT   LIN_RES    RES_NORM    "
              << " RES_U     RES_P      RES_J     NU_NORM     "
              << " NU_U       NU_P       NU_J " << std::endl;

    for (unsigned int i = 0; i < l_width; ++i)
      std::cout << '_';
    std::cout << std::endl;
  }



  template <int dim>
  void Solid<dim>::print_conv_footer()
  {
    static const unsigned int l_width = 150;

    for (unsigned int i = 0; i < l_width; ++i)
      std::cout << '_';
    std::cout << std::endl;

    const std::pair<double, double> error_dil = get_error_dilation();

    std::cout << "Relative errors:" << std::endl
              << "Displacement:\t" << m_error_update.m_u / m_error_update_0.m_u
              << std::endl
              << "Force: \t\t" << m_error_residual.m_u / m_error_residual_0.m_u
              << std::endl
              << "Dilatation:\t" << error_dil.first << std::endl
              << "v / V_0:\t" << error_dil.second * m_vol_reference << " / "
              << m_vol_reference << " = " << error_dil.second << std::endl;

    const double time_ramp = (m_time.current() / m_time.end());
    compute_energy_current(m_energy);

    std::cout << "Strain energy of the body:" << std::endl;
    std::cout << "External magnetic field:\t" << "(" << time_ramp*m_parameters.m_B_applied << ")"<< std::endl;
    std::cout << "Total strain energy:\t\t" << m_energy.m_total_strain_energy << std::endl;
    std::cout << "Volumetric strain energy:\t" << m_energy.m_volumetric_strain_energy << std::endl;
    std::cout << "Isotropic strain energy:\t" << m_energy.m_isotropic_strain_energy << std::endl;
    std::cout << "Mechanical strain energy:\t" << m_energy.m_mechanical_strain_energy << std::endl;
    std::cout << "Magnetic strain energy:\t\t" << m_energy.m_magnetic_strain_energy << std::endl;
  }

  template <int dim>
  double Solid<dim>::compute_vol_current() const
  {
    double vol_current = 0.0;

    FEValues<dim> fe_values(m_fe, m_qf_cell, update_JxW_values);

    for (const auto &cell : m_triangulation.active_cell_iterators())
      {
        fe_values.reinit(cell);

        const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
          m_quadrature_point_history.get_data(cell);
        Assert(lqph.size() == m_n_q_points, ExcInternalError());

        for (const unsigned int q_point : fe_values.quadrature_point_indices())
          {
            const double det_F_qp = lqph[q_point]->get_det_F();
            const double JxW      = fe_values.JxW(q_point);

            vol_current += det_F_qp * JxW;
          }
      }
    Assert(vol_current > 0.0, ExcInternalError());
    return vol_current;
  }

  template <int dim>
  void Solid<dim>::compute_energy_current(Energy &energy) const
  {
    double total_strain_energy = 0.0;
    double volumetric_strain_energy = 0.0;
    double isotropic_strain_energy = 0.0;
    double mechanical_strain_energy = 0.0;
    double magnetic_strain_energy = 0.0;

    FEValues<dim> fe_values(m_fe, m_qf_cell, update_JxW_values);

    for (const auto &cell : m_triangulation.active_cell_iterators())
      {
        fe_values.reinit(cell);

        const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
          m_quadrature_point_history.get_data(cell);
        Assert(lqph.size() == m_n_q_points, ExcInternalError());

        for (const unsigned int q_point : fe_values.quadrature_point_indices())
          {
            const double total_energy_gp = lqph[q_point]->get_total_strain_energy();
            const double volumetric_energy_gp = lqph[q_point]->get_volumetric_strain_energy();
            const double isotropic_energy_gp = lqph[q_point]->get_isotropic_strain_energy();
            const double mechanical_energy_gp = lqph[q_point]->get_mechanical_strain_energy();
            const double magnetic_energy_gp = lqph[q_point]->get_magnetic_strain_energy();
            const double det_F_qp = lqph[q_point]->get_det_F();
            const double JxW      = fe_values.JxW(q_point);

            total_strain_energy += total_energy_gp * det_F_qp * JxW;
            volumetric_strain_energy += volumetric_energy_gp * det_F_qp * JxW;
            isotropic_strain_energy += isotropic_energy_gp * det_F_qp * JxW;
            mechanical_strain_energy += mechanical_energy_gp * det_F_qp * JxW;
            magnetic_strain_energy += magnetic_energy_gp * det_F_qp * JxW;
          }
      }

    energy.m_total_strain_energy = total_strain_energy;
    energy.m_volumetric_strain_energy = volumetric_strain_energy;
    energy.m_isotropic_strain_energy = isotropic_strain_energy;
    energy.m_mechanical_strain_energy = mechanical_strain_energy;
    energy.m_magnetic_strain_energy = magnetic_strain_energy;
  }


  template <int dim>
  std::pair<double, double> Solid<dim>::get_error_dilation() const
  {
    double dil_L2_error = 0.0;

    FEValues<dim> fe_values(m_fe, m_qf_cell, update_JxW_values);

    for (const auto &cell : m_triangulation.active_cell_iterators())
      {
        fe_values.reinit(cell);

        const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
          m_quadrature_point_history.get_data(cell);
        Assert(lqph.size() == m_n_q_points, ExcInternalError());

        for (const unsigned int q_point : fe_values.quadrature_point_indices())
          {
            const double det_F_qp   = lqph[q_point]->get_det_F();
            const double J_tilde_qp = lqph[q_point]->get_J_tilde();
            const double the_error_qp_squared =
              std::pow((det_F_qp - J_tilde_qp), 2);
            const double JxW = fe_values.JxW(q_point);

            dil_L2_error += the_error_qp_squared * JxW;
          }
      }

    return std::make_pair(std::sqrt(dil_L2_error),
                          compute_vol_current() / m_vol_reference);
  }



  template <int dim>
  void Solid<dim>::get_error_residual(Errors &error_residual)
  {
    BlockVector<double> error_res(m_dofs_per_block);

    for (unsigned int i = 0; i < m_dof_handler.n_dofs(); ++i)
      if (!m_constraints.is_constrained(i))
        error_res(i) = m_system_rhs(i);

    error_residual.m_norm = error_res.l2_norm();
    error_residual.m_u    = error_res.block(m_u_dof).l2_norm();
    error_residual.m_p    = error_res.block(m_p_dof).l2_norm();
    error_residual.m_J    = error_res.block(m_J_dof).l2_norm();
  }



  template <int dim>
  void Solid<dim>::get_error_update(const BlockVector<double> &newton_update,
                                    Errors &                   error_update)
  {
    BlockVector<double> error_ud(m_dofs_per_block);
    for (unsigned int i = 0; i < m_dof_handler.n_dofs(); ++i)
      if (!m_constraints.is_constrained(i))
        error_ud(i) = newton_update(i);

    error_update.m_norm = error_ud.l2_norm();
    error_update.m_u    = error_ud.block(m_u_dof).l2_norm();
    error_update.m_p    = error_ud.block(m_p_dof).l2_norm();
    error_update.m_J    = error_ud.block(m_J_dof).l2_norm();
  }




  template <int dim>
  BlockVector<double> Solid<dim>::get_total_solution(
    const BlockVector<double> &solution_delta) const
  {
    BlockVector<double> solution_total(m_solution_n);
    solution_total += solution_delta;
    return solution_total;
  }



  template <int dim>
  void Solid<dim>::assemble_system()
  {
    m_timer.enter_subsection("Assemble system");
    std::cout << " ASM_SYS " << std::flush;

    m_tangent_matrix = 0.0;
    m_system_rhs     = 0.0;

    const UpdateFlags uf_cell(update_values | update_gradients |
			      update_quadrature_points | update_JxW_values);
    const UpdateFlags uf_face(update_values | update_normal_vectors |
                              update_JxW_values);

    PerTaskData_ASM per_task_data(m_dofs_per_cell);
    ScratchData_ASM scratch_data(m_fe, m_qf_cell, uf_cell, m_qf_face, uf_face);

    WorkStream::run(
      m_dof_handler.active_cell_iterators(),
      [this](const typename DoFHandler<dim>::active_cell_iterator &cell,
             ScratchData_ASM &                                     scratch,
             PerTaskData_ASM &                                     data) {
        this->assemble_system_one_cell(cell, scratch, data);
      },
      [this](const PerTaskData_ASM &data) {
        this->m_constraints.distribute_local_to_global(data.m_cell_matrix,
                                                     data.m_cell_rhs,
                                                     data.m_local_dof_indices,
                                                     m_tangent_matrix,
                                                     m_system_rhs);
      },
      scratch_data,
      per_task_data);

    m_timer.leave_subsection();
  }

  template <int dim>
  void Solid<dim>::assemble_system_one_cell(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    ScratchData_ASM &                                     scratch,
    PerTaskData_ASM &                                     data) const
  {
    data.reset();
    scratch.reset();
    scratch.m_fe_values.reinit(cell);
    cell->get_dof_indices(data.m_local_dof_indices);

    const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
      m_quadrature_point_history.get_data(cell);
    Assert(lqph.size() == m_n_q_points, ExcInternalError());

    const double time_ramp = (m_time.current() / m_time.end());
    std::vector<Tensor<1, dim>> rhs_values(m_n_q_points);

    right_hand_side(scratch.m_fe_values.get_quadrature_points(),
		    rhs_values,
		    m_parameters.m_x_component*time_ramp,
		    m_parameters.m_y_component*time_ramp,
		    m_parameters.m_z_component*time_ramp);

    for (const unsigned int q_point :
         scratch.m_fe_values.quadrature_point_indices())
      {
        const Tensor<2, dim> F_inv = lqph[q_point]->get_F_inv();
        for (const unsigned int k : scratch.m_fe_values.dof_indices())
          {
            const unsigned int k_group = m_fe.system_to_base_index(k).first.first;

            if (k_group == m_u_dof)
              {
                scratch.m_Nx_disp[q_point][k] =
                  scratch.m_fe_values[m_u_fe].value(k, q_point);
                scratch.m_grad_Nx[q_point][k] =
                  scratch.m_fe_values[m_u_fe].gradient(k, q_point) * F_inv;
                scratch.m_symm_grad_Nx[q_point][k] =
                  symmetrize(scratch.m_grad_Nx[q_point][k]);
              }
            else if (k_group == m_p_dof)
              scratch.m_Nx[q_point][k] =
                scratch.m_fe_values[m_p_fe].value(k, q_point);
            else if (k_group == m_J_dof)
              scratch.m_Nx[q_point][k] =
                scratch.m_fe_values[m_J_fe].value(k, q_point);
            else
              Assert(k_group <= m_J_dof, ExcInternalError());
          }
      }

    for (const unsigned int q_point :
         scratch.m_fe_values.quadrature_point_indices())
      {
        const SymmetricTensor<2, dim> tau     = lqph[q_point]->get_tau();
        const SymmetricTensor<4, dim> Jc      = lqph[q_point]->get_Jc();
        const double                  det_F   = lqph[q_point]->get_det_F();
        const double                  p_tilde = lqph[q_point]->get_p_tilde();
        const double                  J_tilde = lqph[q_point]->get_J_tilde();
        const double dPsi_vol_dJ   = lqph[q_point]->get_dPsi_vol_dJ();
        const double d2Psi_vol_dJ2 = lqph[q_point]->get_d2Psi_vol_dJ2();
        const SymmetricTensor<2, dim> &I =
          Physics::Elasticity::StandardTensors<dim>::I;

        const Tensor<2, dim> tau_magnetic = lqph[q_point]->get_tau_magnetic() * time_ramp;
	const Tensor<4, dim> Jc_magnetic = lqph[q_point]->get_Jc_magnetic() * time_ramp;

        Tensor<2, dim> grad_Nx_i_x_Jc;
        Tensor<1, dim> grad_Nx_i_comp_i_x_tau;

        const std::vector<double> &                 N = scratch.m_Nx[q_point];
        const std::vector<Tensor<1,dim>> &          N_disp = scratch.m_Nx_disp[q_point];

        const std::vector<SymmetricTensor<2, dim>> &symm_grad_Nx =
          scratch.m_symm_grad_Nx[q_point];
        const std::vector<Tensor<2, dim>> &grad_Nx = scratch.m_grad_Nx[q_point];
        const double                       JxW = scratch.m_fe_values.JxW(q_point);

        for (const unsigned int i : scratch.m_fe_values.dof_indices())
          {
            const unsigned int component_i =
              m_fe.system_to_component_index(i).first;
            const unsigned int i_group = m_fe.system_to_base_index(i).first.first;

            if (i_group == m_u_dof)
              {
                data.m_cell_rhs(i) -= scalar_product(grad_Nx[i], tau + tau_magnetic) * JxW;

		// contributions from the body force to right-hand side
		data.m_cell_rhs(i) += N_disp[i] * rhs_values[q_point] * JxW;
              }
            else if (i_group == m_p_dof)
              data.m_cell_rhs(i) -= N[i] * (det_F - J_tilde) * JxW;
            else if (i_group == m_J_dof)
              data.m_cell_rhs(i) -= N[i] * (dPsi_vol_dJ - p_tilde) * JxW;
            else
              Assert(i_group <= m_J_dof, ExcInternalError());

            if (i_group == m_u_dof)
              {
		grad_Nx_i_x_Jc = double_contract<0,0,1,1>(grad_Nx[i], Jc + Jc_magnetic);
                grad_Nx_i_comp_i_x_tau = grad_Nx[i][component_i] * tau;
              }

            for (const unsigned int j : scratch.m_fe_values.dof_indices())
              {
                const unsigned int component_j =
                  m_fe.system_to_component_index(j).first;
                const unsigned int j_group =
                  m_fe.system_to_base_index(j).first.first;

                if ((i_group == j_group) && (i_group == m_u_dof))
                  {
		    data.m_cell_matrix(i, j) += scalar_product(grad_Nx_i_x_Jc,
		    		                               grad_Nx[j]      ) * JxW;

                    if (component_i == component_j)
                      data.m_cell_matrix(i, j) +=
                        grad_Nx_i_comp_i_x_tau * grad_Nx[j][component_j] * JxW;
                  }
                else if ((i_group == m_p_dof) && (j_group == m_u_dof))
                  {
                    data.m_cell_matrix(i, j) += N[i] * det_F *               //
                                              (symm_grad_Nx[j] * I) * JxW; //
                    data.m_cell_matrix(j, i) = data.m_cell_matrix(i, j);
                  }
                else if ((i_group == m_J_dof) && (j_group == m_p_dof))
                  {
                    data.m_cell_matrix(i, j) -= N[i] * N[j] * JxW;
                    data.m_cell_matrix(j, i) = data.m_cell_matrix(i, j);
                  }
                else if ((i_group == j_group) && (i_group == m_J_dof))
                  {
                    data.m_cell_matrix(i, j) += N[i] * d2Psi_vol_dJ2 * N[j] * JxW;
                  }
                else
                  Assert((i_group <= m_J_dof) && (j_group <= m_J_dof),
                         ExcInternalError());
              }
          }
      }

    // if there is surface pressure
    const unsigned int face_pressure_id = 100;
    const double p0 = 0.0;
    for (const auto &face : cell->face_iterators())
      if (face->at_boundary() && face->boundary_id() == face_pressure_id)
        {
          scratch.m_fe_face_values.reinit(cell, face);

          for (const unsigned int f_q_point :
               scratch.m_fe_face_values.quadrature_point_indices())
            {
              const Tensor<1, dim> &N =
                scratch.m_fe_face_values.normal_vector(f_q_point);

              const double         pressure  = p0 * time_ramp;
              const Tensor<1, dim> traction  = pressure * N;

              for (const unsigned int i : scratch.m_fe_values.dof_indices())
                {
                  const unsigned int i_group =
                    m_fe.system_to_base_index(i).first.first;

                  if (i_group == m_u_dof)
                    {
                      const unsigned int component_i =
                        m_fe.system_to_component_index(i).first;
                      const double Ni =
                        scratch.m_fe_face_values.shape_value(i, f_q_point);
                      const double JxW = scratch.m_fe_face_values.JxW(f_q_point);

                      data.m_cell_rhs(i) += (Ni * traction[component_i]) * JxW;
                    }
                }
            }
        }
  }

  template <int dim>
  void Solid<dim>::make_constraints(const unsigned int it_nr)
  {
    const bool apply_dirichlet_bc = (it_nr == 0);

    if (it_nr > 1)
      {
        std::cout << " --- " << std::flush;
        return;
      }

    std::cout << " CST " << std::flush;

    if (apply_dirichlet_bc)
      {
	m_constraints.clear();
	DoFTools::make_hanging_node_constraints(m_dof_handler, m_constraints);
	const FEValuesExtractors::Scalar x_displacement(0);
	const FEValuesExtractors::Scalar y_displacement(1);
	const FEValuesExtractors::Scalar z_displacement(2);

	if (m_parameters.m_scenario == 1)
	  {
	    // Dirichlet B.C. 1
	    /*
	    int boundary_id = 6;
	    VectorTools::interpolate_boundary_values(m_dof_handler,
						     boundary_id,
						     Functions::ZeroFunction<dim>(dim),
						     m_constraints,
						     m_fe.component_mask(z_displacement));

	    boundary_id = 7;
	    VectorTools::interpolate_boundary_values(m_dof_handler,
						     boundary_id,
						     Functions::ZeroFunction<dim>(dim),
						     m_constraints,
						     m_fe.component_mask(x_displacement));

	    boundary_id = 8;
	    VectorTools::interpolate_boundary_values(m_dof_handler,
						     boundary_id,
						     Functions::ZeroFunction<dim>(dim),
						     m_constraints,
						     m_fe.component_mask(y_displacement));
	    */

	    Triangulation<3>::active_vertex_iterator vertex_itr;
	    vertex_itr = m_triangulation.begin_active_vertex();
	    std::vector<types::global_dof_index> node_lower_left_dofs(m_fe.dofs_per_vertex);
	    std::vector<types::global_dof_index> node_lower_right_dofs(m_fe.dofs_per_vertex);
	    std::vector<types::global_dof_index> node_upper_left_dofs(m_fe.dofs_per_vertex);

	    for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)
	      {
		if (   (std::fabs(vertex_itr->vertex()[0] + 0.5) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] + 0.5) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] - 0.0) < 1.0e-9) )
		  {
		    node_lower_left_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
		if (   (std::fabs(vertex_itr->vertex()[0] - 0.5) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] + 0.5) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] - 0.0) < 1.0e-9) )
		  {
		    node_lower_right_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
		if (   (std::fabs(vertex_itr->vertex()[0] + 0.5) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] - 0.5) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] - 0.0) < 1.0e-9) )
		  {
		    node_upper_left_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
	      }

	    // node at lower left corner is fixed at all three directions
	    m_constraints.add_line(node_lower_left_dofs[0]);
	    m_constraints.set_inhomogeneity(node_lower_left_dofs[0], 0.0);

	    m_constraints.add_line(node_lower_left_dofs[1]);
	    m_constraints.set_inhomogeneity(node_lower_left_dofs[1], 0.0);

	    m_constraints.add_line(node_lower_left_dofs[2]);
	    m_constraints.set_inhomogeneity(node_lower_left_dofs[2], 0.0);

	    // node at lower right corner is fixed at the y and z directions
	    m_constraints.add_line(node_lower_right_dofs[1]);
	    m_constraints.set_inhomogeneity(node_lower_right_dofs[1], 0.0);

	    m_constraints.add_line(node_lower_right_dofs[2]);
	    m_constraints.set_inhomogeneity(node_lower_right_dofs[2], 0.0);

	    // node at upper left corner is fixed at the x and z directions
	    m_constraints.add_line(node_upper_left_dofs[0]);
	    m_constraints.set_inhomogeneity(node_upper_left_dofs[0], 0.0);

	    m_constraints.add_line(node_upper_left_dofs[2]);
	    m_constraints.set_inhomogeneity(node_upper_left_dofs[2], 0.0);

	    m_constraints.close();
	  }
	else if (m_parameters.m_scenario == 2)
	  {
	    int boundary_id = 6;
	    VectorTools::interpolate_boundary_values(m_dof_handler,
						     boundary_id,
						     Functions::ZeroFunction<dim>(m_n_components),
						     m_constraints,
						     (m_fe.component_mask(x_displacement) |
						      m_fe.component_mask(y_displacement) |
						      m_fe.component_mask(z_displacement))  );

	  }
	else if (m_parameters.m_scenario == 3)
	  {
	    int boundary_id = 6;
	    VectorTools::interpolate_boundary_values(m_dof_handler,
						     boundary_id,
						     Functions::ZeroFunction<dim>(dim),
						     m_constraints);
	  }
	else if (m_parameters.m_scenario == 4)
	  {
	    /*
	    int boundary_id = 6;
	    VectorTools::interpolate_boundary_values(m_dof_handler,
						     boundary_id,
						     Functions::ZeroFunction<dim>(dim),
						     m_constraints);
	    */
	    Triangulation<3>::active_vertex_iterator vertex_itr;
	    vertex_itr = m_triangulation.begin_active_vertex();
	    std::vector<types::global_dof_index> node_bottom_center_dofs(m_fe.dofs_per_vertex);
	    std::vector<types::global_dof_index> node_right_center_dofs(m_fe.dofs_per_vertex);
	    std::vector<types::global_dof_index> node_top_center_dofs(m_fe.dofs_per_vertex);

	    for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)
	      {
		if (   (std::fabs(vertex_itr->vertex()[0] - 0.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] - 0.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] - 0.0) < 1.0e-9) )
		  {
		    node_bottom_center_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
		if (   (std::fabs(vertex_itr->vertex()[0] -27.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] -27.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] - 0.0) < 1.0e-9) )
		  {
		    node_right_center_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
		if (   (std::fabs(vertex_itr->vertex()[0] - 0.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] -54.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] - 0.0) < 1.0e-9) )
		  {
		    node_top_center_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
	      }

	    m_constraints.add_line(node_bottom_center_dofs[0]);
	    m_constraints.set_inhomogeneity(node_bottom_center_dofs[0], 0.0);
	    m_constraints.add_line(node_bottom_center_dofs[1]);
	    m_constraints.set_inhomogeneity(node_bottom_center_dofs[1], 0.0);
	    m_constraints.add_line(node_bottom_center_dofs[2]);
	    m_constraints.set_inhomogeneity(node_bottom_center_dofs[2], 0.0);

	    m_constraints.add_line(node_right_center_dofs[2]);
	    m_constraints.set_inhomogeneity(node_right_center_dofs[2], 0.0);

	    m_constraints.add_line(node_top_center_dofs[0]);
	    m_constraints.set_inhomogeneity(node_top_center_dofs[0], 0.0);
	    m_constraints.add_line(node_top_center_dofs[2]);
	    m_constraints.set_inhomogeneity(node_top_center_dofs[2], 0.0);

	    m_constraints.close();

	  }
	else if (m_parameters.m_scenario == 5)
	  {
	    const double pi = std::acos(-1);
	    const double b = 17.0 / std::tan((56 * pi)/180);

	    Triangulation<3>::active_vertex_iterator vertex_itr;
	    vertex_itr = m_triangulation.begin_active_vertex();
	    std::vector<types::global_dof_index> node_center_dofs(m_fe.dofs_per_vertex);
	    std::vector<types::global_dof_index> node_left_corner_1_dofs(m_fe.dofs_per_vertex);
	    std::vector<types::global_dof_index> node_left_corner_2_dofs(m_fe.dofs_per_vertex);
	    std::vector<types::global_dof_index> node_right_corner_1_dofs(m_fe.dofs_per_vertex);
	    std::vector<types::global_dof_index> node_right_corner_2_dofs(m_fe.dofs_per_vertex);
	    std::vector<types::global_dof_index> node_right_center_dofs(m_fe.dofs_per_vertex);


	    for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)
	      {
		if (   (std::fabs(vertex_itr->vertex()[0] - 34.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] - 34.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] -  0.0) < 1.0e-9) )
		  {
		    node_center_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
		if (   (std::fabs(vertex_itr->vertex()[0] +    b) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] - 17.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] -  0.0) < 1.0e-9) )
		  {
		    node_left_corner_1_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
		if (   (std::fabs(vertex_itr->vertex()[0] +    b) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] - 51.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] -  0.0) < 1.0e-9) )
		  {
		    node_left_corner_2_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
		if (   (std::fabs(vertex_itr->vertex()[0] - 68.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] -  0.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] -  0.0) < 1.0e-9) )
		  {
		    node_right_corner_1_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
		if (   (std::fabs(vertex_itr->vertex()[0] - 68.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] - 68.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] -  0.0) < 1.0e-9) )
		  {
		    node_right_corner_2_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
		if (   (std::fabs(vertex_itr->vertex()[0] - 68.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] - 34.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] -  0.0) < 1.0e-9) )
		  {
		    node_right_center_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
	      }

	    m_constraints.add_line(node_center_dofs[0]);
	    m_constraints.set_inhomogeneity(node_center_dofs[0], 0.0);
	    m_constraints.add_line(node_center_dofs[1]);
	    m_constraints.set_inhomogeneity(node_center_dofs[1], 0.0);
	    //m_constraints.add_line(node_center_dofs[2]);
	    //m_constraints.set_inhomogeneity(node_center_dofs[2], 0.0);

	    m_constraints.add_line(node_left_corner_1_dofs[2]);
	    m_constraints.set_inhomogeneity(node_left_corner_1_dofs[2], 0.0);
	    m_constraints.add_line(node_left_corner_2_dofs[2]);
	    m_constraints.set_inhomogeneity(node_left_corner_2_dofs[2], 0.0);
	    m_constraints.add_line(node_right_corner_1_dofs[2]);
	    m_constraints.set_inhomogeneity(node_right_corner_1_dofs[2], 0.0);
	    m_constraints.add_line(node_right_corner_2_dofs[2]);
	    m_constraints.set_inhomogeneity(node_right_corner_2_dofs[2], 0.0);

	    m_constraints.add_line(node_right_center_dofs[1]);
	    m_constraints.set_inhomogeneity(node_right_center_dofs[1], 0.0);
	    m_constraints.add_line(node_right_center_dofs[2]);
	    m_constraints.set_inhomogeneity(node_right_center_dofs[2], 0.0);

	    m_constraints.close();

	  }
	else if (m_parameters.m_scenario == 6)
	  {
	    Triangulation<3>::active_vertex_iterator vertex_itr;
	    vertex_itr = m_triangulation.begin_active_vertex();
	    std::vector<types::global_dof_index> node_bottom_center_1_dofs(m_fe.dofs_per_vertex);
	    std::vector<types::global_dof_index> node_bottom_center_2_dofs(m_fe.dofs_per_vertex);
	    std::vector<types::global_dof_index> node_right_center_1_dofs(m_fe.dofs_per_vertex);
	    std::vector<types::global_dof_index> node_right_center_2_dofs(m_fe.dofs_per_vertex);
	    std::vector<types::global_dof_index> node_left_center_1_dofs(m_fe.dofs_per_vertex);
	    std::vector<types::global_dof_index> node_left_center_2_dofs(m_fe.dofs_per_vertex);
	    std::vector<types::global_dof_index> node_top_center_1_dofs(m_fe.dofs_per_vertex);
	    std::vector<types::global_dof_index> node_top_center_2_dofs(m_fe.dofs_per_vertex);

	    const double outer_radius = 20.0;
	    const double inner_radius = 12.5;

	    for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)
	      {
		if (   (std::fabs(vertex_itr->vertex()[0] -          0.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] + outer_radius) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] -          0.0) < 1.0e-9) )
		  {
		    node_bottom_center_1_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
		if (   (std::fabs(vertex_itr->vertex()[0] -          0.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] + inner_radius) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] -          0.0) < 1.0e-9) )
		  {
		    node_bottom_center_2_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
		if (   (std::fabs(vertex_itr->vertex()[0] - outer_radius) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] -          0.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] -          0.0) < 1.0e-9) )
		  {
		    node_right_center_1_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
		if (   (std::fabs(vertex_itr->vertex()[0] - inner_radius) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] -          0.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] -          0.0) < 1.0e-9) )
		  {
		    node_right_center_2_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
		if (   (std::fabs(vertex_itr->vertex()[0] + outer_radius) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] -          0.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] -          0.0) < 1.0e-9) )
		  {
		    node_left_center_1_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
		if (   (std::fabs(vertex_itr->vertex()[0] + inner_radius) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] -          0.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] -          0.0) < 1.0e-9) )
		  {
		    node_left_center_2_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
		if (   (std::fabs(vertex_itr->vertex()[0] +          0.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] - outer_radius) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] -          0.0) < 1.0e-9) )
		  {
		    node_top_center_1_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
		if (   (std::fabs(vertex_itr->vertex()[0] +          0.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] - inner_radius) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] -          0.0) < 1.0e-9) )
		  {
		    node_top_center_2_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
	      }

	    m_constraints.add_line(node_bottom_center_1_dofs[0]);
	    m_constraints.set_inhomogeneity(node_bottom_center_1_dofs[0], 0.0);
	    m_constraints.add_line(node_bottom_center_1_dofs[1]);
	    m_constraints.set_inhomogeneity(node_bottom_center_1_dofs[1], 0.0);
	    m_constraints.add_line(node_bottom_center_1_dofs[2]);
	    m_constraints.set_inhomogeneity(node_bottom_center_1_dofs[2], 0.0);

	    m_constraints.add_line(node_bottom_center_2_dofs[0]);
	    m_constraints.set_inhomogeneity(node_bottom_center_2_dofs[0], 0.0);
	    m_constraints.add_line(node_bottom_center_2_dofs[2]);
	    m_constraints.set_inhomogeneity(node_bottom_center_2_dofs[2], 0.0);

	    m_constraints.add_line(node_right_center_1_dofs[2]);
	    m_constraints.set_inhomogeneity(node_right_center_1_dofs[2], 0.0);
	    m_constraints.add_line(node_right_center_2_dofs[2]);
	    m_constraints.set_inhomogeneity(node_right_center_2_dofs[2], 0.0);

	    m_constraints.add_line(node_left_center_1_dofs[2]);
	    m_constraints.set_inhomogeneity(node_left_center_1_dofs[2], 0.0);
	    m_constraints.add_line(node_left_center_2_dofs[2]);
	    m_constraints.set_inhomogeneity(node_left_center_2_dofs[2], 0.0);

	    m_constraints.add_line(node_top_center_1_dofs[2]);
	    m_constraints.set_inhomogeneity(node_top_center_1_dofs[2], 0.0);
	    m_constraints.add_line(node_top_center_2_dofs[2]);
	    m_constraints.set_inhomogeneity(node_top_center_2_dofs[2], 0.0);

	    m_constraints.close();

	  }
	else if (m_parameters.m_scenario == 7)
	  {
	    Triangulation<3>::active_vertex_iterator vertex_itr;
	    vertex_itr = m_triangulation.begin_active_vertex();
	    std::vector<types::global_dof_index> node_bottom_left_dofs(m_fe.dofs_per_vertex);
	    std::vector<types::global_dof_index> node_bottom_right_dofs(m_fe.dofs_per_vertex);
	    std::vector<types::global_dof_index> node_top_left_dofs(m_fe.dofs_per_vertex);
	    std::vector<types::global_dof_index> node_top_right_dofs(m_fe.dofs_per_vertex);

	    std::vector<types::global_dof_index> node_left_bottom_dofs(m_fe.dofs_per_vertex);
	    std::vector<types::global_dof_index> node_left_top_dofs(m_fe.dofs_per_vertex);
	    std::vector<types::global_dof_index> node_right_bottom_dofs(m_fe.dofs_per_vertex);
	    std::vector<types::global_dof_index> node_right_top_dofs(m_fe.dofs_per_vertex);

	    std::vector<types::global_dof_index> node_center_dofs(m_fe.dofs_per_vertex);
	    std::vector<types::global_dof_index> node_center_top_dofs(m_fe.dofs_per_vertex);

	    for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)
	      {
		if (   (std::fabs(vertex_itr->vertex()[0] +  3.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] + 27.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] -  0.0) < 1.0e-9) )
		  {
		    node_bottom_left_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
		if (   (std::fabs(vertex_itr->vertex()[0] -  3.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] + 27.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] -  0.0) < 1.0e-9) )
		  {
		    node_bottom_right_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
		if (   (std::fabs(vertex_itr->vertex()[0] +  3.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] - 27.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] -  0.0) < 1.0e-9) )
		  {
		    node_top_left_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
		if (   (std::fabs(vertex_itr->vertex()[0] -  3.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] - 27.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] -  0.0) < 1.0e-9) )
		  {
		    node_top_right_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
		if (   (std::fabs(vertex_itr->vertex()[0] + 27.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] +  3.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] -  0.0) < 1.0e-9) )
		  {
		    node_left_bottom_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
		if (   (std::fabs(vertex_itr->vertex()[0] + 27.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] -  3.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] -  0.0) < 1.0e-9) )
		  {
		    node_left_top_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
		if (   (std::fabs(vertex_itr->vertex()[0] - 27.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] +  3.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] -  0.0) < 1.0e-9) )
		  {
		    node_right_bottom_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
		if (   (std::fabs(vertex_itr->vertex()[0] - 27.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] -  3.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] -  0.0) < 1.0e-9) )
		  {
		    node_right_top_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
		if (   (std::fabs(vertex_itr->vertex()[0] -  0.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] -  0.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] -  0.0) < 1.0e-9) )
		  {
		    node_center_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
		if (   (std::fabs(vertex_itr->vertex()[0] -  0.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[1] -  3.0) < 1.0e-9)
		    && (std::fabs(vertex_itr->vertex()[2] -  0.0) < 1.0e-9) )
		  {
		    node_center_top_dofs = get_vertex_dofs(vertex_itr, m_dof_handler);
		  }
	      }

	    m_constraints.add_line(node_bottom_left_dofs[2]);
	    m_constraints.set_inhomogeneity(node_bottom_left_dofs[2], 0.0);
	    m_constraints.add_line(node_bottom_right_dofs[2]);
	    m_constraints.set_inhomogeneity(node_bottom_right_dofs[2], 0.0);

	    m_constraints.add_line(node_top_left_dofs[2]);
	    m_constraints.set_inhomogeneity(node_top_left_dofs[2], 0.0);
	    m_constraints.add_line(node_top_right_dofs[2]);
	    m_constraints.set_inhomogeneity(node_top_right_dofs[2], 0.0);

	    m_constraints.add_line(node_left_bottom_dofs[2]);
	    m_constraints.set_inhomogeneity(node_left_bottom_dofs[2], 0.0);
	    m_constraints.add_line(node_left_top_dofs[2]);
	    m_constraints.set_inhomogeneity(node_left_top_dofs[2], 0.0);

	    m_constraints.add_line(node_right_bottom_dofs[2]);
	    m_constraints.set_inhomogeneity(node_right_bottom_dofs[2], 0.0);
	    m_constraints.add_line(node_right_top_dofs[2]);
	    m_constraints.set_inhomogeneity(node_right_top_dofs[2], 0.0);

	    m_constraints.add_line(node_center_dofs[0]);
	    m_constraints.set_inhomogeneity(node_center_dofs[0], 0.0);
	    m_constraints.add_line(node_center_dofs[1]);
	    m_constraints.set_inhomogeneity(node_center_dofs[1], 0.0);

	    m_constraints.add_line(node_center_top_dofs[0]);
	    m_constraints.set_inhomogeneity(node_center_top_dofs[0], 0.0);

	    m_constraints.close();
	  }
	else
	  Assert(false, ExcMessage("The scenario has not been implemented!"));


        // Dirichlet B.C. 2 (displacement loading)
        /*
        {
          const int boundary_id = 7;
          const FEValuesExtractors::Scalar z_displacement(2);
          const double time_inc = m_time.get_delta_t();
          const double disp_magnitude = -1.0;

          // For displacement control, prescribe displacement increment
          VectorTools::interpolate_boundary_values(m_dof_handler,
						   boundary_id,
						   Functions::ConstantFunction<dim>(
						       disp_magnitude*time_inc, dim),
						   m_constraints,
						   m_fe.component_mask(z_displacement));

          const FEValuesExtractors::Scalar x_displacement(0);
          const FEValuesExtractors::Scalar y_displacement(1);
          VectorTools::interpolate_boundary_values(m_dof_handler,
						   boundary_id,
						   Functions::ZeroFunction<dim>(dim),
						   m_constraints,
						   m_fe.component_mask(x_displacement) |
						   m_fe.component_mask(y_displacement));
        }
        */
      }
    else  // inhomogeneous constraints
      {
        if (m_constraints.has_inhomogeneities())
          {
            AffineConstraints<double> homogeneous_constraints(m_constraints);
            for (unsigned int dof = 0; dof != m_dof_handler.n_dofs(); ++dof)
              if (homogeneous_constraints.is_inhomogeneously_constrained(dof))
                homogeneous_constraints.set_inhomogeneity(dof, 0.0);

            m_constraints.clear();
            m_constraints.copy_from(homogeneous_constraints);
          }
      }
    m_constraints.close();
  }



  template <int dim>
  void Solid<dim>::assemble_sc()
  {
    m_timer.enter_subsection("Perform static condensation");
    std::cout << " ASM_SC " << std::flush;

    PerTaskData_SC per_task_data(m_dofs_per_cell,
                                 m_element_indices_u.size(),
                                 m_element_indices_p.size(),
                                 m_element_indices_J.size());
    ScratchData_SC scratch_data;

    WorkStream::run(m_dof_handler.active_cell_iterators(),
                    *this,
                    &Solid::assemble_sc_one_cell,
                    &Solid::copy_local_to_global_sc,
                    scratch_data,
                    per_task_data);

    m_timer.leave_subsection();
  }


  template <int dim>
  void Solid<dim>::copy_local_to_global_sc(const PerTaskData_SC &data)
  {
    for (unsigned int i = 0; i < m_dofs_per_cell; ++i)
      for (unsigned int j = 0; j < m_dofs_per_cell; ++j)
        m_tangent_matrix.add(data.m_local_dof_indices[i],
                           data.m_local_dof_indices[j],
                           data.m_cell_matrix(i, j));
  }


  template <int dim>
  void Solid<dim>::assemble_sc_one_cell(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    ScratchData_SC &                                      scratch,
    PerTaskData_SC &                                      data)
  {
    data.reset();
    scratch.reset();
    cell->get_dof_indices(data.m_local_dof_indices);


    data.m_k_orig.extract_submatrix_from(m_tangent_matrix,
                                       data.m_local_dof_indices,
                                       data.m_local_dof_indices);
    data.m_k_pu.extract_submatrix_from(data.m_k_orig,
                                     m_element_indices_p,
                                     m_element_indices_u);
    data.m_k_pJ.extract_submatrix_from(data.m_k_orig,
                                     m_element_indices_p,
                                     m_element_indices_J);
    data.m_k_JJ.extract_submatrix_from(data.m_k_orig,
                                     m_element_indices_J,
                                     m_element_indices_J);

    data.m_k_pJ_inv.invert(data.m_k_pJ);

    data.m_k_pJ_inv.mmult(data.m_A, data.m_k_pu);
    data.m_k_JJ.mmult(data.m_B, data.m_A);
    data.m_k_pJ_inv.Tmmult(data.m_C, data.m_B);
    data.m_k_pu.Tmmult(data.m_k_bbar, data.m_C);
    data.m_k_bbar.scatter_matrix_to(m_element_indices_u,
                                  m_element_indices_u,
                                  data.m_cell_matrix);

    data.m_k_pJ_inv.add(-1.0, data.m_k_pJ);
    data.m_k_pJ_inv.scatter_matrix_to(m_element_indices_p,
                                    m_element_indices_J,
                                    data.m_cell_matrix);
  }

  template <int dim>
  std::pair<unsigned int, double>
  Solid<dim>::solve_linear_system(BlockVector<double> &newton_update)
  {
    unsigned int lin_it  = 0;
    double       lin_res = 0.0;

    if (m_parameters.m_use_static_condensation == true)
      {

        BlockVector<double> A(m_dofs_per_block);
        BlockVector<double> B(m_dofs_per_block);


        {
          assemble_sc();

          m_tangent_matrix.block(m_p_dof, m_J_dof)
            .vmult(A.block(m_J_dof), m_system_rhs.block(m_p_dof));
          m_tangent_matrix.block(m_J_dof, m_J_dof)
            .vmult(B.block(m_J_dof), A.block(m_J_dof));
          A.block(m_J_dof) = m_system_rhs.block(m_J_dof);
          A.block(m_J_dof) -= B.block(m_J_dof);
          m_tangent_matrix.block(m_p_dof, m_J_dof)
            .Tvmult(A.block(m_p_dof), A.block(m_J_dof));
          m_tangent_matrix.block(m_u_dof, m_p_dof)
            .vmult(A.block(m_u_dof), A.block(m_p_dof));
          m_system_rhs.block(m_u_dof) -= A.block(m_u_dof);

          m_timer.enter_subsection("Linear solver");
          std::cout << " SLV " << std::flush;
          if (m_parameters.m_type_lin == "CG")
            {
              const auto solver_its = static_cast<unsigned int>(
                m_tangent_matrix.block(m_u_dof, m_u_dof).m() *
                m_parameters.m_max_iterations_lin);
              const double tol_sol =
                m_parameters.m_tol_lin * m_system_rhs.block(m_u_dof).l2_norm();

              SolverControl solver_control(solver_its, tol_sol);

              GrowingVectorMemory<Vector<double>> GVM;
              SolverCG<Vector<double>> solver_CG(solver_control, GVM);

              PreconditionSelector<SparseMatrix<double>, Vector<double>>
                preconditioner(m_parameters.m_preconditioner_type,
                               m_parameters.m_preconditioner_relaxation);
              preconditioner.use_matrix(m_tangent_matrix.block(m_u_dof, m_u_dof));

              solver_CG.solve(m_tangent_matrix.block(m_u_dof, m_u_dof),
                              newton_update.block(m_u_dof),
                              m_system_rhs.block(m_u_dof),
                              preconditioner);

              lin_it  = solver_control.last_step();
              lin_res = solver_control.last_value();
            }
          else if (m_parameters.m_type_lin == "Direct")
            {
              SparseDirectUMFPACK A_direct;
              A_direct.initialize(m_tangent_matrix.block(m_u_dof, m_u_dof));
              A_direct.vmult(newton_update.block(m_u_dof),
                             m_system_rhs.block(m_u_dof));

              lin_it  = 1;
              lin_res = 0.0;
            }
          else
            Assert(false, ExcMessage("Linear solver type not implemented"));

          m_timer.leave_subsection();
        }

        m_constraints.distribute(newton_update);

        m_timer.enter_subsection("Linear solver postprocessing");
        std::cout << " PP " << std::flush;

        {
          m_tangent_matrix.block(m_p_dof, m_u_dof)
            .vmult(A.block(m_p_dof), newton_update.block(m_u_dof));
          A.block(m_p_dof) *= -1.0;
          A.block(m_p_dof) += m_system_rhs.block(m_p_dof);
          m_tangent_matrix.block(m_p_dof, m_J_dof)
            .vmult(newton_update.block(m_J_dof), A.block(m_p_dof));
        }

        m_constraints.distribute(newton_update);

        {
          m_tangent_matrix.block(m_J_dof, m_J_dof)
            .vmult(A.block(m_J_dof), newton_update.block(m_J_dof));
          A.block(m_J_dof) *= -1.0;
          A.block(m_J_dof) += m_system_rhs.block(m_J_dof);
          m_tangent_matrix.block(m_p_dof, m_J_dof)
            .Tvmult(newton_update.block(m_p_dof), A.block(m_J_dof));
        }

        m_constraints.distribute(newton_update);

        m_timer.leave_subsection();
      }
    else
      {
        std::cout << " ------ " << std::flush;

        m_timer.enter_subsection("Linear solver");
        std::cout << " SLV " << std::flush;

        if (m_parameters.m_type_lin == "CG")
          {

            const Vector<double> &f_u = m_system_rhs.block(m_u_dof);
            const Vector<double> &f_p = m_system_rhs.block(m_p_dof);
            const Vector<double> &f_J = m_system_rhs.block(m_J_dof);

            Vector<double> &d_u = newton_update.block(m_u_dof);
            Vector<double> &d_p = newton_update.block(m_p_dof);
            Vector<double> &d_J = newton_update.block(m_J_dof);

            const auto K_uu =
              linear_operator(m_tangent_matrix.block(m_u_dof, m_u_dof));
            const auto K_up =
              linear_operator(m_tangent_matrix.block(m_u_dof, m_p_dof));
            const auto K_pu =
              linear_operator(m_tangent_matrix.block(m_p_dof, m_u_dof));
            const auto K_Jp =
              linear_operator(m_tangent_matrix.block(m_J_dof, m_p_dof));
            const auto K_JJ =
              linear_operator(m_tangent_matrix.block(m_J_dof, m_J_dof));

            PreconditionSelector<SparseMatrix<double>, Vector<double>>
              preconditioner_K_Jp_inv("jacobi");
            preconditioner_K_Jp_inv.use_matrix(
              m_tangent_matrix.block(m_J_dof, m_p_dof));
            ReductionControl solver_control_K_Jp_inv(
              static_cast<unsigned int>(m_tangent_matrix.block(m_J_dof, m_p_dof).m() *
                                        m_parameters.m_max_iterations_lin),
              1.0e-30,
              m_parameters.m_tol_lin);
            SolverSelector<Vector<double>> solver_K_Jp_inv;
            solver_K_Jp_inv.select("cg");
            solver_K_Jp_inv.set_control(solver_control_K_Jp_inv);
            const auto K_Jp_inv =
              inverse_operator(K_Jp, solver_K_Jp_inv, preconditioner_K_Jp_inv);

            const auto K_pJ_inv     = transpose_operator(K_Jp_inv);
            const auto K_pp_bar     = K_Jp_inv * K_JJ * K_pJ_inv;
            const auto K_uu_bar_bar = K_up * K_pp_bar * K_pu;
            const auto K_uu_con     = K_uu + K_uu_bar_bar;

            PreconditionSelector<SparseMatrix<double>, Vector<double>>
              preconditioner_K_con_inv(m_parameters.m_preconditioner_type,
                                       m_parameters.m_preconditioner_relaxation);
            preconditioner_K_con_inv.use_matrix(
              m_tangent_matrix.block(m_u_dof, m_u_dof));
            ReductionControl solver_control_K_con_inv(
              static_cast<unsigned int>(m_tangent_matrix.block(m_u_dof, m_u_dof).m() *
                                        m_parameters.m_max_iterations_lin),
              1.0e-30,
              m_parameters.m_tol_lin);
            SolverSelector<Vector<double>> solver_K_con_inv;
            solver_K_con_inv.select("cg");
            solver_K_con_inv.set_control(solver_control_K_con_inv);
            const auto K_uu_con_inv =
              inverse_operator(K_uu_con,
                               solver_K_con_inv,
                               preconditioner_K_con_inv);

            d_u =
              K_uu_con_inv * (f_u - K_up * (K_Jp_inv * f_J - K_pp_bar * f_p));

            m_timer.leave_subsection();

            m_timer.enter_subsection("Linear solver postprocessing");
            std::cout << " PP " << std::flush;

            d_J = K_pJ_inv * (f_p - K_pu * d_u);
            d_p = K_Jp_inv * (f_J - K_JJ * d_J);

            lin_it  = solver_control_K_con_inv.last_step();
            lin_res = solver_control_K_con_inv.last_value();
          }
        else if (m_parameters.m_type_lin == "Direct")
          {
            SparseDirectUMFPACK A_direct;
            A_direct.initialize(m_tangent_matrix);
            A_direct.vmult(newton_update, m_system_rhs);

            lin_it  = 1;
            lin_res = 0.0;

            std::cout << " -- " << std::flush;
          }
        else
          Assert(false, ExcMessage("Linear solver type not implemented"));

        m_timer.leave_subsection();

        m_constraints.distribute(newton_update);
      }

    return std::make_pair(lin_it, lin_res);
  }

  template <int dim>
  void Solid<dim>::output_results() const
  {
    DataOut<dim> data_out;
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(
      DataComponentInterpretation::component_is_scalar);
    data_component_interpretation.push_back(
      DataComponentInterpretation::component_is_scalar);

    std::vector<std::string> solution_name(dim, "displacement");
    solution_name.emplace_back("pressure");
    solution_name.emplace_back("dilatation");

    DataOutBase::VtkFlags output_flags;
    output_flags.write_higher_order_cells       = true;
    output_flags.physical_units["displacement"] = "m";
    data_out.set_flags(output_flags);

    data_out.attach_dof_handler(m_dof_handler);
    data_out.add_data_vector(m_solution_n,
                             solution_name,
                             DataOut<dim>::type_dof_data,
                             data_component_interpretation);

    Vector<double> cell_material_id(m_triangulation.n_active_cells());
    // output material ID for each cell
    for (const auto &cell : m_dof_handler.active_cell_iterators())
      {
	cell_material_id(cell->active_cell_index()) = cell->material_id();
      }
    data_out.add_data_vector(cell_material_id, "materialID");

    Vector<double> cell_total_strain_energy(m_triangulation.n_active_cells());
    Vector<double> cell_mechanical_strain_energy(m_triangulation.n_active_cells());
    Vector<double> cell_volumetric_strain_energy(m_triangulation.n_active_cells());
    Vector<double> cell_isotropic_strain_energy(m_triangulation.n_active_cells());
    Vector<double> cell_magnetic_strain_energy(m_triangulation.n_active_cells());
    FEValues<dim> fe_values(m_fe, m_qf_cell, update_JxW_values);

    for (const auto &cell : m_dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);

        const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
          m_quadrature_point_history.get_data(cell);
        Assert(lqph.size() == m_n_q_points, ExcInternalError());

        double element_total_strain_energy = 0.0;
        double element_mechanical_strain_energy = 0.0;
        double element_volumetric_strain_energy = 0.0;
        double element_isotropic_strain_energy = 0.0;
        double element_magnetic_strain_energy = 0.0;
        double element_volume = 0.0;

        for (unsigned int q_point = 0; q_point < m_n_q_points; ++q_point)
          {
            const double det_F_qp = lqph[q_point]->get_det_F();
            const double JxW      = fe_values.JxW(q_point);

            element_total_strain_energy += lqph[q_point]->get_total_strain_energy()
        	                             * det_F_qp * JxW;
            element_mechanical_strain_energy += lqph[q_point]->get_mechanical_strain_energy()
                                                  * det_F_qp * JxW;
            element_volumetric_strain_energy += lqph[q_point]->get_volumetric_strain_energy()
        	                                  * det_F_qp * JxW;
            element_isotropic_strain_energy += lqph[q_point]->get_isotropic_strain_energy()
        	                                 * det_F_qp * JxW;
            element_magnetic_strain_energy += lqph[q_point]->get_magnetic_strain_energy()
        	                                * det_F_qp * JxW;
            element_volume += det_F_qp * JxW;
          }
        cell_total_strain_energy(cell->active_cell_index()) =
            element_total_strain_energy/element_volume;

        cell_mechanical_strain_energy(cell->active_cell_index()) =
            element_mechanical_strain_energy/element_volume;

        cell_volumetric_strain_energy(cell->active_cell_index()) =
            element_volumetric_strain_energy/element_volume;

        cell_isotropic_strain_energy(cell->active_cell_index()) =
            element_isotropic_strain_energy/element_volume;

        cell_magnetic_strain_energy(cell->active_cell_index()) =
            element_magnetic_strain_energy/element_volume;
      }
    data_out.add_data_vector(cell_total_strain_energy, "TotalStrainEnergy");
    data_out.add_data_vector(cell_mechanical_strain_energy, "MechanicalStrainEnergy");
    data_out.add_data_vector(cell_volumetric_strain_energy, "VolumetricStrainEnergy");
    data_out.add_data_vector(cell_isotropic_strain_energy, "IsotropicStrainEnergy");
    data_out.add_data_vector(cell_magnetic_strain_energy, "MagneticStrainEnergy");

    // MappingQEulerian maps the original mesh to the deformed mesh
    Vector<double> soln(m_solution_n.size());
    for (unsigned int i = 0; i < soln.size(); ++i)
      soln(i) = m_solution_n(i);
    MappingQEulerian<dim> q_mapping(m_degree, m_dof_handler, soln);
    data_out.build_patches(q_mapping, m_degree);

    std::ofstream output("solution-" + std::to_string(dim) + "d-" +
                         std::to_string(m_time.get_timestep()) + ".vtu");
    data_out.write_vtu(output);
  }

  template <int dim>
  void Solid<dim>::read_time_data(const std::string &data_file,
				  std::vector<std::array<double, 3>> & time_table) const
  {
    std::ifstream myfile (data_file);

    double t_0, t_1, delta_t;

    if (myfile.is_open())
      {
	std::cout << "Reading time data file ..." << std::endl;

	while ( myfile >> t_0
		       >> t_1
		       >> delta_t)
	  {
	    Assert( t_0 < t_1,
		    ExcMessage("For each time pair, "
			       "the start time should be smaller than the end time"));
	    time_table.push_back({{t_0, t_1, delta_t}});
	  }

	Assert(std::fabs(t_1 - m_parameters.m_end_time) < 1.0e-9,
	       ExcMessage("End time in time table is inconsistent with input data in parameters.prm"))

	Assert(time_table.size() > 0,
	       ExcMessage("Time data file is empty."));
	myfile.close();
      }
    else
      {
        std::cout << "Time data file : " << data_file << " not exist!" << std::endl;
        Assert(false, ExcMessage("Failed to read time data file"));
      }

    for (auto & time_group : time_table)
      {
	std::cout << time_group[0] << ",\t"
	          << time_group[1] << ",\t"
		  << time_group[2] << std::endl;
      }
  }


  template <int dim>
  void Solid<dim>::read_material_data(const std::string &data_file,
				      const unsigned int total_material_regions)
  {
    std::ifstream myfile (data_file);

    double shear_modulus, bulk_modulus, Br_tilde_x, Br_tilde_y, Br_tilde_z;
    int material_region;
    double poisson_ratio;
    if (myfile.is_open())
      {
        std::cout << "Reading material data file ..." << std::endl;

        while ( myfile >> material_region
                       >> shear_modulus
		       >> bulk_modulus
		       >> Br_tilde_x
		       >> Br_tilde_y
		       >> Br_tilde_z )
          {
            m_material_data[material_region] = {shear_modulus,
        	                                bulk_modulus,
						Br_tilde_x,
						Br_tilde_y,
						Br_tilde_z};
            poisson_ratio = (3*bulk_modulus-2*shear_modulus) / (2*(3*bulk_modulus+shear_modulus));
            Assert( (poisson_ratio <= 0.5)&(poisson_ratio >=-1.0) , ExcInternalError());

            std::cout << "Region " << material_region << " : " << std::endl;
            std::cout << "    shear modulus = " << shear_modulus << std::endl;
            std::cout << "    bulk modulus = "  << bulk_modulus << std::endl;
            std::cout << "    poisson ratio = "  << poisson_ratio << std::endl;
            std::cout << "    Residual magnetic flux density in reference configuration = "
        	      << "(" << Br_tilde_x << ", "
		             << Br_tilde_y << ", "
			     << Br_tilde_z
		      << ")" << std::endl;
          }

        if (m_material_data.size() != total_material_regions)
          {
            std::cout << "Material data file has " << m_material_data.size() << " rows. However, "
        	      << "the mesh has " << total_material_regions << " material regions."
		      << std::endl;
            Assert(m_material_data.size() == total_material_regions,
                       ExcDimensionMismatch(m_material_data.size(), total_material_regions));
          }
        myfile.close();
      }
    else
      {
	std::cout << "Material data file : " << data_file << " not exist!" << std::endl;
	Assert(false, ExcMessage("Failed to read material data file"));
      }
  }


} // namespace FiniteDeformationMagnetoElasticityThreeField


int main()
{
  using namespace FiniteDeformationMagnetoElasticityThreeField;

  try
    {
      const unsigned int dim = 3;
      Solid<dim>         solid("parameters.prm");
      solid.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
