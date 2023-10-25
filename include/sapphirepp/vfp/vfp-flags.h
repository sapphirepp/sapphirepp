/**
 * @file vfp-flags.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define VFPFlags and other enums used in the VFP module
 * @version 0.1
 * @date 2023-10-25
 *
 * @copyright Copyright (c) 2023
 *
 * This file defines enums used in the VFP module.
 * The VFPFlags enum is used to specify which terms should be included in the
 * VFP equation, wether they are time dependent, and if they are momentum
 * depdentent.
 */

#ifndef VFP_VFPFLAGS_H
#define VFP_VFPFLAGS_H

namespace Sapphire
{
  namespace VFP
  {
    enum class VFPFlags
    {
      none = 0,

      /** Do we have a spatical advection term,
       * \f$ (\mathbf{u} + \mathbf{v}) \cdot \nabla_x f \f$ ? */
      spatial_advection = 1 << 0,

      /** Do we have a collision term,
       * \f$ \frac{\nu}{2} \Delta_{\theta, \varphi} f \f$ ? */
      collision = 1 << 1,

      /** Do we have a magnectic field,
       * \f$ q \mathbf{v} \cdot \left( \mathbf{B} \times \nabla_{p} f \right)
       * \f$ ? */
      magnetic = 1 << 2,

      /** The fields are time dependent,
       * \f$ \mathbf{u}(\mathbf{x}, \mathbf{p}, t),
       * \mathbf{B}(\mathbf{x}, \mathbf{p}, t) \f$ . */
      time_dependent_fields = 0 << 3,

      /** The fields are time independnet,
       * \f$ \mathbf{u}(\mathbf{x}, \mathbf{p}),
       * \mathbf{B}(\mathbf{x}, \mathbf{p}) \f$ . */
      time_independent_fields = 1 << 3,

      /** Do we have a momentum term
       * \f$ \left( \gamma m \frac{\mathrm{D} \mathbf{u}}{\mathrm{D} t} +
       * \mathbf{p} \cdot\nabla_{x} \mathbf{u} \right) \cdot \nabla_{p} f \f$ ?
       */
      momentum = 1 << 4,

      /** Use a logarithmic momentum variable, \f$ \ln p \f$ . */
      logarithmic_p = 0 << 5,

      /** Use a linear momentum variable, \f$ p \f$ ? */
      linear_p = 1 << 5,

      /** Do we have a source term,
       * \f$ S(\mathbf{x}, \mathbf{p}, t) \f$ ? */
      source = 1 << 6,

      /** The source term is time dependent,
       * \f$ S(\mathbf{x}, \mathbf{p}, t) \f$ . */
      time_dependent_source = 0 << 7,

      /** The source term is time independent,
       * \f$ S(\mathbf{x}, \mathbf{p}) \f$ . */
      time_independent_source = 1 << 7
    };

    constexpr VFPFlags
    operator|(VFPFlags f1, VFPFlags f2)
    {
      return static_cast<VFPFlags>(static_cast<int>(f1) | static_cast<int>(f2));
    }

    constexpr VFPFlags
    operator&(VFPFlags f1, VFPFlags f2)
    {
      return static_cast<VFPFlags>(static_cast<int>(f1) & static_cast<int>(f2));
    }

    template <typename StreamType>
    inline StreamType &
    operator<<(StreamType &os, VFPFlags f)
    {
      os << "VFP flags: \n";
      if ((f & VFPFlags::spatial_advection) != VFPFlags::none)
        os << "	 - Spatial Advection\n";
      if ((f & VFPFlags::collision) != VFPFlags::none)
        os << "	 - Collision\n";
      if ((f & VFPFlags::magnetic) != VFPFlags::none)
        os << "	 - Magnetic\n";
      if ((f & VFPFlags::time_independent_fields) != VFPFlags::none)
        os << "	 - Time Independent Fields\n";
      else
        os << "	 - Time Dependent Fields\n";
      if ((f & VFPFlags::momentum) != VFPFlags::none)
        {
          os << "	 - Momentum";
          if ((f & VFPFlags::linear_p) != VFPFlags::none)
            os << "	(linear)\n";
          else
            os << "	(logarithmic)\n";
        }
      if ((f & VFPFlags::source) != VFPFlags::none)
        {
          os << "	 - Source";
          if ((f & VFPFlags::time_independent_source) != VFPFlags::none)
            os << "	(time independent)\n";
          else
            os << " (time dependent)\n";
        }
      return os;
    }

    enum class BoundaryConditions
    {
      continous_gradients,
      zero_inflow,
      periodic
    };

    enum class TimeSteppingMethod
    {
      crank_nicolson,
      forward_euler,
      backward_euler,
      erk4,
      lserk
    };

    enum class GridType
    {
      hypercube,
      file
    };
  } // namespace VFP
} // namespace Sapphire
#endif
