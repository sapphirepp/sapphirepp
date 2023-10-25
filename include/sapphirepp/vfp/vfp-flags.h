#ifndef VFP_VFPFLAGS_H
#define VFP_VFPFLAGS_H

namespace Sapphire
{
  namespace VFP
  {
    enum class TermFlags
    {
      none              = 0,
      spatial_advection = 1 << 0,
      collision         = 1 << 1,
      magnetic          = 1 << 2,
      momentum          = 1 << 3,
      source            = 1 << 4
    };

    constexpr TermFlags
    operator|(TermFlags f1, TermFlags f2)
    {
      return static_cast<TermFlags>(static_cast<int>(f1) |
                                    static_cast<int>(f2));
    }

    constexpr TermFlags
    operator&(TermFlags f1, TermFlags f2)
    {
      return static_cast<TermFlags>(static_cast<int>(f1) &
                                    static_cast<int>(f2));
    }

    template <typename StreamType>
    inline StreamType &
    operator<<(StreamType &os, TermFlags f)
    {
      os << "Term flags: \n";
      if ((f & TermFlags::spatial_advection) != TermFlags::none)
        os << "	 - Spatial Advection\n";
      if ((f & TermFlags::collision) != TermFlags::none)
        os << "	 - Collision\n";
      if ((f & TermFlags::magnetic) != TermFlags::none)
        os << "	 - Magnetic\n";
      if ((f & TermFlags::momentum) != TermFlags::none)
        os << "	 - Momentum\n";
      if ((f & TermFlags::source) != TermFlags::none)
        os << "	 - Source\n";
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
