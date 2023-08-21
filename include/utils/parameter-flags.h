#ifndef UTILS_PARAMETERFLAGS_H
#define UTILS_PARAMETERFLAGS_H

namespace Sapphire
{
  namespace Utils
  {
    enum class OutputFormat
    {
      vtu,
      pvtu,
      hdf5
    };
    enum class TimeSteppingMethod
    {
      crank_nicolson,
      forward_euler,
      backward_euler,
      erk4,
      lserk
    };
  } // namespace Utils

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
  } // namespace VFP
} // namespace Sapphire
#endif
