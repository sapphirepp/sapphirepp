/**
 * @file sapphire-logstream.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define saplog
 * @version 0.1
 * @date 2023-07-17
 */

#ifndef SAPPHIREUTISL_SAPPHIRELOGSTREAM_H
#define SAPPHIREUTISL_SAPPHIRELOGSTREAM_H

#include <deal.II/base/logstream.h>

#include <string>

namespace Sapphire
{
  namespace Utils
  {

    class SapphireLogStream : public dealii::LogStream
    {
    public:
      SapphireLogStream();
    };

  } // namespace Utils

  extern Sapphire::Utils::SapphireLogStream saplog;
} // namespace Sapphire
#endif
