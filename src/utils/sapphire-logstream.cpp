#include "sapphire-logstream.h"

#include <iostream>

Sapphire::Utils::SapphireLogStream Sapphire::saplog;

Sapphire::Utils::SapphireLogStream::SapphireLogStream()
  : dealii::LogStream()
{
  pop();
  push("Sapphire");
}