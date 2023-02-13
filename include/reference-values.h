#include <deal.II/base/conditional_ostream.h>

#include <ostream>
namespace VFPSolver {
struct ReferenceValues;

std::ostream &operator<<(std::ostream &os,
                         const ReferenceValues &reference_values);

// dealii::ConditionalOStream &operator<<(dealii::ConditionalOStream &pos,
//                                        const ReferenceValues &reference_values);

}  // namespace VFPSolver
