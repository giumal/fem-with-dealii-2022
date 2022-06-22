#include "step-3.h"


int
main()
{
  deallog.depth_console(2);
  // ParameterHandler prm;
  Step3 laplace_problem;
  ParameterAcceptor::initialize("poisson.prm");
  laplace_problem.run();
  return 0;
}