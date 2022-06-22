#include <gtest/gtest.h>

#include <fstream>

#include "poisson.h"

using namespace dealii;

// Test Fixture for Poisson problem
template <class Integral>
class PoissonTester : public ::testing::Test, public Poisson<Integral::value>
{
public:
  PoissonTester() = default;
};

using PoissonTestTypes = ::testing::Types<std::integral_constant<int, 1>,
                                          std::integral_constant<int, 2>,
                                          std::integral_constant<int, 3>>;


using Poisson2DTester = PoissonTester<std::integral_constant<int, 2>>;

TYPED_TEST_CASE(PoissonTester, PoissonTestTypes);

TYPED_TEST(PoissonTester, MakeGrid)
{
  // Output dimension
  std::cout << "Working on dim " << TypeParam::value << std::endl;
  this->make_grid();
}