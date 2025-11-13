#include <GlymphBIE.hpp>

#include "unit_test_framework.hpp"

template <class Real>
TEST_CASE(Test1)
{

}

template <class Real>
TEST_SUITE(TestSuite1)
{
    TEST(Test1<Real>);
}


auto
main() -> int
{
  // Run the unit tests. If a test fails, the program will print failure info
  // and return 1.
  RUN_SUITE(TestSuite1<float>);
  RUN_SUITE(TestSuite1<double>);
  return 0; 
}