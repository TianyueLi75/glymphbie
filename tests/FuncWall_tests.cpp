#include <GlymphBIE.hpp>

#include "unit_test_framework.hpp"

// Test constant functions on inner and outer walls
template <class Real>
TEST_CASE(Test1)
{

}

// Test constant functions on inner and sinusoidal on outer walls

// Test sinusoidal function on inner and constant on outer walls


// Test sinusoidal functions on both inner and outer walls

// Test more complex functions on both inner and outer walls

// Test enforcement of minimum gap and minimum inner radius

// Test updating beyond number of time steps

 

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