#include <GlymphBIE.hpp>

#include "unit_test_framework.hpp"

#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>

template <class Real> 
TEST_CASE(test1) {
    sctl::Comm comm_;
    StokesBIO<Real> biop(1.0, 1.0, comm_);
}

template <class Real>
TEST_SUITE(TestSuite1) {
    TEST(test1<Real>);
}

template <class Real>
void test(const sctl::Comm comm_) {
  StokesBIO<Real> biop(1.0, 1.0, comm_);
}

int main(int argc, char** argv)
{
  sctl::Comm::MPI_Init(&argc, &argv);
  using Real=double;
  // Run the unit tests. If a test fails, the program will print failure info
  // and return 1.
  // RUN_SUITE(TestSuite1<float>);
  // RUN_SUITE(TestSuite1<double>);
  {
    sctl::Comm comm = sctl::Comm::World();
    test<Real>(comm);
  }

  sctl::Comm::MPI_Finalize();
  return 0; 
}