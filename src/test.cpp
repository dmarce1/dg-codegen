#include <gtest/gtest.h>

#include "dgTransforms.hpp"

bool testTransforms() {
	std::array<double, 3> x0, x1, x2;
	x0[0] = 1;
	x0[1] = 1;
	x0[2] = 1;
	x1[0] = -sqrt(0.6);
	x1[1] = 0.0;
	x1[2] = sqrt(0.6);
	x2[0] = 0.4;
	x2[1] = -0.5;
	x2[2] = 0.4;
	auto y0 = dgMassInverse<double, 1, 3>(dgStiffness<double, 1, 3>(0, dgMassInverse<double, 1, 3>(dgAnalyze<double, 1, 3>(x0))));
	auto y1 = dgMassInverse<double, 1, 3>(dgStiffness<double, 1, 3>(0, dgMassInverse<double, 1, 3>(dgAnalyze<double, 1, 3>(x1))));
	auto y2 = dgMassInverse<double, 1, 3>(dgStiffness<double, 1, 3>(0, dgMassInverse<double, 1, 3>(dgAnalyze<double, 1, 3>(x2))));
	printf( "%e %e %e \n", y0[0], y0[1], y0[2]);
	printf( "%e %e %e \n", y1[0], y1[1], y1[2]);
	printf( "%e %e %e \n", y2[0], y2[1], y2[2]);

	return true;
}

TEST(UnitTests, Legendre) {
    EXPECT_EQ(testTransforms(), true);
}
