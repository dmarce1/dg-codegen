#pragma once

#include "Util.hpp"
#include <array>
#include <cmath>

enum class Quadrature : int {
	gaussLegendre, gaussLobatto
};

template<int D, int O>
constexpr int triangleSize = binco(O + D - 1, D);

template<int D, int O>
constexpr int squareSize = ipow(O, D);


template<typename T>
std::array<T, triangleSize<1, 1>> dgAnalyze_1D_O1(std::array<T, squareSize<1, 1>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 1> buffer;
   std::array<T, 1> output;
   static U const C0 = U(2.0000000000000000e+00);
   output[0] = C0 * input[0];
   return output;
}

template<typename T>
std::array<T, squareSize<1, 1>> dgSynthesize_1D_O1(std::array<T, triangleSize<1, 1>> const& input) {
   std::array<T, 1> buffer;
   std::array<T, 1> output;
   output[0] = input[0];
   return output;
}

template<typename T>
std::array<T, triangleSize<1, 1>> dgMassInverse_1D_O1(std::array<T, triangleSize<1, 1>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(5.0000000000000000e-01);
   std::array<T, triangleSize<1, 1>> output;
   output[0] = C0 * input[0];
   return output;
}

template<typename T>
std::array<T, triangleSize<1, 1>> dgStiffness_1D_O1(int dimension, std::array<T, triangleSize<1, 1>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, triangleSize<1, 1>> output;
   output[0] = U(0.0);
   return output;
}

template<typename T>
std::array<T, triangleSize<0, 1>> dgTrace_1D_O1(int face, std::array<T, triangleSize<1, 1>> const& input) {
   std::array<T, triangleSize<0, 1>> output;
   if(face == 0) {
      output[0] = input[0];
   } else /*if(face == 1)*/ {
      output[0] = input[0];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<1, 1>> dgTraceInverse_1D_O1(int face, std::array<T, triangleSize<0, 1>> const& input) {
   std::array<T, triangleSize<1, 1>> output;
   if(face == 0) {
      output[0] = input[0];
   } else /*if(face == 1)*/ {
      output[0] = input[0];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<1, 2>> dgAnalyze_1D_O2(std::array<T, squareSize<1, 2>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 2> buffer;
   std::array<T, 2> output;
   static U const C0 = U(5.7735026918962573e-01);
   output[0] = input[0] + input[1];
   output[1] = C0 * (-input[0] + input[1]);
   return output;
}

template<typename T>
std::array<T, squareSize<1, 2>> dgSynthesize_1D_O2(std::array<T, triangleSize<1, 2>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 2> buffer;
   std::array<T, 2> output;
   static U const C0 = U(5.7735026918962573e-01);
   output[0] = -C0 * input[1] + input[0];
   output[1] = C0 * input[1] + input[0];
   return output;
}

template<typename T>
std::array<T, triangleSize<1, 2>> dgMassInverse_1D_O2(std::array<T, triangleSize<1, 2>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(5.0000000000000000e-01);
   static U const C1 = U(1.5000000000000000e+00);
   std::array<T, triangleSize<1, 2>> output;
   output[0] = C0 * input[0];
   output[1] = C1 * input[1];
   return output;
}

template<typename T>
std::array<T, triangleSize<1, 2>> dgStiffness_1D_O2(int dimension, std::array<T, triangleSize<1, 2>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(2.0000000000000000e+00);
   std::array<T, triangleSize<1, 2>> output;
   output[0] = U(0.0);
   output[1] = C0 * input[0];
   return output;
}

template<typename T>
std::array<T, triangleSize<0, 2>> dgTrace_1D_O2(int face, std::array<T, triangleSize<1, 2>> const& input) {
   std::array<T, triangleSize<0, 2>> output;
   if(face == 0) {
      output[0] = input[0] - input[1];
   } else /*if(face == 1)*/ {
      output[0] = input[0] + input[1];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<1, 2>> dgTraceInverse_1D_O2(int face, std::array<T, triangleSize<0, 2>> const& input) {
   std::array<T, triangleSize<1, 2>> output;
   if(face == 0) {
      output[0] = input[0];
      output[1] = -input[0];
   } else /*if(face == 1)*/ {
      output[0] = input[0];
      output[1] = input[0];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<1, 3>> dgAnalyze_1D_O3(std::array<T, squareSize<1, 3>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 3> buffer;
   std::array<T, 3> output;
   static U const C0 = U(5.5555555555555558e-01);
   static U const C1 = U(8.8888888888888884e-01);
   static U const C2 = U(4.3033148291193524e-01);
   static U const C3 = U(2.2222222222222221e-01);
   static U const C4 = U(4.4444444444444442e-01);
   output[0] = C0 * (input[0] + input[2]) + C1 * input[1];
   output[1] = C2 * (-input[0] + input[2]);
   output[2] = C3 * (input[0] + input[2]) - C4 * input[1];
   return output;
}

template<typename T>
std::array<T, squareSize<1, 3>> dgSynthesize_1D_O3(std::array<T, triangleSize<1, 3>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 3> buffer;
   std::array<T, 3> output;
   static U const C0 = U(4.0000000000000002e-01);
   static U const C1 = U(7.7459666924148340e-01);
   static U const C2 = U(5.0000000000000000e-01);
   output[0] = C0 * input[2] - C1 * input[1] + input[0];
   output[1] = -C2 * input[2] + input[0];
   output[2] = C0 * input[2] + C1 * input[1] + input[0];
   return output;
}

template<typename T>
std::array<T, triangleSize<1, 3>> dgMassInverse_1D_O3(std::array<T, triangleSize<1, 3>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(5.0000000000000000e-01);
   static U const C1 = U(1.5000000000000000e+00);
   static U const C2 = U(2.5000000000000000e+00);
   std::array<T, triangleSize<1, 3>> output;
   output[0] = C0 * input[0];
   output[1] = C1 * input[1];
   output[2] = C2 * input[2];
   return output;
}

template<typename T>
std::array<T, triangleSize<1, 3>> dgStiffness_1D_O3(int dimension, std::array<T, triangleSize<1, 3>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(2.0000000000000000e+00);
   std::array<T, triangleSize<1, 3>> output;
   output[0] = U(0.0);
   output[1] = C0 * input[0];
   output[2] = C0 * input[1];
   return output;
}

template<typename T>
std::array<T, triangleSize<0, 3>> dgTrace_1D_O3(int face, std::array<T, triangleSize<1, 3>> const& input) {
   std::array<T, triangleSize<0, 3>> output;
   if(face == 0) {
      output[0] = input[0] - input[1] + input[2];
   } else /*if(face == 1)*/ {
      output[0] = input[0] + input[1] + input[2];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<1, 3>> dgTraceInverse_1D_O3(int face, std::array<T, triangleSize<0, 3>> const& input) {
   std::array<T, triangleSize<1, 3>> output;
   if(face == 0) {
      output[0] = input[0];
      output[1] = -input[0];
      output[2] = input[0];
   } else /*if(face == 1)*/ {
      output[0] = input[0];
      output[1] = input[0];
      output[2] = input[0];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<1, 4>> dgAnalyze_1D_O4(std::array<T, squareSize<1, 4>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 4> buffer;
   std::array<T, 4> output;
   static U const C0 = U(3.4785484513745385e-01);
   static U const C1 = U(6.5214515486254609e-01);
   static U const C2 = U(2.2171699031897615e-01);
   static U const C3 = U(2.9955043831178735e-01);
   static U const C4 = U(2.1300321680756459e-01);
   static U const C5 = U(1.0600771525769923e-01);
   static U const C6 = U(2.6850642010792947e-01);
   output[0] = C0 * (input[0] + input[3]) + C1 * (input[1] + input[2]);
   output[1] = C2 * (-input[1] + input[2]) + C3 * (-input[0] + input[3]);
   output[2] = C4 * (input[0] - input[1] - input[2] + input[3]);
   output[3] = C5 * (-input[0] + input[3]) + C6 * (input[1] - input[2]);
   return output;
}

template<typename T>
std::array<T, squareSize<1, 4>> dgSynthesize_1D_O4(std::array<T, triangleSize<1, 4>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 4> buffer;
   std::array<T, 4> output;
   static U const C0 = U(3.0474698495520619e-01);
   static U const C1 = U(6.1233362071871378e-01);
   static U const C2 = U(8.6113631159405257e-01);
   static U const C3 = U(3.2661933500442808e-01);
   static U const C4 = U(3.3998104358485626e-01);
   static U const C5 = U(4.1172799967289958e-01);
   output[0] = -C0 * input[3] + C1 * input[2] - C2 * input[1] + input[0];
   output[1] = -C3 * input[2] - C4 * input[1] + C5 * input[3] + input[0];
   output[2] = -C3 * input[2] + C4 * input[1] - C5 * input[3] + input[0];
   output[3] = C0 * input[3] + C1 * input[2] + C2 * input[1] + input[0];
   return output;
}

template<typename T>
std::array<T, triangleSize<1, 4>> dgMassInverse_1D_O4(std::array<T, triangleSize<1, 4>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(5.0000000000000000e-01);
   static U const C1 = U(1.4999999999999998e+00);
   static U const C2 = U(2.5000000000000004e+00);
   static U const C3 = U(3.5000000000000000e+00);
   std::array<T, triangleSize<1, 4>> output;
   output[0] = C0 * input[0];
   output[1] = C1 * input[1];
   output[2] = C2 * input[2];
   output[3] = C3 * input[3];
   return output;
}

template<typename T>
std::array<T, triangleSize<1, 4>> dgStiffness_1D_O4(int dimension, std::array<T, triangleSize<1, 4>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(2.0000000000000000e+00);
   std::array<T, triangleSize<1, 4>> output;
   output[0] = U(0.0);
   output[1] = C0 * input[0];
   output[2] = C0 * input[1];
   output[3] = C0 * (input[0] + input[2]);
   return output;
}

template<typename T>
std::array<T, triangleSize<0, 4>> dgTrace_1D_O4(int face, std::array<T, triangleSize<1, 4>> const& input) {
   std::array<T, triangleSize<0, 4>> output;
   if(face == 0) {
      output[0] = input[0] - input[1] + input[2] - input[3];
   } else /*if(face == 1)*/ {
      output[0] = input[0] + input[1] + input[2] + input[3];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<1, 4>> dgTraceInverse_1D_O4(int face, std::array<T, triangleSize<0, 4>> const& input) {
   std::array<T, triangleSize<1, 4>> output;
   if(face == 0) {
      output[0] = input[0];
      output[1] = -input[0];
      output[2] = input[0];
      output[3] = -input[0];
   } else /*if(face == 1)*/ {
      output[0] = input[0];
      output[1] = input[0];
      output[2] = input[0];
      output[3] = input[0];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<1, 5>> dgAnalyze_1D_O5(std::array<T, squareSize<1, 5>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 5> buffer;
   std::array<T, 5> output;
   static U const C0 = U(2.3692688505618908e-01);
   static U const C1 = U(4.7862867049936647e-01);
   static U const C2 = U(5.6888888888888889e-01);
   static U const C3 = U(2.1469836819894497e-01);
   static U const C4 = U(2.5772685000059420e-01);
   static U const C5 = U(3.1147336576387012e-02);
   static U const C6 = U(1.7336955879860924e-01);
   static U const C7 = U(2.8444444444444444e-01);
   static U const C8 = U(1.1870775467166646e-01);
   static U const C9 = U(1.9977104139692381e-01);
   static U const C10 = U(5.8221336871210075e-02);
   static U const C11 = U(1.6488800353787675e-01);
   static U const C12 = U(2.1333333333333335e-01);
   output[0] = C0 * (input[0] + input[4]) + C1 * (input[1] + input[3]) + C2 * input[2];
   output[1] = C3 * (-input[0] + input[4]) + C4 * (-input[1] + input[3]);
   output[2] = C5 * (-input[1] - input[3]) + C6 * (input[0] + input[4]) - C7 * input[2];
   output[3] = C8 * (-input[0] + input[4]) + C9 * (input[1] - input[3]);
   output[4] = C10 * (input[0] + input[4]) + C11 * (-input[1] - input[3]) + C12 * input[2];
   return output;
}

template<typename T>
std::array<T, squareSize<1, 5>> dgSynthesize_1D_O5(std::array<T, triangleSize<1, 5>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 5> buffer;
   std::array<T, 5> output;
   static U const C0 = U(2.4573545909491201e-01);
   static U const C1 = U(5.0103117104466199e-01);
   static U const C2 = U(7.3174286977813119e-01);
   static U const C3 = U(9.0617984593866396e-01);
   static U const C4 = U(6.5076203111464545e-02);
   static U const C5 = U(3.4450089119367744e-01);
   static U const C6 = U(4.1738210372666812e-01);
   static U const C7 = U(5.3846931010568311e-01);
   static U const C8 = U(3.7500000000000000e-01);
   static U const C9 = U(5.0000000000000000e-01);
   output[0] = C0 * input[4] - C1 * input[3] + C2 * input[2] - C3 * input[1] + input[0];
   output[1] = -C4 * input[2] - C5 * input[4] + C6 * input[3] - C7 * input[1] + input[0];
   output[2] = C8 * input[4] - C9 * input[2] + input[0];
   output[3] = -C4 * input[2] - C5 * input[4] - C6 * input[3] + C7 * input[1] + input[0];
   output[4] = C0 * input[4] + C1 * input[3] + C2 * input[2] + C3 * input[1] + input[0];
   return output;
}

template<typename T>
std::array<T, triangleSize<1, 5>> dgMassInverse_1D_O5(std::array<T, triangleSize<1, 5>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(5.0000000000000000e-01);
   static U const C1 = U(1.5000000000000000e+00);
   static U const C2 = U(2.5000000000000000e+00);
   static U const C3 = U(3.5000000000000000e+00);
   static U const C4 = U(4.5000000000000000e+00);
   std::array<T, triangleSize<1, 5>> output;
   output[0] = C0 * input[0];
   output[1] = C1 * input[1];
   output[2] = C2 * input[2];
   output[3] = C3 * input[3];
   output[4] = C4 * input[4];
   return output;
}

template<typename T>
std::array<T, triangleSize<1, 5>> dgStiffness_1D_O5(int dimension, std::array<T, triangleSize<1, 5>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(2.0000000000000000e+00);
   std::array<T, triangleSize<1, 5>> output;
   output[0] = U(0.0);
   output[1] = C0 * input[0];
   output[2] = C0 * input[1];
   output[3] = C0 * (input[0] + input[2]);
   output[4] = C0 * (input[1] + input[3]);
   return output;
}

template<typename T>
std::array<T, triangleSize<0, 5>> dgTrace_1D_O5(int face, std::array<T, triangleSize<1, 5>> const& input) {
   std::array<T, triangleSize<0, 5>> output;
   if(face == 0) {
      output[0] = input[0] - input[1] + input[2] - input[3] + input[4];
   } else /*if(face == 1)*/ {
      output[0] = input[0] + input[1] + input[2] + input[3] + input[4];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<1, 5>> dgTraceInverse_1D_O5(int face, std::array<T, triangleSize<0, 5>> const& input) {
   std::array<T, triangleSize<1, 5>> output;
   if(face == 0) {
      output[0] = input[0];
      output[1] = -input[0];
      output[2] = input[0];
      output[3] = -input[0];
      output[4] = input[0];
   } else /*if(face == 1)*/ {
      output[0] = input[0];
      output[1] = input[0];
      output[2] = input[0];
      output[3] = input[0];
      output[4] = input[0];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<2, 1>> dgAnalyze_2D_O1(std::array<T, squareSize<2, 1>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 1> buffer;
   std::array<T, 1> output;
   static U const C0 = U(2.0000000000000000e+00);
   buffer[0] = C0 * input[0];
   output[0] = C0 * buffer[0];
   return output;
}

template<typename T>
std::array<T, squareSize<2, 1>> dgSynthesize_2D_O1(std::array<T, triangleSize<2, 1>> const& input) {
   std::array<T, 1> buffer;
   std::array<T, 1> output;
   buffer[0] = input[0];
   output[0] = buffer[0];
   return output;
}

template<typename T>
std::array<T, triangleSize<2, 1>> dgMassInverse_2D_O1(std::array<T, triangleSize<2, 1>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(2.5000000000000000e-01);
   std::array<T, triangleSize<2, 1>> output;
   output[0] = C0 * input[0];
   return output;
}

template<typename T>
std::array<T, triangleSize<2, 1>> dgStiffness_2D_O1(int dimension, std::array<T, triangleSize<2, 1>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, triangleSize<2, 1>> output;
   if(dimension == 0) {
      output[0] = U(0.0);
   } else /*if(dimension == 1)*/ {
      output[0] = U(0.0);
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<1, 1>> dgTrace_2D_O1(int face, std::array<T, triangleSize<2, 1>> const& input) {
   std::array<T, triangleSize<1, 1>> output;
   if(face == 0) {
      output[0] = input[0];
   } else if(face == 1) {
      output[0] = input[0];
   } else if(face == 2) {
      output[0] = input[0];
   } else /*if(face == 3)*/ {
      output[0] = input[0];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<2, 1>> dgTraceInverse_2D_O1(int face, std::array<T, triangleSize<1, 1>> const& input) {
   std::array<T, triangleSize<2, 1>> output;
   if(face == 0) {
      output[0] = input[0];
   } else if(face == 1) {
      output[0] = input[0];
   } else if(face == 2) {
      output[0] = input[0];
   } else /*if(face == 3)*/ {
      output[0] = input[0];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<2, 2>> dgAnalyze_2D_O2(std::array<T, squareSize<2, 2>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 4> buffer;
   std::array<T, 3> output;
   static U const C0 = U(5.7735026918962573e-01);
   buffer[0] = input[0] + input[1];
   buffer[1] = input[2] + input[3];
   buffer[2] = C0 * (-input[0] + input[1]);
   buffer[3] = C0 * (-input[2] + input[3]);
   output[0] = buffer[0] + buffer[1];
   output[1] = buffer[2] + buffer[3];
   output[2] = C0 * (-buffer[0] + buffer[1]);
   return output;
}

template<typename T>
std::array<T, squareSize<2, 2>> dgSynthesize_2D_O2(std::array<T, triangleSize<2, 2>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 4> buffer;
   std::array<T, 4> output;
   static U const C0 = U(5.7735026918962573e-01);
   buffer[0] = -C0 * input[2] + input[0];
   buffer[1] = C0 * input[2] + input[0];
   buffer[2] = input[1];
   buffer[3] = input[1];
   output[0] = -C0 * buffer[2] + buffer[0];
   output[1] = C0 * buffer[2] + buffer[0];
   output[2] = -C0 * buffer[3] + buffer[1];
   output[3] = C0 * buffer[3] + buffer[1];
   return output;
}

template<typename T>
std::array<T, triangleSize<2, 2>> dgMassInverse_2D_O2(std::array<T, triangleSize<2, 2>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(2.5000000000000000e-01);
   static U const C1 = U(7.5000000000000000e-01);
   std::array<T, triangleSize<2, 2>> output;
   output[0] = C0 * input[0];
   output[1] = C1 * input[1];
   output[2] = C1 * input[2];
   return output;
}

template<typename T>
std::array<T, triangleSize<2, 2>> dgStiffness_2D_O2(int dimension, std::array<T, triangleSize<2, 2>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(2.0000000000000000e+00);
   std::array<T, triangleSize<2, 2>> output;
   if(dimension == 0) {
      output[0] = U(0.0);
      output[1] = U(0.0);
      output[2] = C0 * input[0];
   } else /*if(dimension == 1)*/ {
      output[0] = U(0.0);
      output[1] = C0 * input[0];
      output[2] = U(0.0);
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<1, 2>> dgTrace_2D_O2(int face, std::array<T, triangleSize<2, 2>> const& input) {
   std::array<T, triangleSize<1, 2>> output;
   if(face == 0) {
      output[0] = input[0] - input[2];
      output[1] = input[1];
   } else if(face == 1) {
      output[0] = input[0] + input[2];
      output[1] = input[1];
   } else if(face == 2) {
      output[0] = input[0] - input[1];
      output[1] = input[2];
   } else /*if(face == 3)*/ {
      output[0] = input[0] + input[1];
      output[1] = input[2];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<2, 2>> dgTraceInverse_2D_O2(int face, std::array<T, triangleSize<1, 2>> const& input) {
   std::array<T, triangleSize<2, 2>> output;
   if(face == 0) {
      output[0] = input[0];
      output[1] = input[1];
      output[2] = -input[0];
   } else if(face == 1) {
      output[0] = input[0];
      output[1] = input[1];
      output[2] = input[0];
   } else if(face == 2) {
      output[0] = input[0];
      output[1] = -input[0];
      output[2] = input[1];
   } else /*if(face == 3)*/ {
      output[0] = input[0];
      output[1] = input[0];
      output[2] = input[1];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<2, 3>> dgAnalyze_2D_O3(std::array<T, squareSize<2, 3>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 9> buffer;
   std::array<T, 6> output;
   static U const C0 = U(5.5555555555555558e-01);
   static U const C1 = U(8.8888888888888884e-01);
   static U const C2 = U(4.3033148291193524e-01);
   static U const C3 = U(2.2222222222222221e-01);
   static U const C4 = U(4.4444444444444442e-01);
   buffer[0] = C0 * (input[0] + input[2]) + C1 * input[1];
   buffer[1] = C0 * (input[3] + input[5]) + C1 * input[4];
   buffer[2] = C0 * (input[6] + input[8]) + C1 * input[7];
   buffer[3] = C2 * (-input[0] + input[2]);
   buffer[4] = C2 * (-input[3] + input[5]);
   buffer[5] = C2 * (-input[6] + input[8]);
   buffer[6] = C3 * (input[0] + input[2]) - C4 * input[1];
   buffer[7] = C3 * (input[3] + input[5]) - C4 * input[4];
   buffer[8] = C3 * (input[6] + input[8]) - C4 * input[7];
   output[0] = C0 * (buffer[0] + buffer[2]) + C1 * buffer[1];
   output[1] = C0 * (buffer[3] + buffer[5]) + C1 * buffer[4];
   output[2] = C2 * (-buffer[0] + buffer[2]);
   output[3] = C0 * (buffer[6] + buffer[8]) + C1 * buffer[7];
   output[4] = C2 * (-buffer[3] + buffer[5]);
   output[5] = C3 * (buffer[0] + buffer[2]) - C4 * buffer[1];
   return output;
}

template<typename T>
std::array<T, squareSize<2, 3>> dgSynthesize_2D_O3(std::array<T, triangleSize<2, 3>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 9> buffer;
   std::array<T, 9> output;
   static U const C0 = U(4.0000000000000002e-01);
   static U const C1 = U(7.7459666924148340e-01);
   static U const C2 = U(5.0000000000000000e-01);
   buffer[0] = C0 * input[5] - C1 * input[2] + input[0];
   buffer[1] = -C2 * input[5] + input[0];
   buffer[2] = C0 * input[5] + C1 * input[2] + input[0];
   buffer[3] = -C1 * input[4] + input[1];
   buffer[4] = input[1];
   buffer[5] = C1 * input[4] + input[1];
   buffer[6] = input[3];
   buffer[7] = input[3];
   buffer[8] = input[3];
   output[0] = C0 * buffer[6] - C1 * buffer[3] + buffer[0];
   output[1] = -C2 * buffer[6] + buffer[0];
   output[2] = C0 * buffer[6] + C1 * buffer[3] + buffer[0];
   output[3] = C0 * buffer[7] - C1 * buffer[4] + buffer[1];
   output[4] = -C2 * buffer[7] + buffer[1];
   output[5] = C0 * buffer[7] + C1 * buffer[4] + buffer[1];
   output[6] = C0 * buffer[8] - C1 * buffer[5] + buffer[2];
   output[7] = -C2 * buffer[8] + buffer[2];
   output[8] = C0 * buffer[8] + C1 * buffer[5] + buffer[2];
   return output;
}

template<typename T>
std::array<T, triangleSize<2, 3>> dgMassInverse_2D_O3(std::array<T, triangleSize<2, 3>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(2.5000000000000000e-01);
   static U const C1 = U(7.5000000000000000e-01);
   static U const C2 = U(1.2500000000000000e+00);
   static U const C3 = U(2.2500000000000000e+00);
   std::array<T, triangleSize<2, 3>> output;
   output[0] = C0 * input[0];
   output[1] = C1 * input[1];
   output[2] = C1 * input[2];
   output[3] = C2 * input[3];
   output[4] = C3 * input[4];
   output[5] = C2 * input[5];
   return output;
}

template<typename T>
std::array<T, triangleSize<2, 3>> dgStiffness_2D_O3(int dimension, std::array<T, triangleSize<2, 3>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(2.0000000000000000e+00);
   std::array<T, triangleSize<2, 3>> output;
   if(dimension == 0) {
      output[0] = U(0.0);
      output[1] = U(0.0);
      output[2] = C0 * input[0];
      output[3] = U(0.0);
      output[4] = C0 * input[1];
      output[5] = C0 * input[2];
   } else /*if(dimension == 1)*/ {
      output[0] = U(0.0);
      output[1] = C0 * input[0];
      output[2] = U(0.0);
      output[3] = C0 * input[1];
      output[4] = C0 * input[2];
      output[5] = U(0.0);
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<1, 3>> dgTrace_2D_O3(int face, std::array<T, triangleSize<2, 3>> const& input) {
   std::array<T, triangleSize<1, 3>> output;
   if(face == 0) {
      output[0] = input[0] - input[2] + input[5];
      output[1] = input[1] - input[4];
      output[2] = input[3];
   } else if(face == 1) {
      output[0] = input[0] + input[2] + input[5];
      output[1] = input[1] + input[4];
      output[2] = input[3];
   } else if(face == 2) {
      output[0] = input[0] - input[1] + input[3];
      output[1] = input[2] - input[4];
      output[2] = input[5];
   } else /*if(face == 3)*/ {
      output[0] = input[0] + input[1] + input[3];
      output[1] = input[2] + input[4];
      output[2] = input[5];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<2, 3>> dgTraceInverse_2D_O3(int face, std::array<T, triangleSize<1, 3>> const& input) {
   std::array<T, triangleSize<2, 3>> output;
   if(face == 0) {
      output[0] = input[0];
      output[1] = input[1];
      output[2] = -input[0];
      output[3] = input[2];
      output[4] = -input[1];
      output[5] = input[0];
   } else if(face == 1) {
      output[0] = input[0];
      output[1] = input[1];
      output[2] = input[0];
      output[3] = input[2];
      output[4] = input[1];
      output[5] = input[0];
   } else if(face == 2) {
      output[0] = input[0];
      output[1] = -input[0];
      output[2] = input[1];
      output[3] = input[0];
      output[4] = -input[1];
      output[5] = input[2];
   } else /*if(face == 3)*/ {
      output[0] = input[0];
      output[1] = input[0];
      output[2] = input[1];
      output[3] = input[0];
      output[4] = input[1];
      output[5] = input[2];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<2, 4>> dgAnalyze_2D_O4(std::array<T, squareSize<2, 4>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 16> buffer;
   std::array<T, 10> output;
   static U const C0 = U(3.4785484513745385e-01);
   static U const C1 = U(6.5214515486254609e-01);
   static U const C2 = U(2.2171699031897615e-01);
   static U const C3 = U(2.9955043831178735e-01);
   static U const C4 = U(2.1300321680756459e-01);
   static U const C5 = U(1.0600771525769923e-01);
   static U const C6 = U(2.6850642010792947e-01);
   buffer[0] = C0 * (input[0] + input[3]) + C1 * (input[1] + input[2]);
   buffer[1] = C0 * (input[4] + input[7]) + C1 * (input[5] + input[6]);
   buffer[2] = C0 * (input[8] + input[11]) + C1 * (input[9] + input[10]);
   buffer[3] = C0 * (input[12] + input[15]) + C1 * (input[13] + input[14]);
   buffer[4] = C2 * (-input[1] + input[2]) + C3 * (-input[0] + input[3]);
   buffer[5] = C2 * (-input[5] + input[6]) + C3 * (-input[4] + input[7]);
   buffer[6] = C2 * (-input[9] + input[10]) + C3 * (-input[8] + input[11]);
   buffer[7] = C2 * (-input[13] + input[14]) + C3 * (-input[12] + input[15]);
   buffer[8] = C4 * (input[0] - input[1] - input[2] + input[3]);
   buffer[9] = C4 * (input[4] - input[5] - input[6] + input[7]);
   buffer[10] = C4 * (input[8] - input[9] - input[10] + input[11]);
   buffer[11] = C4 * (input[12] - input[13] - input[14] + input[15]);
   buffer[12] = C5 * (-input[0] + input[3]) + C6 * (input[1] - input[2]);
   buffer[13] = C5 * (-input[4] + input[7]) + C6 * (input[5] - input[6]);
   buffer[14] = C5 * (-input[8] + input[11]) + C6 * (input[9] - input[10]);
   buffer[15] = C5 * (-input[12] + input[15]) + C6 * (input[13] - input[14]);
   output[0] = C0 * (buffer[0] + buffer[3]) + C1 * (buffer[1] + buffer[2]);
   output[1] = C0 * (buffer[4] + buffer[7]) + C1 * (buffer[5] + buffer[6]);
   output[2] = C2 * (-buffer[1] + buffer[2]) + C3 * (-buffer[0] + buffer[3]);
   output[3] = C0 * (buffer[8] + buffer[11]) + C1 * (buffer[9] + buffer[10]);
   output[4] = C2 * (-buffer[5] + buffer[6]) + C3 * (-buffer[4] + buffer[7]);
   output[5] = C4 * (buffer[0] - buffer[1] - buffer[2] + buffer[3]);
   output[6] = C0 * (buffer[12] + buffer[15]) + C1 * (buffer[13] + buffer[14]);
   output[7] = C2 * (-buffer[9] + buffer[10]) + C3 * (-buffer[8] + buffer[11]);
   output[8] = C4 * (buffer[4] - buffer[5] - buffer[6] + buffer[7]);
   output[9] = C5 * (-buffer[0] + buffer[3]) + C6 * (buffer[1] - buffer[2]);
   return output;
}

template<typename T>
std::array<T, squareSize<2, 4>> dgSynthesize_2D_O4(std::array<T, triangleSize<2, 4>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 16> buffer;
   std::array<T, 16> output;
   static U const C0 = U(3.0474698495520619e-01);
   static U const C1 = U(6.1233362071871378e-01);
   static U const C2 = U(8.6113631159405257e-01);
   static U const C3 = U(3.2661933500442808e-01);
   static U const C4 = U(3.3998104358485626e-01);
   static U const C5 = U(4.1172799967289958e-01);
   buffer[0] = -C0 * input[9] + C1 * input[5] - C2 * input[2] + input[0];
   buffer[1] = -C3 * input[5] - C4 * input[2] + C5 * input[9] + input[0];
   buffer[2] = -C3 * input[5] + C4 * input[2] - C5 * input[9] + input[0];
   buffer[3] = C0 * input[9] + C1 * input[5] + C2 * input[2] + input[0];
   buffer[4] = C1 * input[8] - C2 * input[4] + input[1];
   buffer[5] = -C3 * input[8] - C4 * input[4] + input[1];
   buffer[6] = -C3 * input[8] + C4 * input[4] + input[1];
   buffer[7] = C1 * input[8] + C2 * input[4] + input[1];
   buffer[8] = -C2 * input[7] + input[3];
   buffer[9] = -C4 * input[7] + input[3];
   buffer[10] = C4 * input[7] + input[3];
   buffer[11] = C2 * input[7] + input[3];
   buffer[12] = input[6];
   buffer[13] = input[6];
   buffer[14] = input[6];
   buffer[15] = input[6];
   output[0] = -C0 * buffer[12] + C1 * buffer[8] - C2 * buffer[4] + buffer[0];
   output[1] = -C3 * buffer[8] - C4 * buffer[4] + C5 * buffer[12] + buffer[0];
   output[2] = -C3 * buffer[8] + C4 * buffer[4] - C5 * buffer[12] + buffer[0];
   output[3] = C0 * buffer[12] + C1 * buffer[8] + C2 * buffer[4] + buffer[0];
   output[4] = -C0 * buffer[13] + C1 * buffer[9] - C2 * buffer[5] + buffer[1];
   output[5] = -C3 * buffer[9] - C4 * buffer[5] + C5 * buffer[13] + buffer[1];
   output[6] = -C3 * buffer[9] + C4 * buffer[5] - C5 * buffer[13] + buffer[1];
   output[7] = C0 * buffer[13] + C1 * buffer[9] + C2 * buffer[5] + buffer[1];
   output[8] = -C0 * buffer[14] + C1 * buffer[10] - C2 * buffer[6] + buffer[2];
   output[9] = -C3 * buffer[10] - C4 * buffer[6] + C5 * buffer[14] + buffer[2];
   output[10] = -C3 * buffer[10] + C4 * buffer[6] - C5 * buffer[14] + buffer[2];
   output[11] = C0 * buffer[14] + C1 * buffer[10] + C2 * buffer[6] + buffer[2];
   output[12] = -C0 * buffer[15] + C1 * buffer[11] - C2 * buffer[7] + buffer[3];
   output[13] = -C3 * buffer[11] - C4 * buffer[7] + C5 * buffer[15] + buffer[3];
   output[14] = -C3 * buffer[11] + C4 * buffer[7] - C5 * buffer[15] + buffer[3];
   output[15] = C0 * buffer[15] + C1 * buffer[11] + C2 * buffer[7] + buffer[3];
   return output;
}

template<typename T>
std::array<T, triangleSize<2, 4>> dgMassInverse_2D_O4(std::array<T, triangleSize<2, 4>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(2.5000000000000000e-01);
   static U const C1 = U(7.5000000000000000e-01);
   static U const C2 = U(1.2500000000000002e+00);
   static U const C3 = U(2.2499999999999996e+00);
   static U const C4 = U(1.7500000000000000e+00);
   static U const C5 = U(3.7500000000000004e+00);
   std::array<T, triangleSize<2, 4>> output;
   output[0] = C0 * input[0];
   output[1] = C1 * input[1];
   output[2] = C1 * input[2];
   output[3] = C2 * input[3];
   output[4] = C3 * input[4];
   output[5] = C2 * input[5];
   output[6] = C4 * input[6];
   output[7] = C5 * input[7];
   output[8] = C5 * input[8];
   output[9] = C4 * input[9];
   return output;
}

template<typename T>
std::array<T, triangleSize<2, 4>> dgStiffness_2D_O4(int dimension, std::array<T, triangleSize<2, 4>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(2.0000000000000000e+00);
   std::array<T, triangleSize<2, 4>> output;
   if(dimension == 0) {
      output[0] = U(0.0);
      output[1] = U(0.0);
      output[2] = C0 * input[0];
      output[3] = U(0.0);
      output[4] = C0 * input[1];
      output[5] = C0 * input[2];
      output[6] = U(0.0);
      output[7] = C0 * input[3];
      output[8] = C0 * input[4];
      output[9] = C0 * (input[0] + input[5]);
   } else /*if(dimension == 1)*/ {
      output[0] = U(0.0);
      output[1] = C0 * input[0];
      output[2] = U(0.0);
      output[3] = C0 * input[1];
      output[4] = C0 * input[2];
      output[5] = U(0.0);
      output[6] = C0 * (input[0] + input[3]);
      output[7] = C0 * input[4];
      output[8] = C0 * input[5];
      output[9] = U(0.0);
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<1, 4>> dgTrace_2D_O4(int face, std::array<T, triangleSize<2, 4>> const& input) {
   std::array<T, triangleSize<1, 4>> output;
   if(face == 0) {
      output[0] = input[0] - input[2] + input[5] - input[9];
      output[1] = input[1] - input[4] + input[8];
      output[2] = input[3] - input[7];
      output[3] = input[6];
   } else if(face == 1) {
      output[0] = input[0] + input[2] + input[5] + input[9];
      output[1] = input[1] + input[4] + input[8];
      output[2] = input[3] + input[7];
      output[3] = input[6];
   } else if(face == 2) {
      output[0] = input[0] - input[1] + input[3] - input[6];
      output[1] = input[2] - input[4] + input[7];
      output[2] = input[5] - input[8];
      output[3] = input[9];
   } else /*if(face == 3)*/ {
      output[0] = input[0] + input[1] + input[3] + input[6];
      output[1] = input[2] + input[4] + input[7];
      output[2] = input[5] + input[8];
      output[3] = input[9];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<2, 4>> dgTraceInverse_2D_O4(int face, std::array<T, triangleSize<1, 4>> const& input) {
   std::array<T, triangleSize<2, 4>> output;
   if(face == 0) {
      output[0] = input[0];
      output[1] = input[1];
      output[2] = -input[0];
      output[3] = input[2];
      output[4] = -input[1];
      output[5] = input[0];
      output[6] = input[3];
      output[7] = -input[2];
      output[8] = input[1];
      output[9] = -input[0];
   } else if(face == 1) {
      output[0] = input[0];
      output[1] = input[1];
      output[2] = input[0];
      output[3] = input[2];
      output[4] = input[1];
      output[5] = input[0];
      output[6] = input[3];
      output[7] = input[2];
      output[8] = input[1];
      output[9] = input[0];
   } else if(face == 2) {
      output[0] = input[0];
      output[1] = -input[0];
      output[2] = input[1];
      output[3] = input[0];
      output[4] = -input[1];
      output[5] = input[2];
      output[6] = -input[0];
      output[7] = input[1];
      output[8] = -input[2];
      output[9] = input[3];
   } else /*if(face == 3)*/ {
      output[0] = input[0];
      output[1] = input[0];
      output[2] = input[1];
      output[3] = input[0];
      output[4] = input[1];
      output[5] = input[2];
      output[6] = input[0];
      output[7] = input[1];
      output[8] = input[2];
      output[9] = input[3];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<2, 5>> dgAnalyze_2D_O5(std::array<T, squareSize<2, 5>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 25> buffer;
   std::array<T, 15> output;
   static U const C0 = U(2.3692688505618908e-01);
   static U const C1 = U(4.7862867049936647e-01);
   static U const C2 = U(5.6888888888888889e-01);
   static U const C3 = U(2.1469836819894497e-01);
   static U const C4 = U(2.5772685000059420e-01);
   static U const C5 = U(3.1147336576387012e-02);
   static U const C6 = U(1.7336955879860924e-01);
   static U const C7 = U(2.8444444444444444e-01);
   static U const C8 = U(1.1870775467166646e-01);
   static U const C9 = U(1.9977104139692381e-01);
   static U const C10 = U(5.8221336871210075e-02);
   static U const C11 = U(1.6488800353787675e-01);
   static U const C12 = U(2.1333333333333335e-01);
   buffer[0] = C0 * (input[0] + input[4]) + C1 * (input[1] + input[3]) + C2 * input[2];
   buffer[1] = C0 * (input[5] + input[9]) + C1 * (input[6] + input[8]) + C2 * input[7];
   buffer[2] = C0 * (input[10] + input[14]) + C1 * (input[11] + input[13]) + C2 * input[12];
   buffer[3] = C0 * (input[15] + input[19]) + C1 * (input[16] + input[18]) + C2 * input[17];
   buffer[4] = C0 * (input[20] + input[24]) + C1 * (input[21] + input[23]) + C2 * input[22];
   buffer[5] = C3 * (-input[0] + input[4]) + C4 * (-input[1] + input[3]);
   buffer[6] = C3 * (-input[5] + input[9]) + C4 * (-input[6] + input[8]);
   buffer[7] = C3 * (-input[10] + input[14]) + C4 * (-input[11] + input[13]);
   buffer[8] = C3 * (-input[15] + input[19]) + C4 * (-input[16] + input[18]);
   buffer[9] = C3 * (-input[20] + input[24]) + C4 * (-input[21] + input[23]);
   buffer[10] = C5 * (-input[1] - input[3]) + C6 * (input[0] + input[4]) - C7 * input[2];
   buffer[11] = C5 * (-input[6] - input[8]) + C6 * (input[5] + input[9]) - C7 * input[7];
   buffer[12] = C5 * (-input[11] - input[13]) + C6 * (input[10] + input[14]) - C7 * input[12];
   buffer[13] = C5 * (-input[16] - input[18]) + C6 * (input[15] + input[19]) - C7 * input[17];
   buffer[14] = C5 * (-input[21] - input[23]) + C6 * (input[20] + input[24]) - C7 * input[22];
   buffer[15] = C8 * (-input[0] + input[4]) + C9 * (input[1] - input[3]);
   buffer[16] = C8 * (-input[5] + input[9]) + C9 * (input[6] - input[8]);
   buffer[17] = C8 * (-input[10] + input[14]) + C9 * (input[11] - input[13]);
   buffer[18] = C8 * (-input[15] + input[19]) + C9 * (input[16] - input[18]);
   buffer[19] = C8 * (-input[20] + input[24]) + C9 * (input[21] - input[23]);
   buffer[20] = C10 * (input[0] + input[4]) + C11 * (-input[1] - input[3]) + C12 * input[2];
   buffer[21] = C10 * (input[5] + input[9]) + C11 * (-input[6] - input[8]) + C12 * input[7];
   buffer[22] = C10 * (input[10] + input[14]) + C11 * (-input[11] - input[13]) + C12 * input[12];
   buffer[23] = C10 * (input[15] + input[19]) + C11 * (-input[16] - input[18]) + C12 * input[17];
   buffer[24] = C10 * (input[20] + input[24]) + C11 * (-input[21] - input[23]) + C12 * input[22];
   output[0] = C0 * (buffer[0] + buffer[4]) + C1 * (buffer[1] + buffer[3]) + C2 * buffer[2];
   output[1] = C0 * (buffer[5] + buffer[9]) + C1 * (buffer[6] + buffer[8]) + C2 * buffer[7];
   output[2] = C3 * (-buffer[0] + buffer[4]) + C4 * (-buffer[1] + buffer[3]);
   output[3] = C0 * (buffer[10] + buffer[14]) + C1 * (buffer[11] + buffer[13]) + C2 * buffer[12];
   output[4] = C3 * (-buffer[5] + buffer[9]) + C4 * (-buffer[6] + buffer[8]);
   output[5] = C5 * (-buffer[1] - buffer[3]) + C6 * (buffer[0] + buffer[4]) - C7 * buffer[2];
   output[6] = C0 * (buffer[15] + buffer[19]) + C1 * (buffer[16] + buffer[18]) + C2 * buffer[17];
   output[7] = C3 * (-buffer[10] + buffer[14]) + C4 * (-buffer[11] + buffer[13]);
   output[8] = C5 * (-buffer[6] - buffer[8]) + C6 * (buffer[5] + buffer[9]) - C7 * buffer[7];
   output[9] = C8 * (-buffer[0] + buffer[4]) + C9 * (buffer[1] - buffer[3]);
   output[10] = C0 * (buffer[20] + buffer[24]) + C1 * (buffer[21] + buffer[23]) + C2 * buffer[22];
   output[11] = C3 * (-buffer[15] + buffer[19]) + C4 * (-buffer[16] + buffer[18]);
   output[12] = C5 * (-buffer[11] - buffer[13]) + C6 * (buffer[10] + buffer[14]) - C7 * buffer[12];
   output[13] = C8 * (-buffer[5] + buffer[9]) + C9 * (buffer[6] - buffer[8]);
   output[14] = C10 * (buffer[0] + buffer[4]) + C11 * (-buffer[1] - buffer[3]) + C12 * buffer[2];
   return output;
}

template<typename T>
std::array<T, squareSize<2, 5>> dgSynthesize_2D_O5(std::array<T, triangleSize<2, 5>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 25> buffer;
   std::array<T, 25> output;
   static U const C0 = U(2.4573545909491201e-01);
   static U const C1 = U(5.0103117104466199e-01);
   static U const C2 = U(7.3174286977813119e-01);
   static U const C3 = U(9.0617984593866396e-01);
   static U const C4 = U(6.5076203111464545e-02);
   static U const C5 = U(3.4450089119367744e-01);
   static U const C6 = U(4.1738210372666812e-01);
   static U const C7 = U(5.3846931010568311e-01);
   static U const C8 = U(3.7500000000000000e-01);
   static U const C9 = U(5.0000000000000000e-01);
   buffer[0] = C0 * input[14] - C1 * input[9] + C2 * input[5] - C3 * input[2] + input[0];
   buffer[1] = -C4 * input[5] - C5 * input[14] + C6 * input[9] - C7 * input[2] + input[0];
   buffer[2] = C8 * input[14] - C9 * input[5] + input[0];
   buffer[3] = -C4 * input[5] - C5 * input[14] - C6 * input[9] + C7 * input[2] + input[0];
   buffer[4] = C0 * input[14] + C1 * input[9] + C2 * input[5] + C3 * input[2] + input[0];
   buffer[5] = -C1 * input[13] + C2 * input[8] - C3 * input[4] + input[1];
   buffer[6] = -C4 * input[8] + C6 * input[13] - C7 * input[4] + input[1];
   buffer[7] = -C9 * input[8] + input[1];
   buffer[8] = -C4 * input[8] - C6 * input[13] + C7 * input[4] + input[1];
   buffer[9] = C1 * input[13] + C2 * input[8] + C3 * input[4] + input[1];
   buffer[10] = C2 * input[12] - C3 * input[7] + input[3];
   buffer[11] = -C4 * input[12] - C7 * input[7] + input[3];
   buffer[12] = -C9 * input[12] + input[3];
   buffer[13] = -C4 * input[12] + C7 * input[7] + input[3];
   buffer[14] = C2 * input[12] + C3 * input[7] + input[3];
   buffer[15] = -C3 * input[11] + input[6];
   buffer[16] = -C7 * input[11] + input[6];
   buffer[17] = input[6];
   buffer[18] = C7 * input[11] + input[6];
   buffer[19] = C3 * input[11] + input[6];
   buffer[20] = input[10];
   buffer[21] = input[10];
   buffer[22] = input[10];
   buffer[23] = input[10];
   buffer[24] = input[10];
   output[0] = C0 * buffer[20] - C1 * buffer[15] + C2 * buffer[10] - C3 * buffer[5] + buffer[0];
   output[1] = -C4 * buffer[10] - C5 * buffer[20] + C6 * buffer[15] - C7 * buffer[5] + buffer[0];
   output[2] = C8 * buffer[20] - C9 * buffer[10] + buffer[0];
   output[3] = -C4 * buffer[10] - C5 * buffer[20] - C6 * buffer[15] + C7 * buffer[5] + buffer[0];
   output[4] = C0 * buffer[20] + C1 * buffer[15] + C2 * buffer[10] + C3 * buffer[5] + buffer[0];
   output[5] = C0 * buffer[21] - C1 * buffer[16] + C2 * buffer[11] - C3 * buffer[6] + buffer[1];
   output[6] = -C4 * buffer[11] - C5 * buffer[21] + C6 * buffer[16] - C7 * buffer[6] + buffer[1];
   output[7] = C8 * buffer[21] - C9 * buffer[11] + buffer[1];
   output[8] = -C4 * buffer[11] - C5 * buffer[21] - C6 * buffer[16] + C7 * buffer[6] + buffer[1];
   output[9] = C0 * buffer[21] + C1 * buffer[16] + C2 * buffer[11] + C3 * buffer[6] + buffer[1];
   output[10] = C0 * buffer[22] - C1 * buffer[17] + C2 * buffer[12] - C3 * buffer[7] + buffer[2];
   output[11] = -C4 * buffer[12] - C5 * buffer[22] + C6 * buffer[17] - C7 * buffer[7] + buffer[2];
   output[12] = C8 * buffer[22] - C9 * buffer[12] + buffer[2];
   output[13] = -C4 * buffer[12] - C5 * buffer[22] - C6 * buffer[17] + C7 * buffer[7] + buffer[2];
   output[14] = C0 * buffer[22] + C1 * buffer[17] + C2 * buffer[12] + C3 * buffer[7] + buffer[2];
   output[15] = C0 * buffer[23] - C1 * buffer[18] + C2 * buffer[13] - C3 * buffer[8] + buffer[3];
   output[16] = -C4 * buffer[13] - C5 * buffer[23] + C6 * buffer[18] - C7 * buffer[8] + buffer[3];
   output[17] = C8 * buffer[23] - C9 * buffer[13] + buffer[3];
   output[18] = -C4 * buffer[13] - C5 * buffer[23] - C6 * buffer[18] + C7 * buffer[8] + buffer[3];
   output[19] = C0 * buffer[23] + C1 * buffer[18] + C2 * buffer[13] + C3 * buffer[8] + buffer[3];
   output[20] = C0 * buffer[24] - C1 * buffer[19] + C2 * buffer[14] - C3 * buffer[9] + buffer[4];
   output[21] = -C4 * buffer[14] - C5 * buffer[24] + C6 * buffer[19] - C7 * buffer[9] + buffer[4];
   output[22] = C8 * buffer[24] - C9 * buffer[14] + buffer[4];
   output[23] = -C4 * buffer[14] - C5 * buffer[24] - C6 * buffer[19] + C7 * buffer[9] + buffer[4];
   output[24] = C0 * buffer[24] + C1 * buffer[19] + C2 * buffer[14] + C3 * buffer[9] + buffer[4];
   return output;
}

template<typename T>
std::array<T, triangleSize<2, 5>> dgMassInverse_2D_O5(std::array<T, triangleSize<2, 5>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(2.5000000000000000e-01);
   static U const C1 = U(7.5000000000000000e-01);
   static U const C2 = U(1.2500000000000000e+00);
   static U const C3 = U(2.2500000000000000e+00);
   static U const C4 = U(1.7500000000000000e+00);
   static U const C5 = U(3.7500000000000000e+00);
   static U const C6 = U(5.2500000000000000e+00);
   static U const C7 = U(6.2500000000000000e+00);
   std::array<T, triangleSize<2, 5>> output;
   output[0] = C0 * input[0];
   output[1] = C1 * input[1];
   output[2] = C1 * input[2];
   output[3] = C2 * input[3];
   output[4] = C3 * input[4];
   output[5] = C2 * input[5];
   output[6] = C4 * input[6];
   output[7] = C5 * input[7];
   output[8] = C5 * input[8];
   output[9] = C4 * input[9];
   output[10] = C3 * input[10];
   output[11] = C6 * input[11];
   output[12] = C7 * input[12];
   output[13] = C6 * input[13];
   output[14] = C3 * input[14];
   return output;
}

template<typename T>
std::array<T, triangleSize<2, 5>> dgStiffness_2D_O5(int dimension, std::array<T, triangleSize<2, 5>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(2.0000000000000000e+00);
   std::array<T, triangleSize<2, 5>> output;
   if(dimension == 0) {
      output[0] = U(0.0);
      output[1] = U(0.0);
      output[2] = C0 * input[0];
      output[3] = U(0.0);
      output[4] = C0 * input[1];
      output[5] = C0 * input[2];
      output[6] = U(0.0);
      output[7] = C0 * input[3];
      output[8] = C0 * input[4];
      output[9] = C0 * (input[0] + input[5]);
      output[10] = U(0.0);
      output[11] = C0 * input[6];
      output[12] = C0 * input[7];
      output[13] = C0 * (input[1] + input[8]);
      output[14] = C0 * (input[2] + input[9]);
   } else /*if(dimension == 1)*/ {
      output[0] = U(0.0);
      output[1] = C0 * input[0];
      output[2] = U(0.0);
      output[3] = C0 * input[1];
      output[4] = C0 * input[2];
      output[5] = U(0.0);
      output[6] = C0 * (input[0] + input[3]);
      output[7] = C0 * input[4];
      output[8] = C0 * input[5];
      output[9] = U(0.0);
      output[10] = C0 * (input[1] + input[6]);
      output[11] = C0 * (input[2] + input[7]);
      output[12] = C0 * input[8];
      output[13] = C0 * input[9];
      output[14] = U(0.0);
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<1, 5>> dgTrace_2D_O5(int face, std::array<T, triangleSize<2, 5>> const& input) {
   std::array<T, triangleSize<1, 5>> output;
   if(face == 0) {
      output[0] = input[0] - input[2] + input[5] - input[9] + input[14];
      output[1] = input[1] - input[4] + input[8] - input[13];
      output[2] = input[3] - input[7] + input[12];
      output[3] = input[6] - input[11];
      output[4] = input[10];
   } else if(face == 1) {
      output[0] = input[0] + input[2] + input[5] + input[9] + input[14];
      output[1] = input[1] + input[4] + input[8] + input[13];
      output[2] = input[3] + input[7] + input[12];
      output[3] = input[6] + input[11];
      output[4] = input[10];
   } else if(face == 2) {
      output[0] = input[0] - input[1] + input[3] - input[6] + input[10];
      output[1] = input[2] - input[4] + input[7] - input[11];
      output[2] = input[5] - input[8] + input[12];
      output[3] = input[9] - input[13];
      output[4] = input[14];
   } else /*if(face == 3)*/ {
      output[0] = input[0] + input[1] + input[3] + input[6] + input[10];
      output[1] = input[2] + input[4] + input[7] + input[11];
      output[2] = input[5] + input[8] + input[12];
      output[3] = input[9] + input[13];
      output[4] = input[14];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<2, 5>> dgTraceInverse_2D_O5(int face, std::array<T, triangleSize<1, 5>> const& input) {
   std::array<T, triangleSize<2, 5>> output;
   if(face == 0) {
      output[0] = input[0];
      output[1] = input[1];
      output[2] = -input[0];
      output[3] = input[2];
      output[4] = -input[1];
      output[5] = input[0];
      output[6] = input[3];
      output[7] = -input[2];
      output[8] = input[1];
      output[9] = -input[0];
      output[10] = input[4];
      output[11] = -input[3];
      output[12] = input[2];
      output[13] = -input[1];
      output[14] = input[0];
   } else if(face == 1) {
      output[0] = input[0];
      output[1] = input[1];
      output[2] = input[0];
      output[3] = input[2];
      output[4] = input[1];
      output[5] = input[0];
      output[6] = input[3];
      output[7] = input[2];
      output[8] = input[1];
      output[9] = input[0];
      output[10] = input[4];
      output[11] = input[3];
      output[12] = input[2];
      output[13] = input[1];
      output[14] = input[0];
   } else if(face == 2) {
      output[0] = input[0];
      output[1] = -input[0];
      output[2] = input[1];
      output[3] = input[0];
      output[4] = -input[1];
      output[5] = input[2];
      output[6] = -input[0];
      output[7] = input[1];
      output[8] = -input[2];
      output[9] = input[3];
      output[10] = input[0];
      output[11] = -input[1];
      output[12] = input[2];
      output[13] = -input[3];
      output[14] = input[4];
   } else /*if(face == 3)*/ {
      output[0] = input[0];
      output[1] = input[0];
      output[2] = input[1];
      output[3] = input[0];
      output[4] = input[1];
      output[5] = input[2];
      output[6] = input[0];
      output[7] = input[1];
      output[8] = input[2];
      output[9] = input[3];
      output[10] = input[0];
      output[11] = input[1];
      output[12] = input[2];
      output[13] = input[3];
      output[14] = input[4];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<3, 1>> dgAnalyze_3D_O1(std::array<T, squareSize<3, 1>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 2> buffer;
   std::array<T, 1> output;
   static U const C0 = U(2.0000000000000000e+00);
   buffer[0] = C0 * input[0];
   buffer[1] = C0 * buffer[0];
   output[0] = C0 * buffer[1];
   return output;
}

template<typename T>
std::array<T, squareSize<3, 1>> dgSynthesize_3D_O1(std::array<T, triangleSize<3, 1>> const& input) {
   std::array<T, 2> buffer;
   std::array<T, 1> output;
   buffer[0] = input[0];
   buffer[1] = buffer[0];
   output[0] = buffer[1];
   return output;
}

template<typename T>
std::array<T, triangleSize<3, 1>> dgMassInverse_3D_O1(std::array<T, triangleSize<3, 1>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(1.2500000000000000e-01);
   std::array<T, triangleSize<3, 1>> output;
   output[0] = C0 * input[0];
   return output;
}

template<typename T>
std::array<T, triangleSize<3, 1>> dgStiffness_3D_O1(int dimension, std::array<T, triangleSize<3, 1>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, triangleSize<3, 1>> output;
   if(dimension == 0) {
      output[0] = U(0.0);
   } else if(dimension == 1) {
      output[0] = U(0.0);
   } else /*if(dimension == 2)*/ {
      output[0] = U(0.0);
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<2, 1>> dgTrace_3D_O1(int face, std::array<T, triangleSize<3, 1>> const& input) {
   std::array<T, triangleSize<2, 1>> output;
   if(face == 0) {
      output[0] = input[0];
   } else if(face == 1) {
      output[0] = input[0];
   } else if(face == 2) {
      output[0] = input[0];
   } else if(face == 3) {
      output[0] = input[0];
   } else if(face == 4) {
      output[0] = input[0];
   } else /*if(face == 5)*/ {
      output[0] = input[0];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<3, 1>> dgTraceInverse_3D_O1(int face, std::array<T, triangleSize<2, 1>> const& input) {
   std::array<T, triangleSize<3, 1>> output;
   if(face == 0) {
      output[0] = input[0];
   } else if(face == 1) {
      output[0] = input[0];
   } else if(face == 2) {
      output[0] = input[0];
   } else if(face == 3) {
      output[0] = input[0];
   } else if(face == 4) {
      output[0] = input[0];
   } else /*if(face == 5)*/ {
      output[0] = input[0];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<3, 2>> dgAnalyze_3D_O2(std::array<T, squareSize<3, 2>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 14> buffer;
   std::array<T, 4> output;
   static U const C0 = U(5.7735026918962573e-01);
   buffer[0] = input[0] + input[1];
   buffer[1] = input[2] + input[3];
   buffer[2] = input[4] + input[5];
   buffer[3] = input[6] + input[7];
   buffer[4] = C0 * (-input[0] + input[1]);
   buffer[5] = C0 * (-input[2] + input[3]);
   buffer[6] = C0 * (-input[4] + input[5]);
   buffer[7] = C0 * (-input[6] + input[7]);
   buffer[8] = buffer[0] + buffer[1];
   buffer[9] = buffer[2] + buffer[3];
   buffer[10] = buffer[4] + buffer[5];
   buffer[11] = C0 * (-buffer[0] + buffer[1]);
   buffer[12] = buffer[6] + buffer[7];
   buffer[13] = C0 * (-buffer[2] + buffer[3]);
   output[0] = buffer[8] + buffer[9];
   output[1] = buffer[10] + buffer[12];
   output[2] = buffer[11] + buffer[13];
   output[3] = C0 * (-buffer[8] + buffer[9]);
   return output;
}

template<typename T>
std::array<T, squareSize<3, 2>> dgSynthesize_3D_O2(std::array<T, triangleSize<3, 2>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 14> buffer;
   std::array<T, 8> output;
   static U const C0 = U(5.7735026918962573e-01);
   buffer[0] = -C0 * input[3] + input[0];
   buffer[1] = C0 * input[3] + input[0];
   buffer[2] = input[1];
   buffer[3] = input[2];
   buffer[4] = input[1];
   buffer[5] = input[2];
   buffer[6] = -C0 * buffer[3] + buffer[0];
   buffer[7] = C0 * buffer[3] + buffer[0];
   buffer[8] = -C0 * buffer[5] + buffer[1];
   buffer[9] = C0 * buffer[5] + buffer[1];
   buffer[10] = buffer[2];
   buffer[11] = buffer[2];
   buffer[12] = buffer[4];
   buffer[13] = buffer[4];
   output[0] = -C0 * buffer[10] + buffer[6];
   output[1] = C0 * buffer[10] + buffer[6];
   output[2] = -C0 * buffer[11] + buffer[7];
   output[3] = C0 * buffer[11] + buffer[7];
   output[4] = -C0 * buffer[12] + buffer[8];
   output[5] = C0 * buffer[12] + buffer[8];
   output[6] = -C0 * buffer[13] + buffer[9];
   output[7] = C0 * buffer[13] + buffer[9];
   return output;
}

template<typename T>
std::array<T, triangleSize<3, 2>> dgMassInverse_3D_O2(std::array<T, triangleSize<3, 2>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(1.2500000000000000e-01);
   static U const C1 = U(3.7500000000000000e-01);
   std::array<T, triangleSize<3, 2>> output;
   output[0] = C0 * input[0];
   output[1] = C1 * input[1];
   output[2] = C1 * input[2];
   output[3] = C1 * input[3];
   return output;
}

template<typename T>
std::array<T, triangleSize<3, 2>> dgStiffness_3D_O2(int dimension, std::array<T, triangleSize<3, 2>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(2.0000000000000000e+00);
   std::array<T, triangleSize<3, 2>> output;
   if(dimension == 0) {
      output[0] = U(0.0);
      output[1] = U(0.0);
      output[2] = U(0.0);
      output[3] = C0 * input[0];
   } else if(dimension == 1) {
      output[0] = U(0.0);
      output[1] = U(0.0);
      output[2] = C0 * input[0];
      output[3] = U(0.0);
   } else /*if(dimension == 2)*/ {
      output[0] = U(0.0);
      output[1] = C0 * input[0];
      output[2] = U(0.0);
      output[3] = U(0.0);
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<2, 2>> dgTrace_3D_O2(int face, std::array<T, triangleSize<3, 2>> const& input) {
   std::array<T, triangleSize<2, 2>> output;
   if(face == 0) {
      output[0] = input[0] - input[3];
      output[1] = input[1];
      output[2] = input[2];
   } else if(face == 1) {
      output[0] = input[0] + input[3];
      output[1] = input[1];
      output[2] = input[2];
   } else if(face == 2) {
      output[0] = input[0] - input[2];
      output[1] = input[1];
      output[2] = input[3];
   } else if(face == 3) {
      output[0] = input[0] + input[2];
      output[1] = input[1];
      output[2] = input[3];
   } else if(face == 4) {
      output[0] = input[0] - input[1];
      output[1] = input[2];
      output[2] = input[3];
   } else /*if(face == 5)*/ {
      output[0] = input[0] + input[1];
      output[1] = input[2];
      output[2] = input[3];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<3, 2>> dgTraceInverse_3D_O2(int face, std::array<T, triangleSize<2, 2>> const& input) {
   std::array<T, triangleSize<3, 2>> output;
   if(face == 0) {
      output[0] = input[0];
      output[1] = input[1];
      output[2] = input[2];
      output[3] = -input[0];
   } else if(face == 1) {
      output[0] = input[0];
      output[1] = input[1];
      output[2] = input[2];
      output[3] = input[0];
   } else if(face == 2) {
      output[0] = input[0];
      output[1] = input[1];
      output[2] = -input[0];
      output[3] = input[2];
   } else if(face == 3) {
      output[0] = input[0];
      output[1] = input[1];
      output[2] = input[0];
      output[3] = input[2];
   } else if(face == 4) {
      output[0] = input[0];
      output[1] = -input[0];
      output[2] = input[1];
      output[3] = input[2];
   } else /*if(face == 5)*/ {
      output[0] = input[0];
      output[1] = input[0];
      output[2] = input[1];
      output[3] = input[2];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<3, 3>> dgAnalyze_3D_O3(std::array<T, squareSize<3, 3>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 45> buffer;
   std::array<T, 10> output;
   static U const C0 = U(5.5555555555555558e-01);
   static U const C1 = U(8.8888888888888884e-01);
   static U const C2 = U(4.3033148291193524e-01);
   static U const C3 = U(2.2222222222222221e-01);
   static U const C4 = U(4.4444444444444442e-01);
   buffer[0] = C0 * (input[0] + input[2]) + C1 * input[1];
   buffer[1] = C0 * (input[3] + input[5]) + C1 * input[4];
   buffer[2] = C0 * (input[6] + input[8]) + C1 * input[7];
   buffer[3] = C0 * (input[9] + input[11]) + C1 * input[10];
   buffer[4] = C0 * (input[12] + input[14]) + C1 * input[13];
   buffer[5] = C0 * (input[15] + input[17]) + C1 * input[16];
   buffer[6] = C0 * (input[18] + input[20]) + C1 * input[19];
   buffer[7] = C0 * (input[21] + input[23]) + C1 * input[22];
   buffer[8] = C0 * (input[24] + input[26]) + C1 * input[25];
   buffer[9] = C2 * (-input[0] + input[2]);
   buffer[10] = C2 * (-input[3] + input[5]);
   buffer[11] = C2 * (-input[6] + input[8]);
   buffer[12] = C2 * (-input[9] + input[11]);
   buffer[13] = C2 * (-input[12] + input[14]);
   buffer[14] = C2 * (-input[15] + input[17]);
   buffer[15] = C2 * (-input[18] + input[20]);
   buffer[16] = C2 * (-input[21] + input[23]);
   buffer[17] = C2 * (-input[24] + input[26]);
   buffer[18] = C3 * (input[0] + input[2]) - C4 * input[1];
   buffer[19] = C3 * (input[3] + input[5]) - C4 * input[4];
   buffer[20] = C3 * (input[6] + input[8]) - C4 * input[7];
   buffer[21] = C3 * (input[9] + input[11]) - C4 * input[10];
   buffer[22] = C3 * (input[12] + input[14]) - C4 * input[13];
   buffer[23] = C3 * (input[15] + input[17]) - C4 * input[16];
   buffer[24] = C3 * (input[18] + input[20]) - C4 * input[19];
   buffer[25] = C3 * (input[21] + input[23]) - C4 * input[22];
   buffer[26] = C3 * (input[24] + input[26]) - C4 * input[25];
   buffer[27] = C0 * (buffer[0] + buffer[2]) + C1 * buffer[1];
   buffer[28] = C0 * (buffer[3] + buffer[5]) + C1 * buffer[4];
   buffer[29] = C0 * (buffer[6] + buffer[8]) + C1 * buffer[7];
   buffer[30] = C0 * (buffer[9] + buffer[11]) + C1 * buffer[10];
   buffer[31] = C2 * (-buffer[0] + buffer[2]);
   buffer[32] = C0 * (buffer[12] + buffer[14]) + C1 * buffer[13];
   buffer[33] = C2 * (-buffer[3] + buffer[5]);
   buffer[34] = C0 * (buffer[15] + buffer[17]) + C1 * buffer[16];
   buffer[35] = C2 * (-buffer[6] + buffer[8]);
   buffer[36] = C0 * (buffer[18] + buffer[20]) + C1 * buffer[19];
   buffer[37] = C2 * (-buffer[9] + buffer[11]);
   buffer[38] = C3 * (buffer[0] + buffer[2]) - C4 * buffer[1];
   buffer[39] = C0 * (buffer[21] + buffer[23]) + C1 * buffer[22];
   buffer[40] = C2 * (-buffer[12] + buffer[14]);
   buffer[41] = C3 * (buffer[3] + buffer[5]) - C4 * buffer[4];
   buffer[42] = C0 * (buffer[24] + buffer[26]) + C1 * buffer[25];
   buffer[43] = C2 * (-buffer[15] + buffer[17]);
   buffer[44] = C3 * (buffer[6] + buffer[8]) - C4 * buffer[7];
   output[0] = C0 * (buffer[27] + buffer[29]) + C1 * buffer[28];
   output[1] = C0 * (buffer[30] + buffer[34]) + C1 * buffer[32];
   output[2] = C0 * (buffer[31] + buffer[35]) + C1 * buffer[33];
   output[3] = C2 * (-buffer[27] + buffer[29]);
   output[4] = C0 * (buffer[36] + buffer[42]) + C1 * buffer[39];
   output[5] = C0 * (buffer[37] + buffer[43]) + C1 * buffer[40];
   output[6] = C0 * (buffer[38] + buffer[44]) + C1 * buffer[41];
   output[7] = C2 * (-buffer[30] + buffer[34]);
   output[8] = C2 * (-buffer[31] + buffer[35]);
   output[9] = C3 * (buffer[27] + buffer[29]) - C4 * buffer[28];
   return output;
}

template<typename T>
std::array<T, squareSize<3, 3>> dgSynthesize_3D_O3(std::array<T, triangleSize<3, 3>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 45> buffer;
   std::array<T, 27> output;
   static U const C0 = U(4.0000000000000002e-01);
   static U const C1 = U(7.7459666924148340e-01);
   static U const C2 = U(5.0000000000000000e-01);
   buffer[0] = C0 * input[9] - C1 * input[3] + input[0];
   buffer[1] = -C2 * input[9] + input[0];
   buffer[2] = C0 * input[9] + C1 * input[3] + input[0];
   buffer[3] = -C1 * input[7] + input[1];
   buffer[4] = -C1 * input[8] + input[2];
   buffer[5] = input[1];
   buffer[6] = input[2];
   buffer[7] = C1 * input[7] + input[1];
   buffer[8] = C1 * input[8] + input[2];
   buffer[9] = input[4];
   buffer[10] = input[5];
   buffer[11] = input[6];
   buffer[12] = input[4];
   buffer[13] = input[5];
   buffer[14] = input[6];
   buffer[15] = input[4];
   buffer[16] = input[5];
   buffer[17] = input[6];
   buffer[18] = C0 * buffer[11] - C1 * buffer[4] + buffer[0];
   buffer[19] = -C2 * buffer[11] + buffer[0];
   buffer[20] = C0 * buffer[11] + C1 * buffer[4] + buffer[0];
   buffer[21] = C0 * buffer[14] - C1 * buffer[6] + buffer[1];
   buffer[22] = -C2 * buffer[14] + buffer[1];
   buffer[23] = C0 * buffer[14] + C1 * buffer[6] + buffer[1];
   buffer[24] = C0 * buffer[17] - C1 * buffer[8] + buffer[2];
   buffer[25] = -C2 * buffer[17] + buffer[2];
   buffer[26] = C0 * buffer[17] + C1 * buffer[8] + buffer[2];
   buffer[27] = -C1 * buffer[10] + buffer[3];
   buffer[28] = buffer[3];
   buffer[29] = C1 * buffer[10] + buffer[3];
   buffer[30] = -C1 * buffer[13] + buffer[5];
   buffer[31] = buffer[5];
   buffer[32] = C1 * buffer[13] + buffer[5];
   buffer[33] = -C1 * buffer[16] + buffer[7];
   buffer[34] = buffer[7];
   buffer[35] = C1 * buffer[16] + buffer[7];
   buffer[36] = buffer[9];
   buffer[37] = buffer[9];
   buffer[38] = buffer[9];
   buffer[39] = buffer[12];
   buffer[40] = buffer[12];
   buffer[41] = buffer[12];
   buffer[42] = buffer[15];
   buffer[43] = buffer[15];
   buffer[44] = buffer[15];
   output[0] = C0 * buffer[36] - C1 * buffer[27] + buffer[18];
   output[1] = -C2 * buffer[36] + buffer[18];
   output[2] = C0 * buffer[36] + C1 * buffer[27] + buffer[18];
   output[3] = C0 * buffer[37] - C1 * buffer[28] + buffer[19];
   output[4] = -C2 * buffer[37] + buffer[19];
   output[5] = C0 * buffer[37] + C1 * buffer[28] + buffer[19];
   output[6] = C0 * buffer[38] - C1 * buffer[29] + buffer[20];
   output[7] = -C2 * buffer[38] + buffer[20];
   output[8] = C0 * buffer[38] + C1 * buffer[29] + buffer[20];
   output[9] = C0 * buffer[39] - C1 * buffer[30] + buffer[21];
   output[10] = -C2 * buffer[39] + buffer[21];
   output[11] = C0 * buffer[39] + C1 * buffer[30] + buffer[21];
   output[12] = C0 * buffer[40] - C1 * buffer[31] + buffer[22];
   output[13] = -C2 * buffer[40] + buffer[22];
   output[14] = C0 * buffer[40] + C1 * buffer[31] + buffer[22];
   output[15] = C0 * buffer[41] - C1 * buffer[32] + buffer[23];
   output[16] = -C2 * buffer[41] + buffer[23];
   output[17] = C0 * buffer[41] + C1 * buffer[32] + buffer[23];
   output[18] = C0 * buffer[42] - C1 * buffer[33] + buffer[24];
   output[19] = -C2 * buffer[42] + buffer[24];
   output[20] = C0 * buffer[42] + C1 * buffer[33] + buffer[24];
   output[21] = C0 * buffer[43] - C1 * buffer[34] + buffer[25];
   output[22] = -C2 * buffer[43] + buffer[25];
   output[23] = C0 * buffer[43] + C1 * buffer[34] + buffer[25];
   output[24] = C0 * buffer[44] - C1 * buffer[35] + buffer[26];
   output[25] = -C2 * buffer[44] + buffer[26];
   output[26] = C0 * buffer[44] + C1 * buffer[35] + buffer[26];
   return output;
}

template<typename T>
std::array<T, triangleSize<3, 3>> dgMassInverse_3D_O3(std::array<T, triangleSize<3, 3>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(1.2500000000000000e-01);
   static U const C1 = U(3.7500000000000000e-01);
   static U const C2 = U(6.2500000000000000e-01);
   static U const C3 = U(1.1250000000000000e+00);
   std::array<T, triangleSize<3, 3>> output;
   output[0] = C0 * input[0];
   output[1] = C1 * input[1];
   output[2] = C1 * input[2];
   output[3] = C1 * input[3];
   output[4] = C2 * input[4];
   output[5] = C3 * input[5];
   output[6] = C2 * input[6];
   output[7] = C3 * input[7];
   output[8] = C3 * input[8];
   output[9] = C2 * input[9];
   return output;
}

template<typename T>
std::array<T, triangleSize<3, 3>> dgStiffness_3D_O3(int dimension, std::array<T, triangleSize<3, 3>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(2.0000000000000000e+00);
   std::array<T, triangleSize<3, 3>> output;
   if(dimension == 0) {
      output[0] = U(0.0);
      output[1] = U(0.0);
      output[2] = U(0.0);
      output[3] = C0 * input[0];
      output[4] = U(0.0);
      output[5] = U(0.0);
      output[6] = U(0.0);
      output[7] = C0 * input[1];
      output[8] = C0 * input[2];
      output[9] = C0 * input[3];
   } else if(dimension == 1) {
      output[0] = U(0.0);
      output[1] = U(0.0);
      output[2] = C0 * input[0];
      output[3] = U(0.0);
      output[4] = U(0.0);
      output[5] = C0 * input[1];
      output[6] = C0 * input[2];
      output[7] = U(0.0);
      output[8] = C0 * input[3];
      output[9] = U(0.0);
   } else /*if(dimension == 2)*/ {
      output[0] = U(0.0);
      output[1] = C0 * input[0];
      output[2] = U(0.0);
      output[3] = U(0.0);
      output[4] = C0 * input[1];
      output[5] = C0 * input[2];
      output[6] = U(0.0);
      output[7] = C0 * input[3];
      output[8] = U(0.0);
      output[9] = U(0.0);
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<2, 3>> dgTrace_3D_O3(int face, std::array<T, triangleSize<3, 3>> const& input) {
   std::array<T, triangleSize<2, 3>> output;
   if(face == 0) {
      output[0] = input[0] - input[3] + input[9];
      output[1] = input[1] - input[7];
      output[2] = input[2] - input[8];
      output[3] = input[4];
      output[4] = input[5];
      output[5] = input[6];
   } else if(face == 1) {
      output[0] = input[0] + input[3] + input[9];
      output[1] = input[1] + input[7];
      output[2] = input[2] + input[8];
      output[3] = input[4];
      output[4] = input[5];
      output[5] = input[6];
   } else if(face == 2) {
      output[0] = input[0] - input[2] + input[6];
      output[1] = input[1] - input[5];
      output[2] = input[3] - input[8];
      output[3] = input[4];
      output[4] = input[7];
      output[5] = input[9];
   } else if(face == 3) {
      output[0] = input[0] + input[2] + input[6];
      output[1] = input[1] + input[5];
      output[2] = input[3] + input[8];
      output[3] = input[4];
      output[4] = input[7];
      output[5] = input[9];
   } else if(face == 4) {
      output[0] = input[0] - input[1] + input[4];
      output[1] = input[2] - input[5];
      output[2] = input[3] - input[7];
      output[3] = input[6];
      output[4] = input[8];
      output[5] = input[9];
   } else /*if(face == 5)*/ {
      output[0] = input[0] + input[1] + input[4];
      output[1] = input[2] + input[5];
      output[2] = input[3] + input[7];
      output[3] = input[6];
      output[4] = input[8];
      output[5] = input[9];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<3, 3>> dgTraceInverse_3D_O3(int face, std::array<T, triangleSize<2, 3>> const& input) {
   std::array<T, triangleSize<3, 3>> output;
   if(face == 0) {
      output[0] = input[0];
      output[1] = input[1];
      output[2] = input[2];
      output[3] = -input[0];
      output[4] = input[3];
      output[5] = input[4];
      output[6] = input[5];
      output[7] = -input[1];
      output[8] = -input[2];
      output[9] = input[0];
   } else if(face == 1) {
      output[0] = input[0];
      output[1] = input[1];
      output[2] = input[2];
      output[3] = input[0];
      output[4] = input[3];
      output[5] = input[4];
      output[6] = input[5];
      output[7] = input[1];
      output[8] = input[2];
      output[9] = input[0];
   } else if(face == 2) {
      output[0] = input[0];
      output[1] = input[1];
      output[2] = -input[0];
      output[3] = input[2];
      output[4] = input[3];
      output[5] = -input[1];
      output[6] = input[0];
      output[7] = input[4];
      output[8] = -input[2];
      output[9] = input[5];
   } else if(face == 3) {
      output[0] = input[0];
      output[1] = input[1];
      output[2] = input[0];
      output[3] = input[2];
      output[4] = input[3];
      output[5] = input[1];
      output[6] = input[0];
      output[7] = input[4];
      output[8] = input[2];
      output[9] = input[5];
   } else if(face == 4) {
      output[0] = input[0];
      output[1] = -input[0];
      output[2] = input[1];
      output[3] = input[2];
      output[4] = input[0];
      output[5] = -input[1];
      output[6] = input[3];
      output[7] = -input[2];
      output[8] = input[4];
      output[9] = input[5];
   } else /*if(face == 5)*/ {
      output[0] = input[0];
      output[1] = input[0];
      output[2] = input[1];
      output[3] = input[2];
      output[4] = input[0];
      output[5] = input[1];
      output[6] = input[3];
      output[7] = input[2];
      output[8] = input[4];
      output[9] = input[5];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<3, 4>> dgAnalyze_3D_O4(std::array<T, squareSize<3, 4>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 104> buffer;
   std::array<T, 20> output;
   static U const C0 = U(3.4785484513745385e-01);
   static U const C1 = U(6.5214515486254609e-01);
   static U const C2 = U(2.2171699031897615e-01);
   static U const C3 = U(2.9955043831178735e-01);
   static U const C4 = U(2.1300321680756459e-01);
   static U const C5 = U(1.0600771525769923e-01);
   static U const C6 = U(2.6850642010792947e-01);
   buffer[0] = C0 * (input[0] + input[3]) + C1 * (input[1] + input[2]);
   buffer[1] = C0 * (input[4] + input[7]) + C1 * (input[5] + input[6]);
   buffer[2] = C0 * (input[8] + input[11]) + C1 * (input[9] + input[10]);
   buffer[3] = C0 * (input[12] + input[15]) + C1 * (input[13] + input[14]);
   buffer[4] = C0 * (input[16] + input[19]) + C1 * (input[17] + input[18]);
   buffer[5] = C0 * (input[20] + input[23]) + C1 * (input[21] + input[22]);
   buffer[6] = C0 * (input[24] + input[27]) + C1 * (input[25] + input[26]);
   buffer[7] = C0 * (input[28] + input[31]) + C1 * (input[29] + input[30]);
   buffer[8] = C0 * (input[32] + input[35]) + C1 * (input[33] + input[34]);
   buffer[9] = C0 * (input[36] + input[39]) + C1 * (input[37] + input[38]);
   buffer[10] = C0 * (input[40] + input[43]) + C1 * (input[41] + input[42]);
   buffer[11] = C0 * (input[44] + input[47]) + C1 * (input[45] + input[46]);
   buffer[12] = C0 * (input[48] + input[51]) + C1 * (input[49] + input[50]);
   buffer[13] = C0 * (input[52] + input[55]) + C1 * (input[53] + input[54]);
   buffer[14] = C0 * (input[56] + input[59]) + C1 * (input[57] + input[58]);
   buffer[15] = C0 * (input[60] + input[63]) + C1 * (input[61] + input[62]);
   buffer[16] = C2 * (-input[1] + input[2]) + C3 * (-input[0] + input[3]);
   buffer[17] = C2 * (-input[5] + input[6]) + C3 * (-input[4] + input[7]);
   buffer[18] = C2 * (-input[9] + input[10]) + C3 * (-input[8] + input[11]);
   buffer[19] = C2 * (-input[13] + input[14]) + C3 * (-input[12] + input[15]);
   buffer[20] = C2 * (-input[17] + input[18]) + C3 * (-input[16] + input[19]);
   buffer[21] = C2 * (-input[21] + input[22]) + C3 * (-input[20] + input[23]);
   buffer[22] = C2 * (-input[25] + input[26]) + C3 * (-input[24] + input[27]);
   buffer[23] = C2 * (-input[29] + input[30]) + C3 * (-input[28] + input[31]);
   buffer[24] = C2 * (-input[33] + input[34]) + C3 * (-input[32] + input[35]);
   buffer[25] = C2 * (-input[37] + input[38]) + C3 * (-input[36] + input[39]);
   buffer[26] = C2 * (-input[41] + input[42]) + C3 * (-input[40] + input[43]);
   buffer[27] = C2 * (-input[45] + input[46]) + C3 * (-input[44] + input[47]);
   buffer[28] = C2 * (-input[49] + input[50]) + C3 * (-input[48] + input[51]);
   buffer[29] = C2 * (-input[53] + input[54]) + C3 * (-input[52] + input[55]);
   buffer[30] = C2 * (-input[57] + input[58]) + C3 * (-input[56] + input[59]);
   buffer[31] = C2 * (-input[61] + input[62]) + C3 * (-input[60] + input[63]);
   buffer[32] = C4 * (input[0] - input[1] - input[2] + input[3]);
   buffer[33] = C4 * (input[4] - input[5] - input[6] + input[7]);
   buffer[34] = C4 * (input[8] - input[9] - input[10] + input[11]);
   buffer[35] = C4 * (input[12] - input[13] - input[14] + input[15]);
   buffer[36] = C4 * (input[16] - input[17] - input[18] + input[19]);
   buffer[37] = C4 * (input[20] - input[21] - input[22] + input[23]);
   buffer[38] = C4 * (input[24] - input[25] - input[26] + input[27]);
   buffer[39] = C4 * (input[28] - input[29] - input[30] + input[31]);
   buffer[40] = C4 * (input[32] - input[33] - input[34] + input[35]);
   buffer[41] = C4 * (input[36] - input[37] - input[38] + input[39]);
   buffer[42] = C4 * (input[40] - input[41] - input[42] + input[43]);
   buffer[43] = C4 * (input[44] - input[45] - input[46] + input[47]);
   buffer[44] = C4 * (input[48] - input[49] - input[50] + input[51]);
   buffer[45] = C4 * (input[52] - input[53] - input[54] + input[55]);
   buffer[46] = C4 * (input[56] - input[57] - input[58] + input[59]);
   buffer[47] = C4 * (input[60] - input[61] - input[62] + input[63]);
   buffer[48] = C5 * (-input[0] + input[3]) + C6 * (input[1] - input[2]);
   buffer[49] = C5 * (-input[4] + input[7]) + C6 * (input[5] - input[6]);
   buffer[50] = C5 * (-input[8] + input[11]) + C6 * (input[9] - input[10]);
   buffer[51] = C5 * (-input[12] + input[15]) + C6 * (input[13] - input[14]);
   buffer[52] = C5 * (-input[16] + input[19]) + C6 * (input[17] - input[18]);
   buffer[53] = C5 * (-input[20] + input[23]) + C6 * (input[21] - input[22]);
   buffer[54] = C5 * (-input[24] + input[27]) + C6 * (input[25] - input[26]);
   buffer[55] = C5 * (-input[28] + input[31]) + C6 * (input[29] - input[30]);
   buffer[56] = C5 * (-input[32] + input[35]) + C6 * (input[33] - input[34]);
   buffer[57] = C5 * (-input[36] + input[39]) + C6 * (input[37] - input[38]);
   buffer[58] = C5 * (-input[40] + input[43]) + C6 * (input[41] - input[42]);
   buffer[59] = C5 * (-input[44] + input[47]) + C6 * (input[45] - input[46]);
   buffer[60] = C5 * (-input[48] + input[51]) + C6 * (input[49] - input[50]);
   buffer[61] = C5 * (-input[52] + input[55]) + C6 * (input[53] - input[54]);
   buffer[62] = C5 * (-input[56] + input[59]) + C6 * (input[57] - input[58]);
   buffer[63] = C5 * (-input[60] + input[63]) + C6 * (input[61] - input[62]);
   buffer[64] = C0 * (buffer[0] + buffer[3]) + C1 * (buffer[1] + buffer[2]);
   buffer[65] = C0 * (buffer[4] + buffer[7]) + C1 * (buffer[5] + buffer[6]);
   buffer[66] = C0 * (buffer[8] + buffer[11]) + C1 * (buffer[9] + buffer[10]);
   buffer[67] = C0 * (buffer[12] + buffer[15]) + C1 * (buffer[13] + buffer[14]);
   buffer[68] = C0 * (buffer[16] + buffer[19]) + C1 * (buffer[17] + buffer[18]);
   buffer[69] = C2 * (-buffer[1] + buffer[2]) + C3 * (-buffer[0] + buffer[3]);
   buffer[70] = C0 * (buffer[20] + buffer[23]) + C1 * (buffer[21] + buffer[22]);
   buffer[71] = C2 * (-buffer[5] + buffer[6]) + C3 * (-buffer[4] + buffer[7]);
   buffer[72] = C0 * (buffer[24] + buffer[27]) + C1 * (buffer[25] + buffer[26]);
   buffer[73] = C2 * (-buffer[9] + buffer[10]) + C3 * (-buffer[8] + buffer[11]);
   buffer[74] = C0 * (buffer[28] + buffer[31]) + C1 * (buffer[29] + buffer[30]);
   buffer[75] = C2 * (-buffer[13] + buffer[14]) + C3 * (-buffer[12] + buffer[15]);
   buffer[76] = C0 * (buffer[32] + buffer[35]) + C1 * (buffer[33] + buffer[34]);
   buffer[77] = C2 * (-buffer[17] + buffer[18]) + C3 * (-buffer[16] + buffer[19]);
   buffer[78] = C4 * (buffer[0] - buffer[1] - buffer[2] + buffer[3]);
   buffer[79] = C0 * (buffer[36] + buffer[39]) + C1 * (buffer[37] + buffer[38]);
   buffer[80] = C2 * (-buffer[21] + buffer[22]) + C3 * (-buffer[20] + buffer[23]);
   buffer[81] = C4 * (buffer[4] - buffer[5] - buffer[6] + buffer[7]);
   buffer[82] = C0 * (buffer[40] + buffer[43]) + C1 * (buffer[41] + buffer[42]);
   buffer[83] = C2 * (-buffer[25] + buffer[26]) + C3 * (-buffer[24] + buffer[27]);
   buffer[84] = C4 * (buffer[8] - buffer[9] - buffer[10] + buffer[11]);
   buffer[85] = C0 * (buffer[44] + buffer[47]) + C1 * (buffer[45] + buffer[46]);
   buffer[86] = C2 * (-buffer[29] + buffer[30]) + C3 * (-buffer[28] + buffer[31]);
   buffer[87] = C4 * (buffer[12] - buffer[13] - buffer[14] + buffer[15]);
   buffer[88] = C0 * (buffer[48] + buffer[51]) + C1 * (buffer[49] + buffer[50]);
   buffer[89] = C2 * (-buffer[33] + buffer[34]) + C3 * (-buffer[32] + buffer[35]);
   buffer[90] = C4 * (buffer[16] - buffer[17] - buffer[18] + buffer[19]);
   buffer[91] = C5 * (-buffer[0] + buffer[3]) + C6 * (buffer[1] - buffer[2]);
   buffer[92] = C0 * (buffer[52] + buffer[55]) + C1 * (buffer[53] + buffer[54]);
   buffer[93] = C2 * (-buffer[37] + buffer[38]) + C3 * (-buffer[36] + buffer[39]);
   buffer[94] = C4 * (buffer[20] - buffer[21] - buffer[22] + buffer[23]);
   buffer[95] = C5 * (-buffer[4] + buffer[7]) + C6 * (buffer[5] - buffer[6]);
   buffer[96] = C0 * (buffer[56] + buffer[59]) + C1 * (buffer[57] + buffer[58]);
   buffer[97] = C2 * (-buffer[41] + buffer[42]) + C3 * (-buffer[40] + buffer[43]);
   buffer[98] = C4 * (buffer[24] - buffer[25] - buffer[26] + buffer[27]);
   buffer[99] = C5 * (-buffer[8] + buffer[11]) + C6 * (buffer[9] - buffer[10]);
   buffer[100] = C0 * (buffer[60] + buffer[63]) + C1 * (buffer[61] + buffer[62]);
   buffer[101] = C2 * (-buffer[45] + buffer[46]) + C3 * (-buffer[44] + buffer[47]);
   buffer[102] = C4 * (buffer[28] - buffer[29] - buffer[30] + buffer[31]);
   buffer[103] = C5 * (-buffer[12] + buffer[15]) + C6 * (buffer[13] - buffer[14]);
   output[0] = C0 * (buffer[64] + buffer[67]) + C1 * (buffer[65] + buffer[66]);
   output[1] = C0 * (buffer[68] + buffer[74]) + C1 * (buffer[70] + buffer[72]);
   output[2] = C0 * (buffer[69] + buffer[75]) + C1 * (buffer[71] + buffer[73]);
   output[3] = C2 * (-buffer[65] + buffer[66]) + C3 * (-buffer[64] + buffer[67]);
   output[4] = C0 * (buffer[76] + buffer[85]) + C1 * (buffer[79] + buffer[82]);
   output[5] = C0 * (buffer[77] + buffer[86]) + C1 * (buffer[80] + buffer[83]);
   output[6] = C0 * (buffer[78] + buffer[87]) + C1 * (buffer[81] + buffer[84]);
   output[7] = C2 * (-buffer[70] + buffer[72]) + C3 * (-buffer[68] + buffer[74]);
   output[8] = C2 * (-buffer[71] + buffer[73]) + C3 * (-buffer[69] + buffer[75]);
   output[9] = C4 * (buffer[64] - buffer[65] - buffer[66] + buffer[67]);
   output[10] = C0 * (buffer[88] + buffer[100]) + C1 * (buffer[92] + buffer[96]);
   output[11] = C0 * (buffer[89] + buffer[101]) + C1 * (buffer[93] + buffer[97]);
   output[12] = C0 * (buffer[90] + buffer[102]) + C1 * (buffer[94] + buffer[98]);
   output[13] = C0 * (buffer[91] + buffer[103]) + C1 * (buffer[95] + buffer[99]);
   output[14] = C2 * (-buffer[79] + buffer[82]) + C3 * (-buffer[76] + buffer[85]);
   output[15] = C2 * (-buffer[80] + buffer[83]) + C3 * (-buffer[77] + buffer[86]);
   output[16] = C2 * (-buffer[81] + buffer[84]) + C3 * (-buffer[78] + buffer[87]);
   output[17] = C4 * (buffer[68] - buffer[70] - buffer[72] + buffer[74]);
   output[18] = C4 * (buffer[69] - buffer[71] - buffer[73] + buffer[75]);
   output[19] = C5 * (-buffer[64] + buffer[67]) + C6 * (buffer[65] - buffer[66]);
   return output;
}

template<typename T>
std::array<T, squareSize<3, 4>> dgSynthesize_3D_O4(std::array<T, triangleSize<3, 4>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 104> buffer;
   std::array<T, 64> output;
   static U const C0 = U(3.0474698495520619e-01);
   static U const C1 = U(6.1233362071871378e-01);
   static U const C2 = U(8.6113631159405257e-01);
   static U const C3 = U(3.2661933500442808e-01);
   static U const C4 = U(3.3998104358485626e-01);
   static U const C5 = U(4.1172799967289958e-01);
   buffer[0] = -C0 * input[19] + C1 * input[9] - C2 * input[3] + input[0];
   buffer[1] = -C3 * input[9] - C4 * input[3] + C5 * input[19] + input[0];
   buffer[2] = -C3 * input[9] + C4 * input[3] - C5 * input[19] + input[0];
   buffer[3] = C0 * input[19] + C1 * input[9] + C2 * input[3] + input[0];
   buffer[4] = C1 * input[17] - C2 * input[7] + input[1];
   buffer[5] = C1 * input[18] - C2 * input[8] + input[2];
   buffer[6] = -C3 * input[17] - C4 * input[7] + input[1];
   buffer[7] = -C3 * input[18] - C4 * input[8] + input[2];
   buffer[8] = -C3 * input[17] + C4 * input[7] + input[1];
   buffer[9] = -C3 * input[18] + C4 * input[8] + input[2];
   buffer[10] = C1 * input[17] + C2 * input[7] + input[1];
   buffer[11] = C1 * input[18] + C2 * input[8] + input[2];
   buffer[12] = -C2 * input[14] + input[4];
   buffer[13] = -C2 * input[15] + input[5];
   buffer[14] = -C2 * input[16] + input[6];
   buffer[15] = -C4 * input[14] + input[4];
   buffer[16] = -C4 * input[15] + input[5];
   buffer[17] = -C4 * input[16] + input[6];
   buffer[18] = C4 * input[14] + input[4];
   buffer[19] = C4 * input[15] + input[5];
   buffer[20] = C4 * input[16] + input[6];
   buffer[21] = C2 * input[14] + input[4];
   buffer[22] = C2 * input[15] + input[5];
   buffer[23] = C2 * input[16] + input[6];
   buffer[24] = input[10];
   buffer[25] = input[11];
   buffer[26] = input[12];
   buffer[27] = input[13];
   buffer[28] = input[10];
   buffer[29] = input[11];
   buffer[30] = input[12];
   buffer[31] = input[13];
   buffer[32] = input[10];
   buffer[33] = input[11];
   buffer[34] = input[12];
   buffer[35] = input[13];
   buffer[36] = input[10];
   buffer[37] = input[11];
   buffer[38] = input[12];
   buffer[39] = input[13];
   buffer[40] = -C0 * buffer[27] + C1 * buffer[14] - C2 * buffer[5] + buffer[0];
   buffer[41] = -C3 * buffer[14] - C4 * buffer[5] + C5 * buffer[27] + buffer[0];
   buffer[42] = -C3 * buffer[14] + C4 * buffer[5] - C5 * buffer[27] + buffer[0];
   buffer[43] = C0 * buffer[27] + C1 * buffer[14] + C2 * buffer[5] + buffer[0];
   buffer[44] = -C0 * buffer[31] + C1 * buffer[17] - C2 * buffer[7] + buffer[1];
   buffer[45] = -C3 * buffer[17] - C4 * buffer[7] + C5 * buffer[31] + buffer[1];
   buffer[46] = -C3 * buffer[17] + C4 * buffer[7] - C5 * buffer[31] + buffer[1];
   buffer[47] = C0 * buffer[31] + C1 * buffer[17] + C2 * buffer[7] + buffer[1];
   buffer[48] = -C0 * buffer[35] + C1 * buffer[20] - C2 * buffer[9] + buffer[2];
   buffer[49] = -C3 * buffer[20] - C4 * buffer[9] + C5 * buffer[35] + buffer[2];
   buffer[50] = -C3 * buffer[20] + C4 * buffer[9] - C5 * buffer[35] + buffer[2];
   buffer[51] = C0 * buffer[35] + C1 * buffer[20] + C2 * buffer[9] + buffer[2];
   buffer[52] = -C0 * buffer[39] + C1 * buffer[23] - C2 * buffer[11] + buffer[3];
   buffer[53] = -C3 * buffer[23] - C4 * buffer[11] + C5 * buffer[39] + buffer[3];
   buffer[54] = -C3 * buffer[23] + C4 * buffer[11] - C5 * buffer[39] + buffer[3];
   buffer[55] = C0 * buffer[39] + C1 * buffer[23] + C2 * buffer[11] + buffer[3];
   buffer[56] = C1 * buffer[26] - C2 * buffer[13] + buffer[4];
   buffer[57] = -C3 * buffer[26] - C4 * buffer[13] + buffer[4];
   buffer[58] = -C3 * buffer[26] + C4 * buffer[13] + buffer[4];
   buffer[59] = C1 * buffer[26] + C2 * buffer[13] + buffer[4];
   buffer[60] = C1 * buffer[30] - C2 * buffer[16] + buffer[6];
   buffer[61] = -C3 * buffer[30] - C4 * buffer[16] + buffer[6];
   buffer[62] = -C3 * buffer[30] + C4 * buffer[16] + buffer[6];
   buffer[63] = C1 * buffer[30] + C2 * buffer[16] + buffer[6];
   buffer[64] = C1 * buffer[34] - C2 * buffer[19] + buffer[8];
   buffer[65] = -C3 * buffer[34] - C4 * buffer[19] + buffer[8];
   buffer[66] = -C3 * buffer[34] + C4 * buffer[19] + buffer[8];
   buffer[67] = C1 * buffer[34] + C2 * buffer[19] + buffer[8];
   buffer[68] = C1 * buffer[38] - C2 * buffer[22] + buffer[10];
   buffer[69] = -C3 * buffer[38] - C4 * buffer[22] + buffer[10];
   buffer[70] = -C3 * buffer[38] + C4 * buffer[22] + buffer[10];
   buffer[71] = C1 * buffer[38] + C2 * buffer[22] + buffer[10];
   buffer[72] = -C2 * buffer[25] + buffer[12];
   buffer[73] = -C4 * buffer[25] + buffer[12];
   buffer[74] = C4 * buffer[25] + buffer[12];
   buffer[75] = C2 * buffer[25] + buffer[12];
   buffer[76] = -C2 * buffer[29] + buffer[15];
   buffer[77] = -C4 * buffer[29] + buffer[15];
   buffer[78] = C4 * buffer[29] + buffer[15];
   buffer[79] = C2 * buffer[29] + buffer[15];
   buffer[80] = -C2 * buffer[33] + buffer[18];
   buffer[81] = -C4 * buffer[33] + buffer[18];
   buffer[82] = C4 * buffer[33] + buffer[18];
   buffer[83] = C2 * buffer[33] + buffer[18];
   buffer[84] = -C2 * buffer[37] + buffer[21];
   buffer[85] = -C4 * buffer[37] + buffer[21];
   buffer[86] = C4 * buffer[37] + buffer[21];
   buffer[87] = C2 * buffer[37] + buffer[21];
   buffer[88] = buffer[24];
   buffer[89] = buffer[24];
   buffer[90] = buffer[24];
   buffer[91] = buffer[24];
   buffer[92] = buffer[28];
   buffer[93] = buffer[28];
   buffer[94] = buffer[28];
   buffer[95] = buffer[28];
   buffer[96] = buffer[32];
   buffer[97] = buffer[32];
   buffer[98] = buffer[32];
   buffer[99] = buffer[32];
   buffer[100] = buffer[36];
   buffer[101] = buffer[36];
   buffer[102] = buffer[36];
   buffer[103] = buffer[36];
   output[0] = -C0 * buffer[88] + C1 * buffer[72] - C2 * buffer[56] + buffer[40];
   output[1] = -C3 * buffer[72] - C4 * buffer[56] + C5 * buffer[88] + buffer[40];
   output[2] = -C3 * buffer[72] + C4 * buffer[56] - C5 * buffer[88] + buffer[40];
   output[3] = C0 * buffer[88] + C1 * buffer[72] + C2 * buffer[56] + buffer[40];
   output[4] = -C0 * buffer[89] + C1 * buffer[73] - C2 * buffer[57] + buffer[41];
   output[5] = -C3 * buffer[73] - C4 * buffer[57] + C5 * buffer[89] + buffer[41];
   output[6] = -C3 * buffer[73] + C4 * buffer[57] - C5 * buffer[89] + buffer[41];
   output[7] = C0 * buffer[89] + C1 * buffer[73] + C2 * buffer[57] + buffer[41];
   output[8] = -C0 * buffer[90] + C1 * buffer[74] - C2 * buffer[58] + buffer[42];
   output[9] = -C3 * buffer[74] - C4 * buffer[58] + C5 * buffer[90] + buffer[42];
   output[10] = -C3 * buffer[74] + C4 * buffer[58] - C5 * buffer[90] + buffer[42];
   output[11] = C0 * buffer[90] + C1 * buffer[74] + C2 * buffer[58] + buffer[42];
   output[12] = -C0 * buffer[91] + C1 * buffer[75] - C2 * buffer[59] + buffer[43];
   output[13] = -C3 * buffer[75] - C4 * buffer[59] + C5 * buffer[91] + buffer[43];
   output[14] = -C3 * buffer[75] + C4 * buffer[59] - C5 * buffer[91] + buffer[43];
   output[15] = C0 * buffer[91] + C1 * buffer[75] + C2 * buffer[59] + buffer[43];
   output[16] = -C0 * buffer[92] + C1 * buffer[76] - C2 * buffer[60] + buffer[44];
   output[17] = -C3 * buffer[76] - C4 * buffer[60] + C5 * buffer[92] + buffer[44];
   output[18] = -C3 * buffer[76] + C4 * buffer[60] - C5 * buffer[92] + buffer[44];
   output[19] = C0 * buffer[92] + C1 * buffer[76] + C2 * buffer[60] + buffer[44];
   output[20] = -C0 * buffer[93] + C1 * buffer[77] - C2 * buffer[61] + buffer[45];
   output[21] = -C3 * buffer[77] - C4 * buffer[61] + C5 * buffer[93] + buffer[45];
   output[22] = -C3 * buffer[77] + C4 * buffer[61] - C5 * buffer[93] + buffer[45];
   output[23] = C0 * buffer[93] + C1 * buffer[77] + C2 * buffer[61] + buffer[45];
   output[24] = -C0 * buffer[94] + C1 * buffer[78] - C2 * buffer[62] + buffer[46];
   output[25] = -C3 * buffer[78] - C4 * buffer[62] + C5 * buffer[94] + buffer[46];
   output[26] = -C3 * buffer[78] + C4 * buffer[62] - C5 * buffer[94] + buffer[46];
   output[27] = C0 * buffer[94] + C1 * buffer[78] + C2 * buffer[62] + buffer[46];
   output[28] = -C0 * buffer[95] + C1 * buffer[79] - C2 * buffer[63] + buffer[47];
   output[29] = -C3 * buffer[79] - C4 * buffer[63] + C5 * buffer[95] + buffer[47];
   output[30] = -C3 * buffer[79] + C4 * buffer[63] - C5 * buffer[95] + buffer[47];
   output[31] = C0 * buffer[95] + C1 * buffer[79] + C2 * buffer[63] + buffer[47];
   output[32] = -C0 * buffer[96] + C1 * buffer[80] - C2 * buffer[64] + buffer[48];
   output[33] = -C3 * buffer[80] - C4 * buffer[64] + C5 * buffer[96] + buffer[48];
   output[34] = -C3 * buffer[80] + C4 * buffer[64] - C5 * buffer[96] + buffer[48];
   output[35] = C0 * buffer[96] + C1 * buffer[80] + C2 * buffer[64] + buffer[48];
   output[36] = -C0 * buffer[97] + C1 * buffer[81] - C2 * buffer[65] + buffer[49];
   output[37] = -C3 * buffer[81] - C4 * buffer[65] + C5 * buffer[97] + buffer[49];
   output[38] = -C3 * buffer[81] + C4 * buffer[65] - C5 * buffer[97] + buffer[49];
   output[39] = C0 * buffer[97] + C1 * buffer[81] + C2 * buffer[65] + buffer[49];
   output[40] = -C0 * buffer[98] + C1 * buffer[82] - C2 * buffer[66] + buffer[50];
   output[41] = -C3 * buffer[82] - C4 * buffer[66] + C5 * buffer[98] + buffer[50];
   output[42] = -C3 * buffer[82] + C4 * buffer[66] - C5 * buffer[98] + buffer[50];
   output[43] = C0 * buffer[98] + C1 * buffer[82] + C2 * buffer[66] + buffer[50];
   output[44] = -C0 * buffer[99] + C1 * buffer[83] - C2 * buffer[67] + buffer[51];
   output[45] = -C3 * buffer[83] - C4 * buffer[67] + C5 * buffer[99] + buffer[51];
   output[46] = -C3 * buffer[83] + C4 * buffer[67] - C5 * buffer[99] + buffer[51];
   output[47] = C0 * buffer[99] + C1 * buffer[83] + C2 * buffer[67] + buffer[51];
   output[48] = -C0 * buffer[100] + C1 * buffer[84] - C2 * buffer[68] + buffer[52];
   output[49] = -C3 * buffer[84] - C4 * buffer[68] + C5 * buffer[100] + buffer[52];
   output[50] = -C3 * buffer[84] + C4 * buffer[68] - C5 * buffer[100] + buffer[52];
   output[51] = C0 * buffer[100] + C1 * buffer[84] + C2 * buffer[68] + buffer[52];
   output[52] = -C0 * buffer[101] + C1 * buffer[85] - C2 * buffer[69] + buffer[53];
   output[53] = -C3 * buffer[85] - C4 * buffer[69] + C5 * buffer[101] + buffer[53];
   output[54] = -C3 * buffer[85] + C4 * buffer[69] - C5 * buffer[101] + buffer[53];
   output[55] = C0 * buffer[101] + C1 * buffer[85] + C2 * buffer[69] + buffer[53];
   output[56] = -C0 * buffer[102] + C1 * buffer[86] - C2 * buffer[70] + buffer[54];
   output[57] = -C3 * buffer[86] - C4 * buffer[70] + C5 * buffer[102] + buffer[54];
   output[58] = -C3 * buffer[86] + C4 * buffer[70] - C5 * buffer[102] + buffer[54];
   output[59] = C0 * buffer[102] + C1 * buffer[86] + C2 * buffer[70] + buffer[54];
   output[60] = -C0 * buffer[103] + C1 * buffer[87] - C2 * buffer[71] + buffer[55];
   output[61] = -C3 * buffer[87] - C4 * buffer[71] + C5 * buffer[103] + buffer[55];
   output[62] = -C3 * buffer[87] + C4 * buffer[71] - C5 * buffer[103] + buffer[55];
   output[63] = C0 * buffer[103] + C1 * buffer[87] + C2 * buffer[71] + buffer[55];
   return output;
}

template<typename T>
std::array<T, triangleSize<3, 4>> dgMassInverse_3D_O4(std::array<T, triangleSize<3, 4>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(1.2500000000000003e-01);
   static U const C1 = U(3.7500000000000000e-01);
   static U const C2 = U(6.2500000000000011e-01);
   static U const C3 = U(1.1249999999999998e+00);
   static U const C4 = U(8.7500000000000011e-01);
   static U const C5 = U(1.8750000000000002e+00);
   static U const C6 = U(3.3749999999999991e+00);
   std::array<T, triangleSize<3, 4>> output;
   output[0] = C0 * input[0];
   output[1] = C1 * input[1];
   output[2] = C1 * input[2];
   output[3] = C1 * input[3];
   output[4] = C2 * input[4];
   output[5] = C3 * input[5];
   output[6] = C2 * input[6];
   output[7] = C3 * input[7];
   output[8] = C3 * input[8];
   output[9] = C2 * input[9];
   output[10] = C4 * input[10];
   output[11] = C5 * input[11];
   output[12] = C5 * input[12];
   output[13] = C4 * input[13];
   output[14] = C5 * input[14];
   output[15] = C6 * input[15];
   output[16] = C5 * input[16];
   output[17] = C5 * input[17];
   output[18] = C5 * input[18];
   output[19] = C4 * input[19];
   return output;
}

template<typename T>
std::array<T, triangleSize<3, 4>> dgStiffness_3D_O4(int dimension, std::array<T, triangleSize<3, 4>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(2.0000000000000000e+00);
   std::array<T, triangleSize<3, 4>> output;
   if(dimension == 0) {
      output[0] = U(0.0);
      output[1] = U(0.0);
      output[2] = U(0.0);
      output[3] = C0 * input[0];
      output[4] = U(0.0);
      output[5] = U(0.0);
      output[6] = U(0.0);
      output[7] = C0 * input[1];
      output[8] = C0 * input[2];
      output[9] = C0 * input[3];
      output[10] = U(0.0);
      output[11] = U(0.0);
      output[12] = U(0.0);
      output[13] = U(0.0);
      output[14] = C0 * input[4];
      output[15] = C0 * input[5];
      output[16] = C0 * input[6];
      output[17] = C0 * input[7];
      output[18] = C0 * input[8];
      output[19] = C0 * (input[0] + input[9]);
   } else if(dimension == 1) {
      output[0] = U(0.0);
      output[1] = U(0.0);
      output[2] = C0 * input[0];
      output[3] = U(0.0);
      output[4] = U(0.0);
      output[5] = C0 * input[1];
      output[6] = C0 * input[2];
      output[7] = U(0.0);
      output[8] = C0 * input[3];
      output[9] = U(0.0);
      output[10] = U(0.0);
      output[11] = C0 * input[4];
      output[12] = C0 * input[5];
      output[13] = C0 * (input[0] + input[6]);
      output[14] = U(0.0);
      output[15] = C0 * input[7];
      output[16] = C0 * input[8];
      output[17] = U(0.0);
      output[18] = C0 * input[9];
      output[19] = U(0.0);
   } else /*if(dimension == 2)*/ {
      output[0] = U(0.0);
      output[1] = C0 * input[0];
      output[2] = U(0.0);
      output[3] = U(0.0);
      output[4] = C0 * input[1];
      output[5] = C0 * input[2];
      output[6] = U(0.0);
      output[7] = C0 * input[3];
      output[8] = U(0.0);
      output[9] = U(0.0);
      output[10] = C0 * (input[0] + input[4]);
      output[11] = C0 * input[5];
      output[12] = C0 * input[6];
      output[13] = U(0.0);
      output[14] = C0 * input[7];
      output[15] = C0 * input[8];
      output[16] = U(0.0);
      output[17] = C0 * input[9];
      output[18] = U(0.0);
      output[19] = U(0.0);
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<2, 4>> dgTrace_3D_O4(int face, std::array<T, triangleSize<3, 4>> const& input) {
   std::array<T, triangleSize<2, 4>> output;
   if(face == 0) {
      output[0] = input[0] - input[3] + input[9] - input[19];
      output[1] = input[1] - input[7] + input[17];
      output[2] = input[2] - input[8] + input[18];
      output[3] = input[4] - input[14];
      output[4] = input[5] - input[15];
      output[5] = input[6] - input[16];
      output[6] = input[10];
      output[7] = input[11];
      output[8] = input[12];
      output[9] = input[13];
   } else if(face == 1) {
      output[0] = input[0] + input[3] + input[9] + input[19];
      output[1] = input[1] + input[7] + input[17];
      output[2] = input[2] + input[8] + input[18];
      output[3] = input[4] + input[14];
      output[4] = input[5] + input[15];
      output[5] = input[6] + input[16];
      output[6] = input[10];
      output[7] = input[11];
      output[8] = input[12];
      output[9] = input[13];
   } else if(face == 2) {
      output[0] = input[0] - input[2] + input[6] - input[13];
      output[1] = input[1] - input[5] + input[12];
      output[2] = input[3] - input[8] + input[16];
      output[3] = input[4] - input[11];
      output[4] = input[7] - input[15];
      output[5] = input[9] - input[18];
      output[6] = input[10];
      output[7] = input[14];
      output[8] = input[17];
      output[9] = input[19];
   } else if(face == 3) {
      output[0] = input[0] + input[2] + input[6] + input[13];
      output[1] = input[1] + input[5] + input[12];
      output[2] = input[3] + input[8] + input[16];
      output[3] = input[4] + input[11];
      output[4] = input[7] + input[15];
      output[5] = input[9] + input[18];
      output[6] = input[10];
      output[7] = input[14];
      output[8] = input[17];
      output[9] = input[19];
   } else if(face == 4) {
      output[0] = input[0] - input[1] + input[4] - input[10];
      output[1] = input[2] - input[5] + input[11];
      output[2] = input[3] - input[7] + input[14];
      output[3] = input[6] - input[12];
      output[4] = input[8] - input[15];
      output[5] = input[9] - input[17];
      output[6] = input[13];
      output[7] = input[16];
      output[8] = input[18];
      output[9] = input[19];
   } else /*if(face == 5)*/ {
      output[0] = input[0] + input[1] + input[4] + input[10];
      output[1] = input[2] + input[5] + input[11];
      output[2] = input[3] + input[7] + input[14];
      output[3] = input[6] + input[12];
      output[4] = input[8] + input[15];
      output[5] = input[9] + input[17];
      output[6] = input[13];
      output[7] = input[16];
      output[8] = input[18];
      output[9] = input[19];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<3, 4>> dgTraceInverse_3D_O4(int face, std::array<T, triangleSize<2, 4>> const& input) {
   std::array<T, triangleSize<3, 4>> output;
   if(face == 0) {
      output[0] = input[0];
      output[1] = input[1];
      output[2] = input[2];
      output[3] = -input[0];
      output[4] = input[3];
      output[5] = input[4];
      output[6] = input[5];
      output[7] = -input[1];
      output[8] = -input[2];
      output[9] = input[0];
      output[10] = input[6];
      output[11] = input[7];
      output[12] = input[8];
      output[13] = input[9];
      output[14] = -input[3];
      output[15] = -input[4];
      output[16] = -input[5];
      output[17] = input[1];
      output[18] = input[2];
      output[19] = -input[0];
   } else if(face == 1) {
      output[0] = input[0];
      output[1] = input[1];
      output[2] = input[2];
      output[3] = input[0];
      output[4] = input[3];
      output[5] = input[4];
      output[6] = input[5];
      output[7] = input[1];
      output[8] = input[2];
      output[9] = input[0];
      output[10] = input[6];
      output[11] = input[7];
      output[12] = input[8];
      output[13] = input[9];
      output[14] = input[3];
      output[15] = input[4];
      output[16] = input[5];
      output[17] = input[1];
      output[18] = input[2];
      output[19] = input[0];
   } else if(face == 2) {
      output[0] = input[0];
      output[1] = input[1];
      output[2] = -input[0];
      output[3] = input[2];
      output[4] = input[3];
      output[5] = -input[1];
      output[6] = input[0];
      output[7] = input[4];
      output[8] = -input[2];
      output[9] = input[5];
      output[10] = input[6];
      output[11] = -input[3];
      output[12] = input[1];
      output[13] = -input[0];
      output[14] = input[7];
      output[15] = -input[4];
      output[16] = input[2];
      output[17] = input[8];
      output[18] = -input[5];
      output[19] = input[9];
   } else if(face == 3) {
      output[0] = input[0];
      output[1] = input[1];
      output[2] = input[0];
      output[3] = input[2];
      output[4] = input[3];
      output[5] = input[1];
      output[6] = input[0];
      output[7] = input[4];
      output[8] = input[2];
      output[9] = input[5];
      output[10] = input[6];
      output[11] = input[3];
      output[12] = input[1];
      output[13] = input[0];
      output[14] = input[7];
      output[15] = input[4];
      output[16] = input[2];
      output[17] = input[8];
      output[18] = input[5];
      output[19] = input[9];
   } else if(face == 4) {
      output[0] = input[0];
      output[1] = -input[0];
      output[2] = input[1];
      output[3] = input[2];
      output[4] = input[0];
      output[5] = -input[1];
      output[6] = input[3];
      output[7] = -input[2];
      output[8] = input[4];
      output[9] = input[5];
      output[10] = -input[0];
      output[11] = input[1];
      output[12] = -input[3];
      output[13] = input[6];
      output[14] = input[2];
      output[15] = -input[4];
      output[16] = input[7];
      output[17] = -input[5];
      output[18] = input[8];
      output[19] = input[9];
   } else /*if(face == 5)*/ {
      output[0] = input[0];
      output[1] = input[0];
      output[2] = input[1];
      output[3] = input[2];
      output[4] = input[0];
      output[5] = input[1];
      output[6] = input[3];
      output[7] = input[2];
      output[8] = input[4];
      output[9] = input[5];
      output[10] = input[0];
      output[11] = input[1];
      output[12] = input[3];
      output[13] = input[6];
      output[14] = input[2];
      output[15] = input[4];
      output[16] = input[7];
      output[17] = input[5];
      output[18] = input[8];
      output[19] = input[9];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<3, 5>> dgAnalyze_3D_O5(std::array<T, squareSize<3, 5>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 200> buffer;
   std::array<T, 35> output;
   static U const C0 = U(2.3692688505618908e-01);
   static U const C1 = U(4.7862867049936647e-01);
   static U const C2 = U(5.6888888888888889e-01);
   static U const C3 = U(2.1469836819894497e-01);
   static U const C4 = U(2.5772685000059420e-01);
   static U const C5 = U(3.1147336576387012e-02);
   static U const C6 = U(1.7336955879860924e-01);
   static U const C7 = U(2.8444444444444444e-01);
   static U const C8 = U(1.1870775467166646e-01);
   static U const C9 = U(1.9977104139692381e-01);
   static U const C10 = U(5.8221336871210075e-02);
   static U const C11 = U(1.6488800353787675e-01);
   static U const C12 = U(2.1333333333333335e-01);
   buffer[0] = C0 * (input[0] + input[4]) + C1 * (input[1] + input[3]) + C2 * input[2];
   buffer[1] = C0 * (input[5] + input[9]) + C1 * (input[6] + input[8]) + C2 * input[7];
   buffer[2] = C0 * (input[10] + input[14]) + C1 * (input[11] + input[13]) + C2 * input[12];
   buffer[3] = C0 * (input[15] + input[19]) + C1 * (input[16] + input[18]) + C2 * input[17];
   buffer[4] = C0 * (input[20] + input[24]) + C1 * (input[21] + input[23]) + C2 * input[22];
   buffer[5] = C0 * (input[25] + input[29]) + C1 * (input[26] + input[28]) + C2 * input[27];
   buffer[6] = C0 * (input[30] + input[34]) + C1 * (input[31] + input[33]) + C2 * input[32];
   buffer[7] = C0 * (input[35] + input[39]) + C1 * (input[36] + input[38]) + C2 * input[37];
   buffer[8] = C0 * (input[40] + input[44]) + C1 * (input[41] + input[43]) + C2 * input[42];
   buffer[9] = C0 * (input[45] + input[49]) + C1 * (input[46] + input[48]) + C2 * input[47];
   buffer[10] = C0 * (input[50] + input[54]) + C1 * (input[51] + input[53]) + C2 * input[52];
   buffer[11] = C0 * (input[55] + input[59]) + C1 * (input[56] + input[58]) + C2 * input[57];
   buffer[12] = C0 * (input[60] + input[64]) + C1 * (input[61] + input[63]) + C2 * input[62];
   buffer[13] = C0 * (input[65] + input[69]) + C1 * (input[66] + input[68]) + C2 * input[67];
   buffer[14] = C0 * (input[70] + input[74]) + C1 * (input[71] + input[73]) + C2 * input[72];
   buffer[15] = C0 * (input[75] + input[79]) + C1 * (input[76] + input[78]) + C2 * input[77];
   buffer[16] = C0 * (input[80] + input[84]) + C1 * (input[81] + input[83]) + C2 * input[82];
   buffer[17] = C0 * (input[85] + input[89]) + C1 * (input[86] + input[88]) + C2 * input[87];
   buffer[18] = C0 * (input[90] + input[94]) + C1 * (input[91] + input[93]) + C2 * input[92];
   buffer[19] = C0 * (input[95] + input[99]) + C1 * (input[96] + input[98]) + C2 * input[97];
   buffer[20] = C0 * (input[100] + input[104]) + C1 * (input[101] + input[103]) + C2 * input[102];
   buffer[21] = C0 * (input[105] + input[109]) + C1 * (input[106] + input[108]) + C2 * input[107];
   buffer[22] = C0 * (input[110] + input[114]) + C1 * (input[111] + input[113]) + C2 * input[112];
   buffer[23] = C0 * (input[115] + input[119]) + C1 * (input[116] + input[118]) + C2 * input[117];
   buffer[24] = C0 * (input[120] + input[124]) + C1 * (input[121] + input[123]) + C2 * input[122];
   buffer[25] = C3 * (-input[0] + input[4]) + C4 * (-input[1] + input[3]);
   buffer[26] = C3 * (-input[5] + input[9]) + C4 * (-input[6] + input[8]);
   buffer[27] = C3 * (-input[10] + input[14]) + C4 * (-input[11] + input[13]);
   buffer[28] = C3 * (-input[15] + input[19]) + C4 * (-input[16] + input[18]);
   buffer[29] = C3 * (-input[20] + input[24]) + C4 * (-input[21] + input[23]);
   buffer[30] = C3 * (-input[25] + input[29]) + C4 * (-input[26] + input[28]);
   buffer[31] = C3 * (-input[30] + input[34]) + C4 * (-input[31] + input[33]);
   buffer[32] = C3 * (-input[35] + input[39]) + C4 * (-input[36] + input[38]);
   buffer[33] = C3 * (-input[40] + input[44]) + C4 * (-input[41] + input[43]);
   buffer[34] = C3 * (-input[45] + input[49]) + C4 * (-input[46] + input[48]);
   buffer[35] = C3 * (-input[50] + input[54]) + C4 * (-input[51] + input[53]);
   buffer[36] = C3 * (-input[55] + input[59]) + C4 * (-input[56] + input[58]);
   buffer[37] = C3 * (-input[60] + input[64]) + C4 * (-input[61] + input[63]);
   buffer[38] = C3 * (-input[65] + input[69]) + C4 * (-input[66] + input[68]);
   buffer[39] = C3 * (-input[70] + input[74]) + C4 * (-input[71] + input[73]);
   buffer[40] = C3 * (-input[75] + input[79]) + C4 * (-input[76] + input[78]);
   buffer[41] = C3 * (-input[80] + input[84]) + C4 * (-input[81] + input[83]);
   buffer[42] = C3 * (-input[85] + input[89]) + C4 * (-input[86] + input[88]);
   buffer[43] = C3 * (-input[90] + input[94]) + C4 * (-input[91] + input[93]);
   buffer[44] = C3 * (-input[95] + input[99]) + C4 * (-input[96] + input[98]);
   buffer[45] = C3 * (-input[100] + input[104]) + C4 * (-input[101] + input[103]);
   buffer[46] = C3 * (-input[105] + input[109]) + C4 * (-input[106] + input[108]);
   buffer[47] = C3 * (-input[110] + input[114]) + C4 * (-input[111] + input[113]);
   buffer[48] = C3 * (-input[115] + input[119]) + C4 * (-input[116] + input[118]);
   buffer[49] = C3 * (-input[120] + input[124]) + C4 * (-input[121] + input[123]);
   buffer[50] = C5 * (-input[1] - input[3]) + C6 * (input[0] + input[4]) - C7 * input[2];
   buffer[51] = C5 * (-input[6] - input[8]) + C6 * (input[5] + input[9]) - C7 * input[7];
   buffer[52] = C5 * (-input[11] - input[13]) + C6 * (input[10] + input[14]) - C7 * input[12];
   buffer[53] = C5 * (-input[16] - input[18]) + C6 * (input[15] + input[19]) - C7 * input[17];
   buffer[54] = C5 * (-input[21] - input[23]) + C6 * (input[20] + input[24]) - C7 * input[22];
   buffer[55] = C5 * (-input[26] - input[28]) + C6 * (input[25] + input[29]) - C7 * input[27];
   buffer[56] = C5 * (-input[31] - input[33]) + C6 * (input[30] + input[34]) - C7 * input[32];
   buffer[57] = C5 * (-input[36] - input[38]) + C6 * (input[35] + input[39]) - C7 * input[37];
   buffer[58] = C5 * (-input[41] - input[43]) + C6 * (input[40] + input[44]) - C7 * input[42];
   buffer[59] = C5 * (-input[46] - input[48]) + C6 * (input[45] + input[49]) - C7 * input[47];
   buffer[60] = C5 * (-input[51] - input[53]) + C6 * (input[50] + input[54]) - C7 * input[52];
   buffer[61] = C5 * (-input[56] - input[58]) + C6 * (input[55] + input[59]) - C7 * input[57];
   buffer[62] = C5 * (-input[61] - input[63]) + C6 * (input[60] + input[64]) - C7 * input[62];
   buffer[63] = C5 * (-input[66] - input[68]) + C6 * (input[65] + input[69]) - C7 * input[67];
   buffer[64] = C5 * (-input[71] - input[73]) + C6 * (input[70] + input[74]) - C7 * input[72];
   buffer[65] = C5 * (-input[76] - input[78]) + C6 * (input[75] + input[79]) - C7 * input[77];
   buffer[66] = C5 * (-input[81] - input[83]) + C6 * (input[80] + input[84]) - C7 * input[82];
   buffer[67] = C5 * (-input[86] - input[88]) + C6 * (input[85] + input[89]) - C7 * input[87];
   buffer[68] = C5 * (-input[91] - input[93]) + C6 * (input[90] + input[94]) - C7 * input[92];
   buffer[69] = C5 * (-input[96] - input[98]) + C6 * (input[95] + input[99]) - C7 * input[97];
   buffer[70] = C5 * (-input[101] - input[103]) + C6 * (input[100] + input[104]) - C7 * input[102];
   buffer[71] = C5 * (-input[106] - input[108]) + C6 * (input[105] + input[109]) - C7 * input[107];
   buffer[72] = C5 * (-input[111] - input[113]) + C6 * (input[110] + input[114]) - C7 * input[112];
   buffer[73] = C5 * (-input[116] - input[118]) + C6 * (input[115] + input[119]) - C7 * input[117];
   buffer[74] = C5 * (-input[121] - input[123]) + C6 * (input[120] + input[124]) - C7 * input[122];
   buffer[75] = C8 * (-input[0] + input[4]) + C9 * (input[1] - input[3]);
   buffer[76] = C8 * (-input[5] + input[9]) + C9 * (input[6] - input[8]);
   buffer[77] = C8 * (-input[10] + input[14]) + C9 * (input[11] - input[13]);
   buffer[78] = C8 * (-input[15] + input[19]) + C9 * (input[16] - input[18]);
   buffer[79] = C8 * (-input[20] + input[24]) + C9 * (input[21] - input[23]);
   buffer[80] = C8 * (-input[25] + input[29]) + C9 * (input[26] - input[28]);
   buffer[81] = C8 * (-input[30] + input[34]) + C9 * (input[31] - input[33]);
   buffer[82] = C8 * (-input[35] + input[39]) + C9 * (input[36] - input[38]);
   buffer[83] = C8 * (-input[40] + input[44]) + C9 * (input[41] - input[43]);
   buffer[84] = C8 * (-input[45] + input[49]) + C9 * (input[46] - input[48]);
   buffer[85] = C8 * (-input[50] + input[54]) + C9 * (input[51] - input[53]);
   buffer[86] = C8 * (-input[55] + input[59]) + C9 * (input[56] - input[58]);
   buffer[87] = C8 * (-input[60] + input[64]) + C9 * (input[61] - input[63]);
   buffer[88] = C8 * (-input[65] + input[69]) + C9 * (input[66] - input[68]);
   buffer[89] = C8 * (-input[70] + input[74]) + C9 * (input[71] - input[73]);
   buffer[90] = C8 * (-input[75] + input[79]) + C9 * (input[76] - input[78]);
   buffer[91] = C8 * (-input[80] + input[84]) + C9 * (input[81] - input[83]);
   buffer[92] = C8 * (-input[85] + input[89]) + C9 * (input[86] - input[88]);
   buffer[93] = C8 * (-input[90] + input[94]) + C9 * (input[91] - input[93]);
   buffer[94] = C8 * (-input[95] + input[99]) + C9 * (input[96] - input[98]);
   buffer[95] = C8 * (-input[100] + input[104]) + C9 * (input[101] - input[103]);
   buffer[96] = C8 * (-input[105] + input[109]) + C9 * (input[106] - input[108]);
   buffer[97] = C8 * (-input[110] + input[114]) + C9 * (input[111] - input[113]);
   buffer[98] = C8 * (-input[115] + input[119]) + C9 * (input[116] - input[118]);
   buffer[99] = C8 * (-input[120] + input[124]) + C9 * (input[121] - input[123]);
   buffer[100] = C10 * (input[0] + input[4]) + C11 * (-input[1] - input[3]) + C12 * input[2];
   buffer[101] = C10 * (input[5] + input[9]) + C11 * (-input[6] - input[8]) + C12 * input[7];
   buffer[102] = C10 * (input[10] + input[14]) + C11 * (-input[11] - input[13]) + C12 * input[12];
   buffer[103] = C10 * (input[15] + input[19]) + C11 * (-input[16] - input[18]) + C12 * input[17];
   buffer[104] = C10 * (input[20] + input[24]) + C11 * (-input[21] - input[23]) + C12 * input[22];
   buffer[105] = C10 * (input[25] + input[29]) + C11 * (-input[26] - input[28]) + C12 * input[27];
   buffer[106] = C10 * (input[30] + input[34]) + C11 * (-input[31] - input[33]) + C12 * input[32];
   buffer[107] = C10 * (input[35] + input[39]) + C11 * (-input[36] - input[38]) + C12 * input[37];
   buffer[108] = C10 * (input[40] + input[44]) + C11 * (-input[41] - input[43]) + C12 * input[42];
   buffer[109] = C10 * (input[45] + input[49]) + C11 * (-input[46] - input[48]) + C12 * input[47];
   buffer[110] = C10 * (input[50] + input[54]) + C11 * (-input[51] - input[53]) + C12 * input[52];
   buffer[111] = C10 * (input[55] + input[59]) + C11 * (-input[56] - input[58]) + C12 * input[57];
   buffer[112] = C10 * (input[60] + input[64]) + C11 * (-input[61] - input[63]) + C12 * input[62];
   buffer[113] = C10 * (input[65] + input[69]) + C11 * (-input[66] - input[68]) + C12 * input[67];
   buffer[114] = C10 * (input[70] + input[74]) + C11 * (-input[71] - input[73]) + C12 * input[72];
   buffer[115] = C10 * (input[75] + input[79]) + C11 * (-input[76] - input[78]) + C12 * input[77];
   buffer[116] = C10 * (input[80] + input[84]) + C11 * (-input[81] - input[83]) + C12 * input[82];
   buffer[117] = C10 * (input[85] + input[89]) + C11 * (-input[86] - input[88]) + C12 * input[87];
   buffer[118] = C10 * (input[90] + input[94]) + C11 * (-input[91] - input[93]) + C12 * input[92];
   buffer[119] = C10 * (input[95] + input[99]) + C11 * (-input[96] - input[98]) + C12 * input[97];
   buffer[120] = C10 * (input[100] + input[104]) + C11 * (-input[101] - input[103]) + C12 * input[102];
   buffer[121] = C10 * (input[105] + input[109]) + C11 * (-input[106] - input[108]) + C12 * input[107];
   buffer[122] = C10 * (input[110] + input[114]) + C11 * (-input[111] - input[113]) + C12 * input[112];
   buffer[123] = C10 * (input[115] + input[119]) + C11 * (-input[116] - input[118]) + C12 * input[117];
   buffer[124] = C10 * (input[120] + input[124]) + C11 * (-input[121] - input[123]) + C12 * input[122];
   buffer[125] = C0 * (buffer[0] + buffer[4]) + C1 * (buffer[1] + buffer[3]) + C2 * buffer[2];
   buffer[126] = C0 * (buffer[5] + buffer[9]) + C1 * (buffer[6] + buffer[8]) + C2 * buffer[7];
   buffer[127] = C0 * (buffer[10] + buffer[14]) + C1 * (buffer[11] + buffer[13]) + C2 * buffer[12];
   buffer[128] = C0 * (buffer[15] + buffer[19]) + C1 * (buffer[16] + buffer[18]) + C2 * buffer[17];
   buffer[129] = C0 * (buffer[20] + buffer[24]) + C1 * (buffer[21] + buffer[23]) + C2 * buffer[22];
   buffer[130] = C0 * (buffer[25] + buffer[29]) + C1 * (buffer[26] + buffer[28]) + C2 * buffer[27];
   buffer[131] = C3 * (-buffer[0] + buffer[4]) + C4 * (-buffer[1] + buffer[3]);
   buffer[132] = C0 * (buffer[30] + buffer[34]) + C1 * (buffer[31] + buffer[33]) + C2 * buffer[32];
   buffer[133] = C3 * (-buffer[5] + buffer[9]) + C4 * (-buffer[6] + buffer[8]);
   buffer[134] = C0 * (buffer[35] + buffer[39]) + C1 * (buffer[36] + buffer[38]) + C2 * buffer[37];
   buffer[135] = C3 * (-buffer[10] + buffer[14]) + C4 * (-buffer[11] + buffer[13]);
   buffer[136] = C0 * (buffer[40] + buffer[44]) + C1 * (buffer[41] + buffer[43]) + C2 * buffer[42];
   buffer[137] = C3 * (-buffer[15] + buffer[19]) + C4 * (-buffer[16] + buffer[18]);
   buffer[138] = C0 * (buffer[45] + buffer[49]) + C1 * (buffer[46] + buffer[48]) + C2 * buffer[47];
   buffer[139] = C3 * (-buffer[20] + buffer[24]) + C4 * (-buffer[21] + buffer[23]);
   buffer[140] = C0 * (buffer[50] + buffer[54]) + C1 * (buffer[51] + buffer[53]) + C2 * buffer[52];
   buffer[141] = C3 * (-buffer[25] + buffer[29]) + C4 * (-buffer[26] + buffer[28]);
   buffer[142] = C5 * (-buffer[1] - buffer[3]) + C6 * (buffer[0] + buffer[4]) - C7 * buffer[2];
   buffer[143] = C0 * (buffer[55] + buffer[59]) + C1 * (buffer[56] + buffer[58]) + C2 * buffer[57];
   buffer[144] = C3 * (-buffer[30] + buffer[34]) + C4 * (-buffer[31] + buffer[33]);
   buffer[145] = C5 * (-buffer[6] - buffer[8]) + C6 * (buffer[5] + buffer[9]) - C7 * buffer[7];
   buffer[146] = C0 * (buffer[60] + buffer[64]) + C1 * (buffer[61] + buffer[63]) + C2 * buffer[62];
   buffer[147] = C3 * (-buffer[35] + buffer[39]) + C4 * (-buffer[36] + buffer[38]);
   buffer[148] = C5 * (-buffer[11] - buffer[13]) + C6 * (buffer[10] + buffer[14]) - C7 * buffer[12];
   buffer[149] = C0 * (buffer[65] + buffer[69]) + C1 * (buffer[66] + buffer[68]) + C2 * buffer[67];
   buffer[150] = C3 * (-buffer[40] + buffer[44]) + C4 * (-buffer[41] + buffer[43]);
   buffer[151] = C5 * (-buffer[16] - buffer[18]) + C6 * (buffer[15] + buffer[19]) - C7 * buffer[17];
   buffer[152] = C0 * (buffer[70] + buffer[74]) + C1 * (buffer[71] + buffer[73]) + C2 * buffer[72];
   buffer[153] = C3 * (-buffer[45] + buffer[49]) + C4 * (-buffer[46] + buffer[48]);
   buffer[154] = C5 * (-buffer[21] - buffer[23]) + C6 * (buffer[20] + buffer[24]) - C7 * buffer[22];
   buffer[155] = C0 * (buffer[75] + buffer[79]) + C1 * (buffer[76] + buffer[78]) + C2 * buffer[77];
   buffer[156] = C3 * (-buffer[50] + buffer[54]) + C4 * (-buffer[51] + buffer[53]);
   buffer[157] = C5 * (-buffer[26] - buffer[28]) + C6 * (buffer[25] + buffer[29]) - C7 * buffer[27];
   buffer[158] = C8 * (-buffer[0] + buffer[4]) + C9 * (buffer[1] - buffer[3]);
   buffer[159] = C0 * (buffer[80] + buffer[84]) + C1 * (buffer[81] + buffer[83]) + C2 * buffer[82];
   buffer[160] = C3 * (-buffer[55] + buffer[59]) + C4 * (-buffer[56] + buffer[58]);
   buffer[161] = C5 * (-buffer[31] - buffer[33]) + C6 * (buffer[30] + buffer[34]) - C7 * buffer[32];
   buffer[162] = C8 * (-buffer[5] + buffer[9]) + C9 * (buffer[6] - buffer[8]);
   buffer[163] = C0 * (buffer[85] + buffer[89]) + C1 * (buffer[86] + buffer[88]) + C2 * buffer[87];
   buffer[164] = C3 * (-buffer[60] + buffer[64]) + C4 * (-buffer[61] + buffer[63]);
   buffer[165] = C5 * (-buffer[36] - buffer[38]) + C6 * (buffer[35] + buffer[39]) - C7 * buffer[37];
   buffer[166] = C8 * (-buffer[10] + buffer[14]) + C9 * (buffer[11] - buffer[13]);
   buffer[167] = C0 * (buffer[90] + buffer[94]) + C1 * (buffer[91] + buffer[93]) + C2 * buffer[92];
   buffer[168] = C3 * (-buffer[65] + buffer[69]) + C4 * (-buffer[66] + buffer[68]);
   buffer[169] = C5 * (-buffer[41] - buffer[43]) + C6 * (buffer[40] + buffer[44]) - C7 * buffer[42];
   buffer[170] = C8 * (-buffer[15] + buffer[19]) + C9 * (buffer[16] - buffer[18]);
   buffer[171] = C0 * (buffer[95] + buffer[99]) + C1 * (buffer[96] + buffer[98]) + C2 * buffer[97];
   buffer[172] = C3 * (-buffer[70] + buffer[74]) + C4 * (-buffer[71] + buffer[73]);
   buffer[173] = C5 * (-buffer[46] - buffer[48]) + C6 * (buffer[45] + buffer[49]) - C7 * buffer[47];
   buffer[174] = C8 * (-buffer[20] + buffer[24]) + C9 * (buffer[21] - buffer[23]);
   buffer[175] = C0 * (buffer[100] + buffer[104]) + C1 * (buffer[101] + buffer[103]) + C2 * buffer[102];
   buffer[176] = C3 * (-buffer[75] + buffer[79]) + C4 * (-buffer[76] + buffer[78]);
   buffer[177] = C5 * (-buffer[51] - buffer[53]) + C6 * (buffer[50] + buffer[54]) - C7 * buffer[52];
   buffer[178] = C8 * (-buffer[25] + buffer[29]) + C9 * (buffer[26] - buffer[28]);
   buffer[179] = C10 * (buffer[0] + buffer[4]) + C11 * (-buffer[1] - buffer[3]) + C12 * buffer[2];
   buffer[180] = C0 * (buffer[105] + buffer[109]) + C1 * (buffer[106] + buffer[108]) + C2 * buffer[107];
   buffer[181] = C3 * (-buffer[80] + buffer[84]) + C4 * (-buffer[81] + buffer[83]);
   buffer[182] = C5 * (-buffer[56] - buffer[58]) + C6 * (buffer[55] + buffer[59]) - C7 * buffer[57];
   buffer[183] = C8 * (-buffer[30] + buffer[34]) + C9 * (buffer[31] - buffer[33]);
   buffer[184] = C10 * (buffer[5] + buffer[9]) + C11 * (-buffer[6] - buffer[8]) + C12 * buffer[7];
   buffer[185] = C0 * (buffer[110] + buffer[114]) + C1 * (buffer[111] + buffer[113]) + C2 * buffer[112];
   buffer[186] = C3 * (-buffer[85] + buffer[89]) + C4 * (-buffer[86] + buffer[88]);
   buffer[187] = C5 * (-buffer[61] - buffer[63]) + C6 * (buffer[60] + buffer[64]) - C7 * buffer[62];
   buffer[188] = C8 * (-buffer[35] + buffer[39]) + C9 * (buffer[36] - buffer[38]);
   buffer[189] = C10 * (buffer[10] + buffer[14]) + C11 * (-buffer[11] - buffer[13]) + C12 * buffer[12];
   buffer[190] = C0 * (buffer[115] + buffer[119]) + C1 * (buffer[116] + buffer[118]) + C2 * buffer[117];
   buffer[191] = C3 * (-buffer[90] + buffer[94]) + C4 * (-buffer[91] + buffer[93]);
   buffer[192] = C5 * (-buffer[66] - buffer[68]) + C6 * (buffer[65] + buffer[69]) - C7 * buffer[67];
   buffer[193] = C8 * (-buffer[40] + buffer[44]) + C9 * (buffer[41] - buffer[43]);
   buffer[194] = C10 * (buffer[15] + buffer[19]) + C11 * (-buffer[16] - buffer[18]) + C12 * buffer[17];
   buffer[195] = C0 * (buffer[120] + buffer[124]) + C1 * (buffer[121] + buffer[123]) + C2 * buffer[122];
   buffer[196] = C3 * (-buffer[95] + buffer[99]) + C4 * (-buffer[96] + buffer[98]);
   buffer[197] = C5 * (-buffer[71] - buffer[73]) + C6 * (buffer[70] + buffer[74]) - C7 * buffer[72];
   buffer[198] = C8 * (-buffer[45] + buffer[49]) + C9 * (buffer[46] - buffer[48]);
   buffer[199] = C10 * (buffer[20] + buffer[24]) + C11 * (-buffer[21] - buffer[23]) + C12 * buffer[22];
   output[0] = C0 * (buffer[125] + buffer[129]) + C1 * (buffer[126] + buffer[128]) + C2 * buffer[127];
   output[1] = C0 * (buffer[130] + buffer[138]) + C1 * (buffer[132] + buffer[136]) + C2 * buffer[134];
   output[2] = C0 * (buffer[131] + buffer[139]) + C1 * (buffer[133] + buffer[137]) + C2 * buffer[135];
   output[3] = C3 * (-buffer[125] + buffer[129]) + C4 * (-buffer[126] + buffer[128]);
   output[4] = C0 * (buffer[140] + buffer[152]) + C1 * (buffer[143] + buffer[149]) + C2 * buffer[146];
   output[5] = C0 * (buffer[141] + buffer[153]) + C1 * (buffer[144] + buffer[150]) + C2 * buffer[147];
   output[6] = C0 * (buffer[142] + buffer[154]) + C1 * (buffer[145] + buffer[151]) + C2 * buffer[148];
   output[7] = C3 * (-buffer[130] + buffer[138]) + C4 * (-buffer[132] + buffer[136]);
   output[8] = C3 * (-buffer[131] + buffer[139]) + C4 * (-buffer[133] + buffer[137]);
   output[9] = C5 * (-buffer[126] - buffer[128]) + C6 * (buffer[125] + buffer[129]) - C7 * buffer[127];
   output[10] = C0 * (buffer[155] + buffer[171]) + C1 * (buffer[159] + buffer[167]) + C2 * buffer[163];
   output[11] = C0 * (buffer[156] + buffer[172]) + C1 * (buffer[160] + buffer[168]) + C2 * buffer[164];
   output[12] = C0 * (buffer[157] + buffer[173]) + C1 * (buffer[161] + buffer[169]) + C2 * buffer[165];
   output[13] = C0 * (buffer[158] + buffer[174]) + C1 * (buffer[162] + buffer[170]) + C2 * buffer[166];
   output[14] = C3 * (-buffer[140] + buffer[152]) + C4 * (-buffer[143] + buffer[149]);
   output[15] = C3 * (-buffer[141] + buffer[153]) + C4 * (-buffer[144] + buffer[150]);
   output[16] = C3 * (-buffer[142] + buffer[154]) + C4 * (-buffer[145] + buffer[151]);
   output[17] = C5 * (-buffer[132] - buffer[136]) + C6 * (buffer[130] + buffer[138]) - C7 * buffer[134];
   output[18] = C5 * (-buffer[133] - buffer[137]) + C6 * (buffer[131] + buffer[139]) - C7 * buffer[135];
   output[19] = C8 * (-buffer[125] + buffer[129]) + C9 * (buffer[126] - buffer[128]);
   output[20] = C0 * (buffer[175] + buffer[195]) + C1 * (buffer[180] + buffer[190]) + C2 * buffer[185];
   output[21] = C0 * (buffer[176] + buffer[196]) + C1 * (buffer[181] + buffer[191]) + C2 * buffer[186];
   output[22] = C0 * (buffer[177] + buffer[197]) + C1 * (buffer[182] + buffer[192]) + C2 * buffer[187];
   output[23] = C0 * (buffer[178] + buffer[198]) + C1 * (buffer[183] + buffer[193]) + C2 * buffer[188];
   output[24] = C0 * (buffer[179] + buffer[199]) + C1 * (buffer[184] + buffer[194]) + C2 * buffer[189];
   output[25] = C3 * (-buffer[155] + buffer[171]) + C4 * (-buffer[159] + buffer[167]);
   output[26] = C3 * (-buffer[156] + buffer[172]) + C4 * (-buffer[160] + buffer[168]);
   output[27] = C3 * (-buffer[157] + buffer[173]) + C4 * (-buffer[161] + buffer[169]);
   output[28] = C3 * (-buffer[158] + buffer[174]) + C4 * (-buffer[162] + buffer[170]);
   output[29] = C5 * (-buffer[143] - buffer[149]) + C6 * (buffer[140] + buffer[152]) - C7 * buffer[146];
   output[30] = C5 * (-buffer[144] - buffer[150]) + C6 * (buffer[141] + buffer[153]) - C7 * buffer[147];
   output[31] = C5 * (-buffer[145] - buffer[151]) + C6 * (buffer[142] + buffer[154]) - C7 * buffer[148];
   output[32] = C8 * (-buffer[130] + buffer[138]) + C9 * (buffer[132] - buffer[136]);
   output[33] = C8 * (-buffer[131] + buffer[139]) + C9 * (buffer[133] - buffer[137]);
   output[34] = C10 * (buffer[125] + buffer[129]) + C11 * (-buffer[126] - buffer[128]) + C12 * buffer[127];
   return output;
}

template<typename T>
std::array<T, squareSize<3, 5>> dgSynthesize_3D_O5(std::array<T, triangleSize<3, 5>> const& input) {
   using U = typename ElementType<T>::type;
   std::array<T, 200> buffer;
   std::array<T, 125> output;
   static U const C0 = U(2.4573545909491201e-01);
   static U const C1 = U(5.0103117104466199e-01);
   static U const C2 = U(7.3174286977813119e-01);
   static U const C3 = U(9.0617984593866396e-01);
   static U const C4 = U(6.5076203111464545e-02);
   static U const C5 = U(3.4450089119367744e-01);
   static U const C6 = U(4.1738210372666812e-01);
   static U const C7 = U(5.3846931010568311e-01);
   static U const C8 = U(3.7500000000000000e-01);
   static U const C9 = U(5.0000000000000000e-01);
   buffer[0] = C0 * input[34] - C1 * input[19] + C2 * input[9] - C3 * input[3] + input[0];
   buffer[1] = -C4 * input[9] - C5 * input[34] + C6 * input[19] - C7 * input[3] + input[0];
   buffer[2] = C8 * input[34] - C9 * input[9] + input[0];
   buffer[3] = -C4 * input[9] - C5 * input[34] - C6 * input[19] + C7 * input[3] + input[0];
   buffer[4] = C0 * input[34] + C1 * input[19] + C2 * input[9] + C3 * input[3] + input[0];
   buffer[5] = -C1 * input[32] + C2 * input[17] - C3 * input[7] + input[1];
   buffer[6] = -C1 * input[33] + C2 * input[18] - C3 * input[8] + input[2];
   buffer[7] = -C4 * input[17] + C6 * input[32] - C7 * input[7] + input[1];
   buffer[8] = -C4 * input[18] + C6 * input[33] - C7 * input[8] + input[2];
   buffer[9] = -C9 * input[17] + input[1];
   buffer[10] = -C9 * input[18] + input[2];
   buffer[11] = -C4 * input[17] - C6 * input[32] + C7 * input[7] + input[1];
   buffer[12] = -C4 * input[18] - C6 * input[33] + C7 * input[8] + input[2];
   buffer[13] = C1 * input[32] + C2 * input[17] + C3 * input[7] + input[1];
   buffer[14] = C1 * input[33] + C2 * input[18] + C3 * input[8] + input[2];
   buffer[15] = C2 * input[29] - C3 * input[14] + input[4];
   buffer[16] = C2 * input[30] - C3 * input[15] + input[5];
   buffer[17] = C2 * input[31] - C3 * input[16] + input[6];
   buffer[18] = -C4 * input[29] - C7 * input[14] + input[4];
   buffer[19] = -C4 * input[30] - C7 * input[15] + input[5];
   buffer[20] = -C4 * input[31] - C7 * input[16] + input[6];
   buffer[21] = -C9 * input[29] + input[4];
   buffer[22] = -C9 * input[30] + input[5];
   buffer[23] = -C9 * input[31] + input[6];
   buffer[24] = -C4 * input[29] + C7 * input[14] + input[4];
   buffer[25] = -C4 * input[30] + C7 * input[15] + input[5];
   buffer[26] = -C4 * input[31] + C7 * input[16] + input[6];
   buffer[27] = C2 * input[29] + C3 * input[14] + input[4];
   buffer[28] = C2 * input[30] + C3 * input[15] + input[5];
   buffer[29] = C2 * input[31] + C3 * input[16] + input[6];
   buffer[30] = -C3 * input[25] + input[10];
   buffer[31] = -C3 * input[26] + input[11];
   buffer[32] = -C3 * input[27] + input[12];
   buffer[33] = -C3 * input[28] + input[13];
   buffer[34] = -C7 * input[25] + input[10];
   buffer[35] = -C7 * input[26] + input[11];
   buffer[36] = -C7 * input[27] + input[12];
   buffer[37] = -C7 * input[28] + input[13];
   buffer[38] = input[10];
   buffer[39] = input[11];
   buffer[40] = input[12];
   buffer[41] = input[13];
   buffer[42] = C7 * input[25] + input[10];
   buffer[43] = C7 * input[26] + input[11];
   buffer[44] = C7 * input[27] + input[12];
   buffer[45] = C7 * input[28] + input[13];
   buffer[46] = C3 * input[25] + input[10];
   buffer[47] = C3 * input[26] + input[11];
   buffer[48] = C3 * input[27] + input[12];
   buffer[49] = C3 * input[28] + input[13];
   buffer[50] = input[20];
   buffer[51] = input[21];
   buffer[52] = input[22];
   buffer[53] = input[23];
   buffer[54] = input[24];
   buffer[55] = input[20];
   buffer[56] = input[21];
   buffer[57] = input[22];
   buffer[58] = input[23];
   buffer[59] = input[24];
   buffer[60] = input[20];
   buffer[61] = input[21];
   buffer[62] = input[22];
   buffer[63] = input[23];
   buffer[64] = input[24];
   buffer[65] = input[20];
   buffer[66] = input[21];
   buffer[67] = input[22];
   buffer[68] = input[23];
   buffer[69] = input[24];
   buffer[70] = input[20];
   buffer[71] = input[21];
   buffer[72] = input[22];
   buffer[73] = input[23];
   buffer[74] = input[24];
   buffer[75] = C0 * buffer[54] - C1 * buffer[33] + C2 * buffer[17] - C3 * buffer[6] + buffer[0];
   buffer[76] = -C4 * buffer[17] - C5 * buffer[54] + C6 * buffer[33] - C7 * buffer[6] + buffer[0];
   buffer[77] = C8 * buffer[54] - C9 * buffer[17] + buffer[0];
   buffer[78] = -C4 * buffer[17] - C5 * buffer[54] - C6 * buffer[33] + C7 * buffer[6] + buffer[0];
   buffer[79] = C0 * buffer[54] + C1 * buffer[33] + C2 * buffer[17] + C3 * buffer[6] + buffer[0];
   buffer[80] = C0 * buffer[59] - C1 * buffer[37] + C2 * buffer[20] - C3 * buffer[8] + buffer[1];
   buffer[81] = -C4 * buffer[20] - C5 * buffer[59] + C6 * buffer[37] - C7 * buffer[8] + buffer[1];
   buffer[82] = C8 * buffer[59] - C9 * buffer[20] + buffer[1];
   buffer[83] = -C4 * buffer[20] - C5 * buffer[59] - C6 * buffer[37] + C7 * buffer[8] + buffer[1];
   buffer[84] = C0 * buffer[59] + C1 * buffer[37] + C2 * buffer[20] + C3 * buffer[8] + buffer[1];
   buffer[85] = C0 * buffer[64] - C1 * buffer[41] + C2 * buffer[23] - C3 * buffer[10] + buffer[2];
   buffer[86] = -C4 * buffer[23] - C5 * buffer[64] + C6 * buffer[41] - C7 * buffer[10] + buffer[2];
   buffer[87] = C8 * buffer[64] - C9 * buffer[23] + buffer[2];
   buffer[88] = -C4 * buffer[23] - C5 * buffer[64] - C6 * buffer[41] + C7 * buffer[10] + buffer[2];
   buffer[89] = C0 * buffer[64] + C1 * buffer[41] + C2 * buffer[23] + C3 * buffer[10] + buffer[2];
   buffer[90] = C0 * buffer[69] - C1 * buffer[45] + C2 * buffer[26] - C3 * buffer[12] + buffer[3];
   buffer[91] = -C4 * buffer[26] - C5 * buffer[69] + C6 * buffer[45] - C7 * buffer[12] + buffer[3];
   buffer[92] = C8 * buffer[69] - C9 * buffer[26] + buffer[3];
   buffer[93] = -C4 * buffer[26] - C5 * buffer[69] - C6 * buffer[45] + C7 * buffer[12] + buffer[3];
   buffer[94] = C0 * buffer[69] + C1 * buffer[45] + C2 * buffer[26] + C3 * buffer[12] + buffer[3];
   buffer[95] = C0 * buffer[74] - C1 * buffer[49] + C2 * buffer[29] - C3 * buffer[14] + buffer[4];
   buffer[96] = -C4 * buffer[29] - C5 * buffer[74] + C6 * buffer[49] - C7 * buffer[14] + buffer[4];
   buffer[97] = C8 * buffer[74] - C9 * buffer[29] + buffer[4];
   buffer[98] = -C4 * buffer[29] - C5 * buffer[74] - C6 * buffer[49] + C7 * buffer[14] + buffer[4];
   buffer[99] = C0 * buffer[74] + C1 * buffer[49] + C2 * buffer[29] + C3 * buffer[14] + buffer[4];
   buffer[100] = -C1 * buffer[53] + C2 * buffer[32] - C3 * buffer[16] + buffer[5];
   buffer[101] = -C4 * buffer[32] + C6 * buffer[53] - C7 * buffer[16] + buffer[5];
   buffer[102] = -C9 * buffer[32] + buffer[5];
   buffer[103] = -C4 * buffer[32] - C6 * buffer[53] + C7 * buffer[16] + buffer[5];
   buffer[104] = C1 * buffer[53] + C2 * buffer[32] + C3 * buffer[16] + buffer[5];
   buffer[105] = -C1 * buffer[58] + C2 * buffer[36] - C3 * buffer[19] + buffer[7];
   buffer[106] = -C4 * buffer[36] + C6 * buffer[58] - C7 * buffer[19] + buffer[7];
   buffer[107] = -C9 * buffer[36] + buffer[7];
   buffer[108] = -C4 * buffer[36] - C6 * buffer[58] + C7 * buffer[19] + buffer[7];
   buffer[109] = C1 * buffer[58] + C2 * buffer[36] + C3 * buffer[19] + buffer[7];
   buffer[110] = -C1 * buffer[63] + C2 * buffer[40] - C3 * buffer[22] + buffer[9];
   buffer[111] = -C4 * buffer[40] + C6 * buffer[63] - C7 * buffer[22] + buffer[9];
   buffer[112] = -C9 * buffer[40] + buffer[9];
   buffer[113] = -C4 * buffer[40] - C6 * buffer[63] + C7 * buffer[22] + buffer[9];
   buffer[114] = C1 * buffer[63] + C2 * buffer[40] + C3 * buffer[22] + buffer[9];
   buffer[115] = -C1 * buffer[68] + C2 * buffer[44] - C3 * buffer[25] + buffer[11];
   buffer[116] = -C4 * buffer[44] + C6 * buffer[68] - C7 * buffer[25] + buffer[11];
   buffer[117] = -C9 * buffer[44] + buffer[11];
   buffer[118] = -C4 * buffer[44] - C6 * buffer[68] + C7 * buffer[25] + buffer[11];
   buffer[119] = C1 * buffer[68] + C2 * buffer[44] + C3 * buffer[25] + buffer[11];
   buffer[120] = -C1 * buffer[73] + C2 * buffer[48] - C3 * buffer[28] + buffer[13];
   buffer[121] = -C4 * buffer[48] + C6 * buffer[73] - C7 * buffer[28] + buffer[13];
   buffer[122] = -C9 * buffer[48] + buffer[13];
   buffer[123] = -C4 * buffer[48] - C6 * buffer[73] + C7 * buffer[28] + buffer[13];
   buffer[124] = C1 * buffer[73] + C2 * buffer[48] + C3 * buffer[28] + buffer[13];
   buffer[125] = C2 * buffer[52] - C3 * buffer[31] + buffer[15];
   buffer[126] = -C4 * buffer[52] - C7 * buffer[31] + buffer[15];
   buffer[127] = -C9 * buffer[52] + buffer[15];
   buffer[128] = -C4 * buffer[52] + C7 * buffer[31] + buffer[15];
   buffer[129] = C2 * buffer[52] + C3 * buffer[31] + buffer[15];
   buffer[130] = C2 * buffer[57] - C3 * buffer[35] + buffer[18];
   buffer[131] = -C4 * buffer[57] - C7 * buffer[35] + buffer[18];
   buffer[132] = -C9 * buffer[57] + buffer[18];
   buffer[133] = -C4 * buffer[57] + C7 * buffer[35] + buffer[18];
   buffer[134] = C2 * buffer[57] + C3 * buffer[35] + buffer[18];
   buffer[135] = C2 * buffer[62] - C3 * buffer[39] + buffer[21];
   buffer[136] = -C4 * buffer[62] - C7 * buffer[39] + buffer[21];
   buffer[137] = -C9 * buffer[62] + buffer[21];
   buffer[138] = -C4 * buffer[62] + C7 * buffer[39] + buffer[21];
   buffer[139] = C2 * buffer[62] + C3 * buffer[39] + buffer[21];
   buffer[140] = C2 * buffer[67] - C3 * buffer[43] + buffer[24];
   buffer[141] = -C4 * buffer[67] - C7 * buffer[43] + buffer[24];
   buffer[142] = -C9 * buffer[67] + buffer[24];
   buffer[143] = -C4 * buffer[67] + C7 * buffer[43] + buffer[24];
   buffer[144] = C2 * buffer[67] + C3 * buffer[43] + buffer[24];
   buffer[145] = C2 * buffer[72] - C3 * buffer[47] + buffer[27];
   buffer[146] = -C4 * buffer[72] - C7 * buffer[47] + buffer[27];
   buffer[147] = -C9 * buffer[72] + buffer[27];
   buffer[148] = -C4 * buffer[72] + C7 * buffer[47] + buffer[27];
   buffer[149] = C2 * buffer[72] + C3 * buffer[47] + buffer[27];
   buffer[150] = -C3 * buffer[51] + buffer[30];
   buffer[151] = -C7 * buffer[51] + buffer[30];
   buffer[152] = buffer[30];
   buffer[153] = C7 * buffer[51] + buffer[30];
   buffer[154] = C3 * buffer[51] + buffer[30];
   buffer[155] = -C3 * buffer[56] + buffer[34];
   buffer[156] = -C7 * buffer[56] + buffer[34];
   buffer[157] = buffer[34];
   buffer[158] = C7 * buffer[56] + buffer[34];
   buffer[159] = C3 * buffer[56] + buffer[34];
   buffer[160] = -C3 * buffer[61] + buffer[38];
   buffer[161] = -C7 * buffer[61] + buffer[38];
   buffer[162] = buffer[38];
   buffer[163] = C7 * buffer[61] + buffer[38];
   buffer[164] = C3 * buffer[61] + buffer[38];
   buffer[165] = -C3 * buffer[66] + buffer[42];
   buffer[166] = -C7 * buffer[66] + buffer[42];
   buffer[167] = buffer[42];
   buffer[168] = C7 * buffer[66] + buffer[42];
   buffer[169] = C3 * buffer[66] + buffer[42];
   buffer[170] = -C3 * buffer[71] + buffer[46];
   buffer[171] = -C7 * buffer[71] + buffer[46];
   buffer[172] = buffer[46];
   buffer[173] = C7 * buffer[71] + buffer[46];
   buffer[174] = C3 * buffer[71] + buffer[46];
   buffer[175] = buffer[50];
   buffer[176] = buffer[50];
   buffer[177] = buffer[50];
   buffer[178] = buffer[50];
   buffer[179] = buffer[50];
   buffer[180] = buffer[55];
   buffer[181] = buffer[55];
   buffer[182] = buffer[55];
   buffer[183] = buffer[55];
   buffer[184] = buffer[55];
   buffer[185] = buffer[60];
   buffer[186] = buffer[60];
   buffer[187] = buffer[60];
   buffer[188] = buffer[60];
   buffer[189] = buffer[60];
   buffer[190] = buffer[65];
   buffer[191] = buffer[65];
   buffer[192] = buffer[65];
   buffer[193] = buffer[65];
   buffer[194] = buffer[65];
   buffer[195] = buffer[70];
   buffer[196] = buffer[70];
   buffer[197] = buffer[70];
   buffer[198] = buffer[70];
   buffer[199] = buffer[70];
   output[0] = C0 * buffer[175] - C1 * buffer[150] + C2 * buffer[125] - C3 * buffer[100] + buffer[75];
   output[1] = -C4 * buffer[125] - C5 * buffer[175] + C6 * buffer[150] - C7 * buffer[100] + buffer[75];
   output[2] = C8 * buffer[175] - C9 * buffer[125] + buffer[75];
   output[3] = -C4 * buffer[125] - C5 * buffer[175] - C6 * buffer[150] + C7 * buffer[100] + buffer[75];
   output[4] = C0 * buffer[175] + C1 * buffer[150] + C2 * buffer[125] + C3 * buffer[100] + buffer[75];
   output[5] = C0 * buffer[176] - C1 * buffer[151] + C2 * buffer[126] - C3 * buffer[101] + buffer[76];
   output[6] = -C4 * buffer[126] - C5 * buffer[176] + C6 * buffer[151] - C7 * buffer[101] + buffer[76];
   output[7] = C8 * buffer[176] - C9 * buffer[126] + buffer[76];
   output[8] = -C4 * buffer[126] - C5 * buffer[176] - C6 * buffer[151] + C7 * buffer[101] + buffer[76];
   output[9] = C0 * buffer[176] + C1 * buffer[151] + C2 * buffer[126] + C3 * buffer[101] + buffer[76];
   output[10] = C0 * buffer[177] - C1 * buffer[152] + C2 * buffer[127] - C3 * buffer[102] + buffer[77];
   output[11] = -C4 * buffer[127] - C5 * buffer[177] + C6 * buffer[152] - C7 * buffer[102] + buffer[77];
   output[12] = C8 * buffer[177] - C9 * buffer[127] + buffer[77];
   output[13] = -C4 * buffer[127] - C5 * buffer[177] - C6 * buffer[152] + C7 * buffer[102] + buffer[77];
   output[14] = C0 * buffer[177] + C1 * buffer[152] + C2 * buffer[127] + C3 * buffer[102] + buffer[77];
   output[15] = C0 * buffer[178] - C1 * buffer[153] + C2 * buffer[128] - C3 * buffer[103] + buffer[78];
   output[16] = -C4 * buffer[128] - C5 * buffer[178] + C6 * buffer[153] - C7 * buffer[103] + buffer[78];
   output[17] = C8 * buffer[178] - C9 * buffer[128] + buffer[78];
   output[18] = -C4 * buffer[128] - C5 * buffer[178] - C6 * buffer[153] + C7 * buffer[103] + buffer[78];
   output[19] = C0 * buffer[178] + C1 * buffer[153] + C2 * buffer[128] + C3 * buffer[103] + buffer[78];
   output[20] = C0 * buffer[179] - C1 * buffer[154] + C2 * buffer[129] - C3 * buffer[104] + buffer[79];
   output[21] = -C4 * buffer[129] - C5 * buffer[179] + C6 * buffer[154] - C7 * buffer[104] + buffer[79];
   output[22] = C8 * buffer[179] - C9 * buffer[129] + buffer[79];
   output[23] = -C4 * buffer[129] - C5 * buffer[179] - C6 * buffer[154] + C7 * buffer[104] + buffer[79];
   output[24] = C0 * buffer[179] + C1 * buffer[154] + C2 * buffer[129] + C3 * buffer[104] + buffer[79];
   output[25] = C0 * buffer[180] - C1 * buffer[155] + C2 * buffer[130] - C3 * buffer[105] + buffer[80];
   output[26] = -C4 * buffer[130] - C5 * buffer[180] + C6 * buffer[155] - C7 * buffer[105] + buffer[80];
   output[27] = C8 * buffer[180] - C9 * buffer[130] + buffer[80];
   output[28] = -C4 * buffer[130] - C5 * buffer[180] - C6 * buffer[155] + C7 * buffer[105] + buffer[80];
   output[29] = C0 * buffer[180] + C1 * buffer[155] + C2 * buffer[130] + C3 * buffer[105] + buffer[80];
   output[30] = C0 * buffer[181] - C1 * buffer[156] + C2 * buffer[131] - C3 * buffer[106] + buffer[81];
   output[31] = -C4 * buffer[131] - C5 * buffer[181] + C6 * buffer[156] - C7 * buffer[106] + buffer[81];
   output[32] = C8 * buffer[181] - C9 * buffer[131] + buffer[81];
   output[33] = -C4 * buffer[131] - C5 * buffer[181] - C6 * buffer[156] + C7 * buffer[106] + buffer[81];
   output[34] = C0 * buffer[181] + C1 * buffer[156] + C2 * buffer[131] + C3 * buffer[106] + buffer[81];
   output[35] = C0 * buffer[182] - C1 * buffer[157] + C2 * buffer[132] - C3 * buffer[107] + buffer[82];
   output[36] = -C4 * buffer[132] - C5 * buffer[182] + C6 * buffer[157] - C7 * buffer[107] + buffer[82];
   output[37] = C8 * buffer[182] - C9 * buffer[132] + buffer[82];
   output[38] = -C4 * buffer[132] - C5 * buffer[182] - C6 * buffer[157] + C7 * buffer[107] + buffer[82];
   output[39] = C0 * buffer[182] + C1 * buffer[157] + C2 * buffer[132] + C3 * buffer[107] + buffer[82];
   output[40] = C0 * buffer[183] - C1 * buffer[158] + C2 * buffer[133] - C3 * buffer[108] + buffer[83];
   output[41] = -C4 * buffer[133] - C5 * buffer[183] + C6 * buffer[158] - C7 * buffer[108] + buffer[83];
   output[42] = C8 * buffer[183] - C9 * buffer[133] + buffer[83];
   output[43] = -C4 * buffer[133] - C5 * buffer[183] - C6 * buffer[158] + C7 * buffer[108] + buffer[83];
   output[44] = C0 * buffer[183] + C1 * buffer[158] + C2 * buffer[133] + C3 * buffer[108] + buffer[83];
   output[45] = C0 * buffer[184] - C1 * buffer[159] + C2 * buffer[134] - C3 * buffer[109] + buffer[84];
   output[46] = -C4 * buffer[134] - C5 * buffer[184] + C6 * buffer[159] - C7 * buffer[109] + buffer[84];
   output[47] = C8 * buffer[184] - C9 * buffer[134] + buffer[84];
   output[48] = -C4 * buffer[134] - C5 * buffer[184] - C6 * buffer[159] + C7 * buffer[109] + buffer[84];
   output[49] = C0 * buffer[184] + C1 * buffer[159] + C2 * buffer[134] + C3 * buffer[109] + buffer[84];
   output[50] = C0 * buffer[185] - C1 * buffer[160] + C2 * buffer[135] - C3 * buffer[110] + buffer[85];
   output[51] = -C4 * buffer[135] - C5 * buffer[185] + C6 * buffer[160] - C7 * buffer[110] + buffer[85];
   output[52] = C8 * buffer[185] - C9 * buffer[135] + buffer[85];
   output[53] = -C4 * buffer[135] - C5 * buffer[185] - C6 * buffer[160] + C7 * buffer[110] + buffer[85];
   output[54] = C0 * buffer[185] + C1 * buffer[160] + C2 * buffer[135] + C3 * buffer[110] + buffer[85];
   output[55] = C0 * buffer[186] - C1 * buffer[161] + C2 * buffer[136] - C3 * buffer[111] + buffer[86];
   output[56] = -C4 * buffer[136] - C5 * buffer[186] + C6 * buffer[161] - C7 * buffer[111] + buffer[86];
   output[57] = C8 * buffer[186] - C9 * buffer[136] + buffer[86];
   output[58] = -C4 * buffer[136] - C5 * buffer[186] - C6 * buffer[161] + C7 * buffer[111] + buffer[86];
   output[59] = C0 * buffer[186] + C1 * buffer[161] + C2 * buffer[136] + C3 * buffer[111] + buffer[86];
   output[60] = C0 * buffer[187] - C1 * buffer[162] + C2 * buffer[137] - C3 * buffer[112] + buffer[87];
   output[61] = -C4 * buffer[137] - C5 * buffer[187] + C6 * buffer[162] - C7 * buffer[112] + buffer[87];
   output[62] = C8 * buffer[187] - C9 * buffer[137] + buffer[87];
   output[63] = -C4 * buffer[137] - C5 * buffer[187] - C6 * buffer[162] + C7 * buffer[112] + buffer[87];
   output[64] = C0 * buffer[187] + C1 * buffer[162] + C2 * buffer[137] + C3 * buffer[112] + buffer[87];
   output[65] = C0 * buffer[188] - C1 * buffer[163] + C2 * buffer[138] - C3 * buffer[113] + buffer[88];
   output[66] = -C4 * buffer[138] - C5 * buffer[188] + C6 * buffer[163] - C7 * buffer[113] + buffer[88];
   output[67] = C8 * buffer[188] - C9 * buffer[138] + buffer[88];
   output[68] = -C4 * buffer[138] - C5 * buffer[188] - C6 * buffer[163] + C7 * buffer[113] + buffer[88];
   output[69] = C0 * buffer[188] + C1 * buffer[163] + C2 * buffer[138] + C3 * buffer[113] + buffer[88];
   output[70] = C0 * buffer[189] - C1 * buffer[164] + C2 * buffer[139] - C3 * buffer[114] + buffer[89];
   output[71] = -C4 * buffer[139] - C5 * buffer[189] + C6 * buffer[164] - C7 * buffer[114] + buffer[89];
   output[72] = C8 * buffer[189] - C9 * buffer[139] + buffer[89];
   output[73] = -C4 * buffer[139] - C5 * buffer[189] - C6 * buffer[164] + C7 * buffer[114] + buffer[89];
   output[74] = C0 * buffer[189] + C1 * buffer[164] + C2 * buffer[139] + C3 * buffer[114] + buffer[89];
   output[75] = C0 * buffer[190] - C1 * buffer[165] + C2 * buffer[140] - C3 * buffer[115] + buffer[90];
   output[76] = -C4 * buffer[140] - C5 * buffer[190] + C6 * buffer[165] - C7 * buffer[115] + buffer[90];
   output[77] = C8 * buffer[190] - C9 * buffer[140] + buffer[90];
   output[78] = -C4 * buffer[140] - C5 * buffer[190] - C6 * buffer[165] + C7 * buffer[115] + buffer[90];
   output[79] = C0 * buffer[190] + C1 * buffer[165] + C2 * buffer[140] + C3 * buffer[115] + buffer[90];
   output[80] = C0 * buffer[191] - C1 * buffer[166] + C2 * buffer[141] - C3 * buffer[116] + buffer[91];
   output[81] = -C4 * buffer[141] - C5 * buffer[191] + C6 * buffer[166] - C7 * buffer[116] + buffer[91];
   output[82] = C8 * buffer[191] - C9 * buffer[141] + buffer[91];
   output[83] = -C4 * buffer[141] - C5 * buffer[191] - C6 * buffer[166] + C7 * buffer[116] + buffer[91];
   output[84] = C0 * buffer[191] + C1 * buffer[166] + C2 * buffer[141] + C3 * buffer[116] + buffer[91];
   output[85] = C0 * buffer[192] - C1 * buffer[167] + C2 * buffer[142] - C3 * buffer[117] + buffer[92];
   output[86] = -C4 * buffer[142] - C5 * buffer[192] + C6 * buffer[167] - C7 * buffer[117] + buffer[92];
   output[87] = C8 * buffer[192] - C9 * buffer[142] + buffer[92];
   output[88] = -C4 * buffer[142] - C5 * buffer[192] - C6 * buffer[167] + C7 * buffer[117] + buffer[92];
   output[89] = C0 * buffer[192] + C1 * buffer[167] + C2 * buffer[142] + C3 * buffer[117] + buffer[92];
   output[90] = C0 * buffer[193] - C1 * buffer[168] + C2 * buffer[143] - C3 * buffer[118] + buffer[93];
   output[91] = -C4 * buffer[143] - C5 * buffer[193] + C6 * buffer[168] - C7 * buffer[118] + buffer[93];
   output[92] = C8 * buffer[193] - C9 * buffer[143] + buffer[93];
   output[93] = -C4 * buffer[143] - C5 * buffer[193] - C6 * buffer[168] + C7 * buffer[118] + buffer[93];
   output[94] = C0 * buffer[193] + C1 * buffer[168] + C2 * buffer[143] + C3 * buffer[118] + buffer[93];
   output[95] = C0 * buffer[194] - C1 * buffer[169] + C2 * buffer[144] - C3 * buffer[119] + buffer[94];
   output[96] = -C4 * buffer[144] - C5 * buffer[194] + C6 * buffer[169] - C7 * buffer[119] + buffer[94];
   output[97] = C8 * buffer[194] - C9 * buffer[144] + buffer[94];
   output[98] = -C4 * buffer[144] - C5 * buffer[194] - C6 * buffer[169] + C7 * buffer[119] + buffer[94];
   output[99] = C0 * buffer[194] + C1 * buffer[169] + C2 * buffer[144] + C3 * buffer[119] + buffer[94];
   output[100] = C0 * buffer[195] - C1 * buffer[170] + C2 * buffer[145] - C3 * buffer[120] + buffer[95];
   output[101] = -C4 * buffer[145] - C5 * buffer[195] + C6 * buffer[170] - C7 * buffer[120] + buffer[95];
   output[102] = C8 * buffer[195] - C9 * buffer[145] + buffer[95];
   output[103] = -C4 * buffer[145] - C5 * buffer[195] - C6 * buffer[170] + C7 * buffer[120] + buffer[95];
   output[104] = C0 * buffer[195] + C1 * buffer[170] + C2 * buffer[145] + C3 * buffer[120] + buffer[95];
   output[105] = C0 * buffer[196] - C1 * buffer[171] + C2 * buffer[146] - C3 * buffer[121] + buffer[96];
   output[106] = -C4 * buffer[146] - C5 * buffer[196] + C6 * buffer[171] - C7 * buffer[121] + buffer[96];
   output[107] = C8 * buffer[196] - C9 * buffer[146] + buffer[96];
   output[108] = -C4 * buffer[146] - C5 * buffer[196] - C6 * buffer[171] + C7 * buffer[121] + buffer[96];
   output[109] = C0 * buffer[196] + C1 * buffer[171] + C2 * buffer[146] + C3 * buffer[121] + buffer[96];
   output[110] = C0 * buffer[197] - C1 * buffer[172] + C2 * buffer[147] - C3 * buffer[122] + buffer[97];
   output[111] = -C4 * buffer[147] - C5 * buffer[197] + C6 * buffer[172] - C7 * buffer[122] + buffer[97];
   output[112] = C8 * buffer[197] - C9 * buffer[147] + buffer[97];
   output[113] = -C4 * buffer[147] - C5 * buffer[197] - C6 * buffer[172] + C7 * buffer[122] + buffer[97];
   output[114] = C0 * buffer[197] + C1 * buffer[172] + C2 * buffer[147] + C3 * buffer[122] + buffer[97];
   output[115] = C0 * buffer[198] - C1 * buffer[173] + C2 * buffer[148] - C3 * buffer[123] + buffer[98];
   output[116] = -C4 * buffer[148] - C5 * buffer[198] + C6 * buffer[173] - C7 * buffer[123] + buffer[98];
   output[117] = C8 * buffer[198] - C9 * buffer[148] + buffer[98];
   output[118] = -C4 * buffer[148] - C5 * buffer[198] - C6 * buffer[173] + C7 * buffer[123] + buffer[98];
   output[119] = C0 * buffer[198] + C1 * buffer[173] + C2 * buffer[148] + C3 * buffer[123] + buffer[98];
   output[120] = C0 * buffer[199] - C1 * buffer[174] + C2 * buffer[149] - C3 * buffer[124] + buffer[99];
   output[121] = -C4 * buffer[149] - C5 * buffer[199] + C6 * buffer[174] - C7 * buffer[124] + buffer[99];
   output[122] = C8 * buffer[199] - C9 * buffer[149] + buffer[99];
   output[123] = -C4 * buffer[149] - C5 * buffer[199] - C6 * buffer[174] + C7 * buffer[124] + buffer[99];
   output[124] = C0 * buffer[199] + C1 * buffer[174] + C2 * buffer[149] + C3 * buffer[124] + buffer[99];
   return output;
}

template<typename T>
std::array<T, triangleSize<3, 5>> dgMassInverse_3D_O5(std::array<T, triangleSize<3, 5>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(1.2500000000000000e-01);
   static U const C1 = U(3.7500000000000000e-01);
   static U const C2 = U(6.2500000000000000e-01);
   static U const C3 = U(1.1250000000000000e+00);
   static U const C4 = U(8.7500000000000000e-01);
   static U const C5 = U(1.8750000000000000e+00);
   static U const C6 = U(3.3750000000000000e+00);
   static U const C7 = U(2.6250000000000000e+00);
   static U const C8 = U(3.1250000000000000e+00);
   static U const C9 = U(5.6250000000000000e+00);
   std::array<T, triangleSize<3, 5>> output;
   output[0] = C0 * input[0];
   output[1] = C1 * input[1];
   output[2] = C1 * input[2];
   output[3] = C1 * input[3];
   output[4] = C2 * input[4];
   output[5] = C3 * input[5];
   output[6] = C2 * input[6];
   output[7] = C3 * input[7];
   output[8] = C3 * input[8];
   output[9] = C2 * input[9];
   output[10] = C4 * input[10];
   output[11] = C5 * input[11];
   output[12] = C5 * input[12];
   output[13] = C4 * input[13];
   output[14] = C5 * input[14];
   output[15] = C6 * input[15];
   output[16] = C5 * input[16];
   output[17] = C5 * input[17];
   output[18] = C5 * input[18];
   output[19] = C4 * input[19];
   output[20] = C3 * input[20];
   output[21] = C7 * input[21];
   output[22] = C8 * input[22];
   output[23] = C7 * input[23];
   output[24] = C3 * input[24];
   output[25] = C7 * input[25];
   output[26] = C9 * input[26];
   output[27] = C9 * input[27];
   output[28] = C7 * input[28];
   output[29] = C8 * input[29];
   output[30] = C9 * input[30];
   output[31] = C8 * input[31];
   output[32] = C7 * input[32];
   output[33] = C7 * input[33];
   output[34] = C3 * input[34];
   return output;
}

template<typename T>
std::array<T, triangleSize<3, 5>> dgStiffness_3D_O5(int dimension, std::array<T, triangleSize<3, 5>> const& input) {
   using U = typename ElementType<T>::type;
   static U const C0 = U(2.0000000000000000e+00);
   std::array<T, triangleSize<3, 5>> output;
   if(dimension == 0) {
      output[0] = U(0.0);
      output[1] = U(0.0);
      output[2] = U(0.0);
      output[3] = C0 * input[0];
      output[4] = U(0.0);
      output[5] = U(0.0);
      output[6] = U(0.0);
      output[7] = C0 * input[1];
      output[8] = C0 * input[2];
      output[9] = C0 * input[3];
      output[10] = U(0.0);
      output[11] = U(0.0);
      output[12] = U(0.0);
      output[13] = U(0.0);
      output[14] = C0 * input[4];
      output[15] = C0 * input[5];
      output[16] = C0 * input[6];
      output[17] = C0 * input[7];
      output[18] = C0 * input[8];
      output[19] = C0 * (input[0] + input[9]);
      output[20] = U(0.0);
      output[21] = U(0.0);
      output[22] = U(0.0);
      output[23] = U(0.0);
      output[24] = U(0.0);
      output[25] = C0 * input[10];
      output[26] = C0 * input[11];
      output[27] = C0 * input[12];
      output[28] = C0 * input[13];
      output[29] = C0 * input[14];
      output[30] = C0 * input[15];
      output[31] = C0 * input[16];
      output[32] = C0 * (input[1] + input[17]);
      output[33] = C0 * (input[2] + input[18]);
      output[34] = C0 * (input[3] + input[19]);
   } else if(dimension == 1) {
      output[0] = U(0.0);
      output[1] = U(0.0);
      output[2] = C0 * input[0];
      output[3] = U(0.0);
      output[4] = U(0.0);
      output[5] = C0 * input[1];
      output[6] = C0 * input[2];
      output[7] = U(0.0);
      output[8] = C0 * input[3];
      output[9] = U(0.0);
      output[10] = U(0.0);
      output[11] = C0 * input[4];
      output[12] = C0 * input[5];
      output[13] = C0 * (input[0] + input[6]);
      output[14] = U(0.0);
      output[15] = C0 * input[7];
      output[16] = C0 * input[8];
      output[17] = U(0.0);
      output[18] = C0 * input[9];
      output[19] = U(0.0);
      output[20] = U(0.0);
      output[21] = C0 * input[10];
      output[22] = C0 * input[11];
      output[23] = C0 * (input[1] + input[12]);
      output[24] = C0 * (input[2] + input[13]);
      output[25] = U(0.0);
      output[26] = C0 * input[14];
      output[27] = C0 * input[15];
      output[28] = C0 * (input[3] + input[16]);
      output[29] = U(0.0);
      output[30] = C0 * input[17];
      output[31] = C0 * input[18];
      output[32] = U(0.0);
      output[33] = C0 * input[19];
      output[34] = U(0.0);
   } else /*if(dimension == 2)*/ {
      output[0] = U(0.0);
      output[1] = C0 * input[0];
      output[2] = U(0.0);
      output[3] = U(0.0);
      output[4] = C0 * input[1];
      output[5] = C0 * input[2];
      output[6] = U(0.0);
      output[7] = C0 * input[3];
      output[8] = U(0.0);
      output[9] = U(0.0);
      output[10] = C0 * (input[0] + input[4]);
      output[11] = C0 * input[5];
      output[12] = C0 * input[6];
      output[13] = U(0.0);
      output[14] = C0 * input[7];
      output[15] = C0 * input[8];
      output[16] = U(0.0);
      output[17] = C0 * input[9];
      output[18] = U(0.0);
      output[19] = U(0.0);
      output[20] = C0 * (input[1] + input[10]);
      output[21] = C0 * (input[2] + input[11]);
      output[22] = C0 * input[12];
      output[23] = C0 * input[13];
      output[24] = U(0.0);
      output[25] = C0 * (input[3] + input[14]);
      output[26] = C0 * input[15];
      output[27] = C0 * input[16];
      output[28] = U(0.0);
      output[29] = C0 * input[17];
      output[30] = C0 * input[18];
      output[31] = U(0.0);
      output[32] = C0 * input[19];
      output[33] = U(0.0);
      output[34] = U(0.0);
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<2, 5>> dgTrace_3D_O5(int face, std::array<T, triangleSize<3, 5>> const& input) {
   std::array<T, triangleSize<2, 5>> output;
   if(face == 0) {
      output[0] = input[0] - input[3] + input[9] - input[19] + input[34];
      output[1] = input[1] - input[7] + input[17] - input[32];
      output[2] = input[2] - input[8] + input[18] - input[33];
      output[3] = input[4] - input[14] + input[29];
      output[4] = input[5] - input[15] + input[30];
      output[5] = input[6] - input[16] + input[31];
      output[6] = input[10] - input[25];
      output[7] = input[11] - input[26];
      output[8] = input[12] - input[27];
      output[9] = input[13] - input[28];
      output[10] = input[20];
      output[11] = input[21];
      output[12] = input[22];
      output[13] = input[23];
      output[14] = input[24];
   } else if(face == 1) {
      output[0] = input[0] + input[3] + input[9] + input[19] + input[34];
      output[1] = input[1] + input[7] + input[17] + input[32];
      output[2] = input[2] + input[8] + input[18] + input[33];
      output[3] = input[4] + input[14] + input[29];
      output[4] = input[5] + input[15] + input[30];
      output[5] = input[6] + input[16] + input[31];
      output[6] = input[10] + input[25];
      output[7] = input[11] + input[26];
      output[8] = input[12] + input[27];
      output[9] = input[13] + input[28];
      output[10] = input[20];
      output[11] = input[21];
      output[12] = input[22];
      output[13] = input[23];
      output[14] = input[24];
   } else if(face == 2) {
      output[0] = input[0] - input[2] + input[6] - input[13] + input[24];
      output[1] = input[1] - input[5] + input[12] - input[23];
      output[2] = input[3] - input[8] + input[16] - input[28];
      output[3] = input[4] - input[11] + input[22];
      output[4] = input[7] - input[15] + input[27];
      output[5] = input[9] - input[18] + input[31];
      output[6] = input[10] - input[21];
      output[7] = input[14] - input[26];
      output[8] = input[17] - input[30];
      output[9] = input[19] - input[33];
      output[10] = input[20];
      output[11] = input[25];
      output[12] = input[29];
      output[13] = input[32];
      output[14] = input[34];
   } else if(face == 3) {
      output[0] = input[0] + input[2] + input[6] + input[13] + input[24];
      output[1] = input[1] + input[5] + input[12] + input[23];
      output[2] = input[3] + input[8] + input[16] + input[28];
      output[3] = input[4] + input[11] + input[22];
      output[4] = input[7] + input[15] + input[27];
      output[5] = input[9] + input[18] + input[31];
      output[6] = input[10] + input[21];
      output[7] = input[14] + input[26];
      output[8] = input[17] + input[30];
      output[9] = input[19] + input[33];
      output[10] = input[20];
      output[11] = input[25];
      output[12] = input[29];
      output[13] = input[32];
      output[14] = input[34];
   } else if(face == 4) {
      output[0] = input[0] - input[1] + input[4] - input[10] + input[20];
      output[1] = input[2] - input[5] + input[11] - input[21];
      output[2] = input[3] - input[7] + input[14] - input[25];
      output[3] = input[6] - input[12] + input[22];
      output[4] = input[8] - input[15] + input[26];
      output[5] = input[9] - input[17] + input[29];
      output[6] = input[13] - input[23];
      output[7] = input[16] - input[27];
      output[8] = input[18] - input[30];
      output[9] = input[19] - input[32];
      output[10] = input[24];
      output[11] = input[28];
      output[12] = input[31];
      output[13] = input[33];
      output[14] = input[34];
   } else /*if(face == 5)*/ {
      output[0] = input[0] + input[1] + input[4] + input[10] + input[20];
      output[1] = input[2] + input[5] + input[11] + input[21];
      output[2] = input[3] + input[7] + input[14] + input[25];
      output[3] = input[6] + input[12] + input[22];
      output[4] = input[8] + input[15] + input[26];
      output[5] = input[9] + input[17] + input[29];
      output[6] = input[13] + input[23];
      output[7] = input[16] + input[27];
      output[8] = input[18] + input[30];
      output[9] = input[19] + input[32];
      output[10] = input[24];
      output[11] = input[28];
      output[12] = input[31];
      output[13] = input[33];
      output[14] = input[34];
   }
   return output;
}

template<typename T>
std::array<T, triangleSize<3, 5>> dgTraceInverse_3D_O5(int face, std::array<T, triangleSize<2, 5>> const& input) {
   std::array<T, triangleSize<3, 5>> output;
   if(face == 0) {
      output[0] = input[0];
      output[1] = input[1];
      output[2] = input[2];
      output[3] = -input[0];
      output[4] = input[3];
      output[5] = input[4];
      output[6] = input[5];
      output[7] = -input[1];
      output[8] = -input[2];
      output[9] = input[0];
      output[10] = input[6];
      output[11] = input[7];
      output[12] = input[8];
      output[13] = input[9];
      output[14] = -input[3];
      output[15] = -input[4];
      output[16] = -input[5];
      output[17] = input[1];
      output[18] = input[2];
      output[19] = -input[0];
      output[20] = input[10];
      output[21] = input[11];
      output[22] = input[12];
      output[23] = input[13];
      output[24] = input[14];
      output[25] = -input[6];
      output[26] = -input[7];
      output[27] = -input[8];
      output[28] = -input[9];
      output[29] = input[3];
      output[30] = input[4];
      output[31] = input[5];
      output[32] = -input[1];
      output[33] = -input[2];
      output[34] = input[0];
   } else if(face == 1) {
      output[0] = input[0];
      output[1] = input[1];
      output[2] = input[2];
      output[3] = input[0];
      output[4] = input[3];
      output[5] = input[4];
      output[6] = input[5];
      output[7] = input[1];
      output[8] = input[2];
      output[9] = input[0];
      output[10] = input[6];
      output[11] = input[7];
      output[12] = input[8];
      output[13] = input[9];
      output[14] = input[3];
      output[15] = input[4];
      output[16] = input[5];
      output[17] = input[1];
      output[18] = input[2];
      output[19] = input[0];
      output[20] = input[10];
      output[21] = input[11];
      output[22] = input[12];
      output[23] = input[13];
      output[24] = input[14];
      output[25] = input[6];
      output[26] = input[7];
      output[27] = input[8];
      output[28] = input[9];
      output[29] = input[3];
      output[30] = input[4];
      output[31] = input[5];
      output[32] = input[1];
      output[33] = input[2];
      output[34] = input[0];
   } else if(face == 2) {
      output[0] = input[0];
      output[1] = input[1];
      output[2] = -input[0];
      output[3] = input[2];
      output[4] = input[3];
      output[5] = -input[1];
      output[6] = input[0];
      output[7] = input[4];
      output[8] = -input[2];
      output[9] = input[5];
      output[10] = input[6];
      output[11] = -input[3];
      output[12] = input[1];
      output[13] = -input[0];
      output[14] = input[7];
      output[15] = -input[4];
      output[16] = input[2];
      output[17] = input[8];
      output[18] = -input[5];
      output[19] = input[9];
      output[20] = input[10];
      output[21] = -input[6];
      output[22] = input[3];
      output[23] = -input[1];
      output[24] = input[0];
      output[25] = input[11];
      output[26] = -input[7];
      output[27] = input[4];
      output[28] = -input[2];
      output[29] = input[12];
      output[30] = -input[8];
      output[31] = input[5];
      output[32] = input[13];
      output[33] = -input[9];
      output[34] = input[14];
   } else if(face == 3) {
      output[0] = input[0];
      output[1] = input[1];
      output[2] = input[0];
      output[3] = input[2];
      output[4] = input[3];
      output[5] = input[1];
      output[6] = input[0];
      output[7] = input[4];
      output[8] = input[2];
      output[9] = input[5];
      output[10] = input[6];
      output[11] = input[3];
      output[12] = input[1];
      output[13] = input[0];
      output[14] = input[7];
      output[15] = input[4];
      output[16] = input[2];
      output[17] = input[8];
      output[18] = input[5];
      output[19] = input[9];
      output[20] = input[10];
      output[21] = input[6];
      output[22] = input[3];
      output[23] = input[1];
      output[24] = input[0];
      output[25] = input[11];
      output[26] = input[7];
      output[27] = input[4];
      output[28] = input[2];
      output[29] = input[12];
      output[30] = input[8];
      output[31] = input[5];
      output[32] = input[13];
      output[33] = input[9];
      output[34] = input[14];
   } else if(face == 4) {
      output[0] = input[0];
      output[1] = -input[0];
      output[2] = input[1];
      output[3] = input[2];
      output[4] = input[0];
      output[5] = -input[1];
      output[6] = input[3];
      output[7] = -input[2];
      output[8] = input[4];
      output[9] = input[5];
      output[10] = -input[0];
      output[11] = input[1];
      output[12] = -input[3];
      output[13] = input[6];
      output[14] = input[2];
      output[15] = -input[4];
      output[16] = input[7];
      output[17] = -input[5];
      output[18] = input[8];
      output[19] = input[9];
      output[20] = input[0];
      output[21] = -input[1];
      output[22] = input[3];
      output[23] = -input[6];
      output[24] = input[10];
      output[25] = -input[2];
      output[26] = input[4];
      output[27] = -input[7];
      output[28] = input[11];
      output[29] = input[5];
      output[30] = -input[8];
      output[31] = input[12];
      output[32] = -input[9];
      output[33] = input[13];
      output[34] = input[14];
   } else /*if(face == 5)*/ {
      output[0] = input[0];
      output[1] = input[0];
      output[2] = input[1];
      output[3] = input[2];
      output[4] = input[0];
      output[5] = input[1];
      output[6] = input[3];
      output[7] = input[2];
      output[8] = input[4];
      output[9] = input[5];
      output[10] = input[0];
      output[11] = input[1];
      output[12] = input[3];
      output[13] = input[6];
      output[14] = input[2];
      output[15] = input[4];
      output[16] = input[7];
      output[17] = input[5];
      output[18] = input[8];
      output[19] = input[9];
      output[20] = input[0];
      output[21] = input[1];
      output[22] = input[3];
      output[23] = input[6];
      output[24] = input[10];
      output[25] = input[2];
      output[26] = input[4];
      output[27] = input[7];
      output[28] = input[11];
      output[29] = input[5];
      output[30] = input[8];
      output[31] = input[12];
      output[32] = input[9];
      output[33] = input[13];
      output[34] = input[14];
   }
   return output;
}

template<typename T, int D, int O>
auto dgAnalyze(std::array<T, squareSize<D, O>> const& input) {
   if constexpr(D == 1) {
      if constexpr(O == 1) {
         return dgAnalyze_1D_O1(input);
      } else if constexpr(O == 2) {
         return dgAnalyze_1D_O2(input);
      } else if constexpr(O == 3) {
         return dgAnalyze_1D_O3(input);
      } else if constexpr(O == 4) {
         return dgAnalyze_1D_O4(input);
      }
   } else if constexpr(D == 2) {
      if constexpr(O == 1) {
         return dgAnalyze_2D_O1(input);
      } else if constexpr(O == 2) {
         return dgAnalyze_2D_O2(input);
      } else if constexpr(O == 3) {
         return dgAnalyze_2D_O3(input);
      } else if constexpr(O == 4) {
         return dgAnalyze_2D_O4(input);
      }
   } else if constexpr(D == 3) {
      if constexpr(O == 1) {
         return dgAnalyze_3D_O1(input);
      } else if constexpr(O == 2) {
         return dgAnalyze_3D_O2(input);
      } else if constexpr(O == 3) {
         return dgAnalyze_3D_O3(input);
      } else if constexpr(O == 4) {
         return dgAnalyze_3D_O4(input);
      }
   }
}

template<typename T, int D, int O>
auto dgSynthesize(std::array<T, triangleSize<D, O>> const& input) {
   if constexpr(D == 1) {
      if constexpr(O == 1) {
         return dgSynthesize_1D_O1(input);
      } else if constexpr(O == 2) {
         return dgSynthesize_1D_O2(input);
      } else if constexpr(O == 3) {
         return dgSynthesize_1D_O3(input);
      } else if constexpr(O == 4) {
         return dgSynthesize_1D_O4(input);
      }
   } else if constexpr(D == 2) {
      if constexpr(O == 1) {
         return dgSynthesize_2D_O1(input);
      } else if constexpr(O == 2) {
         return dgSynthesize_2D_O2(input);
      } else if constexpr(O == 3) {
         return dgSynthesize_2D_O3(input);
      } else if constexpr(O == 4) {
         return dgSynthesize_2D_O4(input);
      }
   } else if constexpr(D == 3) {
      if constexpr(O == 1) {
         return dgSynthesize_3D_O1(input);
      } else if constexpr(O == 2) {
         return dgSynthesize_3D_O2(input);
      } else if constexpr(O == 3) {
         return dgSynthesize_3D_O3(input);
      } else if constexpr(O == 4) {
         return dgSynthesize_3D_O4(input);
      }
   }
}

template<typename T, int D, int O>
auto dgMassInverse(std::array<T, triangleSize<D, O>> const& input) {
   if constexpr(D == 1) {
      if constexpr(O == 1) {
         return dgMassInverse_1D_O1(input);
      } else if constexpr(O == 2) {
         return dgMassInverse_1D_O2(input);
      } else if constexpr(O == 3) {
         return dgMassInverse_1D_O3(input);
      } else if constexpr(O == 4) {
         return dgMassInverse_1D_O4(input);
      }
   } else if constexpr(D == 2) {
      if constexpr(O == 1) {
         return dgMassInverse_2D_O1(input);
      } else if constexpr(O == 2) {
         return dgMassInverse_2D_O2(input);
      } else if constexpr(O == 3) {
         return dgMassInverse_2D_O3(input);
      } else if constexpr(O == 4) {
         return dgMassInverse_2D_O4(input);
      }
   } else if constexpr(D == 3) {
      if constexpr(O == 1) {
         return dgMassInverse_3D_O1(input);
      } else if constexpr(O == 2) {
         return dgMassInverse_3D_O2(input);
      } else if constexpr(O == 3) {
         return dgMassInverse_3D_O3(input);
      } else if constexpr(O == 4) {
         return dgMassInverse_3D_O4(input);
      }
   }
}

template<typename T, int D, int O>
auto dgStiffness(int dimension, std::array<T, triangleSize<D, O>> const& input) {
   if constexpr(D == 1) {
      if constexpr(O == 1) {
         return dgStiffness_1D_O1(dimension, input);
      } else if constexpr(O == 2) {
         return dgStiffness_1D_O2(dimension, input);
      } else if constexpr(O == 3) {
         return dgStiffness_1D_O3(dimension, input);
      } else if constexpr(O == 4) {
         return dgStiffness_1D_O4(dimension, input);
      }
   } else if constexpr(D == 2) {
      if constexpr(O == 1) {
         return dgStiffness_2D_O1(dimension, input);
      } else if constexpr(O == 2) {
         return dgStiffness_2D_O2(dimension, input);
      } else if constexpr(O == 3) {
         return dgStiffness_2D_O3(dimension, input);
      } else if constexpr(O == 4) {
         return dgStiffness_2D_O4(dimension, input);
      }
   } else if constexpr(D == 3) {
      if constexpr(O == 1) {
         return dgStiffness_3D_O1(dimension, input);
      } else if constexpr(O == 2) {
         return dgStiffness_3D_O2(dimension, input);
      } else if constexpr(O == 3) {
         return dgStiffness_3D_O3(dimension, input);
      } else if constexpr(O == 4) {
         return dgStiffness_3D_O4(dimension, input);
      }
   }
}

template<typename T, int D, int O>
auto dgTrace(int face, std::array<T, triangleSize<D, O>> const& input) {
   if constexpr(D == 1) {
      if constexpr(O == 1) {
         return dgTrace_1D_O1(face, input);
      } else if constexpr(O == 2) {
         return dgTrace_1D_O2(face, input);
      } else if constexpr(O == 3) {
         return dgTrace_1D_O3(face, input);
      } else if constexpr(O == 4) {
         return dgTrace_1D_O4(face, input);
      }
   } else if constexpr(D == 2) {
      if constexpr(O == 1) {
         return dgTrace_2D_O1(face, input);
      } else if constexpr(O == 2) {
         return dgTrace_2D_O2(face, input);
      } else if constexpr(O == 3) {
         return dgTrace_2D_O3(face, input);
      } else if constexpr(O == 4) {
         return dgTrace_2D_O4(face, input);
      }
   } else if constexpr(D == 3) {
      if constexpr(O == 1) {
         return dgTrace_3D_O1(face, input);
      } else if constexpr(O == 2) {
         return dgTrace_3D_O2(face, input);
      } else if constexpr(O == 3) {
         return dgTrace_3D_O3(face, input);
      } else if constexpr(O == 4) {
         return dgTrace_3D_O4(face, input);
      }
   }
}

template<typename T, int D, int O>
auto dgTraceInverse(int face, std::array<T, triangleSize<D - 1, O>> const& input) {
   if constexpr(D == 1) {
      if constexpr(O == 1) {
         return dgTraceInverse_1D_O1(face, input);
      } else if constexpr(O == 2) {
         return dgTraceInverse_1D_O2(face, input);
      } else if constexpr(O == 3) {
         return dgTraceInverse_1D_O3(face, input);
      } else if constexpr(O == 4) {
         return dgTraceInverse_1D_O4(face, input);
      }
   } else if constexpr(D == 2) {
      if constexpr(O == 1) {
         return dgTraceInverse_2D_O1(face, input);
      } else if constexpr(O == 2) {
         return dgTraceInverse_2D_O2(face, input);
      } else if constexpr(O == 3) {
         return dgTraceInverse_2D_O3(face, input);
      } else if constexpr(O == 4) {
         return dgTraceInverse_2D_O4(face, input);
      }
   } else if constexpr(D == 3) {
      if constexpr(O == 1) {
         return dgTraceInverse_3D_O1(face, input);
      } else if constexpr(O == 2) {
         return dgTraceInverse_3D_O2(face, input);
      } else if constexpr(O == 3) {
         return dgTraceInverse_3D_O3(face, input);
      } else if constexpr(O == 4) {
         return dgTraceInverse_3D_O4(face, input);
      }
   }
}

template<typename T, int D, int O>
std::array<T, D> getQuadraturePoint(int flatIndex) {
   if constexpr(D == 1) {
      if constexpr(O == 1) {
         constexpr std::array<std::array<T, D>, squareSize<D, O>> map{{
            {T( 0.00000000000000000e+00)}
         }};
         return map[flatIndex];
      } else if constexpr(O == 2) {
         constexpr std::array<std::array<T, D>, squareSize<D, O>> map{{
            {T(-5.77350269189625731e-01)}, 
            {T( 5.77350269189625731e-01)}
         }};
         return map[flatIndex];
      } else if constexpr(O == 3) {
         constexpr std::array<std::array<T, D>, squareSize<D, O>> map{{
            {T(-7.74596669241483404e-01)}, 
            {T( 0.00000000000000000e+00)}, 
            {T( 7.74596669241483404e-01)}
         }};
         return map[flatIndex];
      } else if constexpr(O == 4) {
         constexpr std::array<std::array<T, D>, squareSize<D, O>> map{{
            {T(-8.61136311594052573e-01)}, 
            {T(-3.39981043584856257e-01)}, 
            {T( 3.39981043584856257e-01)}, 
            {T( 8.61136311594052573e-01)}
         }};
         return map[flatIndex];
      }
   } else if constexpr(D == 2) {
      if constexpr(O == 1) {
         constexpr std::array<std::array<T, D>, squareSize<D, O>> map{{
            {T( 0.00000000000000000e+00), T( 0.00000000000000000e+00)}
         }};
         return map[flatIndex];
      } else if constexpr(O == 2) {
         constexpr std::array<std::array<T, D>, squareSize<D, O>> map{{
            {T(-5.77350269189625731e-01), T(-5.77350269189625731e-01)}, 
            {T(-5.77350269189625731e-01), T( 5.77350269189625731e-01)}, 
            {T( 5.77350269189625731e-01), T(-5.77350269189625731e-01)}, 
            {T( 5.77350269189625731e-01), T( 5.77350269189625731e-01)}
         }};
         return map[flatIndex];
      } else if constexpr(O == 3) {
         constexpr std::array<std::array<T, D>, squareSize<D, O>> map{{
            {T(-7.74596669241483404e-01), T(-7.74596669241483404e-01)}, 
            {T(-7.74596669241483404e-01), T( 0.00000000000000000e+00)}, 
            {T(-7.74596669241483404e-01), T( 7.74596669241483404e-01)}, 
            {T( 0.00000000000000000e+00), T(-7.74596669241483404e-01)}, 
            {T( 0.00000000000000000e+00), T( 0.00000000000000000e+00)}, 
            {T( 0.00000000000000000e+00), T( 7.74596669241483404e-01)}, 
            {T( 7.74596669241483404e-01), T(-7.74596669241483404e-01)}, 
            {T( 7.74596669241483404e-01), T( 0.00000000000000000e+00)}, 
            {T( 7.74596669241483404e-01), T( 7.74596669241483404e-01)}
         }};
         return map[flatIndex];
      } else if constexpr(O == 4) {
         constexpr std::array<std::array<T, D>, squareSize<D, O>> map{{
            {T(-8.61136311594052573e-01), T(-8.61136311594052573e-01)}, 
            {T(-8.61136311594052573e-01), T(-3.39981043584856257e-01)}, 
            {T(-8.61136311594052573e-01), T( 3.39981043584856257e-01)}, 
            {T(-8.61136311594052573e-01), T( 8.61136311594052573e-01)}, 
            {T(-3.39981043584856257e-01), T(-8.61136311594052573e-01)}, 
            {T(-3.39981043584856257e-01), T(-3.39981043584856257e-01)}, 
            {T(-3.39981043584856257e-01), T( 3.39981043584856257e-01)}, 
            {T(-3.39981043584856257e-01), T( 8.61136311594052573e-01)}, 
            {T( 3.39981043584856257e-01), T(-8.61136311594052573e-01)}, 
            {T( 3.39981043584856257e-01), T(-3.39981043584856257e-01)}, 
            {T( 3.39981043584856257e-01), T( 3.39981043584856257e-01)}, 
            {T( 3.39981043584856257e-01), T( 8.61136311594052573e-01)}, 
            {T( 8.61136311594052573e-01), T(-8.61136311594052573e-01)}, 
            {T( 8.61136311594052573e-01), T(-3.39981043584856257e-01)}, 
            {T( 8.61136311594052573e-01), T( 3.39981043584856257e-01)}, 
            {T( 8.61136311594052573e-01), T( 8.61136311594052573e-01)}
         }};
         return map[flatIndex];
      }
   } else if constexpr(D == 3) {
      if constexpr(O == 1) {
         constexpr std::array<std::array<T, D>, squareSize<D, O>> map{{
            {T( 0.00000000000000000e+00), T( 0.00000000000000000e+00), T( 0.00000000000000000e+00)}
         }};
         return map[flatIndex];
      } else if constexpr(O == 2) {
         constexpr std::array<std::array<T, D>, squareSize<D, O>> map{{
            {T(-5.77350269189625731e-01), T(-5.77350269189625731e-01), T(-5.77350269189625731e-01)}, 
            {T(-5.77350269189625731e-01), T(-5.77350269189625731e-01), T( 5.77350269189625731e-01)}, 
            {T(-5.77350269189625731e-01), T( 5.77350269189625731e-01), T(-5.77350269189625731e-01)}, 
            {T(-5.77350269189625731e-01), T( 5.77350269189625731e-01), T( 5.77350269189625731e-01)}, 
            {T( 5.77350269189625731e-01), T(-5.77350269189625731e-01), T(-5.77350269189625731e-01)}, 
            {T( 5.77350269189625731e-01), T(-5.77350269189625731e-01), T( 5.77350269189625731e-01)}, 
            {T( 5.77350269189625731e-01), T( 5.77350269189625731e-01), T(-5.77350269189625731e-01)}, 
            {T( 5.77350269189625731e-01), T( 5.77350269189625731e-01), T( 5.77350269189625731e-01)}
         }};
         return map[flatIndex];
      } else if constexpr(O == 3) {
         constexpr std::array<std::array<T, D>, squareSize<D, O>> map{{
            {T(-7.74596669241483404e-01), T(-7.74596669241483404e-01), T(-7.74596669241483404e-01)}, 
            {T(-7.74596669241483404e-01), T(-7.74596669241483404e-01), T( 0.00000000000000000e+00)}, 
            {T(-7.74596669241483404e-01), T(-7.74596669241483404e-01), T( 7.74596669241483404e-01)}, 
            {T(-7.74596669241483404e-01), T( 0.00000000000000000e+00), T(-7.74596669241483404e-01)}, 
            {T(-7.74596669241483404e-01), T( 0.00000000000000000e+00), T( 0.00000000000000000e+00)}, 
            {T(-7.74596669241483404e-01), T( 0.00000000000000000e+00), T( 7.74596669241483404e-01)}, 
            {T(-7.74596669241483404e-01), T( 7.74596669241483404e-01), T(-7.74596669241483404e-01)}, 
            {T(-7.74596669241483404e-01), T( 7.74596669241483404e-01), T( 0.00000000000000000e+00)}, 
            {T(-7.74596669241483404e-01), T( 7.74596669241483404e-01), T( 7.74596669241483404e-01)}, 
            {T( 0.00000000000000000e+00), T(-7.74596669241483404e-01), T(-7.74596669241483404e-01)}, 
            {T( 0.00000000000000000e+00), T(-7.74596669241483404e-01), T( 0.00000000000000000e+00)}, 
            {T( 0.00000000000000000e+00), T(-7.74596669241483404e-01), T( 7.74596669241483404e-01)}, 
            {T( 0.00000000000000000e+00), T( 0.00000000000000000e+00), T(-7.74596669241483404e-01)}, 
            {T( 0.00000000000000000e+00), T( 0.00000000000000000e+00), T( 0.00000000000000000e+00)}, 
            {T( 0.00000000000000000e+00), T( 0.00000000000000000e+00), T( 7.74596669241483404e-01)}, 
            {T( 0.00000000000000000e+00), T( 7.74596669241483404e-01), T(-7.74596669241483404e-01)}, 
            {T( 0.00000000000000000e+00), T( 7.74596669241483404e-01), T( 0.00000000000000000e+00)}, 
            {T( 0.00000000000000000e+00), T( 7.74596669241483404e-01), T( 7.74596669241483404e-01)}, 
            {T( 7.74596669241483404e-01), T(-7.74596669241483404e-01), T(-7.74596669241483404e-01)}, 
            {T( 7.74596669241483404e-01), T(-7.74596669241483404e-01), T( 0.00000000000000000e+00)}, 
            {T( 7.74596669241483404e-01), T(-7.74596669241483404e-01), T( 7.74596669241483404e-01)}, 
            {T( 7.74596669241483404e-01), T( 0.00000000000000000e+00), T(-7.74596669241483404e-01)}, 
            {T( 7.74596669241483404e-01), T( 0.00000000000000000e+00), T( 0.00000000000000000e+00)}, 
            {T( 7.74596669241483404e-01), T( 0.00000000000000000e+00), T( 7.74596669241483404e-01)}, 
            {T( 7.74596669241483404e-01), T( 7.74596669241483404e-01), T(-7.74596669241483404e-01)}, 
            {T( 7.74596669241483404e-01), T( 7.74596669241483404e-01), T( 0.00000000000000000e+00)}, 
            {T( 7.74596669241483404e-01), T( 7.74596669241483404e-01), T( 7.74596669241483404e-01)}
         }};
         return map[flatIndex];
      } else if constexpr(O == 4) {
         constexpr std::array<std::array<T, D>, squareSize<D, O>> map{{
            {T(-8.61136311594052573e-01), T(-8.61136311594052573e-01), T(-8.61136311594052573e-01)}, 
            {T(-8.61136311594052573e-01), T(-8.61136311594052573e-01), T(-3.39981043584856257e-01)}, 
            {T(-8.61136311594052573e-01), T(-8.61136311594052573e-01), T( 3.39981043584856257e-01)}, 
            {T(-8.61136311594052573e-01), T(-8.61136311594052573e-01), T( 8.61136311594052573e-01)}, 
            {T(-8.61136311594052573e-01), T(-3.39981043584856257e-01), T(-8.61136311594052573e-01)}, 
            {T(-8.61136311594052573e-01), T(-3.39981043584856257e-01), T(-3.39981043584856257e-01)}, 
            {T(-8.61136311594052573e-01), T(-3.39981043584856257e-01), T( 3.39981043584856257e-01)}, 
            {T(-8.61136311594052573e-01), T(-3.39981043584856257e-01), T( 8.61136311594052573e-01)}, 
            {T(-8.61136311594052573e-01), T( 3.39981043584856257e-01), T(-8.61136311594052573e-01)}, 
            {T(-8.61136311594052573e-01), T( 3.39981043584856257e-01), T(-3.39981043584856257e-01)}, 
            {T(-8.61136311594052573e-01), T( 3.39981043584856257e-01), T( 3.39981043584856257e-01)}, 
            {T(-8.61136311594052573e-01), T( 3.39981043584856257e-01), T( 8.61136311594052573e-01)}, 
            {T(-8.61136311594052573e-01), T( 8.61136311594052573e-01), T(-8.61136311594052573e-01)}, 
            {T(-8.61136311594052573e-01), T( 8.61136311594052573e-01), T(-3.39981043584856257e-01)}, 
            {T(-8.61136311594052573e-01), T( 8.61136311594052573e-01), T( 3.39981043584856257e-01)}, 
            {T(-8.61136311594052573e-01), T( 8.61136311594052573e-01), T( 8.61136311594052573e-01)}, 
            {T(-3.39981043584856257e-01), T(-8.61136311594052573e-01), T(-8.61136311594052573e-01)}, 
            {T(-3.39981043584856257e-01), T(-8.61136311594052573e-01), T(-3.39981043584856257e-01)}, 
            {T(-3.39981043584856257e-01), T(-8.61136311594052573e-01), T( 3.39981043584856257e-01)}, 
            {T(-3.39981043584856257e-01), T(-8.61136311594052573e-01), T( 8.61136311594052573e-01)}, 
            {T(-3.39981043584856257e-01), T(-3.39981043584856257e-01), T(-8.61136311594052573e-01)}, 
            {T(-3.39981043584856257e-01), T(-3.39981043584856257e-01), T(-3.39981043584856257e-01)}, 
            {T(-3.39981043584856257e-01), T(-3.39981043584856257e-01), T( 3.39981043584856257e-01)}, 
            {T(-3.39981043584856257e-01), T(-3.39981043584856257e-01), T( 8.61136311594052573e-01)}, 
            {T(-3.39981043584856257e-01), T( 3.39981043584856257e-01), T(-8.61136311594052573e-01)}, 
            {T(-3.39981043584856257e-01), T( 3.39981043584856257e-01), T(-3.39981043584856257e-01)}, 
            {T(-3.39981043584856257e-01), T( 3.39981043584856257e-01), T( 3.39981043584856257e-01)}, 
            {T(-3.39981043584856257e-01), T( 3.39981043584856257e-01), T( 8.61136311594052573e-01)}, 
            {T(-3.39981043584856257e-01), T( 8.61136311594052573e-01), T(-8.61136311594052573e-01)}, 
            {T(-3.39981043584856257e-01), T( 8.61136311594052573e-01), T(-3.39981043584856257e-01)}, 
            {T(-3.39981043584856257e-01), T( 8.61136311594052573e-01), T( 3.39981043584856257e-01)}, 
            {T(-3.39981043584856257e-01), T( 8.61136311594052573e-01), T( 8.61136311594052573e-01)}, 
            {T( 3.39981043584856257e-01), T(-8.61136311594052573e-01), T(-8.61136311594052573e-01)}, 
            {T( 3.39981043584856257e-01), T(-8.61136311594052573e-01), T(-3.39981043584856257e-01)}, 
            {T( 3.39981043584856257e-01), T(-8.61136311594052573e-01), T( 3.39981043584856257e-01)}, 
            {T( 3.39981043584856257e-01), T(-8.61136311594052573e-01), T( 8.61136311594052573e-01)}, 
            {T( 3.39981043584856257e-01), T(-3.39981043584856257e-01), T(-8.61136311594052573e-01)}, 
            {T( 3.39981043584856257e-01), T(-3.39981043584856257e-01), T(-3.39981043584856257e-01)}, 
            {T( 3.39981043584856257e-01), T(-3.39981043584856257e-01), T( 3.39981043584856257e-01)}, 
            {T( 3.39981043584856257e-01), T(-3.39981043584856257e-01), T( 8.61136311594052573e-01)}, 
            {T( 3.39981043584856257e-01), T( 3.39981043584856257e-01), T(-8.61136311594052573e-01)}, 
            {T( 3.39981043584856257e-01), T( 3.39981043584856257e-01), T(-3.39981043584856257e-01)}, 
            {T( 3.39981043584856257e-01), T( 3.39981043584856257e-01), T( 3.39981043584856257e-01)}, 
            {T( 3.39981043584856257e-01), T( 3.39981043584856257e-01), T( 8.61136311594052573e-01)}, 
            {T( 3.39981043584856257e-01), T( 8.61136311594052573e-01), T(-8.61136311594052573e-01)}, 
            {T( 3.39981043584856257e-01), T( 8.61136311594052573e-01), T(-3.39981043584856257e-01)}, 
            {T( 3.39981043584856257e-01), T( 8.61136311594052573e-01), T( 3.39981043584856257e-01)}, 
            {T( 3.39981043584856257e-01), T( 8.61136311594052573e-01), T( 8.61136311594052573e-01)}, 
            {T( 8.61136311594052573e-01), T(-8.61136311594052573e-01), T(-8.61136311594052573e-01)}, 
            {T( 8.61136311594052573e-01), T(-8.61136311594052573e-01), T(-3.39981043584856257e-01)}, 
            {T( 8.61136311594052573e-01), T(-8.61136311594052573e-01), T( 3.39981043584856257e-01)}, 
            {T( 8.61136311594052573e-01), T(-8.61136311594052573e-01), T( 8.61136311594052573e-01)}, 
            {T( 8.61136311594052573e-01), T(-3.39981043584856257e-01), T(-8.61136311594052573e-01)}, 
            {T( 8.61136311594052573e-01), T(-3.39981043584856257e-01), T(-3.39981043584856257e-01)}, 
            {T( 8.61136311594052573e-01), T(-3.39981043584856257e-01), T( 3.39981043584856257e-01)}, 
            {T( 8.61136311594052573e-01), T(-3.39981043584856257e-01), T( 8.61136311594052573e-01)}, 
            {T( 8.61136311594052573e-01), T( 3.39981043584856257e-01), T(-8.61136311594052573e-01)}, 
            {T( 8.61136311594052573e-01), T( 3.39981043584856257e-01), T(-3.39981043584856257e-01)}, 
            {T( 8.61136311594052573e-01), T( 3.39981043584856257e-01), T( 3.39981043584856257e-01)}, 
            {T( 8.61136311594052573e-01), T( 3.39981043584856257e-01), T( 8.61136311594052573e-01)}, 
            {T( 8.61136311594052573e-01), T( 8.61136311594052573e-01), T(-8.61136311594052573e-01)}, 
            {T( 8.61136311594052573e-01), T( 8.61136311594052573e-01), T(-3.39981043584856257e-01)}, 
            {T( 8.61136311594052573e-01), T( 8.61136311594052573e-01), T( 3.39981043584856257e-01)}, 
            {T( 8.61136311594052573e-01), T( 8.61136311594052573e-01), T( 8.61136311594052573e-01)}
         }};
         return map[flatIndex];
      }
   }
}

template<int D, int O>
std::array<int, D> flatToTriangular(int flatIndex) {
   if constexpr(D == 1) {
      if constexpr(O == 1) {
         constexpr std::array<std::array<int, D>, triangleSize<D, O>> map{{
            {0}
         }};
         return map[flatIndex];
      } else if constexpr(O == 2) {
         constexpr std::array<std::array<int, D>, triangleSize<D, O>> map{{
            {0}, {1}
         }};
         return map[flatIndex];
      } else if constexpr(O == 3) {
         constexpr std::array<std::array<int, D>, triangleSize<D, O>> map{{
            {0}, {1}, {2}
         }};
         return map[flatIndex];
      } else if constexpr(O == 4) {
         constexpr std::array<std::array<int, D>, triangleSize<D, O>> map{{
            {0}, {1}, {2}, {3}
         }};
         return map[flatIndex];
      }
   } else if constexpr(D == 2) {
      if constexpr(O == 1) {
         constexpr std::array<std::array<int, D>, triangleSize<D, O>> map{{
            {0, 0}
         }};
         return map[flatIndex];
      } else if constexpr(O == 2) {
         constexpr std::array<std::array<int, D>, triangleSize<D, O>> map{{
            {0, 0}, {0, 1}, {1, 0}
         }};
         return map[flatIndex];
      } else if constexpr(O == 3) {
         constexpr std::array<std::array<int, D>, triangleSize<D, O>> map{{
            {0, 0}, {0, 1}, {1, 0}, {0, 2}, 
            {1, 1}, {2, 0}
         }};
         return map[flatIndex];
      } else if constexpr(O == 4) {
         constexpr std::array<std::array<int, D>, triangleSize<D, O>> map{{
            {0, 0}, {0, 1}, {1, 0}, {0, 2}, 
            {1, 1}, {2, 0}, {0, 3}, {1, 2}, 
            {2, 1}, {3, 0}
         }};
         return map[flatIndex];
      }
   } else if constexpr(D == 3) {
      if constexpr(O == 1) {
         constexpr std::array<std::array<int, D>, triangleSize<D, O>> map{{
            {0, 0, 0}
         }};
         return map[flatIndex];
      } else if constexpr(O == 2) {
         constexpr std::array<std::array<int, D>, triangleSize<D, O>> map{{
            {0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {1, 0, 0}
         }};
         return map[flatIndex];
      } else if constexpr(O == 3) {
         constexpr std::array<std::array<int, D>, triangleSize<D, O>> map{{
            {0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {1, 0, 0}, 
            {0, 0, 2}, {0, 1, 1}, {0, 2, 0}, {1, 0, 1}, 
            {1, 1, 0}, {2, 0, 0}
         }};
         return map[flatIndex];
      } else if constexpr(O == 4) {
         constexpr std::array<std::array<int, D>, triangleSize<D, O>> map{{
            {0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {1, 0, 0}, 
            {0, 0, 2}, {0, 1, 1}, {0, 2, 0}, {1, 0, 1}, 
            {1, 1, 0}, {2, 0, 0}, {0, 0, 3}, {0, 1, 2}, 
            {0, 2, 1}, {0, 3, 0}, {1, 0, 2}, {1, 1, 1}, 
            {1, 2, 0}, {2, 0, 1}, {2, 1, 0}, {3, 0, 0}
         }};
         return map[flatIndex];
      }
   }
}

template<int D, int O>
int triangularToFlat(std::array<int, D> ti) {
	std::array<int, D> num, den;
	for (int d = D - 1; d > 0; d--) {
		ti[d - 1] += ti[d];
	}
	num = ti;
	den.fill(1);
	for (int d1 = 1; d1 < D; d1++) {
		for (int d2 = 0; d2 < D - d1; d2++) {
			num[d2] *= ti[d2] + d1;
			den[d2] *= d1 + 1;
		}
	}
	int flat = 0;
	for (int d = 0; d < D; d++) {
		flat += num[d] / den[d];
	}
	return flat;
}
