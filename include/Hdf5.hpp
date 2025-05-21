// writeDGHDF5.hpp

#ifndef WRITE_DG_HDF5_HPP
#define WRITE_DG_HDF5_HPP

#include <string>
#include <vector>
#include <array>
#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <type_traits>

#include <H5Cpp.h>
#include <MultiIndex.hpp>

#include "Matrix.hpp"
#include "Vector.hpp"

namespace detail {

// Compute Base^Exp at compile time
template<int Base, int Exp>
struct Pow {
	static constexpr std::size_t value = Base * Pow<Base, Exp - 1>::value;
};

template<int Base>
struct Pow<Base, 0> {
	static constexpr std::size_t value = 1;
};

} // namespace detail

//------------------------------------------------------------------------------
// Write XDMF 1.0 .xmf for a D-dimensional rectilinear mesh
//------------------------------------------------------------------------------
template<typename T, int D, int N, int P3, int BW>
void writeXdmfRectilinear1_0(std::string const &xmfFilename, std::string const &h5Filename, std::vector<std::string> const &fieldNames, bool includeGhost) {
//    static_assert(D >= 1, "D must be >=1");
//    static_assert(N >= 1, "N must be >=1");
//    static_assert(P3 >= 1, "P3 must be >=1");
//    static_assert(BW >= 0, "BW must be >=0");
//
//    int cellSize = includeGhost ? (N + 2 * BW) : N;
//    std::array<int, D> cellDims, nodeDims;
//    for (int d = 0; d < D; ++d) {
//        cellDims[d] = cellSize;
//        nodeDims[d] = cellSize + 1;
//    }
//
//    int precision = std::is_same_v<T, double> ? 8 : 4;
//
//    std::ofstream out(xmfFilename);
//    if (!out.good()) {
//        throw std::runtime_error("Failed to open " + xmfFilename);
//    }
//
//    // XDMF 1.0 header with properly quoted entity
//    out << "<?xml version=\"1.0\"?>\n"
//        << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" [\n"
//        << "  <!ENTITY HeavyData \"HDF\">\n"
//        << "]>\n"
//        << "<Xdmf Version=\"1.0\">\n"
//        << "  <Domain>\n"
//        << "    <Grid Name=\"DGMesh\" GridType=\"Uniform\">\n";
//
//    // Topology
//    out << "      <Topology TopologyType=\"" << D << "DRectMesh\" Dimensions=\"";
//    for (int d = D - 1; d >= 0; --d) {
//        out << nodeDims[d] << (d > 0 ? " " : "");
//    }
//    out << "\"/>\n";
//
//    // Geometry
//    out << "      <Geometry GeometryType=\"XYZ\">\n";
//    for (int d = 0; d < D; ++d) {
//        char axis = char('x' + d);
//        out << "        <DataItem Dimensions=\""
//            << nodeDims[d]
//            << "\" NumberType=\"Float\" Precision=\"" << precision
//            << "\" Format=\"&HeavyData;\">"
//            << h5Filename << ":/" << axis
//            << "</DataItem>\n";
//    }
//    out << "      </Geometry>\n";
//
//    // Attributes
//    for (auto const &name : fieldNames) {
//        out << "      <Attribute Name=\"" << name
//            << "\" AttributeType=\"Scalar\" Center=\"Cell\">\n"
//            << "        <DataItem Dimensions=\"";
//        for (int d = D - 1; d >= 0; --d) {
//            out << cellDims[d] << (d > 0 ? " " : "");
//        }
//        out << "\" NumberType=\"Float\" Precision=\"" << precision
//            << "\" Format=\"&HeavyData;\">"
//            << h5Filename << ":/" << name
//            << "</DataItem>\n"
//            << "      </Attribute>\n";
//    }
//
//    // Footer
//    out << "    </Grid>\n"
//        << "  </Domain>\n"
//        << "</Xdmf>\n";
//
//    out.close();
}

//------------------------------------------------------------------------------
// Write HDF5 (.h5) and accompanying XDMF (.xmf) for DG data
//------------------------------------------------------------------------------
//template<typename T, int D, int N, int P, int BW>
//void writeHdf5(std::string filename, T const &h,
//		std::vector<std::array<std::array<T, Math::integerPower(N + 2 * BW, D)>, TriangularIndex<P, D>::elementCount()>> const &fieldData,
//		std::vector<std::string> const &fieldNames) {
//	using hindex_t = IndexTuple<D, repeat<D>(-BW), repeat<D>(N + BW)>;
//	hid_t file_id;
//	file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
//	constexpr int M = N + 2 * BW + 1;
//	std::array<T*, D> coordArrays;
//	std::array<hsize_t, D> dims;
//	int N3;
//	for (int d = 0; d < D; d++) {
//		dims[d] = M;
//	}
//	N3 = 1;
//	for (int d = 0; d < D; d++) {
//		N3 *= dims[d];
//	}
//	for (int d = 0; d < D; d++) {
//		coordArrays[d] = (T*) malloc(N3 * sizeof(T));
//		for (auto I = hindex_t::begin(); I != hindex_t::end(); I++) {
//			coordArrays[d][I.flatIndex()] = (I[d] - BW) * h;
//		}
//	}
//	std::string coordNames[NDIM];
//	for (int d = 0; d < D; d++) {
//		coordNames[d] = std::string(1, 'x' + d);
//	}
//	hid_t dataset_id, dataspace_id;
//	auto const H5T_data_type = std::is_same<T, double>::value ? H5T_NATIVE_DOUBLE : H5T_NATIVE_FLOAT;
//	for (int d = 0; d < D; ++d) {
//		dataspace_id = H5Screate_simple(D, dims.data(), NULL);
//		dataset_id = H5Dcreate(file_id, coordNames[d].c_str(), H5T_data_type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//		H5Dwrite(dataset_id, H5T_data_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, coordArrays[d]);
//		H5Dclose(dataset_id);
//		H5Sclose(dataspace_id);
//	}
//	for (int d = 0; d < D; d++) {
//		dims[d]--;
//	}
//	N3 = 1;
//	for (int d = 0; d < D; d++) {
//		N3 *= dims[d];
//	}
//	using index_type = TriangularIndex<P, D>;
//	int const nf = fieldNames.size();
//	for (int fi = 0; fi < nf; fi++) {
//		for (index_type Q = index_type::begin(); Q != index_type::end(); Q++) {
//			std::string const name = fieldNames[fi] + "_" + std::to_string(Q.flatIndex());
//			dataspace_id = H5Screate_simple(D, dims.data(), NULL);
//			dataset_id = H5Dcreate(file_id, name.c_str(), H5T_data_type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//			H5Dwrite(dataset_id, H5T_data_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, fieldData[fi][Q.flatIndex()].data());
//			H5Dclose(dataset_id);
//			H5Sclose(dataspace_id);
//		}
//
//	}
//	for (int d = 0; d < D; d++) {
//		free(coordArrays[d]);
//	}
//	H5Fclose(file_id);
//	std::string xFilename = filename + ".xmf";
//	std::ofstream xmfFile { xFilename };
//	if (!xmfFile) {
//		throw std::runtime_error("Unable to open XDMF file \"" + xFilename + "\" for writing");
//	}
//	std::string const precision = std::is_same<T, double>::value ? "8" : "4";
//	std::string const dataType = std::is_same<T, double>::value ? "double" : "float";
//	std::string const topology = std::to_string(D) + "DCoRectMesh";
//	std::string cDimensions, nDimensions, geometry = "ORIGIN_";
//	for (int d = 0; d < D; d++) {
//		cDimensions += std::to_string(N + 2 * BW);
//		cDimensions += " ";
//		nDimensions += std::to_string(N + 2 * BW + 1);
//		nDimensions += " ";
//		geometry += "D";
//		geometry += std::string(1, 'X' + d);
//	}
//	cDimensions.pop_back();
//	nDimensions.pop_back();
//	xmfFile << "<?xml version=\"1.0\" ?>" << std::endl;
//	xmfFile << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
//	xmfFile << "<Xdmf Version=\"2.0\">" << std::endl;
//	xmfFile << " <Domain>" << std::endl;
//	xmfFile << "  <Grid Name=\"mesh\" GridType=\"Uniform\">" << std::endl;
//	xmfFile << "   <Topology TopologyType=\"" << topology << "\" Dimensions=\"" << nDimensions << "\"/>" << std::endl;
//	xmfFile << "    <Geometry GeometryType=\"" << geometry << "\">" << std::endl;
//	//(I[d] - BW) * h;
//	std::string minStr;
//	std::string maxStr;
//	for (int d = 0; d < D; d++) {
//		minStr += std::to_string(-BW * h) + " ";
//		maxStr += std::to_string(h) + " ";
////		minStr += std::to_string(-BW * h) + " ";
////		maxStr += std::to_string((N + BW) * h) + " ";
//	}
//	minStr.pop_back();
//	maxStr.pop_back();
//	xmfFile << "     <DataItem Dimensions=\"" + std::to_string(D) + "\" Format=\"XML\">" + minStr + "</DataItem>" << std::endl;
//	xmfFile << "     <DataItem Dimensions=\"" + std::to_string(D) + "\" Format=\"XML\">" + maxStr + "</DataItem>" << std::endl;
////	for (int d = 0; d < D; d++) {
////		xmfFile << "     <DataItem Dimensions=\"" << nDimensions << "\" NumberType=\"" << dataType << "\" Precision=\"" << precision << "\" Format=\"HDF\">"
////				<< std::endl;
////		xmfFile << "      " << filename << ":/" << std::string(1, 'X' + d) << std::endl;
////		xmfFile << "     </DataItem>" << std::endl;
////	}
//	xmfFile << "    </Geometry>" << std::endl;
//	for (int fi = 0; fi < nf; fi++) {
//		for (index_type Q = index_type::begin(); Q != index_type::end(); Q++) {
//			std::string const name = fieldNames[fi] + "_" + std::to_string(Q.flatIndex());
//			xmfFile << "   <Attribute Name=\"" << name << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
//			xmfFile << "    <DataItem Dimensions=\"" << cDimensions << "\" NumberType=\"" << dataType << "\" Precision=\"" << precision << "\" Format=\"HDF\">"
//					<< std::endl;
//			xmfFile << "     " << filename << ":/" << name << std::endl;
//			xmfFile << "    </DataItem>" << std::endl;
//			xmfFile << "   </Attribute>" << std::endl;
//		}
//
//	}
//	xmfFile << "  </Grid>" << std::endl;
//	xmfFile << " </Domain>" << std::endl;
//	xmfFile << "</Xdmf>" << std::endl;
//}
//
#endif // WRITE_DG_HDF5_HPP
