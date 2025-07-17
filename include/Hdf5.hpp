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
#include <TriIndex.hpp>

#include "Matrix.hpp"
#include "Util.hpp"
#include <valarray>

template<typename T, int D, int N, int P, template<typename> typename S>
void writeHdf5(std::string filename, T const &h, S<std::array<std::valarray<T>, binco(D + P - 1, D)>> const &fieldData,
		std::vector<std::string> const &fieldNames) {
	constexpr Range<int, D> Box { repeat<D>(0), repeat<D>(N) };
	using hindex_t = MultiIndex<Box>;
	hid_t file_id;
	file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	constexpr int M = N + 1;
	std::array<double*, D> coordArrays;
	std::array<hsize_t, D> dims;
	int N3;
	for (int d = 0; d < D; d++) {
		dims[d] = M;
	}
	N3 = 1;
	for (int d = 0; d < D; d++) {
		N3 *= dims[d];
	}
	for (int d = 0; d < D; d++) {
		coordArrays[d] = (double*) malloc(N3 * sizeof(T));
		for (auto I = hindex_t::begin(); I != hindex_t::end(); I++) {
			coordArrays[d][I] = I[d] * h;
		}
	}
	std::string coordNames[D];
	for (int d = 0; d < D; d++) {
		coordNames[d] = std::string(1, 'x' + d);
	}
	hid_t dataset_id, dataspace_id;
	auto const H5T_data_type = (sizeof(T) == sizeof(double)) ? H5T_NATIVE_DOUBLE : H5T_NATIVE_FLOAT;
	for (int d = 0; d < D; ++d) {
		dataspace_id = H5Screate_simple(D, dims.data(), NULL);
		dataset_id = H5Dcreate(file_id, coordNames[d].c_str(), H5T_data_type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dwrite(dataset_id, H5T_data_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, coordArrays[d]);
		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
	}
	for (int d = 0; d < D; d++) {
		dims[d]--;
	}
	N3 = 1;
	for (int d = 0; d < D; d++) {
		N3 *= dims[d];
	}
	using index_type = TriIndex<P, D>;
	int const nf = fieldNames.size();
	std::vector<T> buffer(N3);
	for (int fi = 0; fi < nf; fi++) {
		for (index_type Q = index_type::begin(); Q != index_type::end(); Q++) {
			std::string const name = fieldNames[fi] + "_" + std::to_string(Q);
			dataspace_id = H5Screate_simple(D, dims.data(), NULL);
			dataset_id = H5Dcreate(file_id, name.c_str(), H5T_data_type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			auto const &field = fieldData[fi][Q];
			reverseMultiIndexData<Box>(std::begin(field), std::end(field), buffer.data());
			H5Dwrite(dataset_id, H5T_data_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data());
			H5Dclose(dataset_id);
			H5Sclose(dataspace_id);
		}

	}
	for (int d = 0; d < D; d++) {
		free(coordArrays[d]);
	}
	H5Fclose(file_id);
	std::string xFilename = filename + ".xmf";
	std::ofstream xmfFile { xFilename };
	if (!xmfFile) {
		throw std::runtime_error("Unable to open XDMF file \"" + xFilename + "\" for writing");
	}
	std::string const precision = (sizeof(T) == sizeof(double)) ? "8" : "4";
	std::string const dataType = (sizeof(T) == sizeof(double)) ? "double" : "float";
	std::string const topology = std::to_string(D) + "DCoRectMesh";
	std::string cDimensions, nDimensions, geometry = "ORIGIN_";
	for (int d = 0; d < D; d++) {
		cDimensions += std::to_string(N);
		cDimensions += " ";
		nDimensions += std::to_string(N + 1);
		nDimensions += " ";
		geometry += "D";
		geometry += std::string(1, 'X' + d);
	}
	cDimensions.pop_back();
	nDimensions.pop_back();
	xmfFile << "<?xml version=\"1.0\" ?>" << std::endl;
	xmfFile << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
	xmfFile << "<Xdmf Version=\"2.0\">" << std::endl;
	xmfFile << " <Domain>" << std::endl;
	xmfFile << "  <Grid Name=\"mesh\" GridType=\"Uniform\">" << std::endl;
	xmfFile << "   <Topology TopologyType=\"" << topology << "\" Dimensions=\"" << nDimensions << "\"/>" << std::endl;
	xmfFile << "    <Geometry GeometryType=\"" << geometry << "\">" << std::endl;
	std::string minStr;
	std::string maxStr;
	for (int d = 0; d < D; d++) {
		minStr += std::to_string(0) + " ";
		maxStr += std::to_string(h) + " ";
	}
	minStr.pop_back();
	maxStr.pop_back();
	xmfFile << "     <DataItem Dimensions=\"" + std::to_string(D) + "\" Format=\"XML\">" + minStr + "</DataItem>" << std::endl;
	xmfFile << "     <DataItem Dimensions=\"" + std::to_string(D) + "\" Format=\"XML\">" + maxStr + "</DataItem>" << std::endl;
	xmfFile << "    </Geometry>" << std::endl;
	for (int fi = 0; fi < nf; fi++) {
		for (index_type Q = index_type::begin(); Q != index_type::end(); Q++) {
			std::string const name = fieldNames[fi] + "_" + std::to_string(Q);
			xmfFile << "   <Attribute Name=\"" << name << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
			xmfFile << "    <DataItem Dimensions=\"" << cDimensions << "\" NumberType=\"" << dataType << "\" Precision=\"" << precision << "\" Format=\"HDF\">"
					<< std::endl;
			xmfFile << "     " << filename << ":/" << name << std::endl;
			xmfFile << "    </DataItem>" << std::endl;
			xmfFile << "   </Attribute>" << std::endl;
		}

	}
	xmfFile << "  </Grid>" << std::endl;
	xmfFile << " </Domain>" << std::endl;
	xmfFile << "</Xdmf>" << std::endl;
}
//
#endif // WRITE_DG_HDF5_HPP
