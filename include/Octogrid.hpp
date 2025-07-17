#include "HyperSubgrid.hpp"

#include <hpx/include/components.hpp>

#include <bit>
#include <memory>

template<typename T>
int mostSignificantBit(T value) {
	constexpr T maxBit = T(CHAR_BIT * sizeof(T)) - T(1);
	return T(maxBit - std::countl_zero(value));
}

template<typename >
struct OctogridClient;

template<int D>
struct AmrLocationId {
	using ChildType = Child<D>;
	using Integer = uint64_t;
	using FaceType = Face<D>;
	static constexpr Integer zero = Integer(0);
	static constexpr Integer one = Integer(1);
	static constexpr Integer dimCount = Integer(D);
	AmrLocationId() :
			bits_(1) {
	}
	AmrLocationId getChild(ChildType child) const {
		AmrLocationId cid = *this;
		cid.bits_ <<= Integer(D);
		cid.bits_ |= Integer(child);
		return cid;
	}
	Integer level() const {
		int const msb = mostSignificantBit(bits_);
		return (msb - one) / dimCount;
	}
	Integer operator[](Integer dim) const {
		Integer bits = bits_ >> dim;
		Integer loc = zero;
		while (bits) {
			loc |= bits & one;
			bits >>= dimCount;
		}
		return loc;
	}
	bool isOuterBound(FaceType face) const {
		Integer const loc = operator[](face.dimension());
		Integer const lev = level();
		if (lev && (face.direction() > 0)) {
			Integer const maxValue = (one << lev) - one;
			return loc == maxValue;
		} else {
			return loc == zero;
		}
	}
private:
	Integer bits_;
};

enum class BoundaryType : int {
	domainDecomposition, adaptiveMeshRefinement, physical
};

template<typename Subgrid>
struct OctogridServer: public hpx::components::component_base<OctogridServer<Subgrid>> {
	template<typename T>
	using State = typename Subgrid::State<T>;
	using Type = typename Subgrid::Type;
	using RungeKutta = typename Subgrid::RungeKutta;
	static constexpr int dimensionCount = Subgrid::dimensionCount;
	static constexpr int interiorWidth = Subgrid::interiorWidth;
	static constexpr int modeCount = Subgrid::modeCount;
	static constexpr int siblingCount = 2 * dimensionCount;
	static constexpr int childCount = 1 << dimensionCount;
	static constexpr int childrenPerFace = 1 << (dimensionCount - 1);
	using FaceType = Face<Subgrid::dimensionCount>;
	using ChildType = Child<Subgrid::dimensionCount>;
	using FaceChildType = Child<dimensionCount - 1>;
	using SiblingVector = std::array<OctogridClient<Subgrid>, siblingCount>;
	using ChildrenVector = std::array<OctogridClient<Subgrid>, childCount>;
	using FaceChildrenVector = std::array<OctogridClient<Subgrid>, childrenPerFace>;
	using OctogridType = OctogridClient<Subgrid>;
	using LocationType = AmrLocationId<dimensionCount>;

	static OctogridType createRoot();

	OctogridServer();
	FaceChildrenVector getFaceChildren(Face<dimensionCount> face) const; //
	void refine(); //
	void set(LocationType, SiblingVector&&); //
	HPX_DEFINE_COMPONENT_ACTION(OctogridServer, getFaceChildren, getFaceChildrenAction);//
	HPX_DEFINE_COMPONENT_ACTION(OctogridServer, refine, refineAction);//
	HPX_DEFINE_COMPONENT_ACTION(OctogridServer, set, setAction);//

private:
	std::array<BoundaryType, siblingCount> boundaryType;
	ChildrenVector children_;
	bool isLeaf_;
	LocationType location_;
	SiblingVector siblings_;
	std::shared_ptr<Subgrid> subgrid_;
};

template<typename Subgrid>
struct OctogridClient: public hpx::components::client_base<OctogridClient<Subgrid>, OctogridServer<Subgrid>> {
	using base_type = hpx::components::client_base<OctogridClient<Subgrid>, OctogridServer<Subgrid>>;
	using server_type = OctogridServer<Subgrid>;
	static constexpr int dimensionCount = Subgrid::dimensionCount;
	hpx::future<typename server_type::FaceChildrenVector> getFaceChildren(Face<dimensionCount> face) const {
		return hpx::async<typename server_type::getFaceChildrenAction>(this->get_id(), face);
	}
	hpx::future<void> refine() const {
		return hpx::async<typename server_type::refineAction>(this->get_id());
	}
	hpx::future<void> set(typename server_type::LocationType location, typename server_type::SiblingVector &&sibs) const {
		return hpx::async<typename server_type::setAction>(this->get_id(), location, std::move(sibs));
	}
	bool operator==(hpx::id_type const &other) const {
		return this->get_id() == other;
	}
	bool operator!=(hpx::id_type const &other) const {
		return this->get_id() != other;
	}
	OctogridClient() = default;
	OctogridClient(hpx::future<hpx::id_type> &&gid) :
			base_type(std::move(gid)) {
	}
	explicit OctogridClient(hpx::id_type &&id) :
			base_type(std::move(id)) {
	}
};

template<typename Subgrid>
OctogridServer<Subgrid>::OctogridServer() {
	boundaryType.fill(BoundaryType::adaptiveMeshRefinement);
	isLeaf_ = true;
}

template<typename Subgrid>
typename OctogridServer<Subgrid>::OctogridType OctogridServer<Subgrid>::createRoot() {
	OctogridType root = hpx::new_ < OctogridServer > (hpx::find_here());
	SiblingVector siblings;
	siblings.fill(root);
	root.set(LocationType(), std::move(siblings)).get();
	return root;
}

template<typename Subgrid>
typename OctogridServer<Subgrid>::FaceChildrenVector OctogridServer<Subgrid>::getFaceChildren(Face<OctogridServer<Subgrid>::dimensionCount> face) const {
	FaceChildrenVector faceChildren;
	for (FaceChildType to = FaceChildType::begin(); to < FaceChildType::end(); to++) {
		Child<dimensionCount> const from(to, face);
		faceChildren[to] = children_[from];
	}
	return faceChildren;
}

template<typename Subgrid>
void OctogridServer<Subgrid>::refine() {
	std::array<std::array<OctogridType, siblingCount>, childCount> siblingsOfChildren;
	std::vector<hpx::future<FaceChildrenVector>> siblingFutures;
	std::vector<hpx::future<void>> childFutures;
	childFutures.reserve(childCount);
	siblingFutures.reserve(siblingCount);
	for (FaceType face = FaceType::begin(); face != FaceType::end(); face++) {
		if (boundaryType[face] == BoundaryType::domainDecomposition) {
			siblingFutures.push_back(siblings_[face].getFaceChildren(face));
		}
	}
	for (ChildType child = ChildType::begin(); child != ChildType::end(); child++) {
		children_[child] = hpx::new_ < OctogridServer > (hpx::find_here());
	}
	for (ChildType child = ChildType::begin(); child != ChildType::end(); child++) {
		for (int dim = 0; dim < dimensionCount; dim++) {
			siblingsOfChildren[child][dim] = children_[child.flip(dim)];
		}
	}
	for (FaceType face = FaceType::count() - 1; face >= 0; face--) {
		if (boundaryType[face] == BoundaryType::domainDecomposition) {
			auto const theseSiblings = siblingFutures.back().get();
			siblingFutures.pop_back();
			for (FaceChildType from = FaceChildType::begin(); from < FaceChildType::end(); from++) {
				Child<dimensionCount> const to(from, face);
				siblingsOfChildren[to][face] = theseSiblings[from];
			}
		}
	}
	for (ChildType childId = ChildType::begin(); childId != ChildType::end(); childId++) {
		childFutures.push_back(children_[childId].set(location_.getChild(childId), std::move(siblingsOfChildren[childId])));
	}
	hpx::wait_all(childFutures.begin(), childFutures.end());
	isLeaf_ = false;
}

template<typename Subgrid>
void OctogridServer<Subgrid>::set(typename OctogridServer<Subgrid>::LocationType location, typename OctogridServer<Subgrid>::SiblingVector &&sibs) {
	static hpx::id_type const nullId { };
	location_ = location;
	siblings_ = std::move(sibs);
	for (FaceType face = FaceType::begin(); face != FaceType::end(); face++) {
		if (siblings_[face] != nullId) {
			boundaryType[face] = BoundaryType::domainDecomposition;
		}
	}
}

