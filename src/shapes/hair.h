#if defined(_MSC_VER)
#pragma once
#endif

#if !defined(PBRT_SHAPES_HAIR_H)
#define PBRT_SHAPES_HAIR_H

// shapes/hair.h*
#include "shape.h"
#include "shapes/trianglemesh.h"
#include "accelerators/kdtreeaccel.h"

class HairKDTree : public KdTreeAccel {
public:
	HairKDTree(std::vector<Point> &vertices,
		std::vector<bool> &vertexStartsFiber, float radius);

private:
	size_t m_indexCount;
};

/**
* \brief Intersection shape structure for cylindrical hair
* segments with miter joints. This class expects an ASCII file containing
* a list of hairs made from segments. Each line should contain an X,
* Y and Z coordinate separated by a space. An empty line indicates
* the start of a new hair.
*/
class HairShape : public Shape {
public:
	// Construct a new HairShape instance
	HairShape(const Transform *o2w, const Transform *w2o, bool ro, 
		const string filename, float radius, float angleThreshold);
	~HairShape();
	BBox ObjectBound() const;
	BBox WorldBound() const;
	void Refine(vector<Reference<Shape> > &refined) const;

	/// Return the list of vertices underlying the hair shape
	const std::vector<Point> &getVertices() const;

	/**
	* Return a boolean list specifying whether a vertex
	* marks the beginning of a new fiber
	*/
	const std::vector<bool> &getStartFiber() const;

	bool Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
		DifferentialGeometry *dg) const;

	bool CanIntersect() const { return false; }
	void Refine(vector<Reference<Shape> > &refined) const;
	bool IntersectP(const Ray &ray) const;
	BBox ObjectBound() const;
	float Area() const;
	Point Sample(float u1, float u2, Normal *Ns) const;

protected:
	vector<Point> &vertices;
	vector<bool> &vertexStartsFiber;
	float radius;
};

HairShape *CreateHairShape(const Transform *o2w, const Transform *w2o,
	bool reverseOrientation, const ParamSet &params);

#endif /* PBRT_SHAPES_HAIR_H */