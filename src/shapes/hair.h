#if defined(_MSC_VER)
#pragma once
#endif

#if !defined(PBRT_SHAPES_HAIR_H)
#define PBRT_SHAPES_HAIR_H

// shapes/hair.h*
#include "shape.h"
#include "shapes/cylinder.h"

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
		vector<Cylinder*> &cylinders, float radius);
	~HairShape();
	BBox ObjectBound() const;
	BBox WorldBound() const;

	bool Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
		DifferentialGeometry *dg) const;
	bool IntersectP(const Ray &ray) const;

	Point Sample(float u1, float u2, Normal *Ns) const;

protected:
	vector<Cylinder*> cylinders;
	float radius;
};

HairShape *CreateHairShape(const Transform *o2w, const Transform *w2o,
	bool reverseOrientation, vector<Cylinder*> cylinders, float radius);

#endif /* PBRT_SHAPES_HAIR_H */