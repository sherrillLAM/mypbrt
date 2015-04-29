#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_SHAPES_HAIR_H
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
		const vector<Reference<Shape> > &cylinders);
	~HairShape();
	BBox ObjectBound() const;
	BBox WorldBound() const;

	bool CanIntersect() const { return false; }
	void Refine(vector<Reference<Shape> > &refined) const;
	void printShape() const;
	friend class Cylinder;

protected:
	vector<Reference<Shape> > cylinders;
};

HairShape *CreateHairShape(const Transform *o2w, const Transform *w2o,
	bool reverseOrientation, const vector<Reference<Shape> > &cylinders);

#endif /* PBRT_SHAPES_HAIR_H */