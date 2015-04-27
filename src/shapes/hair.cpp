#include "stdafx.h"
#include "shapes/hair.h"
#include <fstream>
#include <sstream>

#define MTS_HAIR_USE_FANCY_CLIPPING 1

/**
* \brief Space-efficient acceleration structure for cylindrical hair
* segments with miter joints. This class expects an ASCII file containing
* a list of hairs made from segments. Each line should contain an X,
* Y and Z coordinate separated by a space. An empty line indicates
* the start of a new hair.
*/

HairShape::HairShape(const Transform *o2w, const Transform *w2o, bool ro,
	const vector<Reference<Shape> > &cys)
	: Shape(o2w, w2o, ro) {
	cylinders = cys;
}

HairShape::~HairShape() {
	cylinders.clear();
}

BBox HairShape::ObjectBound() const {
	BBox objectBounds;
	for (int i = 0; i < cylinders.size(); i++)
		objectBounds = Union(objectBounds, cylinders[i]->ObjectBound());
	return objectBounds;
}

BBox HairShape::WorldBound() const {
	BBox worldBounds;
	for (int i = 0; i < cylinders.size(); i++)
		worldBounds = Union(worldBounds, cylinders[i]->WorldBound());
	return worldBounds;
}

void HairShape::Refine(vector<Reference<Shape> > &refined) const {
	refined = cylinders;
}

HairShape *CreateHairShape(const Transform *o2w, const Transform *w2o,
	bool reverseOrientation, const vector<Reference<Shape> > &cys) {
	return new HairShape(o2w, w2o, reverseOrientation, cys);
}

void HairShape::printShape() const {
	printf("shape: Hair\n");
	printf("num of cylinders: %d\n", cylinders.size());
	fflush(stdout);
}