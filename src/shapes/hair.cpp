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
	vector<Cylinder*> &cys, float r)
	: Shape(o2w, w2o, ro) {
	cylinders = cys;
	radius = r;
}

HairShape::~HairShape() {
	for (int i = 0; i < cylinders.size(); i++)
		delete cylinders[i];
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

bool HairShape::Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
	DifferentialGeometry *dg) const {
	return false;
}

bool HairShape::IntersectP(const Ray &ray) const {
	return false;
}

Point HairShape::Sample(float u1, float u2, Normal *Ns) const {
	return Point(0,0,0);
}

HairShape *CreateHairShape(const Transform *o2w, const Transform *w2o,
	bool reverseOrientation, vector<Cylinder *> cys, float r) {
	return new HairShape(o2w, w2o, reverseOrientation, cys, r);
}