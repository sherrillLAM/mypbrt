#include "shapes/hair.h"
#include "stdafx.h"
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
	float r, float angleThreshold)
	: Shape(o2w, w2o, ro) {

}

bool HairShape::Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
	DifferentialGeometry *dg) const {
	bool result = false;
	for (int i = 0; i < vertices.size() - 1; i++) {
		if (vertexStartsFiber[i + 1]) break;

		Point p1 = (*WorldToObject)(vertices[i]);
		Point p2 = (*WorldToObject)(vertices[i+1]);
		
	}
}

bool HairShape::IntersectP(const Ray &ray) const {
	
}

const std::vector<Point> &HairShape::getVertices() const {
	return vertices;
}

const std::vector<bool> &HairShape::getStartFiber() const {
	return vertexStartsFiber;
}

BBox HairShape::ObjectBound() const {

}

float HairShape::Area() const {

}

Point HairShape::Sample(float u1, float u2, Normal *Ns) const {

}