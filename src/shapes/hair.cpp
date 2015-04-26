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


class HairKDTree : public KdTreeAccel {
public:

	HairKDTree(std::vector<Point> &vertices,
		std::vector<bool> &vertexStartsFiber, float radius)
		: m_radius(radius) {
		/* Take the supplied vertex & start fiber arrays (without copying) */
		m_vertices.swap(vertices);
		m_vertexStartsFiber.swap(vertexStartsFiber);
		m_hairCount = 0;

		/* Compute the index of the first vertex in each segment. */
		m_segIndex.reserve(m_vertices.size());
		for (size_t i = 0; i<m_vertices.size() - 1; i++) {
			if (m_vertexStartsFiber[i])
				m_hairCount++;
			if (!m_vertexStartsFiber[i + 1])
				m_segIndex.push_back((size_t)i);
		}
		m_segmentCount = m_segIndex.size();

		/* Optimization: replace all primitive indices by the
		associated vertex indices (this avoids an extra
		indirection during traversal later on) */
		for (size_t i = 0; i<m_indexCount; ++i)
			m_indices[i] = m_segIndex[m_indices[i]];

		/* Free the segIndex array, it is not needed anymore */
		std::vector<size_t>().swap(m_segIndex);
	}

	/// Return the list of vertices underlying the hair kd-tree
	inline const std::vector<Point> &getVertices() const {
		return m_vertices;
	}

	/**
	* Return a boolean list specifying whether a vertex
	* marks the beginning of a new fiber
	*/
	inline const std::vector<bool> &getStartFiber() const {
		return m_vertexStartsFiber;
	}

	/// Return the radius of the hairs stored in the kd-tree
	inline float getRadius() const {
		return m_radius;
	}

	/// Return the total number of segments
	inline size_t getSegmentCount() const {
		return m_segmentCount;
	}

	/// Return the total number of hairs
	inline size_t getHairCount() const {
		return m_hairCount;
	}

	/// Return the total number of vertices
	inline size_t getVertexCount() const {
		return m_vertices.size();
	}

	/// Intersect a ray with all segments stored in the kd-tree
	inline bool rayIntersect(const Ray &ray, float _mint, float _maxt,
		float &t, void *temp) const {
		float tempT = std::numeric_limits<float>::infinity();
		float mint, maxt;

		if (m_aabb.rayIntersect(ray, mint, maxt)) {
			if (_mint > mint) mint = _mint;
			if (_maxt < maxt) maxt = _maxt;

			if (EXPECT_TAKEN(maxt > mint)) {
				if (rayIntersectHavran<false>(ray, mint, maxt, tempT, temp)) {
					t = tempT;
					return true;
				}
			}
		}
		return false;
	}

	/**
	* \brief Intersect a ray with all segments stored in the kd-tree
	* (Visiblity query version)
	*/
	inline bool rayIntersect(const Ray &ray, float _mint, float _maxt) const {
		float tempT = std::numeric_limits<float>::infinity();
		float mint, maxt;

		if (m_aabb.rayIntersect(ray, mint, maxt)) {
			if (_mint > mint) mint = _mint;
			if (_maxt < maxt) maxt = _maxt;

			if (EXPECT_TAKEN(maxt > mint)) {
				if (rayIntersectHavran<true>(ray, mint, maxt, tempT, NULL))
					return true;
			}
		}
		return false;
	}

#if defined(MTS_HAIR_USE_FANCY_CLIPPING)
	/**
	* Compute the ellipse created by the intersection of an infinite
	* cylinder and a plane. Returns false in the degenerate case.
	* Based on:
	* www.geometrictools.com/Documentation/IntersectionCylinderPlane.pdf
	*/
	bool intersectCylPlane(Point planePt, Normal planeNrml,
		Point cylPt, Vector cylD, float radius, Point &center,
		Vector *axes, float *lengths) const {
		if (AbsDot(planeNrml, cylD) < Epsilon)
			return false;

		Assert(std::abs(planeNrml.Length() - 1) <Epsilon);
		Vector B, A = cylD - Dot(cylD, planeNrml)*planeNrml;

		float length = A.Length();
		if (length > Epsilon && planeNrml != cylD) {
			A /= length;
			B = Cross(planeNrml, A);
		}
		else {
			coordinateSystem(planeNrml, A, B);
		}

		Vector delta = planePt - cylPt,
			deltaProj = delta - cylD*Dot(delta, cylD);

		float aDotD = Dot(A, cylD);
		float bDotD = Dot(B, cylD);
		float c0 = 1 - aDotD*aDotD;
		float c1 = 1 - bDotD*bDotD;
		float c2 = 2 * Dot(A, deltaProj);
		float c3 = 2 * Dot(B, deltaProj);
		float c4 = Dot(delta, deltaProj) - radius*radius;

		float lambda = (c2*c2 / (4 * c0) + c3*c3 / (4 * c1) - c4) / (c0*c1);

		float alpha0 = -c2 / (2 * c0),
			beta0 = -c3 / (2 * c1);

		lengths[0] = std::sqrt(c1*lambda),
			lengths[1] = std::sqrt(c0*lambda);

		center = planePt + alpha0 * A + beta0 * B;
		axes[0] = A;
		axes[1] = B;
		return true;
	}

	/**
	* \brief Intersect an infinite cylinder with an
	* AABB face and bound the resulting clipped ellipse
	*/
	BBox intersectCylFace(int axis,
		const Point &min, const Point &max,
		const Point &cylPt, const Vector &cylD) const {
		int axis1 = (axis + 1) % 3;
		int axis2 = (axis + 2) % 3;

		Normal planeNrml(0.0f, 0.0f, 0.0f);
		planeNrml[axis] = 1;

		Point ellipseCenter;
		Vector ellipseAxes[2];
		float ellipseLengths[2];

		BBox aabb;
		if (!intersectCylPlane(min, planeNrml, cylPt, cylD, m_radius,
			ellipseCenter, ellipseAxes, ellipseLengths)) {
			/* Degenerate case -- return an invalid AABB. This is
			not a problem, since one of the other faces will provide
			enough information to arrive at a correct clipped AABB */
			return aabb;
		}

		/* Intersect the ellipse against the sides of the AABB face */
		for (int i = 0; i<4; ++i) {
			Point p1, p2;
			p1[axis] = p2[axis] = min[axis];
			p1[axis1] = ((i + 1) & 2) ? min[axis1] : max[axis1];
			p1[axis2] = ((i + 0) & 2) ? min[axis2] : max[axis2];
			p2[axis1] = ((i + 2) & 2) ? min[axis1] : max[axis1];
			p2[axis2] = ((i + 1) & 2) ? min[axis2] : max[axis2];

			Point p1l(
				Dot(p1 - ellipseCenter, ellipseAxes[0]) / ellipseLengths[0],
				Dot(p1 - ellipseCenter, ellipseAxes[1]) / ellipseLengths[1]);
			Point p2l(
				Dot(p2 - ellipseCenter, ellipseAxes[0]) / ellipseLengths[0],
				Dot(p2 - ellipseCenter, ellipseAxes[1]) / ellipseLengths[1]);

			Vector rel = p2l - p1l;
			float A = Dot(rel, rel);
			float B = 2 * Dot(Vector(p1l), rel);
			float C = Dot(Vector(p1l), Vector(p1l)) - 1;

			float x0, x1;
			if (solveQuadratic(A, B, C, x0, x1)) {
				if (x0 >= 0 && x0 <= 1)
					aabb.expandBy(p1 + (p2 - p1)*x0);
				if (x1 >= 0 && x1 <= 1)
					aabb.expandBy(p1 + (p2 - p1)*x1);
			}
		}

		ellipseAxes[0] *= ellipseLengths[0];
		ellipseAxes[1] *= ellipseLengths[1];
		BBox faceBounds(min, max);

		/* Find the componentwise maxima of the ellipse */
		for (int i = 0; i<2; ++i) {
			int j = (i == 0) ? axis1 : axis2;
			float alpha = ellipseAxes[0][j];
			float beta = ellipseAxes[1][j];
			float ratio = beta / alpha, tmp = std::sqrt(1 + ratio*ratio);
			float cosTheta = 1 / tmp, sinTheta = ratio / tmp;
			Point p1 = ellipseCenter + cosTheta*ellipseAxes[0] + sinTheta*ellipseAxes[1];
			Point p2 = ellipseCenter - cosTheta*ellipseAxes[0] - sinTheta*ellipseAxes[1];

			if (faceBounds.contains(p1))
				aabb.expandBy(p1);
			if (faceBounds.contains(p2))
				aabb.expandBy(p2);
		}

		return aabb;
	}

	BBox getAABB(int index) const {
		int iv = m_segIndex.at(index);
		Point center;
		Vector axes[2];
		float lengths[2];

		bool success = intersectCylPlane(firstVertex(iv), firstMiterNormal(iv),
			firstVertex(iv), tangent(iv), m_radius, center, axes, lengths);
		Assert(success);

		BBox result;
		axes[0] *= lengths[0]; axes[1] *= lengths[1];
		for (int i = 0; i<3; ++i) {
			float range = std::sqrt(axes[0][i] * axes[0][i] + axes[1][i] * axes[1][i]);
			result.min[i] = fmin(result.min[i], center[i] - range);
			result.max[i] = fmax(result.max[i], center[i] + range);
		}

		success = intersectCylPlane(secondVertex(iv), secondMiterNormal(iv),
			secondVertex(iv), tangent(iv), m_radius, center, axes, lengths);
		Assert(success);

		axes[0] *= lengths[0]; axes[1] *= lengths[1];
		for (int i = 0; i<3; ++i) {
			float range = std::sqrt(axes[0][i] * axes[0][i] + axes[1][i] * axes[1][i]);
			result.min[i] = fmin(result.min[i], center[i] - range);
			result.max[i] = fmax(result.max[i], center[i] + range);
		}
		return result;
	}

	BBox getClippedAABB(int index, const BBox &box) const {
		/* Compute a base bounding box */
		BBox base(getAABB(index));
		base.clip(box);

		int iv = m_segIndex.at(index);

		Point cylPt = firstVertex(iv);
		Vector cylD = tangent(iv);

		/* Now forget about the cylinder ends and
		intersect an infinite cylinder with each AABB face */
		BBox clippedAABB;
		clippedAABB.expandBy(intersectCylFace(0,
			Point(base.min.x, base.min.y, base.min.z),
			Point(base.min.x, base.max.y, base.max.z),
			cylPt, cylD));

		clippedAABB.expandBy(intersectCylFace(0,
			Point(base.max.x, base.min.y, base.min.z),
			Point(base.max.x, base.max.y, base.max.z),
			cylPt, cylD));

		clippedAABB.expandBy(intersectCylFace(1,
			Point(base.min.x, base.min.y, base.min.z),
			Point(base.max.x, base.min.y, base.max.z),
			cylPt, cylD));

		clippedAABB.expandBy(intersectCylFace(1,
			Point(base.min.x, base.max.y, base.min.z),
			Point(base.max.x, base.max.y, base.max.z),
			cylPt, cylD));

		clippedAABB.expandBy(intersectCylFace(2,
			Point(base.min.x, base.min.y, base.min.z),
			Point(base.max.x, base.max.y, base.min.z),
			cylPt, cylD));

		clippedAABB.expandBy(intersectCylFace(2,
			Point(base.min.x, base.min.y, base.max.z),
			Point(base.max.x, base.max.y, base.max.z),
			cylPt, cylD));

		clippedAABB.clip(base);
		return clippedAABB;
	}
#else
	/// Compute the AABB of a segment (only used during tree construction)
BBox getAABB(int index) const {
		int iv = m_segIndex.at(index);

		// cosine of steepest miter angle
		const float cos0 = Dot(firstMiterNormal(iv), tangent(iv));
		const float cos1 = Dot(secondMiterNormal(iv), tangent(iv));
		const float maxInvCos = 1.0 / std::min(cos0, cos1);
		const Vector expandVec(m_radius * maxInvCos);

		const Point a = firstVertex(iv);
		const Point b = secondVertex(iv);

		BBox aabb;
		aabb.expandBy(a - expandVec);
		aabb.expandBy(a + expandVec);
		aabb.expandBy(b - expandVec);
		aabb.expandBy(b + expandVec);
		return aabb;
	}

	/// Compute the clipped AABB of a segment (only used during tree construction)
BBox getClippedAABB(int index, const AABB &box) const {
	BBox aabb(getAABB(index));
		aabb.clip(box);
		return aabb;
	}
#endif

	/// Return the total number of segments
	inline int getPrimitiveCount() const {
		return (int)m_segIndex.size();
	}

	inline bool intersect(const Ray &ray, int iv,
		float mint, float maxt, float &t, void *tmp) const {
		/* First compute the intersection with the infinite cylinder */
		float nearT, farT;
		Vector axis = tangent(iv);

		// Projection of ray onto subspace normal to axis
		Vector relOrigin = ray.o - firstVertex(iv);
		Vector projOrigin = relOrigin - Dot(axis, relOrigin) * axis;
		Vector projDirection = ray.d - Dot(axis, ray.d) * axis;

		// Quadratic to intersect circle in projection
		const float A = projDirection.LengthSquared();
		const float B = 2 * Dot(projOrigin, projDirection);
		const float C = projOrigin.LengthSquared() - m_radius*m_radius;

		if (!solveQuadratic(A, B, C, nearT, farT))
			return false;

		if (nearT > maxt || farT < mint)
			return false;

		/* Next check the intersection points against the miter planes */
		Point pointNear = ray(nearT);
		Point pointFar = ray(farT);
		if (Dot(pointNear - firstVertex(iv), firstMiterNormal(iv)) >= 0 &&
			Dot(pointNear - secondVertex(iv), secondMiterNormal(iv)) <= 0 &&
			nearT >= mint) {
			t = nearT;
		}
		else if (Dot(pointFar - firstVertex(iv), firstMiterNormal(iv)) >= 0 &&
			Dot(pointFar - secondVertex(iv), secondMiterNormal(iv)) <= 0) {
			if (farT > maxt)
				return false;
			t = farT;
		}
		else {
			return false;
		}

		int *storage = static_cast<int *>(tmp);
		if (storage)
			*storage = iv;

		return true;
	}

	inline bool intersect(const Ray &ray, int iv,
		float mint, float maxt) const {
		float tempT;
		return intersect(ray, iv, mint, maxt, tempT, NULL);
	}

	/* Some utility functions */
	inline Point firstVertex(int iv) const { return m_vertices[iv]; }
	inline Point secondVertex(int iv) const { return m_vertices[iv + 1]; }
	inline Point prevVertex(int iv) const { return m_vertices[iv - 1]; }
	inline Point nextVertex(int iv) const { return m_vertices[iv + 2]; }

	inline bool prevSegmentExists(int iv) const { return !m_vertexStartsFiber[iv]; }
	inline bool nextSegmentExists(int iv) const { return !m_vertexStartsFiber[iv + 2]; }

	inline Vector tangent(int iv) const { return Normalize(secondVertex(iv) - firstVertex(iv)); }
	inline Vector prevTangent(int iv) const { return Normalize(firstVertex(iv) - prevVertex(iv)); }
	inline Vector nextTangent(int iv) const { return Normalize(nextVertex(iv) - secondVertex(iv)); }

	inline Vector firstMiterNormal(int iv) const {
		if (prevSegmentExists(iv))
			return Normalize(prevTangent(iv) + tangent(iv));
		else
			return tangent(iv);
	}

	inline Vector secondMiterNormal(int iv) const {
		if (nextSegmentExists(iv))
			return Normalize(tangent(iv) + nextTangent(iv));
		else
			return tangent(iv);
	}

protected:
	std::vector<Point> m_vertices;
	std::vector<bool> m_vertexStartsFiber;
	std::vector<int> m_segIndex;
	size_t m_segmentCount;
	size_t m_hairCount;
	float m_radius;
}; 

HairShape::HairShape(const Transform *o2w, const Transform *w2o, bool ro, 
	const string filename, float r, float angleThreshold)
	: Shape(o2w, w2o, ro) {
	/* Skip segments, whose tangent differs by less than one degree
	compared to the previous one */
	float dpThresh = std::cos(angleThreshold);

	std::ifstream is(filename);
	if (is.fail())
		Error("Could not open \"%s\"!", filename);

	radius = r;

	string line;
	bool newFiber = true;
	Point p, lastP(0.0f, 0.0f, 0.0f);
	Vector tangent(0.0f, 0.0f, 0.0f);
	size_t nDegenerate = 0, nSkipped = 0;

	while (is.good()) {
		getline(is, line);
		if (line.length() > 0 && line[0] == '#') {
			newFiber = true;
			continue;
		}
		std::istringstream iss(line);
		iss >> p.x >> p.y >> p.z;
		if (!iss.fail()) {
			p = (*o2w)(p);
			if (newFiber) {
				vertices.push_back(p);
				vertexStartsFiber.push_back(newFiber);
				lastP = p;
				tangent = Vector(0.0f, 0.0f, 0.0f);
			}
			else if (p != lastP) {
				if (tangent == Vector(0.0f, 0.0f, 0.0f)) {
					vertices.push_back(p);
					vertexStartsFiber.push_back(newFiber);
					tangent = Normalize(p - lastP);
					lastP = p;
				}
				else {
					Vector nextTangent = Normalize(p - lastP);
					if (Dot(nextTangent, tangent) > dpThresh) {
						/* Too small of a difference in the tangent value,
						just overwrite the previous vertex by the current one */
						tangent = Normalize(p - vertices[vertices.size() - 2]);
						vertices[vertices.size() - 1] = p;
						++nSkipped;
					}
					else {
						vertices.push_back(p);
						vertexStartsFiber.push_back(newFiber);
						tangent = nextTangent;
					}
					lastP = p;
				}
			}
			else {
				nDegenerate++;
			}
			newFiber = false;
		}
		else {
			newFiber = true;
		}
	}

	if (nDegenerate > 0)
		Error("Encountered %d degenerate segments!", nDegenerate);
	if (nSkipped > 0)
		Error("Skipped %d low-curvature segments.", nSkipped);

	vertexStartsFiber.push_back(true);
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