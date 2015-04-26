
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// shapes/cylinder.cpp*
#include "stdafx.h"
#include "shapes/cylinder.h"
#include "paramset.h"
#include <iostream>

// Cylinder Method Definitions
Cylinder::Cylinder(const Transform *o2w, const Transform *w2o, bool ro,
                   float rad, float z0, float z1, float pm)
    : Shape(o2w, w2o, ro) {
    radius = rad;
    zmin = min(z0, z1);
    zmax = max(z0, z1);
    phiMax = Radians(Clamp(pm, 0.0f, 360.0f));
}

BBox Cylinder::ObjectBound() const {
    Point p1 = Point(-radius, -radius, zmin);
    Point p2 = Point( radius,  radius, zmax);
    return BBox(p1, p2);
}


bool Cylinder::Intersect(const Ray &r, float *tHit, float *rayEpsilon,
                         DifferentialGeometry *dg) const {
    float phi;
    Point phit;
    // Transform _Ray_ to object space
    Ray ray;
    (*WorldToObject)(r, &ray);

    // Compute quadratic cylinder coefficients
    float A = ray.d.x*ray.d.x + ray.d.y*ray.d.y;
    float B = 2 * (ray.d.x*ray.o.x + ray.d.y*ray.o.y);
    float C = ray.o.x*ray.o.x + ray.o.y*ray.o.y - radius*radius;

    // Solve quadratic equation for _t_ values
    float t0, t1;
    if (!Quadratic(A, B, C, &t0, &t1))
        return false;

    // Compute intersection distance along ray
    if (t0 > ray.maxt || t1 < ray.mint)
        return false;
    float thit = t0;
    if (t0 < ray.mint) {
        thit = t1;
        if (thit > ray.maxt) return false;
    }

    // Compute cylinder hit point and $\phi$
    phit = ray(thit);
    phi = atan2f(phit.y, phit.x);
    if (phi < 0.) phi += 2.f*M_PI;

    // Test cylinder intersection against clipping parameters
    if (phit.z < zmin || phit.z > zmax || phi > phiMax) {
        if (thit == t1) return false;
        thit = t1;
        if (t1 > ray.maxt) return false;
        // Compute cylinder hit point and $\phi$
        phit = ray(thit);
        phi = atan2f(phit.y, phit.x);
        if (phi < 0.) phi += 2.f*M_PI;
        if (phit.z < zmin || phit.z > zmax || phi > phiMax)
            return false;
    }

    // Find parametric representation of cylinder hit
    float u = phi / phiMax;
    float v = (phit.z - zmin) / (zmax - zmin);

    // Compute cylinder $\dpdu$ and $\dpdv$
    Vector dpdu(-phiMax * phit.y, phiMax * phit.x, 0);
    Vector dpdv(0, 0, zmax - zmin);

    // Compute cylinder $\dndu$ and $\dndv$
    Vector d2Pduu = -phiMax * phiMax * Vector(phit.x, phit.y, 0);
    Vector d2Pduv(0, 0, 0), d2Pdvv(0, 0, 0);

    // Compute coefficients for fundamental forms
    float E = Dot(dpdu, dpdu);
    float F = Dot(dpdu, dpdv);
    float G = Dot(dpdv, dpdv);
    Vector N = Normalize(Cross(dpdu, dpdv));
    float e = Dot(N, d2Pduu);
    float f = Dot(N, d2Pduv);
    float g = Dot(N, d2Pdvv);

    // Compute $\dndu$ and $\dndv$ from fundamental form coefficients
    float invEGF2 = 1.f / (E*G - F*F);
    Normal dndu = Normal((f*F - e*G) * invEGF2 * dpdu +
                         (e*F - f*E) * invEGF2 * dpdv);
    Normal dndv = Normal((g*F - f*G) * invEGF2 * dpdu +
                         (f*F - g*E) * invEGF2 * dpdv);

    // Initialize _DifferentialGeometry_ from parametric information
    const Transform &o2w = *ObjectToWorld;
    *dg = DifferentialGeometry(o2w(phit), o2w(dpdu), o2w(dpdv),
                               o2w(dndu), o2w(dndv), u, v, this);

    // Update _tHit_ for quadric intersection
    *tHit = thit;

    // Compute _rayEpsilon_ for quadric intersection
    *rayEpsilon = 5e-4f * *tHit;
    return true;
}


bool Cylinder::IntersectP(const Ray &r) const {
    float phi;
    Point phit;
    // Transform _Ray_ to object space
    Ray ray;
    (*WorldToObject)(r, &ray);

    // Compute quadratic cylinder coefficients
    float A = ray.d.x*ray.d.x + ray.d.y*ray.d.y;
    float B = 2 * (ray.d.x*ray.o.x + ray.d.y*ray.o.y);
    float C = ray.o.x*ray.o.x + ray.o.y*ray.o.y - radius*radius;

    // Solve quadratic equation for _t_ values
    float t0, t1;
    if (!Quadratic(A, B, C, &t0, &t1))
        return false;

    // Compute intersection distance along ray
    if (t0 > ray.maxt || t1 < ray.mint)
        return false;
    float thit = t0;
    if (t0 < ray.mint) {
        thit = t1;
        if (thit > ray.maxt) return false;
    }

    // Compute cylinder hit point and $\phi$
    phit = ray(thit);
    phi = atan2f(phit.y, phit.x);
    if (phi < 0.) phi += 2.f*M_PI;

    // Test cylinder intersection against clipping parameters
    if (phit.z < zmin || phit.z > zmax || phi > phiMax) {
        if (thit == t1) return false;
        thit = t1;
        if (t1 > ray.maxt) return false;
        // Compute cylinder hit point and $\phi$
        phit = ray(thit);
        phi = atan2f(phit.y, phit.x);
        if (phi < 0.) phi += 2.f*M_PI;
        if (phit.z < zmin || phit.z > zmax || phi > phiMax)
            return false;
    }
    return true;
}


float Cylinder::Area() const {
    return (zmax-zmin) * phiMax * radius;
}


Cylinder *CreateCylinderShape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params) {
    float radius = params.FindOneFloat("radius", 1);
    float zmin = params.FindOneFloat("zmin", -1);
    float zmax = params.FindOneFloat("zmax", 1);
    float phimax = params.FindOneFloat("phimax", 360);

	printf("%f,%f,%f,%f\n", o2w->GetMatrix().m[0][0], o2w->GetMatrix().m[0][1], o2w->GetMatrix().m[0][2], o2w->GetMatrix().m[0][3]);
	printf("%f,%f,%f,%f\n", o2w->GetMatrix().m[1][0], o2w->GetMatrix().m[1][1], o2w->GetMatrix().m[1][2], o2w->GetMatrix().m[1][3]);
	printf("%f,%f,%f,%f\n", o2w->GetMatrix().m[2][0], o2w->GetMatrix().m[2][1], o2w->GetMatrix().m[2][2], o2w->GetMatrix().m[2][3]);
	printf("%f,%f,%f,%f\n\n", o2w->GetMatrix().m[3][0], o2w->GetMatrix().m[3][1], o2w->GetMatrix().m[3][2], o2w->GetMatrix().m[3][3]);
	
	//printf("%f,%f,%f,%f\n", w2o->GetInverseMatrix().m[0][0], w2o->GetInverseMatrix().m[0][1], w2o->GetInverseMatrix().m[0][2], w2o->GetInverseMatrix().m[0][3]);
	//printf("%f,%f,%f,%f\n", w2o->GetInverseMatrix().m[1][0], w2o->GetInverseMatrix().m[1][1], w2o->GetInverseMatrix().m[1][2], w2o->GetInverseMatrix().m[1][3]);
	//printf("%f,%f,%f,%f\n", w2o->GetInverseMatrix().m[2][0], w2o->GetInverseMatrix().m[2][1], w2o->GetInverseMatrix().m[2][2], w2o->GetInverseMatrix().m[2][3]);
	//printf("%f,%f,%f,%f\n\n", w2o->GetInverseMatrix().m[3][0], w2o->GetInverseMatrix().m[3][1], w2o->GetInverseMatrix().m[3][2], w2o->GetInverseMatrix().m[3][3]);

	int np1;
	const Point *p = params.FindPoint("p", &np1);

	Point p1 = (*w2o)(Point(0.0f, 0.0f, zmin));
	Point p2 = (*w2o)(Point(0.0f, 0.0f, zmax));

	//printf("%f, %f, %f\n", p1.x, p1.y, p1.z);
	//printf("%f, %f, %f\n", p2.x, p2.y, p2.z);

	if (!p || np1 < 2) {
		printf("no points input.\n");
		return new Cylinder(o2w, w2o, reverseOrientation, radius, zmin, zmax, phimax);
	}

	Vector rel = p[1] - p[0];
	float length = rel.Length();

	Transform ObjectToWorld = fromFrame(rel / length) * (*o2w);
	Transform WorldToObject = Inverse(ObjectToWorld);

	//Point p3 = WorldToObject(p[0]);
	//Point p4 = WorldToObject(p[1]);

	//printf("%f, %f, %f\n", p3.x, p3.y, p3.z);
	//printf("%f, %f, %f\n", p4.x, p4.y, p4.z);

	printf("%f,%f,%f,%f\n", ObjectToWorld.GetMatrix().m[0][0], ObjectToWorld.GetMatrix().m[0][1], ObjectToWorld.GetMatrix().m[0][2], ObjectToWorld.GetMatrix().m[0][3]);
	printf("%f,%f,%f,%f\n", ObjectToWorld.GetMatrix().m[1][0], ObjectToWorld.GetMatrix().m[1][1], ObjectToWorld.GetMatrix().m[1][2], ObjectToWorld.GetMatrix().m[1][3]);
	printf("%f,%f,%f,%f\n", ObjectToWorld.GetMatrix().m[2][0], ObjectToWorld.GetMatrix().m[2][1], ObjectToWorld.GetMatrix().m[2][2], ObjectToWorld.GetMatrix().m[2][3]);
	printf("%f,%f,%f,%f\n\n", ObjectToWorld.GetMatrix().m[3][0], ObjectToWorld.GetMatrix().m[3][1], ObjectToWorld.GetMatrix().m[3][2], ObjectToWorld.GetMatrix().m[3][3]);

	//printf("%f,%f,%f,%f\n", WorldToObject.GetMatrix().m[0][0], WorldToObject.GetMatrix().m[0][1], WorldToObject.GetMatrix().m[0][2], WorldToObject.GetMatrix().m[0][3]);
	//printf("%f,%f,%f,%f\n", WorldToObject.GetMatrix().m[1][0], WorldToObject.GetMatrix().m[1][1], WorldToObject.GetMatrix().m[1][2], WorldToObject.GetMatrix().m[1][3]);
	//printf("%f,%f,%f,%f\n", WorldToObject.GetMatrix().m[2][0], WorldToObject.GetMatrix().m[2][1], WorldToObject.GetMatrix().m[2][2], WorldToObject.GetMatrix().m[2][3]);
	//printf("%f,%f,%f,%f\n", WorldToObject.GetMatrix().m[3][0], WorldToObject.GetMatrix().m[3][1], WorldToObject.GetMatrix().m[3][2], WorldToObject.GetMatrix().m[3][3]);

	//printf("%f,%f,%f,%f\n", w2o->GetInverseMatrix().m[0][0], w2o->GetInverseMatrix().m[0][1], w2o->GetInverseMatrix().m[0][2], w2o->GetInverseMatrix().m[0][3]);
	//printf("%f,%f,%f,%f\n", w2o->GetInverseMatrix().m[1][0], w2o->GetInverseMatrix().m[1][1], w2o->GetInverseMatrix().m[1][2], w2o->GetInverseMatrix().m[1][3]);
	//printf("%f,%f,%f,%f\n", w2o->GetInverseMatrix().m[2][0], w2o->GetInverseMatrix().m[2][1], w2o->GetInverseMatrix().m[2][2], w2o->GetInverseMatrix().m[2][3]);
	//printf("%f,%f,%f,%f\n\n", w2o->GetInverseMatrix().m[3][0], w2o->GetInverseMatrix().m[3][1], w2o->GetInverseMatrix().m[3][2], w2o->GetInverseMatrix().m[3][3]);

	//return new Cylinder(&ObjectToWorld, &WorldToObject, reverseOrientation, radius, 0.0f, length, phimax);
	std::cout << (o2w->GetMatrix() == ObjectToWorld.GetMatrix()) << std::endl;
	std::cout << (w2o->GetMatrix() == WorldToObject.GetMatrix()) << std::endl;

	//o2w = new Transform(ObjectToWorld.GetMatrix());
	//w2o = new Transform(WorldToObject.GetMatrix());
	return new Cylinder(o2w, w2o, reverseOrientation, radius, zmin, zmax, phimax);
}


Point Cylinder::Sample(float u1, float u2, Normal *Ns) const {
    float z = Lerp(u1, zmin, zmax);
    float t = u2 * phiMax;
    Point p = Point(radius * cosf(t), radius * sinf(t), z);
    *Ns = Normalize((*ObjectToWorld)(Normal(p.x, p.y, 0.)));
    if (ReverseOrientation) *Ns *= -1.f;
    return (*ObjectToWorld)(p);
}
