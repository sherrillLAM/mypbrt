
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

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_MATERIALS_HAIR_M_H
#define PBRT_MATERIALS_HAIR_M_H

// materials/hair.h*
#include "pbrt.h"
#include "material.h"

// HairMaterial Declarations
class HairMaterial : public Material {
public:
	// HairMaterial Public Methods
	HairMaterial(Reference<Texture<Spectrum> > kd,
		Reference<Texture<Spectrum> > ks,
		Reference<Texture<float> > roughness,
		Reference<Texture<Spectrum> > refl,
		Reference<Texture<Spectrum> > trans,
		Reference<Texture<float> > bump,
		Reference<Texture<float> > k,
		string mod,
		Reference<Texture<float> > refr,
		Reference<Texture<float> > ecc,
		Reference<Texture<float> > aR,
		Reference<Texture<float> > aTT,
		Reference<Texture<float> > aTRT,
		Reference<Texture<float> > bR,
		Reference<Texture<float> > bTT,
		Reference<Texture<float> > bTRT,
		Reference<Texture<Spectrum> > absorb)
		: Kd(kd), Ks(ks), roughness(roughness), reflect(refl),
		transmit(trans), bumpMap(bump), model(mod), K(k),
		refraction(refr), eccentricity(ecc),
		alphaR(aR), alphaTT(aTT), alphaTRT(aTRT),
		betaR(bR), betaTT(bTT), betaTRT(bTRT),
		absorb(absorb) {
	}
	BSDF *GetBSDF(const DifferentialGeometry &dgGeom,
		const DifferentialGeometry &dgShading,
		MemoryArena &arena) const;
private:
	// HairMaterial Private Data
	Reference<Texture<Spectrum> > Kd, Ks;
	Reference<Texture<Spectrum> > reflect, transmit, absorb;
	Reference<Texture<float> > roughness, bumpMap, K,
							refraction, eccentricity,
							alphaR, alphaTT, alphaTRT, betaR, betaTT, betaTRT;
	string model;
};


HairMaterial *CreateHairMaterial(const Transform &xform,
	const TextureParams &mp);

#endif // PBRT_MATERIALS_HAIR_M_H