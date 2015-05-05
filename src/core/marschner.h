//////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2009, Image Engine Design Inc. All rights reserved.
//
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are
//  met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the name of Image Engine Design nor the names of any
//       other contributors to this software may be used to endorse or
//       promote products derived from this software without specific prior
//       written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
//  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
//  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
//  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//////////////////////////////////////////////////////////////////////////

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <cmath>
#include "pbrt.h"
#include "spectrum.h"
#include "geometry.h"

// Texture Inline Functions
inline float SmoothStep(float min, float max, float value) {
	float v = Clamp((value - min) / (max - min), 0.f, 1.f);
	return v * v * (-2.f * v + 3.f);
}

// Solves a * x + b == 0
float ieSolveLinear(float a, float b, float &root)
{
	float rootCount = -1;
	if (a != 0)
	{
		root = -b / a;
		rootCount = 1;
	}
	else if (b != 0)
	{
		rootCount = 0;
	}
	return rootCount;
}

float ieCubicRoot(float v)
{
	float result = pow(abs(v), 1 / 3);
	if (v < 0.f)
		result *= -1.f;

	return result;
}

float ieSolveQuadratic(float a, float b, float c, float *roots)
{
	const float epsilon = 1e-16;
	float rootCount = 0;
	if (abs(a) < epsilon)
	{
		rootCount = ieSolveLinear(b, c, roots[0]);
	}
	else
	{
		float D = b*b - 4 * a*c;

		if (abs(D) < epsilon)
		{
			roots[0] = -b / (2 * a);
			rootCount = 1;
		}
		else if (D > 0)
		{
			float s = sqrt(D);
			roots[0] = (-b + s) / (2 * a);
			roots[1] = (-b - s) / (2 * a);
			rootCount = 2;
		}
	}
	return rootCount;
}

// Computes real roots for a given cubic polynomial (x^3+Ax^2+Bx+C = 0).
// \todo: make sure it returns the same number of roots as in OpenEXR/ImathRoot.h
float ieSolveNormalizedCubic(float A, float B, float C, float *roots)
{
	const float epsilon = 1e-16;
	int rootCount = 0;
	if (abs(C) < epsilon)
	{
		// We're solving x^3 + A x^2 + Bx = 0. That's got a root where x = 0,
		// and potentially two more where x^2 + A x + B = 0:

		// find quadratic roots, if they exist:
		rootCount = ieSolveQuadratic(1, A, B, roots);

		// add x = 0 root:
		roots[rootCount] = 0;
		rootCount = rootCount + 1;
	}
	else
	{
		float Q = (3 * B - A*A) / 9;
		float R = (9 * A*B - 27 * C - 2 * A*A*A) / 54;
		float D = Q*Q*Q + R*R;	// polynomial discriminant

		if (D > 0) // complex or duplicate roots
		{
			float sqrtD = sqrt(D);
			float S = ieCubicRoot(R + sqrtD);
			float T = ieCubicRoot(R - sqrtD);
			roots[0] = (-A / 3 + (S + T));   // one real root
			rootCount = 1;
		}
		else  // 3 real roots
		{
			float th = acos(R / sqrt(-(Q*Q*Q)));
			float sqrtQ = sqrt(-Q);
			roots[0] = (2 * sqrtQ*cos(th / 3) - A / 3);
			roots[1] = (2 * sqrtQ*cos((th + 2 * M_PI) / 3) - A / 3);
			roots[2] = (2 * sqrtQ*cos((th + 4 * M_PI) / 3) - A / 3);
			rootCount = 3;
		}
	}
	return rootCount;
}

float ieSolveCubic(float a, float b, float c, float d, float *roots)
{
	float epsilon = 1e-16;
	float rootCount;
	if (abs(a) < epsilon)
	{
		rootCount = ieSolveQuadratic(b, c, d, roots);
	}
	else
	{
		rootCount = ieSolveNormalizedCubic(b / a, c / a, d / a, roots);
	}
	return rootCount;
}

// converts a given refraction index ( eta ) to work on a 2d plane that is a cross section of the hair.
// the theta parameter is the angle from the incident light to the cross section plane.
float ieBravaisIndex( float theta, float eta )
{
	float sinTheta = sin( theta );
	float result = sqrt( eta*eta - sinTheta*sinTheta ) / cos( theta );
	return result;
}

// Computes reflectance fresnel with different index (eta) of refractions for perpendicular and parallel polarized light.
// Assumes the source media is vaccuum ( n = 1 ). If invert is non-zero, then assumes the target media is vaccuum.
float ieMarschnerFresnel( float incidenceAngle, float etaPerp, float etaParal, const float invert )
{
	float n1, n2;
	float rPerp = 1;
	float rParal = 1;

	float angle = abs(incidenceAngle);
	if ( angle > M_PI/2 )
	{
		angle = M_PI - angle;
	}

	if ( invert )
	{
		n1 = etaPerp;
		n2 = 1;
	}
	else
	{
		n1 = 1;
		n2 = etaPerp;
	}

	// Perpendicular light reflectance
	float a = (n1/n2)*sin(angle);
	a *= a;
	if ( a <= 1 )
	{
		float b = n2*sqrt(1-a);
		float c = n1*cos(angle);
		rPerp =  ( c - b ) / ( c + b );
		rPerp *= rPerp;
		rPerp = fmin( 1, rPerp );
	}
	if ( invert )
	{
		n1 = etaParal;
		n2 = 1;
	}
	else
	{
		n1 = 1;
		n2 = etaParal;
	}
	// Parallel light reflectance
	float d = (n1/n2)*sin(angle);
	d *= d;
	if ( d <= 1 )
	{
		float e = n1*sqrt(1-d);
		float f = n2*cos(angle);
		rParal = ( e - f ) / ( e + f );
		rParal *= rParal;
		rParal = fmin( 1, rParal );
	}
	return 0.5 * (rPerp + rParal);
}

// computes a new refraction index based on the hair eccentricity and the azimuth distance.
float ieMarschnerEccentricityRefraction( float eccentricity, float refraction, float averageAzimuth )
{
	float n1 = 2 * ( refraction - 1 ) * eccentricity * eccentricity - refraction + 2;
	float n2 = 2 * ( refraction - 1 ) / ( eccentricity * eccentricity ) - refraction + 2;
	return ( (n1 + n2) + cos( 2 * averageAzimuth ) * ( n1 - n2 ) ) / 2;
}

float ieMarschnerExitAnglePolynomial( float p, float eta, float h )
{
	// use polynomial that approximates original equation.
	const float pi3 = M_PI*M_PI*M_PI;
	float gamma = asin(h);
	float c = asin(1/eta);
	return (6*p*c/M_PI - 2)*gamma-8*(p*c/pi3)*gamma*gamma*gamma + p*M_PI;
}

float ieMarschnerDExitAnglePolynomial( float p, float eta, float h )
{
	// computes the derivative of the polynomial relative to h.
	float gamma = asin( h );
	const float pi3 = M_PI*M_PI*M_PI;
	float c = asin(1/eta);
	float dGamma = (6*p*c/M_PI-2) - 3*8*(p*c/pi3)*gamma*gamma;
	float denom = sqrt(1-h*h);
	return dGamma/fmax(1e-5,denom);
}

float ieMarschnerDDExitAnglePolynomial( float p, float eta, float h )
{
	// computes the second derivative of the polynomial relative to h.
	float gamma = asin( h );
	const float pi3 = M_PI*M_PI*M_PI;
	float c = asin(1/eta);
	float dGamma = -2*3*8*(p*c/pi3)*gamma;
	float denom = pow( 1-h*h, 3/2 );
	return (dGamma*h)/fmax(1e-5,denom);
}

Spectrum ieMarschnerA(Spectrum absorption, Vector lightVec, const float p, float gammaI, float refraction, float etaPerp, float etaParal )
{
	if ( p == 0 )
	{
		return Spectrum (ieMarschnerFresnel( gammaI, etaPerp, etaParal, 0 ));
	}

	float h = sin( gammaI );			// from [1] right before equation 3.
	float gammaT = asin( Clamp(h / etaPerp,-1.f,1.f) );	// from [1] right before equation 3.
	float thetaT = acos( (etaPerp/refraction)*cos( lightVec[1] ) );	// definition for equation 20 in [2].
	float cosTheta = cos(thetaT);
	float l;
	// equation 20 in [2]
	l = 2*cos(gammaT)/fmax(1e-5,cosTheta);
	float rgb[3];
	absorption.ToRGB(rgb);
	Spectrum segmentAbsorption = Exp(absorption * l*p*-1.f);
	//Spectrum segmentAbsorption = Spectrum(expf(rgb[0] * l*p*-1.f), expf(rgb[1] * l*p*-1.f), expf(rgb[2] * l*p*-1.f));
	float fresnel;
	// equations 24-28 in [2]
	float invFresnel = ieMarschnerFresnel( gammaT, etaPerp, etaParal, 1 );
	fresnel = ( 1 - ieMarschnerFresnel( gammaI, etaPerp, etaParal, 0 ) )*( 1 - invFresnel );
	if ( p > 1 )
	{
		fresnel = fresnel * invFresnel;
	}
	return fresnel * segmentAbsorption;
}

float ieMarschnerTargetAngle( const float p,  float relativeAzimuth )
{
	float targetAngle = abs(relativeAzimuth);

	// set right range to match polynomial representation of the real curve.
	if ( p != 1 )
	{
		// converts angles to range [-PI,PI]
		if ( targetAngle > M_PI )
			targetAngle -= 2* M_PI;

		// offset center
		targetAngle += p* M_PI;
	}
	return targetAngle;
}

float ieMarschnerRoots( const float p,  float eta,  float targetAngle, float *roots )
{
	float rootCount;
	// Computes the roots of: o(p,y) - targetAngle = 0
	// by using the polynomial approximation: o(p,y) = (6pc / PI - 2)y - 8(pc/PI^3)y^3 + pPI where c = asin( 1/eta )^-1
	const float pi3 = M_PI*M_PI*M_PI;
	float c = asin(1/eta);
	rootCount = ieSolveCubic(-8 * (p*c / pi3), 0, (6 * p*c / M_PI - 2), (p*M_PI - targetAngle), roots);
	return rootCount;
}

Spectrum ieMarschnerNP(  Spectrum absorption,  Vector lightVec, const float p,  float refraction,  float etaPerp,  float etaParal,  float targetAngle )
{
	float roots[3] = { 0,0,0 };
	float rootCount = ieMarschnerRoots( p, etaPerp, targetAngle, roots );

	Spectrum result = Spectrum(0.f);
	const float denomMin = 1e-5;
	int rootIndex;
	for ( rootIndex = 0; rootIndex < rootCount; rootIndex+=1 )
	{
		float gammaI = roots[rootIndex];
		if ( abs(gammaI) <= M_PI/2 )
		{
			float h = sin( gammaI );
			Spectrum finalAbsorption = ieMarschnerA( absorption, lightVec, p, gammaI, refraction, etaPerp, etaParal );
			float dexitAngle;
			dexitAngle = ieMarschnerDExitAnglePolynomial( p, etaPerp, h );
			float denom = fmax( denomMin, 2*abs( dexitAngle ) );
			result += (finalAbsorption / denom);
		}
	}
	return result;
}

Spectrum ieMarschnerNTRT(Spectrum absorption, Vector lightVec, float refraction, float etaPerp, float etaParal, float targetAngle, const float causticLimit,const float causticWidth, const float glintScale, const float causticFade )
{
	float dH, t, hc, Oc1, Oc2;
	if ( etaPerp < 2 )
	{
		float ddexitAngle;
		// compute roots of the polynomials derivative
		float c = asin(1/etaPerp);
		const float pi3 = M_PI*M_PI*M_PI;
		float gammac = sqrt((6 * 2 * c / M_PI - 2) / (3 * 8 * (2 * c / pi3)));
		hc = abs(sin(gammac));
		ddexitAngle = ieMarschnerDDExitAnglePolynomial( 2, etaPerp, hc );
		dH = fmin( causticLimit, 2*sqrt( 2*causticWidth/abs( ddexitAngle ) ) );
		t = 1;
	}
	else
	{
		hc = 0;
		dH = causticLimit;
		t = 1 - SmoothStep( 2, 2 + causticFade, etaPerp );
	}

	Oc1 = ieMarschnerExitAnglePolynomial( 2, etaPerp, hc );
	Oc2 = ieMarschnerExitAnglePolynomial( 2, etaPerp, -hc );

	float a = 1 / (causticWidth * sqrt(2 * M_PI));

	boost::variate_generator<boost::mt19937, boost::normal_distribution<float> >
		causticCenter(boost::mt19937(time(0)),
		boost::normal_distribution<float>(0, causticWidth));

	boost::variate_generator<boost::mt19937, boost::normal_distribution<float> >
		causticLeft(boost::mt19937(time(0)),
		boost::normal_distribution<float>(-targetAngle + Oc1, causticWidth));

	boost::variate_generator<boost::mt19937, boost::normal_distribution<float> >
		causticRight(boost::mt19937(time(0)),
		boost::normal_distribution<float>(-targetAngle + Oc2, causticWidth));

	Spectrum glintAbsorption = ieMarschnerA( absorption, lightVec, 2, asin(hc), refraction, etaPerp, etaParal );

	Spectrum result = ieMarschnerNP( absorption, lightVec, 2, refraction, etaPerp, etaParal, targetAngle );
	result *= 1 - t*causticLeft()/causticCenter();
	result *= 1 - t*causticRight()/causticCenter();
	result += glintAbsorption*t*glintScale*dH*(causticLeft() + causticRight());
	return result;
}

float ieGaussian(float a, float b, float c, float x) {
	float o = x - b;
	return a * expf(-o*o / (2*c*c));
}

void ieGaussianPDF(float mu, float sigma, float &a, float &b, float &c) {
	a = 1/(sigma*sqrt(2*M_PI));
	b = mu;
	c = sigma;
}

float ieMarschnerM(float alpha, float beta, float normWidth, float x) {
	float a, b, c;
	float norm = 1/((normWidth/180.0f*M_PI)*sqrt(2*M_PI));
	ieGaussianPDF(alpha, beta, a, b, c);
	return (ieGaussian(a, b, c , x) / a) * norm;
}
