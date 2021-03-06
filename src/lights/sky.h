
/*
    pbrt source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.
    2011 Stefan Sydow (added sky model after Preetham1999)

    This file is part of pbrt.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_LIGHTS_SKY_H
#define PBRT_LIGHTS_SKY_H

// lights/sky.h*
#include "pbrt.h"
#include "light.h"
#include "texture.h"
#include "shape.h"
#include "scene.h"
#include "mipmap.h"

// SkyAreaLight Declarations
class SkyAreaLight : public Light {
public:
    // SkyAreaLight Public Methods
    SkyAreaLight(const Transform &light2world, int ns, 
        float turb, float lat, float hour, int day);
    ~SkyAreaLight();
    Spectrum Power(const Scene *) const;
    bool IsDeltaLight() const { return false; }
    Spectrum Le(const RayDifferential &r) const;
    Spectrum Sample_L(const Point &p, float pEpsilon, const LightSample &ls,
        float time, Vector *wi, float *pdf, VisibilityTester *visibility) const;
    Spectrum Sample_L(const Scene *scene, const LightSample &ls, float u1, float u2,
        float time, Ray *ray, Normal *Ns, float *pdf) const;
    float Pdf(const Point &, const Vector &) const;
    void SHProject(const Point &p, float pEpsilon, int lmax, const Scene *scene,
        bool computeLightVis, float time, RNG &rng, Spectrum *coeffs) const;

	void setEnvConditions(float turbidity, float local_lat);
	void setTime(int day_of_year, float hour_of_day);
	void updateEnv();

private:
    // SkyAreaLight Private Data
	
	float samplePower(float phi, float theta) const;
	Spectrum sample(float phi, float theta) const;

	float turbidity;
	float localLatitude;
	float hourOfDay;
	int dayOfYear;

	float phiSun, thetaSun;
	float declination; 
	float A[3],B[3],C[3],D[3],E[3];
	float xyY_zenith[3];
	float xyY_norm[3];

    Distribution2D *distribution;
};


SkyAreaLight *CreateSkyLight(const Transform &light2world,
        const ParamSet &paramSet);

#endif // PBRT_LIGHTS_SKY_H
