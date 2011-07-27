
/*
    pbrt source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.

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


// lights/sky.cpp*
#include "stdafx.h"
#include "lights/sky.h"
#include "sh.h"
#include "montecarlo.h"
#include "paramset.h"
#include "imageio.h"

float A_T[3] = { 0.1787, -0.0193, -0.0167};
float A_0[3] = {-1.4630, -0.2592, -0.2608};
float B_T[3] = {-0.3554, -0.0665, -0.0950};
float B_0[3] = { 0.4275,  0.0008,  0.0092};
float C_T[3] = {-0.0227, -0.0004, -0.0079};
float C_0[3] = { 5.3251,  0.2125,  0.2102};
float D_T[3] = { 0.1206, -0.0641, -0.0441};
float D_0[3] = {-2.5771, -0.8989, -1.6537};
float E_T[3] = {-0.0670, -0.0033, -0.0109};
float E_0[3] = { 0.3703,  0.0452,  0.0529};

// SkyAreaLight Utility Classes
struct SkyAreaCube {
    // SkyAreaCube Public Methods
    SkyAreaCube(const SkyAreaLight *l, const Scene *s,
                     float t, bool cv, float pe)
        : light(l), scene(s), time(t), pEpsilon(pe), computeVis(cv) { }
    Spectrum operator()(int, int, const Point &p, const Vector &w) {
        Ray ray(p, w, pEpsilon, INFINITY, time);
        if (!computeVis || !scene->IntersectP(ray))
            return light->Le(RayDifferential(ray));
        return 0.f;
    }
    const SkyAreaLight *light;
    const Scene *scene;
    float time, pEpsilon;
    bool computeVis;
};


// SkyAreaLight Method Definitions
SkyAreaLight::~SkyAreaLight() {
    delete distribution;
}

SkyAreaLight::SkyAreaLight(const Transform &light2world, int ns)
    : Light(light2world, ns) {
    int width = 1024, height = 1024;

	this->turbidity = 6; // clear sky;
	this->localLatitude = (52.0)*M_PI/180; //berlin
	this->hourOfDay = 4.00;
	this->dayOfYear = 160;
	updateEnv();

    // Initialize sampling PDFs for sky area light
    // Compute scalar-valued image _img_ from environment map
    float filter = 1.f / max(width, height);
    float *img = new float[width*height];
    for (int v = 0; v < height; ++v) {
        //float vp = (float)v / (float)height;
        float sinTheta = sinf(M_PI * float(v+.5f)/float(height));
        for (int u = 0; u < width; ++u) {
            //float up = (float)u / (float)width;
            img[u+v*width] = filter * samplePower(
					2*M_PI * float(u+.5f)/float(width), 
					M_PI * float(v+.5f)/float(height));
            img[u+v*width] *= sinTheta;
        }
    }

    // Compute sampling distributions for rows and columns of image
    distribution = new Distribution2D(img, width, height);
    delete[] img;
}

void SkyAreaLight::setEnvConditions(float turbidity, float local_lat){
	this->turbidity = turbidity;
	this->localLatitude = local_lat;
	updateEnv();
}

void SkyAreaLight::setTime(int day_of_year, float hour_of_day){
	this->hourOfDay = hour_of_day;
	this->dayOfYear = day_of_year;
	updateEnv();
}

void SkyAreaLight::updateEnv(){
	// calulate sun position 
	declination = 0.4093*sin((2*M_PI*(dayOfYear - 81))/368);
	float cos_l = cos(localLatitude);
	float cos_d = cos(declination);
	float sin_l_d = sin(localLatitude)*sin(declination);
	float time_angle = M_PI*hourOfDay/12;
	float cos_time_angle = cos(time_angle);
	thetaSun = M_PI/2 - asin(sin_l_d - cos_l*cos_d * cos_time_angle);
	phiSun = atan((-cos_d*sin(time_angle))/
			(cos_l*sin(declination) - sin(localLatitude)*cos_d * cos_time_angle));

	// compute coeff according to turbidity
	for(int i = 0; i < 3; i++){
		A[i] = A_0[i]+turbidity*A_T[i];
		B[i] = B_0[i]+turbidity*B_T[i];
		C[i] = C_0[i]+turbidity*C_T[i];
		D[i] = D_0[i]+turbidity*D_T[i];
		E[i] = E_0[i]+turbidity*E_T[i];
	}

	// compote zenith color
	float xi = (4.0/9.0 - turbidity/120.0)*(M_PI - 2*thetaSun);
	float thetaSun_2 = thetaSun*thetaSun;
	float thetaSun_3 = thetaSun*thetaSun_2;
	xyY_zenith[0] = (4.0453*turbidity - 4.9710)*tan(xi) - 0.2155*turbidity + 2.4192;
	xyY_zenith[1] = 
		turbidity*turbidity*
			( 0.0017*thetaSun_3 - 0.0037*thetaSun_2 + 0.0021*thetaSun + 0.000)
		+ turbidity*
			(-0.0290*thetaSun_3 + 0.0638*thetaSun_2 - 0.0320*thetaSun + 0.0039)
		+   ( 0.1169*thetaSun_3 - 0.2120*thetaSun_2 + 0.0605*thetaSun + 0.2589);
	xyY_zenith[2] = 
		turbidity*turbidity*
			( 0.0028*thetaSun_3 - 0.0061*thetaSun_2 + 0.0032*thetaSun + 0.000)
		+ turbidity*
			(-0.0421*thetaSun_3 + 0.0897*thetaSun_2 - 0.0415*thetaSun + 0.0052)
		+   ( 0.1535*thetaSun_3 - 0.2676*thetaSun_2 + 0.0667*thetaSun + 0.2669);

	float cos_thetaSun = cos(thetaSun);
	for (int i = 0; i < 3; i++){
		float Y_0 = (1 + A[i] * exp(B[i]))*(1 + C[i]*exp(D[i]*thetaSun) + E[i]*cos_thetaSun*cos_thetaSun);
		xyY_norm[i] = xyY_zenith[i]/Y_0;
	}
}

float SkyAreaLight::samplePower(float phi, float theta) const
{
	double cos_gamma = cos(theta)*cos(thetaSun) + sin(theta)*sin(thetaSun)*cos(phi-phiSun);
	double gamma = acos(cos_gamma);
	const double gamma_sun = 0.5*1.392e6/1.496e8;
	if(gamma < gamma_sun){//direct sunlight
		//gamma *= gamma/gamma_sun; //its a hack
		gamma = 0; //its a hack
		cos_gamma = cos(gamma);
	}

	float result = xyY_norm[0] * (1 + A[0] * exp(B[0]/fabs(cos(theta))))*(1 + C[0]*exp(D[0]*gamma) + E[0]*cos_gamma*cos_gamma);
	return result;
}

Spectrum SkyAreaLight::sample(float phi, float theta) const
{
	float xyY[3];
	float XYZ[3];
	RGBSpectrum result;
	double cos_gamma = cos(theta)*cos(thetaSun) + sin(theta)*sin(thetaSun)*cos(phi-phiSun);
	double gamma = acos(cos_gamma);
	float abs_cos_theta = fabs(cos(theta));
	const double gamma_sun = 0.5*1.392e6/1.496e8;
	if(gamma < gamma_sun){//direct sunlight
		//gamma *= gamma/gamma_sun; //its a hack
		gamma = 0; //its a hack
		cos_gamma = cos(gamma);
	}
	for(int i = 0; i < 3; i++){
		xyY[i] = xyY_norm[i] *(1 + A[i] * exp(B[i]/abs_cos_theta))*(1 + C[i]*exp(D[i]*gamma) + E[i]*cos_gamma*cos_gamma);
	}

	XYZ[0] = xyY[0]*xyY[1]/xyY[2];
	XYZ[1] = xyY[0];
	XYZ[2] = xyY[0]*(1-xyY[1]-xyY[2])/xyY[2];
	//assert(XYZ[2] >= 0);
	assert(xyY[0] < 1e6);
	return RGBSpectrum::FromXYZ(XYZ,SPECTRUM_ILLUMINANT);
	//return Spectrum(result,SPECTRUM_ILLUMINANT);
}

Spectrum SkyAreaLight::Power(const Scene *scene) const {
    Point worldCenter;
    float worldRadius;
	float XYZ[3];
	RGBSpectrum result;
	XYZ[0] = xyY_zenith[0]*xyY_zenith[1]/xyY_zenith[2];
	XYZ[1] = xyY_zenith[0];
	XYZ[2] = xyY_zenith[0]*(1-xyY_zenith[1]-xyY_zenith[2])/xyY_zenith[2];
	result.FromXYZ(XYZ);

    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    return M_PI * worldRadius * worldRadius *
		Spectrum(result,SPECTRUM_ILLUMINANT);
}

Spectrum SkyAreaLight::Le(const RayDifferential &r) const {
    Vector wh = Normalize(WorldToLight(r.d));
    float s = SphericalPhi(wh);
    float t = SphericalTheta(wh);
    return sample(s, t);
}


void SkyAreaLight::SHProject(const Point &p, float pEpsilon,
        int lmax, const Scene *scene, bool computeLightVis,
        float time, RNG &rng, Spectrum *coeffs) const {
    // Project _SkyAreaLight_ to SH using Monte Carlo if visibility needed
    if (computeLightVis) {
        Light::SHProject(p, pEpsilon, lmax, scene, computeLightVis,
                         time, rng, coeffs);
        return;
    }
    for (int i = 0; i < SHTerms(lmax); ++i)
        coeffs[i] = 0.f;
    /* TODO
	int ntheta = radianceMap->Height(), nphi = radianceMap->Width();
    if (min(ntheta, nphi) > 50) {
        // Project _SkyAreaLight_ to SH from lat-long representation

        // Precompute $\theta$ and $\phi$ values for lat-long map projection
        float *buf = new float[2*ntheta + 2*nphi];
        float *bufp = buf;
        float *sintheta = bufp;  bufp += ntheta;
        float *costheta = bufp;  bufp += ntheta;
        float *sinphi = bufp;    bufp += nphi;
        float *cosphi = bufp;
        for (int theta = 0; theta < ntheta; ++theta) {
            sintheta[theta] = sinf((theta + .5f)/ntheta * M_PI);
            costheta[theta] = cosf((theta + .5f)/ntheta * M_PI);
        }
        for (int phi = 0; phi < nphi; ++phi) {
            sinphi[phi] = sinf((phi + .5f)/nphi * 2.f * M_PI);
            cosphi[phi] = cosf((phi + .5f)/nphi * 2.f * M_PI);
        }
        float *Ylm = ALLOCA(float, SHTerms(lmax));
        for (int theta = 0; theta < ntheta; ++theta) {
            for (int phi = 0; phi < nphi; ++phi) {
                // Add _SkyAreaLight_ texel's contribution to SH coefficients
                Vector w = Vector(sintheta[theta] * cosphi[phi],
                                  sintheta[theta] * sinphi[phi],
                                  costheta[theta]);
                w = Normalize(LightToWorld(w));
                Spectrum Le = Spectrum(radianceMap->Texel(0, phi, theta),
                                       SPECTRUM_ILLUMINANT);
                SHEvaluate(w, lmax, Ylm);
                for (int i = 0; i < SHTerms(lmax); ++i)
                    coeffs[i] += Le * Ylm[i] * sintheta[theta] *
                        (M_PI / ntheta) * (2.f * M_PI / nphi);
            }
        }

        // Free memory used for lat-long theta and phi values
        delete[] buf;
    }
    else {
	*/
        // Project _SkyAreaLight_ to SH from cube map sampling
        SHProjectCube(SkyAreaCube(this, scene, time, computeLightVis,
                                       pEpsilon),
                      p, 200, lmax, coeffs);
    //TODO}
}


SkyAreaLight *CreateSkyLight(const Transform &light2world,
        const ParamSet &paramSet) {
    // TODO
	//Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
    //Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    int nSamples = paramSet.FindOneInt("nsamples", 1);
    if (PbrtOptions.quickRender) nSamples = max(1, nSamples / 4);
    return new SkyAreaLight(light2world, nSamples);
}


Spectrum SkyAreaLight::Sample_L(const Point &p, float pEpsilon,
        const LightSample &ls, float time, Vector *wi, float *pdf,
        VisibilityTester *visibility) const {
    PBRT_INFINITE_LIGHT_STARTED_SAMPLE();
    // Find $(u,v)$ sample coordinates in sky light texture
    float uv[2], mapPdf;
    distribution->SampleContinuous(ls.uPos[0], ls.uPos[1], uv, &mapPdf);
    if (mapPdf == 0.f) return 0.f;
	
	//float power = samplePower(2*M_PI * ls.uPos[0], M_PI *ls.uPos[1]);

    // Convert sky light sample point to direction
    float theta = uv[1] * M_PI, phi = uv[0] * 2.f * M_PI;
    float costheta = cosf(theta), sintheta = sinf(theta);
    float sinphi = sinf(phi), cosphi = cosf(phi);
    *wi = LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
                              costheta));

    // Compute PDF for sampled sky light direction
    *pdf = mapPdf / (2.f * M_PI * M_PI * sintheta);
    if (sintheta == 0.f) *pdf = 0.f;

    // Return radiance value for sky light direction
    visibility->SetRay(p, pEpsilon, *wi, time);
    Spectrum Ls = sample(phi, theta);
    PBRT_INFINITE_LIGHT_FINISHED_SAMPLE();
    return Ls;
}


float SkyAreaLight::Pdf(const Point &, const Vector &w) const {
    PBRT_INFINITE_LIGHT_STARTED_PDF();
    Vector wi = WorldToLight(w);
    float theta = SphericalTheta(wi), phi = SphericalPhi(wi);
    float sintheta = sinf(theta);
    if (sintheta == 0.f) return 0.f;
    float p = distribution->Pdf(phi * INV_TWOPI, theta * INV_PI) /
           (2.f * M_PI * M_PI * sintheta);
    PBRT_INFINITE_LIGHT_FINISHED_PDF();
    return p;
}


Spectrum SkyAreaLight::Sample_L(const Scene *scene,
        const LightSample &ls, float u1, float u2, float time,
        Ray *ray, Normal *Ns, float *pdf) const {
    PBRT_INFINITE_LIGHT_STARTED_SAMPLE();
    // Compute direction for sky light sample ray

    // Find $(u,v)$ sample coordinates in sky light texture
    float uv[2], mapPdf;
    distribution->SampleContinuous(ls.uPos[0], ls.uPos[1], uv, &mapPdf);
    if (mapPdf == 0.f) return Spectrum(0.f);

    float theta = uv[1] * M_PI, phi = uv[0] * 2.f * M_PI;
    float costheta = cosf(theta), sintheta = sinf(theta);
    float sinphi = sinf(phi), cosphi = cosf(phi);
    Vector d = -LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
                                    costheta));
    *Ns = (Normal)d;

    // Compute origin for sky light sample ray
    Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    Vector v1, v2;
    CoordinateSystem(-d, &v1, &v2);
    float d1, d2;
    ConcentricSampleDisk(u1, u2, &d1, &d2);
    Point Pdisk = worldCenter + worldRadius * (d1 * v1 + d2 * v2);
    *ray = Ray(Pdisk + worldRadius * -d, d, 0., INFINITY, time);

    // Compute _SkyAreaLight_ ray PDF
    float directionPdf = mapPdf / (2.f * M_PI * M_PI * sintheta);
    float areaPdf = 1.f / (M_PI * worldRadius * worldRadius);
    *pdf = directionPdf * areaPdf;
    if (sintheta == 0.f) *pdf = 0.f;
    Spectrum Ls = sample(phi, theta);
    PBRT_INFINITE_LIGHT_FINISHED_SAMPLE();
    return Ls;
}


