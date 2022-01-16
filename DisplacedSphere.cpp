
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

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


// shapes/dsphere.cpp*
#include "shapes/DisplacedSphere.h"
#include "sampling.h"
#include "paramset.h"
#include "efloat.h"
#include "stats.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>

// GLOBALS:
const float minDist = 0.000003; // minimum distance of raymarch before stop (hit)
const float minDisty = 0.03; // minimum distance of raymarch before stop (hit)
const float esix    = 0.000001;  // 1e^-6;

// PPM GLOBALS:
std::vector<unsigned int> r; // red
std::vector<unsigned int> g; // green
std::vector<unsigned int> b; // blue
std::vector<float> rgb;      // greyscale ((r+g+b)/3)

unsigned int width;
unsigned int height;
unsigned int maxColVal = 255;
int size;

bool analytical; // use DisplacementAnalytical if true else DisplacementMap

namespace pbrt 
{

    ////////////////////////////////////
    // dsphere specific added functions:
    ////////////////////////////////////

    // Normalize but it prevents divide by zero:
    Vector3f Normalise(Vector3f abnormal)
    {
        if (abnormal == Vector3f({0, 0, 0})) 
        {
            return abnormal;
        }
        Vector3f normal = Normalize(abnormal);
        return normal;
    }

    // https://en.wikipedia.org/wiki/UV_mapping
    // https://glslsandbox.com/e#37633.02
    Float* Dsphere::VectorToUV(Vector3f input) const
    {
        Float uv[2];
        input = Normalise(input);

        uv[0] = 0.5 + (atan2(input.x, input.z) / (2 * Pi));
        uv[1] = 0.5 - (asin(input.y) / (Pi));

        return uv;
    }

    // Displacement as given in coursework document:
    Float Dsphere::DisplacementAnalytical(Vector3f vector) const 
    {
        Float *uv;

        uv = VectorToUV(-Normalise(vector));

        return (maxDisplacement / 2) * (1 + cos(50 * Pi * uv[0]) * sin(25 * Pi * uv[1]));
    }

    // gives magnitude of displacement at the given sphere coords
    // according to the displacement map
    Float Dsphere::DisplacementMap(Vector3f vector) const
    {
        Float *uv;

        uv = VectorToUV(-Normalise(vector));

        float x = uv[0] * width;
        float y = uv[1] * height;

        int u = ceil(x);  // round to nearest pixel
        int v = ceil(y);

        // Find difference between rounded value
        // and unrounded value
        float du  = u - x;
        float dv  = v - y;
        float duv = std::abs(du - dv);

        const int indexi = (u + v * width) % size; // index in the image 1d array
        const int indexj = (indexi + 1) % size;

        // take weighted sum of the 2 pixels that x and y fall between:
        float rgbi =  ((rgb[indexi] * duv) + ((1 - duv) * rgb[indexj])) * (maxDisplacement / maxColVal);

        // displacement at this index:
        return rgbi;

    }

    // Determines which displacement function to use:
    Float Dsphere::Displacement(Vector3f vector) const 
    {
        return analytical ? DisplacementAnalytical(vector) : DisplacementMap(vector);
    }

    // gives distance to the displacement sphere:
    Float Dsphere::DistanceToDsphere(Vector3f vector) const
    {
        return ((vector.Length() - radius) - Displacement(vector));
    }

    // Implements the ray marching algorithm.
    // Returns an approximate distance travelled by ray to the shape.
    // https://www.shadertoy.com/view/4ldGzB
    Float Dsphere::RayMarch(Ray ray, bool* intersect) const
    {

        Vector3f step = Vector3f(ray(0));

        float max          = step.Length(); // 
        float distance     = 0;             // distance travelled along ray
        float stepDistance = max;           // distance travelled in this step

        // while stepDistance is not too great and not too small
        while (stepDistance <= max * 2 && stepDistance > minDist)
        {
            step = Vector3f(ray(distance)); // take as step along the ray
            
            stepDistance = step.Length() - radius - Displacement(step);  // find distance to dsphere
            
            distance += stepDistance; // move along the ray
        }

        *intersect = false;
        if (stepDistance <= minDist) *intersect = true;

        return distance;
    }

    // https://github.com/sol-prog/Perlin_Noise/blob/master/ppm.cpp
    void ReadPPM(const std::string &fname) 
    {
        std::ifstream inp(fname.c_str(), std::ios::in | std::ios::binary);
        if (inp.is_open()) 
        {
            std::string line;
            std::getline(inp, line);
            if (line != "P6") 
            {
                std::cout << "Error. Unrecognized file format." << std::endl;
                return;
            }
            std::getline(inp, line);
            while (line[0] == '#') 
            {
                std::getline(inp, line);
            }
            std::stringstream dimensions(line);

            try 
            {
                dimensions >> width;
                dimensions >> height;
            } 
            catch (std::exception &e) 
            {
                std::cout << "Header file format error. " << e.what() << std::endl;
                return;
            }

            std::getline(inp, line);
            std::stringstream max_val(line);
            try 
            {
                max_val >> maxColVal;
            } 
            catch (std::exception &e) 
            {
                std::cout << "Header file format error. " << e.what()
                          << std::endl;
                return;
            }
            
            size = width * height;
            r.reserve(size);
            g.reserve(size);
            b.reserve(size);
            rgb.reserve(size);

            char aux;
            for (unsigned int i = 0; i < size; i++) 
            {
                inp.read(&aux, 1);
                r.push_back((unsigned char)aux);
                inp.read(&aux, 1);
                g.push_back((unsigned char)aux);
                inp.read(&aux, 1);
                b.push_back((unsigned char)aux);

                rgb.push_back((r[i] + g[i] + b[i]) / 3);
            }
            
        } 

        else 
        {
            std::cout << "Error. Unable to open " << fname << std::endl;
        }

        inp.close();

        
        // convulution matrix:
        int gaussV[25] = {1,  4, 6, 4,  1,  4,  16, 24, 16, 4, 6, 24, 36,
                          24, 6, 4, 16, 24, 16, 4,  1,  4,  6, 4, 1};

        int gaussVII[49] = {1,  4,  7,  10, 7,  4,  1,  4,  12, 26, 33, 26, 12,
                            4,  7,  26, 55, 71, 55, 26, 7,  10, 33, 71, 91, 71,
                            33, 10, 7,  26, 55, 71, 55, 26, 7,  4,  12, 26, 33,
                            26, 12, 4,  1,  4,  7,  10, 7,  4,  1};

        int normaliser = 0;
        for (auto &num : gaussVII) normaliser += num;

        for (int i = 0; i < size; i++) 
        {
            int result = 0;
            for (int j = -3; j <= 3; j++) 
            {
                for (int k = -3; k <= 3; k++) 
                {
                    int index = (i + j * width + k) % size;
                    result += rgb[index] * gaussVII[(j + 3) * 7 + k + 3];
                }
            }

            rgb[i] = result / normaliser;
        }

    }

    ////////////////////////////////////
    // modified functions:
    ////////////////////////////////////

    bool Dsphere::Intersect(const Ray &r, Float *tHit,
                            SurfaceInteraction *isect,
                            bool testAlphaTexture) const 
    {

        ProfilePhase p(Prof::ShapeIntersect);
        // Transform _Ray_ to object space
        Vector3f oErr, dErr;
        Ray ray = (*WorldToObject)(r, &oErr, &dErr);

        // As in sphere.cpp, we check to see if a ray hits sphere
        // with radius = radius + maxDisplacement
        EFloat ox(ray.o.x, oErr.x), oy(ray.o.y, oErr.y), oz(ray.o.z, oErr.z);
        EFloat dx(ray.d.x, dErr.x), dy(ray.d.y, dErr.y), dz(ray.d.z, dErr.z);
        EFloat a = dx * dx + dy * dy + dz * dz;
        EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
        EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius + maxDisplacement) * EFloat(radius + maxDisplacement);

        // Solve quadratic equation for _t_ values
        EFloat t0, t1;
        if (!Quadratic(a, b, c, &t0, &t1)) return false;

        EFloat tShapeHit = t0;
        if (tShapeHit.LowerBound() <= 0) 
        {
            tShapeHit = t1;
            if (tShapeHit.UpperBound() > ray.tMax) return false;
        }

        // find if ray intersects, and take distance along 
        // ray that intersection occurs
        bool intersect;
        *tHit = RayMarch(ray, &intersect);
        if (intersect) 
        {

            // Hit variables:
            Point3f  pHit = ray(*tHit);     // point of intersection
            Vector3f vHit = Vector3f(pHit); // vector respresentation of pHit
            Vector3f rHit = vHit - (Displacement(vHit) * Normalise(vHit)); // project hit into underlying sphere
            
            // Computing normals:
            Float phi = std::atan2(rHit.y, rHit.x);
            if (phi < 0) phi += 2 * Pi;

            Float theta = std::acos(Clamp(rHit.z / radius, -1, 1));

            Float u = phi / phiMax;
            Float v = (theta - thetaMin) / (thetaMax - thetaMin);

            Float du = u + esix; // small change in u
            Float dv = v + esix; // small change in v

            Float dPhi   = du * phiMax;
            Float dTheta = thetaMin + dv * (thetaMax - thetaMin);

            float r = radius;

            // find xyz of small change in u
            float dudx = r * sin(theta) * cos(dPhi);
            float dudy = r * sin(theta) * sin(dPhi);
            float dudz = r * cos(theta);
            Vector3f duxyz = Normalise(Vector3f(dudx, dudy, dudz));

            Float d = Displacement(duxyz);
            duxyz *= d + radius;

            // find xyz of small change in v
            float dvdx = r * sin(dTheta) * cos(phi);
            float dvdy = r * sin(dTheta) * sin(phi);
            float dvdz = r * cos(dTheta);
            Vector3f dvxyz = Normalise(Vector3f(dvdx, dvdy, dvdz));

            d = Displacement(dvxyz);
            dvxyz *= d + radius;

            Vector3f dpdu = (vHit - duxyz) / esix; // Finite difference normals
            Vector3f dpdv = (vHit - dvxyz) / esix;

            Vector3f N = Normalise(Cross(dpdu, dpdv));
            Normal3f dndu = Normal3f(N);
            Normal3f dndv = Normal3f(N);
            
            Vector3f pError = gamma(5) * Abs(vHit);
            if (Normalise(dpdu) == Normalise(dpdv)) dpdu.x *= 1.05;
            *isect = (*ObjectToWorld)(
                SurfaceInteraction(pHit, pError, Point2f(u, v), -ray.d,
                                   Normalise(dpdu), Normalise(dpdv), dndu, dndv, ray.time, this));
            
        }

        return intersect;

    }
    

    bool Dsphere::IntersectP(const Ray &r, bool testAlphaTexture) const 
    {
        Vector3f oErr, dErr;
        Ray ray = (*WorldToObject)(r, &oErr, &dErr);

        bool intersect;
        RayMarch(ray, &intersect);

        return intersect;
        
    }

    Bounds3f Dsphere::ObjectBound() const 
    {
        return Bounds3f(Point3f(-(radius + maxDisplacement), -(radius + maxDisplacement), zMin),
                        Point3f((radius + maxDisplacement), (radius + maxDisplacement), zMax));
    }

    ////////////////////////////////////
    // sphere functions:
    ////////////////////////////////////

    Float Dsphere::Area() const { return phiMax * radius * (zMax - zMin); }

    Interaction Dsphere::Sample(const Point2f &u, Float *pdf) const {
        Point3f pObj = Point3f(0, 0, 0) + radius * UniformSampleSphere(u);
        Interaction it;
        it.n = Normalize((*ObjectToWorld)(Normal3f(pObj.x, pObj.y, pObj.z)));
        if (reverseOrientation) it.n *= -1;
        // Reproject _pObj_ to sphere surface and compute _pObjError_
        pObj *= radius / Distance(pObj, Point3f(0, 0, 0));
        Vector3f pObjError = gamma(5) * Abs((Vector3f)pObj);
        it.p = (*ObjectToWorld)(pObj, pObjError, &it.pError);
        *pdf = 1 / Area();
        return it;
    }

    Interaction Dsphere::Sample(const Interaction &ref, const Point2f &u,
                               Float *pdf) const {
        Point3f pCenter = (*ObjectToWorld)(Point3f(0, 0, 0));

        // Sample uniformly on sphere if $\pt{}$ is inside it
        Point3f pOrigin =
            OffsetRayOrigin(ref.p, ref.pError, ref.n, pCenter - ref.p);
        if (DistanceSquared(pOrigin, pCenter) <= radius * radius) {
            Interaction intr = Sample(u, pdf);
            Vector3f wi = intr.p - ref.p;
            if (wi.LengthSquared() == 0)
                *pdf = 0;
            else {
                // Convert from area measure returned by Sample() call above to
                // solid angle measure.
                wi = Normalize(wi);
                *pdf *= DistanceSquared(ref.p, intr.p) / AbsDot(intr.n, -wi);
            }
            if (std::isinf(*pdf)) *pdf = 0.f;
            return intr;
        }

        // Sample sphere uniformly inside subtended cone

        // Compute coordinate system for sphere sampling
        Float dc = Distance(ref.p, pCenter);
        Float invDc = 1 / dc;
        Vector3f wc = (pCenter - ref.p) * invDc;
        Vector3f wcX, wcY;
        CoordinateSystem(wc, &wcX, &wcY);

        // Compute $\theta$ and $\phi$ values for sample in cone
        Float sinThetaMax = radius * invDc;
        Float sinThetaMax2 = sinThetaMax * sinThetaMax;
        Float invSinThetaMax = 1 / sinThetaMax;
        Float cosThetaMax = std::sqrt(std::max((Float)0.f, 1 - sinThetaMax2));

        Float cosTheta  = (cosThetaMax - 1) * u[0] + 1;
        Float sinTheta2 = 1 - cosTheta * cosTheta;

        if (sinThetaMax2 < 0.00068523f /* sin^2(1.5 deg) */) {
            /* Fall back to a Taylor series expansion for small angles, where
               the standard approach suffers from severe cancellation errors */
            sinTheta2 = sinThetaMax2 * u[0];
            cosTheta = std::sqrt(1 - sinTheta2);
        }

        // Compute angle $\alpha$ from center of sphere to sampled point on surface
        Float cosAlpha = sinTheta2 * invSinThetaMax +
            cosTheta * std::sqrt(std::max((Float)0.f, 1.f - sinTheta2 * invSinThetaMax * invSinThetaMax));
        Float sinAlpha = std::sqrt(std::max((Float)0.f, 1.f - cosAlpha*cosAlpha));
        Float phi = u[1] * 2 * Pi;

        // Compute surface normal and sampled point on sphere
        Vector3f nWorld =
            SphericalDirection(sinAlpha, cosAlpha, phi, -wcX, -wcY, -wc);
        Point3f pWorld = pCenter + radius * Point3f(nWorld.x, nWorld.y, nWorld.z);

        // Return _Interaction_ for sampled point on sphere
        Interaction it;
        it.p = pWorld;
        it.pError = gamma(5) * Abs((Vector3f)pWorld);
        it.n = Normal3f(nWorld);
        if (reverseOrientation) it.n *= -1;

        // Uniform cone PDF.
        *pdf = 1 / (2 * Pi * (1 - cosThetaMax));

        return it;
    }

    Float Dsphere::Pdf(const Interaction &ref, const Vector3f &wi) const {
        Point3f pCenter = (*ObjectToWorld)(Point3f(0, 0, 0));
        // Return uniform PDF if point is inside sphere
        Point3f pOrigin =
            OffsetRayOrigin(ref.p, ref.pError, ref.n, pCenter - ref.p);
        if (DistanceSquared(pOrigin, pCenter) <= radius * radius)
            return Shape::Pdf(ref, wi);

        // Compute general sphere PDF
        Float sinThetaMax2 = radius * radius / DistanceSquared(ref.p, pCenter);
        Float cosThetaMax = std::sqrt(std::max((Float)0, 1 - sinThetaMax2));
        return UniformConePdf(cosThetaMax);
    }

    Float Dsphere::SolidAngle(const Point3f &p, int nSamples) const {
        Point3f pCenter = (*ObjectToWorld)(Point3f(0, 0, 0));
        if (DistanceSquared(p, pCenter) <= radius * radius)
            return 4 * Pi;
        Float sinTheta2 = radius * radius / DistanceSquared(p, pCenter);
        Float cosTheta = std::sqrt(std::max((Float)0, 1 - sinTheta2));
        return (2 * Pi * (1 - cosTheta));
    }


    // Constructor:
    std::shared_ptr<Shape> CreateDsphereShape(const Transform *o2w,
                                             const Transform *w2o,
                                             bool reverseOrientation,
                                             const ParamSet &params) 
    {
        Float radius = params.FindOneFloat("radius", 1.f);
        Float phimax = params.FindOneFloat("phimax", 360.f);

        // Displacement Sphere specific parameters:
        Float maxDisplacement = radius * params.FindOneFloat("maxdispl", 1.f) / 100; // maximum percentage of radius map can apply 

        Float zmin = params.FindOneFloat("zmin", -radius - maxDisplacement);
        Float zmax = params.FindOneFloat("zmax", radius + maxDisplacement);

        std::string mapString = params.FindOneString("displacementmap", "");
        mapString.empty() ? analytical = true : analytical = false;

        if (!analytical) ReadPPM(mapString);
        else std::cout << "\nNo displacement map detected. Using analytical displacement map.\n";

        return std::make_shared<Dsphere>(o2w, w2o, reverseOrientation, radius, zmin,
                                        zmax, phimax, maxDisplacement);
    }

}  // namespace pbrt
