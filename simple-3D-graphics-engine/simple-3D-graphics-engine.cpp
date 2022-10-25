// simple-3D-graphics-engine.cpp : Este archivo contiene la función "main". La ejecución del programa comienza y termina ahí.
//
#include <iostream>
#include "olcConsoleGameEngine.h"
#include <fstream>
#include <strstream>
#include <algorithm>

using namespace std;

struct vec3d {
    float x = 0;
    float y = 0;
    float z = 0;
    float w = 1;
};

struct triangle {
    vec3d p[3];
    wchar_t sym;
    short col;
};

struct mesh {
    vector<triangle> tris;

    bool loadFromObjectFile(string filePath) {
        ifstream f(filePath);
        if (!f.is_open()) return false;

        vector<vec3d> verts;

        while (!f.eof()) {
            char line[128];
            f.getline(line, 128);

            strstream s;
            s << line;

            char junk;
            if (line[0] == 'v') {
                vec3d v;
                s >> junk >> v.x >> v.y >> v.z;
                verts.push_back(v);
            }

            if (line[0] == 'f') {
                int f[3];
                s >> junk >> f[0] >> f[1] >> f[2];
                tris.push_back({verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1]});
            }
        }
        return true;
    }
};

struct matrix4x4 {
    float m[4][4] = {0};
};

class olcEngine3D : public olcConsoleGameEngine {
public:
    olcEngine3D() {
        m_sAppName = L"3D demo";
    }

private:
    mesh meshObject;
    matrix4x4 matProj;

    vec3d vCamera;
    vec3d vLookDir;

    // Direction which the player is facing (rotation in the Y axis).
    float fYaw;

    float fTheta;

    vec3d matrixMultiplyVector(matrix4x4& m, vec3d& i) {
        vec3d v;

        v.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + i.w * m.m[3][0];
        v.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + i.w * m.m[3][1];
        v.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + i.w * m.m[3][2];
        v.w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + i.w * m.m[3][3];

        return v;
    }

    matrix4x4 matrixMakeIdentity() {
        matrix4x4 matrix;

        matrix.m[0][0] = 1.0f;
        matrix.m[1][1] = 1.0f;
        matrix.m[2][2] = 1.0f;
        matrix.m[3][3] = 1.0f;

        return matrix;
    }

    matrix4x4 matrixMakeRotationX(float fAngleRad) {
        matrix4x4 matrix;

        matrix.m[0][0] = 1.0f;
        matrix.m[1][1] = cosf(fAngleRad);
        matrix.m[1][2] = sinf(fAngleRad);
        matrix.m[2][1] = -sinf(fAngleRad);
        matrix.m[2][2] = cosf(fAngleRad);
        matrix.m[3][3] = 1.0f;

        return matrix;
    }

    matrix4x4 matrixMakeRotationY(float fAngleRad) {
        matrix4x4 matrix;

        matrix.m[0][0] = cosf(fAngleRad);
        matrix.m[0][2] = sinf(fAngleRad);
        matrix.m[2][0] = -sinf(fAngleRad);
        matrix.m[1][1] = 1.0f;
        matrix.m[2][2] = cosf(fAngleRad);
        matrix.m[3][3] = 1.0f;

        return matrix;
    }

    matrix4x4 matrixMakeRotationZ(float fAngleRad) {
        matrix4x4 matrix;

        matrix.m[0][0] = cosf(fAngleRad);
        matrix.m[0][1] = sinf(fAngleRad);
        matrix.m[1][0] = -sinf(fAngleRad);
        matrix.m[1][1] = cosf(fAngleRad);
        matrix.m[2][2] = 1.0f;
        matrix.m[3][3] = 1.0f;

        return matrix;
    }

    matrix4x4 matrixMakeTranslation(float x, float y, float z) {
        matrix4x4 matrix;

        matrix.m[0][0] = 1.0f;
        matrix.m[1][1] = 1.0f;
        matrix.m[2][2] = 1.0f;
        matrix.m[3][3] = 1.0f;
        matrix.m[3][0] = x;
        matrix.m[3][1] = y;
        matrix.m[3][2] = z;

        return matrix;
    }

    matrix4x4 matrixMakeProjection(float fFovDeg, float fAspectRatio, float fNear, float fFar) {
        float fFovRad = 1.0f / tanf(fFovDeg * 0.5f / 180.0f * 3.1416f);
        matrix4x4 matrix;

        matrix.m[0][0] = fAspectRatio * fFovRad;
        matrix.m[1][1] = fFovRad;
        matrix.m[2][2] = fFar / (fFar - fNear);
        matrix.m[3][2] = (-fFar * fNear) / (fFar - fNear);
        matrix.m[2][3] = 1.0f;
        matrix.m[3][3] = 0.0f;

        return matrix;
    }

    matrix4x4 matrixMultiplyMatrix(matrix4x4& m1, matrix4x4& m2) {
        matrix4x4 matrix;

        for (int c = 0; c < 4; c++) {
            for (int r = 0; r < 4; r++) {
                matrix.m[r][c] = m1.m[r][0] * m2.m[0][c] + m1.m[r][1] * m2.m[1][c] + m1.m[r][2] * m2.m[2][c] + m1.m[r][3] * m2.m[3][c];
            }
        }

        return matrix;
    }

    matrix4x4 matrixPointAt(vec3d& pos, vec3d& target, vec3d& up) {
        // Get new forward vector from pos to target.
        vec3d newForward = vectorSub(target, pos);
        newForward = vectorNormalise(newForward);

        // Get new up direction.
        vec3d a = vectorMul(newForward, vectorDotProduct(up, newForward));
        vec3d newUp = vectorSub(up, a);
        newUp = vectorNormalise(newUp);

        // Get new right direction vector, just doing Cross Product of up and forward.
        vec3d newRight = vectorCrossProduct(newUp, newForward);

        // Dimensioning and translation matrix.
        matrix4x4 matrix;
        matrix.m[0][0] = newRight.x;	matrix.m[0][1] = newRight.y;	matrix.m[0][2] = newRight.z;	matrix.m[0][3] = 0.0f;
        matrix.m[1][0] = newUp.x;		matrix.m[1][1] = newUp.y;		matrix.m[1][2] = newUp.z;		matrix.m[1][3] = 0.0f;
        matrix.m[2][0] = newForward.x;	matrix.m[2][1] = newForward.y;	matrix.m[2][2] = newForward.z;	matrix.m[2][3] = 0.0f;
        matrix.m[3][0] = pos.x;			matrix.m[3][1] = pos.y;			matrix.m[3][2] = pos.z;			matrix.m[3][3] = 1.0f;

        return matrix;
    }

    matrix4x4 matrixQuickInverse(matrix4x4& m) {
        // Invert matrix, only for Rotation and Translation matrices.
        matrix4x4 matrix;

        matrix.m[0][0] = m.m[0][0];     matrix.m[0][1] = m.m[1][0];     matrix.m[0][2] = m.m[2][0];     matrix.m[0][3] = 0.0f;
        matrix.m[1][0] = m.m[0][1];     matrix.m[1][1] = m.m[1][1];     matrix.m[1][2] = m.m[2][1];     matrix.m[1][3] = 0.0f;
        matrix.m[2][0] = m.m[0][2];     matrix.m[2][1] = m.m[1][2];     matrix.m[2][2] = m.m[2][2];     matrix.m[2][3] = 0.0f;

        matrix.m[3][0] = -(m.m[3][0] * matrix.m[0][0] + m.m[3][1] * matrix.m[1][0] + m.m[3][2] * matrix.m[2][0]);
        matrix.m[3][1] = -(m.m[3][0] * matrix.m[0][1] + m.m[3][1] * matrix.m[1][1] + m.m[3][2] * matrix.m[2][1]);
        matrix.m[3][2] = -(m.m[3][0] * matrix.m[0][2] + m.m[3][1] * matrix.m[1][2] + m.m[3][2] * matrix.m[2][2]);
        matrix.m[3][3] = 1.0f;

        return matrix;
    }

    vec3d vectorAdd(vec3d& v1, vec3d& v2) {
        return { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
    }

    vec3d vectorSub(vec3d& v1, vec3d& v2) {
        return { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
    }

    vec3d vectorMul(vec3d& v1, float k) {
        return { v1.x * k, v1.y * k, v1.z * k };
    }
    
    vec3d vectorDiv(vec3d& v1, float k) {
        return { v1.x / k, v1.y / k, v1.z / k };
    }

    float vectorDotProduct(vec3d& v1, vec3d& v2) {
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    }

    float vectorLength(vec3d& v) {
        return sqrtf(vectorDotProduct(v, v));
    }

    vec3d vectorNormalise(vec3d& v) {
        float len = vectorLength(v);
        return { v.x / len, v.y / len, v.z / len };
    }

    vec3d vectorCrossProduct(vec3d& v1, vec3d& v2) {
        vec3d v;

        v.x = v1.y * v2.z - v1.z * v2.y;
        v.y = v1.z * v2.x - v1.x * v2.z;
        v.z = v1.x * v2.y - v1.y * v2.x;
        
        return v;
    }

    vec3d vectorIntersectPlane(vec3d& plane_p, vec3d& plane_n, vec3d& lineStart, vec3d& lineEnd) {
        plane_n = vectorNormalise(plane_n);
        float plane_d = -vectorDotProduct(plane_n, plane_p);
        float ad = vectorDotProduct(lineStart, plane_n);
        float bd = vectorDotProduct(lineEnd, plane_n);
        float t = (-plane_d - ad) / (bd - ad);
        vec3d lineStartToEnd = vectorSub(lineEnd, lineStart);
        vec3d lineToIntersect = vectorMul(lineStartToEnd, t);

        return vectorAdd(lineStart, lineToIntersect);
    }

    int triangleClipAgainstPlane(vec3d plane_p, vec3d plane_n, triangle& in_tri, triangle& out_tri1, triangle& out_tri2) {
        plane_n = vectorNormalise(plane_n);

        auto dist = [&](vec3d& p) {
            vec3d n = vectorNormalise(p);
            return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z - vectorDotProduct(plane_n, plane_p));
        };

        vec3d* inside_points[3]; int nInsidePointCount = 0;
        vec3d* outside_points[3]; int nOutsidePointCount = 0;

        float d0 = dist(in_tri.p[0]);
        float d1 = dist(in_tri.p[1]);
        float d2 = dist(in_tri.p[2]);

        if (d0 >= 0) inside_points[nInsidePointCount++] = &in_tri.p[0];
        else outside_points[nOutsidePointCount++] = &in_tri.p[0];
        if (d1 >= 0) inside_points[nInsidePointCount++] = &in_tri.p[1];
        else outside_points[nOutsidePointCount++] = &in_tri.p[1];
        if (d2 >= 0) inside_points[nInsidePointCount++] = &in_tri.p[2];
        else outside_points[nOutsidePointCount++] = &in_tri.p[2];

        // Check how many points are outside and do the corresponding process.
        if (nInsidePointCount == 0) return 0;   // All the points are outside so is no valid returned triangle.

        if (nInsidePointCount == 3) {  // All the points are inside, so the triangle is already valid.
            out_tri1 = in_tri;
            return 1;
        }

        if (nInsidePointCount == 1 && nOutsidePointCount == 2) {    // Will return a smaller triangle.
            out_tri1.col = in_tri.col;
            out_tri1.sym = in_tri.sym;

            // Set valid point.
            out_tri1.p[0] = *inside_points[0];

            // The two new points are in the locations where the sides of the triangle intersect with the plane.
            out_tri1.p[1] = vectorIntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0]);
            out_tri1.p[2] = vectorIntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[1]);

            return 1;
        }

        if (nInsidePointCount == 2 && nOutsidePointCount == 1) {    // Will return two new triangles, maded by the quad former o¿by the two points inside and the points of intersection.
            out_tri1.col = in_tri.col;
            out_tri1.sym = in_tri.sym;

            out_tri2.col = in_tri.col;
            out_tri2.sym = in_tri.sym;

            // First triangle is composed of the two triangles inside and the intersection of one side of the triangle and the outside point.
            out_tri1.p[0] = *inside_points[0];
            out_tri1.p[1] = *inside_points[1];
            out_tri1.p[2] = vectorIntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0]);

            // Second triangle is composed of the two intersections of the two side of the triangle with the outside point and one inside point.
            out_tri2.p[0] = *inside_points[1];
            out_tri2.p[1] = out_tri1.p[2];  // We already have one intersection point.
            out_tri2.p[2] = vectorIntersectPlane(plane_p, plane_n, *inside_points[1], *outside_points[0]);

            return 2;
        }
    }

    CHAR_INFO getColour(float lum) {
        short bg_col, fg_col;
        wchar_t sym;
        int pixel_bw = (int)(13.0f * lum);
        switch (pixel_bw) {
            case 0: bg_col = BG_BLACK; fg_col = FG_BLACK; sym = PIXEL_SOLID; break;

            case 1: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_QUARTER; break;
            case 2: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_HALF; break;
            case 3: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_THREEQUARTERS; break;
            case 4: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_SOLID; break;

            case 5: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_QUARTER; break;
            case 6: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_HALF; break;
            case 7: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_THREEQUARTERS; break;
            case 8: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_SOLID; break;

            case 9:  bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_QUARTER; break;
            case 10: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_HALF; break;
            case 11: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_THREEQUARTERS; break;
            case 12: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_SOLID; break;
            default: bg_col = BG_BLACK; fg_col = FG_BLACK; sym = PIXEL_SOLID;
        }

        CHAR_INFO c;
        c.Attributes = bg_col | fg_col;
        c.Char.UnicodeChar = sym;
        return c;
    }

public:
    bool OnUserCreate() override {
        meshObject.loadFromObjectFile("objects_src/axis.obj");

        matProj = matrixMakeProjection(90.0f, (float)ScreenHeight() / (float)ScreenWidth(), 0.1f, 1000.0f);

        return true;
    }

    bool OnUserUpdate(float fElapsedTime) override {
        
        // Move forward and backwards.
        vec3d vForward = vectorMul(vLookDir, 8.0f * fElapsedTime);

        if (GetKey(L'W').bHeld) vCamera = vectorAdd(vCamera, vForward);
        if (GetKey(L'S').bHeld) vCamera = vectorSub(vCamera, vForward);

        if (GetKey(L'A').bHeld) vCamera.x -= 8.0f * fElapsedTime;
        if (GetKey(L'D').bHeld) vCamera.x += 8.0f * fElapsedTime;

        // Rotate and go up and down.
        if (GetKey(VK_UP).bHeld) vCamera.y -= 8.0f * fElapsedTime;
        if (GetKey(VK_DOWN).bHeld) vCamera.y += 8.0f * fElapsedTime;
        if (GetKey(VK_LEFT).bHeld) fYaw += 2.0f * fElapsedTime;
        if (GetKey(VK_RIGHT).bHeld) fYaw -= 2.0f * fElapsedTime;


        Fill(0, 0, ScreenWidth(), ScreenHeight(), PIXEL_SOLID, BG_BLACK);

        matrix4x4 matRotZ, matRotX;
        // fTheta += 1.0f * fElapsedTime;

        matRotZ = matrixMakeRotationZ(fTheta * 0.5f);
        matRotX = matrixMakeRotationX(fTheta);

        matrix4x4 matTrans;
        matTrans = matrixMakeTranslation(0.0f, 0.0f, 5.0f);

        matrix4x4 matWorld;
        matWorld = matrixMakeIdentity();
        matWorld = matrixMultiplyMatrix(matRotZ, matRotX);
        matWorld = matrixMultiplyMatrix(matWorld, matTrans);

        // Transformation for rotation and tranlation matrices (Camera).
        vec3d vUp = { 0, 1, 0 };
        vec3d vTarget = { 0, 0, 1 };
        matrix4x4 matCameraRot = matrixMakeRotationY(fYaw);
        vLookDir = matrixMultiplyVector(matCameraRot, vTarget);
        vTarget = vectorAdd(vCamera, vLookDir);

        matrix4x4 matCamera = matrixPointAt(vCamera, vTarget, vUp);

        // Make view matrix from camera.
        matrix4x4 matView = matrixQuickInverse(matCamera);

        // Store triangles
        vector<triangle> vecTrianglesToRaster;

        // Draw triangles
        triangle triProjected, triTransformed, triViewed;
        for (auto &tri : meshObject.tris) {

            triTransformed.p[0] = matrixMultiplyVector(matWorld, tri.p[0]);
            triTransformed.p[1] = matrixMultiplyVector(matWorld, tri.p[1]);
            triTransformed.p[2] = matrixMultiplyVector(matWorld, tri.p[2]);

            // Identify which faces are facing the camera
            // We will calculate the Cross Product to obtain the Normal vector, 
            // and then calculate the Dot Product between the normal vector and a camera vector to the triangle,
            // and if it is smaller than 0, that means that is not facing the camera.
            
            // Calculate triangle Normal.
            vec3d normal, line1, line2;

            // Get lines either side of triangle.
            line1 = vectorSub(triTransformed.p[1], triTransformed.p[0]);
            line2 = vectorSub(triTransformed.p[2], triTransformed.p[0]);

            // Make Cross Product of the triangle to get Normal vector.
            normal = vectorCrossProduct(line1, line2);

            // Normalise the normal vector.
            normal = vectorNormalise(normal);

            // Get Ray from the triangle to the camera.
            vec3d vCameraRay = vectorSub(triTransformed.p[0], vCamera);

            if (vectorDotProduct(normal, vCameraRay) < 0.0f) {
                // If the normal vector aligns better with the light vector, more shiny it will be.
                vec3d light_direction = { 0.0f, 1.0f, -1.0f };
                light_direction = vectorNormalise(light_direction);

                // Check how similar are the light_direction vector, and the normal vector.
                float dp = max(0.1f, vectorDotProduct(light_direction, normal));

                // Get console colors.
                CHAR_INFO c = getColour(dp);
                triTransformed.col = c.Attributes;
                triTransformed.sym = c.Char.UnicodeChar;

                // Convert world space to view space.
                triViewed.p[0] = matrixMultiplyVector(matView, triTransformed.p[0]);
                triViewed.p[1] = matrixMultiplyVector(matView, triTransformed.p[1]);
                triViewed.p[2] = matrixMultiplyVector(matView, triTransformed.p[2]);
                triViewed.sym = triTransformed.sym;
                triViewed.col = triTransformed.col;

                // Clip the triangles.
                int nClippedTriangles = 0;
                triangle clipped[2];    // Size = 2 because we can have maximum 2 new triangles.
                nClippedTriangles = triangleClipAgainstPlane({ 0.0f, 0.0f, 0.1f }, { 0.0f, 0.0f, 1.0f }, triViewed, clipped[0], clipped[1]);

                for (int n = 0; n < nClippedTriangles; n++) {

                    // Project triangles into 2d
                    triProjected.p[0] = matrixMultiplyVector(matProj, clipped[n].p[0]);
                    triProjected.p[1] = matrixMultiplyVector(matProj, clipped[n].p[1]);
                    triProjected.p[2] = matrixMultiplyVector(matProj, clipped[n].p[2]);
                    triProjected.col = clipped[n].col;
                    triProjected.sym = clipped[n].sym;

                    // Normalise vectors.
                    triProjected.p[0] = vectorDiv(triProjected.p[0], triProjected.p[0].w);
                    triProjected.p[1] = vectorDiv(triProjected.p[1], triProjected.p[1].w);
                    triProjected.p[2] = vectorDiv(triProjected.p[2], triProjected.p[2].w);

                    // Offset verts into visible normalised space.
                    vec3d vOffsetView = { 1, 1, 0 };
                    triProjected.p[0] = vectorAdd(triProjected.p[0], vOffsetView);
                    triProjected.p[1] = vectorAdd(triProjected.p[1], vOffsetView);
                    triProjected.p[2] = vectorAdd(triProjected.p[2], vOffsetView);

                    triProjected.p[0].x *= 0.5f * (float)ScreenWidth();
                    triProjected.p[0].y *= 0.5f * (float)ScreenHeight();
                    triProjected.p[1].x *= 0.5f * (float)ScreenWidth();
                    triProjected.p[1].y *= 0.5f * (float)ScreenHeight();
                    triProjected.p[2].x *= 0.5f * (float)ScreenWidth();
                    triProjected.p[2].y *= 0.5f * (float)ScreenHeight();

                    // Store triangles for sorting
                    vecTrianglesToRaster.push_back(triProjected);
                }
            }
        }

        // Sort triangles from back to front.
        sort(vecTrianglesToRaster.begin(), vecTrianglesToRaster.end(), 
            [](triangle& t1, triangle& t2) {
                float z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0f;
                float z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0f;

                return z1 > z2;
            });

        for (auto &triToRaster : vecTrianglesToRaster) {
            // Clip triangles against screen edges.
            triangle clipped[2];
            list<triangle> listTriangles;

            listTriangles.push_back(triToRaster);
            int nNewTriangles = 1;

            for (int p = 0; p < 4; p++) {

                int nTrisToAdd = 0;
                while (nNewTriangles > 0) {
                    triangle test = listTriangles.front();
                    listTriangles.pop_front();
                    nNewTriangles--;

                    // Check the four screen edges.
                    switch (p) {
                    case 0: nTrisToAdd = triangleClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
                    case 1: nTrisToAdd = triangleClipAgainstPlane({ 0.0f, (float)ScreenHeight() - 1, 0.0f }, { 0.0f, -1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
                    case 2: nTrisToAdd = triangleClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
                    case 3: nTrisToAdd = triangleClipAgainstPlane({ (float)ScreenWidth() - 1, 0.0f, 0.0f }, { -1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
                    }

                    for (int w = 0; w < nTrisToAdd; w++) 
                        listTriangles.push_back(clipped[w]);
                    
                }
                nNewTriangles = listTriangles.size();

                for (auto& t : listTriangles) {
                    FillTriangle(
                        t.p[0].x, t.p[0].y,
                        t.p[1].x, t.p[1].y,
                        t.p[2].x, t.p[2].y,
                        t.sym, t.col
                    );

                    /*DrawTriangle(
                        t.p[0].x, t.p[0].y,
                        t.p[1].x, t.p[1].y,
                        t.p[2].x, t.p[2].y,
                        PIXEL_SOLID, FG_BLACK
                    );*/
                }
            }
        }

        return true;
    }
};

int main()
{
    olcEngine3D demo;
    if (demo.ConstructConsole(256, 240, 2, 2)) demo.Start();
    else cout << "Couldn't start the app." << "\n";
}

// Ejecutar programa: Ctrl + F5 o menú Depurar > Iniciar sin depurar
// Depurar programa: F5 o menú Depurar > Iniciar depuración

// Sugerencias para primeros pasos: 1. Use la ventana del Explorador de soluciones para agregar y administrar archivos
//   2. Use la ventana de Team Explorer para conectar con el control de código fuente
//   3. Use la ventana de salida para ver la salida de compilación y otros mensajes
//   4. Use la ventana Lista de errores para ver los errores
//   5. Vaya a Proyecto > Agregar nuevo elemento para crear nuevos archivos de código, o a Proyecto > Agregar elemento existente para agregar archivos de código existentes al proyecto
//   6. En el futuro, para volver a abrir este proyecto, vaya a Archivo > Abrir > Proyecto y seleccione el archivo .sln
