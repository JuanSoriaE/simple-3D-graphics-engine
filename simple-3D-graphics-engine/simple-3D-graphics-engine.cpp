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
    mesh meshCube;
    matrix4x4 matProj;
    vec3d vCamera;
    float fTheta = 0.0f;

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
        meshCube.loadFromObjectFile("objects_src/ship.obj");

        matProj = matrixMakeProjection(90.0f, (float)ScreenHeight() / (float)ScreenWidth(), 0.1f, 1000.0f);

        return true;
    }

    bool OnUserUpdate(float fElapsedTime) override {
        Fill(0, 0, ScreenWidth(), ScreenHeight(), PIXEL_SOLID, BG_BLACK);

        matrix4x4 matRotZ, matRotX;
        fTheta += 1.0f * fElapsedTime;

        matRotZ = matrixMakeRotationZ(fTheta * 0.5f);
        matRotX = matrixMakeRotationX(fTheta);

        matrix4x4 matTrans;
        matTrans = matrixMakeTranslation(0.0f, 0.0f, 16.0f);

        matrix4x4 matWorld;
        matWorld = matrixMakeIdentity();
        matWorld = matrixMultiplyMatrix(matRotZ, matRotX);
        matWorld = matrixMultiplyMatrix(matWorld, matTrans);

        // Store triangles
        vector<triangle> vecTrianglesToRaster;

        // Draw triangles
        triangle triProjected, triTransformed;
        for (auto &tri : meshCube.tris) {

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

                // Project triangles into 2d
                triProjected.p[0] = matrixMultiplyVector(matProj, triTransformed.p[0]);
                triProjected.p[1] = matrixMultiplyVector(matProj, triTransformed.p[1]);
                triProjected.p[2] = matrixMultiplyVector(matProj, triTransformed.p[2]);
                triProjected.col = triTransformed.col;
                triProjected.sym = triTransformed.sym;

                // Normalise vectors.
                triProjected.p[0] = vectorDiv(triProjected.p[0], triProjected.p[0].w);
                triProjected.p[1] = vectorDiv(triProjected.p[1], triProjected.p[1].w);
                triProjected.p[2] = vectorDiv(triProjected.p[2], triProjected.p[2].w);

                // Offset verts into visible normalised space.
                vec3d vOffsetView = {1, 1, 0};
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

        // Sort triangles from back to front.
        sort(vecTrianglesToRaster.begin(), vecTrianglesToRaster.end(), 
            [](triangle& t1, triangle& t2) {
                float z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0f;
                float z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0f;

                return z1 > z2;
            });

        for (auto &triProjected : vecTrianglesToRaster) {
            // Rasterize triangles
            FillTriangle(
                triProjected.p[0].x, triProjected.p[0].y,
                triProjected.p[1].x, triProjected.p[1].y,
                triProjected.p[2].x, triProjected.p[2].y,
                triProjected.sym, triProjected.col
            );
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
