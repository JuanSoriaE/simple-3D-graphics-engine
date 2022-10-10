// simple-3D-graphics-engine.cpp : Este archivo contiene la función "main". La ejecución del programa comienza y termina ahí.
//
#include <iostream>
#include "olcConsoleGameEngine.h"
#include <fstream>
#include <strstream>
#include <algorithm>

using namespace std;

struct vec3d {
    float x, y, z;
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
    float f_theta = 0.0f;
    mesh meshCube, meshTest;
    matrix4x4 matProj;
    vec3d vCamera;

    void multiplyMatrixVector(vec3d &i, vec3d &o, matrix4x4 &m) {
        o.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + m.m[3][0];
        o.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + m.m[3][1];
        o.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + m.m[3][2];
        float w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + m.m[3][3];

        if (w != 0.0f) {
            o.x /= w;
            o.y /= w;
            o.z /= w;
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
        /*
        meshTest.tris = {
            {0.0f, 0.0f, 0.0f,  1.0f, 1.0f, 0.0f,   1.0f, 0.0f, 0.0f},

            // Square
            {1.0f, 0.0f, 0.0f,  1.0f, 1.0f, 0.0f,   1.0f, 1.0f, 1.0f},
            {1.0f, 0.0f, 0.0f,  1.0f, 1.0f, 1.0f,   1.0f, 0.0f, 1.0f},

            {1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f,   0.0f, 0.0f, 1.0f},

            {0.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f,   1.0f, 1.0f, 0.0f},
            {0.0f, 0.0f, 1.0f,  1.0f, 1.0f, 0.0f,   0.0f, 0.0f, 0.0f},

            {0.0f, 0.0f, 0.0f,  0.0f, 0.0f, 1.0f,   1.0f, 0.0f, 1.0f},
            {0.0f, 0.0f, 0.0f,  1.0f, 0.0f, 1.0f,   1.0f, 0.0f, 0.0f},

        };

        meshCube.tris = {
            // SOUTH
            {0.0f, 0.0f, 0.0f,  0.0f, 1.0f, 0.0f,   1.0f, 1.0f, 0.0f},
            {0.0f, 0.0f, 0.0f,  1.0f, 1.0f, 0.0f,   1.0f, 0.0f, 0.0f},

            // EAST
            {1.0f, 0.0f, 0.0f,  1.0f, 1.0f, 0.0f,   1.0f, 1.0f, 1.0f},
            {1.0f, 0.0f, 0.0f,  1.0f, 1.0f, 1.0f,   1.0f, 0.0f, 1.0f},

            // NORTH
            {1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f,   0.0f, 1.0f, 1.0f},
            {1.0f, 0.0f, 1.0f,  0.0f, 1.0f, 1.0f,   0.0f, 0.0f, 1.0f},

            // WEST
            {0.0f, 0.0f, 1.0f,  0.0f, 1.0f, 1.0f,   0.0f, 1.0f, 0.0f},
            {0.0f, 0.0f, 1.0f,  0.0f, 1.0f, 0.0f,   0.0f, 0.0f, 0.0f},

            // TOP
            {0.0f, 1.0f, 0.0f,  0.0f, 1.0f, 1.0f,   1.0f, 1.0f, 1.0f},
            {0.0f, 1.0f, 0.0f,  1.0f, 1.0f, 1.0f,   1.0f, 1.0f, 0.0f},

            // BOTTOM
            {1.0f, 0.0f, 1.0f,  0.0f, 0.0f, 1.0f,   0.0f, 0.0f, 0.0f},
            {1.0f, 0.0f, 1.0f,  0.0f, 0.0f, 0.0f,   1.0f, 0.0f, 0.0f
        };
        */

        meshCube.loadFromObjectFile("objects_src/ship.obj");

        // Projection Matrix
        float f_near = 0.1f;
        float f_far = 1000.0f;
        float f_fov = 90.0f;
        float f_aspect_ratio = (float)ScreenHeight() / (float)ScreenWidth();
        float f_fov_rad = 1.0f / tanf(f_fov * 0.5f / 180.0f * 3.1416f);

        matProj.m[0][0] = f_aspect_ratio * f_fov_rad;
        matProj.m[1][1] = f_fov_rad;
        matProj.m[2][2] = f_far / (f_far - f_near);
        matProj.m[3][2] = (-f_far * f_near) / (f_far - f_near);
        matProj.m[2][3] = 1.0f;
        matProj.m[3][3] = 0.0f;

        return true;
    }

    bool OnUserUpdate(float fElapsedTime) override {
        Fill(0, 0, ScreenWidth(), ScreenHeight(), PIXEL_SOLID, BG_BLACK);

        matrix4x4 matRotZ, matRotX;
        f_theta += 1.0f * fElapsedTime;

        // Z Rotation
        matRotZ.m[0][0] = cosf(f_theta);
        matRotZ.m[0][1] = sinf(f_theta);
        matRotZ.m[1][0] = -sinf(f_theta);
        matRotZ.m[1][1] = cosf(f_theta);
        matRotZ.m[2][2] = 1;
        matRotZ.m[3][3] = 1;

        // X Rotation
        matRotX.m[0][0] = 1;
        matRotX.m[1][1] = cosf(f_theta * 0.5f);
        matRotX.m[1][2] = sinf(f_theta * 0.5f);
        matRotX.m[2][1] = -sinf(f_theta * 0.5f);
        matRotX.m[2][2] = cosf(f_theta * 0.5f);
        matRotX.m[3][3] = 1;

        vector<triangle> vecTrianglesToRaster;

        // Draw triangles
        triangle triProjected, triTranslated, triRotatedZ, triRotatedZX;
        for (auto &tri : meshCube.tris) {

            multiplyMatrixVector(tri.p[0], triRotatedZ.p[0], matRotZ);
            multiplyMatrixVector(tri.p[1], triRotatedZ.p[1], matRotZ);
            multiplyMatrixVector(tri.p[2], triRotatedZ.p[2], matRotZ);

            // Rotate in X axis
            multiplyMatrixVector(triRotatedZ.p[0], triRotatedZX.p[0], matRotX);
            multiplyMatrixVector(triRotatedZ.p[1], triRotatedZX.p[1], matRotX);
            multiplyMatrixVector(triRotatedZ.p[2], triRotatedZX.p[2], matRotX);

            // Screen offset
            triTranslated = triRotatedZX;
            triTranslated.p[0].z = triRotatedZX.p[0].z + 7.0f;
            triTranslated.p[1].z = triRotatedZX.p[1].z + 7.0f;
            triTranslated.p[2].z = triRotatedZX.p[2].z + 7.0f;

            // Identify which faces are facing the camera
            // We will calculate the Cross Product to obtain the Normal vector, 
            // and then calculate the Dot Product between the normal vector and a camera vector to the triangle,
            // and if it is smaller than 0, that means that is not facing the camera.
            vec3d normal, line1, line2;
            line1.x = triTranslated.p[1].x - triTranslated.p[0].x;
            line1.y = triTranslated.p[1].y - triTranslated.p[0].y;
            line1.z = triTranslated.p[1].z - triTranslated.p[0].z;

            line2.x = triTranslated.p[2].x - triTranslated.p[0].x;
            line2.y = triTranslated.p[2].y - triTranslated.p[0].y;
            line2.z = triTranslated.p[2].z - triTranslated.p[0].z;

            normal.x = line1.y * line2.z - line1.z * line2.y;
            normal.y = line1.z * line2.x - line1.x * line2.z;
            normal.z = line1.x * line2.y - line1.y * line2.x;

            float normalLen = sqrtf(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
            normal.x /= normalLen;
            normal.y /= normalLen;
            normal.z /= normalLen;

            if (normal.x * (triTranslated.p[0].x - vCamera.x) +
                normal.y * (triTranslated.p[0].y - vCamera.y) +
                normal.z * (triTranslated.p[0].z - vCamera.z) < 0.0) {
                // Illumination test
                // If the normal vector aligns better with the light vector, more shiny it will be.
                vec3d lightDirection = {0.0f, 0.0f, -1.0f};
                float lightLen = sqrtf(lightDirection.x* lightDirection.x + lightDirection.y* lightDirection.y + lightDirection.z* lightDirection.z);
                lightDirection.x /= lightLen;
                lightDirection.y /= lightLen;
                lightDirection.z /= lightLen;

                float dp = normal.x * lightDirection.x + normal.y * lightDirection.y + normal.z * lightDirection.z;

                CHAR_INFO c = getColour(dp);
                triTranslated.col = c.Attributes;
                triTranslated.sym = c.Char.UnicodeChar;

                // Project triangles into 2d
                multiplyMatrixVector(triTranslated.p[0], triProjected.p[0], matProj);
                multiplyMatrixVector(triTranslated.p[1], triProjected.p[1], matProj);
                multiplyMatrixVector(triTranslated.p[2], triProjected.p[2], matProj);
                triProjected.col = triTranslated.col;
                triProjected.sym = triTranslated.sym;

                // Scale
                triProjected.p[0].x += 1.0f; triProjected.p[0].y += 1.0f;
                triProjected.p[1].x += 1.0f; triProjected.p[1].y += 1.0f;
                triProjected.p[2].x += 1.0f; triProjected.p[2].y += 1.0f;

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
