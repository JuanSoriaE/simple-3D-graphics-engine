// simple-3D-graphics-engine.cpp : Este archivo contiene la función "main". La ejecución del programa comienza y termina ahí.
//
#include <iostream>
#include "olcConsoleGameEngine.h"

using namespace std;

struct vec3d {
    float x, y, z;
};

struct triangle {
    vec3d p[3];
};

struct mesh {
    vector<triangle> tris;
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
    mesh meshCube;
    matrix4x4 matProj;

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

public:
    bool OnUserCreate() override {

        meshCube.tris = {
            // SOUTH
            {0.0f, 0.0f, 0.0f,  0.0f, 1.0f, 0.0f,   1.0f, 1.0f, 0.0f},
            {0.0f, 0.0f, 0.0f,  1.0f, 1.0f, 0.0f,   1.0f, 1.0f, 0.0f},

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
            {1.0f, 0.0f, 1.0f,  0.0f, 0.0f, 0.0f,   1.0f, 0.0f, 0.0f},

        };

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

        // Draw triangles
        triangle triProjected, triTranslated, triRotatedZ, triRotatedZX;
        for (auto &tri : meshCube.tris) {

            multiplyMatrixVector(tri.p[0], triRotatedZ.p[0], matRotZ);
            multiplyMatrixVector(tri.p[1], triRotatedZ.p[1], matRotZ);
            multiplyMatrixVector(tri.p[2], triRotatedZ.p[2], matRotZ);

            multiplyMatrixVector(triRotatedZ.p[0], triRotatedZX.p[0], matRotX);
            multiplyMatrixVector(triRotatedZ.p[1], triRotatedZX.p[1], matRotX);
            multiplyMatrixVector(triRotatedZ.p[2], triRotatedZX.p[2], matRotX);

            triTranslated = triRotatedZX;
            triTranslated.p[0].z = triRotatedZX.p[0].z + 3.0f;
            triTranslated.p[1].z = triRotatedZX.p[1].z + 3.0f;
            triTranslated.p[2].z = triRotatedZX.p[2].z + 3.0f;

            multiplyMatrixVector(triTranslated.p[0], triProjected.p[0], matProj);
            multiplyMatrixVector(triTranslated.p[1], triProjected.p[1], matProj);
            multiplyMatrixVector(triTranslated.p[2], triProjected.p[2], matProj);

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

            DrawTriangle(
                triProjected.p[0].x, triProjected.p[0].y, 
                triProjected.p[1].x, triProjected.p[1].y, 
                triProjected.p[2].x, triProjected.p[2].y,
                PIXEL_SOLID, FG_WHITE
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
