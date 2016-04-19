#include <iostream>
#include "ParticleSystem.h"
#include <fstream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include "shader.hpp"
#include "controls.hpp"

GLFWwindow* window;

using namespace std;
using namespace glm;

ParticleSystem *PBF = new ParticleSystem();

void init_liquid() {
    for (double x=5;x<15;x+=0.5) {
        for (double y=5;y<15;y+=1) {
            for (double z=5;z<15;z+=1) {
                PBF->add_particle(x,y,z);
            }
        }
    }
}

std::vector <GLfloat> init_color(size_t length) {
    std::vector <GLfloat> colors;
    for (int i = 0; i < length; i++) {
        colors.push_back(glm::linearRand(0.0f,1.0f));
        colors.push_back(glm::linearRand(0.0f,1.0f));
        colors.push_back(glm::linearRand(0.0f,1.0f));
    }
    return colors;
}

std::vector <GLfloat> get_positions() {
    std::vector <GLfloat> positions;
    for (auto i : PBF->particles) {
        positions.push_back(float(i->x.x));
        positions.push_back(float(i->x.z));
        positions.push_back(float(i->x.y));
    }
    return positions;
}

std::vector <GLfloat> marching_cubes() {
    std::vector <GLfloat> points;
    for (size_t k = 0; k < PBF->kmax; k++) {
        for (size_t j = 0; j < PBF->jmax; j++) {
            for (size_t i = 0; i < PBF->imax; i++) {
                std::vector<double> scalars;
                std::vector<glm::dvec3> grid;
                scalars.push_back(PBF->scalar_field[k * (PBF->jmax * PBF->imax) + j * (PBF->imax) + i]);
                grid.push_back(glm::dvec3(i * PBF->h, j * PBF->h, k * PBF->h));
                scalars.push_back(PBF->scalar_field[k * (PBF->jmax * PBF->imax) + j * (PBF->imax) + (i + 1)]);
                grid.push_back(glm::dvec3((i + 1) * PBF->h, j * PBF->h, k * PBF->h));
                scalars.push_back(PBF->scalar_field[k * (PBF->jmax * PBF->imax) + (j + 1) * (PBF->imax) + (i + 1)]);
                grid.push_back(glm::dvec3((i + 1) * PBF->h, (j + 1) * PBF->h, k * PBF->h));
                scalars.push_back(PBF->scalar_field[k * (PBF->jmax * PBF->imax) + (j + 1) * (PBF->imax) + i]);
                grid.push_back(glm::dvec3(i * PBF->h, (j + 1) * PBF->h, k * PBF->h));
                scalars.push_back(PBF->scalar_field[(k + 1) * (PBF->jmax * PBF->imax) + j * (PBF->imax) + i]);
                grid.push_back(glm::dvec3(i * PBF->h, j * PBF->h, (k + 1) * PBF->h));
                scalars.push_back(PBF->scalar_field[(k + 1) * (PBF->jmax * PBF->imax) + j * (PBF->imax) + (i + 1)]);
                grid.push_back(glm::dvec3((i + 1) * PBF->h, j * PBF->h, (k + 1) * PBF->h));
                scalars.push_back(PBF->scalar_field[(k + 1) * (PBF->jmax * PBF->imax) + (j + 1) * (PBF->imax) + (i + 1)]);
                grid.push_back(glm::dvec3((i + 1) * PBF->h, (j + 1) * PBF->h, (k + 1) * PBF->h));
                scalars.push_back(PBF->scalar_field[(k + 1) * (PBF->jmax * PBF->imax) + (j + 1) * (PBF->imax) + i]);
                grid.push_back(glm::dvec3(i * PBF->h, (j + 1) * PBF->h, (k + 1) * PBF->h));
                std::vector<glm::dvec3> vertices = PBF->polygonise(grid,scalars,PBF->isolevel);
                for (auto v : vertices) {
                    points.push_back(float(v.x));
                    points.push_back(float(v.z));
                    points.push_back(float(v.y));
                }
            }
        }
    }
    return points;
}

int main() {
    // Initialise GLFW
    if( !glfwInit() )
    {
        fprintf( stderr, "Failed to initialize GLFW\n" );
        getchar();
        return -1;
    }

    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    // Open a window and create its OpenGL context
    window = glfwCreateWindow( 1024, 768, "Tutorial 04 - Colored Cube", NULL, NULL);
    if( window == NULL ){
        fprintf( stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n" );
        getchar();
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

    // Initialize GLEW
    glewExperimental = true; // Needed for core profile
    if (glewInit() != GLEW_OK) {
        fprintf(stderr, "Failed to initialize GLEW\n");
        getchar();
        glfwTerminate();
        return -1;
    }

    // Ensure we can capture the escape key being pressed below
    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    glfwPollEvents();
    glfwSetCursorPos(window, 1024/2, 768/2);

    // Dark blue background
    glClearColor(0.0f, 0.0f, 0.4f, 0.0f);

    // Enable depth test
    glEnable(GL_DEPTH_TEST);
    // Accept fragment if it closer to the camera than the former one
    glDepthFunc(GL_LESS);

    glEnable(GL_PROGRAM_POINT_SIZE);

    GLuint VertexArrayID;
    glGenVertexArrays(1, &VertexArrayID);
    glBindVertexArray(VertexArrayID);

    // Create and compile our GLSL program from the shaders
    GLuint programID = LoadShaders( "TransformVertexShader.vertexshader", "ColorFragmentShader.fragmentshader" );

    // Get a handle for our "MVP" uniform
    GLuint MatrixID = glGetUniformLocation(programID, "MVP");

    /*
    // Projection matrix : 45° Field of View, 4:3 ratio, display range : 0.1 unit <-> 100 units
    glm::mat4 Projection = glm::perspective(glm::radians(45.0f), 4.0f / 3.0f, 0.1f, 100.0f);
    // Camera matrix
    glm::mat4 View       = glm::lookAt(
            glm::vec3(50,50,52), // Camera is at (4,3,-3), in World Space
            glm::vec3(0,0,0), // and looks at the origin
            glm::vec3(0,1,0)  // Head is up (set to 0,-1,0 to look upside-down)
    );
    // Model matrix : an identity matrix (model will be at the origin)
    glm::mat4 Model      = glm::mat4(1.0f);
    // Our ModelViewProjection : multiplication of our 3 matrices
    glm::mat4 MVP        = Projection * View * Model; // Remember, matrix multiplication is the other way around
    */
    // Our vertices. Three consecutive floats give a 3D vertex; Three consecutive vertices give a triangle.
    // A cube has 6 faces with 2 triangles each, so this makes 6*2=12 triangles, and 12*3 vertices

    init_liquid();

    PBF->step();

    std::vector <GLfloat> vertices = marching_cubes();

    std::vector <GLfloat> colors = init_color(PBF->particles.size());

    GLuint vertexbuffer;
    glGenBuffers(1, &vertexbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    glBufferData(GL_ARRAY_BUFFER, vertices.size()*sizeof(GLfloat), &vertices.front(), GL_STATIC_DRAW);

    GLuint colorbuffer;
    glGenBuffers(1, &colorbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, colorbuffer);
    glBufferData(GL_ARRAY_BUFFER, colors.size()* sizeof(GLfloat), &colors.front(), GL_STATIC_DRAW);
    glFrontFace(GL_CW);
    do{

        vertices = marching_cubes();
        glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
        glBufferData(GL_ARRAY_BUFFER, vertices.size()*sizeof(GLfloat), &vertices.front(), GL_STATIC_DRAW);

        // Clear the screen
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Use our shader
        glUseProgram(programID);

        computeMatricesFromInputs();
        glm::mat4 ProjectionMatrix = getProjectionMatrix();
        glm::mat4 ViewMatrix = getViewMatrix();
        glm::mat4 ModelMatrix = glm::mat4(1.0);
        glm::mat4 MVP = ProjectionMatrix * ViewMatrix * ModelMatrix;

        // Send our transformation to the currently bound shader,
        // in the "MVP" uniform
        glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);

        // 1rst attribute buffer : vertices
        glEnableVertexAttribArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
        glVertexAttribPointer(
                0,                  // attribute. No particular reason for 0, but must match the layout in the shader.
                3,                  // size
                GL_FLOAT,           // type
                GL_FALSE,           // normalized?
                0,                  // stride
                (void*)0            // array buffer offset
        );

        // 2nd attribute buffer : colors
        glEnableVertexAttribArray(1);
        glBindBuffer(GL_ARRAY_BUFFER, colorbuffer);
        glVertexAttribPointer(
                1,                                // attribute. No particular reason for 1, but must match the layout in the shader.
                3,                                // size
                GL_FLOAT,                         // type
                GL_FALSE,                         // normalized?
                0,                                // stride
                (void*)0                          // array buffer offset
        );

        // Draw the triangle !
        glDrawArrays(GL_TRIANGLES, 0, vertices.size()/3); // 12*3 indices starting at 0 -> 12 triangles

        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);

        PBF->step();

        // Swap buffers
        glfwSwapBuffers(window);
        glfwPollEvents();

    } // Check if the ESC key was pressed or the window was closed
    while( glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
           glfwWindowShouldClose(window) == 0 );

    // Cleanup VBO and shader
    glDeleteBuffers(1, &vertexbuffer);
    glDeleteBuffers(1, &colorbuffer);
    glDeleteProgram(programID);
    glDeleteVertexArrays(1, &VertexArrayID);

    // Close OpenGL window and terminate GLFW
    glfwTerminate();

    return 0;

/*
    ParticleSystem *system = new ParticleSystem();
    //ofstream file;
    for (double x=5;x<15;x+=0.5) {
        for (double y=5;y<15;y+=1) {
            for (double z=5;z<15;z+=1) {
                system->add_particle(x,y,z);
            }
        }
    }

    //file.open("output.csv");

    for (int steps=0;steps<200;steps++) {

        std::vector<glm::dvec3> positions = system->get_positions();
        system->step();
        //std::vector<glm::dvec3> velocities = system->get_velocity();
        for (size_t i=0;i<positions.size();i++) {
            glm::dvec3 pos = positions[i];
            //glm::dvec3 vel = velocities[i];
            //file << pos.x << "," << pos.y << "," << pos.z << endl;
        }
        //cout << positions[0].x << "," << positions[0].y << "," << positions[0].z << endl;


    }
    //file.close();
    */

}