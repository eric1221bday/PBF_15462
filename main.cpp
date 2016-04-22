#include <iostream>
#include "ParticleSystem.h"
#include <fstream>

GLFWwindow* window;

using namespace std;
using namespace glm;

ParticleSystem *PBF = new ParticleSystem();

void init_liquid() {
    for (double x=0.1;x<20;x+=1) {
        for (double y=0.1;y<20;y+=1) {
            for (double z=0.1;z<12;z+=1) {
                PBF->add_particle(x,y,z);
            }
        }
    }
}

std::vector <GLfloat> get_particles() {
    std::vector <GLfloat> positions;
    //cout << PBF->scalar_field[1 * (PBF->jmax * PBF->imax) + 5 * (PBF->imax) + 5] << endl;
    for (auto i : PBF->particles) {
        positions.push_back(float(i->x.x));
        positions.push_back(float(i->x.z));
        positions.push_back(float(i->x.y));
    }
    return positions;
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

    int windowWidth = 1280;
    int windowHeight = 720;

    // Open a window and create its OpenGL context
    window = glfwCreateWindow( windowWidth, windowHeight, "Tutorial 04 - Colored Cube", NULL, NULL);
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
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

    // Enable depth test
    glEnable(GL_DEPTH_TEST);
    // Accept fragment if it closer to the camera than the former one
    glDepthFunc(GL_LESS);

    glEnable(GL_PROGRAM_POINT_SIZE);

    //glPointParameteri(GL_POINT_SPRITE_COORD_ORIGIN,GL_LOWER_LEFT);

    GLuint VertexArrayID;
    glGenVertexArrays(1, &VertexArrayID);
    glBindVertexArray(VertexArrayID);

    // Create and compile our GLSL program from the shaders
    GLuint programID = LoadShaders( "DepthVertex.glsl", "DepthFragment.glsl" );

    // Get a handle for our "MVP" uniform
    GLuint MatrixID = glGetUniformLocation(programID, "MVP");

    GLuint MVID = glGetUniformLocation(programID, "MV");

    GLuint PID = glGetUniformLocation(programID, "P");

    // Our vertices. Three consecutive floats give a 3D vertex; Three consecutive vertices give a triangle.
    // A cube has 6 faces with 2 triangles each, so this makes 6*2=12 triangles, and 12*3 vertices

    init_liquid();

    PBF->step();

    std::vector <GLfloat> vertices = get_particles();

    GLuint vertexbuffer;
    glGenBuffers(1, &vertexbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    glBufferData(GL_ARRAY_BUFFER, vertices.size()*sizeof(GLfloat), NULL, GL_STREAM_DRAW);

    GLuint FramebufferName = 0;
    glGenFramebuffers(1, &FramebufferName);
    glBindFramebuffer(GL_FRAMEBUFFER, FramebufferName);

    // The texture we're going to render to
    GLuint depthTexture;
    glGenTextures(1, &depthTexture);
    glBindTexture(GL_TEXTURE_2D, depthTexture);
    glTexImage2D(GL_TEXTURE_2D, 0,GL_DEPTH_COMPONENT, windowWidth, windowHeight, 0,GL_DEPTH_COMPONENT, GL_FLOAT, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_R_TO_TEXTURE);

    glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, depthTexture, 0);

    // No color output in the bound framebuffer, only depth.
    glDrawBuffer(GL_NONE);

    // Always check that our framebuffer is ok
    if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
        return false;

    GLuint BlurFramebufferName = 0;
    glGenFramebuffers(1, &BlurFramebufferName);
    glBindFramebuffer(GL_FRAMEBUFFER, BlurFramebufferName);
    GLuint blurTexture;
    glGenTextures(1, &blurTexture);
    glBindTexture(GL_TEXTURE_2D, blurTexture);
    glTexImage2D(GL_TEXTURE_2D, 0,GL_DEPTH_COMPONENT, windowWidth, windowHeight, 0,GL_DEPTH_COMPONENT, GL_FLOAT, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_R_TO_TEXTURE);
    glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, blurTexture, 0);
    glDrawBuffer(GL_NONE);

    // Always check that our framebuffer is ok
    if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
        return false;

    // The fullscreen quad's FBO
    static const GLfloat g_quad_vertex_buffer_data[] = {
            -1.0f, -1.0f, 0.0f,
            1.0f, -1.0f, 0.0f,
            -1.0f,  1.0f, 0.0f,
            -1.0f,  1.0f, 0.0f,
            1.0f, -1.0f, 0.0f,
            1.0f,  1.0f, 0.0f,
    };

    GLuint quad_vertexbuffer;
    glGenBuffers(1, &quad_vertexbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, quad_vertexbuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(g_quad_vertex_buffer_data), g_quad_vertex_buffer_data, GL_STATIC_DRAW);

    // Create and compile our GLSL program from the shaders
    GLuint quad_programID = LoadShaders( "Passthrough.glsl", "FinalTexture.glsl" );
    GLuint texID = glGetUniformLocation(quad_programID, "renderedTexture");
    GLuint quad_MatrixID = glGetUniformLocation(quad_programID, "MVP");
    GLuint quad_MVID = glGetUniformLocation(quad_programID, "MV");
    GLuint quad_PID = glGetUniformLocation(quad_programID, "P");

    GLuint blur_programID = LoadShaders( "Passthrough.glsl", "BlurTexture.glsl" );
    GLuint blur_texID = glGetUniformLocation(blur_programID, "renderedTexture");
    GLuint blur_dirID = glGetUniformLocation(blur_programID, "direction");

    do{

        glBindFramebuffer(GL_FRAMEBUFFER, FramebufferName);
        glViewport(0,0,windowWidth,windowHeight);

        vertices = get_particles();
        glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
        glBufferData(GL_ARRAY_BUFFER, vertices.size()*sizeof(GLfloat), NULL, GL_STREAM_DRAW);
        glBufferSubData(GL_ARRAY_BUFFER, 0, vertices.size()*sizeof(GLfloat), &vertices.front());

        // Clear the screen
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Use our shader
        glUseProgram(programID);

        computeMatricesFromInputs();
        glm::mat4 ProjectionMatrix = getProjectionMatrix();
        glm::mat4 ViewMatrix = getViewMatrix();
        glm::mat4 ModelMatrix = glm::mat4(1.0);
        glm::mat4 MVP = ProjectionMatrix * ViewMatrix * ModelMatrix;
        glm::mat4 MV = ViewMatrix * ModelMatrix;

        // Send our transformation to the currently bound shader,
        // in the "MVP" uniform
        glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);
        glUniformMatrix4fv(MVID, 1, GL_FALSE, &MV[0][0]);
        glUniformMatrix4fv(PID, 1, GL_FALSE, &ProjectionMatrix[0][0]);

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

        // Draw the triangle !
        glDrawArrays(GL_POINTS, 0, vertices.size()/3); // 12*3 indices starting at 0 -> 12 triangles

        glDisableVertexAttribArray(0);



        glBindFramebuffer(GL_FRAMEBUFFER, BlurFramebufferName);
        glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, blurTexture, 0);
        // Render on the whole framebuffer, complete from the lower left corner to the upper right
        glViewport(0,0,windowWidth,windowHeight);

        // Clear the screen
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Use our shader
        glUseProgram(blur_programID);

        // Bind our texture in Texture Unit 0
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, depthTexture);
        // Set our "renderedTexture" sampler to user Texture Unit 0
        glUniform1i(blur_texID, 0);
        glUniform1i(blur_dirID, 0);

        // 1rst attribute buffer : vertices
        glEnableVertexAttribArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, quad_vertexbuffer);
        glVertexAttribPointer(
                0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
                3,                  // size
                GL_FLOAT,           // type
                GL_FALSE,           // normalized?
                0,                  // stride
                (void*)0            // array buffer offset
        );

        // Draw the triangles !
        glDrawArrays(GL_TRIANGLES, 0, 6); // 2*3 indices starting at 0 -> 2 triangles

        glDisableVertexAttribArray(0);



        glBindFramebuffer(GL_FRAMEBUFFER, FramebufferName);
        glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, depthTexture, 0);
        // Render on the whole framebuffer, complete from the lower left corner to the upper right
        glViewport(0,0,windowWidth,windowHeight);

        // Clear the screen
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Use our shader
        glUseProgram(blur_programID);

        // Bind our texture in Texture Unit 0
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, blurTexture);
        // Set our "renderedTexture" sampler to user Texture Unit 0
        glUniform1i(blur_texID, 0);
        glUniform1i(blur_dirID, 1);

        // 1rst attribute buffer : vertices
        glEnableVertexAttribArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, quad_vertexbuffer);
        glVertexAttribPointer(
                0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
                3,                  // size
                GL_FLOAT,           // type
                GL_FALSE,           // normalized?
                0,                  // stride
                (void*)0            // array buffer offset
        );

        // Draw the triangles !
        glDrawArrays(GL_TRIANGLES, 0, 6); // 2*3 indices starting at 0 -> 2 triangles

        glDisableVertexAttribArray(0);



        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        // Render on the whole framebuffer, complete from the lower left corner to the upper right
        glViewport(0,0,windowWidth*2,windowHeight*2);

        // Clear the screen
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Use our shader
        glUseProgram(quad_programID);

        // Bind our texture in Texture Unit 0
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, depthTexture);
        // Set our "renderedTexture" sampler to user Texture Unit 0
        glUniform1i(texID, 0);
        glUniformMatrix4fv(quad_MatrixID, 1, GL_FALSE, &MVP[0][0]);
        glUniformMatrix4fv(quad_MVID, 1, GL_FALSE, &MV[0][0]);
        glUniformMatrix4fv(quad_PID, 1, GL_FALSE, &ProjectionMatrix[0][0]);


        // 1rst attribute buffer : vertices
        glEnableVertexAttribArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, quad_vertexbuffer);
        glVertexAttribPointer(
                0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
                3,                  // size
                GL_FLOAT,           // type
                GL_FALSE,           // normalized?
                0,                  // stride
                (void*)0            // array buffer offset
        );

        // Draw the triangles !
        glDrawArrays(GL_TRIANGLES, 0, 6); // 2*3 indices starting at 0 -> 2 triangles

        glDisableVertexAttribArray(0);
        if (glfwGetKey( window, GLFW_KEY_P ) == GLFW_PRESS){
            PBF->step();
        }

        // Swap buffers
        glfwSwapBuffers(window);
        glfwPollEvents();

    } // Check if the ESC key was pressed or the window was closed
    while( glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
           glfwWindowShouldClose(window) == 0 );

    // Cleanup VBO and shader
    glDeleteBuffers(1, &vertexbuffer);
    glDeleteProgram(programID);
    glDeleteVertexArrays(1, &VertexArrayID);
    glDeleteFramebuffers(1, &FramebufferName);
    glDeleteFramebuffers(1, &BlurFramebufferName);
    glDeleteTextures(1, &depthTexture);
    glDeleteTextures(1, &blurTexture);

    glDeleteBuffers(1, &quad_vertexbuffer);

    // Close OpenGL window and terminate GLFW
    glfwTerminate();

    return 0;

}