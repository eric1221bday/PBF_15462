#version 400 core

// Input vertex data, different for all executions of this shader.
layout(location = 0) in vec3 vertexPosition_modelspace;
layout(location = 1) in float density;

// Output data ; will be interpolated for each fragment.
out vec3 posEye;
out float radius;
out mat4 projection;
out float dist;
out float density_out;

// Values that stay constant for the whole mesh.
uniform mat4 MVP;
uniform mat4 MV;
uniform mat4 P;

void main(){	

	// Output position of the vertex, in clip space : MVP * position
	gl_Position =  MVP * vec4(vertexPosition_modelspace,1);

	posEye = (MV * vec4(vertexPosition_modelspace,1)).xyz;

	dist = length(posEye);

    radius = 1.f;

    gl_PointSize = 1000.f/gl_Position.w;

    projection = P;

    density_out = density;

	// The color of each vertex will be interpolated
	// to produce the color of each fragment
}

