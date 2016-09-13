#version 400 core

// Interpolated values from the vertex shaders
in vec3 posEye;
in float radius;
in mat4 projection;
in float dist;
in float density_out;

// Ouput data
out vec3 color;

void main(){

	// Output color = color specified in the vertex shader, 
	// interpolated between all 3 surrounding vertices
    vec3 N;
	N.xy = gl_PointCoord*2.0-1.0;
	float r2 = dot(N.xy,N.xy);
	if (r2 > 1.0) discard;
	N.z = sqrt(1.0 - r2);

	vec4 pixelPos = vec4(posEye + N*radius, 1.0);
	vec4 clipSpacePos = projection * pixelPos;
	float depth = (clipSpacePos.z / clipSpacePos.w)*0.5f+0.5f;
    float z = (2 * 0.1) / (100 + 0.1 - gl_FragDepth * (100 - 0.1));
    gl_FragDepth = depth;
	//float diffuse = max(0.0, dot(N, vec3(1,1,1)));
	color = .2 * vec3(1 - r2);


}