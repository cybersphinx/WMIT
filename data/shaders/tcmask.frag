#version 120
#pragma debug(on)

varying float vertexDistance;
varying vec3 normal, lightDir, eyeVec;

uniform sampler2D Texture0;
uniform sampler2D Texture1;
uniform sampler2D Texture2;
uniform vec4 teamcolour;
uniform int tcmask, normalmap;
uniform int fogEnabled;
uniform bool ecmEffect;
uniform float graphicsCycle;

void main(void)
{
	vec4 mask, colour;

	// Get color from texture unit 0, merge with lighting
	colour = texture2D(Texture0, gl_TexCoord[0].st);

	if (tcmask == 1)
	{
		// Get tcmask information from texture unit 1
		mask = texture2D(Texture1, gl_TexCoord[0].st);
	
		// Apply color using grain merge with tcmask
		gl_FragColor = (colour + (teamcolour - 0.5) * mask.a) * gl_Color;
	}
	else
	{
		gl_FragColor = colour * gl_Color;
	}

	if (ecmEffect)
	{
		gl_FragColor.a = 0.45 + 0.225 * graphicsCycle;
	}

	if (fogEnabled > 0)
	{
		// Calculate linear fog
		float fogFactor = (gl_Fog.end - vertexDistance) / (gl_Fog.end - gl_Fog.start);
		fogFactor = clamp(fogFactor, 0.0, 1.0);
	
		// Return fragment color
		gl_FragColor = mix(gl_Fog.color, gl_FragColor, fogFactor);
	}
}
