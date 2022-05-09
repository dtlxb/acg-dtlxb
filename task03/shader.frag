#version 120

// see the GLSL 1.2 specification:
// https://www.khronos.org/registry/OpenGL/specs/gl/GLSLangSpec.1.20.pdf

uniform float time; // current time given from CPU

/**
 * return: signed distance function of a axis aligned box
 * pos: position to evaluate SDF
 * hsize: half size in the XYZ axis
 */
float sdf_box( vec3 pos, vec3 hsize )
{
  vec3 q = abs(pos) - hsize;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

// Here, I tried to analytically find the nearest *small shpere* from pos. It was too complicated and I chose another way.
#define PI 3.14159265358
float round(float x){ 
  float tmp = fract(x);
  if (tmp > 0.5f){
    return ceil(x);
  }
  return floor(x);
}

// Ctrl+CV from IQ's website!!!!
float sdDeathStar( in vec3 p2, in float ra, float rb, in float d ) // d is on x axis.
{ 
  // sampling independent computations (only depend on shape)
  float a = (ra*ra - rb*rb + d*d)/(2.0*d);
  float b = sqrt(max(ra*ra-a*a,0.0));
	
  // sampling dependant computations
  vec2 p = vec2( p2.x, length(p2.yz) );
  if( p.x*b-p.y*a > d*max(b-p.y,0.0) )
    return length(p-vec2(a,b));
  else
    return max( (length(p          )-ra),
               -(length(p-vec2(d,0))-rb));
}

// Definition of singed distance funtion called from
float SDF(vec3 pos)
{
  // for "problem2", write some code below to return the SDF where
  // many small spheres are caved out from a big sphere.
  // The radius of big sphere is `0.8` and its center is at the origin
  // The radius of the small spheres is `0.12` and it's repeaded at intrval of `0.2` in the grid pattern
  // Look Inigo Quilez's article for hints:
  // https://iquilezles.org/articles/distfunctions/

  float dis = length(pos)-0.8;
  dis = -dis;

  float stepp = 0.2;
  float xx = -0.81;
  float yy = -0.82;
  float zz = -0.83;

  // test: inverse + union = intersect
  /*
  dis = -dis;
  float dis2 = length(pos-vec3(0.8,0.0,0.0))-0.12;
  dis2 = -dis2;
  float dis3 = min(dis,dis2);
  return -dis3; */

  /////////////////// MY IDEA: /////////////////////
  // final result = Inverse(
  //                    Union(
  //                          Inverse(deathstar1),Inverse(deathstar2),Inverse(deathstar3)...
  //                          )
  //                       ).
  // To realize these 2 operations,
  // Union(Shape1, Shape2)=min(SDF1,SDF2), and Inverse(Shape1)= -1*SDF(Shape1).
  
  xx = 0.01;
  zz = 0.01;
  yy = 0.01;
  while (xx <= 0.81){
    zz = 0.01;
    yy = 0.01;
    while (zz <= 0.81){
      yy = 0.01;
      while (yy <= 0.81){
        float theta = atan(-zz/yy);
        float a = atan((sin(theta)*zz-cos(theta)*yy) / xx);

        // rotate to x axis, because deathstar function was implemented in x axis for simplicity.
        vec3 v1 = vec3(cos(a),              sin(a),               0);
        vec3 v2 = vec3(-sin(a)*cos(theta),  cos(a)*cos(theta),    sin(theta));
        vec3 v3 = vec3(sin(a)*sin(theta),   -cos(a)*sin(theta),   cos(theta));
        mat3 toxaxis = mat3(v1,v2,v3); // columns

        float tmpdis1 = -sdDeathStar(toxaxis*abs(pos), 0.8, 0.12, length(vec3(xx,yy,zz))); // one small sphere
        dis = min(dis, tmpdis1); // 
        yy = yy + 0.2;
      }
      zz = zz + 0.2;
    }
    xx = xx + 0.2;
  }
  return -dis;

}

void main()
{
  // camera position
  vec3 cam_pos = normalize( vec3(sin(time),cos(time),sin(3*time)) );

  // local frame defined on the cameera
  vec3 frame_z = cam_pos;
  vec3 frame_x = normalize(cross(vec3(0,0,1),frame_z));
  vec3 frame_y = cross(frame_z,frame_x);

  // gl_FragCoord: the coordinate of the pixel
  // left-bottom is (0,0), right-top is (W,H)
  // https://www.khronos.org/registry/OpenGL-Refpages/gl4/html/gl_FragCoord.xhtml
  vec2 scr_xy = gl_FragCoord.xy / vec2(500,500) * 2.0 - vec2(1,1); // canonical screen position [-1,+1] x [-1,+1]
  vec3 src = frame_x * scr_xy.x + frame_y * scr_xy.y + frame_z * 1;  // source of ray from pixel
  vec3 dir = -frame_z;  // direction of ray (looking at the origin)

  vec3 pos_cur = src; // the current ray position
  for(int itr=0;itr<30;++itr){
    float s0 = SDF(pos_cur);
    if( s0 < 0.0 ){ // ray starting from inside the object
      gl_FragColor = vec4(1, 0, 0, 1); // paint red
      return;
    }
    if( s0 < 1.0e-3 ){ // the ray hit the implicit surfacee
      float eps = 1.0e-3;
      float sx = SDF(pos_cur+vec3(eps,0,0))-s0;
      float sy = SDF(pos_cur+vec3(0,eps,0))-s0;
      float sz = SDF(pos_cur+vec3(0,0,eps))-s0;
      vec3 nrm = normalize(vec3(sx,sy,sz)); // normal direction
      float c = -dot(nrm, dir); // Lambersian reflection. The light is at the camera position.
      gl_FragColor = vec4(c, c, c, 1);
      return;
    }
    pos_cur += s0 * dir; // advance ray
  }
  gl_FragColor = vec4(0.5, 0.7, 0.9, 1); // ray doesn't hit the object
}
