#pragma OPENCL EXTENSION cl_khr_fp64 : enable

// Define C++ Classes as OpenCL structs

typedef struct __attribute__ ((packed)) {
    uchar4 color;
    float4 points[3];
    float16 mvp;
} cl_polygon;

float3 __cross__(float3 v1, float3 v2) {

    return (float3)(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}

float3 barycentric(float4 A, float4 B, float4 C, float3 P) {
    float3 s[2];
    for (int i=2; i--; ) {
        s[i][0] = C[i]-A[i];
        s[i][1] = B[i]-A[i];
        s[i][2] = A[i]-P[i];
    }

    float3 u = cross(s[0], s[1]);

    float uz = u.z;
    if (uz < 0) {
        uz = -uz;
    } 
    
    if ( uz > 1e-2 ) 
    {
        return (float3)(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
    }
    return (float3)( -1, 1, 1);
}


void _swap_float4(float4 *a, float4 *b) {
    float4 tmp = *b;
    *b = *a;
    *a = tmp;
} 

// void triangle(const float4 *pts, global float *zbuffer, global uchar4 *image, uint w, uint h, uchar4 color) {
//     float4 t0 = pts[0];
//     float4 t1 = pts[1];
//     float4 t2 = pts[2];

//     if (t0.y==t1.y && t0.y==t2.y) return; // i dont care about degenerate triangles
    
//     if (t0.y>t1.y) _swap_float4(&t0, &t1);
//     if (t0.y>t2.y) _swap_float4(&t0, &t2);
//     if (t1.y>t2.y) _swap_float4(&t1, &t2);

//     int total_height = t2.y-t0.y;
//     for (int i=0; i < total_height; i++) {
//         bool second_half = i>t1.y-t0.y || t1.y==t0.y;
//         int segment_height = second_half ? t2.y-t1.y : t1.y-t0.y;
//         float alpha = (float)i/total_height;
//         float beta  = (float)(i-(second_half ? t1.y-t0.y : 0))/segment_height; // be careful: with above conditions no division by zero here
//         float4 A =               t0 + (float4)(t2-t0)*alpha;
//         float4 B = second_half ? t1 + (float4)(t2-t1)*beta : t0 + (float4)(t1-t0)*beta;
        
//         if (A.x>B.x) _swap_float4(&A, &B);
//         for (int j=A.x; j<=B.x; j++) {
//             float phi = B.x==A.x ? 1. : (float)(j-A.x)/(float)(B.x-A.x);
//             float4 P = (float4)(A) + (float4)(B-A)*phi;

//             if (P.x >= w || P.y >= h || P.x < 0 || P.y < 0) continue;

//             int idx = P.x + P.y * w;
//             if (zbuffer[idx] < P.z) {
//                 zbuffer[idx] = P.z;
//                 image[idx] = color;
//             }
//         }
//     }
// }

float4 normal(const float4 *pts){
    float4 a = pts[0] - pts[1];    
    float4 b = pts[0] - pts[2];    
    
    float4 res = cross(a,b);
    res = normalize(res);

    float4 i = -pts[0];
    if (acos(dot(res, i) / length(res) / length(i)) > 0)
        res = -res;

    return res;
}

void triangle(const float4 *pts, global float *zbuffer, global uchar4 *image, uint w, uint h, uchar4 color) {
    float3 LightVec = (float3)(3 / 5.0f, -4 / 5.0f, 0);

    float4 n = normal(pts);
  
    float intens =
            0.2f + 0.8f - cos((
                    n.x * LightVec.x +
                    n.y * LightVec.y +
                    n.z * LightVec.z
            ) / M_PI) * 0.8f;


   
    // color = (uchar4)( color.x * intens, color.y * intens, color.z * intens, 255);

    float2 bboxmin = (float2)( INFINITY,  INFINITY);
    float2 bboxmax = (float2)( -INFINITY, -INFINITY);
    float2   clamp   = (float2)(w - 1, h - 1);

    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            bboxmin[j] = max(0.f,      min(bboxmin[j], pts[i][j]));
            bboxmax[j] = min(clamp[j], max(bboxmax[j], pts[i][j]));
        }
    }

    // printf("bboxmin: %f %f, bboxmax: %f %f\n", bboxmin.x, bboxmin.y, bboxmax.x, bboxmax.y);

    // printf(" %f %f %f\n %f %f %f\n %f %f %f\n\n",
    //      pts[0].x, pts[0].y, pts[0].z,
    //      pts[1].x, pts[1].y, pts[1].z,
    //      pts[2].x, pts[2].y, pts[2].z
    //     );

    float3 P;
    for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
        for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
            float3 bc_screen  = barycentric(pts[0], pts[1], pts[2], P);
            if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) {
                continue;
            }

            P.z = 0;
            for (int i=0; i<3; i++) {
                P.z += pts[i][2] * bc_screen[i];
            }
            uint pos = (uint)(P.x) + (uint)(P.y) * w;
            
            // printf("P: %f %f %f, pos: %d, z: %f \n", P.x, P.y, P.z, pos, zbuffer[pos]);

            if (zbuffer[pos] < P.z) {
                // printf("FILL P: %f %f %f, pos: %d, z: %f Color: %u %u %u %u \n", P.x, P.y, P.z, pos, zbuffer[pos],  color.x, color.y, color.z, color.w);

                zbuffer[pos] = P.z;

                // printf("After FILL P: %f %f %f, pos: %d, z: %f \n", P.x, P.y, P.z, pos, zbuffer[pos]);

                image[pos] = color;

            }
        }
    }
}


float4 applyMatrix(float4 p, float16 mvp) {
    float4 res = (float4)(
        p.x * mvp[0] +  p.y * mvp[4] + p.z * mvp[8 ] + p.w * mvp[12],
        p.x * mvp[1] +  p.y * mvp[5] + p.z * mvp[9 ] + p.w * mvp[13], 
        p.x * mvp[2] +  p.y * mvp[6] + p.z * mvp[10] + p.w * mvp[14], 
        p.x * mvp[3] +  p.y * mvp[7] + p.z * mvp[11] + p.w * mvp[15]
    );

    res.x /= res.w;
    res.y /= res.w;
    // res.z /= res.w;

    res.w = 1;
    
    return res;
}


// float4 applyMatrix(float4 p, float16 mvp) {
//     float4 res = (float4)(
//         p.x * mvp[0] +  p.y * mvp[1] + p.z * mvp[2 ] + p.w * mvp[3],
//         p.x * mvp[4] +  p.y * mvp[5] + p.z * mvp[6 ] + p.w * mvp[7], 
//         p.x * mvp[8] +  p.y * mvp[9] + p.z * mvp[10] + p.w * mvp[11], 
//         p.x * mvp[12] +  p.y * mvp[13] + p.z * mvp[14] + p.w * mvp[15]
//     );

//     res.x /= res.w;
//     res.y /= res.w;
//     res.z /= res.w;

//     res.w = 1;
    
//     return res;
// }

float4 world2screen(float4 v, uint2 wh) {
    v.x += wh.x / 2;
    v.y += wh.y / 2;
    return v;
}


kernel void render(global uchar4 *pixels,
                   global float *zbuffer,
                   global const uint2 *wh,
                   global const cl_polygon *polygons
                   ) {

    int i = get_global_id(0);
    
    uchar4 colors[] = {
        (uchar4)(255, 0, 0,255),
        (uchar4)(0, 255, 0, 255),
        (uchar4)(0,0,255, 255),
        (uchar4)(255,0,255, 255),
    };


    // pixels[i] = (uchar4)(i % 255, (i + 50) % 255, i % 255, 255);

    cl_polygon P = polygons[i];
    const global float4 *p = polygons[i].points;

    
    // printf("%f %f %f %f \n%f %f %f %f \n%f %f %f %f \n%f %f %f %f \n\n", 
    //     P.mvp[0], P.mvp[1], P.mvp[2], P.mvp[3],
    //     P.mvp[4], P.mvp[5], P.mvp[6], P.mvp[7],
    //     P.mvp[8], P.mvp[9], P.mvp[10], P.mvp[11],
    //     P.mvp[12], P.mvp[13], P.mvp[14], P.mvp[15]
    // );



    float4 points[3] = {
        world2screen(P.points[0], *wh),
        world2screen(P.points[1], *wh),
        world2screen(P.points[2], *wh),
    };

    float4 pointsApplyed[3] = {
        applyMatrix(points[0], P.mvp),
        applyMatrix(points[1], P.mvp),
        applyMatrix(points[2], P.mvp)
    };

    


    triangle(pointsApplyed, zbuffer, pixels, wh->x, wh->y, colors[i % 4]);

}