#pragma OPENCL EXTENSION cl_khr_fp64 : enable

// Define C++ Classes as OpenCL structs

typedef struct __attribute__ ((packed)) {
    uchar4 color;
    float4 points[3];
    float16 mvp;
    float16 mvp_old;
} cl_polygon;


float16 inverce(float16 m) {
    float SubFactor00 = m[10] * m[15] - m[14] * m[11];
	float SubFactor01 = m[9] * m[15] - m[13] * m[11];
	float SubFactor02 = m[9] * m[14] - m[13] * m[10];
	float SubFactor03 = m[8] * m[15] - m[12] * m[11];
	float SubFactor04 = m[8] * m[14] - m[12] * m[10];
	float SubFactor05 = m[8] * m[13] - m[12] * m[9];
	float SubFactor06 = m[6] * m[15] - m[14] * m[7];
	float SubFactor07 = m[5] * m[15] - m[13] * m[7];
	float SubFactor08 = m[5] * m[14] - m[13] * m[6];
	float SubFactor09 = m[4] * m[15] - m[12] * m[7];
	float SubFactor10 = m[4] * m[14] - m[12] * m[6];
	float SubFactor11 = m[4] * m[13] - m[12] * m[5];
	float SubFactor12 = m[6] * m[11] - m[10] * m[7];
	float SubFactor13 = m[5] * m[11] - m[9] * m[7];
	float SubFactor14 = m[5] * m[10] - m[9] * m[6];
	float SubFactor15 = m[4] * m[11] - m[8] * m[7];
	float SubFactor16 = m[4] * m[10] - m[8] * m[6];
	float SubFactor17 = m[4] * m[9] - m[8] * m[5];

    float16 Inverse;
    Inverse[0] = + (m[5] * SubFactor00 - m[6] * SubFactor01 + m[7] * SubFactor02);
    Inverse[1] = - (m[4] * SubFactor00 - m[6] * SubFactor03 + m[7] * SubFactor04);
    Inverse[2] = + (m[4] * SubFactor01 - m[5] * SubFactor03 + m[7] * SubFactor05);
    Inverse[3] = - (m[4] * SubFactor02 - m[5] * SubFactor04 + m[6] * SubFactor05);

    Inverse[4] = - (m[1] * SubFactor00 - m[2] * SubFactor01 + m[3] * SubFactor02);
    Inverse[5] = + (m[0] * SubFactor00 - m[2] * SubFactor03 + m[3] * SubFactor04);
    Inverse[6] = - (m[0] * SubFactor01 - m[1] * SubFactor03 + m[3] * SubFactor05);
    Inverse[7] = + (m[0] * SubFactor02 - m[1] * SubFactor04 + m[2] * SubFactor05);

    Inverse[8] = + (m[1] * SubFactor06 - m[2] * SubFactor07 + m[3] * SubFactor08);
    Inverse[9] = - (m[0] * SubFactor06 - m[2] * SubFactor09 + m[3] * SubFactor10);
    Inverse[10] = + (m[0] * SubFactor07 - m[1] * SubFactor09 + m[3] * SubFactor11);
    Inverse[11] = - (m[0] * SubFactor08 - m[1] * SubFactor10 + m[2] * SubFactor11);

    Inverse[12] = - (m[1] * SubFactor12 - m[2] * SubFactor13 + m[3] * SubFactor14);
    Inverse[13] = + (m[0] * SubFactor12 - m[2] * SubFactor15 + m[3] * SubFactor16);
    Inverse[14] = - (m[0] * SubFactor13 - m[1] * SubFactor15 + m[3] * SubFactor17);
    Inverse[15] = + (m[0] * SubFactor14 - m[1] * SubFactor16 + m[2] * SubFactor17);

    float Determinant =
        + m[0] * Inverse[0]
        + m[1] * Inverse[1]
        + m[2] * Inverse[2]
        + m[3] * Inverse[3];

    Inverse /= Determinant;

    float16 InverseT;
    for(uint i = 0; i < 4; i++) {
        for(uint j = 0; j < 4; j++) {
            InverseT[j + i*4] = Inverse[i + j *4];
        }
    }
	return InverseT;
}

// float16 inverce(float16 m) {
//     float16 inv;
//     float det;
//     int i;

//     inv[0] = m[5]  * m[10] * m[15] - 
//              m[5]  * m[11] * m[14] - 
//              m[9]  * m[6]  * m[15] + 
//              m[9]  * m[7]  * m[14] +
//              m[13] * m[6]  * m[11] - 
//              m[13] * m[7]  * m[10];

//     inv[4] = -m[4]  * m[10] * m[15] + 
//               m[4]  * m[11] * m[14] + 
//               m[8]  * m[6]  * m[15] - 
//               m[8]  * m[7]  * m[14] - 
//               m[12] * m[6]  * m[11] + 
//               m[12] * m[7]  * m[10];

//     inv[8] = m[4]  * m[9] * m[15] - 
//              m[4]  * m[11] * m[13] - 
//              m[8]  * m[5] * m[15] + 
//              m[8]  * m[7] * m[13] + 
//              m[12] * m[5] * m[11] - 
//              m[12] * m[7] * m[9];

//     inv[12] = -m[4]  * m[9] * m[14] + 
//                m[4]  * m[10] * m[13] +
//                m[8]  * m[5] * m[14] - 
//                m[8]  * m[6] * m[13] - 
//                m[12] * m[5] * m[10] + 
//                m[12] * m[6] * m[9];

//     inv[1] = -m[1]  * m[10] * m[15] + 
//               m[1]  * m[11] * m[14] + 
//               m[9]  * m[2] * m[15] - 
//               m[9]  * m[3] * m[14] - 
//               m[13] * m[2] * m[11] + 
//               m[13] * m[3] * m[10];

//     inv[5] = m[0]  * m[10] * m[15] - 
//              m[0]  * m[11] * m[14] - 
//              m[8]  * m[2] * m[15] + 
//              m[8]  * m[3] * m[14] + 
//              m[12] * m[2] * m[11] - 
//              m[12] * m[3] * m[10];

//     inv[9] = -m[0]  * m[9] * m[15] + 
//               m[0]  * m[11] * m[13] + 
//               m[8]  * m[1] * m[15] - 
//               m[8]  * m[3] * m[13] - 
//               m[12] * m[1] * m[11] + 
//               m[12] * m[3] * m[9];

//     inv[13] = m[0]  * m[9] * m[14] - 
//               m[0]  * m[10] * m[13] - 
//               m[8]  * m[1] * m[14] + 
//               m[8]  * m[2] * m[13] + 
//               m[12] * m[1] * m[10] - 
//               m[12] * m[2] * m[9];

//     inv[2] = m[1]  * m[6] * m[15] - 
//              m[1]  * m[7] * m[14] - 
//              m[5]  * m[2] * m[15] + 
//              m[5]  * m[3] * m[14] + 
//              m[13] * m[2] * m[7] - 
//              m[13] * m[3] * m[6];

//     inv[6] = -m[0]  * m[6] * m[15] + 
//               m[0]  * m[7] * m[14] + 
//               m[4]  * m[2] * m[15] - 
//               m[4]  * m[3] * m[14] - 
//               m[12] * m[2] * m[7] + 
//               m[12] * m[3] * m[6];

//     inv[10] = m[0]  * m[5] * m[15] - 
//               m[0]  * m[7] * m[13] - 
//               m[4]  * m[1] * m[15] + 
//               m[4]  * m[3] * m[13] + 
//               m[12] * m[1] * m[7] - 
//               m[12] * m[3] * m[5];

//     inv[14] = -m[0]  * m[5] * m[14] + 
//                m[0]  * m[6] * m[13] + 
//                m[4]  * m[1] * m[14] - 
//                m[4]  * m[2] * m[13] - 
//                m[12] * m[1] * m[6] + 
//                m[12] * m[2] * m[5];

//     inv[3] = -m[1] * m[6] * m[11] + 
//               m[1] * m[7] * m[10] + 
//               m[5] * m[2] * m[11] - 
//               m[5] * m[3] * m[10] - 
//               m[9] * m[2] * m[7] + 
//               m[9] * m[3] * m[6];

//     inv[7] = m[0] * m[6] * m[11] - 
//              m[0] * m[7] * m[10] - 
//              m[4] * m[2] * m[11] + 
//              m[4] * m[3] * m[10] + 
//              m[8] * m[2] * m[7] - 
//              m[8] * m[3] * m[6];

//     inv[11] = -m[0] * m[5] * m[11] + 
//                m[0] * m[7] * m[9] + 
//                m[4] * m[1] * m[11] - 
//                m[4] * m[3] * m[9] - 
//                m[8] * m[1] * m[7] + 
//                m[8] * m[3] * m[5];

//     inv[15] = m[0] * m[5] * m[10] - 
//               m[0] * m[6] * m[9] - 
//               m[4] * m[1] * m[10] + 
//               m[4] * m[2] * m[9] + 
//               m[8] * m[1] * m[6] - 
//               m[8] * m[2] * m[5];

//     det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

//     if (det == 0)
//         return false;

//     det = 1.0 / det;

//     for (i = 0; i < 16; i++)
//         inv[i] = inv[i] * det;

//     return true;
// }

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

void triangle(const float4 *pts, global float *zbuffer, global float2 *velocity, global uchar4 *image, uint w, uint h, uchar4 color, float16 MVP, float16 MVP_old) {
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

    float16 iMVP = inverce(MVP);
    // float16 iMVP = MVP;

    // printf("MVP\n%f %f %f %f \n%f %f %f %f \n%f %f %f %f \n%f %f %f %f \n\n"
    // "iMVP\n%f %f %f %f \n%f %f %f %f \n%f %f %f %f \n%f %f %f %f \n\n", 
    //     MVP[0], MVP[1], MVP[2], MVP[3],
    //     MVP[4], MVP[5], MVP[6], MVP[7],
    //     MVP[8], MVP[9], MVP[10],MVP[11],
    //     MVP[12],MVP[13],MVP[14],MVP[15], 
    //     iMVP[0], iMVP[1], iMVP[2], iMVP[3],
    //     iMVP[4], iMVP[5], iMVP[6], iMVP[7],
    //     iMVP[8], iMVP[9], iMVP[10], iMVP[11],
    //     iMVP[12],iMVP[13], iMVP[14], iMVP[15]
    // );

    // float4 testP = (float4)(5,5,5, 1);
    // float4 testPA = applyMatrix(testP, MVP);
    // float4 testPAR = applyMatrix(testPA, iMVP);


    // printf("%f %f %f", testP.x, testPA.x, testPAR.x);

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

            if (zbuffer[pos] < P.z && P.z < 10.f) {
                // printf("FILL P: %f %f %f, pos: %d, z: %f Color: %u %u %u %u \n", P.x, P.y, P.z, pos, zbuffer[pos],  color.x, color.y, color.z, color.w);

                zbuffer[pos] = P.z;

                // printf("After FILL P: %f %f %f, pos: %d, z: %f \n", P.x, P.y, P.z, pos, zbuffer[pos]);

                image[pos] = color;

                float4 nP = (float4)(P, 1);

                float4 oldP = applyMatrix(applyMatrix(nP, iMVP), MVP_old);

                velocity[pos] = (float2)(oldP.x - P.x, oldP.y - P.y);

            }
        }
    }
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
                   global const cl_polygon *polygons,
                   global float2 *velocity
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

    


    triangle(pointsApplyed, zbuffer, velocity, pixels, wh->x, wh->y, colors[i % 4], P.mvp, P.mvp_old);

}