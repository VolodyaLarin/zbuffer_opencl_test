
// Define C++ Classes as OpenCL structs

typedef struct __attribute__ ((packed)) {
    uchar4 color;
    float3 points[3];
} cl_polygon;

float3 barycentric(float3 A, float3 B, float3 C, float3 P) {
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


void triangle(float3 *pts, global float *zbuffer, global uchar4 *image, int w, int h, uchar4 color) {
    float2 bboxmin = (float2)( INFINITY,  INFINITY);
    float2 bboxmax = (float2)( -INFINITY, -INFINITY);
    float2   clamp   = (float2)(w - 1, h - 1);

    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            bboxmin[j] = max(0.f,      min(bboxmin[j], pts[i][j]));
            bboxmax[j] = min(clamp[j], max(bboxmax[j], pts[i][j]));
        }
    }

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
            int pos = P.x+P.y*w;
            

            if (zbuffer[pos] < P.z) {
                zbuffer[pos] = P.z;

                image[pos] = color;
            }
        }
    }
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
        (uchar4)(255),
    };


    // pixels[i] = (uchar4)(i % 255, (i + 50) % 255, i % 255, 255);
    

    cl_polygon p = polygons[i];

    float3 P = p.points[0];
    uchar4 C = p.color;

    printf("%d. %f %f %f \n", i, P.x, P.y, P.z);
    printf("%d. %u %u %u %u \n", i, C.x, C.y, C.z, C.w);



    triangle(p.points, zbuffer, pixels, wh->x, wh->y, colors[i % 4]);

}