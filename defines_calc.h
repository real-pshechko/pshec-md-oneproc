#ifndef __defines_calc__
#define __defines_calc__
//a = (t *) malloc ((n) * sizeof (t))

#define AllocMem(a, n, t)   a = new t [(n)]  //a = (t *) malloc ((n) * sizeof (t))

#define VAdd(v1, v2, v3)									\
	(v1).x = (v2).x + (v3).x,								\
	(v1).y = (v2).y + (v3).x,								\
	(v1).z = (v2).z + (v3).z

#define VSub(v1, v2, v3)									\
	(v1).x = (v2).x - (v3).x,								\
	(v1).y = (v2).y - (v3).x,								\
	(v1).z = (v2).z - (v3).z

#define VDot(v1, v2)                                        \
   ((v1).x * (v2).x + (v1).y * (v2).y + (v1).z * (v2).z)

#define VSAdd(v1, v2, s3, v3)                               \
   (v1).x = (v2).x + (s3) * (v3).x,                         \
   (v1).y = (v2).y + (s3) * (v3).y,                         \
   (v1).z = (v2).z + (s3) * (v3).z

#define VSet(v, sx, sy, sz)                                 \
   (v).x = sx,                                              \
   (v).y = sy,                                              \
   (v).z = sz

#define VSetAll(v, s) VSet (v, s, s, s)
#define VZero(v)  VSetAll (v, 0)
#define VLenSq(v)  VDot (v, v)
#define VVSAdd(v1, s2, v2) VSAdd (v1, v1, s2, v2)

#define Sqr(x) ((x) * (x))
#define Cube(x) ((x) * (x) * (x))
#define DO_MOL for (n = 0; n < nMol; n ++)

#define VWrap(v, t)                                         \
   if (v.t >= 0.5 * region.t)      v.t -= region.t;         \
   else if (v.t < -0.5 * region.t) v.t += region.t
#define VWrapAll(v)                                         \
   {VWrap (v, x);                                           \
   VWrap (v, y);                                            \
   VWrap (v, z);}

#endif