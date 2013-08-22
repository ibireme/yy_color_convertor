//
//  yy_color_converter.c
//
//  Created by ibireme on 13-8-17.
//


#include "yy_color_converter.h"


#define CLAMP_COLOR_VALUE(v) (v) = (v) < 0 ? 0 : (v) > 1 ? 1 : (v)


//////////////////////////////////////////////////////////////////////////////// RGB <-> HSL

void RGB2HSL(CGFloat r, CGFloat g, CGFloat b,
             CGFloat *h, CGFloat *s, CGFloat *l) {
    CLAMP_COLOR_VALUE(r);
    CLAMP_COLOR_VALUE(g);
    CLAMP_COLOR_VALUE(b);
    
    CGFloat max, min, delta, sum;
    max = fmaxf(r, fmaxf(g, b));
    min = fminf(r, fminf(g, b));
    delta = max - min;
    sum = max + min;
    
    *l = sum / 2;           // Lightness
    if (delta == 0) {       // No Saturation, so Hue is undefined (achromatic)
        *h = *s = 0;
        return;
    }
    *s = delta / (sum < 1 ? sum : 2 - sum);             // Saturation
    if      (r == max) *h = (g - b) / delta / 6;        // color between y & m
    else if (g == max) *h = (2 + (b - r) / delta) / 6;  // color between c & y
    else               *h = (4 + (r - g) / delta) / 6;  // color between m & y
    if (*h < 0) *h += 1;
}


void HSL2RGB(CGFloat h, CGFloat s, CGFloat l,
             CGFloat *r, CGFloat *g, CGFloat *b) {
    CLAMP_COLOR_VALUE(h);
    CLAMP_COLOR_VALUE(s);
    CLAMP_COLOR_VALUE(l);
    
    if (s == 0) { // No Saturation, Hue is undefined (achromatic)
        *r = *g = *b = l;
        return;
    }
    
    CGFloat q;
    q = (l <= 0.5) ? (l * (1 + s)) : (l + s - (l * s));
    if (q <= 0) {
        *r = *g = *b = 0.0;
    } else {
        int sextant;
        CGFloat m, sv, fract, vsf, mid1, mid2;
        m = l + l - q;
        sv = (q - m) / q;
        if (h == 1) h = 0;
        h *= 6.0;
        sextant = h;
        fract = h - sextant;
        vsf = q * sv * fract;
        mid1 = m + vsf;
        mid2 = q - vsf;
        switch (sextant) {
            case 0: *r = q; *g = mid1; *b = m; break;
            case 1: *r = mid2; *g = q; *b = m; break;
            case 2: *r = m; *g = q; *b = mid1; break;
            case 3: *r = m; *g = mid2; *b = q; break;
            case 4: *r = mid1; *g = m; *b = q; break;
            case 5: *r = q; *g = m; *b = mid2; break;
        }
    }
}


//////////////////////////////////////////////////////////////////////////////// RGB <-> HSB


void RGB2HSB(CGFloat r, CGFloat g, CGFloat b,
             CGFloat *h, CGFloat *s, CGFloat *v) {
    CLAMP_COLOR_VALUE(r);
    CLAMP_COLOR_VALUE(g);
    CLAMP_COLOR_VALUE(b);
    
    CGFloat max, min, delta;
    max = fmaxf(r, fmaxf(g, b));
    min = fminf(r, fminf(g, b));
    delta = max - min;
    
    *v = max;               // Brightness
    if (delta == 0) {       // No Saturation, so Hue is undefined (achromatic)
        *h = *s = 0;
        return;
    }
    *s = delta / max;       // Saturation
    
    if (r == max) *h = (g - b) / delta / 6;             // color between y & m
    else if (g == max) *h = (2 + (b - r) / delta) / 6;  // color between c & y
    else *h = (4 + (r - g) / delta) / 6;                // color between m & c
    if (*h < 0) *h += 1;
}


void HSB2RGB(CGFloat h, CGFloat s, CGFloat v,
             CGFloat *r, CGFloat *g, CGFloat *b) {
    CLAMP_COLOR_VALUE(h);
    CLAMP_COLOR_VALUE(s);
    CLAMP_COLOR_VALUE(v);
    
    if (s == 0) {
        *r = *g = *b = v; // No Saturation, so Hue is undefined (Achromatic)
    } else {
        int sextant;
        CGFloat f, p, q, t;
        if (h == 1) h = 0;
        h *= 6;
        sextant = floorf(h);
        f = h - sextant;
        p = v * (1 - s);
        q = v * (1 - s * f);
        t = v * (1 - s * (1 - f) );
        switch (sextant) {
            case 0: *r = v; *g = t; *b = p; break;
            case 1: *r = q; *g = v; *b = p; break;
            case 2: *r = p; *g = v; *b = t; break;
            case 3: *r = p; *g = q; *b = v; break;
            case 4: *r = t; *g = p; *b = v; break;
            case 5: *r = v; *g = p; *b = q; break;
        }
    }
}


//////////////////////////////////////////////////////////////////////////////// RGB <-> CMYK


void RGB2CMYK(CGFloat r, CGFloat g, CGFloat b,
              CGFloat *c, CGFloat *m, CGFloat *y, CGFloat *k) {
    CLAMP_COLOR_VALUE(r);
    CLAMP_COLOR_VALUE(g);
    CLAMP_COLOR_VALUE(b);
    
    *c = 1 - r;
    *m = 1 - g;
    *y = 1 - b;
    *k = fminf(*c, fminf(*m, *y));
    
    if (*k == 1) {
        *c = *m = *y = 0;   // Pure black
    } else {
        *c = (*c - *k) / (1 - *k);
        *m = (*m - *k) / (1 - *k);
        *y = (*y - *k) / (1 - *k);
    }
}


void CMYK2RGB(CGFloat c, CGFloat m, CGFloat y, CGFloat k,
              CGFloat *r, CGFloat *g, CGFloat *b) {
    CLAMP_COLOR_VALUE(c);
    CLAMP_COLOR_VALUE(m);
    CLAMP_COLOR_VALUE(y);
    CLAMP_COLOR_VALUE(k);
    
    *r = (1 - c) * (1 - k);
    *g = (1 - m) * (1 - k);
    *b = (1 - y) * (1 - k);
}




//////////////////////////////////////////////////////////////////////////////// RGB <-> CMY
void RGB2CMY(CGFloat r, CGFloat g, CGFloat b,
             CGFloat *c, CGFloat *m, CGFloat *y) {
    CLAMP_COLOR_VALUE(r);
    CLAMP_COLOR_VALUE(g);
    CLAMP_COLOR_VALUE(b);
    
    *c = 1 - r;
    *m = 1 - g;
    *y = 1 - b;
}

void CMY2RGB(CGFloat c, CGFloat m, CGFloat y,
             CGFloat *r, CGFloat *g, CGFloat *b) {
    CLAMP_COLOR_VALUE(c);
    CLAMP_COLOR_VALUE(m);
    CLAMP_COLOR_VALUE(y);
    
    *r = 1 - c;
    *g = 1 - m;
    *b = 1 - y;
}



//////////////////////////////////////////////////////////////////////////////// CMYK <-> CMY
void CMY2CMYK(CGFloat c, CGFloat m, CGFloat y,
              CGFloat *cc, CGFloat *mm, CGFloat *yy, CGFloat *kk) {
    CLAMP_COLOR_VALUE(c);
    CLAMP_COLOR_VALUE(m);
    CLAMP_COLOR_VALUE(y);
    
    *kk = fminf(c, fminf(m, y));
    if (*kk == 1) {
        *cc = *mm = *yy = 0;   // Pure black
    } else {
        *cc = (c - *kk) / (1 - *kk);
        *mm = (m - *kk) / (1 - *kk);
        *yy = (y - *kk) / (1 - *kk);
    }
}

void CMYK2CMY(CGFloat c, CGFloat m, CGFloat y, CGFloat k,
              CGFloat *cc, CGFloat *mm, CGFloat *yy) {
    CLAMP_COLOR_VALUE(c);
    CLAMP_COLOR_VALUE(m);
    CLAMP_COLOR_VALUE(y);
    CLAMP_COLOR_VALUE(k);
    
    *cc = c * (1 - k) + k;
    *mm = m * (1 - k) + k;
    *yy = y * (1 - k) + k;
}



//////////////////////////////////////////////////////////////////////////////// HSB <-> HSL


void HSB2HSL(CGFloat h, CGFloat s, CGFloat b,
             CGFloat *hh, CGFloat *ss, CGFloat *ll) {
    CLAMP_COLOR_VALUE(h);
    CLAMP_COLOR_VALUE(s);
    CLAMP_COLOR_VALUE(b);
    
    *hh = h;
    *ll = (2 - s) * b / 2;
    if (*ll <= 0.5) {
        *ss = (s) / ((2 - s));
    } else {
        *ss = (s * b) / (2 - (2 - s) * b);
    }
}


void HSL2HSB(CGFloat h, CGFloat s, CGFloat l,
             CGFloat *hh, CGFloat *ss, CGFloat *bb) {
    CLAMP_COLOR_VALUE(h);
    CLAMP_COLOR_VALUE(s);
    CLAMP_COLOR_VALUE(l);
    
    *hh = h;
    if (l <= 0.5) {
        *bb = (s + 1) * l;
        *ss = (2 * s) / (s + 1);
    } else {
        *bb = l + s * (1 - l);
        *ss = (2 * s * (1 - l)) / *bb;
    }
}






//////////////////////////////////////////////////////////////////////////////// CHEXYZ <-> CIExyY


/// XYZ:0~1 ** xyY:0~1
void CIEXYZ2CIExyY(CGFloat X,CGFloat Y, CGFloat Z,
                   CGFloat *xx,CGFloat *yy, CGFloat *YY) {
    *YY = Y;
    CGFloat sum = X + Y + Z;
    if (sum > 0.0) {
        *xx = X / sum;
        *yy = Y / sum;
    } else {
        *xx = *yy = 0;
    }
}

void CIExyY2CIEXYZ(CGFloat x, CGFloat y, CGFloat Y,
                   CGFloat *XX, CGFloat *YY, CGFloat *ZZ) {
    *YY = Y;
    if (y == 0) {
        *XX = *YY = *ZZ = 0;
    } else {
        *XX = x * Y / y;
        *ZZ =  (1 - x - y) * Y / y;
    }
}





void CIEXYZ2HunterLab(CGFloat X, CGFloat Y, CGFloat Z,
                      CGFloat *L, CGFloat *a, CGFloat *b) {
    *L = 10 * sqrt(Y);
    *a = 17.5 * ((1.02 * X) - Y) / sqrt(Y);
    *b = 7 * (Y - (0.847 * Z)) / sqrt(Y);
}

void HunterLab2CIEXYZ(CGFloat L, CGFloat a, CGFloat b,
                      CGFloat *X, CGFloat *Y, CGFloat *Z) {
    CGFloat x = L / 10;
    CGFloat y = a / 17.5 * L / 10;
    CGFloat z = b / 7 * L / 10;
    *Y = y * y;
    *X = (x + *Y) / 1.02;
    *Z = -(z - *Y) / 0.847;
}













////////////////////////////////////////////////////////////////////////////////Define


typedef struct {
    CGFloat x, y, z;
} yy_triple_t;


typedef struct{
    CGFloat m00,m01,m02;
    CGFloat m10,m11,m12;
    CGFloat m20,m21,m22;
}yy_matrix_t;


typedef struct {
    CGFloat xr,yr,xg,yg,xb,yb;
} yy_xyrgb_t;


static yy_triple_t illum_define [] = {
    { 1.09850, 1.0, 0.35585 }, //A (ASTM E308-01)
    { 0.99072, 1.0, 0.85223 }, //B (Wyszecki & Stiles, p. 769)
    { 0.98074, 1.0, 1.18232 }, //C (ASTM E308-01)
    { 0.96422, 1.0, 0.82521 }, //D50 (ASTM E308-01)
    { 0.95682, 1.0, 0.92149 }, //D55 (ASTM E308-01)
    { 0.95047, 1.0, 1.08883 }, //D65 (ASTM E308-01)
    { 0.94972, 1.0, 1.22638 }, //D75 (ASTM E308-01)
    { 1.00000, 1.0, 1.00000 }, //E (ASTM E308-01)
    { 0.99186, 1.0, 0.67393 }, //F2 (ASTM E308-01)
    { 0.95041, 1.0, 1.08747 }, //F7 (ASTM E308-01)
    { 1.00962, 1.0, 0.64350 }  //F11 (ASTM E308-01)
};


static float gamma_define[] = {
    1.0,
    1.8,
    2.2,
   -2.2,
    0
};


static yy_xyrgb_t xyrgb_define[] = {
    { 0.64,   0.33,   0.21,   0.71,   0.15,   0.06   },
    { 0.625,  0.340,  0.280,  0.595,  0.155,  0.070  },
    { 0.7347, 0.2653, 0.2150, 0.7750, 0.1300, 0.0350 },
    { 0.6888, 0.3112, 0.1986, 0.7551, 0.1265, 0.0352 },
    { 0.64,   0.33,   0.28,   0.65,   0.15,   0.06   },
    { 0.735,  0.265,  0.274,  0.717,  0.167,  0.009  },
    { 0.630,  0.340,  0.295,  0.605,  0.150,  0.075  },
    { 0.696,  0.300,  0.215,  0.765,  0.130,  0.035  },
    { 0.67,   0.33,   0.21,   0.71,   0.14,   0.08   },
    { 0.695,  0.305,  0.260,  0.700,  0.110,  0.005  },
    { 0.67,   0.33,   0.21,   0.71,   0.14,   0.08   },
    { 0.64,   0.33,   0.29,   0.60,   0.15,   0.06   },
    { 0.7347, 0.2653, 0.1596, 0.8404, 0.0366, 0.0001 },
    { 0.630,  0.340,  0.310,  0.595,  0.155,  0.070  },
    { 0.64,   0.33,   0.30,   0.60,   0.15,   0.06   }, ///default
    { 0.735,  0.265,  0.115,  0.826,  0.157,  0.018  }
};

static yy_triple_t illum_rgb_define[] = {
    { 0.95047, 1.0, 1.08883 },
    { 0.95047, 1.0, 1.08883 },
    { 0.96422, 1.0, 0.82521 },
    { 0.96422, 1.0, 0.82521 },
    { 0.95047, 1.0, 1.08883 },
    { 1.00000, 1.0, 1.00000 },
    { 0.96422, 1.0, 0.82521 },
    { 0.96422, 1.0, 0.82521 },
    { 0.96422, 1.0, 0.82521 },
    { 0.96422, 1.0, 0.82521 },
    { 0.98074, 1.0, 1.18232 },
    { 0.95047, 1.0, 1.08883 },
    { 0.96422, 1.0, 0.82521 },
    { 0.95047, 1.0, 1.08883 },
    { 0.95047, 1.0, 1.08883 },
    { 0.96422, 1.0, 0.82521 }
};




static yy_matrix_t matrix_adapt_ma[] = {
    {
        0.8951,-0.7502, 0.0389,
        0.2664, 1.7135,-0.0685,
       -0.1614, 0.0367, 1.0296
    },
    {
        0.40024,-0.22630, 0.00000,
        0.70760, 1.16532, 0.00000,
       -0.08081, 0.04570, 0.91822
    },
    {
        1,0,0,
        0,1,0,
        0,0,1
    },
    {
        1,0,0,
        0,1,0,
        0,0,1
    }
};






//////////////////////////////////////////////////////////////////////////////// CHEXYZ <-> CIELab

void CIEXYZ2CIELab(yy_illuminant illum,
                   CGFloat X,CGFloat Y,CGFloat Z,
                   CGFloat *L,CGFloat *a,CGFloat *b) {
    CGFloat xr = X / illum_define[illum].x;
    CGFloat yr = Y / illum_define[illum].y;
    CGFloat zr = Z / illum_define[illum].z;
    
    static CGFloat m_kE = 216.0 / 24389.0;
    static CGFloat m_kK = 24389.0 / 27.0;
    
    CGFloat fx = (xr > m_kE)
               ? pow(xr, 1.0 / 3.0)
               : ((m_kK * xr + 16.0) / 116.0);
    CGFloat fy = (yr > m_kE)
               ? pow(yr, 1.0 / 3.0)
               : ((m_kK * yr + 16.0) / 116.0);
    CGFloat fz = (zr > m_kE)
               ? pow(zr, 1.0 / 3.0)
               : ((m_kK * zr + 16.0) / 116.0);
    
    *L = 116.0 * fy - 16.0;
    *a = 500.0 * (fx - fy);
    *b = 200.0 * (fy - fz);
}


void CIELab2CIEXYZ(yy_illuminant illum,
                   CGFloat L,CGFloat a,CGFloat b,
                   CGFloat *X,CGFloat *Y,CGFloat *Z) {
    CGFloat fy = (L + 16.0) / 116.0;
    CGFloat fx = 0.002 * a + fy;
    CGFloat fz = fy - 0.005 * b;
    
    CGFloat fx3 = fx * fx * fx;
    CGFloat fz3 = fz * fz * fz;
    
    static CGFloat m_kE = 216.0 / 24389.0;
    static CGFloat m_kK = 24389.0 / 27.0;
    static CGFloat m_kKE = 8.0;
    
    CGFloat xr = (fx3 > m_kE) ? fx3 : ((116.0 * fx - 16.0) / m_kK);
    CGFloat yr = (L > m_kKE) ? pow((L + 16.0) / 116.0, 3.0) : (L / m_kK);
    CGFloat zr = (fz3 > m_kE) ? fz3 : ((116.0 * fz - 16.0) / m_kK);
    
    *X = xr * illum_define[illum].x;
    *Y = yr * illum_define[illum].y;
    *Z = zr * illum_define[illum].z;
}








//////////////////////////////////////////////////////////////////////////////// CIEXYZ <-> CIELuv

void CIEXYZ2CIELuv(yy_illuminant illum,
                   CGFloat X,CGFloat Y,CGFloat Z,
                   CGFloat *L,CGFloat *u,CGFloat *v) {
    CGFloat Den = X + 15.0 * Y + 3.0 * Z;
    CGFloat up = (Den > 0.0) ? ((4.0 * X) / (X + 15.0 * Y + 3.0 * Z)) : 0.0;
    CGFloat vp = (Den > 0.0) ? ((9.0 * Y) / (X + 15.0 * Y + 3.0 * Z)) : 0.0;
    
    yy_triple_t triple = illum_define[illum];
    CGFloat urp = (4.0 * triple.x)
        / (triple.x + 15.0 * triple.y + 3.0 * triple.z);
    CGFloat vrp = (9.0 * triple.y)
        / (triple.x + 15.0 * triple.y + 3.0 * triple.z);
    
    CGFloat yr = Y / triple.y;
    
    static CGFloat m_kE = 216.0 / 24389.0;
    static CGFloat m_kK = 24389.0 / 27.0;
    *L = (yr > m_kE) ? (116.0 * pow(yr, 1.0 / 3.0) - 16.0) : (m_kK * yr);
    *u = 13.0 * *L * (up - urp);
    *v = 13.0 * *L * (vp - vrp);
}


void CIELuv2CIEXYZ(yy_illuminant illum,
                   CGFloat L,CGFloat u,CGFloat v,
                   CGFloat *X,CGFloat *Y,CGFloat *Z) {
    static CGFloat m_kK = 24389.0 / 27.0;
    static CGFloat m_kKE = 8.0;
    *Y = (L > m_kKE) ? pow((L + 16.0) / 116.0, 3.0) : (L / m_kK);
    yy_triple_t triple = illum_define[illum];
    CGFloat u0 = (4.0 * triple.x)
        / (triple.x + 15.0 * triple.y + 3.0 * triple.z);
    CGFloat v0 = (9.0 * triple.y)
        / (triple.x + 15.0 * triple.y + 3.0 * triple.z);
    
    CGFloat a = (((52.0 * L) / (u + 13.0 * L * u0)) - 1.0) / 3.0;
    CGFloat b = -5.0 * *Y;
    CGFloat c = -1.0 / 3.0;
    CGFloat d = *Y * (((39.0 * L) / (v + 13.0 * L * v0)) - 5.0);
    
    *X = (d - b) / (a - c);
    *Z = *X * a + b;
}






////////////////////////////////////////////////////////////////////////////////Lab <-> LCHab

void CIELab2CIELCHab(CGFloat L,CGFloat a,CGFloat b,
                     CGFloat *LL,CGFloat *CC,CGFloat *HH) {
    *LL = L;
    *CC = sqrt(a * a + b * b);
    *HH = 180.0 * atan2(b, a) / M_PI;
    if (*HH < 0.0) {
        *HH += 360.0;
    }
}

void CIELCHab2CIELab(CGFloat L,CGFloat C,CGFloat H,
                     CGFloat *LL,CGFloat *aa,CGFloat *bb) {
    *LL = L;
    *aa = C * cos(H * M_PI / 180.0);
    *bb = C * sin(H * M_PI / 180.0);
}

////////////////////////////////////////////////////////////////////////////////Luv <-> LCHuv
// same as Lab <-> LCHab
void CIELuv2CIELCHuv(CGFloat L,CGFloat u,CGFloat v,
                     CGFloat *LL,CGFloat *CC,CGFloat *HH) {
    *LL = L;
    *CC = sqrt(u * u + v * v);
    *HH = 180.0 * atan2(v, u) / M_PI;
    if (*HH < 0.0) {
        *HH += 360.0;
    }
}

void CIELCHuv2CIELuv(CGFloat L,CGFloat C,CGFloat H,
                     CGFloat *LL,CGFloat *aa,CGFloat *bb) {
    *LL = L;
    *aa = C * cos(H * M_PI / 180.0);
    *bb = C * sin(H * M_PI / 180.0);
}



////////////////////////////////////////////////////////////////////////////////Util




/// |M|
static inline CGFloat Determinant3x3(const yy_matrix_t m) {
    CGFloat det = m.m00 * (m.m22 * m.m11 - m.m21 * m.m12)
                - m.m10 * (m.m22 * m.m01 - m.m21 * m.m02)
                + m.m20 * (m.m12 * m.m01 - m.m11 * m.m02);
    return det;
}

/// AB = BA
static inline void MatrixInvert3x3(const yy_matrix_t m, yy_matrix_t *i) {
    CGFloat scale = 1.0 / Determinant3x3(m);
    
    i->m00 =  scale * (m.m22 * m.m11 - m.m21 * m.m12);
    i->m01 = -scale * (m.m22 * m.m01 - m.m21 * m.m02);
    i->m02 =  scale * (m.m12 * m.m01 - m.m11 * m.m02);
    
    i->m10 = -scale * (m.m22 * m.m10 - m.m20 * m.m12);
    i->m11 =  scale * (m.m22 * m.m00 - m.m20 * m.m02);
    i->m12 = -scale * (m.m12 * m.m00 - m.m10 * m.m02);
    
    i->m20 =  scale * (m.m21 * m.m10 - m.m20 * m.m11);
    i->m21 = -scale * (m.m21 * m.m00 - m.m20 * m.m01);
    i->m22 =  scale * (m.m11 * m.m00 - m.m10 * m.m01);
}


static inline void MatrixTranspose3x3(yy_matrix_t *m) {
    CGFloat v;
    v = m->m01;
    m->m01 = m->m10;
    m->m10 = v;
    
    v = m->m02;
    m->m02 = m->m20;
    m->m20 = v;
    
    v = m->m12;
    m->m12 = m->m21;
    m->m21 = v;
}


////////////////////////////////////////////////////////////////////////////////RGB XYZ


static inline void getMatrixRGB2XYZ(yy_rgbspace rgbspace,
                                    yy_matrix_t *matrixRGB2XYZ){
    CGFloat xr, yr, xg, yg, xb, yb;
    
    yy_triple_t t_RefWhiteRGB = illum_rgb_define[rgbspace];
    xr = xyrgb_define[rgbspace].xr;
    yr = xyrgb_define[rgbspace].yr;
    xg = xyrgb_define[rgbspace].xg;
    yg = xyrgb_define[rgbspace].yg;
    xb = xyrgb_define[rgbspace].xb;
    yb = xyrgb_define[rgbspace].yb;
    
    yy_matrix_t m = {
        .m00 = xr/yr,        .m01 = xg/yg,        .m02 = xb/yb,
        .m10 = 1.0,          .m11 = 1.0,          .m12 = 1.0,
        .m20=(1.0-xr-yr)/yr, .m21=(1.0-xg-yg)/yg, .m22 = (1.0-xb-yb)/yb
    };
    
    yy_matrix_t mi = {1,0,0, 0,1,0, 0,0,1};
    
    MatrixInvert3x3(m, &mi);
    
    CGFloat sr = t_RefWhiteRGB.x * mi.m00
        + t_RefWhiteRGB.y * mi.m01 + t_RefWhiteRGB.z * mi.m02;
    CGFloat sg = t_RefWhiteRGB.x * mi.m10
        + t_RefWhiteRGB.y * mi.m11 + t_RefWhiteRGB.z * mi.m12;
    CGFloat sb = t_RefWhiteRGB.x * mi.m20
        + t_RefWhiteRGB.y * mi.m21 + t_RefWhiteRGB.z * mi.m22;
    
    matrixRGB2XYZ->m00 = sr * m.m00;
    matrixRGB2XYZ->m01 = sg * m.m01;
    matrixRGB2XYZ->m02 = sb * m.m02;
    matrixRGB2XYZ->m10 = sr * m.m10;
    matrixRGB2XYZ->m11 = sg * m.m11;
    matrixRGB2XYZ->m12 = sb * m.m12;
    matrixRGB2XYZ->m20 = sr * m.m20;
    matrixRGB2XYZ->m21 = sg * m.m21;
    matrixRGB2XYZ->m22 = sb * m.m22;
    
    MatrixTranspose3x3(matrixRGB2XYZ);
}



static inline void getMatrixXYZ2RGB(yy_rgbspace rgbspace,
                                      yy_matrix_t *matrixXYZ2RGB){
    yy_matrix_t matrixRGB2XYZ;
    getMatrixRGB2XYZ(rgbspace, &matrixRGB2XYZ);
    MatrixInvert3x3(matrixRGB2XYZ, matrixXYZ2RGB);
}





////////////////////////////////////////////////////////////////////////////////Gamma



static inline CGFloat EncodeGammasRGBFast(CGFloat value) {
    if (value <= (0.0031306684425005883)) return(12.92f * value);
    
    div_t quotient;
    CGFloat p, term[9];
    int exponent;
    
    static const CGFloat coeff[] = { // Chebychevi poly: x^(5/12), x=1.5
        1.1758200232996901923,
        0.16665763094889061230,
       -0.0083154894939042125035,
        0.00075187976780420279038,
       -0.000083240178519391795367,
        0.000010229209410070008679,
       -1.3400466409860246e-06,
        1.8333422241635376682e-07,
       -2.5878596761348859722e-08
    };
    
    static const CGFloat powers_of_two[] = { // (2^N)^(5/12)
        1.0,
        1.3348398541700343678,
        1.7817974362806785482,
        2.3784142300054420538,
        3.1748021039363991669,
        4.2378523774371812394,
        5.6568542494923805819,
        7.5509945014535482244,
        1.0079368399158985525e1,
        1.3454342644059433809e1,
        1.7959392772949968275e1,
        2.3972913230026907883e1
    };
    
    // Compute x^(1/2.4) == x^(5/12) == pow(x,1.0/2.4).
    term[0] = 1.0;
    term[1] = 4.0 * frexp(value, &exponent) - 3.0;
    term[2] = 2.0 * term[1] * term[1] - term[0];
    term[3] = 2.0 * term[1] * term[2] - term[1];
    term[4] = 2.0 * term[1] * term[3] - term[2];
    term[5] = 2.0 * term[1] * term[4] - term[3];
    term[6] = 2.0 * term[1] * term[5] - term[4];
    term[7] = 2.0 * term[1] * term[6] - term[5];
    term[8] = 2.0 * term[1] * term[7] - term[6];
    p = coeff[0] * term[0] + coeff[1] * term[1] + coeff[2] * term[2]
      + coeff[3] * term[3] + coeff[4] * term[4] + coeff[5] * term[5]
      + coeff[6] * term[6] + coeff[7] * term[7] + coeff[8] * term[8];
    quotient = div(exponent - 1, 12);
    if (quotient.rem < 0) {
        quotient.quot -= 1;
        quotient.rem += 12;
    }
    value = (ldexp(powers_of_two[quotient.rem] * p, 5 * quotient.quot));
    
    return (1.055 * value - 0.055);
}


static inline CGFloat DecodeGammasRGBFast(CGFloat value) {
    if (value <= (0.0404482362771076)) return(value / 12.92f);
    value = (value + 0.055) / 1.055;
    
    div_t quotient;
    CGFloat p, term[9];
    int exponent;
    
    static const CGFloat coeff[] = { /* terms for x^(7/5), x=1.5 */
        1.7917488588043277509,
        0.82045614371976854984,
        0.027694100686325412819,
       -0.00094244335181762134018,
        0.000064355540911469709545,
       -5.7224404636060757485e-06,
        5.8767669437311184313e-07,
       -6.6139920053589721168e-08,
        7.9323242696227458163e-09
    };
    
    static const CGFloat powers_of_two[] = { /* (2^x)^(7/5) */
        1.0,
        2.6390158215457883983,
        6.9644045063689921093,
        1.8379173679952558018e+01,
        4.8502930128332728543e+01
    };
    
    /*
     * Compute x^2.4 == x*x^(7/5) == pow(x,2.4).
     */
    term[0] = 1.0;
    term[1] = 4.0 * frexp(value, &exponent) - 3.0;
    term[2] = 2.0 * term[1] * term[1] - term[0];
    term[3] = 2.0 * term[1] * term[2] - term[1];
    term[4] = 2.0 * term[1] * term[3] - term[2];
    term[5] = 2.0 * term[1] * term[4] - term[3];
    term[6] = 2.0 * term[1] * term[5] - term[4];
    term[7] = 2.0 * term[1] * term[6] - term[5];
    term[8] = 2.0 * term[1] * term[7] - term[6];
    p = coeff[0] * term[0] + coeff[1] * term[1] + coeff[2] * term[2]
      + coeff[3] * term[3] + coeff[4] * term[4] + coeff[5] * term[5]
      + coeff[6] * term[6] + coeff[7] * term[7] + coeff[8] * term[8];
    quotient = div(exponent - 1, 5);
    if (quotient.rem < 0) {
        quotient.quot -= 1;
        quotient.rem += 5;
    }
    return (value * ldexp(powers_of_two[quotient.rem] * p, 7 * quotient.quot));
}


static inline CGFloat EncodeGammasRGB(CGFloat value) {
    CGFloat sign = 1.0;
    if (value < 0.0) {
        sign = -1.0;
        value = -value;
    }
    value = value <= 0.0031308
          ? (value * 12.92)
          : (1.055 * pow(value, 1.0 / 2.4) - 0.055);
    value *= sign;
    return value;
}


static inline CGFloat DecodeGammasRGB(CGFloat value) {
    CGFloat sign = 1.0;
    if (value < 0.0) {
        sign = -1.0;
        value = -value;
    }
    value = value <= 0.04045
          ? (value / 12.92)
          : pow((value + 0.055) / 1.055, 2.4);
    value *= sign;
    return value;
}



///gamma > 0 : normal, gamma < 0 : sRGB,  0 : L*
static inline CGFloat EncodeGamma(CGFloat value, CGFloat gamma) {
    if (gamma > 0.0) {
        return (value >= 0.0)
              ? pow(value, 1.0 / gamma)
              : -pow(-value, 1.0 / gamma);
    } else if (gamma < 0.0) {
        return EncodeGammasRGB(value);
    } else {
        /* L* */
        CGFloat sign = 1.0;
        if (value < 0.0) {
            sign = -1.0;
            value = -value;
        }
        value = (value <= (216.0 / 24389.0))
              ? (value * 24389.0 / 2700.0)
              : (1.16 * pow(value, 1.0 / 3.0) - 0.16);
        value *= sign;
        return value;
    }
}


///gamma > 0 : normal, gamma < 0 : sRGB,  0 : L*
static inline CGFloat DecodeGamma(CGFloat value, CGFloat gamma) {
    if (gamma > 0) {
        return (value >= 0.0) ? pow(value, gamma) : -pow(-value, gamma);
    } else if (gamma < 0) {
        return DecodeGammasRGB(value);
    } else {
        CGFloat linear;
        CGFloat sign = 1.0;
        if (value < 0.0) {
            sign = -1.0;
            value = -value;
        }
        linear = (value <= 0.08)
               ? (2700.0 * value / 24389.0)
               : ((((1000000.0 * value + 480000.0) * value + 76800.0)
                 * value + 4096.0) / 1560896.0);
        linear *= sign;
        return linear;
    }
}










// Observer = 2°, Illuminant = D65, sRGB sRGB Adapt=NONE
// X:[0,0.95047], Y:[0,1.0], Z:[0,1.08883]

void RGB2CIEXYZDefault(CGFloat R, CGFloat G, CGFloat B,
                       CGFloat *X, CGFloat *Y,CGFloat *Z) {
    R = DecodeGammasRGBFast(R);
    G = DecodeGammasRGBFast(G);
    B = DecodeGammasRGBFast(B);
    *X = 0.41239558896741421610 * R + 0.35758343076371481710 * G
       + 0.18049264738170157350 * B;
    *Y = 0.21258623078559555160 * R + 0.71517030370341084990 * G
       + 0.07220049864333622685 * B;
    *Z = 0.01929721549174694484 * R + 0.11918386458084853180 * G
       + 0.95049712513157976600 * B;
}

// Observer = 2°, Illuminant = D65, sRGB sRGB Adapt=NONE
void CIEXYZ2RGBDefault(CGFloat X, CGFloat Y, CGFloat Z,
                       CGFloat *R, CGFloat *G, CGFloat *B) {
    *R =  3.2406 * X - 1.5372 * Y - 0.4986 * Z;
    *G = -0.9689 * X + 1.8758 * Y + 0.0415 * Z;
    *B =  0.0557 * X - 0.2040 * Y + 1.0570 * Z;
    *R = EncodeGammasRGBFast(*R);
    *G = EncodeGammasRGBFast(*G);
    *B = EncodeGammasRGBFast(*B);
}







void RGB2CIEXYZ(yy_illuminant illum, yy_gamma gammaEnum,
                yy_rgbspace rgbspace, yy_color_adaptation adapt,
                CGFloat R, CGFloat G, CGFloat B,
                CGFloat *X, CGFloat *Y, CGFloat *Z) {
    
    float gamma = gamma_define[gammaEnum];
    R = DecodeGamma(R, gamma);
    G = DecodeGamma(G, gamma);
    B = DecodeGamma(B, gamma);
    
    yy_matrix_t matrixRGB2XYZ;
    getMatrixRGB2XYZ(rgbspace, &matrixRGB2XYZ);
    
    /// XYZ = |RGB| x |MatrixRGB2XYZ|
    *X = R * matrixRGB2XYZ.m00 + G * matrixRGB2XYZ.m10 + B * matrixRGB2XYZ.m20;
    *Y = R * matrixRGB2XYZ.m01 + G * matrixRGB2XYZ.m11 + B * matrixRGB2XYZ.m21;
    *Z = R * matrixRGB2XYZ.m02 + G * matrixRGB2XYZ.m12 + B * matrixRGB2XYZ.m22;
    

    if (adapt != yy_adaption_NONE) {
        yy_triple_t t_illum = illum_define[illum];
        yy_triple_t t_illum_rgb = illum_rgb_define[rgbspace];
        yy_matrix_t m_adapt_ma = matrix_adapt_ma[adapt];
        yy_matrix_t m_adapt_mai;
        MatrixInvert3x3(m_adapt_ma, &m_adapt_mai);
        
        CGFloat Ad = t_illum.x * m_adapt_ma.m00
                   + t_illum.y * m_adapt_ma.m10
                   + t_illum.z * m_adapt_ma.m20;
        
        CGFloat Bd = t_illum.x * m_adapt_ma.m01
                   + t_illum.y * m_adapt_ma.m11
                   + t_illum.z * m_adapt_ma.m21;
        
        CGFloat Cd = t_illum.x * m_adapt_ma.m02
                   + t_illum.y * m_adapt_ma.m12
                   + t_illum.z * m_adapt_ma.m22;
        
        CGFloat As = t_illum_rgb.x * m_adapt_ma.m00
                   + t_illum_rgb.y * m_adapt_ma.m10
                   + t_illum_rgb.z * m_adapt_ma.m20;
        
        CGFloat Bs = t_illum_rgb.x * m_adapt_ma.m01
                   + t_illum_rgb.y * m_adapt_ma.m11
                   + t_illum_rgb.z * m_adapt_ma.m21;
        
        CGFloat Cs = t_illum_rgb.x * m_adapt_ma.m02
                   + t_illum_rgb.y * m_adapt_ma.m12
                   + t_illum_rgb.z * m_adapt_ma.m22;
        
        CGFloat XX = *X * m_adapt_ma.m00
                   + *Y * m_adapt_ma.m10
                   + *Z * m_adapt_ma.m20;
        
        CGFloat YY = *X * m_adapt_ma.m01
                   + *Y * m_adapt_ma.m11
                   + *Z * m_adapt_ma.m21;
        
        CGFloat ZZ = *X * m_adapt_ma.m02
                   + *Y * m_adapt_ma.m12
                   + *Z * m_adapt_ma.m22;
        
        XX *= (Ad / As);
        YY *= (Bd / Bs);
        ZZ *= (Cd / Cs);
        
        *X = XX * m_adapt_mai.m00
           + YY * m_adapt_mai.m10 + ZZ * m_adapt_mai.m20;
        *Y = XX * m_adapt_mai.m01
           + YY * m_adapt_mai.m11 + ZZ * m_adapt_mai.m21;
        *Z = XX * m_adapt_mai.m02
           + YY * m_adapt_mai.m12 + ZZ * m_adapt_mai.m22;
    }
}






void CIEXYZ2RGB(yy_illuminant illum, yy_gamma gammaEnum,
                yy_rgbspace colorspace, yy_color_adaptation adapt,
                CGFloat X, CGFloat Y, CGFloat Z,
                CGFloat *R, CGFloat *G, CGFloat *B) {
    CGFloat X2 = X;
    CGFloat Y2 = Y;
    CGFloat Z2 = Z;
    
    if (adapt != yy_adaption_NONE) {
        yy_triple_t t_illum = illum_define[illum];
        yy_triple_t t_illum_rgb = illum_rgb_define[colorspace];
        yy_matrix_t m_adapt_ma = matrix_adapt_ma[adapt];
        yy_matrix_t m_adapt_mai;
        MatrixInvert3x3(m_adapt_ma, &m_adapt_mai);
        
        CGFloat As = t_illum.x * m_adapt_ma.m00
                   + t_illum.y * m_adapt_ma.m10
                   + t_illum.z * m_adapt_ma.m20;
        
        CGFloat Bs = t_illum.x * m_adapt_ma.m01
                   + t_illum.y * m_adapt_ma.m11
                   + t_illum.z * m_adapt_ma.m21;
        
        CGFloat Cs = t_illum.x * m_adapt_ma.m02
                   + t_illum.y * m_adapt_ma.m12
                   + t_illum.z * m_adapt_ma.m22;
        
        CGFloat Ad = t_illum_rgb.x * m_adapt_ma.m00
                   + t_illum_rgb.y * m_adapt_ma.m10
                   + t_illum_rgb.z * m_adapt_ma.m20;
        
        CGFloat Bd = t_illum_rgb.x * m_adapt_ma.m01
                   + t_illum_rgb.y * m_adapt_ma.m11
                   + t_illum_rgb.z * m_adapt_ma.m21;
        
        CGFloat Cd = t_illum_rgb.x * m_adapt_ma.m02
                   + t_illum_rgb.y * m_adapt_ma.m12
                   + t_illum_rgb.z * m_adapt_ma.m22;
        
        CGFloat X1 = X * m_adapt_ma.m00
                   + Y * m_adapt_ma.m10
                   + Z * m_adapt_ma.m20;
        
        CGFloat Y1 = X * m_adapt_ma.m01
                   + Y * m_adapt_ma.m11
                   + Z * m_adapt_ma.m21;
        
        CGFloat Z1 = X * m_adapt_ma.m02
                   + Y * m_adapt_ma.m12
                   + Z * m_adapt_ma.m22;
        
        X1 *= (Ad / As);
        Y1 *= (Bd / Bs);
        Z1 *= (Cd / Cs);
        
        X2 = X1 * m_adapt_mai.m00
           + Y1 * m_adapt_mai.m10
           + Z1 * m_adapt_mai.m20;
        
        Y2 = X1 * m_adapt_mai.m01
           + Y1 * m_adapt_mai.m11
           + Z1 * m_adapt_mai.m21;
        
        Z2 = X1 * m_adapt_mai.m02
           + Y1 * m_adapt_mai.m12
           + Z1 * m_adapt_mai.m22;
    }
    
    yy_matrix_t t_MatrixXYZ2RGB;
    getMatrixXYZ2RGB(colorspace, &t_MatrixXYZ2RGB);
    float gamma = gamma_define[gammaEnum];
    *R = EncodeGamma(X2 * t_MatrixXYZ2RGB.m00
                   + Y2 * t_MatrixXYZ2RGB.m10
                   + Z2 * t_MatrixXYZ2RGB.m20, gamma);
    *G = EncodeGamma(X2 * t_MatrixXYZ2RGB.m01
                   + Y2 * t_MatrixXYZ2RGB.m11
                   + Z2 * t_MatrixXYZ2RGB.m21, gamma);
    *B = EncodeGamma(X2 * t_MatrixXYZ2RGB.m02
                   + Y2 * t_MatrixXYZ2RGB.m12
                   + Z2 * t_MatrixXYZ2RGB.m22, gamma);
}

































//////////////////////////////////////////////////////////////////////////////// XYZ <-> LMS

void CIEXYZ2LMS(CGFloat x, CGFloat y, CGFloat z,
                CGFloat *L, CGFloat *M, CGFloat *S) {
    *L =  0.7328 * x + 0.4296 * y - 0.1624 * z;
    *M = -0.7036 * x + 1.6975 * y + 0.0061 * z;
    *S =  0.0030 * x + 0.0136 * y + 0.9834 * z;
}

void CIELMSToXYZ(CGFloat L, CGFloat M, CGFloat S,
                 CGFloat *X, CGFloat *Y, CGFloat *Z) {
    *X = 1.096123820835514 * L
       - 0.278869000218287 * M + 0.182745179382773 * S;
    *Y = 0.454369041975359 * L
       + 0.473533154307412 * M + 0.072097803717229 * S;
    *Z =-0.009627608738429 * L
       - 0.005698031216113 * M + 1.015325639954543 * S;
}



//////////////////////////////////////////////////////////////////////////////// RGB <-> YUV


///RGB [0,1] Y[0,1] U[-0.436,0.436] V[-0.615,0.615]
void RGB2YUV(CGFloat R, CGFloat G, CGFloat B,
             CGFloat *Y, CGFloat *U, CGFloat *V) {
    CLAMP_COLOR_VALUE(R);
    CLAMP_COLOR_VALUE(G);
    CLAMP_COLOR_VALUE(B);
    *Y =  0.298839 * R + 0.586811 * G + 0.114350 * B;
    *U = -0.147    * R - 0.289    * G + 0.436    * B + 0.5;
    *V =  0.615    * R - 0.515    * G - 0.100    * B + 0.5;
    
}

///RGB [0,1] Y[0,1] U[-0.436,0.436] V[-0.615,0.615]
void YUV2RGB(CGFloat Y, CGFloat U, CGFloat V,
             CGFloat *R, CGFloat *G, CGFloat *B) {
    U -= 0.5;
    V -= 0.5;
    *R = Y - 3.945707070708279e-05 * U + 1.1398279671717170825 * V;
    *G = Y - 0.3946101641414141437 * U - 0.5805003156565656797 * V;
    *B = Y + 2.0319996843434342537 * U - 4.813762626262513e-04 * V;
    CLAMP_COLOR_VALUE(*R);
    CLAMP_COLOR_VALUE(*G);
    CLAMP_COLOR_VALUE(*B);
}



//////////////////////////////////////////////////////////////////////////////// RGB <-> YCbCr

void RGB2YCbCr(CGFloat R, CGFloat G, CGFloat B,
               CGFloat *Y, CGFloat *Cb, CGFloat *Cr) {
    CLAMP_COLOR_VALUE(R);
    CLAMP_COLOR_VALUE(G);
    CLAMP_COLOR_VALUE(B);
    *Y  =  0.298839 * R + 0.586811 * G + 0.114350 * B;
    *Cb = -0.168737 * R - 0.331264 * G + 0.500000 * B + 0.5;
    *Cr =  0.500000 * R - 0.418688 * G - 0.081312 * B + 0.5;
}

void YCbCr2RGB(CGFloat Y, CGFloat Cb, CGFloat Cr,
               CGFloat *R, CGFloat *G, CGFloat *B) {
    Cb -= 0.5;
    Cr -= 0.5;
    *R = 0.99999999999914679361 * Y
       - 1.2188941887145875e-06 * Cb + 1.4019995886561440468 * Cr;
    *G = 0.99999975910502514331 * Y
       - 0.34413567816504303521 * Cb - 0.71413649331646789076 * Cr;
    *B = 1.00000124040004623180 * Y
       + 1.77200006607230409200 * Cb + 2.1453384174593273e-06 * Cr;
    CLAMP_COLOR_VALUE(*R);
    CLAMP_COLOR_VALUE(*G);
    CLAMP_COLOR_VALUE(*B);
}

//////////////////////////////////////////////////////////////////////////////// RGB <-> YPbPr

void RGB2YPbPr(CGFloat R, CGFloat G, CGFloat B,
               CGFloat *Y, CGFloat *Pb, CGFloat *Pr) {
    CLAMP_COLOR_VALUE(R);
    CLAMP_COLOR_VALUE(G);
    CLAMP_COLOR_VALUE(B);
    *Y  =  0.298839 * R + 0.586811 * G + 0.114350 * B;
    *Pb = -0.168737 * R - 0.331264 * G + 0.500000 * B + 0.5;
    *Pr =  0.500000 * R - 0.418688 * G - 0.081312 * B + 0.5;
}

void YPbPr2RGB(CGFloat Y, CGFloat Pb, CGFloat Pr,
               CGFloat *R, CGFloat *G, CGFloat *B) {
    Pb -= 0.5;
    Pr -= 0.5;
    *R = 0.99999999999914679361 * Y
       - 1.2188941887145875e-06 * Pb + 1.4019995886561440468 * Pr;
    *G = 0.99999975910502514331 * Y
       - 0.34413567816504303521 * Pb - 0.71413649331646789076 * Pr;
    *B = 1.00000124040004623180 * Y
       + 1.77200006607230409200 * Pb + 2.1453384174593273e-06 * Pr;
    CLAMP_COLOR_VALUE(*R);
    CLAMP_COLOR_VALUE(*G);
    CLAMP_COLOR_VALUE(*B);
}

//////////////////////////////////////////////////////////////////////////////// RGB <-> YDbDr

void RGB2YDbDr(CGFloat R, CGFloat G, CGFloat B,
               CGFloat *Y, CGFloat *Db, CGFloat *Dr) {
    CLAMP_COLOR_VALUE(R);
    CLAMP_COLOR_VALUE(G);
    CLAMP_COLOR_VALUE(B);
    *Y  =  0.298839 * R + 0.586811 * G + 0.114350 * B;
    *Db = -0.450    * R - 0.883    * G + 1.333    * B + 0.5;
    *Dr = -1.333    * R + 1.116    * G + 0.217    * B + 0.5;
}

void YDbDr2RGB(CGFloat Y, CGFloat Db, CGFloat Dr,
               CGFloat *R, CGFloat *G, CGFloat *B) {
    Dr -= 0.5;
    Db -= 0.5;
    *R =  Y + 9.2303716147657e-05 * Db - 0.52591263066186533 * Dr;
    *G =  Y - 0.12913289889050927 * Db + 0.26789932820759876 * Dr;
    *B =  Y + 0.66467905997895482 * Db - 7.9202543533108e-05 * Dr;
    CLAMP_COLOR_VALUE(*R);
    CLAMP_COLOR_VALUE(*G);
    CLAMP_COLOR_VALUE(*B);
}

//////////////////////////////////////////////////////////////////////////////// RGB <-> YIQ

void RGB2YIQ(CGFloat R, CGFloat G, CGFloat B,
             CGFloat *Y, CGFloat *I, CGFloat *Q) {
    CLAMP_COLOR_VALUE(R);
    CLAMP_COLOR_VALUE(G);
    CLAMP_COLOR_VALUE(B);
    *Y = 0.298839 * R + 0.586811 * G + 0.114350 * B;
    *I = 0.595716 * R - 0.274453 * G - 0.321263 * B + 0.5;
    *Q = 0.211456 * R - 0.522591 * G + 0.311135 * B + 0.5;
}

void YIQ2RGB(CGFloat Y, CGFloat I, CGFloat Q,
             CGFloat *R, CGFloat *G, CGFloat *B) {
    I -= 0.5;
    Q -= 0.5;
    *R = Y + 0.9562957197589482261 * I + 0.6210244164652610754 * Q;
    *G = Y - 0.2721220993185104464 * I - 0.6473805968256950427 * Q;
    *B = Y - 1.1069890167364901945 * I + 1.7046149983646481374 * Q;
    CLAMP_COLOR_VALUE(*R);
    CLAMP_COLOR_VALUE(*G);
    CLAMP_COLOR_VALUE(*B);
}






