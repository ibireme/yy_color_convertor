//
//  yy_color_converter.h
//
//  Created by ibireme on 13-8-17.
//

#ifndef yy_color_h
#define yy_color_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef CGFLOAT_DEFINED
#define CGFLOAT_DEFINED 1
typedef double   CGFloat;
#endif





typedef enum  {
    yy_illuminant_A = 0,
    yy_illuminant_B,
    yy_illuminant_C,
    yy_illuminant_D50,  ///default
    yy_illuminant_D55,
    yy_illuminant_D65,
    yy_illuminant_D75,
    yy_illuminant_E,
    yy_illuminant_F2,
    yy_illuminant_F7,
    yy_illuminant_F11
}yy_illuminant;


typedef enum {
    yy_gamma_1_0 = 0,
    yy_gamma_1_8 = 1,
    yy_gamma_2_2 = 2,
    yy_gamma_sRGB = 3,   ///default
    yy_gamma_L = 4
} yy_gamma;


typedef enum {
    yy_rgbspace_AdobeRGB = 0,
    yy_rgbspace_AppleRGB,
    yy_rgbspace_BestRGB,
    yy_rgbspace_BetaRGB,
    yy_rgbspace_BruceRGB,
    yy_rgbspace_CIE_RGB,
    yy_rgbspace_ColorMatchRGB,
    yy_rgbspace_DonRGB4,
    yy_rgbspace_ECI_RGBv2,
    yy_rgbspace_EktaSpacePS5,
    yy_rgbspace_NTSC_RGB,
    yy_rgbspace_PAL_RGB,
    yy_rgbspace_ProPhotoRGB,
    yy_rgbspace_SMPTE_C_RGB,
    yy_rgbspace_sRGB,       ///default
    yy_rgbspace_WideGamutRGB
} yy_rgbspace;



typedef enum {
    yy_adaption_Bradform = 0,
    yy_adaption_VonKries = 1,
    yy_adaption_XYZ_Saling = 2,
    yy_adaption_NONE = 3     ///default
} yy_color_adaptation;



/**
 * RGB <-> HSL
 * RGB <-> HSB
 * RGB <-> CMYK
 * RGB <-> CMY
 * RGB <-> YUV
 * RGB <-> YCbCr
 * RGB <-> YDbDr
 * RGB <-> YPbPr
 * RGB <-> YIQ
 * CMYK <-> CMY
 * HSL <-> HSB
 *
 * RGB <-> CIEXYZ
 *
 * CIEXYZ <-> CIExyY
 * CIEXYZ <-> HunterLab
 * CIEXYZ <-> LMS
 * CIEXYZ <-> CIELab
 * CIEXYZ <-> CIELuv
 * CIELab <-> CIELCHab
 * CIELuv <-> CIELCHuv
 *
 * RGB,HSL,HSB,CMYK,YUV,YIQ,XYZ,xyY,Lab,Luv,LCH ...
 * mostly,input and output value is in range 0 to 1
 */



extern void RGB2HSL(CGFloat r, CGFloat g, CGFloat b,
                    CGFloat *h, CGFloat *s, CGFloat *l);

extern void HSL2RGB(CGFloat h, CGFloat s, CGFloat l,
                    CGFloat *r, CGFloat *g, CGFloat *b);



extern void RGB2HSB(CGFloat r, CGFloat g, CGFloat b,
                    CGFloat *h, CGFloat *s, CGFloat *v);

extern void HSB2RGB(CGFloat h, CGFloat s, CGFloat v,
                    CGFloat *r, CGFloat *g, CGFloat *b);



extern void HSB2HSL(CGFloat h, CGFloat s, CGFloat b,
                    CGFloat *hh, CGFloat *ss, CGFloat *ll);

extern void HSL2HSB(CGFloat h, CGFloat s, CGFloat l,
                    CGFloat *hh, CGFloat *ss, CGFloat *bb);



extern void RGB2CMYK(CGFloat r, CGFloat g, CGFloat b,
                     CGFloat *c, CGFloat *m, CGFloat *y, CGFloat *k);

extern void CMYK2RGB(CGFloat c, CGFloat m, CGFloat y, CGFloat k,
                     CGFloat *r, CGFloat *g, CGFloat *b);



extern void RGB2CMY(CGFloat r, CGFloat g, CGFloat b,
                    CGFloat *c, CGFloat *m, CGFloat *y);

extern void CMY2RGB(CGFloat c, CGFloat m, CGFloat y,
                    CGFloat *r, CGFloat *g, CGFloat *b);



extern void CMY2CMYK(CGFloat c, CGFloat m, CGFloat y,
                     CGFloat *cc, CGFloat *mm, CGFloat *yy, CGFloat *kk);

extern void CMYK2CMY(CGFloat c, CGFloat m, CGFloat y, CGFloat k,
                     CGFloat *cc, CGFloat *mm, CGFloat *yy);



extern void RGB2YUV(CGFloat R, CGFloat G, CGFloat B,
               CGFloat *Y, CGFloat *U, CGFloat *V);

extern void YUV2RGB(CGFloat Y, CGFloat U, CGFloat V,
               CGFloat *R, CGFloat *G, CGFloat *B);



extern void RGB2YCbCr(CGFloat R, CGFloat G, CGFloat B,
               CGFloat *Y, CGFloat *Cb, CGFloat *Cr);

extern void YCbCr2RGB(CGFloat Y, CGFloat Cb, CGFloat Cr,
               CGFloat *R, CGFloat *G, CGFloat *B);



extern void RGB2YPbPr(CGFloat R, CGFloat G, CGFloat B,
               CGFloat *Y, CGFloat *Pb, CGFloat *Pr);

extern void YPbPr2RGB(CGFloat Y, CGFloat Pb, CGFloat Pr,
               CGFloat *R, CGFloat *G, CGFloat *B);



extern void RGB2YDbDr(CGFloat R, CGFloat G, CGFloat B,
               CGFloat *Y, CGFloat *Db, CGFloat *Dr);

extern void YDbDr2RGB(CGFloat Y, CGFloat Db, CGFloat Dr,
               CGFloat *R, CGFloat *G, CGFloat *B);



extern void RGB2YIQ(CGFloat R, CGFloat G, CGFloat B,
             CGFloat *Y, CGFloat *I, CGFloat *Q);

extern void YIQ2RGB(CGFloat Y, CGFloat I, CGFloat Q,
             CGFloat *R, CGFloat *G, CGFloat *B);



extern void RGB2CIEXYZDefault(CGFloat R, CGFloat G, CGFloat B,
                       CGFloat *X, CGFloat *Y, CGFloat *Z);

extern void CIEXYZ2RGBDefault(CGFloat X, CGFloat Y, CGFloat Z,
                       CGFloat *R, CGFloat *G, CGFloat *B);



extern void RGB2CIEXYZ(yy_illuminant illum, yy_gamma gammaEnum,
                yy_rgbspace rgbspace, yy_color_adaptation adapt,
                CGFloat R, CGFloat G, CGFloat B,
                CGFloat *X, CGFloat *Y, CGFloat *Z);

extern void CIEXYZ2RGB(yy_illuminant illum, yy_gamma gammaEnum,
                yy_rgbspace colorspace, yy_color_adaptation adapt,
                CGFloat X, CGFloat Y, CGFloat Z,
                CGFloat *R, CGFloat *G, CGFloat *B);



extern void CIEXYZ2CIExyY(CGFloat X, CGFloat Y, CGFloat Z,
                   CGFloat *xx, CGFloat *yy, CGFloat *YY);

extern void CIExyY2CIEXYZ(CGFloat x, CGFloat y, CGFloat Y,
                   CGFloat *XX, CGFloat *YY, CGFloat *ZZ);



extern void CIEXYZ2HunterLab(CGFloat X, CGFloat Y, CGFloat Z,
                      CGFloat *L, CGFloat *a, CGFloat *b);

extern void HunterLab2CIEXYZ(CGFloat L, CGFloat a, CGFloat b,
                      CGFloat *X, CGFloat *Y, CGFloat *Z);



extern void CIEXYZ2CIELab(yy_illuminant illum,
                   CGFloat X, CGFloat Y, CGFloat Z,
                   CGFloat *L, CGFloat *a, CGFloat *b);

extern void CIELab2CIEXYZ(yy_illuminant illum,
                   CGFloat L, CGFloat a, CGFloat b,
                   CGFloat *X, CGFloat *Y, CGFloat *Z);



extern void CIEXYZ2CIELuv(yy_illuminant illum,
                   CGFloat X, CGFloat Y, CGFloat Z,
                   CGFloat *L, CGFloat *u, CGFloat *v);

extern void CIELuv2CIEXYZ(yy_illuminant illum,
                   CGFloat L, CGFloat u, CGFloat v,
                   CGFloat *X, CGFloat *Y, CGFloat *Z);



extern void CIELab2CIELCHab(CGFloat L, CGFloat a, CGFloat b,
                     CGFloat *LL, CGFloat *CC, CGFloat *HH);

extern void CIELCHab2CIELab(CGFloat L, CGFloat C, CGFloat H,
                     CGFloat *LL, CGFloat *aa, CGFloat *bb);



extern void CIELuv2CIELCHuv(CGFloat L, CGFloat u, CGFloat v,
                     CGFloat *LL, CGFloat *CC, CGFloat *HH);

extern void CIELCHuv2CIELuv(CGFloat L, CGFloat C, CGFloat H,
                     CGFloat *LL, CGFloat *aa, CGFloat *bb);



extern void CIEXYZ2LMS(CGFloat x, CGFloat y, CGFloat z,
             CGFloat *L, CGFloat *M, CGFloat *S);

extern void LMSToCIEXYZ(CGFloat L, CGFloat M, CGFloat S,
              CGFloat *X, CGFloat *Y, CGFloat *Z);


#endif
