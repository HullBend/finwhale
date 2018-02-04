/*
 * Copyright 2013 SPZ
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
/*
 * ====================================================
 * Notice of SunSoft copyrighted software this class is derived from:
 *
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunSoft, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */
package math;


/**
 * This is a utility class that provides computation methods related to the
 * error functions.
 */
@Deprecated
public final class Erf {

    private static final double erx = 8.45062911510467529297e-01;
    private static final double efx = 1.28379167095512586316e-01;

    // Coefficients for approximation to erf in [0,0.84375]

    private static final double pp0 = 1.28379167095512558561e-01;
    private static final double pp1 = -3.25042107247001499370e-01;
    private static final double pp2 = -2.84817495755985104766e-02;
    private static final double pp3 = -5.77027029648944159157e-03;
    private static final double pp4 = -2.37630166566501626084e-05;

    private static final double qq0 = 3.97917223959155352819e-01;
    private static final double qq1 = 6.50222499887672944485e-02;
    private static final double qq2 = 5.08130628187576562776e-03;
    private static final double qq3 = 1.32494738004321644526e-04;
    private static final double qq4 = -3.96022827877536812320e-06;

    // Coefficients for approximation to erf in [0.84375,1.25]

    private static final double pa0 = -2.36211856075265944077e-03;
    private static final double pa1 = 4.14856118683748331666e-01;
    private static final double pa2 = -3.72207876035701323847e-01;
    private static final double pa3 = 3.18346619901161753674e-01;
    private static final double pa4 = -1.10894694282396677476e-01;
    private static final double pa5 = 3.54783043256182359371e-02;
    private static final double pa6 = -2.16637559486879084300e-03;

    private static final double qa0 = 1.06420880400844228286e-01;
    private static final double qa1 = 5.40397917702171048937e-01;
    private static final double qa2 = 7.18286544141962662868e-02;
    private static final double qa3 = 1.26171219808761642112e-01;
    private static final double qa4 = 1.36370839120290507362e-02;
    private static final double qa5 = 1.19844998467991074170e-02;

    // Coefficients for approximation to erfc in [1.25,1/.35]

    private static final double ra0 = -9.86494403484714822705e-03;
    private static final double ra1 = -6.93858572707181764372e-01;
    private static final double ra2 = -1.05586262253232909814e01;
    private static final double ra3 = -6.23753324503260060396e01;
    private static final double ra4 = -1.62396669462573470355e02;
    private static final double ra5 = -1.84605092906711035994e02;
    private static final double ra6 = -8.12874355063065934246e01;
    private static final double ra7 = -9.81432934416914548592e00;

    private static final double sa0 = 1.96512716674392571292e01;
    private static final double sa1 = 1.37657754143519042600e02;
    private static final double sa2 = 4.34565877475229228821e02;
    private static final double sa3 = 6.45387271733267880336e02;
    private static final double sa4 = 4.29008140027567833386e02;
    private static final double sa5 = 1.08635005541779435134e02;
    private static final double sa6 = 6.57024977031928170135e00;
    private static final double sa7 = -6.04244152148580987438e-02;

    // Coefficients for approximation to erfc in [1/.35,28]

    private static final double rb0 = -9.86494292470009928597e-03;
    private static final double rb1 = -7.99283237680523006574e-01;
    private static final double rb2 = -1.77579549177547519889e01;
    private static final double rb3 = -1.60636384855821916062e02;
    private static final double rb4 = -6.37566443368389627722e02;
    private static final double rb5 = -1.02509513161107724954e03;
    private static final double rb6 = -4.83519191608651397019e02;

    private static final double sb0 = 3.03380607434824582924e01;
    private static final double sb1 = 3.25792512996573918826e02;
    private static final double sb2 = 1.53672958608443695994e03;
    private static final double sb3 = 3.19985821950859553908e03;
    private static final double sb4 = 2.55305040643316442583e03;
    private static final double sb5 = 4.74528541206955367215e02;
    private static final double sb6 = -2.24409524465858183362e01;

    /**
     * erf(). Derived from C code for the error function by SunSoft (1993).
     */
    public static double erf(final double x) {
        final double absX = (x >= 0.0 ? x : -x);
        double ret;
        if (absX < 0.84375) { // 0 < |x| < 0.84375
            if (absX < 3.7252902984619141e-9) { // |x| < 2**-28
                ret = (absX + absX) * efx;
            } else {
                final double s = x * x;
                final double P = pp0 + s
                        * (pp1 + s * (pp2 + s * (pp3 + s * pp4)));
                final double Q = 1.0 + s
                        * (qq0 + s * (qq1 + s * (qq2 + s * (qq3 + s * qq4))));
                ret = (absX + absX) * (P / Q);
            }
        } else if (absX < 1.25) { // 0.84375 < |x| < 1.25
            final double s = absX - 1.0;
            final double P = pa0
                    + s
                    * (pa1 + s
                            * (pa2 + s
                                    * (pa3 + s * (pa4 + s * (pa5 + s * pa6)))));
            final double Q = 1.0
                    + s
                    * (qa0 + s
                            * (qa1 + s
                                    * (qa2 + s * (qa3 + s * (qa4 + s * qa5)))));
            ret = erx + (P / Q);
        } else if (absX >= 6.0) {
            ret = 1.0;
        } else {
            // 1.25 < |x| < 6.0
            ret = 1.0 - erfc(absX);
        }
        return (x >= 0.0) ? ret : -ret;
    }

    /**
     * erfc(). Derived from C code for the complementary error function by
     * SunSoft (1993).
     */
    public static double erfc(final double x) {
        final double absX = (x >= 0.0 ? x : -x);
        double ret;
        if (absX < 1.25) {
            ret = 1.0 - erf(absX);
        } else if (absX > 28.0) {
            ret = 0.0;
        } else { // 1.25 < |x| < 28
            final double s = 1.0 / (absX * absX);
            double R;
            double S;
            if (absX < 2.8571428) { // ( |x| < 1/0.35 )
                R = ra0
                        + s
                        * (ra1 + s
                                * (ra2 + s
                                        * (ra3 + s
                                                * (ra4 + s
                                                        * (ra5 + s
                                                                * (ra6 + s
                                                                        * ra7))))));
                S = 1.0
                        + s
                        * (sa0 + s
                                * (sa1 + s
                                        * (sa2 + s
                                                * (sa3 + s
                                                        * (sa4 + s
                                                                * (sa5 + s
                                                                        * (sa6 + s
                                                                                * sa7)))))));
            } else { // ( |x| > 1/0.35 )
                R = rb0
                        + s
                        * (rb1 + s
                                * (rb2 + s
                                        * (rb3 + s
                                                * (rb4 + s * (rb5 + s * rb6)))));
                S = 1.0
                        + s
                        * (sb0 + s
                                * (sb1 + s
                                        * (sb2 + s
                                                * (sb3 + s
                                                        * (sb4 + s
                                                                * (sb5 + s
                                                                        * sb6))))));
            }
            ret = FastMath.exp((-x * x) - 0.5625 + (R / S)) / absX;
        }
        return (x >= 0.0) ? ret : 2.0 - ret;
    }

    private Erf() {
    }
}
