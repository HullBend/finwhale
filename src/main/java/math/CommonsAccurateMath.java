/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package math;

/**
 * Faster, more accurate, portable alternative to {@link Math} and
 * {@link StrictMath} for large scale computation.
 * <p>
 * CommonsAccurateMath is a drop-in replacement for some methods from Math and
 * StrictMath. This means that for any method in Math (say {@code Math.exp(x)}
 * or {@code Math.cbrt(y)}) that is defined in CommonsAccurateMath, user can
 * directly change the class and use the methods as is (using
 * {@code CommonsAccurateMath.exp(x)} or {@code CommonsAccurateMath.cbrt(y)} in
 * the previous example).
 * </p>
 * <p>
 * CommonsAccurateMath speed is achieved by relying heavily on optimizing
 * compilers to native code present in many JVMs today and use of large tables.
 * The larger tables are lazily initialised on first use, so that the setup time
 * does not penalise methods that don't need them.
 * </p>
 * <p>
 * CommonsAccurateMath accuracy should be mostly independent of the JVM as it
 * relies only on IEEE-754 basic operations and on embedded tables. Almost all
 * operations are accurate to about 0.5 ulp throughout the domain range. This
 * statement, of course is only a rough global observed behavior, it is
 * <em>not</em> a guarantee for <em>every</em> double numbers input (see William
 * Kahan's <a
 * href="http://en.wikipedia.org/wiki/Rounding#The_table-maker.27s_dilemma"
 * >Table Maker's Dilemma</a>).
 * </p>
 * 
 * @version $Id: FastMath.java 1462503 2013-03-29 15:48:27Z luc $
 * @since 2.2
 */
final class CommonsAccurateMath {
    /*
     * There are 52 bits in the mantissa of a double. For additional precision,
     * the code splits double numbers into two parts, by clearing the low order
     * 30 bits if possible, and then performs the arithmetic on each half
     * separately.
     */

    /**
     * 0x40000000 - used to split a double into two parts, both with the low
     * order bits cleared. Equivalent to 2^30.
     */
    private static final long HEX_40000000 = 0x40000000L; // 1073741824L

    /**
     * 2^52 - double numbers this large must be integral (no fraction) or NaN or
     * Infinite
     */
    private static final double TWO_POWER_52 = 4503599627370496.0;
    /** 2^53 - double numbers this large must be even. */
    private static final double TWO_POWER_53 = 2 * TWO_POWER_52;

    private static final boolean RECOMPUTE_TABLES_AT_RUNTIME = false;

    /** Index of exp(0) in the array of integer exponentials. */
    static final int EXP_INT_TABLE_MAX_INDEX = 750;

    /** Length of the array of integer exponentials. */
    static final int EXP_INT_TABLE_LEN = EXP_INT_TABLE_MAX_INDEX * 2;
    /** Logarithm table length. */
    static final int LN_MANT_LEN = 1024;

    /** log(2) (high bits). */
    private static final double LN_2_A = 0.693147063255310059;

    /** log(2) (low bits). */
    private static final double LN_2_B = 1.17304635250823482e-7;

    /** Coefficients for log, when input 0.99 < x < 1.01. */
    private static final double LN_QUICK_COEF[][] = {
            { 1.0, 5.669184079525E-24 }, { -0.25, -0.25 },
            { 0.3333333134651184, 1.986821492305628E-8 },
            { -0.25, -6.663542893624021E-14 },
            { 0.19999998807907104, 1.1921056801463227E-8 },
            { -0.1666666567325592, -7.800414592973399E-9 },
            { 0.1428571343421936, 5.650007086920087E-9 },
            { -0.12502530217170715, -7.44321345601866E-11 },
            { 0.11113807559013367, 9.219544613762692E-9 }, };

    /** Coefficients for log in the range of 1.0 < x < 1.0 + 2^-10. */
    private static final double LN_HI_PREC_COEF[][] = {
            { 1.0, -6.032174644509064E-23 }, { -0.25, -0.25 },
            { 0.3333333134651184, 1.9868161777724352E-8 },
            { -0.2499999701976776, -2.957007209750105E-8 },
            { 0.19999954104423523, 1.5830993332061267E-10 },
            { -0.16624879837036133, -2.6033824355191673E-8 } };

    /** Exponential fractions table length. */
    static final int EXP_FRAC_TABLE_LEN = 1025; // 0, 1/1024, ... 1024/1024

    /** StrictMath.log(Double.MAX_VALUE): {@value} */
    private static final double LOG_MAX_VALUE = StrictMath
            .log(Double.MAX_VALUE);

    /** Table of 2^((n+2)/3) */
    private static final double CBRTTWO[] = { 0.6299605249474366,
            0.7937005259840998, 1.0, 1.2599210498948732, 1.5874010519681994 };

    /**
     * Enclose large data table in nested static class so it's only loaded on
     * first access.
     */
    private static final class ExpIntTable {
        /**
         * Exponential evaluated at integer values, exp(x) = expIntTableA[x +
         * EXP_INT_TABLE_MAX_INDEX] + expIntTableB[x+EXP_INT_TABLE_MAX_INDEX].
         */
        private static final double[] EXP_INT_TABLE_A;
        /**
         * Exponential evaluated at integer values, exp(x) = expIntTableA[x +
         * EXP_INT_TABLE_MAX_INDEX] + expIntTableB[x+EXP_INT_TABLE_MAX_INDEX]
         */
        private static final double[] EXP_INT_TABLE_B;

        static {
            if (RECOMPUTE_TABLES_AT_RUNTIME) {
                EXP_INT_TABLE_A = new double[CommonsAccurateMath.EXP_INT_TABLE_LEN];
                EXP_INT_TABLE_B = new double[CommonsAccurateMath.EXP_INT_TABLE_LEN];

                final double tmp[] = new double[2];
                final double recip[] = new double[2];

                // Populate expIntTable
                for (int i = 0; i < CommonsAccurateMath.EXP_INT_TABLE_MAX_INDEX; i++) {
                    CommonsCalc.expint(i, tmp);
                    EXP_INT_TABLE_A[i
                            + CommonsAccurateMath.EXP_INT_TABLE_MAX_INDEX] = tmp[0];
                    EXP_INT_TABLE_B[i
                            + CommonsAccurateMath.EXP_INT_TABLE_MAX_INDEX] = tmp[1];

                    if (i != 0) {
                        // Negative integer powers
                        CommonsCalc.splitReciprocal(tmp, recip);
                        EXP_INT_TABLE_A[CommonsAccurateMath.EXP_INT_TABLE_MAX_INDEX
                                - i] = recip[0];
                        EXP_INT_TABLE_B[CommonsAccurateMath.EXP_INT_TABLE_MAX_INDEX
                                - i] = recip[1];
                    }
                }
            } else {
                EXP_INT_TABLE_A = CommonsMathLiterals.loadExpIntA();
                EXP_INT_TABLE_B = CommonsMathLiterals.loadExpIntB();
            }
        }
    }

    /**
     * Enclose large data table in nested static class so it's only loaded on
     * first access.
     */
    private static final class ExpFracTable {
        /**
         * Exponential over the range of 0 - 1 in increments of 2^-10
         * exp(x/1024) = expFracTableA[x] + expFracTableB[x]. 1024 = 2^10
         */
        private static final double[] EXP_FRAC_TABLE_A;
        /**
         * Exponential over the range of 0 - 1 in increments of 2^-10
         * exp(x/1024) = expFracTableA[x] + expFracTableB[x].
         */
        private static final double[] EXP_FRAC_TABLE_B;

        static {
            if (RECOMPUTE_TABLES_AT_RUNTIME) {
                EXP_FRAC_TABLE_A = new double[CommonsAccurateMath.EXP_FRAC_TABLE_LEN];
                EXP_FRAC_TABLE_B = new double[CommonsAccurateMath.EXP_FRAC_TABLE_LEN];

                final double tmp[] = new double[2];

                // Populate expFracTable
                final double factor = 1d / (EXP_FRAC_TABLE_LEN - 1);
                for (int i = 0; i < EXP_FRAC_TABLE_A.length; i++) {
                    CommonsCalc.slowexp(i * factor, tmp);
                    EXP_FRAC_TABLE_A[i] = tmp[0];
                    EXP_FRAC_TABLE_B[i] = tmp[1];
                }
            } else {
                EXP_FRAC_TABLE_A = CommonsMathLiterals.loadExpFracA();
                EXP_FRAC_TABLE_B = CommonsMathLiterals.loadExpFracB();
            }
        }
    }

    /**
     * Enclose large data table in nested static class so it's only loaded on
     * first access.
     */
    private static final class lnMant {
        /**
         * Extended precision logarithm table over the range 1 - 2 in increments
         * of 2^-10.
         */
        private static final double[][] LN_MANT;

        static {
            if (RECOMPUTE_TABLES_AT_RUNTIME) {
                LN_MANT = new double[CommonsAccurateMath.LN_MANT_LEN][];

                // Populate lnMant table
                for (int i = 0; i < LN_MANT.length; i++) {
                    final double d = Double
                            .longBitsToDouble((((long) i) << 42) | 0x3ff0000000000000L);
                    LN_MANT[i] = CommonsCalc.slowLog(d);
                }
            } else {
                LN_MANT = CommonsMathLiterals.loadLnMant();
            }
        }
    }

    /**
     * Power function. Compute x^y.
     * 
     * @param x
     *            a double
     * @param y
     *            a double
     * @return double
     */
    static double pow(double x, double y) {
        final double lns[] = new double[2];

        if (y == 0.0) {
            return 1.0;
        }

        if (x != x) { // X is NaN
            return x;
        }

        if (x == 0) {
            long bits = Double.doubleToRawLongBits(x);
            if ((bits & 0x8000000000000000L) != 0) {
                // -zero
                long yi = (long) y;

                if (y < 0 && y == yi && (yi & 1) == 1) {
                    return Double.NEGATIVE_INFINITY;
                }

                if (y > 0 && y == yi && (yi & 1) == 1) {
                    return -0.0;
                }
            }

            if (y < 0) {
                return Double.POSITIVE_INFINITY;
            }
            if (y > 0) {
                return 0.0;
            }

            return Double.NaN;
        }

        if (x == Double.POSITIVE_INFINITY) {
            if (y != y) { // y is NaN
                return y;
            }
            if (y < 0.0) {
                return 0.0;
            } else {
                return Double.POSITIVE_INFINITY;
            }
        }

        if (y == Double.POSITIVE_INFINITY) {
            if (x * x == 1.0) {
                return Double.NaN;
            }

            if (x * x > 1.0) {
                return Double.POSITIVE_INFINITY;
            } else {
                return 0.0;
            }
        }

        if (x == Double.NEGATIVE_INFINITY) {
            if (y != y) { // y is NaN
                return y;
            }

            if (y < 0) {
                long yi = (long) y;
                if (y == yi && (yi & 1) == 1) {
                    return -0.0;
                }

                return 0.0;
            }

            if (y > 0) {
                long yi = (long) y;
                if (y == yi && (yi & 1) == 1) {
                    return Double.NEGATIVE_INFINITY;
                }

                return Double.POSITIVE_INFINITY;
            }
        }

        if (y == Double.NEGATIVE_INFINITY) {

            if (x * x == 1.0) {
                return Double.NaN;
            }

            if (x * x < 1.0) {
                return Double.POSITIVE_INFINITY;
            } else {
                return 0.0;
            }
        }

        /* Handle special case x<0 */
        if (x < 0) {
            // y is an even integer in this case
            if (y >= TWO_POWER_53 || y <= -TWO_POWER_53) {
                return pow(-x, y);
            }

            if (y == (long) y) {
                // If y is an integer
                return ((long) y & 1) == 0 ? pow(-x, y) : -pow(-x, y);
            } else {
                return Double.NaN;
            }
        }

        /* Split y into ya and yb such that y = ya+yb */
        double ya;
        double yb;
        if (y < 8e298 && y > -8e298) {
            double tmp1 = y * HEX_40000000;
            ya = y + tmp1 - tmp1;
            yb = y - ya;
        } else {
            double tmp1 = y * 9.31322574615478515625E-10;
            double tmp2 = tmp1 * 9.31322574615478515625E-10;
            ya = (tmp1 + tmp2 - tmp1) * HEX_40000000 * HEX_40000000;
            yb = y - ya;
        }

        /* Compute ln(x) */
        final double lores = log_(x, lns);
        if (Double.isInfinite(lores)) { // don't allow this to be converted to
                                        // NaN
            return lores;
        }

        double lna = lns[0];
        double lnb = lns[1];

        /* resplit lns */
        double tmp1 = lna * HEX_40000000;
        double tmp2 = lna + tmp1 - tmp1;
        lnb += lna - tmp2;
        lna = tmp2;

        // y*ln(x) = (aa+ab)
        final double aa = lna * ya;
        final double ab = lna * yb + lnb * ya + lnb * yb;

        lna = aa + ab;
        lnb = -(lna - aa - ab);

        double z = 1.0 / 120.0;
        z = z * lnb + (1.0 / 24.0);
        z = z * lnb + (1.0 / 6.0);
        z = z * lnb + 0.5;
        z = z * lnb + 1.0;
        z = z * lnb;

        final double result = exp_(lna, z, null);
        // result = result + result * z;
        return result;
    }

    /**
     * Raise a double to an int power.
     * 
     * @param d
     *            Number to raise.
     * @param e
     *            Exponent.
     * @return d<sup>e</sup>
     * @since 3.1
     */
    @SuppressWarnings("unused")
    private static double pow_(double d, int e) {

        if (e == 0) {
            return 1.0;
        } else if (e < 0) {
            e = -e;
            d = 1.0 / d;
        }

        // split d as two 26 bits numbers
        // beware the following expressions must NOT be simplified, they rely on
        // floating point arithmetic properties
        final int splitFactor = 0x8000001;
        final double cd = splitFactor * d;
        final double d1High = cd - (cd - d);
        final double d1Low = d - d1High;

        // prepare result
        double resultHigh = 1;
        double resultLow = 0;

        // d^(2p)
        double d2p = d;
        double d2pHigh = d1High;
        double d2pLow = d1Low;

        while (e != 0) {

            if ((e & 0x1) != 0) {
                // accurate multiplication result = result * d^(2p) using
                // Veltkamp TwoProduct algorithm
                // beware the following expressions must NOT be simplified, they
                // rely on floating point arithmetic properties
                final double tmpHigh = resultHigh * d2p;
                final double cRH = splitFactor * resultHigh;
                final double rHH = cRH - (cRH - resultHigh);
                final double rHL = resultHigh - rHH;
                final double tmpLow = rHL
                        * d2pLow
                        - (((tmpHigh - rHH * d2pHigh) - rHL * d2pHigh) - rHH
                                * d2pLow);
                resultHigh = tmpHigh;
                resultLow = resultLow * d2p + tmpLow;
            }

            // accurate squaring d^(2(p+1)) = d^(2p) * d^(2p) using Veltkamp
            // TwoProduct algorithm
            // beware the following expressions must NOT be simplified, they
            // rely on floating point arithmetic properties
            final double tmpHigh = d2pHigh * d2p;
            final double cD2pH = splitFactor * d2pHigh;
            final double d2pHH = cD2pH - (cD2pH - d2pHigh);
            final double d2pHL = d2pHigh - d2pHH;
            final double tmpLow = d2pHL
                    * d2pLow
                    - (((tmpHigh - d2pHH * d2pHigh) - d2pHL * d2pHigh) - d2pHH
                            * d2pLow);
            final double cTmpH = splitFactor * tmpHigh;
            d2pHigh = cTmpH - (cTmpH - tmpHigh);
            d2pLow = d2pLow * d2p + tmpLow + (tmpHigh - d2pHigh);
            d2p = d2pHigh + d2pLow;

            e = e >> 1;

        }

        return resultHigh + resultLow;

    }

    /**
     * Exponential function.
     * 
     * Computes exp(x), function result is nearly rounded. It will be correctly
     * rounded to the theoretical value for 99.9% of input values, otherwise it
     * will have a 1 UPL error.
     * 
     * Method: Lookup intVal = exp(int(x)) Lookup fracVal = exp(int(x-int(x) /
     * 1024.0) * 1024.0 ); Compute z as the exponential of the remaining bits by
     * a polynomial minus one exp(x) = intVal * fracVal * (1 + z)
     * 
     * Accuracy: Calculation is done with 63 bits of precision, so result should
     * be correctly rounded for 99.9% of input values, with less than 1 ULP
     * error otherwise.
     * 
     * @param x
     *            a double
     * @return double e<sup>x</sup>
     */
    static double exp(double x) {
        return exp_(x, 0.0, null);
    }

    /**
     * Internal helper method for exponential function.
     * 
     * @param x
     *            original argument of the exponential function
     * @param extra
     *            extra bits of precision on input (To Be Confirmed)
     * @param hiPrec
     *            extra bits of precision on output (To Be Confirmed)
     * @return exp(x)
     */
    private static double exp_(double x, double extra, double[] hiPrec) {
        double intPartA;
        double intPartB;
        int intVal;

        /*
         * Lookup exp(floor(x)). intPartA will have the upper 22 bits, intPartB
         * will have the lower 52 bits.
         */
        if (x < 0.0) {
            intVal = (int) -x;

            if (intVal > 746) {
                if (hiPrec != null) {
                    hiPrec[0] = 0.0;
                    hiPrec[1] = 0.0;
                }
                return 0.0;
            }

            if (intVal > 709) {
                /* This will produce a subnormal output */
                final double result = exp_(x + 40.19140625, extra, hiPrec) / 285040095144011776.0;
                if (hiPrec != null) {
                    hiPrec[0] /= 285040095144011776.0;
                    hiPrec[1] /= 285040095144011776.0;
                }
                return result;
            }

            if (intVal == 709) {
                /* exp(1.494140625) is nearly a machine number... */
                final double result = exp_(x + 1.494140625, extra, hiPrec) / 4.455505956692756620;
                if (hiPrec != null) {
                    hiPrec[0] /= 4.455505956692756620;
                    hiPrec[1] /= 4.455505956692756620;
                }
                return result;
            }

            intVal++;

            intPartA = ExpIntTable.EXP_INT_TABLE_A[EXP_INT_TABLE_MAX_INDEX
                    - intVal];
            intPartB = ExpIntTable.EXP_INT_TABLE_B[EXP_INT_TABLE_MAX_INDEX
                    - intVal];

            intVal = -intVal;
        } else {
            intVal = (int) x;

            if (intVal > 709) {
                if (hiPrec != null) {
                    hiPrec[0] = Double.POSITIVE_INFINITY;
                    hiPrec[1] = 0.0;
                }
                return Double.POSITIVE_INFINITY;
            }

            intPartA = ExpIntTable.EXP_INT_TABLE_A[EXP_INT_TABLE_MAX_INDEX
                    + intVal];
            intPartB = ExpIntTable.EXP_INT_TABLE_B[EXP_INT_TABLE_MAX_INDEX
                    + intVal];
        }

        /*
         * Get the fractional part of x, find the greatest multiple of 2^-10
         * less than x and look up the exp function of it. fracPartA will have
         * the upper 22 bits, fracPartB the lower 52 bits.
         */
        final int intFrac = (int) ((x - intVal) * 1024.0);
        final double fracPartA = ExpFracTable.EXP_FRAC_TABLE_A[intFrac];
        final double fracPartB = ExpFracTable.EXP_FRAC_TABLE_B[intFrac];

        /*
         * epsilon is the difference in x from the nearest multiple of 2^-10. It
         * has a value in the range 0 <= epsilon < 2^-10. Do the subtraction
         * from x as the last step to avoid possible loss of precision.
         */
        final double epsilon = x - (intVal + intFrac / 1024.0);

        /*
         * Compute z = exp(epsilon) - 1.0 via a minimax polynomial. z has full
         * double precision (52 bits). Since z < 2^-10, we will have 62 bits of
         * precision when combined with the constant 1. This will be used in the
         * last addition below to get proper rounding.
         */

        /*
         * Remez generated polynomial. Converges on the interval [0, 2^-10],
         * error is less than 0.5 ULP
         */
        double z = 0.04168701738764507;
        z = z * epsilon + 0.1666666505023083;
        z = z * epsilon + 0.5000000000042687;
        z = z * epsilon + 1.0;
        z = z * epsilon + -3.940510424527919E-20;

        /*
         * Compute (intPartA+intPartB) * (fracPartA+fracPartB) by binomial
         * expansion. tempA is exact since intPartA and intPartB only have 22
         * bits each. tempB will have 52 bits of precision.
         */
        double tempA = intPartA * fracPartA;
        double tempB = intPartA * fracPartB + intPartB * fracPartA + intPartB
                * fracPartB;

        /*
         * Compute the result. (1+z)(tempA+tempB). Order of operations is
         * important. For accuracy add by increasing size. tempA is exact and
         * much larger than the others. If there are extra bits specified from
         * the pow() function, use them.
         */
        final double tempC = tempB + tempA;
        final double result;
        if (extra != 0.0) {
            result = tempC * extra * z + tempC * extra + tempC * z + tempB
                    + tempA;
        } else {
            result = tempC * z + tempB + tempA;
        }

        if (hiPrec != null) {
            // If requesting high precision
            hiPrec[0] = tempA;
            hiPrec[1] = tempC * extra * z + tempC * extra + tempC * z + tempB;
        }

        return result;
    }

    /**
     * Internal helper method for expm1
     * 
     * @param x
     *            number to compute shifted exponential
     * @param hiPrecOut
     *            receive high precision result for -1.0 < x < 1.0
     * @return exp(x) - 1
     */
    private static double expm1_(double x, double hiPrecOut[]) {
        if (x != x || x == 0.0) { // NaN or zero
            return x;
        }

        if (x <= -1.0 || x >= 1.0) {
            // If not between +/- 1.0
            // return exp(x) - 1.0;
            double hiPrec[] = new double[2];
            exp_(x, 0.0, hiPrec);
            if (x > 0.0) {
                return -1.0 + hiPrec[0] + hiPrec[1];
            } else {
                final double ra = -1.0 + hiPrec[0];
                double rb = -(ra + 1.0 - hiPrec[0]);
                rb += hiPrec[1];
                return ra + rb;
            }
        }

        double baseA;
        double baseB;
        double epsilon;
        boolean negative = false;

        if (x < 0.0) {
            x = -x;
            negative = true;
        }

        {
            int intFrac = (int) (x * 1024.0);
            double tempA = ExpFracTable.EXP_FRAC_TABLE_A[intFrac] - 1.0;
            double tempB = ExpFracTable.EXP_FRAC_TABLE_B[intFrac];

            double temp = tempA + tempB;
            tempB = -(temp - tempA - tempB);
            tempA = temp;

            temp = tempA * HEX_40000000;
            baseA = tempA + temp - temp;
            baseB = tempB + (tempA - baseA);

            epsilon = x - intFrac / 1024.0;
        }

        /* Compute expm1(epsilon) */
        double zb = 0.008336750013465571;
        zb = zb * epsilon + 0.041666663879186654;
        zb = zb * epsilon + 0.16666666666745392;
        zb = zb * epsilon + 0.49999999999999994;
        zb = zb * epsilon;
        zb = zb * epsilon;

        double za = epsilon;
        double temp = za + zb;
        zb = -(temp - za - zb);
        za = temp;

        temp = za * HEX_40000000;
        temp = za + temp - temp;
        zb += za - temp;
        za = temp;

        /*
         * Combine the parts. expm1(a+b) = expm1(a) + expm1(b) +
         * expm1(a)*expm1(b)
         */
        double ya = za * baseA;
        // double yb = za*baseB + zb*baseA + zb*baseB;
        temp = ya + za * baseB;
        double yb = -(temp - ya - za * baseB);
        ya = temp;

        temp = ya + zb * baseA;
        yb += -(temp - ya - zb * baseA);
        ya = temp;

        temp = ya + zb * baseB;
        yb += -(temp - ya - zb * baseB);
        ya = temp;

        // ya = ya + za + baseA;
        // yb = yb + zb + baseB;
        temp = ya + baseA;
        yb += -(temp - baseA - ya);
        ya = temp;

        temp = ya + za;
        // yb += (ya > za) ? -(temp - ya - za) : -(temp - za - ya);
        yb += -(temp - ya - za);
        ya = temp;

        temp = ya + baseB;
        // yb += (ya > baseB) ? -(temp - ya - baseB) : -(temp - baseB - ya);
        yb += -(temp - ya - baseB);
        ya = temp;

        temp = ya + zb;
        // yb += (ya > zb) ? -(temp - ya - zb) : -(temp - zb - ya);
        yb += -(temp - ya - zb);
        ya = temp;

        if (negative) {
            /* Compute expm1(-x) = -expm1(x) / (expm1(x) + 1) */
            double denom = 1.0 + ya;
            double denomr = 1.0 / denom;
            double denomb = -(denom - 1.0 - ya) + yb;
            double ratio = ya * denomr;
            temp = ratio * HEX_40000000;
            final double ra = ratio + temp - temp;
            double rb = ratio - ra;

            temp = denom * HEX_40000000;
            za = denom + temp - temp;
            zb = denom - za;

            rb += (ya - za * ra - za * rb - zb * ra - zb * rb) * denomr;

            // f(x) = x/1+x
            // Compute f'(x)
            // Product rule: d(uv) = du*v + u*dv
            // Chain rule: d(f(g(x)) = f'(g(x))*f(g'(x))
            // d(1/x) = -1/(x*x)
            // d(1/1+x) = -1/( (1+x)^2) * 1 = -1/((1+x)*(1+x))
            // d(x/1+x) = -x/((1+x)(1+x)) + 1/1+x = 1 / ((1+x)(1+x))

            // Adjust for yb
            rb += yb * denomr; // numerator
            rb += -ya * denomb * denomr * denomr; // denominator

            // negate
            ya = -ra;
            yb = -rb;
        }

        if (hiPrecOut != null) {
            hiPrecOut[0] = ya;
            hiPrecOut[1] = yb;
        }

        return ya + yb;
    }

    /**
     * Compute the hyperbolic cosine of a number.
     * 
     * @param x
     *            number on which evaluation is done
     * @return hyperbolic cosine of x
     */
    static double cosh(double x) {
        if (x != x) {
            return x;
        }

        // cosh[z] = (exp(z) + exp(-z))/2

        // for numbers with magnitude 20 or so,
        // exp(-z) can be ignored in comparison with exp(z)

        if (x > 20) {
            if (x >= LOG_MAX_VALUE) {
                // Avoid overflow (MATH-905).
                final double t = exp(0.5 * x);
                return (0.5 * t) * t;
            } else {
                return 0.5 * exp(x);
            }
        } else if (x < -20) {
            if (x <= -LOG_MAX_VALUE) {
                // Avoid overflow (MATH-905).
                final double t = exp(-0.5 * x);
                return (0.5 * t) * t;
            } else {
                return 0.5 * exp(-x);
            }
        }

        final double hiPrec[] = new double[2];
        if (x < 0.0) {
            x = -x;
        }
        exp_(x, 0.0, hiPrec);

        double ya = hiPrec[0] + hiPrec[1];
        double yb = -(ya - hiPrec[0] - hiPrec[1]);

        double temp = ya * HEX_40000000;
        double yaa = ya + temp - temp;
        double yab = ya - yaa;

        // recip = 1/y
        double recip = 1.0 / ya;
        temp = recip * HEX_40000000;
        double recipa = recip + temp - temp;
        double recipb = recip - recipa;

        // Correct for rounding in division
        recipb += (1.0 - yaa * recipa - yaa * recipb - yab * recipa - yab
                * recipb)
                * recip;
        // Account for yb
        recipb += -yb * recip * recip;

        // y = y + 1/y
        temp = ya + recipa;
        yb += -(temp - ya - recipa);
        ya = temp;
        temp = ya + recipb;
        yb += -(temp - ya - recipb);
        ya = temp;

        double result = ya + yb;
        result *= 0.5;
        return result;
    }

    /**
     * Compute the hyperbolic sine of a number.
     * 
     * @param x
     *            number on which evaluation is done
     * @return hyperbolic sine of x
     */
    static double sinh(double x) {
        boolean negate = false;
        if (x != x) {
            return x;
        }

        // sinh[z] = (exp(z) - exp(-z) / 2

        // for values of z larger than about 20,
        // exp(-z) can be ignored in comparison with exp(z)

        if (x > 20) {
            if (x >= LOG_MAX_VALUE) {
                // Avoid overflow (MATH-905).
                final double t = exp(0.5 * x);
                return (0.5 * t) * t;
            } else {
                return 0.5 * exp(x);
            }
        } else if (x < -20) {
            if (x <= -LOG_MAX_VALUE) {
                // Avoid overflow (MATH-905).
                final double t = exp(-0.5 * x);
                return (-0.5 * t) * t;
            } else {
                return -0.5 * exp(-x);
            }
        }

        if (x == 0) {
            return x;
        }

        if (x < 0.0) {
            x = -x;
            negate = true;
        }

        double result;

        if (x > 0.25) {
            double hiPrec[] = new double[2];
            exp_(x, 0.0, hiPrec);

            double ya = hiPrec[0] + hiPrec[1];
            double yb = -(ya - hiPrec[0] - hiPrec[1]);

            double temp = ya * HEX_40000000;
            double yaa = ya + temp - temp;
            double yab = ya - yaa;

            // recip = 1/y
            double recip = 1.0 / ya;
            temp = recip * HEX_40000000;
            double recipa = recip + temp - temp;
            double recipb = recip - recipa;

            // Correct for rounding in division
            recipb += (1.0 - yaa * recipa - yaa * recipb - yab * recipa - yab
                    * recipb)
                    * recip;
            // Account for yb
            recipb += -yb * recip * recip;

            recipa = -recipa;
            recipb = -recipb;

            // y = y + 1/y
            temp = ya + recipa;
            yb += -(temp - ya - recipa);
            ya = temp;
            temp = ya + recipb;
            yb += -(temp - ya - recipb);
            ya = temp;

            result = ya + yb;
            result *= 0.5;
        } else {
            double hiPrec[] = new double[2];
            expm1_(x, hiPrec);

            double ya = hiPrec[0] + hiPrec[1];
            double yb = -(ya - hiPrec[0] - hiPrec[1]);

            /* Compute expm1(-x) = -expm1(x) / (expm1(x) + 1) */
            double denom = 1.0 + ya;
            double denomr = 1.0 / denom;
            double denomb = -(denom - 1.0 - ya) + yb;
            double ratio = ya * denomr;
            double temp = ratio * HEX_40000000;
            double ra = ratio + temp - temp;
            double rb = ratio - ra;

            temp = denom * HEX_40000000;
            double za = denom + temp - temp;
            double zb = denom - za;

            rb += (ya - za * ra - za * rb - zb * ra - zb * rb) * denomr;

            // Adjust for yb
            rb += yb * denomr; // numerator
            rb += -ya * denomb * denomr * denomr; // denominator

            // y = y - 1/y
            temp = ya + ra;
            yb += -(temp - ya - ra);
            ya = temp;
            temp = ya + rb;
            yb += -(temp - ya - rb);
            ya = temp;

            result = ya + yb;
            result *= 0.5;
        }

        if (negate) {
            result = -result;
        }

        return result;
    }

    /**
     * Compute the hyperbolic tangent of a number.
     * 
     * @param x
     *            number on which evaluation is done
     * @return hyperbolic tangent of x
     */
    static double tanh(double x) {
        boolean negate = false;

        if (x != x) {
            return x;
        }

        // tanh[z] = sinh[z] / cosh[z]
        // = (exp(z) - exp(-z)) / (exp(z) + exp(-z))
        // = (exp(2x) - 1) / (exp(2x) + 1)

        // for magnitude > 20, sinh[z] == cosh[z] in double precision

        if (x > 20.0) {
            return 1.0;
        }

        if (x < -20) {
            return -1.0;
        }

        if (x == 0) {
            return x;
        }

        if (x < 0.0) {
            x = -x;
            negate = true;
        }

        double result;
        if (x >= 0.5) {
            double hiPrec[] = new double[2];
            // tanh(x) = (exp(2x) - 1) / (exp(2x) + 1)
            exp_(x * 2.0, 0.0, hiPrec);

            double ya = hiPrec[0] + hiPrec[1];
            double yb = -(ya - hiPrec[0] - hiPrec[1]);

            /* Numerator */
            double na = -1.0 + ya;
            double nb = -(na + 1.0 - ya);
            double temp = na + yb;
            nb += -(temp - na - yb);
            na = temp;

            /* Denominator */
            double da = 1.0 + ya;
            double db = -(da - 1.0 - ya);
            temp = da + yb;
            db += -(temp - da - yb);
            da = temp;

            temp = da * HEX_40000000;
            double daa = da + temp - temp;
            double dab = da - daa;

            // ratio = na/da
            double ratio = na / da;
            temp = ratio * HEX_40000000;
            double ratioa = ratio + temp - temp;
            double ratiob = ratio - ratioa;

            // Correct for rounding in division
            ratiob += (na - daa * ratioa - daa * ratiob - dab * ratioa - dab
                    * ratiob)
                    / da;

            // Account for nb
            ratiob += nb / da;
            // Account for db
            ratiob += -db * na / da / da;

            result = ratioa + ratiob;
        } else {
            double hiPrec[] = new double[2];
            // tanh(x) = expm1(2x) / (expm1(2x) + 2)
            expm1_(x * 2.0, hiPrec);

            double ya = hiPrec[0] + hiPrec[1];
            double yb = -(ya - hiPrec[0] - hiPrec[1]);

            /* Numerator */
            double na = ya;
            double nb = yb;

            /* Denominator */
            double da = 2.0 + ya;
            double db = -(da - 2.0 - ya);
            double temp = da + yb;
            db += -(temp - da - yb);
            da = temp;

            temp = da * HEX_40000000;
            double daa = da + temp - temp;
            double dab = da - daa;

            // ratio = na/da
            double ratio = na / da;
            temp = ratio * HEX_40000000;
            double ratioa = ratio + temp - temp;
            double ratiob = ratio - ratioa;

            // Correct for rounding in division
            ratiob += (na - daa * ratioa - daa * ratiob - dab * ratioa - dab
                    * ratiob)
                    / da;

            // Account for nb
            ratiob += nb / da;
            // Account for db
            ratiob += -db * na / da / da;

            result = ratioa + ratiob;
        }

        if (negate) {
            result = -result;
        }

        return result;
    }

    /**
     * Compute the cubic root of a number.
     * 
     * @param x
     *            number on which evaluation is done
     * @return cubic root of x
     */
    static double cbrt(double x) {
        /* Convert input double to bits */
        long inbits = Double.doubleToRawLongBits(x);
        int exponent = (int) ((inbits >> 52) & 0x7ff) - 1023;
        boolean subnormal = false;

        if (exponent == -1023) {
            if (x == 0) {
                return x;
            }

            /* Subnormal, so normalize */
            subnormal = true;
            x *= 1.8014398509481984E16; // 2^54
            inbits = Double.doubleToRawLongBits(x);
            exponent = (int) ((inbits >> 52) & 0x7ff) - 1023;
        }

        if (exponent == 1024) {
            // Nan or infinity. Don't care which.
            return x;
        }

        /* Divide the exponent by 3 */
        int exp3 = exponent / 3;

        /* p2 will be the nearest power of 2 to x with its exponent divided by 3 */
        double p2 = Double.longBitsToDouble((inbits & 0x8000000000000000L)
                | (long) (((exp3 + 1023) & 0x7ff)) << 52);

        /* This will be a number between 1 and 2 */
        final double mant = Double
                .longBitsToDouble((inbits & 0x000fffffffffffffL) | 0x3ff0000000000000L);

        /* Estimate the cube root of mant by polynomial */
        double est = -0.010714690733195933;
        est = est * mant + 0.0875862700108075;
        est = est * mant + -0.3058015757857271;
        est = est * mant + 0.7249995199969751;
        est = est * mant + 0.5039018405998233;

        est *= CBRTTWO[exponent % 3 + 2];

        // est should now be good to about 15 bits of precision. Do 2 rounds of
        // Newton's method to get closer, this should get us full double
        // precision
        // Scale down x for the purpose of doing newtons method. This avoids
        // over/under flows.
        final double xs = x / (p2 * p2 * p2);
        est += (xs - est * est * est) / (3 * est * est);
        est += (xs - est * est * est) / (3 * est * est);

        // Do one round of Newton's method in extended precision to get the last
        // bit right.
        double temp = est * HEX_40000000;
        double ya = est + temp - temp;
        double yb = est - ya;

        double za = ya * ya;
        double zb = ya * yb * 2.0 + yb * yb;
        temp = za * HEX_40000000;
        double temp2 = za + temp - temp;
        zb += za - temp2;
        za = temp2;

        zb = za * yb + ya * zb + zb * yb;
        za = za * ya;

        double na = xs - za;
        double nb = -(na - xs + za);
        nb -= zb;

        est += (na + nb) / (3 * est * est);

        /* Scale by a power of two, so this is exact. */
        est *= p2;

        if (subnormal) {
            est *= 3.814697265625E-6; // 2^-18
        }

        return est;
    }

    /**
     * Internal helper method for natural logarithm function.
     * 
     * @param x
     *            original argument of the natural logarithm function
     * @param hiPrec
     *            extra bits of precision on output (To Be Confirmed)
     * @return log(x)
     */
    private static double log_(final double x, final double[] hiPrec) {
        if (x == 0) { // Handle special case of +0/-0
            return Double.NEGATIVE_INFINITY;
        }
        long bits = Double.doubleToRawLongBits(x);

        /* Handle special cases of negative input, and NaN */
        if (((bits & 0x8000000000000000L) != 0 || x != x) && x != 0.0) {
            if (hiPrec != null) {
                hiPrec[0] = Double.NaN;
            }

            return Double.NaN;
        }

        /* Handle special cases of Positive infinity. */
        if (x == Double.POSITIVE_INFINITY) {
            if (hiPrec != null) {
                hiPrec[0] = Double.POSITIVE_INFINITY;
            }

            return Double.POSITIVE_INFINITY;
        }

        /* Extract the exponent */
        int exp = (int) (bits >> 52) - 1023;

        if ((bits & 0x7ff0000000000000L) == 0) {
            // Subnormal!
            if (x == 0) {
                // Zero
                if (hiPrec != null) {
                    hiPrec[0] = Double.NEGATIVE_INFINITY;
                }

                return Double.NEGATIVE_INFINITY;
            }

            /* Normalize the subnormal number. */
            bits <<= 1;
            while ((bits & 0x0010000000000000L) == 0) {
                --exp;
                bits <<= 1;
            }
        }

        if ((exp == -1 || exp == 0) && x < 1.01 && x > 0.99 && hiPrec == null) {
            /*
             * The normal method doesn't work well in the range [0.99, 1.01], so
             * call do a straight polynomial expansion in higer precision.
             */

            /* Compute x - 1.0 and split it */
            double xa = x - 1.0;
            double xb = xa - x + 1.0;
            double tmp = xa * HEX_40000000;
            double aa = xa + tmp - tmp;
            double ab = xa - aa;
            xa = aa;
            xb = ab;

            final double[] lnCoef_last = LN_QUICK_COEF[LN_QUICK_COEF.length - 1];
            double ya = lnCoef_last[0];
            double yb = lnCoef_last[1];

            for (int i = LN_QUICK_COEF.length - 2; i >= 0; i--) {
                /* Multiply a = y * x */
                aa = ya * xa;
                ab = ya * xb + yb * xa + yb * xb;
                /* split, so now y = a */
                tmp = aa * HEX_40000000;
                ya = aa + tmp - tmp;
                yb = aa - ya + ab;

                /* Add a = y + lnQuickCoef */
                final double[] lnCoef_i = LN_QUICK_COEF[i];
                aa = ya + lnCoef_i[0];
                ab = yb + lnCoef_i[1];
                /* Split y = a */
                tmp = aa * HEX_40000000;
                ya = aa + tmp - tmp;
                yb = aa - ya + ab;
            }

            /* Multiply a = y * x */
            aa = ya * xa;
            ab = ya * xb + yb * xa + yb * xb;
            /* split, so now y = a */
            tmp = aa * HEX_40000000;
            ya = aa + tmp - tmp;
            yb = aa - ya + ab;

            return ya + yb;
        }

        // lnm is a log of a number in the range of 1.0 - 2.0, so 0 <= lnm <
        // ln(2)
        final double[] lnm = lnMant.LN_MANT[(int) ((bits & 0x000ffc0000000000L) >> 42)];

        /*
         * double epsilon = x / Double.longBitsToDouble(bits &
         * 0xfffffc0000000000L);
         * 
         * epsilon -= 1.0;
         */

        // y is the most significant 10 bits of the mantissa
        // double y = Double.longBitsToDouble(bits & 0xfffffc0000000000L);
        // double epsilon = (x - y) / y;
        final double epsilon = (bits & 0x3ffffffffffL)
                / (TWO_POWER_52 + (bits & 0x000ffc0000000000L));

        double lnza = 0.0;
        double lnzb = 0.0;

        if (hiPrec != null) {
            /* split epsilon -> x */
            double tmp = epsilon * HEX_40000000;
            double aa = epsilon + tmp - tmp;
            double ab = epsilon - aa;
            double xa = aa;
            double xb = ab;

            /* Need a more accurate epsilon, so adjust the division. */
            final double numer = bits & 0x3ffffffffffL;
            final double denom = TWO_POWER_52 + (bits & 0x000ffc0000000000L);
            aa = numer - xa * denom - xb * denom;
            xb += aa / denom;

            /* Remez polynomial evaluation */
            final double[] lnCoef_last = LN_HI_PREC_COEF[LN_HI_PREC_COEF.length - 1];
            double ya = lnCoef_last[0];
            double yb = lnCoef_last[1];

            for (int i = LN_HI_PREC_COEF.length - 2; i >= 0; i--) {
                /* Multiply a = y * x */
                aa = ya * xa;
                ab = ya * xb + yb * xa + yb * xb;
                /* split, so now y = a */
                tmp = aa * HEX_40000000;
                ya = aa + tmp - tmp;
                yb = aa - ya + ab;

                /* Add a = y + lnHiPrecCoef */
                final double[] lnCoef_i = LN_HI_PREC_COEF[i];
                aa = ya + lnCoef_i[0];
                ab = yb + lnCoef_i[1];
                /* Split y = a */
                tmp = aa * HEX_40000000;
                ya = aa + tmp - tmp;
                yb = aa - ya + ab;
            }

            /* Multiply a = y * x */
            aa = ya * xa;
            ab = ya * xb + yb * xa + yb * xb;

            /* split, so now lnz = a */
            /*
             * tmp = aa * 1073741824.0; lnza = aa + tmp - tmp; lnzb = aa - lnza
             * + ab;
             */
            lnza = aa + ab;
            lnzb = -(lnza - aa - ab);
        } else {
            /*
             * High precision not required. Eval Remez polynomial using standard
             * double precision
             */
            lnza = -0.16624882440418567;
            lnza = lnza * epsilon + 0.19999954120254515;
            lnza = lnza * epsilon + -0.2499999997677497;
            lnza = lnza * epsilon + 0.3333333333332802;
            lnza = lnza * epsilon + -0.5;
            lnza = lnza * epsilon + 1.0;
            lnza = lnza * epsilon;
        }

        /*
         * Relative sizes: lnzb [0, 2.33E-10] lnm[1] [0, 1.17E-7] ln2B*exp [0,
         * 1.12E-4] lnza [0, 9.7E-4] lnm[0] [0, 0.692] ln2A*exp [0, 709]
         */

        /*
         * Compute the following sum: lnzb + lnm[1] + ln2B*exp + lnza + lnm[0] +
         * ln2A*exp;
         */

        // return lnzb + lnm[1] + ln2B*exp + lnza + lnm[0] + ln2A*exp;
        double a = LN_2_A * exp;
        double b = 0.0;
        double c = a + lnm[0];
        double d = -(c - a - lnm[0]);
        a = c;
        b = b + d;

        c = a + lnza;
        d = -(c - a - lnza);
        a = c;
        b = b + d;

        c = a + LN_2_B * exp;
        d = -(c - a - LN_2_B * exp);
        a = c;
        b = b + d;

        c = a + lnm[1];
        d = -(c - a - lnm[1]);
        a = c;
        b = b + d;

        c = a + lnzb;
        d = -(c - a - lnzb);
        a = c;
        b = b + d;

        if (hiPrec != null) {
            hiPrec[0] = a;
            hiPrec[1] = b;
        }

        return a + b;
    }
}
