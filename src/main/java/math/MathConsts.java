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
package math;

/**
 * Some numerical constants from the Cephes library.
 */
public final class MathConsts {

    /** The IEEE 754 machine epsilon from Cephes (2^-53). */
    public static final double MACH_EPS = 1.11022302462515654042E-16; /* 2^-53 */

    public static final double MAX_LOG = 7.09782712893383996732E2;

    public static final double MIN_LOG = -7.451332191019412076235E2;

    public static final double MAX_GAMMA = 171.624376956302725;

    /** sqrt(2) */
    public static final double SQRT_TWO = 1.41421356237309504880E0;

    /** sqrt(2*PI) */
    public static final double SQRT_TWO_PI = 2.50662827463100050242E0;

    /** sqrt(2)/2 */
    public static final double SQRT_TWO_HALF = 7.07106781186547524401E-1;

    /** ln(PI) */
    public static final double LN_PI = 1.14472988584940017414; /* ln(PI) */

    /** ln(10) */
    public static final double LN_10 = 2.302585092994046; /* ln(10) */

    public static final double BIG = 4.503599627370496e15;

    public static final double BIG_INV = 2.22044604925031308085e-16;

    /** 4.450147717014403e-308 */
    public static final double MIN_VAL = 2.0 * Double.MIN_NORMAL;

    /** 5.218048215738236e-15 */
    public static final double MIN_TOL = (45.0 * MACH_EPS) + BIG_INV;

    /** Largest int x such that 10^x is representable (approximately) as double */
    public static final int MAX_X_FOR_10_EXP_X_AS_DOUBLE = 308;

    private MathConsts() {
    }
}
