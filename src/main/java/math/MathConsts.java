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

    public static final double BIG = 4.503599627370496e15;

    public static final double BIG_INV = 2.22044604925031308085e-16;

    private MathConsts() {
    }
}
