package math.stats.distribution;

import math.rng.PseudoRandom;

public abstract class AbstractContinuousDistribution implements
        ContinuousDistribution {

    protected final PseudoRandom prng;

    protected AbstractContinuousDistribution(final PseudoRandom prng) {
        this.prng = prng;
    }

    @Override
    public double probability(final double x0, final double x1) {
        if (x0 > x1) {
            throw new IllegalArgumentException("Lower endpoint (" + x0
                    + ") must be less than or equal to upper endpoint (" + x1
                    + ")");
        }
        return cdf(x1) - cdf(x0);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] sample(final int sampleSize) {
        final double[] samples = new double[sampleSize];
        for (int i = 0; i < sampleSize; ++i) {
            samples[i] = sample();
        }
        return samples;
    }

    private static final double FINDROOT_ACCURACY = 1.0e-15;
    private static final int FINDROOT_MAX_ITERATIONS = 150;

    /**
     * This method approximates the value of {@code x} for which
     * {@code P(X <= x) = p} where {@code p} is a given probability.
     * <p>
     * It applies a combination of the Newton-Raphson algorithm and the
     * bisection method to the value {@code start} as a starting point.
     * <p>
     * Furthermore, to ensure convergence and stability, the caller must supply
     * an interval {@code [xMin, xMax]} in which the probability distribution
     * reaches the value {@code p}.
     * <p>
     * Caution: this method does not check its arguments! It will produce wrong
     * results if bad values for the parameters are supplied. To be used with
     * care!
     * 
     * @param p
     *            the given probability for which we want to find the
     *            corresponding value of {@code x} such that
     *            {@code P(X <= x) = p}
     * @param start
     *            an initial guess that must lie in the interval
     *            {@code [xMin, xMax]} as a starting point for the search for
     *            {@code x}
     * @param xMin
     *            lower bound for an interval that must contain the searched
     *            {@code x}
     * @param xMax
     *            upper bound for an interval that must contain the searched
     *            {@code x}
     * @return an approximation for the value of {@code x} for which
     *         {@code P(X <= x) = p}
     */
    protected final double findRoot(double p, double start, double xMin,
            double xMax) {
        double x = start;
        double xNew = start;
        double dx = 1.0;
        int i = 0;
        while (Math.abs(dx) > FINDROOT_ACCURACY
                && i++ < FINDROOT_MAX_ITERATIONS) {
            // apply Newton-Raphson step
            double error = cdf(x) - p;
            if (error < 0.0) {
                xMin = x;
            } else {
                xMax = x;
            }
            double density = pdf(x);
            if (density != 0.0) { // avoid division by zero
                dx = error / density;
                xNew = x - dx;
            }
            // If Newton-Raphson fails to converge (which, for example, may be
            // the case if the initial guess is too rough) we apply a bisection
            // step to determine a more narrow interval around the root
            if (xNew < xMin || xNew > xMax || density == 0.0) {
                xNew = (xMin + xMax) / 2.0;
                dx = xNew - x;
            }
            x = xNew;
        }
        return x;
    }
}
