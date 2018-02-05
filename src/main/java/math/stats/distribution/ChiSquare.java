package math.stats.distribution;

import math.rng.DefaultRng;
import math.rng.PseudoRandom;

/**
 * TODO
 */
public class ChiSquare extends AbstractContinuousDistribution {

    private static final double BIG = 100.0;

    private final int degreesOfFreedom;
    private final Gamma gamma;

    public ChiSquare(final int degreesOfFreedom) {
        this(DefaultRng.newPseudoRandom(), degreesOfFreedom);
    }

    public ChiSquare(final PseudoRandom prng, final int degreesOfFreedom) {
        super(prng);
        if (degreesOfFreedom < 1) {
            throw new IllegalArgumentException("degreesOfFreedom < 1 : "
                    + degreesOfFreedom);
        }
        this.degreesOfFreedom = degreesOfFreedom;
        this.gamma = new Gamma(this.prng, (this.degreesOfFreedom / 2.0), 2.0);
    }

    @Override
    public double pdf(final double x) {
        if (x <= 0.0) {
            return 0.0;
        }
        return gamma.pdf(x);
    }

    @Override
    public double cdf(final double x) {
        if (x <= 0.0) {
            return 0.0;
        }
        if (x >= BIG * degreesOfFreedom) {
            return 1.0;
        }
        return gamma.cdf(x);
    }

    @Override
    public double sample() {
        return gamma.sample();
    }

    /**
     * Inverse of the cumulative chi-squared distribution function.
     * 
     * @return the value X for which P(x&lt;=X).
     */
    public double inverse(double probability) {
        if (probability <= 0.0) {
            return 0.0; // TODO is this correct?
        }
        if (probability >= 1.0) {
            return Double.MAX_VALUE;
        }
        return gamma.inverse(probability);
    }

    @Override
    public double mean() {
        return degreesOfFreedom;
    }

    @Override
    public double variance() {
        return 2.0 * degreesOfFreedom;
    }

    /**
     * @return the degreesOfFreedom
     */
    public int getDegreesOfFreedom() {
        return degreesOfFreedom;
    }
}
