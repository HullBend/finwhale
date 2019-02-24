package math.stats.descriptive;

import java.util.Objects;

import math.function.DoubleConsumer;

/**
 * A state object for collecting statistics such as count, min, max, sum,
 * average and variance (or standard deviation).
 * <p>
 * <b>Implementation Note:</b><br>
 * This implementation is <b>not</b> thread-safe.
 */
public class DoubleSummaryStatistics implements DoubleConsumer {
    private long count;
    private double sum;
    private double sumCompensation; // Negative low order bits of sum
    private double min = Double.POSITIVE_INFINITY;
    private double max = Double.NEGATIVE_INFINITY;
    /**
     * the sum of squares of differences from the (current) mean
     * http://en.wikipedia
     * .org/wiki/Algorithms_for_calculating_variance#On-line_algorithm
     * (Welford's algorithm)
     */
    private double sumDiffFromCurrMeanSquared;
    private double sumDiffFromCurrMeanSquaredCompensation;
    /**
     * the variance - recursively calculated via Welford's algorithm
     */
    private double variance;

    /**
     * Construct an empty instance with zero count, zero sum,
     * {@code Double.POSITIVE_INFINITY} min, {@code Double.NEGATIVE_INFINITY}
     * max and zero average.
     */
    public DoubleSummaryStatistics() {
    }

    /**
     * Constructs a non-empty instance with an initial state that corresponds to
     * the current state of the specified {@code other} DoubleSummaryStatistics
     * instance.
     *
     * <p>
     * If {@code other.count} is zero then the remaining arguments are ignored
     * and an empty instance is constructed.
     *
     * <p>
     * If the state of {@code other} is inconsistent then an
     * {@code IllegalArgumentException} is thrown. The necessary conditions for
     * a consistent state are:
     * <ul>
     * <li>{@code other.count >= 0}</li>
     * <li>{@code (other.min <= other.max && !isNaN(other.sum)) || (isNaN(other.min) && isNaN(other.max) && isNaN(other.sum))}</li>
     * </ul>
     * <p>
     * <b>API Note:</b><br>
     * The enforcement of state correctness means that the retrieved set of
     * recorded values obtained from a {@code DoubleSummaryStatistics} source
     * instance may not be a legal state for this constructor due to arithmetic
     * overflow of the source's recorded count of values. The consistency
     * conditions are not sufficient to prevent the creation of an internally
     * inconsistent instance. An example of such a state would be an instance
     * with: {@code other.count} = 2, {@code other.min} = 1, {@code other.max} =
     * 2, and {@code other.sum} = 0.
     *
     * @param other
     *            the DoubleSummaryStatistics instance whose state should be
     *            replicated
     * @throws NullPointerException
     *             if {@code other} is null
     * @throws IllegalArgumentException
     *             if the internal state of the {@code other} object is
     *             inconsistent
     */
    public DoubleSummaryStatistics(DoubleSummaryStatistics other) throws IllegalArgumentException {
        Objects.requireNonNull(other);
        if (other.count < 0L) {
            throw new IllegalArgumentException("Negative count value");
        } else if (other.count > 0L) {
            if (other.min > other.max) {
                throw new IllegalArgumentException("Minimum greater than maximum");
            }
            // All NaN or non NaN
            int ncount = 0;
            if (Double.isNaN(other.min)) {
                ++ncount;
            }
            if (Double.isNaN(other.max)) {
                ++ncount;
            }
            if (Double.isNaN(other.sum)) {
                ++ncount;
            }
            if (ncount > 0 && ncount < 3) {
                throw new IllegalArgumentException("Some, not all, of the minimum, maximum, or sum is NaN");
            }

            this.count = other.count;
            this.sum = other.sum;
            this.sumCompensation = 0.0d;
            this.min = other.min;
            this.max = other.max;
            this.sumDiffFromCurrMeanSquared = other.sumDiffFromCurrMeanSquared;
            this.sumDiffFromCurrMeanSquaredCompensation = other.sumDiffFromCurrMeanSquaredCompensation;
            this.variance = other.variance;
        }
        // Use default field values if count == 0
    }

    /**
     * Records another value into the summary information.
     * 
     * @param value
     *            the input value
     */
    @Override
    public void accept(double value) {
        long countSoFar = count;
        double average = getAverage();
        ++count;
        sumWithCompensation(value);
        min = Math.min(min, value);
        max = Math.max(max, value);
        double delta = value - average;
        average = ((countSoFar * average) / count) + (value / count);
        sumDiffFromCurrMeanSquaredWithCompensation(delta * (value - average));
        if (count > 1L) {
            variance = (sumDiffFromCurrMeanSquared - sumDiffFromCurrMeanSquaredCompensation) / countSoFar;
        }
    }

    /**
     * Incorporate a new double value using Kahan summation / compensated
     * summation.
     */
    private void sumWithCompensation(double value) {
        // https://en.wikipedia.org/wiki/Kahan_summation_algorithm
        double tmp = value - sumCompensation;
        double velvel = sum + tmp; // Little wolf of rounding error
        sumCompensation = (velvel - sum) - tmp;
        sum = velvel;
    }

    /**
     * Incorporate a new double value using Kahan summation / compensated
     * summation.
     */
    private void sumDiffFromCurrMeanSquaredWithCompensation(double value) {
        double tmp = value - sumDiffFromCurrMeanSquaredCompensation;
        double velvel = sumDiffFromCurrMeanSquared + tmp;
        sumDiffFromCurrMeanSquaredCompensation = (velvel - sumDiffFromCurrMeanSquared) - tmp;
        sumDiffFromCurrMeanSquared = velvel;
    }

    /**
     * Return the count of values recorded.
     * 
     * @return the count of values
     */
    public final long getCount() {
        return count;
    }

    /**
     * Returns the sum of values recorded, or zero if no values have been
     * recorded.
     * 
     * <p>
     * The value of a floating-point sum is a function both of the input values
     * as well as the order of addition operations. The order of addition
     * operations of this method is intentionally not defined to allow for
     * implementation flexibility to improve the speed and accuracy of the
     * computed result.
     * 
     * In particular, this method may be implemented using compensated summation
     * or other technique to reduce the error bound in the numerical sum
     * compared to a simple summation of {@code double} values. Because of the
     * unspecified order of operations and the possibility of using differing
     * summation schemes, the output of this method may vary on the same input
     * values.
     * 
     * <p>
     * Various conditions can result in a non-finite sum being computed. This
     * can occur even if the all the recorded values being summed are finite. If
     * any recorded value is non-finite, the sum will be non-finite:
     * 
     * <ul>
     * 
     * <li>If any recorded value is a NaN, then the final sum will be NaN.
     * 
     * <li>If the recorded values contain one or more infinities, the sum will
     * be infinite or NaN.
     * 
     * <ul>
     * 
     * <li>If the recorded values contain infinities of opposite sign, the sum
     * will be NaN.
     * 
     * <li>If the recorded values contain infinities of one sign and an
     * intermediate sum overflows to an infinity of the opposite sign, the sum
     * may be NaN.
     * 
     * </ul>
     * 
     * </ul>
     * 
     * It is possible for intermediate sums of finite values to overflow into
     * opposite-signed infinities; if that occurs, the final sum will be NaN
     * even if the recorded values are all finite.
     * 
     * If all the recorded values are zero, the sign of zero is <em>not</em>
     * guaranteed to be preserved in the final sum.
     * 
     * <p>
     * <b>API Note:</b><br>
     * Values sorted by increasing absolute magnitude tend to yield more
     * accurate results.
     * 
     * @return the sum of values, or zero if none
     */
    public final double getSum() {
        // Better error bounds to add both terms as the final sum
        return sum - sumCompensation;
    }

    /**
     * Returns the minimum recorded value, {@code Double.NaN} if any recorded
     * value was NaN or {@code Double.POSITIVE_INFINITY} if no values were
     * recorded. Unlike the numerical comparison operators, this method
     * considers negative zero to be strictly smaller than positive zero.
     * 
     * @return the minimum recorded value, {@code Double.NaN} if any recorded
     *         value was NaN or {@code Double.POSITIVE_INFINITY} if no values
     *         were recorded
     */
    public final double getMin() {
        return min;
    }

    /**
     * Returns the maximum recorded value, {@code Double.NaN} if any recorded
     * value was NaN or {@code Double.NEGATIVE_INFINITY} if no values were
     * recorded. Unlike the numerical comparison operators, this method
     * considers negative zero to be strictly smaller than positive zero.
     * 
     * @return the maximum recorded value, {@code Double.NaN} if any recorded
     *         value was NaN or {@code Double.NEGATIVE_INFINITY} if no values
     *         were recorded
     */
    public final double getMax() {
        return max;
    }

    /**
     * Returns the arithmetic mean of values recorded, or zero if no values have
     * been recorded.
     * 
     * <p>
     * The computed average can vary numerically and have the special case
     * behavior as computing the sum; see {@link #getSum} for details.
     * 
     * <p>
     * <b>API Note:</b><br>
     * Values sorted by increasing absolute magnitude tend to yield more
     * accurate results.
     * 
     * @return the arithmetic mean of values, or zero if none
     */
    public final double getAverage() {
        return getCount() > 0 ? getSum() / getCount() : 0.0d;
    }

    public final double getVariance() {
        double var = variance;
        if (var < 0.0) {
            var = 0.0;
        }
        return var;
    }

    public final double getStandardDeviation() {
        return Math.sqrt(getVariance());
    }

    /**
     * Returns a non-empty string representation of this object suitable for
     * debugging. The exact presentation format is unspecified and may vary
     * between implementations and versions.
     * 
     * @return a string representation of this object
     */
    @Override
    public String toString() {
        return String.format("%s{count=%d, sum=%f, min=%f, average=%f, max=%f, stddev=%f}",
                this.getClass().getSimpleName(), getCount(), getSum(), getMin(), getAverage(), getMax(),
                getStandardDeviation());
    }
}
