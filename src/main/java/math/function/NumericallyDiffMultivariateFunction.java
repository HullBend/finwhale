package math.function;

import math.MathConsts;

/**
 * Finite difference numerical gradient calculation using the forward difference
 * approximation {@code f'(x) = (f(x + h) - f(x)) / h}.
 * <p>
 * Scaling of {@code h} is taken into account by each individual {@code h} being
 * based upon the absolute magnitude of the corresponding element in the vector
 * {@code x}.
 */
public abstract class NumericallyDiffMultivariateFunction implements
        DiffMultivariateFunction {

    protected final double diffScale;

    public NumericallyDiffMultivariateFunction() {
        this(1.5 * Math.sqrt(MathConsts.MACH_EPS));
    }

    public NumericallyDiffMultivariateFunction(double diffScale) {
        this.diffScale = diffScale;
    }

    @Override
    public final void derivativeAt(double[] x, double[] grad) {
        double fx = this.valueAt(x);

        for (int i = 0; i < x.length; ++i) {
            double xi = x[i];
            double hi = (xi != 0) ? diffScale * Math.abs(xi) : diffScale;

            double xi_plus_hi = xi + hi;

            // account for potential round-off errors
            hi = xi_plus_hi - xi;

            x[i] = xi_plus_hi;
            // new function value for advance in variable i
            double fx_plus_hi = this.valueAt(x);
            // estimated gradient component for variable i
            grad[i] = (fx_plus_hi - fx) / hi;

            // restore the old value for variable i
            x[i] = xi;
        }
    }
}
