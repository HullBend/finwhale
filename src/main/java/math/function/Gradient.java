package math.function;

/**
 * Interface for the gradient of once-differentiable double-valued functions
 * over double[] arrays.
 */
public interface Gradient {
    /**
     * The first-derivative vector (a.k.a. gradient) of a double-valued function
     * over a double[] array evaluated at the input location {@code x} gets
     * stored into the output vector {@code grad}.
     * 
     * @param x
     *            a <code>double[]</code> input vector (not modified)
     * @param grad
     *            a <code>double[]</code> output vector containing the gradient
     *            at location {@code x} (modified)
     */
    void derivativeAt(double[] x, double[] grad);
}
