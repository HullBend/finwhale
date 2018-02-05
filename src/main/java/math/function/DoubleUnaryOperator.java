package math.function;

/**
 * Represents an operation on a single {@code double}-valued operand that produces
 * a {@code double}-valued result.
 */
public interface DoubleUnaryOperator {
    /**
     * Applies this operator to the given operand.
     *
     * @param x the operand
     * @return the operator result
     */
    double applyAsDouble(double x);
}
