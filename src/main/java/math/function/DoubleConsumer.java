package math.function;

/**
 * Represents an operation that accepts a single {@code double}-valued argument
 * and returns no result. This is the primitive type specialization of
 * {@code Consumer} for {@code double}. Unlike most other functional interfaces,
 * {@code DoubleConsumer} is expected to operate via side-effects.
 */
public interface DoubleConsumer {

    /**
     * Performs this operation on the given argument.
     * 
     * @param value
     *            the input argument
     */
    void accept(double value);
}
