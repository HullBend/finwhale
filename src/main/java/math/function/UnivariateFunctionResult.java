package math.function;

public final class UnivariateFunctionResult {

    public final double point;
    public final double value;
    public final int iterations;
    public final boolean converged;

    public UnivariateFunctionResult(double point, double value, int iterations,
            boolean converged) {
        this.point = point;
        this.value = value;
        this.iterations = iterations;
        this.converged = converged;
    }

    public UnivariateFunctionResult(double point, double value, int iterations) {
        this(point, value, iterations, true);
    }

    public UnivariateFunctionResult(double point, double value) {
        this(point, value, 0);
    }
}
