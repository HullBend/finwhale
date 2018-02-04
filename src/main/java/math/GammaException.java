package math;

/**
 * Runtime exception thrown from {@link Gamma} methods.
 */
public final class GammaException extends ArithmeticException {

    private static final long serialVersionUID = 6738428636131168196L;

    public GammaException(final String msg) {
        super(msg);
    }

    public GammaException(final String msg, final int location) {
        super(msg + " : loc:" + location);
    }

    public GammaException(final String msg, final String reason, final int location) {
        super(msg + " : reason = \"" + reason + "\", loc:" + location);
    }

    public GammaException(final String msg, final double value, final String reason, final int location) {
        super(msg + " : [x = " + Double.toString(value) + " ], reason = \"" + reason + "\", loc:" + location);
    }
}
