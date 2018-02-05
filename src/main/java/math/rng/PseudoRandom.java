package math.rng;

public interface PseudoRandom {

    long nextLong();

    // TODO: explain: is this [0, 1] or [0, 1)? { -> rather [0, 1)}
    double nextDouble();

    double nextGaussian();

    float nextFloat();

    int nextInt();

    // TODO: testing!
    void nextBytes(byte[] bytes);

    boolean nextBoolean();

    long nextLong(long n);

    int nextInt(int n);

    int nextInt(int min, int max);

    long nextLong(long min, long max);

    int next(int bits);

    void nextLongs(long[] longs);
}
