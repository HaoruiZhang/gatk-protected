package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;
import org.apache.commons.math3.random.Well19937a;
import org.broadinstitute.hellbender.tools.exome.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.Genome;
import org.broadinstitute.hellbender.tools.exome.SegmentedModel;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;

/**
 * Provide simulated AlleleFractionData objects and stores the "truth" data from which they were created
 *
 * @author David Benjamin
 */
public final class AlleleFractionSimulatedData {
    private static final int MIN_HETS_PER_SEGMENT = 3;
    private static final int RANDOM_SEED = 13;
    private static final int SEGMENT_LENGTHS_SEED = 17;
    private static final int DEPTHS_SEED = 19;
    private static final int MINOR_FRACTION_SEED = 23;
    private static final int BIAS_SEED = 29;
    private static final Random rng = new Random(RANDOM_SEED);

    private static PoissonDistribution makePoisson(final int seed, final double mean) {
        return new PoissonDistribution(new Well19937a(seed), mean, PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);
    }

    private final AlleleFractionState trueState;
    private final SegmentedModel segmentedModel;
    private final int numSegments;

    public AlleleFractionSimulatedData(final double averageHetsPerSegment, final int numSegments,
            final double averageDepth, final double biasMean, final double biasVariance, final double outlierProbability) {
        this.numSegments = numSegments;
        final AlleleFractionState.MinorFractions minorFractions = new AlleleFractionState.MinorFractions();
        final List<AllelicCount> alleleCounts = new ArrayList<>();
        final List<SimpleInterval> segments = new ArrayList<>();

        final PoissonDistribution segmentLengthGenerator = makePoisson(SEGMENT_LENGTHS_SEED, averageHetsPerSegment);
        final PoissonDistribution readDepthGenerator = makePoisson(DEPTHS_SEED, averageDepth);
        final UniformRealDistribution minorFractionGenerator = new UniformRealDistribution(new Well19937a(MINOR_FRACTION_SEED), 0.0, 0.5);

        //translate to ApacheCommons' parametrization of the gamma distribution
        final double gammaShape = biasMean * biasMean / biasVariance;
        final double gammaScale = biasVariance / biasMean;
        final GammaDistribution biasGenerator = new GammaDistribution(new Well19937a(BIAS_SEED), gammaShape, gammaScale);

        for (int segment = 0; segment < numSegments; segment++) {
            // calculate the range of het indices for this segment
            final int numHetsInSegment = Math.max(MIN_HETS_PER_SEGMENT, segmentLengthGenerator.sample());

            final double minorFraction = minorFractionGenerator.sample();
            minorFractions.add(minorFraction);

            //we will put all the hets in this segment on chromosome <segment> at loci 1, 2, 3 etc
            segments.add(new SimpleInterval(Integer.toString(segment), 1, numHetsInSegment + 1));
            for (int het = 1; het < numHetsInSegment + 1; het++) {
                final double bias = biasGenerator.sample();

                //flip a coin to decide alt minor (alt fraction = minor fraction) or ref minor (alt fraction = 1 - minor fraction)
                final double altFraction =  rng.nextDouble() < 0.5 ? minorFraction : 1 - minorFraction;

                //the probability of an alt read is the alt fraction modified by the bias or, in the case of an outlier, random
                final double pAlt = rng.nextDouble() < outlierProbability ? rng.nextDouble() : altFraction / (altFraction + (1 - altFraction) * bias);

                final int numReads = readDepthGenerator.sample();
                final int numAltReads = new BinomialDistribution(numReads, pAlt).sample();
                final int numRefReads = numReads - numAltReads;
                alleleCounts.add(new AllelicCount(new SimpleInterval(Integer.toString(segment), het, het), numRefReads, numAltReads));
            }
        }

        segmentedModel = new SegmentedModel(segments, new Genome(new ArrayList<>(), alleleCounts, "SAMPLE"));
        trueState = new AlleleFractionState(biasMean, biasVariance, outlierProbability, minorFractions);
    };

    public AlleleFractionState getTrueState() { return trueState; }
    public SegmentedModel getSegmentedModel() { return segmentedModel; }

    public AlleleFractionStateError error(final AlleleFractionState state) {
        final double averageMinorFractionError = IntStream.range(0, numSegments)
                .mapToDouble(s -> Math.abs(trueState.minorFractionInSegment(s) - state.minorFractionInSegment(s)))
                .average().getAsDouble();
        return new AlleleFractionStateError(averageMinorFractionError, trueState.meanBias() - state.meanBias(),
                trueState.biasVariance() - state.biasVariance(), trueState.outlierProbability() - state.outlierProbability());
    }

    public static final class AlleleFractionStateError {
        public final double averageMinorFractionError;
        public final double biasMeanError;
        public final double biasVarianceError;
        public final double outlierProbabilityError;

        public AlleleFractionStateError(final double averageMinorFractionError, final double biasMeanError,
                                        final double biasVarianceError, final double outlierProbabilityError) {
            this.averageMinorFractionError = averageMinorFractionError;
            this.biasMeanError = biasMeanError;
            this.biasVarianceError = biasVarianceError;
            this.outlierProbabilityError = outlierProbabilityError;
        }
    }
}
