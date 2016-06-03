package org.broadinstitute.hellbender.tools.exome;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionModeller;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.mcmc.DecileCollection;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Represents an ACNV segmented model for copy ratio and allele fraction.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class ACNVModeller {
    //use 95% HPD interval to construct {@link PosteriorSummary} for segment means and minor allele fractions
    private static final double CREDIBLE_INTERVAL_ALPHA = 0.05;
    private static final DecileCollection NAN_DECILE_COLLECTION = new DecileCollection(Collections.singletonList(Double.NaN));

    public static final Logger logger = LogManager.getLogger(ACNVModeller.class);

    private SegmentedModel segmentedModel;
    private final AllelicPanelOfNormals allelicPON;
    private final List<ACNVModeledSegment> segments = new ArrayList<>();

    //similar-segment merging may leave model in a state where it is not completely fit (i.e., deciles will be unspecified)
    private boolean isModelFit;

    private final int numSamplesCopyRatio;
    private final int numBurnInCopyRatio;
    private final int numSamplesAlleleFraction;
    private final int numBurnInAlleleFraction;
    private final JavaSparkContext ctx;

    public List<ACNVModeledSegment> getACNVModeledSegments() {
        return Collections.unmodifiableList(segments);
    }

    /**
     * Constructs a copy-ratio and allele-fraction modeller for a {@link SegmentedModel},
     * specifying number of total samples and number of burn-in samples for Markov-Chain Monte Carlo model fitting.
     * An initial model fit is performed.
     *
     * @param segmentedModel            contains segments, target coverages, and SNP counts to model
     * @param numSamplesCopyRatio       number of total samples for copy-ratio model MCMC
     * @param numBurnInCopyRatio        number of burn-in samples to discard for copy-ratio model MCMC
     * @param numSamplesAlleleFraction  number of total samples for allele-fraction model MCMC
     * @param numBurnInAlleleFraction   number of burn-in samples to discard for allele-fraction model MCMC
     * @param ctx                       JavaSparkContext, used for kernel density estimation in {@link PosteriorSummary}
     */
    public ACNVModeller(final SegmentedModel segmentedModel,
                        final int numSamplesCopyRatio, final int numBurnInCopyRatio,
                        final int numSamplesAlleleFraction, final int numBurnInAlleleFraction,
                        final JavaSparkContext ctx) {
        this(segmentedModel, AllelicPanelOfNormals.EMPTY_PON, numSamplesCopyRatio, numBurnInCopyRatio, numSamplesAlleleFraction, numBurnInAlleleFraction, ctx);
    }

    /**
     * Constructs a copy-ratio and allele-fraction modeller for a {@link SegmentedModel},
     * specifying number of total samples and number of burn-in samples for Markov-Chain Monte Carlo model fitting.
     * An initial model fit is performed.
     *
     * @param segmentedModel            contains segments, target coverages, and SNP counts to model
     * @param allelicPON                allelic-bias panel of normals
     * @param numSamplesCopyRatio       number of total samples for copy-ratio model MCMC
     * @param numBurnInCopyRatio        number of burn-in samples to discard for copy-ratio model MCMC
     * @param numSamplesAlleleFraction  number of total samples for allele-fraction model MCMC
     * @param numBurnInAlleleFraction   number of burn-in samples to discard for allele-fraction model MCMC
     * @param ctx                       JavaSparkContext, used for kernel density estimation in {@link PosteriorSummary}
     */
    public ACNVModeller(final SegmentedModel segmentedModel, final AllelicPanelOfNormals allelicPON,
                        final int numSamplesCopyRatio, final int numBurnInCopyRatio,
                        final int numSamplesAlleleFraction, final int numBurnInAlleleFraction,
                        final JavaSparkContext ctx) {
        this.segmentedModel = segmentedModel;
        this.allelicPON = allelicPON;
        this.numSamplesCopyRatio = numSamplesCopyRatio;
        this.numBurnInCopyRatio = numBurnInCopyRatio;
        this.numSamplesAlleleFraction = numSamplesAlleleFraction;
        this.numBurnInAlleleFraction = numBurnInAlleleFraction;
        this.ctx = ctx;
        logger.info("Fitting initial model...");
        fitModel();
    }

    /**
     * Performs one iteration of similar-segment merging on the list of {@link ACNVModeledSegment} held internally.
     * Markov-Chain Monte Carlo model fitting is optionally performed after each iteration using the
     * number of total samples and number of burn-in samples specified at construction.
     * @param intervalThresholdSegmentMean         threshold number of credible intervals for segment-mean similarity
     * @param intervalThresholdMinorAlleleFraction threshold number of credible intervals for minor-allele-fraction similarity
     * @param doModelFit                           if true, refit MCMC model after merging
     */
    public void performSimilarSegmentMergingIteration(final double intervalThresholdSegmentMean,
                                                      final double intervalThresholdMinorAlleleFraction,
                                                      final boolean doModelFit) {
        logger.info("Number of segments before similar-segment merging iteration: " + segments.size());
        final List<ACNVModeledSegment> mergedSegments =
                SegmentMergeUtils.mergeSimilarSegments(segments, intervalThresholdSegmentMean, intervalThresholdMinorAlleleFraction);
        logger.info("Number of segments after similar-segment merging iteration: " + mergedSegments.size());
        segmentedModel = new SegmentedModel(toUnmodeledSegments(mergedSegments), segmentedModel.getGenome());
        if (doModelFit) {
            fitModel();
        } else {
            //deciles will be set to NaN
            segments.clear();
            segments.addAll(mergedSegments);
            segments.forEach(s -> s.getSegmentMeanPosteriorSummary().setDeciles(NAN_DECILE_COLLECTION));
            segments.forEach(s -> s.getMinorAlleleFractionPosteriorSummary().setDeciles(NAN_DECILE_COLLECTION));
            isModelFit = false;
        }
    }

    /**
     * Performs Markov-Chain Monte Carlo model fitting using the
     * number of total samples and number of burn-in samples pecified at construction.
     */
    public void fitModel() {
        //perform MCMC to generate posterior samples
        logger.info("Fitting copy-ratio model...");
        final ACNVCopyRatioModeller copyRatioModeller = new ACNVCopyRatioModeller(segmentedModel);
        copyRatioModeller.fitMCMC(numSamplesCopyRatio, numBurnInCopyRatio);
        logger.info("Fitting allele-fraction model...");
        final AlleleFractionModeller alleleFractionModeller = new AlleleFractionModeller(segmentedModel, allelicPON);
        alleleFractionModeller.fitMCMC(numSamplesAlleleFraction, numBurnInAlleleFraction);

        //update list of ACNVModeledSegment with new PosteriorSummaries
        segments.clear();
        final List<SimpleInterval> unmodeledSegments = segmentedModel.getSegments();
        final List<PosteriorSummary> segmentMeansPosteriorSummaries =
                copyRatioModeller.getSegmentMeansPosteriorSummaries(CREDIBLE_INTERVAL_ALPHA, ctx);
        final List<PosteriorSummary> minorAlleleFractionsPosteriorSummaries =
                alleleFractionModeller.getMinorAlleleFractionsPosteriorSummaries(CREDIBLE_INTERVAL_ALPHA, ctx);
        for (int segment = 0; segment < unmodeledSegments.size(); segment++) {
            segments.add(new ACNVModeledSegment(unmodeledSegments.get(segment),
                    segmentMeansPosteriorSummaries.get(segment), minorAlleleFractionsPosteriorSummaries.get(segment)));
        }
        isModelFit = true;
    }

    /**
     * Return the state of the model fit, which we trust is accurately recorded by {@link ACNVModeller#isModelFit}.
     * This should be true if the model is completely fit (i.e., the deciles are specified);
     * however, similar-segment merging may leave the model in a state where it is not completely fit
     * (i.e., the deciles are unspecified).
     */
    public boolean isModelFit() {
        return isModelFit;
    }

    /**
     * Writes the list of {@link ACNVModeledSegment} held internally to file.
     * See {@link SegmentUtils#writeACNVModeledSegmentFile}.
     * @param outFile   output file
     */
    public void writeACNVModeledSegmentFile(final File outFile) {
        if (!isModelFit) {
            logger.warn("Attempted to write ACNV segment file when model was not completely fit. Performing model fit now.");
            fitModel();
        }
        SegmentUtils.writeACNVModeledSegmentFile(outFile, segments, segmentedModel.getGenome());
    }

    //converts list of ACNVModeledSegments to list of SimpleIntervals
    private List<SimpleInterval> toUnmodeledSegments(final List<ACNVModeledSegment> segments) {
        final List<SimpleInterval> unmodeledSegments = new ArrayList<>();
        for (final ACNVModeledSegment segment : segments) {
            unmodeledSegments.add(segment.getInterval());
        }
        return unmodeledSegments;
    }
}
