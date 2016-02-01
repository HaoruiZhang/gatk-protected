package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Log;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Collections;
import java.util.List;

/**
 * Unit tests for {@link ACNVModeller}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class ACNVModellerUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/";

    private static final String SAMPLE_NAME = "test";
    private static final File COVERAGES_FILE = new File(TEST_SUB_DIR
            + "coverages-for-similar-segment-merging.tsv");
    private static final File SNP_COUNTS_FILE = new File(TEST_SUB_DIR
            + "snps-for-similar-segment-merging.tsv");
    private static final File SEGMENT_FILE =
            new File(TEST_SUB_DIR + "segments-for-similar-segment-merging.seg");
    private static final File SEGMENTS_TRUTH_FILE = new File(TEST_SUB_DIR
            + "segments-truth-for-similar-segment-merging.seg");

    private static final int NUM_SAMPLES = 100;
    private static final int NUM_BURN_IN = 50;

    private static final double INTERVAL_THRESHOLD = 1.; //merge when parameter modes are within 95% HPD interval

    /**
     * Test of similar-segment merging using only copy-ratio data (simulated coverages and segments).
     * Spurious breakpoints have been introduced into the list of true segments; similar-segment merging should
     * remerge segments broken by these breakpoints and reproduce the original list of true segments.
     */
    @Test
    public void testMergeSimilarSegmentsCopyRatio() {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);

        final String tempDir = publicTestDir + "similar-segment-copy-ratio-test";
        final File tempDirFile = createTempDir(tempDir);

        //load data (coverages and segments)
        final List<TargetCoverage> targetCoverages = TargetCoverageUtils.readTargetsWithCoverage(COVERAGES_FILE);
        final List<AllelicCount> snpCountsDummy =
                Collections.singletonList(new AllelicCount(new SimpleInterval("1", 1, 1), 0, 1));
        final Genome genome = new Genome(targetCoverages, snpCountsDummy, SAMPLE_NAME);
        final SegmentedModel segmentedModel = new SegmentedModel(SEGMENT_FILE, genome);

        //run MCMC
        final ACNVModeller modeller =
                new ACNVModeller(segmentedModel, tempDirFile.getAbsolutePath() + "/test",
                        NUM_SAMPLES, NUM_BURN_IN, 10, 0, ctx);
        modeller.mergeSimilarSegments(INTERVAL_THRESHOLD, Double.POSITIVE_INFINITY,
                NUM_SAMPLES, NUM_BURN_IN, 10, 0);

        //check equality of segments
        final List<SimpleInterval> segmentsResult =
                SegmentUtils.readIntervalsFromSegmentFile(new File(tempDirFile.getAbsolutePath() +
                        "/test-" + ACNVModeller.FINAL_SEG_FILE_TAG + ".seg"));
        final List<SimpleInterval> segmentsTruth =
                SegmentUtils.readIntervalsFromSegmentFile(SEGMENTS_TRUTH_FILE);
        Assert.assertEquals(segmentsResult, segmentsTruth);
    }

    /**
     * Test of similar-segment merging using simulated data (coverages, SNP counts, and segments).
     * Spurious breakpoints have been introduced into the list of true segments; similar-segment merging should
     * remerge segments broken by these breakpoints and reproduce the original list of true segments.
     */
    @Test
    public void testMergeSimilarSegments() {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);

        final String tempDir = publicTestDir + "similar-segment-test";
        final File tempDirFile = createTempDir(tempDir);

        //load data (coverages, SNP counts, and segments)
        final Genome genome = new Genome(COVERAGES_FILE, SNP_COUNTS_FILE, SAMPLE_NAME);
        final SegmentedModel segmentedModel = new SegmentedModel(SEGMENT_FILE, genome);

        //run MCMC
        final ACNVModeller modeller =
                new ACNVModeller(segmentedModel, tempDirFile.getAbsolutePath() + "/test",
                        NUM_SAMPLES, NUM_BURN_IN, NUM_SAMPLES, NUM_BURN_IN, ctx);
        modeller.mergeSimilarSegments(INTERVAL_THRESHOLD, INTERVAL_THRESHOLD,
                NUM_SAMPLES, NUM_BURN_IN, NUM_SAMPLES, NUM_BURN_IN);

        //check equality of segments
        final List<SimpleInterval> segmentsResult =
                SegmentUtils.readIntervalsFromSegmentFile(new File(tempDirFile.getAbsolutePath() +
                        "/test-" + ACNVModeller.FINAL_SEG_FILE_TAG + ".seg"));
        final List<SimpleInterval> segmentsTruth =
                SegmentUtils.readIntervalsFromSegmentFile(SEGMENTS_TRUTH_FILE);
        Assert.assertEquals(segmentsResult, segmentsTruth);
    }
}