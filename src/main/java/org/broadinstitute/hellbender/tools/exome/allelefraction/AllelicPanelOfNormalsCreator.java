package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.tools.exome.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.AllelicCountCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Creates the panel of normals used for allele-bias correction.  See docs/CNVs/CNV-methods.pdf.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AllelicPanelOfNormalsCreator {
    private static final Logger logger = LogManager.getLogger(AllelicPanelOfNormalsCreator.class);

    private final JavaSparkContext ctx;
    private final List<File> pulldownFiles;

    public AllelicPanelOfNormalsCreator(final JavaSparkContext ctx, final List<File> pulldownFiles) {
        Utils.nonNull(ctx, "JavaSparkContext cannot be null when creating allelic panel of normals.");
        Utils.nonNull(pulldownFiles, "List of pulldown files cannot be null.");
        pulldownFiles.stream().forEach(Utils::regularReadableUserFile);

        this.ctx = ctx;
        this.pulldownFiles = new ArrayList<>(pulldownFiles);
    }

    public AllelicPanelOfNormals create(final double siteFrequency) {
        logger.info("Creating allelic panel of normals...");
        final Map<SimpleInterval, MutableInt> numberOfSamplesMap = new HashMap<>();
        final Map<SimpleInterval, AllelicCount> totalCountsMap = new HashMap<>();
        int pulldownFileCounter = 1;
        final int totalNumberOfSamples = pulldownFiles.size();
        for (final File pulldownFile : pulldownFiles) {
            logger.info("Processing pulldown file " + pulldownFileCounter++ + "/" + totalNumberOfSamples + " (" + pulldownFile + ")...");
            final AllelicCountCollection pulldownCounts = new AllelicCountCollection(pulldownFile);
            for (final AllelicCount count : pulldownCounts.getCounts()) {
                final SimpleInterval site = count.getInterval();
                final AllelicCount currentCountAtSite = totalCountsMap.getOrDefault(site, new AllelicCount(site, 0, 0));
                final AllelicCount updatedCountAtSite = new AllelicCount(
                        site,
                        currentCountAtSite.getRefReadCount() + count.getRefReadCount(),
                        currentCountAtSite.getAltReadCount() + count.getAltReadCount());
                totalCountsMap.put(site, updatedCountAtSite);
                final MutableInt numberOfSamplesAtSite = numberOfSamplesMap.get(site);
                if (numberOfSamplesAtSite == null) {
                    numberOfSamplesMap.put(site, new MutableInt(1));
                } else {
                    numberOfSamplesAtSite.increment();
                }
            }
        }

        logger.info("Total number of unique sites present in samples: " + totalCountsMap.size());
        final AllelicCountCollection totalCounts = new AllelicCountCollection();
        numberOfSamplesMap.entrySet().stream()
                .filter(e -> e.getValue().doubleValue() / totalNumberOfSamples >= siteFrequency)
                .map(e -> totalCountsMap.get(e.getKey()))
                .forEach(totalCounts::add);
        logger.info(String.format("Number of unique sites present in samples above site frequency = %4.3f: %d", siteFrequency, totalCounts.getCounts().size()));
        return new AllelicPanelOfNormals(totalCounts);
    }
}
