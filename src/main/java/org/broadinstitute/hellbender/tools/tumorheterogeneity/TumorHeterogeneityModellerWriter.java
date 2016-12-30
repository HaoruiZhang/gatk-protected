package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.SegmentTableColumn;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.posteriorsummary.PosteriorSummary;
import org.broadinstitute.hellbender.utils.mcmc.posteriorsummary.PosteriorSummaryUtils;
import org.broadinstitute.hellbender.utils.mcmc.posteriorsummary.PosteriorSummaryWriter;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Writes results generated by the {@link TumorHeterogeneityModeller} to file.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityModellerWriter {
    private static final Logger logger = LogManager.getLogger(TumorHeterogeneityModellerWriter.class);

    private static final double CREDIBLE_INTERVAL_ALPHA = 0.05;
    private static final String TUMOR_HETEROGENEITY_DOUBLE_FORMAT = "%4.6f";

    private static final String POPULATION_FRACTION_NAME_PREFIX = "POPULATION_FRACTION_";
    private static final String NORMAL_FRACTION_NAME = "NORMAL_FRACTION";
    private static final String POPULATION_INDEX_NAME = "POPULATION_INDEX";
    private static final String SEGMENT_INDEX_NAME = "SEGMENT_INDEX";
    private static final String LINE_SEPARATOR = System.getProperty("line.separator");

    private final TumorHeterogeneityModeller modeller;

    public TumorHeterogeneityModellerWriter(final TumorHeterogeneityModeller modeller) {
        this.modeller = Utils.nonNull(modeller);
    }

    public void writePopulationFractionAndPloidySamples(final File outputFile) {
        Utils.nonNull(outputFile);
        if (modeller.getConcentrationSamples().size() == 0) {
            throw new IllegalStateException("Cannot output modeller result before samples have been generated.");
        }
        try (final FileWriter writer = new FileWriter(outputFile)) {
            final List<PopulationMixture.PopulationFractions> populationFractionsSamples = modeller.getPopulationFractionsSamples();

            final int numPopulations = populationFractionsSamples.get(0).size();

            //column headers
            writer.write(TumorHeterogeneityParameter.CONCENTRATION.name + "\t");
            writer.write(TumorHeterogeneityParameter.COPY_RATIO_NORMALIZATION.name + "\t");
            writer.write(TumorHeterogeneityParameter.COPY_RATIO_NOISE_CONSTANT.name + "\t");
            writer.write(TumorHeterogeneityParameter.OUTLIER_PROBABILITY.name + "\t");
            writer.write(TumorHeterogeneityParameter.PLOIDY.name + "\t");
            for (int populationIndex = 0; populationIndex < numPopulations - 1; populationIndex++) {
                writer.write(POPULATION_FRACTION_NAME_PREFIX + populationIndex + "\t");
            }
            writer.write(NORMAL_FRACTION_NAME);
            writer.write(LINE_SEPARATOR);

            //rows
            for (int sampleIndex = 0; sampleIndex < populationFractionsSamples.size(); sampleIndex++) {
                final double concentration = modeller.getConcentrationSamples().get(sampleIndex);
                final double copyRatioNormalization = modeller.getCopyRatioNormalizationSamples().get(sampleIndex);
                final double copyRatioNoiseConstant = modeller.getCopyRatioNoiseConstantSamples().get(sampleIndex);
                final double outlierProbability = modeller.getOutlierProbabilitySamples().get(sampleIndex);
                final double ploidy = modeller.getPloidySamples().get(sampleIndex);

                writer.write(String.format(TUMOR_HETEROGENEITY_DOUBLE_FORMAT + "\t", concentration));
                writer.write(String.format(TUMOR_HETEROGENEITY_DOUBLE_FORMAT + "\t", copyRatioNormalization));
                writer.write(String.format(TUMOR_HETEROGENEITY_DOUBLE_FORMAT + "\t", copyRatioNoiseConstant));
                writer.write(String.format(TUMOR_HETEROGENEITY_DOUBLE_FORMAT + "\t", outlierProbability));
                writer.write(String.format(TUMOR_HETEROGENEITY_DOUBLE_FORMAT + "\t", ploidy));
                for (int populationIndex = 0; populationIndex < numPopulations - 1; populationIndex++) {
                    final double populationFraction = populationFractionsSamples.get(sampleIndex).get(populationIndex);
                    writer.write(String.format(TUMOR_HETEROGENEITY_DOUBLE_FORMAT + "\t", populationFraction));
                }
                writer.write(String.format(TUMOR_HETEROGENEITY_DOUBLE_FORMAT + LINE_SEPARATOR, populationFractionsSamples.get(sampleIndex).normalFraction()));
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
    }

    public void writeAveragedProfiles(final File outputFile,
                                      final List<Integer> sampleIndices) {
        Utils.nonNull(outputFile);
        Utils.nonNull(sampleIndices);

        final List<PopulationMixture.VariantProfileCollection> variantProfileCollectionSamples = modeller.getVariantProfileCollectionSamples();
        final List<PopulationMixture.VariantProfileCollection> selectedVariantProfileCollectionSamples =
                sampleIndices.stream().map(variantProfileCollectionSamples::get).collect(Collectors.toList());

        if (sampleIndices.size() == 0) {
            logger.warn("Sample of posterior mode was discarded with burn-in samples and no other samples were present in purity-ploidy bin. " +
                    "Adjust burn-in and total number of samples accordingly. Using posterior mode only...");
            final PopulationMixture.VariantProfileCollection variantProfileCollectionMode =
                    modeller.getPosteriorMode().populationMixture().variantProfileCollection();
            selectedVariantProfileCollectionSamples.add(variantProfileCollectionMode);
        }
        final int numSamples = selectedVariantProfileCollectionSamples.size();

        try (final FileWriter writer = new FileWriter(outputFile)) {
            final int numVariantPopulations = selectedVariantProfileCollectionSamples.get(0).numVariantPopulations();
            final TumorHeterogeneityData data = modeller.getData();

            //column headers
            writer.write(POPULATION_INDEX_NAME + "\t");
            writer.write(SEGMENT_INDEX_NAME + "\t");
            writer.write(SegmentTableColumn.CONTIG.toString() + "\t");
            writer.write(SegmentTableColumn.START.toString() + "\t");
            writer.write(SegmentTableColumn.END.toString() + "\t");
            final List<PloidyState> ploidyStates = data.priors().ploidyStatePrior().ploidyStates();
            final int numPloidyStates = ploidyStates.size();
            for (int ploidyStateIndex = 0; ploidyStateIndex < numPloidyStates - 1; ploidyStateIndex++) {
                final PloidyState ploidyState = ploidyStates.get(ploidyStateIndex);
                writer.write(ploidyState.toString() + "\t");
            }
            final PloidyState ploidyState = ploidyStates.get(numPloidyStates - 1);
            writer.write(ploidyState.toString());
            writer.write(LINE_SEPARATOR);

            //rows
            for (int populationIndex = 0; populationIndex < numVariantPopulations; populationIndex++) {
                final int pi = populationIndex;
                for (int segmentIndex = 0; segmentIndex < data.numSegments(); segmentIndex++) {
                    final SimpleInterval segment = data.segments().get(segmentIndex).getInterval();
                    writer.write(populationIndex + "\t");
                    writer.write(segmentIndex + "\t");
                    writer.write(segment.getContig() + "\t");
                    writer.write(segment.getStart() + "\t");
                    writer.write(segment.getEnd() + "\t");

                    final int si = segmentIndex;
                    for (int ploidyStateIndex = 0; ploidyStateIndex < numPloidyStates; ploidyStateIndex++) {
                        final int vpsi = ploidyStateIndex;
                        final double numPloidyStateSamples = selectedVariantProfileCollectionSamples.stream()
                                .filter(vpc -> vpc.get(pi).ploidyState(si).equals(ploidyStates.get(vpsi)))
                                .count();
                        final double ploidyStatePosteriorMean = numPloidyStateSamples / numSamples;
                        writer.write(String.format(TUMOR_HETEROGENEITY_DOUBLE_FORMAT, ploidyStatePosteriorMean));
                        if (ploidyStateIndex != numPloidyStates - 1) {
                            writer.write("\t");
                        }
                    }
                    if (!(segmentIndex == data.numSegments() - 1 && populationIndex == numVariantPopulations - 1)) {
                        writer.write(LINE_SEPARATOR);
                    }
                }
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
    }

    public void writePosteriorSummaries(final File outputFile,
                                        final List<Integer> sampleIndices,
                                        final JavaSparkContext ctx) {
        Utils.nonNull(outputFile);
        Utils.nonNull(sampleIndices);
        Utils.nonNull(ctx);

        try (final PosteriorSummaryWriter<String> writer = new PosteriorSummaryWriter<>(outputFile, TUMOR_HETEROGENEITY_DOUBLE_FORMAT)) {
            writer.writeAllRecords(getGlobalParameterPosteriorSummaries(sampleIndices, CREDIBLE_INTERVAL_ALPHA, ctx).entrySet());
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
    }

    private Map<String, PosteriorSummary> getGlobalParameterPosteriorSummaries(final List<Integer> sampleIndices,
                                                                               final double credibleIntervalAlpha,
                                                                               final JavaSparkContext ctx) {
        //get samples for all indices
        final List<Double> concentrationSamples = modeller.getConcentrationSamples();
        final List<Double> copyRatioNormalizationSamples = modeller.getCopyRatioNormalizationSamples();
        final List<Double> copyRatioNoiseConstantSamples = modeller.getCopyRatioNoiseConstantSamples();
        final List<Double> outlierProbabilitySamples = modeller.getOutlierProbabilitySamples();
        final List<PopulationMixture.PopulationFractions> populationFractionsSamples = modeller.getPopulationFractionsSamples();
        final List<Double> ploidySamples = modeller.getPloidySamples();

        //collect samples for selected indices
        final int numVariantPopulations = modeller.getVariantProfileCollectionSamples().get(0).numVariantPopulations();
        final List<Double> selectedConcentrationSamples = sampleIndices.stream().map(concentrationSamples::get).collect(Collectors.toList());
        final List<Double> selectedCopyRatioNormalizationSamples = sampleIndices.stream().map(copyRatioNormalizationSamples::get).collect(Collectors.toList());
        final List<Double> selectedCopyRatioNoiseConstantSamples = sampleIndices.stream().map(copyRatioNoiseConstantSamples::get).collect(Collectors.toList());
        final List<Double> selectedOutlierProbabilitySamples = sampleIndices.stream().map(outlierProbabilitySamples::get).collect(Collectors.toList());
        final List<List<Double>> selectedVariantPopulationFractionsSamples = IntStream.range(0, numVariantPopulations).boxed()
                .map(populationIndex -> sampleIndices.stream().map(sampleIndex -> populationFractionsSamples.get(sampleIndex).get(populationIndex)).collect(Collectors.toList()))
                .collect(Collectors.toList());
        final List<Double> selectedNormalPopulationFractionSamples = sampleIndices.stream()
                .map(sampleIndex -> populationFractionsSamples.get(sampleIndex).normalFraction()).collect(Collectors.toList());
        final List<Double> selectedPloidySamples = sampleIndices.stream().map(ploidySamples::get).collect(Collectors.toList());

        //if there are no selected indices, then use only the posterior mode and log a warning
        if (sampleIndices.size() == 0) {
            logger.warn("Sample of posterior mode was discarded with burn-in samples and no other samples were present in purity-ploidy bin. " +
                    "Adjust burn-in and total number of samples accordingly. Using posterior mode only...");
            final double concentrationMode = modeller.getPosteriorMode().concentration();
            final double copyRatioNormalizationMode = modeller.getPosteriorMode().copyRatioNormalization();
            final double copyRatioNoiseConstantMode = modeller.getPosteriorMode().copyRatioNoiseConstant();
            final double outlierProbabilityMode = modeller.getPosteriorMode().outlierProbability();
            final List<Double> variantPopulationFractionsMode = IntStream.range(0, numVariantPopulations).boxed()
                    .map(populationIndex -> modeller.getPosteriorMode().populationMixture().populationFraction(populationIndex))
                    .collect(Collectors.toList());
            final double normalPopulationFractionMode = modeller.getPosteriorMode().populationMixture().populationFractions().normalFraction();
            final double ploidyMode = modeller.getPosteriorMode().ploidy();
            selectedConcentrationSamples.add(concentrationMode);
            selectedCopyRatioNormalizationSamples.add(copyRatioNormalizationMode);
            selectedCopyRatioNoiseConstantSamples.add(copyRatioNoiseConstantMode);
            selectedOutlierProbabilitySamples.add(outlierProbabilityMode);
            IntStream.range(0, numVariantPopulations)
                    .forEach(populationIndex -> selectedVariantPopulationFractionsSamples.get(populationIndex).add(variantPopulationFractionsMode.get(populationIndex)));
            selectedNormalPopulationFractionSamples.add(normalPopulationFractionMode);
            selectedPloidySamples.add(ploidyMode);
        }

        //construct a map containing parameter name -> posterior summary for each parameter
        final Map<String, PosteriorSummary> posteriorSummaries = new LinkedHashMap<>();
        posteriorSummaries.put(TumorHeterogeneityParameter.CONCENTRATION.name,
                PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(selectedConcentrationSamples, credibleIntervalAlpha, ctx));
        posteriorSummaries.put(TumorHeterogeneityParameter.COPY_RATIO_NORMALIZATION.name,
                PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(selectedCopyRatioNormalizationSamples, credibleIntervalAlpha, ctx));
        posteriorSummaries.put(TumorHeterogeneityParameter.COPY_RATIO_NOISE_CONSTANT.name,
                PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(selectedCopyRatioNoiseConstantSamples, credibleIntervalAlpha, ctx));
        posteriorSummaries.put(TumorHeterogeneityParameter.OUTLIER_PROBABILITY.name,
                PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(selectedOutlierProbabilitySamples, credibleIntervalAlpha, ctx));
        IntStream.range(0, numVariantPopulations).forEach(populationIndex ->
                posteriorSummaries.put(POPULATION_FRACTION_NAME_PREFIX + populationIndex,
                        PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(selectedVariantPopulationFractionsSamples.get(populationIndex), credibleIntervalAlpha, ctx)));
        posteriorSummaries.put(NORMAL_FRACTION_NAME,
                PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(selectedNormalPopulationFractionSamples, credibleIntervalAlpha, ctx));
        posteriorSummaries.put(TumorHeterogeneityParameter.PLOIDY.name,
                PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(selectedPloidySamples, credibleIntervalAlpha, ctx));
        return posteriorSummaries;
    }
}
