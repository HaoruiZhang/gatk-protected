package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.*;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityModeller {
    private static final double EPSILON = 1E-10;

    private static final double CONCENTRATION_MIN = EPSILON;
    private static final double CONCENTRATION_MAX = 10.;
    private static final double CONCENTRATION_SLICE_SAMPLING_WIDTH = 0.005;

    private final ParameterizedModel<TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> model;

    private final List<Double> concentrationSamples = new ArrayList<>();
    private final List<TumorHeterogeneityState.PopulationFractions> populationFractionsSamples = new ArrayList<>();
    private final List<TumorHeterogeneityState.PopulationStates> populationStatesSamples = new ArrayList<>();

    /**
     */
    public TumorHeterogeneityModeller(final List<ACNVModeledSegment> segments,
                                      final int numPopulations,
                                      final int numCells,
                                      final double concentrationPriorAlpha,
                                      final double concentrationPriorBeta,
                                      final double variantSegmentFractionPriorAlpha,
                                      final double variantSegmentFractionPriorBeta,
                                      final RandomGenerator rng) {
        Utils.nonNull(segments);
        Utils.validateArg(numPopulations > 0, "Maximum number of populations must be positive.");
        Utils.validateArg(numCells > 0, "Number of auxiliary cells must be positive.");

        //create TumorHeterogeneityData from ACNV segments
        final TumorHeterogeneityData data = new TumorHeterogeneityData(segments);

        //initialize population fractions to be evenly distributed
        final TumorHeterogeneityState.PopulationFractions initialPopulationFractions =
                new TumorHeterogeneityState.PopulationFractions(Collections.nCopies(numPopulations, 1. / numPopulations));
        //randomly initialize population indicators for each cell
        final TumorHeterogeneityState.PopulationIndicators initialPopulationIndicators =
                initializePopulationIndicators(numPopulations, numCells, rng);
        //initialize population states to non-variant
        final TumorHeterogeneityState.PopulationStates initialPopulationStates =
                initializePopulationStates(numPopulations, data.numSegments());

        //initialize TumorHeterogeneityState
        final double initialConcentration = concentrationPriorAlpha / concentrationPriorBeta;
        final TumorHeterogeneityState.HyperparameterValues initialVariantSegmentFractionHyperparameters =
                new TumorHeterogeneityState.HyperparameterValues(variantSegmentFractionPriorAlpha, variantSegmentFractionPriorBeta);
        final TumorHeterogeneityState initialState = new TumorHeterogeneityState(
                initialConcentration, initialPopulationFractions, initialPopulationIndicators,
                initialVariantSegmentFractionHyperparameters, initialPopulationStates);

        //define samplers
        final TumorHeterogeneitySamplers.ConcentrationSampler concentrationSampler =
                new TumorHeterogeneitySamplers.ConcentrationSampler(CONCENTRATION_MIN, CONCENTRATION_MAX,
                        CONCENTRATION_SLICE_SAMPLING_WIDTH, concentrationPriorAlpha, concentrationPriorBeta);
        final TumorHeterogeneitySamplers.PopulationFractionsSampler populationFractionsSampler =
                new TumorHeterogeneitySamplers.PopulationFractionsSampler();
        final TumorHeterogeneitySamplers.PopulationIndicatorsSampler populationIndicatorsSampler =
                new TumorHeterogeneitySamplers.PopulationIndicatorsSampler();
        final TumorHeterogeneitySamplers.VariantSegmentFractionHyperparametersSampler variantSegmentFractionHyperparametersSampler =
                new TumorHeterogeneitySamplers.VariantSegmentFractionHyperparametersSampler();
        final TumorHeterogeneitySamplers.PopulationStatesSampler populationStatesSampler =
                new TumorHeterogeneitySamplers.PopulationStatesSampler();


        model = new ParameterizedModel.GibbsBuilder<>(initialState, data)
                .addParameterSampler(TumorHeterogeneityParameter.CONCENTRATION, concentrationSampler, Double.class)
                .addParameterSampler(TumorHeterogeneityParameter.POPULATION_INDICATORS, populationIndicatorsSampler, TumorHeterogeneityState.PopulationIndicators.class)
                .addParameterSampler(TumorHeterogeneityParameter.POPULATION_FRACTIONS, populationFractionsSampler, TumorHeterogeneityState.PopulationFractions.class)
                .addParameterSampler(TumorHeterogeneityParameter.VARIANT_SEGMENT_FRACTION_HYPERPARAMETERS, variantSegmentFractionHyperparametersSampler, TumorHeterogeneityState.HyperparameterValues.class)
                .addParameterSampler(TumorHeterogeneityParameter.POPULATION_STATES, populationStatesSampler, TumorHeterogeneityState.PopulationStates.class)
                .build();
    }

    /**
     * Adds {@code numSamples - numBurnIn} Markov-Chain Monte-Carlo samples of the parameter posteriors (generated using
     * Gibbs sampling) to the collections held internally.  The current {@link TumorHeterogeneityState} held internally is used
     * to initialize the Markov Chain.
     * @param numSamples    total number of samples per posterior
     * @param numBurnIn     number of burn-in samples to discard
     */
    public void fitMCMC(final int numSamples, final int numBurnIn) {
        //run MCMC
        final GibbsSampler<TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> gibbsSampler
                = new GibbsSampler<>(numSamples, model);
        gibbsSampler.runMCMC();
        //update posterior samples
        concentrationSamples.addAll(gibbsSampler.getSamples(TumorHeterogeneityParameter.CONCENTRATION,
                Double.class, numBurnIn));
        populationFractionsSamples.addAll(gibbsSampler.getSamples(TumorHeterogeneityParameter.POPULATION_FRACTIONS,
                TumorHeterogeneityState.PopulationFractions.class, numBurnIn));
        populationStatesSamples.addAll(gibbsSampler.getSamples(TumorHeterogeneityParameter.POPULATION_STATES,
                TumorHeterogeneityState.PopulationStates.class, numBurnIn));
    }

    /**
     * Returns an unmodifiable view of the list of samples of the concentration posterior.
     * @return  unmodifiable view of the list of samples of the concentration posterior
     */
    public List<Double> getConcentrationSamples() {
        return Collections.unmodifiableList(concentrationSamples);
    }

    /**
     * Returns an unmodifiable view of the list of samples of the population-fractions posterior, represented as a list of
     * {@link TumorHeterogeneityState.PopulationFractions} objects.
     * @return  unmodifiable view of the list of samples of the population-fractions posterior
     */
    public List<TumorHeterogeneityState.PopulationFractions> getPopulationFractionsSamples() {
        return Collections.unmodifiableList(populationFractionsSamples);
    }

    /**
     * Returns an unmodifiable view of the list of samples of the population-states posterior, represented as a list of
     * {@link TumorHeterogeneityState.PopulationFractions} objects.
     * @return  unmodifiable view of the list of samples of the population-states posterior
     */
    public List<TumorHeterogeneityState.PopulationStates> getPopulationStatesSamples() {
        return Collections.unmodifiableList(populationStatesSamples);
    }

    /**
     * Returns a Map of {@link PosteriorSummary} elements summarizing the global parameters.
     * Should only be called after {@link TumorHeterogeneityModeller#fitMCMC(int, int)} has been called.
     * @param credibleIntervalAlpha credible-interval alpha, must be in (0, 1)
     * @param ctx                   {@link JavaSparkContext} used for mllib kernel density estimation
     * @return                      list of {@link PosteriorSummary} elements summarizing the global parameters
     */
    public Map<TumorHeterogeneityParameter, PosteriorSummary> getGlobalParameterPosteriorSummaries(final double credibleIntervalAlpha, final JavaSparkContext ctx) {
        final Map<TumorHeterogeneityParameter, PosteriorSummary> posteriorSummaries = new LinkedHashMap<>();
        posteriorSummaries.put(TumorHeterogeneityParameter.CONCENTRATION, PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(concentrationSamples, credibleIntervalAlpha, ctx));
        return posteriorSummaries;
    }

    private TumorHeterogeneityState.PopulationIndicators initializePopulationIndicators(final int numPopulations,
                                                                                        final int numCells,
                                                                                        final RandomGenerator rng) {
        final List<Integer> populationIndices = IntStream.range(0, numPopulations).boxed().collect(Collectors.toList());
        final Function<Integer, Double> probabilityFunction = j -> 1. / numPopulations;
        return new TumorHeterogeneityState.PopulationIndicators(IntStream.range(0, numCells).boxed()
                .map(p -> GATKProtectedMathUtils.randomSelect(populationIndices, probabilityFunction, rng))
                .collect(Collectors.toList()));
    }

    private TumorHeterogeneityState.PopulationStates initializePopulationStates(final int numPopulations,
                                                                                final int numSegments) {
        return new TumorHeterogeneityState.PopulationStates(Collections.nCopies(numPopulations, initializePopulationState(numSegments)));
    }

    private TumorHeterogeneityState.PopulationState initializePopulationState(final int numSegments) {
        final double variantSegmentFraction = 0.;
        final TumorHeterogeneityState.PopulationState.VariantIndicators variantIndicators =
                new TumorHeterogeneityState.PopulationState.VariantIndicators(Collections.nCopies(numSegments, false));
        final TumorHeterogeneityState.PopulationState.VariantPloidyStateIndicators variantPloidyStateIndicators =
                new TumorHeterogeneityState.PopulationState.VariantPloidyStateIndicators(Collections.nCopies(numSegments, 0));
        return new TumorHeterogeneityState.PopulationState(variantSegmentFraction, variantIndicators, variantPloidyStateIndicators);
    }
}
