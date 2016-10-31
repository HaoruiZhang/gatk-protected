package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.apache.commons.lang3.mutable.MutableInt;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.Parameter;
import org.broadinstitute.hellbender.utils.mcmc.ParameterizedState;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityState extends ParameterizedState<TumorHeterogeneityParameter> {
    public static final class PopulationFractions extends ArrayList<Double> {
        //list of doubles, size = number of populations (including normal), i-th element = fraction of population i by cell number
        //normal population is last element
        private static final long serialVersionUID = 79454L;
        private final int numPopulations;
        public PopulationFractions(final List<Double> populationFractions) {
            super(Utils.nonNull(populationFractions));
            Utils.validateArg(populationFractions.size() > 1, "Number of populations must be strictly greater than 1.");
            final double populationFractionNormalization = populationFractions.stream().mapToDouble(Double::doubleValue).sum();
            Utils.validateArg(Math.abs(1. - populationFractionNormalization) <= POPULATION_FRACTION_NORMALIZATION_EPSILON,
                    "Population fractions must sum to unity.");
            numPopulations = populationFractions.size();
        }
    }

    public static final class PopulationIndicators extends ArrayList<Integer> {
        //list of integers, size = number of cells, i-th element = population index for cell i
        private static final long serialVersionUID = 81915L;
        private final int numCells;
        public PopulationIndicators(final List<Integer> populationIndicators) {
            super(Utils.nonNull(populationIndicators));
            Utils.validateArg(populationIndicators.size() > 0, "Number of cells must be positive.");
            Utils.validateArg(populationIndicators.stream().allMatch(i -> i >= 0), "Population indicators must be non-negative.");
            numCells = populationIndicators.size();
        }
    }

    public static final class VariantProfileCollection extends ArrayList<VariantProfile> {
        //list of VariantProfiles, size = number of populations (excluding normal), i-th element = VariantProfile for population i
        private static final long serialVersionUID = 76498L;
        private final int numSegments;
        private final int numVariantPopulations;
        public VariantProfileCollection(final List<VariantProfile> variantProfiles) {
            super(Utils.nonNull(variantProfiles));
            Utils.validateArg(variantProfiles.size() > 0, "Number of variants must be positive.");
            final int numSegmentsForFirstVariant = variantProfiles.get(0).numSegments;
            Utils.validateArg(variantProfiles.stream().map(s -> s.numSegments).allMatch(n -> n == numSegmentsForFirstVariant),
                    "Number of segments must be same for all variants.");
            Utils.validateArg(numSegmentsForFirstVariant > 0, "Number of segments must be positive.");
            numSegments = numSegmentsForFirstVariant;
            numVariantPopulations = variantProfiles.size();
        }
    }

    private static final double POPULATION_FRACTION_NORMALIZATION_EPSILON = 1E-4;

    private final int numPopulations;   //variant populations + normal population
    private final int numCells;
    private final int numSegments;
    private final int numVariantPloidyStates;
    private final List<MutableInt> populationCounts;

    private final TumorHeterogeneityPriorCollection priors;

    public TumorHeterogeneityState(final double concentration,
                                   final PopulationFractions populationFractions,
                                   final PopulationIndicators populationIndicators,
                                   final VariantProfileCollection variantProfileCollection,
                                   final TumorHeterogeneityPriorCollection priors) {
        super(Arrays.asList(
                new Parameter<>(TumorHeterogeneityParameter.CONCENTRATION, concentration),
                new Parameter<>(TumorHeterogeneityParameter.POPULATION_FRACTIONS, populationFractions),
                new Parameter<>(TumorHeterogeneityParameter.POPULATION_INDICATORS, populationIndicators),
                new Parameter<>(TumorHeterogeneityParameter.VARIANT_PROFILES, variantProfileCollection)));
        Utils.validateArg(concentration > 0, "Concentration must be positive.");
        Utils.nonNull(populationFractions);
        Utils.nonNull(populationIndicators);
        Utils.nonNull(variantProfileCollection);
        Utils.validateArg(populationFractions.numPopulations == variantProfileCollection.numVariantPopulations + 1,
                "Number of populations must be equal to number of variant populations + 1.");
        Utils.validateArg(Collections.max(populationIndicators) < populationFractions.numPopulations,
                "Population indicators must be consistent with number of populations.");
        Utils.validateArg(variantProfileCollection.stream().map(s -> Collections.max(s.variantPloidyStateIndicators)).allMatch(i -> i < priors.variantPloidyStatePrior().numPloidyStates()),
                "Variant ploidy-state indicators must be consistent with number of variant ploidy states.");
        Utils.nonNull(priors);
        numPopulations = populationFractions.numPopulations;
        numCells = populationIndicators.numCells;
        numSegments = variantProfileCollection.numSegments;
        numVariantPloidyStates = priors.variantPloidyStatePrior().numPloidyStates();
        populationCounts = sumPopulationCounts();
        this.priors = priors;
    }

    /*===============================================================================================================*
     * GETTERS                                                                                                       *
     *===============================================================================================================*/

    public int numPopulations() {
        return numPopulations;
    }

    public int numCells() {
        return numCells;
    }

    public int numSegments() {
        return numSegments;
    }

    public double concentration() {
        return get(TumorHeterogeneityParameter.CONCENTRATION, Double.class);
    }

    public double populationFraction(final int populationIndex) {
        validatePopulationIndex(populationIndex);
        return get(TumorHeterogeneityParameter.POPULATION_FRACTIONS, PopulationFractions.class).get(populationIndex);
    }

    public int populationCount(final int populationIndex) {
        validatePopulationIndex(populationIndex);
        return populationCounts.get(populationIndex).intValue();
    }

    public int populationIndex(final int cellIndex) {
        validateCellIndex(cellIndex);
        return get(TumorHeterogeneityParameter.POPULATION_INDICATORS, PopulationIndicators.class).get(cellIndex);
    }

    public double variantSegmentFraction(final int populationIndex) {
        validatePopulationIndex(populationIndex);
        if (populationIndex == numPopulations - 1) {
            throw new IllegalStateException("Attempted to get variant-segment fraction for normal population.");
        }
        return get(TumorHeterogeneityParameter.VARIANT_PROFILES, VariantProfileCollection.class).get(populationIndex).variantSegmentFraction;
    }

    public boolean isVariant(final int populationIndex, final int segmentIndex) {
        validatePopulationIndex(populationIndex);
        validateSegmentIndex(segmentIndex);
        if (populationIndex == numPopulations - 1) {
            //last population is normal
            return false;
        }
        return get(TumorHeterogeneityParameter.VARIANT_PROFILES, VariantProfileCollection.class).get(populationIndex).variantIndicators.get(segmentIndex);
    }

    public int variantPloidyStateIndex(final int populationIndex, final int segmentIndex) {
        validatePopulationIndex(populationIndex);
        validateSegmentIndex(segmentIndex);
        if (populationIndex == numPopulations - 1) {
            throw new IllegalStateException("Attempted to get variant-ploidy-state index for normal population.");
        }
        return get(TumorHeterogeneityParameter.VARIANT_PROFILES, VariantProfileCollection.class).get(populationIndex).variantPloidyStateIndicators.get(segmentIndex);
    }

    public TumorHeterogeneityPriorCollection priors() {
        return priors;
    }

    /*===============================================================================================================*
     * METHODS                                                                                                       *
     *===============================================================================================================*/

    public double calculatePopulationWeightedTotalCopyNumber(final TumorHeterogeneityData data, final int segmentIndex) {
        return calculatePopulationWeightedCopyNumber(data, segmentIndex, PloidyState::total, false);
    }

    public double calculatePopulationWeightedMAlleleCopyNumber(final TumorHeterogeneityData data, final int segmentIndex) {
        return calculatePopulationWeightedCopyNumber(data, segmentIndex, PloidyState::m, false);
    }

    public double calculatePopulationWeightedNAlleleCopyNumber(final TumorHeterogeneityData data, final int segmentIndex) {
        return calculatePopulationWeightedCopyNumber(data, segmentIndex, PloidyState::n, false);
    }

    public double calculatePopulationWeightedGenomicAveragedPloidy(final TumorHeterogeneityData data) {
        return IntStream.range(0, numSegments).mapToDouble(i -> calculatePopulationWeightedCopyNumber(data, i, PloidyState::total, true)).sum();
    }

    private double calculatePopulationWeightedCopyNumber(final TumorHeterogeneityData data, final int segmentIndex,
                                                         final Function<PloidyState, Integer> copyNumber,
                                                         final boolean isWeightedByFractionalSegmentLength) {
        Utils.nonNull(data);
        Utils.validateArg(data.numSegments() == numSegments,
                "Tumor-heterogeneity state and data collection must have same number of segments.");
        validateSegmentIndex(segmentIndex);

        final double numCopiesInSegment = IntStream.range(0, numPopulations)
                .mapToDouble(populationIndex ->
                        isVariant(populationIndex, segmentIndex) ?
                        populationCount(populationIndex) * copyNumber.apply(priors.variantPloidyStatePrior().ploidyStates().get(variantPloidyStateIndex(populationIndex, segmentIndex))) :
                        populationCount(populationIndex) * copyNumber.apply(priors.normalPloidyState()))
                .sum();
        final double populationWeightedCopyNumber = numCopiesInSegment / numCells;
        return isWeightedByFractionalSegmentLength ?
                populationWeightedCopyNumber * data.segmentLength(segmentIndex) / data.totalLength() :
                populationWeightedCopyNumber;
    }

    /*===============================================================================================================*
     * SETTERS (SHOULD ONLY BE USED BY SAMPLERS TO MODIFY STATE FOR SAMPLING OF INDICATOR VARIABLES                  *
     *===============================================================================================================*/

    void setPopulationIndex(final int cellIndex, final int populationIndex) {
        validateCellIndex(cellIndex);
        validatePopulationIndex(populationIndex);
        final int oldPopulationIndex = get(TumorHeterogeneityParameter.POPULATION_INDICATORS, PopulationIndicators.class).get(cellIndex);
        populationCounts.get(oldPopulationIndex).decrement();
        populationCounts.get(populationIndex).increment();
        get(TumorHeterogeneityParameter.POPULATION_INDICATORS, PopulationIndicators.class).set(cellIndex, populationIndex);
    }

    void setIsVariant(final int populationIndex, final int segmentIndex, final boolean isVariant) {
        validatePopulationIndex(populationIndex);
        validateSegmentIndex(segmentIndex);
        if (populationIndex == numPopulations - 1) {
            throw new IllegalStateException("Attempted to set variant-segment fraction for normal population.");
        }
        get(TumorHeterogeneityParameter.VARIANT_PROFILES, VariantProfileCollection.class).get(populationIndex).variantIndicators.set(segmentIndex, isVariant);
    }

    void setVariantPloidyStateIndex(final int populationIndex, final int segmentIndex, final int variantPloidyStateIndex) {
        validatePopulationIndex(populationIndex);
        validateSegmentIndex(segmentIndex);
        validateVariantPloidyStateIndex(variantPloidyStateIndex);
        if (populationIndex == numPopulations - 1) {
            throw new IllegalStateException("Attempted to set variant-ploidy-state index for normal population.");
        }
        get(TumorHeterogeneityParameter.VARIANT_PROFILES, VariantProfileCollection.class).get(populationIndex).variantPloidyStateIndicators.set(segmentIndex, variantPloidyStateIndex);
    }

    /**
     * For each variant population, represents variant-segment fraction and per-segment variant and ploidy indicators.
     */
    static class VariantProfile {
        public static final class VariantIndicators extends ArrayList<Boolean> {
            //list of booleans, size = number of segments, i-th element = true if segment i is variant
            private static final long serialVersionUID = 35746L;
            private final int numSegments;
            public VariantIndicators(final List<Boolean> variantIndicators) {
                super(Utils.nonNull(variantIndicators));
                Utils.validateArg(variantIndicators.size() > 0, "Number of segments must be positive.");
                numSegments = variantIndicators.size();
            }
        }

        public static final class VariantPloidyStateIndicators extends ArrayList<Integer> {
            //list of integers, size = number of segments, i-th element = variant-ploidy-state index of segment i
            private static final long serialVersionUID = 78476L;
            private final int numSegments;
            public VariantPloidyStateIndicators(final List<Integer> variantPloidyStateIndicators) {
                super(Utils.nonNull(variantPloidyStateIndicators));
                Utils.validateArg(variantPloidyStateIndicators.size() > 0, "Number of segments must be positive.");
                Utils.validateArg(variantPloidyStateIndicators.stream().allMatch(i -> i >= 0), "Population indicators must be non-negative.");
                numSegments = variantPloidyStateIndicators.size();
            }
        }

        private final int numSegments;
        private final double variantSegmentFraction;
        private final VariantIndicators variantIndicators;
        private final VariantPloidyStateIndicators variantPloidyStateIndicators;

        VariantProfile(final double variantSegmentFraction,
                       final VariantIndicators variantIndicators,
                       final VariantPloidyStateIndicators variantPloidyStateIndicators) {
            Utils.nonNull(variantIndicators);
            Utils.nonNull(variantPloidyStateIndicators);
            Utils.validateArg(0. <= variantSegmentFraction && variantSegmentFraction <= 1.,
                    "Variant-segment fraction must be in [0, 1].");
            Utils.validateArg(variantIndicators.numSegments == variantPloidyStateIndicators.numSegments,
                    "Number of segments must be same for variant and ploidy-state indicators.");
            numSegments = variantIndicators.numSegments;
            this.variantSegmentFraction = variantSegmentFraction;
            this.variantIndicators = variantIndicators;
            this.variantPloidyStateIndicators = variantPloidyStateIndicators;
        }
    }

    private void validatePopulationIndex(final int populationIndex) {
        Utils.validateArg(0 <= populationIndex && populationIndex < numPopulations, "Population index out of range.");
    }

    private void validateCellIndex(final int cellIndex) {
        Utils.validateArg(0 <= cellIndex && cellIndex < numCells, "Cell index out of range.");
    }

    private void validateSegmentIndex(int segmentIndex) {
        Utils.validateArg(0 <= segmentIndex && segmentIndex < numSegments, "Segment index out of range.");
    }

    private void validateVariantPloidyStateIndex(final int variantPloidyStateIndex) {
        Utils.validateArg(0 <= variantPloidyStateIndex && variantPloidyStateIndex < numVariantPloidyStates, "Variant-ploidy-state index out of range.");
    }

    private List<MutableInt> sumPopulationCounts() {
        final List<MutableInt> populationCounts = IntStream.range(0, numPopulations).boxed()
                .map(j -> new MutableInt(0)).collect(Collectors.toList());
        for (int cellIndex = 0; cellIndex < numCells; cellIndex++) {
            final int populationIndex = populationIndex(cellIndex);
            populationCounts.get(populationIndex).increment();
        }
        return populationCounts;
    }
}
