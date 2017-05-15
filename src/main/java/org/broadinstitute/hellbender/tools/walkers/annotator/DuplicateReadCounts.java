package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by tsato on 5/10/17.
 */
public class DuplicateReadCounts extends GenotypeAnnotation implements StandardSomaticAnnotation {
    public static final String UNIQUE_ALT_READ_SET_COUNT_KEY = "UNIQ_ALT_READ_SET_COUNT";
    public static final String DUPLICATE_READ_COUNT = "DUP_READ_COUNT";

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(UNIQUE_ALT_READ_SET_COUNT_KEY);
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Arrays.asList(
                new VCFFormatHeaderLine(UNIQUE_ALT_READ_SET_COUNT_KEY, 1, VCFHeaderLineType.Integer, "Number of reads with unique start and end positions at a variant site"),
                new VCFFormatHeaderLine(DUPLICATE_READ_COUNT, 3, VCFHeaderLineType.Integer, "Counts of reads in each duplicate set")); // replace 3 with VCFHeaderLineCount.UNBOUNDED
    }

    @Override
    public void annotate(final ReferenceContext ref,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final ReadLikelihoods<Allele> likelihoods) {
        if (g.isHomRef()){
            // skip the normal sample
            return;
        }

        final Allele altAllele = vc.getAlternateAllele(0); // assume single-allelic
        final String tumorSampleName = g.getSampleName();
        Collection<ReadLikelihoods<Allele>.BestAllele> tumorBestAlleles = likelihoods.bestAlleles(tumorSampleName);
        Map<ImmutablePair<Integer, Integer>, Long> duplicateReadMap = tumorBestAlleles.stream()
                .filter(ba -> ba.allele.equals(altAllele) && ba.isInformative())
                .map(ba -> new ImmutablePair<Integer, Integer>(ba.read.getStart(), ba.read.getEnd())) // build fails without the <Integer, Integer>...why?
                .collect(Collectors.groupingBy(x -> x, Collectors.counting()));

        final List<Long> duplicateCounts = new ArrayList<>(duplicateReadMap.values());
        Collections.reverse(duplicateCounts); // sort duplicateCounts in descending order
        final int topN = 3;
        final List<Long> topNDuplicateCounts = duplicateCounts.size() <= 3 ? duplicateCounts
                : duplicateCounts.subList(0, topN);

        gb.attribute(UNIQUE_ALT_READ_SET_COUNT_KEY, duplicateReadMap.size());
        gb.attribute(DUPLICATE_READ_COUNT, topNDuplicateCounts);
    }


}
