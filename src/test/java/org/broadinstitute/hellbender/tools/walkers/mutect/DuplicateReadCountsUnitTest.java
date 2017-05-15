package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.DuplicateReadCounts;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.internal.junit.ArrayAsserts;

import java.io.IOException;
import java.util.*;

/**
 * Created by tsato on 5/15/17.
 */
public class DuplicateReadCountsUnitTest {
    final int numChromosomes = 2;
    final int startingChromosome = 1;
    final int chromosomeSize = 1000;
    final SAMFileHeader SAM_HEADER = ArtificialReadUtils.createArtificialSamHeader(numChromosomes, startingChromosome, chromosomeSize);

    final int chromosomeIndex = 1;
    final long variantSite = 105;
    final int alignmentStart = 100;
    final int readLength = 9;

    final List<Allele> alleles = Arrays.asList(Allele.create((byte) 'C', true), Allele.create((byte) 'A', false));
    final AlleleList<Allele> alleleList = new IndexedAlleleList<>(alleles);
    final VariantContext vc = new VariantContextBuilder("source", Integer.toString(chromosomeIndex), variantSite, variantSite, alleles)
            .attribute(GATKVCFConstants.TUMOR_LOD_KEY, new double[]{ 6.4 }).make();

    final String sampleName = "SamThreetree";
    final SampleList sampleList = new IndexedSampleList(sampleName);

    final byte baseq = 30;
    final int mapq = 60;


    @Test
    public void testDuplicates() throws IOException {
        final int numAltReads = 10;
        final int numRefReads = 10;
        final int numReads = numAltReads + numRefReads;

        final Map<String, List<GATKRead>> readMap = new LinkedHashMap<>();
        final List<GATKRead> reads = new ArrayList<>(numReads);
        final byte[] quals = new byte[readLength];
        Arrays.fill(quals, baseq);

        for (int i = 0; i < numAltReads; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(SAM_HEADER, "Read" + i, chromosomeIndex, alignmentStart,
                    "CCCCACCCC".getBytes(), quals, "9M");
            reads.add(read);
        }

        for (int j = numAltReads; j < numReads; j++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(SAM_HEADER, "Read" + j, chromosomeIndex, alignmentStart,
                    "CCCCCCCCC".getBytes(), quals, "9M");
            reads.add(read);
        }

        readMap.put(sampleName, reads);

        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readMap);
        final int sampleIndex = 0;
        final LikelihoodMatrix<Allele> matrix = likelihoods.sampleMatrix(sampleIndex);

        final double logLikelihoodOfBestAllele = 10.0;
        final int refAlleleIndex = 0;
        final int altAlleleIndex = 1;

        // Log likelihoods are initialized to 0 for all alleles. For each read_i, set its log likelihood for ALT allele to a positive value.
        // This makes read_i an ALT read
        for (int i = 0; i < numAltReads; i++) {
            matrix.set(altAlleleIndex, i, logLikelihoodOfBestAllele);
        }

        // Analogously, make read_j a REF read
        for (int j = numAltReads; j < numReads; j++) {
            matrix.set(refAlleleIndex, j, logLikelihoodOfBestAllele);
        }

        final DuplicateReadCounts duplicateReadCountsAnnotation = new DuplicateReadCounts();
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sampleName);
        duplicateReadCountsAnnotation.annotate(null, vc, genotypeBuilder.make(), genotypeBuilder, likelihoods);

        final Genotype genotype = genotypeBuilder.make();

        final int uniqueReadSetCount = GATKProtectedVariantContextUtils.getAttributeAsInt(genotype, DuplicateReadCounts.UNIQUE_ALT_READ_SET_COUNT_KEY, -1);
        final int[] duplicateReadCount = GATKProtectedVariantContextUtils.getAttributeAsIntArray(genotype, DuplicateReadCounts.DUPLICATE_READ_COUNT, () -> null, -1);

        Assert.assertEquals(uniqueReadSetCount, 1);
        ArrayAsserts.assertArrayEquals(duplicateReadCount, new int[]{ 10 });
    }


    @Test
    public void testMultipleDuplicateSets() throws IOException {
        // TODO: REFACTOR code
        final int numAltReads = 10;
        final int numRefReads = 10;
        final int numReads = numAltReads + numRefReads;

        final Map<String, List<GATKRead>> readMap = new LinkedHashMap<>();
        final List<GATKRead> reads = new ArrayList<>(numReads);
        final byte[] quals = new byte[readLength];
        Arrays.fill(quals, baseq);

        for (int i = 0; i < numAltReads; i++) {
            final int startPosition = i % 2 == 0 ? alignmentStart : alignmentStart + 1;
            final GATKRead read = ArtificialReadUtils.createArtificialRead(SAM_HEADER, "Read" + i, chromosomeIndex, startPosition,
                        "CCCCACCCC".getBytes(), quals, "9M");
            reads.add(read);
        }

        for (int j = numAltReads; j < numReads; j++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(SAM_HEADER, "Read" + j, chromosomeIndex, alignmentStart,
                    "CCCCCCCCC".getBytes(), quals, "9M");
            reads.add(read);
        }

        readMap.put(sampleName, reads);

        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readMap);
        final int sampleIndex = 0;
        final LikelihoodMatrix<Allele> matrix = likelihoods.sampleMatrix(sampleIndex);

        final double logLikelihoodOfBestAllele = 10.0;
        final int refAlleleIndex = 0;
        final int altAlleleIndex = 1;

        // Log likelihoods are initialized to 0 for all alleles. For each read_i, set its log likelihood for ALT allele to a positive value.
        // This makes read_i an ALT read
        for (int i = 0; i < numAltReads; i++) {
            matrix.set(altAlleleIndex, i, logLikelihoodOfBestAllele);
        }

        // Analogously, make read_j a REF read
        for (int j = numAltReads; j < numReads; j++) {
            matrix.set(refAlleleIndex, j, logLikelihoodOfBestAllele);
        }

        final DuplicateReadCounts duplicateReadCountsAnnotation = new DuplicateReadCounts();
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sampleName);
        duplicateReadCountsAnnotation.annotate(null, vc, genotypeBuilder.make(), genotypeBuilder, likelihoods);

        final Genotype genotype = genotypeBuilder.make();

        final int uniqueReadSetCount = GATKProtectedVariantContextUtils.getAttributeAsInt(genotype, DuplicateReadCounts.UNIQUE_ALT_READ_SET_COUNT_KEY, -1);
        final int[] duplicateReadCount = GATKProtectedVariantContextUtils.getAttributeAsIntArray(genotype, DuplicateReadCounts.DUPLICATE_READ_COUNT, () -> null, -1);

        Assert.assertEquals(uniqueReadSetCount, 2);
        ArrayAsserts.assertArrayEquals(duplicateReadCount, new int[]{ 5, 5 });
    }
}