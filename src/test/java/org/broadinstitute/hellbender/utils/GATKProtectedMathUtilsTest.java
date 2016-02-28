package org.broadinstitute.hellbender.utils;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 * Created by davidben on 1/22/16.
 */
public class GATKProtectedMathUtilsTest {
    @Test
    public void naturalLogSumNaturalLogTest() {
        final double[] realSpaceValues = {0.1, 0.2, 1.4, 5.9};
        final double[] logSpaceValues= Arrays.stream(realSpaceValues).map(Math::log).toArray();

        final double expected = Math.log(Arrays.stream(realSpaceValues).sum());
        final double actual = GATKProtectedMathUtils.naturalLogSumExp(logSpaceValues);
        Assert.assertEquals(actual, expected, 1e-10);
    }

    @Test
    public void interquartileRangeTest() {
        Assert.assertEquals(GATKProtectedMathUtils.interquartileRange(new double[] {1,2,3,4,5}), 3, 0.001);
        Assert.assertEquals(GATKProtectedMathUtils.interquartileRange(new double[] {4,2,1,5,3}), 3, 0.001);
        Assert.assertEquals(GATKProtectedMathUtils.interquartileRange(new double[] {1,2,3,4,5,6,7,8,9}), 5, 0.001);
    }

}