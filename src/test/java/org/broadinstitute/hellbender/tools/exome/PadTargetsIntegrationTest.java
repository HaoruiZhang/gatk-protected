package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.codec.digest.DigestUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.List;


public class PadTargetsIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_FILE_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome");
    private static final File TEST_FILE = new File(TEST_FILE_DIR, "targets.tsv");

    @Test
    public void testZeroPadding() {
        final File outputFile = createTempFile("test", ".tsv");

        final String[] arguments = {
                "-" + PadTargets.PADDING_SHORT_NAME, "0",
                "-" + PadTargets.TARGET_FILE_SHORT_NAME, TEST_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
        };
        runCommandLine(arguments);

        String testMd5 = "NEVER POPULATED: TEST";
        String gtMd5 = "NEVER POPULATED: GT";
        try (final FileInputStream fp = new FileInputStream(outputFile)) {
            testMd5 = DigestUtils.md5Hex(fp);

        } catch (final IOException ioe) {
            Assert.fail("Could not open output file", ioe);
        }
        try (final FileInputStream fp = new FileInputStream(TEST_FILE)) {
            gtMd5 = DigestUtils.md5Hex(fp);

        } catch (final IOException ioe) {
            Assert.fail("Could not open test file", ioe);
        }

        Assert.assertEquals(testMd5, gtMd5, "Files should have been exactly the same, but weren't.");

        final List<Target> tc = TargetTableReader.readTargetFile(outputFile);
        Assert.assertNotNull(tc);
    }

    @Test
    public void testSimplePadding() {
        final File outputFile = createTempFile("test", ".tsv");
        final int testPadding = 250;
        final String[] arguments = {
                "-" + PadTargets.PADDING_SHORT_NAME, String.valueOf(testPadding),
                "-" + PadTargets.TARGET_FILE_SHORT_NAME, TEST_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
        };
        runCommandLine(arguments);

        final List<Target> tc = TargetTableReader.readTargetFile(outputFile);
        Assert.assertNotNull(tc);

        Assert.assertEquals(tc.get(0).getStart(), 27003849 - testPadding);
        Assert.assertEquals(tc.get(0).getEnd(), 27003987 + testPadding);
        Assert.assertEquals(tc.get(0).getName(), "target_179698_CRYBB1");

        // index: 93 22	29754711	29754962	target_179839_AP1B1
        Assert.assertEquals(tc.get(93).getStart(), 29754711 - testPadding);
        Assert.assertEquals(tc.get(93).getEnd(), 29754962 + testPadding);
        Assert.assertEquals(tc.get(93).getName(), "target_179839_AP1B1");
    }

    @Test
    public void testBigPadding() {
        final File outputFile = createTempFile("test", ".tsv");
        final int testPadding = 25000;
        final String[] arguments = {
                "-" + PadTargets.PADDING_SHORT_NAME, String.valueOf(testPadding),
                "-" + PadTargets.TARGET_FILE_SHORT_NAME, TEST_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
        };
        runCommandLine(arguments);

        final List<Target> tc = TargetTableReader.readTargetFile(outputFile);
        Assert.assertNotNull(tc);

        Assert.assertEquals(tc.get(0).getStart(), 27003849 - testPadding);
        Assert.assertEquals(tc.get(0).getEnd(), ((27003987 + 27008032)/2));
        Assert.assertEquals(tc.get(0).getName(), "target_179698_CRYBB1");

        // index: 93 (line 95) 22	29754711	29754962	target_179839_AP1B1
        Assert.assertEquals(tc.get(93).getStart(), (29752607 + 29754711)/2 + 1);
        Assert.assertEquals(tc.get(93).getEnd(), ((29754962 + 29755809)/2));
        Assert.assertEquals(tc.get(93).getName(), "target_179839_AP1B1");
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testInvalidPadding() {
        final File outputFile = createTempFile("test", ".tsv");
        final int testPadding = -25000;
        final String[] arguments = {
                "-" + PadTargets.PADDING_SHORT_NAME, String.valueOf(testPadding),
                "-" + PadTargets.TARGET_FILE_SHORT_NAME, TEST_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
        };
        runCommandLine(arguments);
    }
}
