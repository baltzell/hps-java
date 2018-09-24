package org.hps.analysis.alignment.straighttrack;

import java.io.File;
import java.net.URL;
import junit.framework.TestCase;
import org.lcsim.util.cache.FileCache;
import org.lcsim.util.loop.LCSimLoop;

/**
 *
 * @author ngraf
 */
public class StraightTrackGeometryTest extends TestCase {

    public void testIt() throws Exception {
        FileCache cache = new FileCache();
        int nEvents = 1;
        LCSimLoop loop = new LCSimLoop();
        loop.add(new StraightTrackGeometryDriver());
        String fileName = "hpsForwardFullEnergyElectrons_z-2338_bottom_0_SLIC-v06-00-00_QGSP_BERT_HPS-PhysicsRun2016-Nominal-v5-0-nofield_nomsc_recon.slcio";
//        String fileName = "mu-_1.056GeV_slic-3.1.5_geant4-v9r6p1_QGSP_BERT_HPS-EngRun2015-Nominal-v1_fieldOff_++_reco.slcio";
        File inputFile = cache.getCachedFile(new URL("http://www.lcsim.org/test/hps-java/fieldoff/" + fileName));
        loop.setLCIORecordSource(inputFile);
        loop.loop(nEvents);

        System.out.println("Loop processed " + loop.getTotalSupplied() + " events.");
        System.out.println("Done!");
    }

}
