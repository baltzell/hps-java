package org.hps.monitoring.svt;

import hep.aida.IAnalysisFactory;
import hep.aida.IHistogram2D;
import hep.aida.IPlotter;
import hep.aida.IProfile1D;

import java.util.List;

import org.lcsim.detector.tracker.silicon.SiSensor;
import org.lcsim.event.EventHeader;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.geometry.Detector;
import org.lcsim.hps.recon.tracking.HPSFittedRawTrackerHit;
import org.lcsim.hps.recon.tracking.HPSSVTCalibrationConstants;
import org.lcsim.hps.recon.tracking.SvtUtils;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 * This Driver makes plots of sensor occupancies across a run. It is intended to
 * be used with the monitoring system. It will currently only run on a test run
 * detector, as the number of sensors, and hence plotter regions, is hardcoded
 * to 20.
 *
 * @author Jeremy McCormick <jeremym@slac.stanford.edu>
 * @version $Id: SensorOccupancyPlotsDriver.java,v 1.12 2012/04/13 00:06:55
 * jeremy Exp $
 *
 */
public class SVTMonitoringPlots extends Driver {

    private String inputCollection = "SVTRawTrackerHits";
    private String fitCollection = "SVTFittedRawTrackerHits";
    private AIDA aida = AIDA.defaultInstance();
    //private AIDAFrame plotterFrame;
    private IPlotter plotter, plotter2, plotter3;
//    private IHistogram1D[][] plots = new IHistogram1D[2][10];
    private IProfile1D[][] pedestalShifts = new IProfile1D[2][10];
//    private IHistogram2D[][] pedestalShifts = new IHistogram2D[2][10];
    private IHistogram2D[][] t0s = new IHistogram2D[2][10];
    private IHistogram2D[][] amps = new IHistogram2D[2][10];

    public SVTMonitoringPlots() {
    }

    public void setInputCollection(String inputCollection) {
        this.inputCollection = inputCollection;
    }

    private int computePlotterRegion(int layer, int module) {
        // Compute the sensor's x and y grid coordinates and then translate to region number.
        int iy = (layer - 1) / 2;
        int ix = 0;
        if (module > 0) {
            ix += 2;
        }
        if (layer % 2 == 0) {
            ix += 1;
        }
        int region = ix * 5 + iy;
        //System.out.println(sensor.getName() + "; lyr=" + layer + "; mod=" + module + " -> xy[" + ix + "][" + iy + "] -> reg="+region);
        return region;
    }

    protected void detectorChanged(Detector detector) {
        //plotterFrame = new AIDAFrame();
        //plotterFrame.setTitle("HPS SVT Monitoring");

        // Setup the plotter.
        IAnalysisFactory fac = aida.analysisFactory();

        plotter = fac.createPlotterFactory().create("Pedestal Shifts");
        plotter.setTitle("Pedestal Shifts");
        //plotterFrame.addPlotter(plotter);
        plotter.style().dataStyle().errorBarStyle().setVisible(false);
        plotter.style().statisticsBoxStyle().setVisible(false);
        plotter.style().dataStyle().fillStyle().setColor("red");
        plotter.style().dataStyle().markerStyle().setShape("dot");
        plotter.style().dataStyle().markerStyle().setSize(1);
        plotter.createRegions(4, 5);

        plotter2 = fac.createPlotterFactory().create("Fitted T0");
        plotter2.setTitle("Fitted T0");
        //plotterFrame.addPlotter(plotter2);
        plotter2.style().dataStyle().errorBarStyle().setVisible(false);
        plotter2.style().statisticsBoxStyle().setVisible(false);
        plotter2.style().setParameter("hist2DStyle", "colorMap");
        plotter2.style().dataStyle().fillStyle().setParameter("colorMapScheme", "rainbow");
        plotter2.style().zAxisStyle().setParameter("scale", "log");
        plotter2.style().zAxisStyle().setVisible(false);
        plotter2.createRegions(4, 5);

        plotter3 = fac.createPlotterFactory().create("Fitted Amplitude");
        plotter3.setTitle("Fitted Amplitude");
        //plotterFrame.addPlotter(plotter3);
        plotter3.style().dataStyle().errorBarStyle().setVisible(false);
        plotter3.style().statisticsBoxStyle().setVisible(false);
        plotter3.style().setParameter("hist2DStyle", "colorMap");
        plotter3.style().dataStyle().fillStyle().setParameter("colorMapScheme", "rainbow");
        plotter3.style().zAxisStyle().setParameter("scale", "log");
        plotter3.style().zAxisStyle().setVisible(false);
        plotter3.createRegions(4, 5);

        // Setup the occupancy plots.
        aida.tree().cd("/");
        for (int module = 0; module < 2; module++) {
            for (int layer = 0; layer < 10; layer++) {
                int region = computePlotterRegion(layer + 1, module);
                SiSensor sensor = SvtUtils.getInstance().getSensor(module, layer);
                pedestalShifts[module][layer] = aida.profile1D(sensor.getName() + " Pedestal Shifts", 640, -0.5, 639.5);
//                pedestalShifts[module][layer] = aida.histogram2D(HPSSVTDAQMaps.sensorArray[module][layer].getName() + " Pedestal Shifts", 640, -0.5, 639.5, 100, -500, 500);
                plotter.region(region).plot(pedestalShifts[module][layer]);
                t0s[module][layer] = aida.histogram2D(sensor.getName() + " Fitted T0", 640, -0.5, 639.5, 100, -24, 72);
                plotter2.region(region).plot(t0s[module][layer]);
                amps[module][layer] = aida.histogram2D(sensor.getName() + " Fitted Amplitude", 640, -0.5, 639.5, 100, 0, 2000);
                plotter3.region(region).plot(amps[module][layer]);
            }
        }
        //plotterFrame.pack();
        //plotterFrame.setVisible(true);
    }

    public void process(EventHeader event) {
        if (event.hasCollection(RawTrackerHit.class, inputCollection)) {

            // Get RawTrackerHit collection from event.
            List<RawTrackerHit> rawTrackerHits = event.get(RawTrackerHit.class, inputCollection);

            // Clear histograms.
            for (int module = 0; module < 2; module++) {
                for (int layer = 1; layer < 11; layer++) {
//                    plotter.region(computePlotterRegion(layer, module)).setYLimits(pedestalShifts[module][layer - 1].minBinHeight()-100, pedestalShifts[module][layer - 1].maxBinHeight()+100);
//                    pedestalShifts[module][layer - 1].reset();
                }
            }

            // Increment strip hit count.
            for (RawTrackerHit hit : rawTrackerHits) {
                int layer = hit.getIdentifierFieldValue("layer"); // 1-10; axial layers are odd layers; stereo layers are even
                int module = hit.getIdentifierFieldValue("module"); // 0-1; module number is top or bottom
                double pedestal = HPSSVTCalibrationConstants.getPedestal((SiSensor) hit.getDetectorElement(), hit.getIdentifierFieldValue("strip"));
                pedestalShifts[module][layer - 1].fill(hit.getIdentifierFieldValue("strip"), hit.getADCValues()[0] - pedestal);
            }
        }

        if (event.hasCollection(HPSFittedRawTrackerHit.class, fitCollection)) {
            List<HPSFittedRawTrackerHit> fits = event.get(HPSFittedRawTrackerHit.class, fitCollection);
            for (HPSFittedRawTrackerHit fit : fits) {
                int layer = fit.getRawTrackerHit().getIdentifierFieldValue("layer"); // 1-10; axial layers are odd layers; stereo layers are even
                int module = fit.getRawTrackerHit().getIdentifierFieldValue("module"); // 0-1; module number is top or bottom
                int strip = fit.getRawTrackerHit().getIdentifierFieldValue("strip");
                if (fit.getShapeFitParameters().getChiSq() < 5) {
                    double noise = HPSSVTCalibrationConstants.getNoise((SiSensor) fit.getRawTrackerHit().getDetectorElement(), strip);
                    if (fit.getAmp() > 4 * noise) {
                        t0s[module][layer - 1].fill(strip, fit.getT0());
                    }
                    amps[module][layer - 1].fill(strip, fit.getAmp());
                }
            }
        }
    }

    public void endOfData() {
    	//plotterFrame.dispose();
    }
}