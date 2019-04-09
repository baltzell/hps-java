package org.hps.analysis.hodoscope;

import java.util.List;

import org.lcsim.event.EventHeader;
import org.lcsim.event.EventHeader.LCMetaData;
import org.lcsim.event.SimTrackerHit;
import org.lcsim.geometry.Detector;
import org.lcsim.geometry.IDDecoder;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

import hep.aida.ICloud1D;
import hep.aida.ICloud2D;
import hep.physics.vec.Hep3Vector;

public class SimpleHodoscopeAnalysisDriver extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    private String hitCollName = "HodoscopeHits";

    private ICloud1D simHitXPlot = aida.cloud1D("Sim Hit X");
    private ICloud1D simHitYPlot = aida.cloud1D("Sim Hit Y");
    private ICloud1D simHitZPlot = aida.cloud1D("Sim Hit Z");
    private ICloud2D simHitXYPlot = aida.cloud2D("Sim X vs Y");
    private ICloud1D simHitLenPlot = aida.cloud1D("Sim Hit Len");
    private ICloud1D simHitLayerPlot = aida.cloud1D("Sim Hit Layer");
    private ICloud1D simHitIdXPlot = aida.cloud1D("Sim Hit X ID");
    private ICloud1D simHitIdYPlot = aida.cloud1D("Sim Hit Y ID");
    private ICloud1D simHitEnergyPlot = aida.cloud1D("Sim Hit Energy");
    private ICloud1D simHitPathLenPlot = aida.cloud1D("Sim Hit Path Len");
    private ICloud1D simHitTimePlot = aida.cloud1D("Sim Hit Time");    
    
    protected void detectorChanged(Detector det) {
    }
    
    protected void startOfData() {
    }    

    protected void endOfData() {
    }

    protected void process(EventHeader event) {
        List<SimTrackerHit> simHits = event.get(SimTrackerHit.class, hitCollName);
        LCMetaData meta = event.getMetaData(simHits);
        IDDecoder dec = meta.getIDDecoder();
        for (SimTrackerHit simHit : simHits) {
            double[] pos = simHit.getPosition();
            simHitXPlot.fill(pos[0]);
            simHitYPlot.fill(pos[1]);
            simHitZPlot.fill(pos[2]);
            simHitXYPlot.fill(pos[0], pos[1]);
            Hep3Vector posVec = simHit.getPositionVec();
            simHitLenPlot.fill(posVec.magnitude());
            
            dec.setID(simHit.getCellID64());
            int layer = dec.getValue("layer");
            int ix = dec.getValue("ix");
            int iy = dec.getValue("iy");
            simHitLayerPlot.fill(layer);
            simHitIdXPlot.fill(ix);
            simHitIdYPlot.fill(iy);
            
            simHitEnergyPlot.fill(simHit.getdEdx());
            
            simHitPathLenPlot.fill(simHit.getPathLength());
            
            simHitTimePlot.fill(simHit.getTime());
        }
    }
}