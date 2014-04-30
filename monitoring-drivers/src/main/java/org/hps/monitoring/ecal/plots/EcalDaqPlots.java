package org.hps.monitoring.ecal.plots;

import hep.aida.IHistogram1D;
import hep.aida.IPlotter;
import hep.aida.IPlotterFactory;
import hep.aida.IPlotterStyle;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.hps.util.Resettable;
import org.lcsim.conditions.ConditionsManager;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.EventHeader;
import org.lcsim.geometry.Detector;
import org.lcsim.geometry.compact.Subdetector;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;
/*Conditions system imports*/
import org.hps.conditions.DatabaseConditionsManager;
import org.hps.conditions.DefaultTestSetup;
import org.hps.conditions.ecal.EcalChannel;
import org.hps.conditions.ecal.EcalChannelConstants;
import org.lcsim.detector.converter.compact.EcalCrystal;


/**
 * The driver <code>EcalDaqPlots</code> implements the histogram shown to the user 
 * in the fourth tab of the Monitoring Application, when using the Ecal monitoring lcsim file.
 * It contains only a sub-tab, showing the number of hits recorded by the different FADC channels.
 * It is a very preliminary driver to monitor the DAQ status.
 * These plots are updated continuosly.
 * @author Andrea Celentano
 * @TODO: integrate with the new conditions system.
 * 
 *  */

public class EcalDaqPlots extends Driver implements Resettable {

    private String subdetectorName = "Ecal";
    private String inputCollection = "EcalReadoutHits";
    private IPlotter plotter;
    private AIDA aida;
    private Detector detector;
    private List<IHistogram1D> plots;
    
    private List<Integer> slotsT,slotsB,crates;
    
    private EcalChannel.EcalChannelCollection channels;
    //private DatabaseConditionsManager manager;
    private ConditionsManager manager;
    public EcalDaqPlots() {
    	//manager = DatabaseConditionsManager.getInstance();
    	manager = ConditionsManager.defaultInstance();
    }

    public void setSubdetectorName(String subdetectorName) {
        this.subdetectorName = subdetectorName;
    }

    public void setInputCollection(String inputCollection) {
        this.inputCollection = inputCollection;
    }

    public void detectorChanged(Detector detector) {

    	
    	 if (subdetectorName == null) {
             throw new RuntimeException("The subdetectorName parameter was not set.");
         }

         if (inputCollection == null) {
             throw new RuntimeException("The inputCollection parameter was not set.");
         }
    	    	
        this.detector = detector;
        
        Subdetector subdetector = detector.getSubdetector(subdetectorName);
        
        
    	/*Setup the conditions system*/
        
    	 // Get the channel information from the database.                
        channels = manager.getCachedConditions(EcalChannel.EcalChannelCollection.class, "ecal_channels").getCachedData();

        List<EcalCrystal> crystals = detector.getDetectorElement().findDescendants(EcalCrystal.class);
        /*I do not want the ECAL Crates and Slots to be hard-coded. 
         * It is fine to assume that the channels are from 0 to 15:
         * This is determined by JLAB FADC architecture
         * It is also fine to say that there are 14 slots occupied by FADCs in each crate: 14*16=224, number of channel in each ecal sector (a part from the hole)
         *          * 
         */
        
        slotsT=new ArrayList<Integer>();
        slotsB=new ArrayList<Integer>();
        crates=new ArrayList<Integer>();
        
        // Loop over crystals and get the list of slots-crates
        for (EcalCrystal crystal : crystals) {
        
        	//y>0 means TOP, y<0 means BOTTOM      	       	
        	int y = crystal.getY(); 
            int slot = crystal.getSlot();
            int crate = crystal.getCrate();
            
            if (y>0){
            	if (!slotsT.contains(slot)) slotsT.add(slot);
            }
            else {
            	if (!slotsB.contains(slot)) slotsB.add(slot);
            }
            if (!crates.contains(crate)) crates.add(crate);
        }
        
        System.out.println("These DAQ slots found:");
        /*Order the slots in increasing order*/
     
        Collections.sort(slotsB);
        Collections.sort(slotsT);
        
        System.out.println("TOP: ");
       	for (int slot : slotsT){
            System.out.println(slot);
        }
        System.out.println("BOTTOM: ");
       	for (int slot : slotsT){
            System.out.println(slot);
        }
           
        
        aida = AIDA.defaultInstance();
        aida.tree().cd("/");
        plots = new ArrayList<IHistogram1D>();
        
        
    
        for (int j = 0; j < 14; j++) { // TOP slot           
                plots.add(aida.histogram1D("ECAL: Top Crate Slot " + slotsT.get(j), 16, 0, 16));   
        }
        
        for (int j = 0; j < 14; j++) { // BOTTOM slot           
            plots.add(aida.histogram1D("ECAL: Bottom Crate Slot " + slotsB.get(j), 16, 0, 16));   
        }

        IPlotterFactory factory= aida.analysisFactory().createPlotterFactory("ECAL DAQ Plots");
        plotter =factory.create("DAQ Plots");
        IPlotterStyle pstyle = plotter.style();
        pstyle.dataStyle().fillStyle().setColor("orange");
        pstyle.dataStyle().markerStyle().setColor("orange");
        pstyle.dataStyle().errorBarStyle().setVisible(false);
        plotter.createRegions(7, 4);
        
        int id,plot_id;
        for (int i = 0; i < 2; i++) { // crate
            for (int j = 0; j < 14; j++) { // slot               
                id=i*14+j;
            	plot_id = 0;
                if (i==0){
                    if (j%2==0) plot_id=j*2;
                    else plot_id=(j-1)*2+1;
                }
                else if (i==1){
                    if (j%2==0) plot_id=j*2+2;
                    else plot_id=(j-1)*2+3;
                }
                        
                plotter.region(plot_id).plot(plots.get(id));
            }
        }
        plotter.show();
    }
   
    @Override
    public void reset() {
        if (plotter != null) {
            for (IHistogram1D plot : plots) {
                plot.reset();
            }
        }
    }

    public void process(EventHeader event) {
        if (event.hasCollection(CalorimeterHit.class, inputCollection)) {
            List<CalorimeterHit> hits = event.get(CalorimeterHit.class, inputCollection);
            for (CalorimeterHit hit : hits) {
            	
            // Make an ID to find.
            //  	EcalChannel daq = new EcalChannel.DaqId();
            //	daq.crate = 1;
            //	daq.slot = 2;
            // 	daq.channel = 3;

            	// Find the matching channel.                                     	
            /*
            	EcalChannel.GeometryId geomID=new EcalChannel.GeometryId();
            	geomID.x=hit.getIdentifierFieldValue("ix");
            	geomID.y=hit.getIdentifierFieldValue("iy");
            	EcalChannel channel=channels.findChannel(geomID);

            	int crateN=channel.getCrate();
            	int slotN=channel.getSlot();
            	int channelN=channel.getChannel();
              
            	System.out.println("found channel at " + geomID.x + " " + geomID.y + " corresponding to DAQ crate/slot/channel " + crateN + " "+slotN+" "+channelN);
                
            	//Top CRATE
            	if (geomID.y>0){
            		 int index = slotsT.indexOf(slotN);
            		    plots.get(index).fill(channelN);
            	}
            	else if (geomID.y<0){
            		 int index = slotsB.indexOf(slotN);
            		 plots.get(index+14).fill(channelN);
            		}	*/	          
                }
        }
    }         
}
