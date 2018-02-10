package org.hps.readout.ecal.updated;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.hps.conditions.database.DatabaseConditionsManager;
import org.hps.conditions.ecal.EcalChannelConstants;
import org.hps.conditions.ecal.EcalConditions;
import org.hps.readout.ReadoutDataManager;
import org.hps.readout.ReadoutDriver;
import org.hps.readout.TempOutputWriter;
import org.hps.readout.util.LcsimCollection;
import org.hps.recon.ecal.EcalRawConverter;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.EventHeader;
import org.lcsim.event.RawCalorimeterHit;
import org.lcsim.geometry.Detector;
import org.lcsim.lcio.LCIOConstants;

/**
 * Class <code>EcalReadoutRawConverterDriver</code> serves as an
 * interface to the {@link org.hps.recon.ecal.EcalRawConverter
 * EcalRawConverter} object that is responsible for converting raw
 * hits generated by the driver {@link
 * org.hps.readout.ecal.updated.EcalReadoutDriver EcalReadoutDriver}
 * into proper {@link org.lcsim.event.CalorimeterHit CalorimeterHit}
 * objects to be used in clustering.
 * 
 * @author Kyle McCarty <mccarty@jlab.org>
 */
public class EcalReadoutRawConverterDriver extends ReadoutDriver {
    
    // ==============================================================
    // ==== LCIO Collections ========================================
    // ==============================================================
    
    /**
     * Sets the name of the input {@link
     * org.lcsim.event.RawCalorimeterHit RawCalorimeterHit}
     * collection.
     */
    private String inputCollectionName = "EcalRawHits";
    /**
     * Sets the name of the output {@link
     * org.lcsim.event.CalorimeterHit CalorimeterHit} collection.
     */
    private String outputCollectionName = "EcalCorrectedHits";
    
    // ==============================================================
    // ==== Driver Options ==========================================
    // ==============================================================
    
    /**
     * Indicates whether channels that are marked as "bad" in the
     * conditions database should be skipped when producing hits.
     */
    private boolean skipBadChannels = true;
    
    // ==============================================================
    // ==== Driver Parameters =======================================
    // ==============================================================
    
    /**
     * The converter object responsible for processing raw hits into
     * proper {@link org.lcsim.event.CalorimeterHit CalorimeterHit}
     * objects.
     */
    private EcalRawConverter converter = new EcalRawConverter();
    /**
     * Cached copy of the calorimeter conditions. All calorimeter
     * conditions should be called from here, rather than by directly
     * accessing the database manager.
     */
    private EcalConditions ecalConditions = null;
    /**
     * Tracks the current local time in nanoseconds for this driver.
     */
    private double localTime = 2.0;
    
    // ==============================================================
    // ==== Debug Output Writers ====================================
    // ==============================================================
    
    /**
     * Outputs debug comparison data for both input hits and output
     * converted hits to a text file.
     */
    private final TempOutputWriter writer = new TempOutputWriter("converted_hits_new.log");
    
    @Override
    public void detectorChanged(Detector detector) {
        // Reset the converter calorimeter conditions.
        converter.setDetector(detector);
        
        // Cache the calorimeter conditions object.
        ecalConditions = DatabaseConditionsManager.getInstance().getEcalConditions();
    }
    
    @Override
    public void process(EventHeader event) {
        writer.write("> Event " + event.getEventNumber() + " - " + ReadoutDataManager.getCurrentTime() + " (Current) - "
                + (ReadoutDataManager.getCurrentTime() - ReadoutDataManager.getTotalTimeDisplacement(outputCollectionName)) + " (Local)");
        writer.write("Input");
        
        // Check the data management driver to determine whether the
        // input collection is available or not.
        if(!ReadoutDataManager.checkCollectionStatus(inputCollectionName, localTime)) {
            return;
        }
        
        // Get all of the raw hits in the current clock-cycle.
        Collection<RawCalorimeterHit> rawHits = ReadoutDataManager.getData(localTime, localTime + 4.0, inputCollectionName, RawCalorimeterHit.class);
        
        // DEBUG :: Write the raw hits seen.
        for(RawCalorimeterHit hit : rawHits) {
            writer.write(String.format("%d;%d;%d", hit.getAmplitude(), hit.getTimeStamp(), hit.getCellID()));
        }
        
        // Increment the local time.
        localTime += 4.0;
        
        // Pass the raw hits to the raw converter to obtain proper
        // calorimeter hits. In readout, raw hits are always Mode-3,
        // so there is no need to check the form.
        List<CalorimeterHit> newHits = new ArrayList<CalorimeterHit>();
        for(RawCalorimeterHit hit : rawHits) {
            // Convert the raw hit.
            CalorimeterHit newHit = converter.HitDtoA(event, hit, 0.0);
            
            // If the hit is on a bad channel, and these are set to
            // be skipped, ignore the hit. Otherwise, add it to the
            // output list.
            if(!(skipBadChannels && isBadChannel(newHit))) {
                newHits.add(newHit);
            }
        }
        
        // Add the calorimeter hit collection to the data manager.
        ReadoutDataManager.addData(outputCollectionName, newHits, CalorimeterHit.class);
        
        // DEBUG :: Write the converted hits seen.
        writer.write("Output");
        for(CalorimeterHit hit : newHits) {
            writer.write(String.format("%f;%f;%d", hit.getRawEnergy(), hit.getTime(), hit.getCellID()));
        }
    }
    
    @Override
    public void startOfData() {
        // Set the LCIO flags for the output collection. Flags are
        // set to store the hit time and hit position respectively.
        int flags = 0;
        flags += 1 << LCIOConstants.RCHBIT_TIME;
        flags += 1 << LCIOConstants.RCHBIT_LONG;
        
        // Define the LCSim collection parameters for this driver's
        // output.
        LcsimCollection<CalorimeterHit> hitCollectionParams = new LcsimCollection<CalorimeterHit>(outputCollectionName,
                this, CalorimeterHit.class, getTimeDisplacement());
        hitCollectionParams.setFlags(flags);
        hitCollectionParams.setPersistent(false);
        
        // Set the dependencies for the driver and register its
        // output collections with the data management driver.
        addDependency(inputCollectionName);
        ReadoutDataManager.registerCollection(hitCollectionParams);
        
        // DEBUG :: Pass the writer to the superclass writer list.
        writers.add(writer);
        
        // Run the superclass method.
        super.startOfData();
    }
    
    @Override
    protected double getTimeDisplacement() {
        return 0;
    }

    @Override
    protected double getTimeNeededForLocalOutput() {
        return 0;
    }
    
    /**
     * Gets the channel parameters for a given channel ID.
     * @param cellID - The <code>long</code> ID value that represents
     * the channel. This is typically acquired from the method {@link
     * org.lcsim.event.CalorimeterHit#getCellID() getCellID()} in a
     * {@link org.lcsim.event.CalorimeterHit CalorimeterHit} object.
     * @return Returns the channel parameters for the channel as an
     * {@link org.hps.conditions.ecal.EcalChannelConstants
     * EcalChannelConstants} object.
     */
    private EcalChannelConstants findChannel(long cellID) {
        return ecalConditions.getChannelConstants(ecalConditions.getChannelCollection().findGeometric(cellID));
    }
    
    /**
     * Indicates whether or not the channel on which a hit occurs is
     * a "bad" channel according to the conditions database.
     * @param hit - The hit to check.
     * @return Returns <code>true</code> if the hit channel is
     * flagged as "bad" and <code>false</code> otherwise.
     */
    private boolean isBadChannel(CalorimeterHit hit) {
        return findChannel(hit.getCellID()).isBadChannel();
    }
    
    /**
     * Indicates whether or not data from channels flagged as "bad"
     * in the conditions system should be ignored. <code>true</code>
     * indicates that they should be ignored, and <code>false</code>
     * that they should not.
     * @param apply - <code>true</code> indicates that "bad" channels
     * will be ignored and <code>false</code> that they will not.
     */
    public void setSkipBadChannels(boolean state) {
        skipBadChannels = state;
    }
    
    /**
     * Sets the name of the input collection containing the objects
     * of type {@link org.lcsim.event.RawCalorimeterHit
     * RawCalorimeterHit} that are output by the digitization driver.
     * This is <code>"EcalRawHits"</code> by default.
     * @param collection - The name of the input raw hit collection.
     */
    public void setInputCollectionName(String collection) {
        inputCollectionName = collection;
    }
    
    /**
     * Sets the number of integration samples that should be included
     * in a pulse integral after the threshold-crossing event.
     * @param samples - The number of samples, where a sample is a
     * 4 ns clock-cycle.
     */
    public void setNsa(int samples) {
        converter.setNSA(4 * samples);
    }
    
    /**
     * Sets the number of integration samples that should be included
     * in a pulse integral before the threshold-crossing event.
     * @param samples - The number of samples, where a sample is a
     * 4 ns clock-cycle.
     */
    public void setNsb(int samples) {
        converter.setNSB(4 * samples);
    }
    
    /**
     * Sets the name of the output collection containing the objects
     * of type {@link org.lcsim.event.CalorimeterHit CalorimeterHit}
     * that are output by this driver. This is
     * <code>"EcalCorrectedHits"</code> by default.
     * @param collection - The name of the output hit collection.
     */
    public void setOutputCollectionName(String collection) {
        outputCollectionName = collection;
    }
}