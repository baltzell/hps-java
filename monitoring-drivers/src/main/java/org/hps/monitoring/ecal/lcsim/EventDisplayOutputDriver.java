package org.hps.monitoring.ecal.lcsim;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.EventHeader;
import org.hps.recon.ecal.HPSEcalCluster;
import org.hps.recon.ecal.HPSCalorimeterHit;
import org.lcsim.util.Driver;

/**
 * <code>EventDisplayOutputDriver</code> writes the results from clustering
 * and hit reconstruction into a text format that can be read offline by the
 * event display.
 *
 * @author Kyle McCarty
 */
public class EventDisplayOutputDriver extends Driver {
    private FileWriter writer;
    private int eventNum = 0;
    String ecalCollectionName = "EcalHits";
    String clusterCollectionName = "EcalClusters";
    String outputFileName = "cluster-hit.txt";
    
    public void setEcalCollectionName(String ecalCollectionName) {
        this.ecalCollectionName = ecalCollectionName;
    }
    
    public void setClusterCollectionName(String clusterCollectionName) {
        this.clusterCollectionName = clusterCollectionName;
    }
    
    public void setOutputFileName(String outputFileName) {
        this.outputFileName = outputFileName;
    }
    
    public void startOfData() {
        try {
            // Initialize the writer.
            writer = new FileWriter(outputFileName);
            
            // Clear the file.
            writer.write("");
        } catch(IOException e) {
            System.err.println("Error initializing output file for event display.");
        }
    }
    
    public void endOfData() {
        try {
            // Close the file writer.
            writer.close();
        } catch (IOException e) {
            System.err.println("Error closing output file for event display.");
        }
    }
    
    public void process(EventHeader event) {
        // Get the list of clusters.
        List<HPSEcalCluster> clusters;

        // If no cluster collection is present, then make an
        // empty list instead to avoid crashes.
        try {
            clusters = event.get(HPSEcalCluster.class, clusterCollectionName);
            if (clusters == null) {
                throw new RuntimeException("Missing cluster collection!");
            }
        }
        catch(IllegalArgumentException e) {
            clusters = new ArrayList<HPSEcalCluster>(0);
        }
        
        // Get the list of calorimeter hits.
        List<CalorimeterHit> hits;
        
        // If no hit collection is present, then make an empty
        // list instead to avoid crahses.
        try {
            hits = event.get(CalorimeterHit.class, ecalCollectionName);
            if (hits == null) {
                throw new RuntimeException("Missing hit collection!");
            }
        }
        catch(IllegalArgumentException e) {
            hits = new ArrayList<CalorimeterHit>(0);
        }
        
        try {
            if(hits.size() != 0) {//if(clusters.size() != 0) {
                // Increment the event number.
                eventNum++;
                
                // Write the event header.
                writer.append("Event\t" + eventNum + "\n");
                
                // Process the calorimeter hits.
                for (CalorimeterHit hit : hits) {
                    // Get the x/y coordinates for the current hit.
                    int ix = hit.getIdentifierFieldValue("ix");
                    int iy = hit.getIdentifierFieldValue("iy");
                    double energy = hit.getRawEnergy();
                    double time = hit.getTime();
                    
                    // Write the hit to the output file.
                    writer.append(String.format("EcalHit\t%d\t%d\t%f\t%f%n", ix, iy, energy, time));
                }
                
                // Process the clusters.
                for (HPSEcalCluster cluster : clusters) {
                    // Get the seed hit for the cluster.
                    HPSCalorimeterHit seedHit = (HPSCalorimeterHit)cluster.getSeedHit();
                    int ix = seedHit.getIdentifierFieldValue("ix");
                    int iy = seedHit.getIdentifierFieldValue("iy");
                    double time = seedHit.getTime();
                    
                    // Get the cluster's total energy.
                    double energy = cluster.getEnergy();
                    
                    // Write the seed hit to start a cluster.
                    writer.append(String.format("Cluster\t%d\t%d\t%f\t%f%n", ix, iy, energy, time));
                    
                    // Write the component hits to the cluster.
                    for (CalorimeterHit hit : cluster.getCalorimeterHits()) {
                        // Get each component hit's x/y coordinates.
                        ix = hit.getIdentifierFieldValue("ix");
                        iy = hit.getIdentifierFieldValue("iy");
                        
                        // Write them as component hits.
                        writer.append(String.format("CompHit\t%d\t%d%n", ix, iy));
                    }
                }
                
                // Append the end of event indicator.
                writer.append("EndEvent\n");
            }
        } catch(IOException e) {
            System.err.println("Error writing to output for event display.");
        }
    }
}
