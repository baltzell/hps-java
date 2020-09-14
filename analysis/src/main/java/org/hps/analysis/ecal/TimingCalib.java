package org.hps.analysis.ecal;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Scanner;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.lcsim.event.Cluster;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.EventHeader;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

import org.hps.evio.RfHit;

/**
 * 2019 ECAL Timing calibration
 * 
 * @author baltzell
 */
public class TimingCalib extends Driver {
   
    AIDA aida = AIDA.defaultInstance();
  
    String inputTimeShiftFile = null;
   
    // prior to 2019, this was 2.004;
    final static double RFPERIOD = 4.008; // ns
   
    final static int NX = 46;
    final static int NY = 10;
    final static int NCHAN = 442;
    final static int XHOLE[] = {-10,-2};
    final static int IXHOLE[] = {13,21};
    
    final static String RFHITS_NAME = "RFHits";
    final static String RAWHITS_NAME = "EcalUncalHits";
    final static String HITS_NAME = "EcalCalHits";
    final static String CLUSTERS_NAME = "EcalClusters";

    final static Mapper MAPPER = new Mapper();
    
    private final TimeShifts timeShifts = new TimeShifts();
   
    public void setTimeShiftFile(String filename) {
        this.inputTimeShiftFile=filename;
    }
    
    // class just to store time shifts and read/write files
    // formatted for conditions database, where the map's key
    // is the database's "ecal_channel_id":
    public static final class TimeShifts extends TreeMap<Integer,Double> {
        private String filename = null;
        public TimeShifts(){
            super();
            for (int ii=0; ii<NCHAN; ii++) this.put(ii,0.0);
        }
        public TimeShifts(String filename) throws IOException {
            super();
            for (int ii=0; ii<NCHAN; ii++) this.put(ii,0.0);
            this.readFile(filename);
        }
        public void readFile(String filename) throws IOException {
            try (FileReader freader = new FileReader(filename)) {
                try (BufferedReader breader = new BufferedReader(freader)) {
                    String line;
                    while ((line=breader.readLine()) != null) {
                        Scanner scanner = new Scanner(line);
                        scanner.useDelimiter(" *");
                        final int id = scanner.nextInt();
                        final double timeShift = scanner.nextDouble();
                        this.put(id,timeShift);
                    }
                    breader.close();
                    freader.close();
                }
            }
        }
        @Override
        public String toString() {
            String s="";
            for (Integer id : this.keySet()) {
                s += String.format("%3d %12.4f\n",id,this.get(id));
            }
            return s;
        }
        public void toFile(String filename) throws IOException {
            try (FileWriter fwriter = new FileWriter(filename)) {
                fwriter.write(this.toString());
            }
        }
        public TimeShifts add(TimeShifts t1,TimeShifts t2) {
            TimeShifts tnew=new TimeShifts();
            for (int ii : this.keySet()) {
                tnew.put(ii,t1.get(ii)+t2.get(ii));
            }
            return tnew;
        }
    }
  
    // class just to get ecal_channel_id for conditions database
    // from ix/iy, presumably this is already available somewhere ...
    public static final class Mapper extends HashMap<Long,Integer> {
        public Mapper() {
            super();
            for (int iy=NY/2; iy>=-NY/2; iy--) {
                for (int ix=-NX/2; ix<=NX/2; ix++) {
                    if (Math.abs(iy)==1 && ix>=XHOLE[0] && ix<=XHOLE[1]) {
                        continue;
                    }
                    this.put(Mapper.hash(ix,iy),Mapper.ixy2id(ix,iy));
                }
            }
        }
        public final int getID(int ix,int iy) {
            return this.get(Mapper.hash(ix,iy));
        }
        public final int getID(CalorimeterHit hit) {
            return this.getID(hit.getIdentifierFieldValue("ix"),
                              hit.getIdentifierFieldValue("iy"));
        }
        private static Long hash(int ix,int iy) {
            return 1000L*(iy+NY/2) + (ix+NX/2);
        }
        private static int ixy2id(int ix,int iy) {
            if (Math.abs(iy)==1 && ix>=XHOLE[0] && ix<=XHOLE[1]) return -1;
            // stitch fake gaps:
            if (ix>0) ix -= 1;
            if (iy>0) iy -= 1;
            // indices from zero, top to bottom, left to right:
            ix += NX/2;
            iy += NY/2;
            iy = (NY-1) - iy;
            // lazy, just loop, slow but correct:
            int id = 0;
            for (int iiy=0; iiy<=iy; iiy++) {
                for (int iix=0; iix<NX; iix++) {
                    if (iy==iiy && iix>=ix) break;
                    if (iiy==NY/2-1 || iiy==NY/2) {
                        if (iix >= IXHOLE[0] && iix <= IXHOLE[1]) continue;
                    }
                    id++;
                }
            }
            return id+1;
        }
    }
   
    private double timeWalk(double energy) {
        final double e = energy>2.3 ? 2.3 : energy;
        final double p[]={1.2,-8.902,1.08,-0.8422,0.1731};
        final double a = p[0]+p[1]*e;
        final double b = p[2]+p[3]*e+p[4]*e*e;
        return Math.exp(a)+b;
    }
    
    @Override
    protected void detectorChanged(Detector detector) {

        // load time offsets from file:
        // (same file format that can be used for conditions database)
        if (inputTimeShiftFile != null) {
            try {
                timeShifts.readFile(inputTimeShiftFile);
            } catch (IOException ex) {
                Logger.getLogger(TimingCalib.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        // intialize histograms:
        aida.tree().cd("/");
        for (int ii=0; ii<NCHAN; ii++) {
            aida.histogram1D(String.format("t_%03d",ii),400,0,200);
            aida.histogram1D(String.format("dt_%03d",ii),100,-20,20);
            aida.histogram1D(String.format("rfdt_%03d",ii),100,-2,2);
        }
        aida.histogram2D("dt_e",200,0,3,200,0,2);
        aida.histogram2D("dtc_e",200,0,3,200,0,2);
        aida.histogram1D("rfdt",1000,-100,100);
        aida.histogram1D("rf",1000,0,400);
        aida.histogram1D("cdt",200,-20,20);
    }
    
    public void printCollections(EventHeader event) {
        List<List<Object>> lol = event.get(Object.class);
        for (List<Object> list : lol) {
            String name = event.getMetaData(list).getName();
            Class type = event.getMetaData(list).getType();
            System.out.println(name+" : "+type);
        }
    } 

    @Override
    public void process(EventHeader event) {

        // get the RF time:
        double rf=-1;
        for (RfHit rfhit : event.get(RfHit.class,RFHITS_NAME)) {
            for (int ii=0; ii<rfhit.getNDouble(); ii++) {
                if (rfhit.getDoubleVal(ii) > 0) {
                    rf=rfhit.getDoubleVal(ii);
                    break;
                }
            }
        }
        if (rf<0) {
            Logger.getLogger(TimingCalib.class.getName()).log(Level.SEVERE, null, "No RF");
            return;
        }
        
        aida.histogram1D("rf").fill(rf);

        List <Cluster> clusters = event.get(Cluster.class,CLUSTERS_NAME);

        // cluster-cluster:
        for (int ii=0; ii<clusters.size(); ii++) {
            Cluster c1 = event.get(Cluster.class,CLUSTERS_NAME).get(ii);
            for (int jj=ii+1; jj<clusters.size(); jj++) {
                Cluster c2 = event.get(Cluster.class,CLUSTERS_NAME).get(jj);
                // top+bottom oly:
                if (c1.getPosition()[1]*c2.getPosition()[1]>=0) continue;
                aida.histogram1D("cdt").fill(c1.getCalorimeterHits().get(0).getTime()
                                            -c2.getCalorimeterHits().get(0).getTime());
            }
        }

        // hits within a cluster:
        for (Cluster cluster : event.get(Cluster.class,CLUSTERS_NAME)) {

            List<CalorimeterHit> hits = cluster.getCalorimeterHits();

            // hits are ordered by energy, the first one is the seed:
            CalorimeterHit seed = hits.get(0);
            final int idSeed = MAPPER.getID(seed);

            // ignore clusters with small seed energies::
            if (seed.getCorrectedEnergy() < 1.0) continue;

            // ignore out-of-time clusters:
            if (seed.getTime()<20 || seed.getTime()>80) continue;
            
            final double rfdt = (seed.getTime()-rf+1e5+RFPERIOD/2)%RFPERIOD-RFPERIOD/2;

            aida.histogram1D(String.format("rfdt_%03d",idSeed)).fill(rfdt);
            aida.histogram1D("rfdt").fill(seed.getTime()-rf);
            
            for (int ihit=1; ihit<hits.size(); ihit++) {

                CalorimeterHit hit = hits.get(ihit);
                final int id = MAPPER.getID(hit);

                if (hit.getCorrectedEnergy() < 0.2) continue;

                aida.histogram1D(String.format("t_%03d",id)).fill(
                        hit.getTime());
                aida.histogram1D(String.format("dt_%03d",id)).fill(
                        hit.getTime()-seed.getTime());
                aida.histogram2D("dt_e").fill(
                        hit.getCorrectedEnergy(),
                        hit.getTime()-seed.getTime());
                aida.histogram2D("dtc_e").fill(
                        hit.getCorrectedEnergy(),
                        hit.getTime()-this.timeWalk(hit.getCorrectedEnergy())
                        -(seed.getTime()-this.timeWalk(seed.getCorrectedEnergy())));

            }
           
        }

    }
}
