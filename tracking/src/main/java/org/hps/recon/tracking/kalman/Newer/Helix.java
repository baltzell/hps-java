package kalman;

import java.util.Random;

class Helix { // Create a simple helix oriented along the B field axis for testing the Kalman
              // fit
    Vec p; // Helix parameters drho, phi0, K, dz, tanl
    Vec X0; // Pivot point in the B field reference frame
    private double alpha;
    private double B; // Magnetic field magnitude
    private Vec tB; // Magnetic field direction in the global system
    private Vec uB;
    private Vec vB;
    private double rho;
    private double Q;
    private double radLen;
    private HelixPlaneIntersect hpi;
    RotMatrix R; // Rotation from global coordinates to the B field frame
    Vec origin; // Origin of the B-field reference frame in global coordinates
    private FieldMap fM;
    private Random rndm;

    Helix(Vec HelixParams, Vec pivotGlobal, Vec origin, FieldMap fM) {
        this.origin = origin.copy();
        this.fM = fM;
        p = HelixParams.copy();
        Vec Bf = fM.getField(pivotGlobal);
        B = Bf.mag();
        double c = 2.99793e8;
        alpha = 1000.0 * 1.0E9 / (c * B); // Units are Tesla, mm, GeV
        rho = 2.329; // Density of silicon in g/cm^2
        Q = Math.signum(p.v[2]);
        radLen = (21.82 / rho) * 10.0; // Radiation length of silicon in millimeters
        // System.out.format("helix.java: new helix with radius = %10.2f mm.\n",
        // alpha/p[2]);

        // HelixParams.print("helix constructor helix parameters");
        // X0.print("Helix constructor: pivot in B-field frame");
        // Bf.print("Helix constructor: B field");
        tB = Bf.unitVec(B);
        Vec yhat = new Vec(0., 1.0, 0.);
        uB = yhat.cross(tB).unitVec();
        vB = tB.cross(uB);
        R = new RotMatrix(uB, vB, tB);
        //R.print("new helix rotation matrix");
        X0 = R.rotate(pivotGlobal.dif(origin));
        hpi = new HelixPlaneIntersect();
        //origin.print("new helix origin");
        //pivotGlobal.print("new helix pivot global");
        //X0.print("new helix pivot");
        long seedLong = 123445439;
        rndm = new Random();
        rndm.setSeed(seedLong);
    }

    Helix copy() {
        return new Helix(p, R.inverseRotate(X0).sum(origin), origin, fM);
    }

    void print(String s) {
        System.out.format("Helix parameters for %s:", s);
        System.out.format(" drho=%10.5f", p.v[0]);
        System.out.format(" phi0=%10.5f", p.v[1]);
        System.out.format(" K=%10.5f", p.v[2]);
        System.out.format(" dz=%10.5f", p.v[3]);
        System.out.format(" tanlL=%10.5f\n", p.v[4]);
        System.out.format("         Pivot in B-field frame=%10.5f, %10.5f, %10.5f\n", X0.v[0], X0.v[1], X0.v[2]);
        Vec pivotGlobal = R.inverseRotate(X0).sum(origin);
        System.out.format("         Pivot in global frame=%10.5f, %10.5f, %10.5f\n", pivotGlobal.v[0], pivotGlobal.v[1], pivotGlobal.v[2]);
        System.out.format("         Helix radius=%10.5f, with field B=%10.7f\n", alpha / p.v[2], B);
        R.print("from global frame to B-field frame");
        origin.print("origin of the B-field frame in global coordinates");
    }

    double drho() {
        return p.v[0];
    }

    double phi0() {
        return p.v[1];
    }

    double K() {
        return p.v[2];
    }

    double dz() {
        return p.v[3];
    }

    double tanl() {
        return p.v[4];
    }

    Vec atPhi(double phi) { // return the local coordinates on the helix at a particular phi value
        double x = X0.v[0] + (p.v[0] + (alpha / (p.v[2]))) * Math.cos(p.v[1]) - (alpha / (p.v[2])) * Math.cos(p.v[1] + phi);
        double y = X0.v[1] + (p.v[0] + (alpha / (p.v[2]))) * Math.sin(p.v[1]) - (alpha / (p.v[2])) * Math.sin(p.v[1] + phi);
        double z = X0.v[2] + p.v[3] - (alpha / (p.v[2])) * phi * p.v[4];
        return new Vec(x, y, z);
    }

    Vec atPhiGlobal(double phi) { // return the global coordinates on the helix at a particular phi value
        double x = X0.v[0] + (p.v[0] + (alpha / (p.v[2]))) * Math.cos(p.v[1]) - (alpha / (p.v[2])) * Math.cos(p.v[1] + phi);
        double y = X0.v[1] + (p.v[0] + (alpha / (p.v[2]))) * Math.sin(p.v[1]) - (alpha / (p.v[2])) * Math.sin(p.v[1] + phi);
        double z = X0.v[2] + p.v[3] - (alpha / (p.v[2])) * phi * p.v[4];
        return R.inverseRotate(new Vec(x, y, z)).sum(origin);
    }

    double planeIntersect(Plane Pin) { // phi value where the helix intersects the plane P (given in global
                                       // coordinates)
        Plane P = Pin.toLocal(R, origin);
        double phi = hpi.planeIntersect(p, X0, alpha, P);
        // System.out.format("Helix:planeIntersect: phi = %12.10f\n", phi);
        return phi;
    }

    Vec getMom(double phi) { // get the particle momentum vector at a particular phi value
        double px = -Math.sin(p.v[1] + phi) / Math.abs(p.v[2]);
        double py = Math.cos(p.v[1] + phi) / Math.abs(p.v[2]);
        double pz = p.v[4] / Math.abs(p.v[2]);
        return new Vec(px, py, pz);
    }

    Vec getMomGlobal(double phi) { // get the particle momentum vector at a particular phi value
        double px = -Math.sin(p.v[1] + phi) / Math.abs(p.v[2]);
        double py = Math.cos(p.v[1] + phi) / Math.abs(p.v[2]);
        double pz = p.v[4] / Math.abs(p.v[2]);
        return R.inverseRotate(new Vec(px, py, pz));
    }

    double[] gausRan() { // Return two gaussian random numbers

        double x1, x2, w;
        double[] gran = new double[2];
        do {
            x1 = 2.0 * Math.random() - 1.0;
            x2 = 2.0 * Math.random() - 1.0;
            //x1 = 2.0 * rndm.nextDouble() - 1.0;
            //x2 = 2.0 * rndm.nextDouble() - 1.0;
            w = x1 * x1 + x2 * x2;
        } while (w >= 1.0);
        w = Math.sqrt((-2.0 * Math.log(w)) / w);
        gran[0] = x1 * w;
        gran[1] = x2 * w;

        return gran;
    }    

    Helix randomScat(Plane P, Vec r, Vec pmom, double X) { // Produce a new helix scattered randomly in a given plane P
        // X is the thickness of the silicon material in meters
        // r is the intersection point and pmom the momentum at that point
        double phi = this.planeIntersect(P); // Here the plane P is assumed to be given in global coordinates
        // p.print("randomScat: helix parameters before scatter");

        //Vec r = this.atPhiGlobal(phi);
        double tst = r.dif(P.X()).dot(P.T());
        //System.out.format("randomScat: test dot product %12.6e should be zero\n", tst);
        Vec Bf = fM.getField(r);
        double Bnew = Bf.mag();
        Vec tBnew = Bf.unitVec(Bnew);
        // r.print("randomScat: r global");
        //Vec pmom = getMomGlobal(phi);
        //pmom.print("randomScat: p global");
        Vec t = pmom.unitVec();
        //System.out.format("randomScat: original direction in global coordinates=%10.7f, %10.7f, %10.7f\n", t.v[0],t.v[1],t.v[2]);
        Vec zhat = new Vec(0., 0., 1.);
        Vec uhat = t.cross(zhat).unitVec(); // A unit vector u perpendicular to the helix direction
        Vec vhat = t.cross(uhat);
        RotMatrix R = new RotMatrix(uhat, vhat, t);
        //t.print("initial helix direction in Helix.randomScat");
        // R.print("rotation matrix in Helix.randomScat");
        double ct = Math.abs(P.T().dot(t));
        double theta0;
        // Get the scattering angle
        if (X == 0.)
            theta0 = 0.;
        else
            theta0 = Math.sqrt((X / radLen) / ct) * (0.0136 / pmom.mag()) * (1.0 + 0.038 * Math.log((X / radLen) / ct));
        double[] gran = gausRan();
        double thetaX = gran[0]*theta0;
        double thetaY = gran[1]*theta0;
        //double thetaX = rndm.nextGaussian() * theta0;
        //double thetaY = rndm.nextGaussian() * theta0;
        //System.out.format("Helix.randomScat: X=%12.5e, ct=%12.5e, theta0=%12.5e, thetaX=%12.5e, thetaY=%12.5e\n",X,ct,theta0,thetaX,thetaY);
        double tx = Math.sin(thetaX);
        double ty = Math.sin(thetaY);
        Vec tLoc = new Vec(tx, ty, Math.sqrt(1.0 - tx * tx - ty * ty));
        // tLoc.print("tLoc in Helix.randomScat");
        Vec tnew = R.inverseRotate(tLoc);
        //tnew.print("tnew in Helix.randomScat");
        //System.out.format("tnew dot tnew= %14.10f\n", tnew.dot(tnew));
        //System.out.format("t dot tnew= %14.10f\n", t.dot(tnew));
        double check = Math.acos(Math.min(t.dot(tnew),1.));
        //System.out.format("recalculated scattered angle=%10.7f\n", check);
        
        // Rotate the direction into the frame of the new field (evaluated at the new pivot)
        Vec yhat = new Vec(0.,1.,0.);
        Vec uBnew = yhat.cross(tBnew).unitVec();
        Vec vBnew = tBnew.cross(uBnew);
        RotMatrix RB = new RotMatrix(uBnew, vBnew, tBnew);
        //RB.print("randomscat: field rotation matrix");
        tnew = RB.rotate(tnew);

        double E = pmom.mag(); // Everything is assumed electron
        double sp = 0.0; // 0.002; // Estar collision stopping power for electrons in silicon at about a
                         // GeV, in GeV cm2/g
        double dEdx = 0.1 * sp * rho; // in GeV/mm
        double eLoss = dEdx * X / ct;
        // System.out.format("randomScat: energy=%10.7f, energy loss=%10.7f\n", E, eLoss);
        E = E - eLoss;
        double tanl = tnew.v[2] / Math.sqrt(1.0 - tnew.v[2] * tnew.v[2]);
        double pt = E / Math.sqrt(1.0 + tanl * tanl);
        double K = Q / pt;
        double phi0 = Math.atan2(-tnew.v[0], tnew.v[1]);
        Vec H = new Vec(0., phi0, K, 0., tanl); // Pivot point is on the helix, at the plane intersection point, so drho and dz are zero
        // H.print("scattered helix parameters");  
 
        return new Helix(H, r, P.X(), fM); // Create the new helix with new origin and pivot point
    }

    private SquareMatrix makeF(double[] aP, double[] a, double B) { // For testing purposes only
        double alpha = 1. / B;
        double[][] f = new double[5][5];
        f[0][0] = Math.cos(aP[1] - a[1]);
        f[0][1] = (a[0] + alpha / a[2]) * Math.sin(aP[1] - a[1]);
        f[0][2] = (alpha / (a[2] * a[2])) * (1.0 - Math.cos(aP[1] - a[1]));
        f[0][3] = 0.;
        f[0][4] = 0.;
        f[1][0] = -Math.sin(aP[1] - a[1]) / (aP[0] + alpha / a[2]);
        f[1][1] = (a[0] + alpha / a[2]) * Math.cos(aP[1] - a[1]) / (aP[0] + alpha / a[2]);
        f[1][2] = (alpha / (a[2] * a[2])) * Math.sin(aP[1] - a[1]) / (aP[0] + alpha / a[2]);
        f[1][3] = 0.;
        f[1][4] = 0.;
        f[2][0] = 0.;
        f[2][1] = 0.;
        f[2][2] = 1.0;
        f[2][3] = 0.;
        f[2][4] = 0.;
        f[3][0] = (alpha / a[2]) * a[4] * Math.sin(aP[1] - a[1]) / (aP[0] + alpha / a[2]);
        f[3][1] = (alpha / a[2]) * a[4] * (1.0 - (a[0] + alpha / a[2]) * Math.cos(aP[1] - a[1]) / (aP[0] + alpha / a[2]));
        f[3][2] = (alpha / (a[2] * a[2])) * a[4] * (aP[1] - a[1] - (alpha / a[2]) * Math.sin(aP[1] - a[1]) / (aP[0] + alpha / a[2]));
        f[3][3] = 1.0;
        f[3][4] = -(alpha / a[2]) * (aP[1] - a[1]);
        f[4][0] = 0.;
        f[4][1] = 0.;
        f[4][2] = 0.;
        f[4][3] = 0.;
        f[4][4] = 1.0;
        return new SquareMatrix(5, f);
    }

    Vec pivotTransform(Vec pivot) {
        double xC = X0.v[0] + (p.v[0] + alpha / p.v[2]) * Math.cos(p.v[1]); // Center of the helix circle
        double yC = X0.v[1] + (p.v[0] + alpha / p.v[2]) * Math.sin(p.v[1]);
        // System.out.format("pivotTransform center=%10.6f, %10.6f\n", xC, yC);

        // Predicted state vector
        double[] aP = new double[5];
        aP[2] = p.v[2];
        aP[4] = p.v[4];
        if (p.v[2] > 0) {
            aP[1] = Math.atan2(yC - pivot.v[1], xC - pivot.v[0]);
        } else {
            aP[1] = Math.atan2(pivot.v[1] - yC, pivot.v[0] - xC);
        }
        aP[0] = (xC - pivot.v[0]) * Math.cos(aP[1]) + (yC - pivot.v[1]) * Math.sin(aP[1]) - alpha / p.v[2];
        aP[3] = X0.v[2] - pivot.v[2] + p.v[3] - (alpha / p.v[2]) * (aP[1] - p.v[1]) * p.v[4];

        return new Vec(5, aP);
    }
    
    private double[] pivotTransform(double[] pivot, double[] a, double B) { // For testing purposes only
        double alpha = 1 / B;
        double xC = X0.v[0] + (a[0] + alpha / a[2]) * Math.cos(a[1]); // Center of the helix circle
        double yC = X0.v[1] + (a[0] + alpha / a[2]) * Math.sin(a[1]);

        double[] aP = new double[5];
        aP[2] = a[2];
        aP[4] = a[4];
        if (a[2] > 0) {
            aP[1] = Math.atan2(yC - pivot[1], xC - pivot[0]);
        } else {
            aP[1] = Math.atan2(pivot[1] - yC, pivot[0] - xC);
        }
        aP[0] = (xC - pivot[0]) * Math.cos(aP[1]) + (yC - pivot[1]) * Math.sin(aP[1]) - alpha / a[2];
        aP[3] = X0.v[2] - pivot[2] + a[3] - (alpha / a[2]) * (aP[1] - a[1]) * a[4];

        return aP;
    }
}
