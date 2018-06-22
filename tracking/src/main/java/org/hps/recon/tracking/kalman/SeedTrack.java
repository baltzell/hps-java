package org.hps.recon.tracking.kalman;
//package kalman;

import java.util.ArrayList;
import java.util.Comparator;

// Fit a line/parabola approximation of a helix to a set of measurement points.
// The line and parabola are fit simultaneously in order to handle properly the stereo layers.
// The pivot of the helix is assumed to be the origin of the coordinate system, but the
// coordinate system in which the helix is described may be rotated slightly with respect to
// the global coordinates, in order to optimize its alignment with the field.
class SeedTrack {
    boolean success;
    int[] hits; // Save information for the hit used in each layer
    private double drho, phi0, K, dz, tanl; // Helix parameters derived from the line/parabola fits
    private Vec hParm; // Final helix parameters rotated into the field frame
    //private RotMatrix Rot; // Orthogonal transformation from global to helix (B-field) coordinates
    private SquareMatrix C; // Covariance of helix parameters
    private boolean verbose; // Set true to generate lots of debug printout
    private double alpha; // Conversion constant from 1/pt to helix radius
    private double R, xc, yc; // Radius and center of the "helix"
    private Vec sol; // Fitted polynomial coefficients
    private SquareMatrix Csol; // Covariance matrix of the fitted polynomial coefficients
    private double Bavg; // Average B field
    private double yOrigin;
    private static Plane p0; // x,z plane at y=0
    private static double minDistXZ; // Minimum difference in distance to origin for it to be used in sorting

    double getAlpha() {
        return alpha;
    }

    void print(String s) {
        if (success) {
            System.out.format("Seed track %s: B=%10.7f helix= %10.6f, %10.6f, %10.6f, %10.6f, %10.6f\n", s, Bavg, drho, phi0, K, dz, tanl);
            System.out.format("  seed track hits:");
            for (int j = 0; j < hits.length; j++) {
                System.out.format(" %d ", hits[j]);
            }
            System.out.format("\n");
            Vec pInt = planeIntersection(p0);
            System.out.format("  Distance from origin in X,Z at y=0 is %10.5f\n", pInt.mag());
            C.print("  seed track covariance");
        } else {
            System.out.format("Seed track %s fit unsuccessful.\n", s);
        }
    }

    SeedTrack(ArrayList<SiModule> data, // List of Si modules with data
            double yOrigin, // New origin along beam to use for the fit
            ArrayList<int[]> hitList, // Element 0= index of Si module; Element 1= hit number
            boolean verbose // Set true for lots of debug printout
    ) {
        minDistXZ = 0.25;
        p0 = new Plane(new Vec(0., 0., 0.), new Vec(0., 1., 0.));
        this.verbose = verbose;
        this.yOrigin = yOrigin;
        hits = new int[hitList.size()];
        for (int i = 0; i < hitList.size(); i++) {
            hits[i] = hitList.get(i)[1];
        }

        // Fit a straight line in the non-bending plane and a parabola in the bending plane

        double[] xMC = null;
        double[] y = new double[hitList.size()]; // Global y coordinates of measurements (along beam direction)
        double[] zMC = null;
        double[] t = null;
        double[] yMC = null;
        double[] mTrue = null;
        if (verbose) {
            xMC = new double[hitList.size()]; // Global x coordinates of measurements (bending plane)
            zMC = new double[hitList.size()]; // Global z coordinates of measurements (along B field direction)
            t = new double[hitList.size()]; // Stereo angle
            yMC = new double[hitList.size()];
            mTrue = new double[hitList.size()];
        }
        double[] v = new double[hitList.size()]; // Measurement value (i.e. position of the hit strip)
        double[] s = new double[hitList.size()]; // Uncertainty in the measurement (spatial resolution of the SSD)
        double[][] delta = new double[hitList.size()][3];
        double[][] R2 = new double[hitList.size()][3];
        int N = 0;
        int Nbending = 0;
        int Nnonbending = 0;

        // First find the average field
        Vec Bvec = new Vec(0., 0., 0.);
        for (int[] pnt : hitList) {
            SiModule thisSi = data.get(pnt[0]);
            Vec thisB = thisSi.Bfield.getField(thisSi.toGlobal(new Vec(0., 0., 0.))); // Taking the field from the center of the module
            Bvec = Bvec.sum(thisB);
        }
        int Npnt = hitList.size();
        double sF = 1.0 / ((double) Npnt);
        Bvec = Bvec.scale(sF);
        Bavg = Bvec.mag();
        if (verbose) {
            System.out.format("Entering SeedTrack: Npnt=%d, Bavg=%10.5e\n", Npnt, Bavg);
        }
        double c = 2.99793e8; // Speed of light in m/s
        alpha = 1000.0 * 1.0e9 / (c * Bavg); // Convert from pt in GeV to curvature in mm

        N = 0;
        for (int[] itr : hitList) {
            SiModule thisSi = data.get(itr[0]);
            if (!thisSi.isStereo) {
                Nnonbending++;
            } else {
                Nbending++;
            }
            Measurement m = thisSi.hits.get(itr[1]);
            Vec pnt = new Vec(0., m.v, 0.);
            pnt = thisSi.toGlobal(pnt);
            if (verbose) {
                System.out.format("itr=%d %d, Measurement %d = %10.7f, stereo=%10.7f\n", itr[0], itr[1], N, m.v, thisSi.stereo);
                pnt.print("point global");
                t[N] = thisSi.stereo;
            }
            y[N] = pnt.v[1] - yOrigin;
            if (verbose) {
                xMC[N] = m.rGlobal.v[0]; 
                yMC[N] = m.rGlobal.v[1] - yOrigin;
                zMC[N] = m.rGlobal.v[2]; 
                mTrue[N] = m.vTrue;
            }
            v[N] = m.v;
            s[N] = m.sigma;
            for (int j = 0; j < 3; j++) {
                if (j == 1) {
                    delta[N][j] = thisSi.p.X().v[j] - yOrigin;
                } else {
                    delta[N][j] = thisSi.p.X().v[j];
                }
                R2[N][j] = thisSi.Rinv.M[1][j];
            }
            N++;
        }
        if (Nnonbending < 2) {
            System.out.format("SeedTrack: not enough points in the non-bending plane; N=%d\n", Nnonbending);
            success = false;
            return;
        }
        if (Nbending < 3) {
            System.out.format("SeedTrack: not enough points in the bending plane; N=%d\n", Nbending);
            success = false;
            return;
        }
        if (verbose) {
            System.out.format("SeedTrack: data in global coordinates: y, yMC, zMC, xMC, m, mTrue, check, sigma, theta\n");
            for (int i = 0; i < N; i++) {
                double vcheck = R2[i][0]*(xMC[i]-delta[i][0]) + R2[i][2]*(zMC[i]-delta[i][2]) + R2[i][1]*(yMC[i]-delta[i][1]);
                System.out.format("%d  %10.6f  %10.6f  %10.6f  %10.6f   %10.6f   %10.6f   %10.6f   %10.6f   %8.5f\n", i, y[i], yMC[i], zMC[i], xMC[i], v[i], mTrue[i], vcheck, s[i], t[i]);
            }
        }

        // Here we do the 5-parameter linear fit:
        LinearHelixFit fit = new LinearHelixFit(N, y, v, s, delta, R2, verbose);
        if (verbose) {
            fit.print(N, xMC, y, zMC, v, s);
        }

        // Now, derive the helix parameters and covariance from the two fits
        double sgn = -1.0; // Careful with this sign--->selects the correct half of the circle!
        sol = fit.solution();
        Vec coef = new Vec(sol.v[2], sol.v[3], sol.v[4]); // Parabola coefficients

        double[] circleParams = parabolaToCircle(sgn, coef);
        phi0 = circleParams[1];
        K = circleParams[2];
        drho = circleParams[0];
        double dphi0da = (2.0 * coef.v[2] / coef.v[1]) * square(Math.sin(phi0));
        double dphi0db = Math.sin(phi0) * Math.cos(phi0) / coef.v[1] - square(Math.sin(phi0));
        double dphi0dc = (2.0 * coef.v[0] / coef.v[1]) * square(Math.sin(phi0));
        SquareMatrix D = new SquareMatrix(5);
        double temp = xc * Math.tan(phi0) / Math.cos(phi0);

        double slope = sol.v[1];
        double intercept = sol.v[0];
        tanl = slope * Math.cos(phi0);
        dz = intercept + drho * tanl * Math.tan(phi0);
        if (verbose) {
            System.out.format("SeedTrack: dz=%12.5e   tanl=%12.5e\n", dz, tanl);
            System.out.format("     R=%10.6f, xc=%10.6f, yc=%10.6f\n", R, xc, yc);
            System.out.format("     phi0=%10.7f,  K=%10.6f,   drho=%10.6f\n", phi0, K, drho);
        }
        // Some derivatives to transform the covariance from line and parabola
        // coefficients to helix parameters
        D.M[0][2] = 1.0 / Math.cos(phi0) + temp * dphi0da;
        D.M[0][3] = -(coef.v[1] / (2.0 * coef.v[2])) / Math.cos(phi0) + temp * dphi0db;
        D.M[0][4] = -(sgn + (1.0 - 0.5 * coef.v[1] * coef.v[1]) / Math.cos(phi0)) / (2.0 * coef.v[2] * coef.v[2]) + temp * dphi0dc;
        D.M[1][2] = dphi0da;
        D.M[1][3] = dphi0db;
        D.M[1][4] = dphi0dc;
        D.M[2][4] = -2.0 * alpha * sgn;
        D.M[3][0] = 1.0;
        D.M[3][1] = drho * Math.sin(phi0);
        D.M[4][1] = Math.cos(phi0);
        Csol = fit.covariance();
        C = Csol.similarity(D); // Covariance of helix parameters
        if (verbose) {
            D.print("line/parabola to helix derivatives");
        }

        // Note that the non-bending plane is assumed to be y,z (B field in z
        // direction), and the track is assumed to start out more-or-less
        // in the y direction, so that phi0 should be close to zero. phi0 at 90 degrees
        // will give trouble here!

        // Rotate the result into the frame of the B field at the specified origin
        int[] itr = hitList.get(0);
        SiModule firstSi = data.get(itr[0]);
        Vec firstB = firstSi.Bfield.getField(new Vec(0., yOrigin, 0.));
        Vec zhat = firstB.unitVec();
        Vec yhat = new Vec(0., 1., 0.);
        Vec xhat = (yhat.cross(zhat)).unitVec();
        yhat = zhat.cross(xhat);
        RotMatrix Rot = new RotMatrix(xhat, yhat, zhat);

        hParm = rotateHelix(new Vec(drho, phi0, K, dz, tanl), Rot);
        if (verbose) {
            hParm.print("Seedtrack, rotated helix");
        }

        success = true;
    }

    private double square(double x) {
        return x * x;
    }

    Vec helixParams() { // Return the fitted helix parameters
        return hParm;
    }

    SquareMatrix covariance() { // Return covariance matrix of the fitted helix parameters
        return C;
    }

    Vec solution() { // Return the 5 polynomial coefficients
        return sol;
    }

    double B() { // Return the average field
        return Bavg;
    }

    SquareMatrix solutionCovariance() { // Return covariance of the polynomial coefficients
        return Csol;
    }

    Vec solutionErrors() { // Return errors on the polynomial coefficients (for testing)
        return new Vec(Math.sqrt(Csol.M[0][0]), Math.sqrt(Csol.M[1][1]), Math.sqrt(Csol.M[2][2]), Math.sqrt(Csol.M[3][3]), Math.sqrt(Csol.M[4][4]));
    }

    Vec errors() { // Return errors on the helix parameters
        return new Vec(Math.sqrt(C.M[0][0]), Math.sqrt(C.M[1][1]), Math.sqrt(C.M[2][2]), Math.sqrt(C.M[3][3]), Math.sqrt(C.M[4][4]));
    }

    private double[] parabolaToCircle(double sgn, Vec coef) { // Utility to convert from parabola coefficients to circle (i.e. helix)
                                                              // parameters drho, phi0, and K
        R = -sgn / (2.0 * coef.v[2]);
        yc = sgn * R * coef.v[1];
        xc = coef.v[0] - sgn * R * (1.0 - 0.5 * coef.v[1] * coef.v[1]);
        double[] r = new double[3];
        r[1] = Math.atan2(yc, xc);
        if (R < 0.) {
            r[1] += Math.PI;
        }
        r[2] = alpha / R;
        r[0] = xc / Math.cos(r[1]) - R;
        if (verbose) {
            System.out.format("parabolaToCircle:     R=%10.6f, xc=%10.6f, yc=%10.6f, drho=%10.7f, phi0=%10.7f, K=%10.7f\n", R, xc, yc, r[0], r[1], r[2]);
            coef.print("parabola coefficients");
            double phi02 = Math.atan(-coef.v[1] / (2.0 * coef.v[0] * coef.v[2] + (1.0 - coef.v[1] * coef.v[1])));
            System.out.format("phi02 = %10.7f\n", phi02);
        }
        return r;
    }

    // Transformation of a helix from one B-field frame to another, by rotation R
    Vec rotateHelix(Vec a, RotMatrix R) {
        // a = 5 helix parameters
        // R = 3 by 3 rotation matrix
        // The rotation is easily applied to the momentum vector, so first we transform from helix parameters
        // to momentum, apply the rotation, and then transform back to helix parameters.
        Vec p_prime = R.rotate(aTOp(a));
        double Q = Math.signum(a.v[2]);
        Vec a_prime = pTOa(p_prime, Q);
        return new Vec(a.v[0], a_prime.v[0], a_prime.v[1], a.v[3], a_prime.v[2]);
    }

    // Momentum at the start of the given helix (point closest to the pivot)
    Vec aTOp(Vec a) {
        double px = -Math.sin(a.v[1]) / Math.abs(a.v[2]);
        double py = Math.cos(a.v[1]) / Math.abs(a.v[2]);
        double pz = a.v[4] / Math.abs(a.v[2]);
        if (verbose) {
            a.print("helix parameters in SeedTrack.aTOp");
            System.out.format("SeedTrack.aTOp: p=%10.5f %10.5f %10.5f\n", px, py, pz);
        }
        return new Vec(px, py, pz);
    }

    // Transform from momentum at helix starting point back to the helix parameters
    Vec pTOa(Vec p, Double Q) {
        double phi0 = Math.atan2(-p.v[0], p.v[1]);
        double K = Q / Math.sqrt(p.v[0] * p.v[0] + p.v[1] * p.v[1]);
        double tanl = p.v[2] / Math.sqrt(p.v[0] * p.v[0] + p.v[1] * p.v[1]);
        if (verbose) {
            System.out.format("SeedTrack pTOa: Q=%5.1f phi0=%10.7f K=%10.6f tanl=%10.7f\n", Q, phi0, K, tanl);
            p.print("input momentum vector in SeedTrack.pTOa");
        }
        // Note: the following only makes sense when a.v[0] and a.v[3] (drho and dz) are
        // both zero, i.e. pivot is on the helix
        return new Vec(phi0, K, tanl);
    }

    // Comparator function for sorting seed tracks by curvature
    static Comparator<SeedTrack> curvatureComparator = new Comparator<SeedTrack>() {
        public int compare(SeedTrack t1, SeedTrack t2) {
            double K1 = Math.abs(t1.helixParams().v[2]);
            double K2 = Math.abs(t2.helixParams().v[2]);
            if (K1 > K2) {
                return 1;
            } else {
                return -1;
            }
        }
    };

    // Comparator function for sorting seeds by distance from origin in x,z plane
    static Comparator<SeedTrack> dRhoComparator = new Comparator<SeedTrack>() {
        public int compare(SeedTrack t1, SeedTrack t2) {
            Vec pInt1 = t1.planeIntersection(p0);
            Vec pInt2 = t2.planeIntersection(p0);
            double diff = pInt1.mag() - pInt2.mag();
            if (Math.abs(diff) > minDistXZ) {
                if (diff > 0.) {
                    return 1;
                } else {
                    return -1;
                }
            } else { // The distances aren't very different, so use curvature to break the tie
                double K1 = Math.abs(t1.helixParams().v[2]);
                double K2 = Math.abs(t2.helixParams().v[2]);
                if (K1 > K2) {
                    return 1;
                } else {
                    return -1;
                }
            }
        }
    };

    Vec planeIntersection(Plane p) {
        double arg = (K / alpha) * ((drho + (alpha / K)) * Math.sin(phi0) - (p.X().v[1] - yOrigin));
        double phiInt = -phi0 + Math.asin(arg);
        return atPhi(phiInt);
    }

    private Vec atPhi(double phi) { // point on the helix at the angle phi
        double x = (drho + (alpha / K)) * Math.cos(phi0) - (alpha / K) * Math.cos(phi0 + phi);
        double y = yOrigin + (drho + (alpha / K)) * Math.sin(phi0) - (alpha / K) * Math.sin(phi0 + phi);
        double z = dz - (alpha / K) * phi * tanl;
        return new Vec(x, y, z);
    }
}
