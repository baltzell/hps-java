package org.hps.recon.tracking.lit;

/**
 *
 * @author Norman A Graf
 *
 * @version $Id:
 */
public class CbmLitMaterialInfo implements Comparable {

    //TODO might want to abstract out the specific (length, position) from the general
    double fLength; // Length of the material [cm]
    double fRL; // Radiation length [cm]
    double fRho; // Density [g/cm^3]
    double fZ; // Atomic number
    double fA; // Atomic mass
    double fZpos; // Z position of the material
    String fName; // Name of material

    // default constructor
    public CbmLitMaterialInfo() {

    }

    // fully qualified constructor
    public CbmLitMaterialInfo(String Name, double Zpos, double Length, double Rho, double RL, double Z, double A) {
        fName = Name;
        fZpos = Zpos;
        fA = A;
        fZ = Z;
        fRho = Rho;
        fRL = RL;
        fLength = Length;
    }

    public static CbmLitMaterialInfo getSilicon() {
        return new CbmLitMaterialInfo("Silicon", 0., 0., 2.329, 9.370, 14, 28.0855);
    }


    /*
     * @return Length of the material
     */
    double GetLength() {
        return fLength;
    }

    /*
     * @return Radiation length
     */
    double GetRL() {
        return fRL;
    }

    /*
     * @return Density
     */
    double GetRho() {
        return fRho;
    }

    /*
     * @return Atomic number
     */
    double GetZ() {
        return fZ;
    }

    /*
     * @return Atomic mass
     */
    double GetA() {
        return fA;
    }

    /*
     * @return Z position of the material
     */
    double GetZpos() {
        return fZpos;
    }

    String GetName() {
        return fName;
    }

    /*
     * Sets length of the material
     */
    void SetLength(double length) {
        fLength = length;
    }

    /*
     * Sets radiation length of the material
     */
    void SetRL(double rl) {
        fRL = rl;
    }

    /*
     * Sets density
     */
    void SetRho(double rho) {
        fRho = rho;
    }

    /*
     * Sets atomic number
     */
    void SetZ(double Z) {
        fZ = Z;
    }

    /*
     * Sets atomic mass
     */
    void SetA(double A) {
        fA = A;
    }

    /*
     * Sets Z position of the material
     */
    void SetZpos(double zpos) {
        fZpos = zpos;
    }

    void SetName(String name) {
        fName = name;
    }

    /*
     * @return String representation of the class
     */
    public String toString() {
        StringBuffer ss = new StringBuffer();
        ss.append("MaterialInfo: " + fName + " length=" + fLength + " cm, X0=" + fRL
                + " cm, rho=" + fRho + " g/cm^3 Z=" + fZ + " A=" + fA + " g/mole zpos=" + fZpos
                + " cm \n");
        return ss.toString();
    }

    public int compareTo(Object o) {
        if (o == this) {
            return 0;
        }
        CbmLitMaterialInfo that = (CbmLitMaterialInfo) o;
        if (this.GetZpos() < that.GetZpos()) {
            return -1;
        }
        if (this.GetZpos() > that.GetZpos()) {
            return 1;
        }
        return 0;
    }

}
