package org.hps.users.baltzell;

import hep.aida.ref.function.AbstractIFunction;

/*
 * Function for fitting the leading edge of the RF waveform.
 */
public class RfFitFunction extends AbstractIFunction {
	protected double pedestal=0;
	protected double time=0;
	protected double slope=0;
	public RfFitFunction() {
		this("");
	}
	public RfFitFunction(String title) {
		super();
		this.variableNames=new String[]{"time"};
		this.parameterNames=new String[]{"pedestal","time","slope"};
		init(title);
	}
	public double value(double [] v) {
		return  pedestal + (v[0]-time)*slope;
	}
	public void setParameters(double[] pars) throws IllegalArgumentException {
		super.setParameters(pars);
		pedestal=pars[0];
		time=pars[1];
		slope=pars[2];
	}
	public void setParameter(String key,double value) throws IllegalArgumentException{
		super.setParameter(key,value);
		if      (key.equals("pedestal")) pedestal=value;
		else if (key.equals("time"))     time=value;
		else if (key.equals("slope"))    slope=value;
	}
}