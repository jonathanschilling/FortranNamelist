package de.labathome;

import de.labathome.namelist_variable;

/**
 * Prototype class as target for parsing. The default is to look for variables
 * named the same as the variable name in the namelist. Optionally, you can
 * supply a name in the annotation to define a name used in the namelist.
 */
public class MGridParameters {

	@namelist_variable
	String mgrid_ext;

	@namelist_variable
	String mgrid_mode;

	@namelist_variable(name = "lstell_sym")
	boolean stellaratorSymmetric;

	@namelist_variable
	double Rmin;

	@namelist_variable
	double Rmax;

	@namelist_variable
	double Zmin;

	@namelist_variable
	double Zmax;

	@namelist_variable(name = "ir")
	int numR;

	@namelist_variable(name = "jz")
	int numZ;

	@namelist_variable(name = "kp")
	int numPhi;

	// Init variables (especially arrays!) to default values in the constructor.
	public MGridParameters() {
		mgrid_ext = "";
		mgrid_mode = "";
	}
}
