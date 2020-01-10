package examples;

import java.util.Arrays;

import de.labathome.FortranNamelist;
import de.labathome.namelist_variable;

public class FortranNamelistExamples {

	public static void main(String[] args) {
		run_example_1();
		run_example_2();
	}

	/**
	 * This basic example parses Strings, integers and doubles.
	 */
	public static void run_example_1() {

		/**
		 * Prototype class as target for parsing.
		 * The default is to look for variables named the same as the variable name in the namelist.
		 * Optionally, you can supply a name in the annotation to define a name used in the namelist.
		 */
		class MgridParameters {
			@namelist_variable
			String mgrid_ext;

			@namelist_variable
			String mgrid_mode;

			@namelist_variable(name="lstell_sym")
			boolean stellaratorSymmetric;

			@namelist_variable
			double Rmin;

			@namelist_variable
			double Rmax;

			@namelist_variable
			double Zmin;

			@namelist_variable
			double Zmax;

			@namelist_variable(name="ir")
			int numR;

			@namelist_variable(name="jz")
			int numZ;

			@namelist_variable(name="kp")
			int numPhi;

			// Init variables (especially arrays!) to default values in the constructor.
			public MgridParameters() {
				mgrid_ext  = "";
				mgrid_mode = "";
			}
		};

		MgridParameters testClass = new MgridParameters();

		String namelist="&MGRID_NLI\r\n" + 
				"   MGRID_EXT = 'w7x_conf17_rev'\r\n" + 
				"   MGRID_MODE = 'R'\r\n" + 
				"   LSTELL_SYM = .TRUE.\r\n" + 
				"   RMIN = 4.30\r\n" + 
				"   RMAX = 6.30\r\n" + 
				"   ZMIN = -1.20\r\n" + 
				"   ZMAX = 1.20\r\n" + 
				"   IR = 211\r\n" + 
				"   JZ = 241\r\n" + 
				"   KP = 36\r\n" + 
				"/";


		FortranNamelist parser = new FortranNamelist(namelist, "mgrid_nli", testClass);
		testClass = (MgridParameters)parser.getParsed();
	}

	/**
	 * This example parses a 2d array.
	 */
	public static void run_example_2() {

		/**
		 * Prototype class as target for parsing.
		 * The default is to look for variables named the same as the variable name in the namelist.
		 * Optionally, you can supply a name in the annotation to define a name used in the namelist.
		 */
		class VmecInputNamelist {
			/** maximum number of poloidal harmonics (in r,z,lam fourier series) */
			final public static int mpold = 101;
			
			/** maximum number of toroidal harmonics */
			final public static int ntord = 101;

			@namelist_variable
			String mgrid_file;

			@namelist_variable(dim0min=0)
			double[] ac;                  

			@namelist_variable(dim0min=-ntord, dim1min=0)
			double[][] rbc;

			// Init variables (especially arrays!) to default values/sizes in the constructor.
			public VmecInputNamelist() {
				mgrid_file = "NONE";
				ac = new double[21];
				rbc = new double[2*ntord+1][mpold+1]; // actual indices are in [-ntord:ntord][0:mpold]
			}
		}

		VmecInputNamelist vmecInput = new VmecInputNamelist();
		
		String namelist="&indata\r\n" + 
				"   mgrid_file = 'w7x_conf17_rev.nc'\n" + 
				"   AC = 1.0, 2.0, 3.0, 4.0, 42.0\n" +
				"   RBC( 0,0) =  5.52    ZBS(0,0) = -0.0\n" + 
				"   RBC( 1,0) =  0.28143 ZBS(1,0) = -0.23796\n" + 
				"   RBC( 0,1) =  0.48552 ZBS(0,1) =  0.62135\n" + 
				"   RBC( 1,1) = -0.23402 ZBS(1,1) =  0.1838\n" + 
				"/";
		
		FortranNamelist parser = new FortranNamelist(namelist, "indata", vmecInput);
		vmecInput = (VmecInputNamelist)parser.getParsed();
	}

}
