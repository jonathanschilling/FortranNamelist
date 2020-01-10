package examples;

import java.nio.file.Files;
import java.nio.file.Paths;

import de.labathome.FortranNamelist;
import de.labathome.namelist_variable;

public class FortranNamelistExamples {

	public static void main(String[] args) {
		
		run_example_1();
		run_example_2();
	}
	
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
	
	public static void run_example_2() {
		
		
		
	}

}
