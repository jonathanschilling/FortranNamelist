# FortranNamelist
Reading class for Fortran namelists

## Access
Maven central:

```
<dependency>
	<groupId>de.labathome</groupId>
	<artifactId>FortranNamelist</artifactId>
	<version>1.0.6</version>
</dependency>
```

## Useage
1. Declare a class and decorate its members with the `@namelist_variable` annotation.

	```java
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
	```
	
	The name of the class does not matter yet.
	Every member, no matter if public or not, that is decorated with `@namelist_variable`, will be looked for in the given text input.
	If you want to use a different name between the Java source that you just wrote and the Fortran namelist,
	you can specify the corresponding name in the namelist via `(name="lstell_sym")` for example.
	
	Parsing a namelist only works reliably if you know the dimensions of the target arrays.
	In the Java constructor of the class, initialize all dynamically-allocated members, e.g. `String`s and especially all arrays.

2. Instantiate the class.

	```java
	MgridParameters testClass = new MgridParameters();
	```

3. Put the namelist you want to parse into a `String`.

	```java
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
	```
	
	You can also read the contents of the namelist from a text file, which might be the more usual way.

	```java
	String filename = "/home/jonathan/someNamelist.nml";
	String inputFile = "";
	try {
		inputFile = new String(Files.readAllBytes(Paths.get(filename)));
	} catch (Exception e) {
		e.printStackTrace();
	}
	```
	
4. Instantiate the namelist parser and parse the namelist into the previously-defined class structure.

	```java
	FortranNamelist parser = new FortranNamelist(namelist, "mgrid_nli", testClass);
	testClass = (MgridParameters)parser.getParsed();
	```
	
	Here, you need to specify the name of the namelist to parse as the second parameter to the constructor
	of the `FortranNamelist` class.
	The method `getParsed()` then parses the given textual namelist into the class given as third parameter to
	the constructor and returns the instance with all specified fields filled with the values from the namelist.
	
5. If you are parsing one- or two-dimensional arrays, you might need to specify the starting indices of the arrays.
	In Fortran, array indices usually start at 1 and this is assumed by default in the `FortranNamelist` parser as well.
	However, you can also have array indices starting at 0 or even a negative number.
	
	Here is an example:
	
	```java
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
	```
	
	The `indata` namelist contains a `String`, a one-dimensional array `ac` and a two-dimensional array `rbc`.
	You already learned above how getting the String value into the parsed class, so let's concentrate on the arrays here.
	
	You can think of the parsing process in Fortran like the following.
	Target arrays for namelists to read by Fortran have to be already initialized at the time of parsing.
	Usually, people initialize them to some default size that is large enough to accomodate any forseeable user input.
	In the given example, the default (and therefore maximum) size of the `ac` array is 21.
	The default size of the `rbc` array is 203x102.
	When using the `FortranNamelist`parser, this behaviour is mimiced by initializing the dynamically-allocated
	members in the constructor to the *same sizes* as they would be initialized in the Fortran equivalent.
	
	It is valid input to the Fortran parser to specify the array contents in either of the following ways:
	1. with explicit indices or index ranges, e.g. `rbc(0,0) = 1.0` or `rbc(0,1:4) = 1.0, 2.0, 3.0, 4.0`
	2. column-wise in a flattened format, e.g. `rbc = 1280*0.0, 1.0, 2.0, ... 475*0.0 ...`
	
	Especially the latter way of specifying the array contents from the top left corner to the bottom right corner
	one value after another creates ambigiuities when parsing the array without knowledge about the intended
	array dimensions. The `FortranNamelist` parser was implemented to specifically address this issue.
	
	Another issue (which is easier to account for) is that starting indices can be any integer.
	The `ac` array starts at an index of 0 and this is told to the parser
	by specifying `dim0min=0` in the `namelist_variable` 	annotation.
	The `rbc` array has starting indices `dim0min = -ntord` and `dim1min = 0`.
	
	During parsing, the indices are simply shifted into valid index ranges of Java.
	Therefore, you can think of then `dimXmin` parameters as an index offset that is subtracted
	from the Fortran indices	during parsing.

