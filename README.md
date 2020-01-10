# FortranNamelist
Reading class for Fortran namelists

## Access
Maven central:

```
<dependency>
	<groupId>de.labathome</groupId>
	<artifactId>FortranNamelist</artifactId>
	<version>1.0.0</version>
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
	String inputFile = "";
	try {
		inputFile = new String(Files.readAllBytes(Paths.get(getClass().getResource(testVmecInput).toURI())));
	} catch (Exception e) {
		e.printStackTrace();
	}
	```
	
4. Instantiate the namelist parser and parse the namelist into the previously-defined class structure.

	```java
	FortranNamelist parser = new FortranNamelist(namelist, "mgrid_nli", testClass);
	testClass = (MgridParameters)parser.getParsed();
	```

	




