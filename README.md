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

2. Instantiate it.