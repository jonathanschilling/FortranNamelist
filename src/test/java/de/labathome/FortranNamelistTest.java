package de.labathome;

import static org.junit.jupiter.api.Assertions.assertNotNull;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

/**
 * Tests for a class to parse Fortran namelists into Java classes. This uses the
 * information given in the class definition to constrain the parser.
 */
public class FortranNamelistTest {

	@Test
	public final void basicTest() {

		/**
		 * Prototype class as target for parsing. The default is to look for variables
		 * named the same as the variable name in the namelist. Optionally, you can
		 * supply a name in the annotation to define a name used in the namelist.
		 */
		class MgridParameters {

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
			public MgridParameters() {
				mgrid_ext = "";
				mgrid_mode = "";
			}
		}
		;

		MgridParameters testClass = new MgridParameters();

		String namelist = "&MGRID_NLI\r\n" + "   MGRID_EXT = 'w7x_conf17_rev'\r\n" + "   MGRID_MODE = 'R'\r\n"
				+ "   LSTELL_SYM = .TRUE.\r\n" + "   RMIN = 4.30\r\n" + "   RMAX = 6.30\r\n" + "   ZMIN = -1.20\r\n"
				+ "   ZMAX = 1.20\r\n" + "   IR = 211\r\n" + "   JZ = 241\r\n" + "   KP = 36\r\n" + "/";

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(namelist, "mgrid_nli", testClass);
		testClass = (MgridParameters) parser.getParsed();
		assertNotNull(testClass);

		// test for values
		Assertions.assertEquals(testClass.mgrid_ext, "w7x_conf17_rev");
		Assertions.assertEquals(testClass.mgrid_mode, "R");
		Assertions.assertEquals(testClass.stellaratorSymmetric, true);
		Assertions.assertEquals(testClass.Rmin, 4.3);
		Assertions.assertEquals(testClass.Rmax, 6.3);
		Assertions.assertEquals(testClass.Zmin, -1.2);
		Assertions.assertEquals(testClass.Zmax, 1.2);
		Assertions.assertEquals(testClass.numR, 211);
		Assertions.assertEquals(testClass.numZ, 241);
		Assertions.assertEquals(testClass.numPhi, 36);
	}

	/**
	 * Prototype class as target for parsing. The default is to look for variables
	 * named the same as the variable name in the namelist. Optionally, you can
	 * supply a name in the annotation to define a name used in the namelist.
	 */
	class VmecInputNamelist {

		/** maximum number of poloidal harmonics (in r,z,lam fourier series) */
		final public static int mpold = 101;

		/** maximum number of toroidal harmonics */
		final public static int ntord = 101;

		/** 4000; // 1001 // some of these is the actual limit... */
		final public static int nsd = 10001;

		final public static int ndatafmax = nsd;

		final public static int mpol_default = 6;
		final public static int ntor_default = 0;
		final public static int ns_default = 31;

		// from LIBSTELL/Sources/Modules/vsvd0.f

		/** number of external current groups */
		final public static int nigroup = 300;

		/** number of external b-field loop sets allowed */
		final public static int nbsetsp = 5;

		/** number of external poloidal flux loops */
		final public static int nfloops = 100;

		/** number of external b-field coils per set */
		final public static int nbcoilsp = 100;

		/** number of mse measurements */
		final public static int nmse = 100;

		/** number of thompson scattering measurements */
		final public static int ntse = 100;

		// These input variables are mentioned in readin.f, but not in the namelist
		// reading routine anymore.
		// Maybe they have been abandoned in the meantime?
		// ??? phifac: factor scaling toroidal flux to match apres or limiter
		// ??? pknots: array of pressure knot values in SQRT(s) space
		// ??? sknots: array of iota knot values in SQRT(s) space
		// ??? mseprof: offset array from NAMELIST MSEPROFIL \
		// ??? so that the total offset on the i-th MSE data point is | different
		// namelist:
		// ??? taken to be | mseprofile --> ignore here
		// ??? = mseangle_offset+mseangle_offsetM*mseprof(i) /
		// ??? nobd: number of connected flux loop measurements
		// ??? nobser: number of individual flux loop positions
		// ??? nbsets: number of B-coil sets defined in mgrid file
		// ??? nbcoils(n): number of bfield coils in each set defined in mgrid file
		// ??? nbcoilsn: total number of bfield coils defined in mgrid file
		// ??? rbcoil(m,n): R position of the m-th coil in the n-th set from mgrid file
		// ??? zbcoil(m,n): Z position of the m-th coil in the n-th set from mgrid file
		// ??? abcoil(m,n): orientation (surface normal wrt R axis; in radians)
		// ??? of the m-th coil in the n-th set from mgrid file.
		// ??? nflxs: number of flux loop measurements used in matching
		// ??? nbfld(n): number of selected EXTERNAL bfield measurements in set n from
		// nml file
		// ??? nbfldn: total number of EXTERNAL bfield measurements used in matching
		// ??? sigma_current: standard deviation (A) in toroidal current
		// ??? sigma_delphid: standard deviation (Wb) for diamagnetic match

		// NOTE: FOR STANDARD DEVIATIONS (sigma''s) < 0, INTERPRET
		// AS PERCENT OF RESPECTIVE MEASUREMENT
		// THE INPUT DATA FILE WILL ACCEPT BOTH POSITIVE AND NEGATIVE
		// SIGMAS, WHICH IT INTERPRETS DIFFERENTLY. FOR SIGMA > 0, IT
		// TAKES SIGMA TO BE THE STANDARD DEVIATION FOR THAT MEASUREMENT
		// AS DESCRIBED ABOVE. FOR SIGMA < 0, SIGMA IS INTERPRETED AS
		// THE FRACTION OF THE MEASURED DATA NEEDED TO COMPUTE THE ABSOLUTE
		// SIGMA, I.E., (-SIGMA * DATA) = ACTUAL SIGMA USED IN CODE.

		// Commends are based on trunk/VMEC2000/Input_Output/readin.f of the ORNL
		// stellinstall repository.
		// Namelist variables and names are from
		// trunk/LIBSTELL/Sources/Modules/vmec_input.f (same repo).

		/** full path for vacuum Green's function data */
		@namelist_variable
		String mgrid_file;

		/** ??? */
		@namelist_variable
		double time_slice;

		/** number of toroidal field periods ( =1 for Tokamak) */
		@namelist_variable
		int nfp;

		/** flux conserving (=0) or prescribed toroidal current (=1) */
		@namelist_variable
		int ncurr;

		/** ??? */
		@namelist_variable
		int nsin;

		/** ??? probably a single limit for number of iterations */
		@namelist_variable
		int niter;

		/** number of timesteps between printouts on screen */
		@namelist_variable
		int nstep;

		/**
		 * number of iteration steps between accurate calculation of vacuum response;
		 * use fast interpolation scheme in between
		 */
		@namelist_variable
		int nvacskip;

		/** time derivative in force iterations */
		@namelist_variable
		double delt;

		/** ??? probably a single convergence criterion */
		@namelist_variable
		double ftol;

		/** value of compressibility index (gamma=0 => pressure prescribed) */
		@namelist_variable
		double gamma;

		/**
		 * mass or pressure (gamma=0) expansion coefficients (series in s) in MKS units
		 * [NWT/M**2]; Interpretation changes with pmass_type
		 */
		@namelist_variable(dim0min = 0)
		double[] am;

		/**
		 * expansion coefficients for iota (power series in s) used when ncurr=0;
		 * Interpretation changes with piota_type
		 */
		@namelist_variable(dim0min = 0)
		double[] ai;

		/**
		 * expansion coefficients for the normalized (pcurr(s=1) = 1) radial derivative
		 * of the flux-averaged toroidal current density (power series in s) used when
		 * ncurr=1; Interpretation changes with pcurr_type
		 */
		@namelist_variable(dim0min = 0)
		double[] ac;

		/**
		 * expansion coefficients for radial coordinate; can be used to refine grid in
		 * center or at LCFS for example
		 */
		@namelist_variable
		double[] aphi;

		/**
		 * Specifies parameterization type of pcurr function: 'power_series' -
		 * I'(s)=Sum[ ac(j) s ** j] - Default; 'gauss_trunc' - I'(s)=ac(0)
		 * (exp(-(s/ac(1)) ** 2) - exp(-(1/ac(1)) ** 2)); others - see function pcurr
		 */
		@namelist_variable
		String pcurr_type;

		/**
		 * Specifies parameterization type of pmass function: 'power_series' - p(s)=Sum[
		 * am(j) s ** j] - Default; 'gauss_trunc' - p(s)=am(0) (exp(-(s/am(1)) ** 2) -
		 * exp(-(1/am(1)) ** 2)); others - see function pmass
		 */
		@namelist_variable
		String pmass_type;

		/**
		 * Specifies parameterization type of piota function: 'power_series' - p(s)=Sum[
		 * am(j) s ** j] - Default; others - see function piota
		 */
		@namelist_variable
		String piota_type;

		/** Auxiliary array for mass profile. Used for splines, s values */
		@namelist_variable
		double[] am_aux_s;

		/** Auxiliary array for mass profile. Used for splines, function values */
		@namelist_variable
		double[] am_aux_f;

		/** Auxiliary array for iota profile. Used for splines, s values */
		@namelist_variable
		double[] ai_aux_s;

		/** Auxiliary array for iota profile. Used for splines, function values */
		@namelist_variable
		double[] ai_aux_f;

		/** Auxiliary array for current profile. Used for splines, s values */
		@namelist_variable
		double[] ac_aux_s;

		/** Auxiliary array for current profile. Used for splines, function values */
		@namelist_variable
		double[] ac_aux_f;

		/** WAC (anisotropic pres) */
		@namelist_variable(dim0min = 0)
		double[] ah;

		/** WAC (anisotropic pres) */
		@namelist_variable(dim0min = 0)
		double[] at;

		/** WAC (anisotropic pres) */
		@namelist_variable
		double bcrit;

		/** boundary coefficients of COS(m*theta-n*zeta) for R [m] */
		@namelist_variable(dim0min = -ntord, dim1min = 0)
		double[][] rbc;

		/** boundary coefficients of SIN(m*theta-n*zeta) for Z [m] */
		@namelist_variable(dim0min = -ntord, dim1min = 0)
		double[][] zbs;

		/** boundary coefficients of SIN(m*theta-n*zeta) for R [m] */
		@namelist_variable(dim0min = -ntord, dim1min = 0)
		double[][] rbs;

		/** boundary coefficients of COS(m*theta-n*zeta) for Z [m] */
		@namelist_variable(dim0min = -ntord, dim1min = 0)
		double[][] zbc;

		/**
		 * pressure value of pedestal (finite value of pressure at and beyond the LCFS)
		 */
		@namelist_variable
		double spres_ped;

		/**
		 * factor used to scale pressure profile (default value = 1) useful so user can
		 * fix profile and change beta without having to change all AM coefficients
		 * separately
		 */
		@namelist_variable
		double pres_scale;

		/** stellarator-symmetric coefficients of magnetic axis position r */
		@namelist_variable(dim0min = 0)
		double[] raxis_cc;

		/** stellarator-symmetric coefficients of magnetic axis position z */
		@namelist_variable(dim0min = 0)
		double[] zaxis_cs;

		/** stellarator-asymmetric coefficients of magnetic axis position r */
		@namelist_variable(dim0min = 0)
		double[] raxis_cs;

		/** stellarator-asymmetric coefficients of magnetic axis position z */
		@namelist_variable(dim0min = 0)
		double[] zaxis_cc;

		/** upper limit for poloidal mode numbers */
		@namelist_variable
		int mpol;

		/** upper limit for toroidal mode number range */
		@namelist_variable
		int ntor;

		/** number of theta grid points (>=2*mpol+6) */
		@namelist_variable
		int ntheta;

		/** number of zeta grid points (=1 IF ntor=0) */
		@namelist_variable
		int nzeta;

		/** ??? */
		@namelist_variable
		int mfilter_fbdy;

		/** ??? */
		@namelist_variable
		int nfilter_fbdy;

		/**
		 * array of number of iterations (used to terminate run) at each multigrid
		 * iteration
		 */
		@namelist_variable
		int[] niter_array;

		/** array of radial mesh sizes to be used in multigrid sequence */
		@namelist_variable
		int[] ns_array;

		/** array of value of residual(s) at which each multigrid iteration ends */
		@namelist_variable
		double[] ftol_array;

		/** weight factor for constraint force (=1 by DEFAULT) */
		@namelist_variable
		double tcon0;

		/**
		 * specifies type of 2D preconditioner to use ('default', diagonal in m,n,
		 * tri-diagonal in s; 'conjugate-gradient', block tri-di, evolve using cg
		 * method; 'gmres', block tri-di, generalized minimal residual method; 'tfqmr',
		 * block tri-di, transpose-free quasi minimum residual
		 */
		@namelist_variable
		String precon_type;

		/**
		 * value of preconditioned force residuals at which block (2d) tri-di solver is
		 * turned on, if requested via type_prec2d
		 */
		@namelist_variable
		double prec2d_threshold;

		/**
		 * value of toroidal current [A]. Used if ncurr = 1 to specify current profile,
		 * or IF in data reconstruction mode.
		 */
		@namelist_variable
		double curtor;

		/**
		 * ???
		 *
		 * @namelist_variable double sigma_current;
		 *
		 *                    /** array of currents in each external current group. Used
		 *                    to multiply Green's function for fields and loops read in
		 *                    from MGRID file. Should use real current units (A).
		 */
		@namelist_variable
		double[] extcur;

		/** ??? */
		@namelist_variable
		int omp_num_threads;

		/** toroidal flux enclosed by plasma at edge (in Wb) */
		@namelist_variable
		double phiedge;

		/** ??? */
		@namelist_variable
		double[] psa;

		/** ??? */
		@namelist_variable
		double[] pfa;

		/** ??? */
		@namelist_variable
		double[] isa;

		/** ??? */
		@namelist_variable
		double[] ifa;

		/**
		 * = 1 (default), match value of PHIEDGE in input file; = 0, USE pressure
		 * profile width to determine PHIEDGE; = 2, USE LIMPOS data (in mgrid file) to
		 * find PHIEDGE; = 3, USE Ip to find PHIEDGE (fixed-boundary only)
		 */
		@namelist_variable
		int imatch_phiedge;

		/** ??? */
		@namelist_variable
		int opt_raxis;

		/** spline tension for iota */
		@namelist_variable
		double tensi;

		/** spline tension for pressure profile */
		@namelist_variable
		double tensp;

		/** uniform EXPerimental offset of MSE data (calibration offset) ... PLUS ... */
		@namelist_variable
		double mseangle_offset;

		/** multiplier on mseprof offset array (calibration offset) */
		@namelist_variable
		double mseangle_offsetm;

		/**
		 * number of Motional Stark effect data points: >0, USE mse data to find iota;
		 * <=0, fixed iota profile ai
		 */
		@namelist_variable
		int imse;

		/**
		 * number of iota spline points (computed internally unless specified
		 * explicitly)
		 */
		@namelist_variable
		int isnodes;

		/** ??? */
		@namelist_variable
		double[] rstark;

		/** pitch angle data from stark measurement */
		@namelist_variable
		double[] datastark;

		/** standard deviation (degrees) in MSE data */
		@namelist_variable
		double[] sigma_stark;

		/**
		 * number of pressure profile data points: = 0, no thompson scattering data to
		 * READ
		 */
		@namelist_variable
		int itse;

		/**
		 * number of pressure spline points (computed internally unless specified
		 * explicitly)
		 */
		@namelist_variable
		int ipnodes;

		/** number by which Thomson scattering data is scaled to get actual pressure */
		@namelist_variable
		double presfac;

		/** uniform arbitrary radial offset of pressure data */
		@namelist_variable
		double pres_offset;

		/** ??? */
		@namelist_variable
		double[] rthom;

		/** pressure data from Thompson, CHEERS (Pa) */
		@namelist_variable
		double[] datathom;

		/** standard deviation (Pa) for pressure profile data */
		@namelist_variable
		double[] sigma_thom;

		/** diamagnetic toroidal flux (Wb) */
		@namelist_variable
		double phidiam;

		/** ??? */
		@namelist_variable
		double sigma_delphid;

		/** vbl spline tension for iota */
		@namelist_variable
		double tensi2;

		/**
		 * vbl spline tension form factor (note: IF tensi!=tensi2 THEN tension(i-th
		 * point) = tensi+(tensi2-tensi)*(i/n-1))**fpolyi
		 */
		@namelist_variable
		double fpolyi;

		/** number of flux loop measurements used in matching */
		@namelist_variable
		int nflxs;

		/** array giving INDEX of flux measurement in iconnect array */
		@namelist_variable
		int[] indxflx;

		/**
		 * measured flux loop signals corresponding to the combination of signals in
		 * iconnect array
		 */
		@namelist_variable
		double[] dsiobt;

		/** standard deviaton (Wb) for EXTERNAL poloidal flux data */
		@namelist_variable
		double[] sigma_flux;

		/** number of selected EXTERNAL bfield measurements in set n from nml file */
		@namelist_variable
		int[] nbfld;

		/** array giving INDEX of bfield measurement used in matching */
		@namelist_variable
		int[][] indxbfld;

		/** scaling factor for phiedge: bloat plasma (bloat>1) */
		@namelist_variable
		double bloat;

		/** Backwards compatibility: Obsolete --> see raxis_cc, raxis_cs */
		@namelist_variable(dim0min = 0)
		double[] raxis;

		/** Backwards compatibility: Obsolete --> see raxis_cc, raxis_cs */
		@namelist_variable(dim0min = 0)
		double[] zaxis;

		/**
		 * measured magnetic field at rbcoil(m,n),zbcoil(m,n) at the orientation
		 * br*COS(abcoil) + bz*SIN(abcoil)
		 */
		@namelist_variable
		double[][] bbc;

		/** standard deviation (T) for EXTERNAL magnetic field data */
		@namelist_variable
		double[][] sigma_b;

		/**
		 * LOGICAL variable. =.true. IF pressure data are prescribed in REAL space.
		 * =.false. IF data in flux space.
		 */
		@namelist_variable
		boolean lpofr;

		/**
		 * =T, use non-variational forces to ensure <EQUIF> = 0; =F, use variational
		 * form of forces, <EQUIF> ~ 0
		 */
		@namelist_variable
		boolean lforbal;

		/** =T, run in free boundary mode if mgrid_file exists */
		@namelist_variable
		boolean lfreeb;

		/** ??? */
		@namelist_variable
		boolean lmove_axis;

		/** ??? reconstruction mode ??? */
		@namelist_variable
		boolean lrecon;

		/** ??? */
		@namelist_variable
		boolean lmac;

		/** =T, run in asymmetric mode; =F, run in stellarator symmetry mode */
		@namelist_variable
		boolean lasym;

		/** ??? */
		@namelist_variable
		boolean ledge_dump;

		/** Obsolete */
		@namelist_variable
		boolean lspectrum_dump;

		/** Obsolete */
		@namelist_variable
		boolean loptim;

		/** ??? reverse-field pinch ??? */
		@namelist_variable
		boolean lrfp;

		/** to obtain old fort.8 file */
		@namelist_variable
		boolean loldout;

		/** J.Geiger: for txt- and diagno-output */
		@namelist_variable
		boolean lwouttxt;

		/** J.Geiger: for txt- and diagno-output */
		@namelist_variable
		boolean ldiagno;

		/** J.Geiger: to force full 3D1-output */
		@namelist_variable
		boolean lfull3d1out;

		/**
		 * maximum number of iterations of the main loop, i.e. of all values in ns_array
		 */
		@namelist_variable
		int max_main_iterations;

		/** inserted M.Drevlak */
		@namelist_variable
		boolean lgiveup;

		/** inserted M.Drevlak, giveup-factor for ftolv */
		@namelist_variable
		double fgiveup;

		/** J Hanson See jxbforce coding */
		@namelist_variable
		boolean lbsubs;

		/** which NESTOR version to use */
		@namelist_variable
		int vac_1_2;

		// Init variables (especially arrays!) to default values/sizes in the
		// constructor.
		// Updated based on trunk/LIBSTELL/Sources/Modules/{vmec_input.f, vparams.f} and
		// trunk/VMEC2000/Sources/Input_Output/readin.f
		// from the stellinstall repository of ORNL at
		// https://github.com/ORNL-Fusion/stellinstall
		public VmecInputNamelist() {

			// actual indices are in [-ntord:ntord][0:mpold]
			rbc = new double[2 * ntord + 1][mpold + 1];
			zbs = new double[2 * ntord + 1][mpold + 1];
			rbs = new double[2 * ntord + 1][mpold + 1];
			zbc = new double[2 * ntord + 1][mpold + 1];

			am = new double[21];
			ai = new double[21];
			ac = new double[21];
			aphi = new double[20];

			am_aux_s = new double[ndatafmax];
			am_aux_f = new double[ndatafmax];
			ai_aux_s = new double[ndatafmax];
			ai_aux_f = new double[ndatafmax];
			ac_aux_s = new double[ndatafmax];
			ac_aux_f = new double[ndatafmax];

			ah = new double[21];
			at = new double[21];

			// actual indices are in [0:ntord]
			raxis_cc = new double[ntord + 1];
			zaxis_cs = new double[ntord + 1];
			raxis_cs = new double[ntord + 1];
			zaxis_cc = new double[ntord + 1];
			raxis = new double[ntord + 1];
			zaxis = new double[ntord + 1];

			ns_array = new int[100];
			ftol_array = new double[100];
			niter_array = new int[100];

			psa = new double[ndatafmax];
			pfa = new double[ndatafmax];
			isa = new double[ndatafmax];
			ifa = new double[ndatafmax];

			extcur = new double[nigroup];

			nbfld = new int[nbsetsp];
			indxflx = new int[nfloops];
			indxbfld = new int[nbcoilsp][nbsetsp];

			rstark = new double[nmse];
			datastark = new double[nmse];
			sigma_stark = new double[nmse];

			rthom = new double[ntse];
			datathom = new double[ntse];
			sigma_thom = new double[ntse];

			dsiobt = new double[nfloops];
			sigma_flux = new double[nfloops];
			bbc = new double[nbcoilsp][nbsetsp];
			sigma_b = new double[nbcoilsp][nbsetsp];

			omp_num_threads = 8;
			gamma = 0;
			spres_ped = 1;
			mpol = mpol_default;
			ntor = ntor_default;
			ntheta = 0;
			nzeta = 0;
			ns_array[0] = ns_default;
			Arrays.fill(niter_array, -1);
			bloat = 1;
			time_slice = 0;
			nfp = 1;
			ncurr = 0;
			nsin = ns_default;
			niter = 100;
			nstep = 10;
			nvacskip = 1;
			delt = 1;
			ftol = 1.e-10;
			ftol_array[0] = ftol;
			aphi[0] = 1;
			pres_scale = 1;
			mfilter_fbdy = -1;
			nfilter_fbdy = -1;
			tcon0 = 1;
			precon_type = "NONE";
			prec2d_threshold = 1.e-30;
			curtor = 0;
			phiedge = 1;
			mgrid_file = "NONE";
			lfreeb = true;
			lmove_axis = true;
			lmac = false;

			// SPH: changed 05-14-14
			lforbal = false;

			lasym = false;
			lrfp = false;

			// obtain old fort.8 file
			loldout = false;

			// J Geiger 2010-05-04: for txt- and diagno-output
			ldiagno = false;

			// inserted M.Drevlak
			lgiveup = false;

			// inserted M.Drevlak
			fgiveup = 3.E+01;

			// J Hanson. See jxbforce coding
			lbsubs = false;

			// J Geiger & SPH (5-21-15): to force full 3D1-output
			lfull3d1out = true;

			// to keep a presumably expected standard behavior.
			max_main_iterations = 1;

			// #if defined(NETCDF)
			// to keep functionality as expected with netcdf
			lwouttxt = false;
			// J.Geiger: for txt- and diagno-output
			// #else
			// lwouttxt = .true. ! and without netcdf
			// #endif
			// J Geiger 2010-05-04 end

			pcurr_type = "power_series";
			piota_type = "power_series";
			pmass_type = "power_series";

			// ANISTROPY PARAMETERS
			bcrit = 1;
			at[0] = 1;

			vac_1_2 = 1;

			//
			// BACKWARDS COMPATIBILITY
			//
			// Work around a bug in gfortran. When performing an optimized build, the WHERE
			// statement would produce incorrect results. Work around this bug by expanding
			// the full WHERE statment. This should have no adverse effects on any other
			// compiler since these statements are equivalent to the older code statement.
			//
			// WHERE (raxis .ne. 0.0_dp) raxis_cc = raxis
			// WHERE (zaxis .ne. 0.0_dp) zaxis_cs = zaxis
			//
			// The effect of this bug optimized to code to effectively ignore the WHERE
			// statement and assign all value values of the r/zaxis to the r/zaxis_cc/s
			// arrays. Explicitly adding the r/zaxis .ne. 0.0_dp section prevents this. This
			// bug is known to exist in gfortran 4.9. It may manifest in other versions.
			// WHERE (raxis .ne. 0.0_dp)
			// raxis_cc = raxis
			// ELSEWHERE
			// raxis_cc = raxis_cc
			// ENDWHERE
			// WHERE (zaxis .ne. 0.0_dp)
			// zaxis_cs = zaxis
			// ELSEWHERE
			// zaxis_cs = zaxis_cs
			// ENDWHERE
			//
			// ==> If something is given in raxis or zaxis, put it in the corresponding
			// positions of raxis_cc and zaxis_cs.
		}
	}

	@Test
	void vmecInputTest1() throws IOException {

		String inputFilename = "/exampleVmecInput/input.v3fit_jsch_20171207_006_2200ms_pm_10ms_v0_Bprecheck";
		InputStream inputFileStream = getClass().getResourceAsStream(inputFilename);

		String inputFile = new String(inputFileStream.readAllBytes(), StandardCharsets.UTF_8);

		VmecInputNamelist vmecInput = new VmecInputNamelist();

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
		vmecInput = (VmecInputNamelist) parser.getParsed();
		assertNotNull(vmecInput);

		// Numerical Resolution and Symmetry Assumptions
		Assertions.assertEquals(false, vmecInput.lasym);
		Assertions.assertEquals(5, vmecInput.nfp);
		Assertions.assertEquals(8, vmecInput.mpol);
		Assertions.assertEquals(8, vmecInput.ntor);
		Assertions.assertEquals(22, vmecInput.ntheta);
		Assertions.assertEquals(36, vmecInput.nzeta);

		// Multi-Grid Steps
		int[] ns_array = new int[100];
		ns_array[0] = 51;
		Assertions.assertArrayEquals(ns_array, vmecInput.ns_array);

		double[] ftol_array = new double[100];
		ftol_array[0] = 1.0e-19;
		Assertions.assertArrayEquals(ftol_array, vmecInput.ftol_array);

		int[] niter_array = new int[100];
		Arrays.fill(niter_array, 100000);
		Assertions.assertArrayEquals(niter_array, vmecInput.niter_array);

		// Global Physics Parameters
		Assertions.assertEquals(1.86086390731511, vmecInput.phiedge);
		Assertions.assertEquals(1, vmecInput.ncurr);

		// Profile of Mass or Pressure
		Assertions.assertEquals("two_power", vmecInput.pmass_type.trim());

		double[] am = new double[21];
		am[0] = 66978.0973204654;
		am[1] = 1.90820658424883;
		am[2] = 4.46347452430851;
		Assertions.assertArrayEquals(am, vmecInput.am);

		Assertions.assertArrayEquals(new double[VmecInputNamelist.nsd], vmecInput.am_aux_s);
		Assertions.assertArrayEquals(new double[VmecInputNamelist.nsd], vmecInput.am_aux_f);

		Assertions.assertEquals(1.0, vmecInput.pres_scale);
		Assertions.assertEquals(0.0, vmecInput.gamma);
		Assertions.assertEquals(1.0, vmecInput.spres_ped);

		// (Initial Guess for) Rotational Transform Profile
		Assertions.assertEquals("power_series", vmecInput.piota_type.trim());

		Assertions.assertArrayEquals(new double[21], vmecInput.ai);

		Assertions.assertArrayEquals(new double[VmecInputNamelist.nsd], vmecInput.ai_aux_s);
		Assertions.assertArrayEquals(new double[VmecInputNamelist.nsd], vmecInput.ai_aux_f);

		// (Initial Guess for) Toroidal Current Profile
		Assertions.assertEquals("sum_atan", vmecInput.pcurr_type.trim());

		double[] ac = new double[21];
		ac[1] = 1.0;
		ac[2] = 4.65367221956161;
		ac[3] = 2.27775992394378;
		ac[4] = 1.0;
		ac[5] = -0.943974780239535;
		ac[6] = 4.97719715801010;
		ac[7] = 1.0;
		ac[8] = 1.0;
		Assertions.assertArrayEquals(ac, vmecInput.ac);

		Assertions.assertArrayEquals(new double[VmecInputNamelist.nsd], vmecInput.ac_aux_s);
		Assertions.assertArrayEquals(new double[VmecInputNamelist.nsd], vmecInput.ac_aux_f);

		Assertions.assertEquals(-913.791680825795, vmecInput.curtor);
		Assertions.assertEquals(1.0, vmecInput.bloat);

		// Free-Boundary Parameters
		Assertions.assertEquals(true, vmecInput.lfreeb);
		Assertions.assertEquals("/u/jsch/v3fit_runs/MGRIDS/mgrid_w7x_v2_10mm_grid.nc", vmecInput.mgrid_file.trim());

		double[] extcur = new double[VmecInputNamelist.nigroup];
		extcur[0] = -12989.0;
		extcur[1] = -12989.0;
		extcur[2] = -12989.0;
		extcur[3] = -12989.0;
		extcur[4] = -12989.0;
		Assertions.assertArrayEquals(extcur, vmecInput.extcur);

		Assertions.assertEquals(10, vmecInput.nvacskip);

		// Tweaking Parameters
		Assertions.assertEquals(200, vmecInput.nstep);

		double[] aphi = new double[20];
		aphi[0] = 1.0;
		Assertions.assertArrayEquals(aphi, vmecInput.aphi);

		Assertions.assertEquals(0.9, vmecInput.delt);
		Assertions.assertEquals(1.0, vmecInput.tcon0);
		Assertions.assertEquals(false, vmecInput.lforbal);

		// Initial Guess for Magnetic Axis Geometry
		double[] raxis_cc = new double[VmecInputNamelist.ntord + 1];
		raxis_cc[0] = 5.60519782116720;
		raxis_cc[1] = 0.359696693960237;
		raxis_cc[2] = 1.335514497262088E-002;
		raxis_cc[3] = 7.259351998922190E-004;
		raxis_cc[4] = -6.817192552674108E-005;
		raxis_cc[5] = -4.866604255356447E-005;
		raxis_cc[6] = -9.461962764212337E-005;
		raxis_cc[7] = -2.828730415733434E-005;
		raxis_cc[8] = -1.554463199957404E-004;
		Assertions.assertArrayEquals(raxis_cc, vmecInput.raxis_cc);

		double[] zaxis_cs = new double[VmecInputNamelist.ntord + 1];
		zaxis_cs[1] = -0.301269608014063;
		zaxis_cs[2] = -1.618262861498088E-002;
		zaxis_cs[3] = 1.319614987094394E-004;
		zaxis_cs[4] = 7.237442086550272E-005;
		zaxis_cs[5] = -7.887293984398275E-005;
		zaxis_cs[6] = -5.489994367001566E-006;
		zaxis_cs[7] = 2.463216529452668E-005;
		zaxis_cs[8] = -1.340415732215003E-004;
		Assertions.assertArrayEquals(zaxis_cs, vmecInput.zaxis_cs);

		double[] raxis_cs = new double[VmecInputNamelist.ntord + 1];
		raxis_cs[1] = -2.093456611578366E-016;
		raxis_cs[2] = -4.361367940788264E-016;
		raxis_cs[3] = -2.267911329209897E-015;
		raxis_cs[4] = -1.500310571631163E-015;
		raxis_cs[5] = -1.395637741052244E-016;
		raxis_cs[6] = 1.290964910473326E-015;
		raxis_cs[7] = -1.116510192841796E-015;
		raxis_cs[8] = 2.233020385683591E-015;
		Assertions.assertArrayEquals(raxis_cs, vmecInput.raxis_cs);

		double[] zaxis_cc = new double[VmecInputNamelist.ntord + 1];
		zaxis_cc[1] = 3.679904200040097E-017;
		zaxis_cc[2] = -4.906538933386797E-018;
		zaxis_cc[3] = 6.814637407481662E-018;
		zaxis_cc[4] = -2.807630611882445E-017;
		zaxis_cc[5] = -6.187690765993350E-017;
		zaxis_cc[6] = 8.450150385277261E-018;
		zaxis_cc[7] = -5.642519773394817E-017;
		zaxis_cc[8] = 4.497660688937897E-018;
		Assertions.assertArrayEquals(zaxis_cc, vmecInput.zaxis_cc);

		// (Initial Guess for) Boundary Geometry

		double[][] rbc = new double[2 * VmecInputNamelist.ntord + 1][VmecInputNamelist.mpold + 1];
		rbc[VmecInputNamelist.ntord + 0][0] = 5.53767206572562;
		rbc[VmecInputNamelist.ntord + 1][0] = 0.274510550203623;
		rbc[VmecInputNamelist.ntord + 2][0] = -6.981826227607880E-003;
		rbc[VmecInputNamelist.ntord + 3][0] = 1.813036409387450E-003;
		rbc[VmecInputNamelist.ntord + 4][0] = -1.816182024820829E-003;
		rbc[VmecInputNamelist.ntord + 5][0] = -9.457063618548155E-007;
		rbc[VmecInputNamelist.ntord + 6][0] = -3.598476797219176E-004;
		rbc[VmecInputNamelist.ntord + 7][0] = 5.201212576617179E-005;
		rbc[VmecInputNamelist.ntord + 8][0] = -3.252552333940336E-004;
		//
		rbc[VmecInputNamelist.ntord - 8][1] = 1.896612909948546E-004;
		rbc[VmecInputNamelist.ntord - 7][1] = -2.456907136966591E-004;
		rbc[VmecInputNamelist.ntord - 6][1] = 2.171005552378400E-004;
		rbc[VmecInputNamelist.ntord - 5][1] = -3.651594436424253E-004;
		rbc[VmecInputNamelist.ntord - 4][1] = -5.338074727927220E-004;
		rbc[VmecInputNamelist.ntord - 3][1] = 2.506984615739400E-004;
		rbc[VmecInputNamelist.ntord - 2][1] = -8.846005854412054E-004;
		rbc[VmecInputNamelist.ntord - 1][1] = 1.678418664958850E-002;
		rbc[VmecInputNamelist.ntord + 0][1] = 0.495662386868802;
		rbc[VmecInputNamelist.ntord + 1][1] = -0.206910682038032;
		rbc[VmecInputNamelist.ntord + 2][1] = -1.680826411787574E-002;
		rbc[VmecInputNamelist.ntord + 3][1] = 7.103573636538324E-004;
		rbc[VmecInputNamelist.ntord + 4][1] = 4.245820612023922E-004;
		rbc[VmecInputNamelist.ntord + 5][1] = 6.969490097295169E-004;
		rbc[VmecInputNamelist.ntord + 6][1] = -7.039468898728505E-004;
		rbc[VmecInputNamelist.ntord + 7][1] = 5.333996925958589E-004;
		rbc[VmecInputNamelist.ntord + 8][1] = 6.689061049649439E-005;
		//
		rbc[VmecInputNamelist.ntord - 8][2] = 1.151620243875582E-004;
		rbc[VmecInputNamelist.ntord - 7][2] = -3.184094551899957E-004;
		rbc[VmecInputNamelist.ntord - 6][2] = 4.300662767725242E-004;
		rbc[VmecInputNamelist.ntord - 5][2] = -3.210600174477836E-004;
		rbc[VmecInputNamelist.ntord - 4][2] = 2.624531379976409E-004;
		rbc[VmecInputNamelist.ntord - 3][2] = 5.677789104917515E-005;
		rbc[VmecInputNamelist.ntord - 2][2] = 2.317403085550396E-003;
		rbc[VmecInputNamelist.ntord - 1][2] = 9.514153491031787E-003;
		rbc[VmecInputNamelist.ntord + 0][2] = 3.302699802057374E-002;
		rbc[VmecInputNamelist.ntord + 1][2] = 4.021410132068875E-002;
		rbc[VmecInputNamelist.ntord + 2][2] = 6.948281839736793E-002;
		rbc[VmecInputNamelist.ntord + 3][2] = -8.442177055766123E-004;
		rbc[VmecInputNamelist.ntord + 4][2] = -1.591791378279287E-004;
		rbc[VmecInputNamelist.ntord + 5][2] = 2.486332537387373E-005;
		rbc[VmecInputNamelist.ntord + 6][2] = -7.528422072335360E-005;
		rbc[VmecInputNamelist.ntord + 7][2] = -8.021064997689406E-005;
		rbc[VmecInputNamelist.ntord + 8][2] = -5.237201084004603E-005;
		//
		rbc[VmecInputNamelist.ntord - 8][3] = 1.714830561905995E-004;
		rbc[VmecInputNamelist.ntord - 7][3] = 7.083574967044192E-005;
		rbc[VmecInputNamelist.ntord - 6][3] = -1.030190624061461E-004;
		rbc[VmecInputNamelist.ntord - 5][3] = 7.464723877566652E-005;
		rbc[VmecInputNamelist.ntord - 4][3] = 1.594997910406184E-004;
		rbc[VmecInputNamelist.ntord - 3][3] = -1.315594838027009E-004;
		rbc[VmecInputNamelist.ntord - 2][3] = 1.896710074727297E-004;
		rbc[VmecInputNamelist.ntord - 1][3] = -4.186302609439121E-004;
		rbc[VmecInputNamelist.ntord + 0][3] = 3.311501954692046E-004;
		rbc[VmecInputNamelist.ntord + 1][3] = -8.268632528457621E-003;
		rbc[VmecInputNamelist.ntord + 2][3] = -2.035417064634164E-002;
		rbc[VmecInputNamelist.ntord + 3][3] = -1.448421989297009E-002;
		rbc[VmecInputNamelist.ntord + 4][3] = 7.509697962325974E-004;
		rbc[VmecInputNamelist.ntord + 5][3] = -1.328389193189970E-004;
		rbc[VmecInputNamelist.ntord + 6][3] = 1.474276373936834E-004;
		rbc[VmecInputNamelist.ntord + 7][3] = 1.042690595884966E-004;
		rbc[VmecInputNamelist.ntord + 8][3] = 6.708871806270065E-005;
		//
		rbc[VmecInputNamelist.ntord - 8][4] = -1.440375584040593E-004;
		rbc[VmecInputNamelist.ntord - 7][4] = -1.704752360936465E-004;
		rbc[VmecInputNamelist.ntord - 6][4] = 9.011920036675072E-005;
		rbc[VmecInputNamelist.ntord - 5][4] = -1.499685240577500E-004;
		rbc[VmecInputNamelist.ntord - 4][4] = 5.756550703305372E-005;
		rbc[VmecInputNamelist.ntord - 3][4] = 1.900144565741371E-004;
		rbc[VmecInputNamelist.ntord - 2][4] = 4.923760064924362E-005;
		rbc[VmecInputNamelist.ntord - 1][4] = -5.103652305253303E-004;
		rbc[VmecInputNamelist.ntord + 0][4] = 2.486409772891553E-003;
		rbc[VmecInputNamelist.ntord + 1][4] = 3.763697963319822E-003;
		rbc[VmecInputNamelist.ntord + 2][4] = 9.220272047894581E-003;
		rbc[VmecInputNamelist.ntord + 3][4] = 4.017321543601945E-003;
		rbc[VmecInputNamelist.ntord + 4][4] = 9.476230338947471E-004;
		rbc[VmecInputNamelist.ntord + 5][4] = -7.056521343060718E-004;
		rbc[VmecInputNamelist.ntord + 6][4] = -6.013036923002932E-005;
		rbc[VmecInputNamelist.ntord + 7][4] = -8.827308310929046E-005;
		rbc[VmecInputNamelist.ntord + 8][4] = -7.602682766245118E-005;
		//
		rbc[VmecInputNamelist.ntord - 8][5] = 2.153814980670956E-005;
		rbc[VmecInputNamelist.ntord - 7][5] = 1.183952494315462E-005;
		rbc[VmecInputNamelist.ntord - 6][5] = -7.154739133251302E-005;
		rbc[VmecInputNamelist.ntord - 5][5] = -1.003063968491245E-004;
		rbc[VmecInputNamelist.ntord - 4][5] =  1.646390554518521E-004;
		rbc[VmecInputNamelist.ntord - 3][5] = -5.492928459569241E-005;
		rbc[VmecInputNamelist.ntord - 2][5] = -1.008690713410809E-004;
		rbc[VmecInputNamelist.ntord - 1][5] = 5.052245395750897E-004;
		rbc[VmecInputNamelist.ntord + 0][5] = 6.070892885457450E-004;
		rbc[VmecInputNamelist.ntord + 1][5] = 7.467412314074342E-004;
		rbc[VmecInputNamelist.ntord + 2][5] = -1.665044640784653E-003;
		rbc[VmecInputNamelist.ntord + 3][5] = -1.332894131049668E-003;
		rbc[VmecInputNamelist.ntord + 4][5] = 8.621325208882533E-004;
		rbc[VmecInputNamelist.ntord + 5][5] = 2.267011310380025E-004;
		rbc[VmecInputNamelist.ntord + 6][5] = 1.753152728586062E-004;
		rbc[VmecInputNamelist.ntord + 7][5] = 2.258512907309003E-005;
		rbc[VmecInputNamelist.ntord + 8][5] = 4.622663427651913E-005;
		//
		rbc[VmecInputNamelist.ntord - 8][6] = 2.261179964816800E-005;
		rbc[VmecInputNamelist.ntord - 7][6] = -1.811117336074824E-005;
		rbc[VmecInputNamelist.ntord - 6][6] = 7.536505110908006E-005;
		rbc[VmecInputNamelist.ntord - 5][6] = -3.456196403903005E-006;
		rbc[VmecInputNamelist.ntord - 4][6] = -5.139872953188579E-005;
		rbc[VmecInputNamelist.ntord - 3][6] = 3.979845546091016E-005;
		rbc[VmecInputNamelist.ntord - 2][6] = 1.243950000533798E-004;
		rbc[VmecInputNamelist.ntord - 1][6] = -1.106699933734676E-004;
		rbc[VmecInputNamelist.ntord + 0][6] = -1.223947247241717E-004;
		rbc[VmecInputNamelist.ntord + 1][6] = -7.731238353022811E-005;
		rbc[VmecInputNamelist.ntord + 2][6] = -1.569941375256721E-003;
		rbc[VmecInputNamelist.ntord + 3][6] = -1.222191707450016E-003;
		rbc[VmecInputNamelist.ntord + 4][6] = -1.140339566197411E-003;
		rbc[VmecInputNamelist.ntord + 5][6] = -9.518279470355112E-004;
		rbc[VmecInputNamelist.ntord + 6][6] = 2.558431811535456E-004;
		rbc[VmecInputNamelist.ntord + 7][6] = 2.595101629939122E-005;
		rbc[VmecInputNamelist.ntord + 8][6] = 1.620499278874425E-005;
		//
		rbc[VmecInputNamelist.ntord - 8][7] = -1.398163883777544E-005;
		rbc[VmecInputNamelist.ntord - 7][7] = -3.923123999457625E-005;
		rbc[VmecInputNamelist.ntord - 6][7] = -2.306491449291393E-005;
		rbc[VmecInputNamelist.ntord - 5][7] = -1.158175514942671E-004;
		rbc[VmecInputNamelist.ntord - 4][7] = 9.504349582371491E-005;
		rbc[VmecInputNamelist.ntord - 3][7] = -7.789561325605551E-005;
		rbc[VmecInputNamelist.ntord - 2][7] = 3.045433428931733E-005;
		rbc[VmecInputNamelist.ntord - 1][7] = -1.527742119656307E-005;
		rbc[VmecInputNamelist.ntord + 0][7] = -1.998021858917175E-005;
		rbc[VmecInputNamelist.ntord + 1][7] = 4.953044093333177E-005;
		rbc[VmecInputNamelist.ntord + 2][7] = 4.714962252031772E-004;
		rbc[VmecInputNamelist.ntord + 3][7] = 5.784116518161847E-004;
		rbc[VmecInputNamelist.ntord + 4][7] = 7.049240981285555E-004;
		rbc[VmecInputNamelist.ntord + 5][7] = -1.675739529392045E-004;
		rbc[VmecInputNamelist.ntord + 6][7] = -5.586250774777524E-005;
		rbc[VmecInputNamelist.ntord + 7][7] = -2.029459227327149E-004;
		rbc[VmecInputNamelist.ntord + 8][7] = 5.906682990243203E-006;
		//
		Assertions.assertArrayEquals(rbc, vmecInput.rbc);

		double[][] zbs = new double[2 * VmecInputNamelist.ntord + 1][VmecInputNamelist.mpold + 1];
		zbs[VmecInputNamelist.ntord + 0][0] = 0.0;
		zbs[VmecInputNamelist.ntord + 1][0] = -0.231215947095178;
		zbs[VmecInputNamelist.ntord + 2][0] = 3.000448318252753E-004;
		zbs[VmecInputNamelist.ntord + 3][0] = 2.451759201610496E-004;
		zbs[VmecInputNamelist.ntord + 4][0] = 2.002740558876558E-003;
		zbs[VmecInputNamelist.ntord + 5][0] = 2.042500448056660E-004;
		zbs[VmecInputNamelist.ntord + 6][0] = 1.110653375943014E-004;
		zbs[VmecInputNamelist.ntord + 7][0] = 1.233633688773586E-004;
		zbs[VmecInputNamelist.ntord + 8][0] = -5.454884672029991E-004;
		//
		zbs[VmecInputNamelist.ntord - 8][1] = -6.512120712505940E-004;
		zbs[VmecInputNamelist.ntord - 7][1] = -2.473300842408858E-004;
		zbs[VmecInputNamelist.ntord - 6][1] = 1.226721407577620E-003;
		zbs[VmecInputNamelist.ntord - 5][1] = 1.157564964099680E-005;
		zbs[VmecInputNamelist.ntord - 4][1] = 1.948945874541883E-004;
		zbs[VmecInputNamelist.ntord - 3][1] = 1.584625523874076E-003;
		zbs[VmecInputNamelist.ntord - 2][1] = 5.725438140714447E-003;
		zbs[VmecInputNamelist.ntord - 1][1] = 6.183781833751985E-003;
		zbs[VmecInputNamelist.ntord + 0][1] = 0.628787962252855;
		zbs[VmecInputNamelist.ntord + 1][1] = 0.223524141570276;
		zbs[VmecInputNamelist.ntord + 2][1] = 1.323272993170185E-002;
		zbs[VmecInputNamelist.ntord + 3][1] = -6.316519355198045E-004;
		zbs[VmecInputNamelist.ntord + 4][1] = -4.160839711727137E-004;
		zbs[VmecInputNamelist.ntord + 5][1] = -1.431427743599354E-004;
		zbs[VmecInputNamelist.ntord + 6][1] = -3.965264842806391E-004;
		zbs[VmecInputNamelist.ntord + 7][1] = -2.485908172021094E-005;
		zbs[VmecInputNamelist.ntord + 8][1] = 3.549072371917396E-004;
		//
		zbs[VmecInputNamelist.ntord - 8][2] = -3.300731983353361E-004;
		zbs[VmecInputNamelist.ntord - 7][2] = 4.473218222210373E-004;
		zbs[VmecInputNamelist.ntord - 6][2] = -6.065256605704031E-004;
		zbs[VmecInputNamelist.ntord - 5][2] = 3.241308593613608E-004;
		zbs[VmecInputNamelist.ntord - 4][2] = 6.840233808173711E-005;
		zbs[VmecInputNamelist.ntord - 3][2] = -9.056185242554235E-005;
		zbs[VmecInputNamelist.ntord - 2][2] = 1.955219044278804E-004;
		zbs[VmecInputNamelist.ntord - 1][2] = 7.523213388325822E-003;
		zbs[VmecInputNamelist.ntord + 0][2] = -5.215130122068161E-003;
		zbs[VmecInputNamelist.ntord + 1][2] = 2.382604615534300E-002;
		zbs[VmecInputNamelist.ntord + 2][2] = -5.184435337334851E-002;
		zbs[VmecInputNamelist.ntord + 3][2] = 2.258335290509644E-003;
		zbs[VmecInputNamelist.ntord + 4][2] = 6.590696228657478E-004;
		zbs[VmecInputNamelist.ntord + 5][2] = 9.616300338029773E-005;
		zbs[VmecInputNamelist.ntord + 6][2] = 1.018064657002039E-004;
		zbs[VmecInputNamelist.ntord + 7][2] = -6.953157311656576E-005;
		zbs[VmecInputNamelist.ntord + 8][2] = 9.909633701842818E-005;
		//
		zbs[VmecInputNamelist.ntord - 8][3] = -4.519578534113128E-005;
		zbs[VmecInputNamelist.ntord - 7][3] = -8.527405237368145E-005;
		zbs[VmecInputNamelist.ntord - 6][3] = 1.753730817444603E-004;
		zbs[VmecInputNamelist.ntord - 5][3] = -2.326348438307375E-005;
		zbs[VmecInputNamelist.ntord - 4][3] = -2.624975041422071E-004;
		zbs[VmecInputNamelist.ntord - 3][3] = -1.390802867156250E-005;
		zbs[VmecInputNamelist.ntord - 2][3] = 8.074618587295024E-004;
		zbs[VmecInputNamelist.ntord - 1][3] = -1.531636429051088E-003;
		zbs[VmecInputNamelist.ntord + 0][3] = -1.645890760289522E-003;
		zbs[VmecInputNamelist.ntord + 1][3] = -8.199761340258446E-003;
		zbs[VmecInputNamelist.ntord + 2][3] = 5.766395955724331E-003;
		zbs[VmecInputNamelist.ntord + 3][3] = 1.163122094676383E-002;
		zbs[VmecInputNamelist.ntord + 4][3] = -5.968673861601821E-004;
		zbs[VmecInputNamelist.ntord + 5][3] = 1.622756876768911E-004;
		zbs[VmecInputNamelist.ntord + 6][3] = -3.916688872863889E-005;
		zbs[VmecInputNamelist.ntord + 7][3] = -6.313722893293626E-005;
		zbs[VmecInputNamelist.ntord + 8][3] = -3.394782233208609E-005;
		//
		zbs[VmecInputNamelist.ntord - 8][4] = -2.405330978649966E-005;
		zbs[VmecInputNamelist.ntord - 7][4] = 2.084090191207986E-004;
		zbs[VmecInputNamelist.ntord - 6][4] = -6.937710034163128E-005;
		zbs[VmecInputNamelist.ntord - 5][4] = 3.866322821755758E-005;
		zbs[VmecInputNamelist.ntord - 4][4] = 1.296929167587253E-004;
		zbs[VmecInputNamelist.ntord - 3][4] = -1.292084644773919E-004;
		zbs[VmecInputNamelist.ntord - 2][4] = -9.511070777900972E-005;
		zbs[VmecInputNamelist.ntord - 1][4] = -2.460373856685125E-004;
		zbs[VmecInputNamelist.ntord + 0][4] = 4.198477990498741E-003;
		zbs[VmecInputNamelist.ntord + 1][4] = -2.836264233928369E-003;
		zbs[VmecInputNamelist.ntord + 2][4] = 9.239409675099113E-003;
		zbs[VmecInputNamelist.ntord + 3][4] = -4.334921025609129E-003;
		zbs[VmecInputNamelist.ntord + 4][4] = -1.784089978532210E-003;
		zbs[VmecInputNamelist.ntord + 5][4] = 4.675626087533873E-004;
		zbs[VmecInputNamelist.ntord + 6][4] = -3.540448326425835E-005;
		zbs[VmecInputNamelist.ntord + 7][4] = 5.212974445405975E-005;
		zbs[VmecInputNamelist.ntord + 8][4] = 2.676132432009841E-005;
		//
		zbs[VmecInputNamelist.ntord - 8][5] = 3.279177584182842E-005;
		zbs[VmecInputNamelist.ntord - 7][5] = 2.724539949718756E-005;
		zbs[VmecInputNamelist.ntord - 6][5] = -1.593085585738485E-004;
		zbs[VmecInputNamelist.ntord - 5][5] = 2.410957443918596E-004;
		zbs[VmecInputNamelist.ntord - 4][5] = -1.044020785453555E-004;
		zbs[VmecInputNamelist.ntord - 3][5] = 2.056825728184085E-005;
		zbs[VmecInputNamelist.ntord - 2][5] = -3.741888572862706E-005;
		zbs[VmecInputNamelist.ntord - 1][5] = 3.578584050879604E-004;
		zbs[VmecInputNamelist.ntord + 0][5] = 6.752128022657002E-004;
		zbs[VmecInputNamelist.ntord + 1][5] = 1.854756872143505E-003;
		zbs[VmecInputNamelist.ntord + 2][5] = -1.507318481326176E-003;
		zbs[VmecInputNamelist.ntord + 3][5] = -9.597505773106182E-004;
		zbs[VmecInputNamelist.ntord + 4][5] = 9.885927576869936E-004;
		zbs[VmecInputNamelist.ntord + 5][5] = 1.133619037566474E-005;
		zbs[VmecInputNamelist.ntord + 6][5] = -1.325292640508833E-004;
		zbs[VmecInputNamelist.ntord + 7][5] = 2.726109026309649E-005;
		zbs[VmecInputNamelist.ntord + 8][5] = -1.850542092438669E-005;
		//
		zbs[VmecInputNamelist.ntord - 8][6] = 7.386763623731291E-007;
		zbs[VmecInputNamelist.ntord - 7][6] = -2.508309149435569E-005;
		zbs[VmecInputNamelist.ntord - 6][6] = 2.117983834955785E-005;
		zbs[VmecInputNamelist.ntord - 5][6] = -5.753278329735706E-005;
		zbs[VmecInputNamelist.ntord - 4][6] = -1.155249029943377E-005;
		zbs[VmecInputNamelist.ntord - 3][6] = 1.153115758377326E-004;
		zbs[VmecInputNamelist.ntord - 2][6] = -2.830337070611086E-005;
		zbs[VmecInputNamelist.ntord - 1][6] = 1.324540090482543E-005;
		zbs[VmecInputNamelist.ntord + 0][6] = -7.800374644189332E-004;
		zbs[VmecInputNamelist.ntord + 1][6] = 3.838347842484833E-004;
		zbs[VmecInputNamelist.ntord + 2][6] = -6.137053321851298E-004;
		zbs[VmecInputNamelist.ntord + 3][6] = -1.301407081337896E-003;
		zbs[VmecInputNamelist.ntord + 4][6] = -4.955356310734652E-004;
		zbs[VmecInputNamelist.ntord + 5][6] = 3.462627709485832E-004;
		zbs[VmecInputNamelist.ntord + 6][6] = -9.919741284083595E-005;
		zbs[VmecInputNamelist.ntord + 7][6] = -1.784091935229791E-005;
		zbs[VmecInputNamelist.ntord + 8][6] = -2.302543207412257E-005;
		//
		zbs[VmecInputNamelist.ntord - 8][7] = -1.484141014465394E-005;
		zbs[VmecInputNamelist.ntord - 7][7] = -1.816080256001114E-005;
		zbs[VmecInputNamelist.ntord - 6][7] = -3.833112693126000E-006;
		zbs[VmecInputNamelist.ntord - 5][7] = -1.865844947458311E-006;
		zbs[VmecInputNamelist.ntord - 4][7] = -9.801486449695234E-007;
		zbs[VmecInputNamelist.ntord - 3][7] = 8.061675676629167E-006;
		zbs[VmecInputNamelist.ntord - 2][7] = -1.094928404347744E-006;
		zbs[VmecInputNamelist.ntord - 1][7] = -3.601316979928921E-005;
		zbs[VmecInputNamelist.ntord + 0][7] = 2.940524142344392E-004;
		zbs[VmecInputNamelist.ntord + 1][7] = -8.496033274994724E-004;
		zbs[VmecInputNamelist.ntord + 2][7] = 6.255838033712031E-004;
		zbs[VmecInputNamelist.ntord + 3][7] = 3.900613621829464E-004;
		zbs[VmecInputNamelist.ntord + 4][7] = 1.845971251755716E-003;
		zbs[VmecInputNamelist.ntord + 5][7] = 2.917639239820137E-004;
		zbs[VmecInputNamelist.ntord + 6][7] = -5.665554773695777E-005;
		zbs[VmecInputNamelist.ntord + 7][7] = 1.038283071268155E-004;
		zbs[VmecInputNamelist.ntord + 8][7] = -2.979906810937261E-006;
		//
		Assertions.assertArrayEquals(zbs, vmecInput.zbs);

		double[][] rbs = new double[2 * VmecInputNamelist.ntord + 1][VmecInputNamelist.mpold + 1];
		Assertions.assertArrayEquals(rbs, vmecInput.rbs);

		double[][] zbc = new double[2 * VmecInputNamelist.ntord + 1][VmecInputNamelist.mpold + 1];
		Assertions.assertArrayEquals(zbc, vmecInput.zbc);
	}

	@Test
	void vmecInputTest2() throws IOException {

		String inputFilename = "/exampleVmecInput/input.minerva-vmec-9aea06f0c9b27276d9a55c2817eda280";
		InputStream inputFileStream = getClass().getResourceAsStream(inputFilename);

		String inputFile = new String(inputFileStream.readAllBytes(), StandardCharsets.UTF_8);

		VmecInputNamelist vmecInput = new VmecInputNamelist();

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
		vmecInput = (VmecInputNamelist) parser.getParsed();
		assertNotNull(vmecInput);

		Assertions.assertEquals(vmecInput.lasym, false);
		// TODO implement tests for other quantities
	}

	@Test
	void vmecInputTest3() throws IOException {
		InputStream inputFileStream = getClass().getResourceAsStream("/exampleVmecInput/input.demo_vmec6_90");

		String inputFile = new String(inputFileStream.readAllBytes(), StandardCharsets.UTF_8);

		VmecInputNamelist vmecInput = new VmecInputNamelist();

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
		vmecInput = (VmecInputNamelist) parser.getParsed();
		assertNotNull(vmecInput);

		Assertions.assertEquals(vmecInput.lasym, false);
		// TODO implement tests for other quantities
	}

	@Test
	void vmecInputTest4() throws IOException {
		InputStream inputFileStream = getClass().getResourceAsStream("/exampleVmecInput/input.demo3d_vmec6_90");

		String inputFile = new String(inputFileStream.readAllBytes(), StandardCharsets.UTF_8);

		VmecInputNamelist vmecInput = new VmecInputNamelist();

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
		vmecInput = (VmecInputNamelist) parser.getParsed();
		assertNotNull(vmecInput);

		Assertions.assertEquals(vmecInput.lasym, false);
		// TODO implement tests for other quantities
	}

	@Test
	void vmecInputTest5() throws IOException {
		InputStream inputFileStream = getClass().getResourceAsStream("/exampleVmecInput/input.AsMessyAsPossible");

		String inputFile = new String(inputFileStream.readAllBytes(), StandardCharsets.UTF_8);

		VmecInputNamelist vmecInput = new VmecInputNamelist();

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
		vmecInput = (VmecInputNamelist) parser.getParsed();
		assertNotNull(vmecInput);

		Assertions.assertEquals(vmecInput.lasym, false);
		// TODO implement tests for other quantities
	}

	@Test
	void vmecInputTest6() throws IOException {
		InputStream inputFileStream = getClass().getResourceAsStream("/exampleVmecInput/input.ffhr_d1_nflux_v1.1");

		String inputFile = new String(inputFileStream.readAllBytes(), StandardCharsets.UTF_8);

		VmecInputNamelist vmecInput = new VmecInputNamelist();

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
		vmecInput = (VmecInputNamelist) parser.getParsed();
		assertNotNull(vmecInput);

		Assertions.assertEquals(vmecInput.lasym, false);
		// TODO implement tests for other quantities
	}

	@Test
	void vmecInputTest7() throws IOException {
		InputStream inputFileStream = getClass()
				.getResourceAsStream("/exampleVmecInput/input.ITER_nflux_verification_0001.0001");

		String inputFile = new String(inputFileStream.readAllBytes(), StandardCharsets.UTF_8);

		VmecInputNamelist vmecInput = new VmecInputNamelist();

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
		vmecInput = (VmecInputNamelist) parser.getParsed();
		assertNotNull(vmecInput);

		Assertions.assertEquals(vmecInput.lasym, false);
		// TODO implement tests for other quantities
	}

	// => forgot value of am
	@Test
	void vmecInputTest8() throws IOException {

		String inputFilename = "/exampleVmecInput/input.dboe_id_1000_1000_1000_1000_+0000_+0000_v_00_pres_00_it_6";
		InputStream inputFileStream = getClass().getResourceAsStream(inputFilename);

		String inputFile = new String(inputFileStream.readAllBytes(), StandardCharsets.UTF_8);

		VmecInputNamelist vmecInput = new VmecInputNamelist();

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
		vmecInput = (VmecInputNamelist) parser.getParsed();
		assertNotNull(vmecInput);

		Assertions.assertEquals(vmecInput.lasym, false);
		// TODO implement tests for other quantities
	}

	// EXTCUR has spaces in indices
	@Test
	void vmecInputTest9() throws IOException {
		InputStream inputFileStream = getClass().getResourceAsStream("/exampleVmecInput/input.BETA_5_ICUR_5K");

		String inputFile = new String(inputFileStream.readAllBytes(), StandardCharsets.UTF_8);

		VmecInputNamelist vmecInput = new VmecInputNamelist();

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
		vmecInput = (VmecInputNamelist) parser.getParsed();
		assertNotNull(vmecInput);

		Assertions.assertEquals(vmecInput.lasym, false);
		// TODO implement tests for other quantities
	}

	@Test
	void vmecInputTest10() throws IOException {
		InputStream inputFileStream = getClass().getResourceAsStream("/exampleVmecInput/input.JDHtest7");

		String inputFile = new String(inputFileStream.readAllBytes(), StandardCharsets.UTF_8);

		VmecInputNamelist vmecInput = new VmecInputNamelist();

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
		vmecInput = (VmecInputNamelist) parser.getParsed();
		assertNotNull(vmecInput);

		Assertions.assertEquals(vmecInput.lasym, false);
		// TODO implement tests for other quantities
	}

	@Test
	void vmecInputTest11() throws IOException {
		InputStream inputFileStream = getClass().getResourceAsStream("/exampleVmecInput/input.test.vmec_asym");

		String inputFile = new String(inputFileStream.readAllBytes(), StandardCharsets.UTF_8);

		VmecInputNamelist vmecInput = new VmecInputNamelist();

		// parse the namelist into a Java object
		FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
		vmecInput = (VmecInputNamelist) parser.getParsed();
		assertNotNull(vmecInput);

		Assertions.assertEquals(vmecInput.lasym, true);
		// TODO implement tests for other quantities
	}
}
