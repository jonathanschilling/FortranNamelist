package de.labathome;

import java.util.Arrays;

/**
 * INDATA namelist class for VMEC MHD equilibrium solver.
 * This class is the most exact representation of the variables in the INDATA namelist,
 * including array shapes and initial values outside the original sources.
 */
public class VmecIndataNamelist {

	/** maximum number of poloidal harmonics (in r,z,lam fourier series) */
	final public static int mpold = 101;

	/** maximum number of toroidal harmonics */
	final public static int ntord = 101;

	/** 4000 / 1001 / ... some of these is the actual limit... */
	final public static int nsd = 10001;

	/** 4000 / 1001 / ... some of these is the actual limit... */
	final public static int ndatafmax = nsd;

	/** default poloidal Fourier resolution */
	final public static int mpol_default = 6;

	/** default toroidal Fourier resolution */
	final public static int ntor_default = 0;

	/** default number of flux surfaces */
	final public static int ns_default = 31;

	/** number of external current groups; from LIBSTELL/Sources/Modules/vsvd0.f */
	final public static int nigroup = 300;

	/**
	 * number of external b-field loop sets allowed; from
	 * LIBSTELL/Sources/Modules/vsvd0.f
	 */
	final public static int nbsetsp = 5;

	/**
	 * number of external poloidal flux loops; from LIBSTELL/Sources/Modules/vsvd0.f
	 */
	final public static int nfloops = 100;

	/**
	 * number of external b-field coils per set; from
	 * LIBSTELL/Sources/Modules/vsvd0.f
	 */
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




	// Numerical Resolution and Symmetry Assumptions

	/** =T, run in asymmetric mode; =F, run in stellarator symmetry mode */
	@namelist_variable
	public boolean lasym;

	/** number of toroidal field periods ( =1 for Tokamak) */
	@namelist_variable
	public int nfp;

	/** upper limit for poloidal mode numbers */
	@namelist_variable
	public int mpol;

	/** upper limit for toroidal mode number range */
	@namelist_variable
	public int ntor;

	/** number of theta grid points (>=2*mpol+6) */
	@namelist_variable
	public int ntheta;

	/** number of zeta grid points (=1 IF ntor=0) */
	@namelist_variable
	public int nzeta;


	// Multi-Grid Steps

	/** array of radial mesh sizes to be used in multigrid sequence */
	@namelist_variable
	public int[] ns_array;

	/** array of value of residual(s) at which each multigrid iteration ends */
	@namelist_variable
	public double[] ftol_array;

	/** array of number of iterations (used to terminate run) at each multigrid step */
	@namelist_variable
	public int[] niter_array;


	// Global Physics Parameters

	/** toroidal flux enclosed by plasma at edge (in Wb) */
	@namelist_variable
	public double phiedge;

	/** flux conserving (=0) or prescribed toroidal current (=1) */
	@namelist_variable
	public int ncurr;


	// Profile of Mass or Pressure

	/**
	 * Specifies parameterization type of pmass function.
	 * <table style="border-collapse: collapse; border: 1px solid black;">
	 * <thead>
	 * <tr>
	 * <th style="border: 1px solid black;">pmass_type</th>
	 * <th style="border: 1px solid black;">description</th>
	 * </tr>
	 * </thead> <tbody>
	 * <tr>
	 * <td style="border: 1px solid black;">gauss_trunc</td>
	 * <td style="border: 1px solid black;">Truncated gaussian. p(s)=am(0)
	 * (exp(-(s/am(1)) ** 2) - exp(-(1/am(1)) ** 2))</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">two_power</td>
	 * <td style="border: 1px solid black;">am(0)*(1 - s^am(1))^am(2)</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">two_power_gs</td>
	 * <td style="border: 1px solid black;">am(0)*((1 - s^am(1))^am(2)*(1 +
	 * Sum_{i}am(i)*exp(-((s - am(i + 1))/am(i + 2))^2))</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">two_Lorentz</td>
	 * <td style="border: 1px solid black;">two Lorentz function</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">Akima_spline</td>
	 * <td style="border: 1px solid black;">Akima splines.</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">cublic_spline</td>
	 * <td style="border: 1px solid black;">Cublic splines.</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">pedestal</td>
	 * <td style="border: 1px solid black;">Pedestal profile.</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">rational</td>
	 * <td style="border: 1px solid black;">Rational function (ratio of
	 * polynomials).</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">line_segment</td>
	 * <td style="border: 1px solid black;">Line segments.</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">power_series</td>
	 * <td style="border: 1px solid black;">Power series. p(s)=Sum[ am(j) s ** j]
	 * (Default)</td>
	 * </tr>
	 * </tbody>
	 * </table>
	 */
	@namelist_variable
	public String pmass_type;

	/**
	 * mass or pressure (gamma=0) expansion coefficients (series in s) in MKS units
	 * [NWT/M**2]; Interpretation changes with pmass_type
	 */
	@namelist_variable(dim0min = 0)
	public double[] am;

	/** Auxiliary array for mass profile. Used for splines, s values */
	@namelist_variable
	public double[] am_aux_s;

	/** Auxiliary array for mass profile. Used for splines, function values */
	@namelist_variable
	public double[] am_aux_f;

	/**
	 * actor used to scale pressure profile (default value = 1) useful so user can
	 * fix profile and change beta without having to change all AM coefficients
	 * separately
	 */
	@namelist_variable
	public double pres_scale;

	/** value of compressibility index (gamma=0 => pressure prescribed) */
	@namelist_variable
	public double gamma;

	/**
	 * pressure value of pedestal (finite value of pressure at and beyond the LCFS)
	 */
	@namelist_variable
	public double spres_ped;


	// (Initial Guess for) Rotational Transform Profile

	/**
	 * Specifies parameterization type of piota function.
	 * <table style="border-collapse: collapse; border: 1px solid black;">
	 * <thead>
	 * <tr>
	 * <th style="border: 1px solid black;">piota_type</th>
	 * <th style="border: 1px solid black;">description</th>
	 * </tr>
	 * </thead> <tbody>
	 * <tr>
	 * <td style="border: 1px solid black;">sum_atan</td>
	 * <td style="border: 1px solid black;">Sum of arctangents.</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">Akima_spline</td>
	 * <td style="border: 1px solid black;">Akima splines.</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">cublic_spline</td>
	 * <td style="border: 1px solid black;">Cublic splines.</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">rational</td>
	 * <td style="border: 1px solid black;">Rational function (ratio of
	 * polynomials).</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">nice_quadratic</td>
	 * <td style="border: 1px solid black;">Nice quadratic function.</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">line_segment</td>
	 * <td style="border: 1px solid black;">Line segments.</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">power_series</td>
	 * <td style="border: 1px solid black;">Power series. iota(s)=Sum[ ai(j) s ** j]
	 * (Default)</td>
	 * </tr>
	 * </tbody>
	 * </table>
	 */
	@namelist_variable
	public String piota_type;

	/**
	 * expansion coefficients for iota (power series in s) used when ncurr=0;
	 * Interpretation changes with piota_type
	 */
	@namelist_variable(dim0min = 0)
	public double[] ai;

	/** Auxiliary array for iota profile. Used for splines, s values */
	@namelist_variable
	public double[] ai_aux_s;

	/** Auxiliary array for iota profile. Used for splines, function values */
	@namelist_variable
	public double[] ai_aux_f;


	// (Initial Guess for) Toroidal Current Profile

	/**
	 * Specifies parameterization type of pcurr function.
	 * <table style="border-collapse: collapse; border: 1px solid black;">
	 * <thead>
	 * <tr>
	 * <th style="border: 1px solid black;">pcurr_type</th>
	 * <th style="border: 1px solid black;">description</th>
	 * </tr>
	 * </thead> <tbody>
	 * <tr>
	 * <td style="border: 1px solid black;">sum_cossq_s</td>
	 * <td style="border: 1px solid black;">Sum of cos^2 waves. I'</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">sum_cossq_sqrts</td>
	 * <td style="border: 1px solid black;">Sum of cos^2 with respect to sqrt(s). I'
	 * </td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">sum_cossq_s_free</td>
	 * <td style="border: 1px solid black;">Sum of cos^2 waves up to 7. Free
	 * position and width. I'</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">gauss_trunc</td>
	 * <td style="border: 1px solid black;">Truncated gaussian. I'(s)=ac(0)
	 * (exp(-(s/ac(1)) ** 2) - exp(-(1/ac(1)) ** 2))</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">two_power</td>
	 * <td style="border: 1px solid black;">ac(0)*(1 - s^ac(1))^ac(2)</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">two_power_gs</td>
	 * <td style="border: 1px solid black;">ac(0)*((1 - s^ac(1))^ac(2)*(1 +
	 * Sum_{i}ac(i)*exp(-((s - ac(i + 1))/ac(i + 2))^2))</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">sum_atan</td>
	 * <td style="border: 1px solid black;">Sum of arctangents.</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">power_series_I</td>
	 * <td style="border: 1px solid black;">Power series.</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">Akima_spline_I</td>
	 * <td style="border: 1px solid black;">Akima splines.</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">Akima_spline_Ip</td>
	 * <td style="border: 1px solid black;">Akima splines. I'</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">cublic_splines_I</td>
	 * <td style="border: 1px solid black;">Cublic splines.</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">cublic_splines_Ip</td>
	 * <td style="border: 1px solid black;">Cublic splines. I'</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">pedestal</td>
	 * <td style="border: 1px solid black;">Pedestal profile.</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">rational</td>
	 * <td style="border: 1px solid black;">Rational function (ratio of
	 * polynomials).</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">line_segment_Ip</td>
	 * <td style="border: 1px solid black;">Line segments. I'</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">line_segment_I</td>
	 * <td style="border: 1px solid black;">Line segments.</td>
	 * </tr>
	 * <tr>
	 * <td style="border: 1px solid black;">power_series</td>
	 * <td style="border: 1px solid black;">Power series. I'(s)=Sum[ ac(j) s ** j]
	 * (Default)</td>
	 * </tr>
	 * </tbody>
	 * </table>
	 */
	@namelist_variable
	public String pcurr_type;

	/**
	 * expansion coefficients for the normalized (pcurr(s=1) = 1) radial derivative
	 * of the flux-averaged toroidal current density (power series in s) used when
	 * ncurr=1; Interpretation changes with pcurr_type
	 */
	@namelist_variable(dim0min = 0)
	public double[] ac;

	/** Auxiliary array for current profile. Used for splines, s values */
	@namelist_variable
	public double[] ac_aux_s;

	/** Auxiliary array for current profile. Used for splines, function values */
	@namelist_variable
	public double[] ac_aux_f;

	/**
	 * value of toroidal current [A]. Used if ncurr = 1 to specify current profile,
	 * or IF in data reconstruction mode.
	 */
	@namelist_variable
	public double curtor;

	/** scaling factor for phiedge: bloat plasma (bloat>=1) */
	@namelist_variable
	public double bloat;


	// Free-Boundary Parameters

	/** =T, run in free boundary mode if mgrid_file exists */
	@namelist_variable
	public boolean lfreeb;

	/** full path for vacuum Green's function data */
	@namelist_variable
	public String mgrid_file;

	/**
	 * array of currents in each external current group. Used to multiply Green's
	 * function for fields and loops read in from MGRID file. Should use real
	 * current units (A).
	 */
	@namelist_variable
	public double[] extcur;

	/**
	 * number of iteration steps between accurate calculation of vacuum response;
	 * use fast interpolation scheme in between
	 */
	@namelist_variable
	public int nvacskip;


	// Tweaking Parameters

	/** number of timesteps between printouts on screen */
	@namelist_variable
	public int nstep;

	/**
	 * Expansion coefficients (power series) for the redistribution of the radial
	 * positions of the flux surfaces. This can be used to refine grid in center or
	 * at LCFS for example.
	 */
	@namelist_variable
	public double[] aphi;

	/** time derivative in force iterations */
	@namelist_variable
	public double delt;

	/** weight factor for constraint force (=1 by DEFAULT) */
	@namelist_variable
	public double tcon0;

	/**
	 * =T, use non-variational forces to ensure <EQUIF> = 0; =F, use variational
	 * form of forces, <EQUIF> ~ 0
	 */
	@namelist_variable
	public boolean lforbal;


	// Initial Guess for Magnetic Axis Geometry

	/** stellarator-symmetric coefficients of magnetic axis position r */
	@namelist_variable(dim0min = 0)
	public double[] raxis_cc; //

	/** stellarator-symmetric coefficients of magnetic axis position z */
	@namelist_variable(dim0min = 0)
	public double[] zaxis_cs; //

	/** stellarator-asymmetric coefficients of magnetic axis position r */
	@namelist_variable(dim0min = 0)
	public double[] raxis_cs; //

	/** stellarator-asymmetric coefficients of magnetic axis position z */
	@namelist_variable(dim0min = 0)
	public double[] zaxis_cc; //


	// (Initial Guess for) Boundary Geometry


	/** boundary coefficients of COS(m*theta-n*zeta) for R [m] */
	@namelist_variable(dim0min = -ntord, dim1min = 0)
	public double[][] rbc;

	/** boundary coefficients of SIN(m*theta-n*zeta) for Z [m] */
	@namelist_variable(dim0min = -ntord, dim1min = 0)
	public double[][] zbs;

	/** boundary coefficients of SIN(m*theta-n*zeta) for R [m] */
	@namelist_variable(dim0min = -ntord, dim1min = 0)
	public double[][] rbs;

	/** boundary coefficients of COS(m*theta-n*zeta) for Z [m] */
	@namelist_variable(dim0min = -ntord, dim1min = 0)
	public double[][] zbc;




	// ----------------------------------
	// ----------------------------------
	// ----------------------------------
	// Bring these back eventually

	/** inserted M.Drevlak */
	@namelist_variable
	@Deprecated
	public boolean lgiveup;

	/** inserted M.Drevlak, giveup-factor for ftolv */
	@namelist_variable
	@Deprecated
	public double fgiveup;

	/** J Hanson See jxbforce coding */
	@namelist_variable
	@Deprecated
	public boolean lbsubs;

	/**
	 * specifies type of 2D preconditioner to use: <br>
	 * 'default', diagonal in m,n, tri-diagonal in s <br>
	 * 'conjugate-gradient', block tri-di, evolve using cg method <br>
	 * 'gmres', block tri-di, generalized minimal residual method <br>
	 * 'tfqmr' block tri-di, transpose-free quasi minimum residual
	 */
	@namelist_variable
	@Deprecated
	public String precon_type;

	/**
	 * value of preconditioned force residuals at which block (2d) tri-di solver is
	 * turned on, if requested via type_prec2d
	 */
	@namelist_variable
	@Deprecated
	public double prec2d_threshold;

	/** ??? reverse-field pinch ??? */
	@namelist_variable
	@Deprecated
	public boolean lrfp;






	// ----------------------------------
	// ----------------------------------
	// ----------------------------------
	// Below variables are gone for good.

	/** ??? */
	@namelist_variable
	@Deprecated
	public double time_slice;

	/** *** */
	@namelist_variable
	@Deprecated
	public int nsin;

	/** ??? probably a single limit for number of iterations */
	@namelist_variable
	@Deprecated
	public int niter;

	/** ??? probably a single convergence criterion */
	@namelist_variable
	@Deprecated
	public double ftol;

	/** integer specifier for current profile (deprecated) */
	@namelist_variable
	@Deprecated
	public int ipcurr;

	/** integer specifier for mass/pressure profile (deprecated) */
	@namelist_variable
	@Deprecated
	public int ipmass;

	/** integer specifier for iota profile (deprecated) */
	@namelist_variable
	@Deprecated
	public int ipiota;

	/** WAC (anisotropic pres) */
	@namelist_variable(dim0min = 0)
	@Deprecated
	public double[] ah;

	/** WAC (anisotropic pres) */
	@namelist_variable(dim0min = 0)
	@Deprecated
	public double[] at;

	/** WAC (anisotropic pres) */
	@namelist_variable
	@Deprecated
	public double bcrit;

	/** ??? */
	@namelist_variable
	@Deprecated
	public int mfilter_fbdy;

	/** ??? */
	@namelist_variable
	@Deprecated
	public int nfilter_fbdy;

	/** ??? */
	@namelist_variable
	@Deprecated
	public double sigma_current;

	/** ??? */
	@namelist_variable
	@Deprecated
	public int omp_num_threads;

	/** ??? */
	@namelist_variable
	@Deprecated
	public double[] psa;

	/** ??? */
	@namelist_variable
	@Deprecated
	public double[] pfa;

	/** ??? */
	@namelist_variable
	@Deprecated
	public double[] isa;

	/** ??? */
	@namelist_variable
	@Deprecated
	public double[] ifa;

	/**
	 * parameter for reconstruction: <br>
	 * 0: USE pressure profile width to determine PHIEDGE <br>
	 * 1: match value of PHIEDGE in input file (default) <br>
	 * 2: USE LIMPOS data (in mgrid file) to find PHIEDGE <br>
	 * 3: USE Ip to find PHIEDGE (fixed-boundary only)
	 */
	@namelist_variable
	@Deprecated
	public int imatch_phiedge;

	/** ??? */
	@namelist_variable
	@Deprecated
	public int opt_raxis;

	/** spline tension for iota */
	@namelist_variable
	@Deprecated
	public double tensi;

	/** spline tension for pressure profile */
	@namelist_variable
	@Deprecated
	public double tensp;

	/** uniform EXPerimental offset of MSE data (calibration offset) ... PLUS ... */
	@namelist_variable
	@Deprecated
	public double mseangle_offset;

	/** multiplier on mseprof offset array (calibration offset) */
	@namelist_variable
	@Deprecated
	public double mseangle_offsetm;

	/**
	 * number of Motional Stark effect data points: <br>
	 * >0, USE mse data to find iota <br>
	 * <=0, fixed iota profile ai
	 */
	@namelist_variable
	@Deprecated
	public int imse;

	/**
	 * number of iota spline points (computed internally unless specified
	 * explicitly)
	 */
	@namelist_variable
	@Deprecated
	public int isnodes;

	/** ??? */
	@namelist_variable
	@Deprecated
	public double[] rstark;

	/** pitch angle data from stark measurement */
	@namelist_variable
	@Deprecated
	public double[] datastark;

	/** standard deviation (degrees) in MSE data */
	@namelist_variable
	@Deprecated
	public double[] sigma_stark;

	/**
	 * number of pressure profile data points: = 0, no thompson scattering data to
	 * READ
	 */
	@namelist_variable
	@Deprecated
	public int itse;

	/**
	 * number of pressure spline points (computed internally unless specified
	 * explicitly)
	 */
	@namelist_variable
	@Deprecated
	public int ipnodes;

	/** number by which Thomson scattering data is scaled to get actual pressure */
	@namelist_variable
	@Deprecated
	public double presfac;

	/** uniform arbitrary radial offset of pressure data */
	@namelist_variable
	@Deprecated
	public double pres_offset;

	/** ??? */
	@namelist_variable
	@Deprecated
	public double[] rthom;

	/** pressure data from Thompson, CHEERS (Pa) */
	@namelist_variable
	@Deprecated
	public double[] datathom;

	/** standard deviation (Pa) for pressure profile data */
	@namelist_variable
	@Deprecated
	public double[] sigma_thom;

	/** diamagnetic toroidal flux (Wb) */
	@namelist_variable
	@Deprecated
	public double phidiam;

	/** ??? */
	@namelist_variable
	@Deprecated
	public double sigma_delphid;

	/** vbl spline tension for iota */
	@namelist_variable
	@Deprecated
	public double tensi2;

	/**
	 * vbl spline tension form factor (note: IF tensi!=tensi2 THEN tension(i-th
	 * point) = tensi+(tensi2-tensi)*(i/n-1))**fpolyi
	 */
	@namelist_variable
	@Deprecated
	public double fpolyi;

	/** number of flux loop measurements used in matching */
	@namelist_variable
	@Deprecated
	public int nflxs;

	/** array giving INDEX of flux measurement in iconnect array */
	@namelist_variable
	@Deprecated
	public int[] indxflx;

	/**
	 * measured flux loop signals corresponding to the combination of signals in
	 * iconnect array
	 */
	@namelist_variable
	@Deprecated
	public double[] dsiobt;

	/** standard deviaton (Wb) for EXTERNAL poloidal flux data */
	@namelist_variable
	@Deprecated
	public double[] sigma_flux;

	/** number of selected EXTERNAL bfield measurements in set n from nml file */
	@namelist_variable
	@Deprecated
	public int[] nbfld;

	/** array giving INDEX of bfield measurement used in matching */
	@namelist_variable
	@Deprecated
	public int[][] indxbfld;

	/** Backwards compatibility: Obsolete --> see raxis_cc, raxis_cs */
	@namelist_variable(dim0min = 0)
	@Deprecated
	public double[] raxis;

	/** Backwards compatibility: Obsolete --> see raxis_cc, raxis_cs */
	@namelist_variable(dim0min = 0)
	@Deprecated
	public double[] zaxis;

	/**
	 * measured magnetic field at rbcoil(m,n),zbcoil(m,n) at the orientation
	 * br*COS(abcoil) + bz*SIN(abcoil)
	 */
	@namelist_variable
	@Deprecated
	public double[][] bbc;

	/** standard deviation (T) for EXTERNAL magnetic field data */
	@namelist_variable
	@Deprecated
	public double[][] sigma_b;

	/**
	 * =.true. IF pressure data are prescribed in REAL space; =.false. IF data in
	 * flux space.
	 */
	@namelist_variable
	@Deprecated
	public boolean lpofr;



	/** ??? */
	@namelist_variable
	@Deprecated
	public boolean lmove_axis;

	/** ??? reconstruction mode ??? */
	@namelist_variable
	@Deprecated
	public boolean lrecon;

	/** ??? */
	@namelist_variable
	@Deprecated
	public boolean lmac;

	/** ??? */
	@namelist_variable
	@Deprecated
	public boolean ledge_dump;

	/** Obsolete */
	@namelist_variable
	@Deprecated
	public boolean lspectrum_dump;

	/** Obsolete */
	@namelist_variable
	@Deprecated
	public boolean loptim;

	/** to obtain old fort.8 file */
	@namelist_variable
	@Deprecated
	public boolean loldout;

	/** J.Geiger: for txt- and diagno-output */
	@namelist_variable
	@Deprecated
	public boolean lwouttxt;

	/** J.Geiger: for txt- and diagno-output */
	@namelist_variable
	@Deprecated
	public boolean ldiagno;

	/** J.Geiger: to force full 3D1-output */
	@namelist_variable
	@Deprecated
	public boolean lfull3d1out;

	/**
	 * maximum number of iterations of the main loop, i.e. of all values in ns_array
	 */
	@namelist_variable
	@Deprecated
	public int max_main_iterations;

	/**
	 * Init variables (especially arrays!) to default values/sizes in the constructor.
	 * Updated based on trunk/LIBSTELL/Sources/Modules/{vmec_input.f, vparams.f} and
	 * trunk/VMEC2000/Sources/Input_Output/readin.f
	 * from the stellinstall repository of ORNL at
	 * https://github.com/ORNL-Fusion/stellinstall
	 */
	public VmecIndataNamelist() {

		// Numerical Resolution and Symmetry Assumptions
		lasym = false;
		nfp = 1;
		mpol = mpol_default;
		ntor = ntor_default;
		ntheta = 0;
		nzeta = 0;

		// Multi-Grid Steps
		ns_array = new int[100];
		ns_array[0] = ns_default;

		ftol_array = new double[100];
		ftol_array[0] = 1.e-10;

		niter_array = new int[100];
		Arrays.fill(niter_array, -1);
		niter_array[0] = 100;

		// Global Physics Parameters
		phiedge = 1;
		ncurr = 0;

		// Profile of Mass or Pressure
		pmass_type = "power_series";
		am = new double[21];
		am_aux_s = new double[ndatafmax];
		am_aux_f = new double[ndatafmax];
		pres_scale = 1;
		gamma = 0;
		spres_ped = 1;

		// (Initial Guess for) Rotational Transform Profile
		piota_type = "power_series";
		ai = new double[21];
		ai_aux_s = new double[ndatafmax];
		ai_aux_f = new double[ndatafmax];

		// (Initial Guess for) Toroidal Current Profile
		pcurr_type = "power_series";
		ac = new double[21];
		ac_aux_s = new double[ndatafmax];
		ac_aux_f = new double[ndatafmax];
		curtor = 0;
		bloat = 1;

		// Free-Boundary Parameters
		lfreeb = false;
		mgrid_file = "NONE";
		extcur = new double[nigroup];
		nvacskip = 1;

		// Tweaking Parameters
		nstep = 10;
		aphi = new double[20];
		aphi[0] = 1;
		delt = 1;
		tcon0 = 1;
		lforbal = false; // SPH: changed 05-14-14

		// Initial Guess for Magnetic Axis
		raxis_cc = new double[ntord + 1];
		zaxis_cs = new double[ntord + 1];
		raxis_cs = new double[ntord + 1];
		zaxis_cc = new double[ntord + 1];

		// (Initial Guess for) Boundary Geometry
		rbc = new double[2 * ntord + 1][mpold + 1]; // actual max indices are [-ntord:ntord][0:mpold]
		zbs = new double[2 * ntord + 1][mpold + 1];
		rbs = new double[2 * ntord + 1][mpold + 1];
		zbc = new double[2 * ntord + 1][mpold + 1];

		// -------------------------------
		// Bring these back eventually

		// inserted M.Drevlak
		lgiveup = false;

		// inserted M.Drevlak
		fgiveup = 3.E+01;

		// J Hanson. See jxbforce coding
		lbsubs = false;

		precon_type = "NONE";
		prec2d_threshold = 1.e-30;

		lrfp = false;

		// -------------------------------
		// Below are gone for good

		raxis = new double[ntord + 1];
		zaxis = new double[ntord + 1];

		ipcurr = 0;
		ipiota = 0;
		ipmass = 0;

		nsin = ns_default;
		niter = 100;
		ftol = 1.e-10;


		ah = new double[21];
		at = new double[21];

		psa = new double[ndatafmax];
		pfa = new double[ndatafmax];
		isa = new double[ndatafmax];
		ifa = new double[ndatafmax];

		nbfld = new int[nbsetsp];
		indxflx = new int[nfloops];

		rstark = new double[nmse];
		datastark = new double[nmse];
		sigma_stark = new double[nmse];

		rthom = new double[ntse];
		datathom = new double[ntse];
		sigma_thom = new double[ntse];

		dsiobt = new double[nfloops];
		sigma_flux = new double[nfloops];

		indxbfld = new int[nbcoilsp][nbsetsp];
		bbc = new double[nbcoilsp][nbsetsp];
		sigma_b = new double[nbcoilsp][nbsetsp];

		omp_num_threads = 8;

		time_slice = 0;
		mfilter_fbdy = -1;
		nfilter_fbdy = -1;
		lmove_axis = true;
		lmac = false;

		// obtain old fort.8 file
		loldout = false;

		// J Geiger 2010-05-04: for txt- and diagno-output
		ldiagno = false;

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

		// ANISTROPY PARAMETERS
		bcrit = 1;
		at[0] = 1;
	}

	/**
	 * Merge non-zero values given for raxis and zaxis into raxis_cc and zaxis_cs.
	 * This is needed, since raxis and zaxis are deprecated but still might be
	 * around in some old input files. You will probably need this only in case you
	 * work with some old input files...
	 */
	public void sanitize() {

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
		// ==> If something is given in raxis or zaxis, merge them with raxis_cc and zaxis_cs.
		// ==> Call VmecIndataNamelist.sanitize() to accomplish this.

		for (int i = 0; i < raxis.length; ++i) {
			if (raxis[i] != 0.0) {
				raxis_cc[i] = raxis[i];
			}
			if (zaxis[i] != 0.0) {
				zaxis_cs[i] = zaxis[i];
			}
		}

		// TODO: use nsin, niter, ftol if ns_array, ftol_array, niter_array are not set

		// TODO: use ipmass, ipiota, ipcurr if p****_type are not set
	}


	/**
	 * Parse a namelist given in text form.
	 * @param namelist text form of namelist to parse
	 * @return parsed object
	 */
	public static VmecIndataNamelist fromString(String namelist) {

		VmecIndataNamelist indataNamelist = new VmecIndataNamelist();
		FortranNamelist parser = new FortranNamelist(namelist, "indata", indataNamelist);
		indataNamelist = (VmecIndataNamelist) parser.getParsed();

		indataNamelist.sanitize();

		return indataNamelist;
	}
}
