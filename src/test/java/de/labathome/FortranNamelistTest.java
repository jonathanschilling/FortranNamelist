package de.labathome;

import static org.junit.jupiter.api.Assertions.assertNotNull;
import java.lang.reflect.Field;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import org.junit.jupiter.api.Test;

/**
 * Tests for a class to parse Fortran namelists into Java classes.
 * This uses the information given in the class definition to constrain the parser.
 * @author Jonathan Schilling (jonathan.schilling@mail.de)
 * @version 0.9.0 2018-07-25 intial implementation
 * @version 1.0.1 2018-09-11 1d and 2d arrays with fancy specifiers
 */
public class FortranNamelistTest {

    @Test
    public final void basicTest() {


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

        System.out.println("before parsing:");
        dumpObject(testClass);
        System.out.println("\n\n");


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

        // parse the namelist into a Java object
        FortranNamelist parser = new FortranNamelist(namelist, "mgrid_nli", testClass);
        testClass = (MgridParameters)parser.getParsed();

        System.out.println("\n\n");
        System.out.println("after parsing:");
        dumpObject(testClass);

        assertNotNull(testClass);
    }

    @Test
    public final void vmecInputTest() {

        /**
         * Prototype class as target for parsing.
         * The default is to look for variables named the same as the variable name in the namelist.
         * Optionally, you can supply a name in the annotation to define a name used in the namelist.
         */
        class VmecInputNamelist {
            final public static int mpold = 101;      // maximum number of poloidal harmonics (in r,z,lam fourier series)
            final public static int ntord = 101;      // maximum number of toroidal harmonics
            final public static int nsd = 10001; // 4000; // 1001 // some of these is the actual limit...
            final public static int ndatafmax  = nsd;

            final public static int mpol_default = 6;
            final public static int ntor_default = 0;
            final public static int ns_default   = 31;

            // from LIBSTELL/Sources/Modules/vsvd0.f
            final public static int nigroup = 300;    // number of external current groups
            final public static int nbsetsp = 5;      // number of external b-field loop sets allowed
            final public static int nfloops = 100;    // number of external poloidal flux loops
            final public static int nbcoilsp = 100;   // number of external b-field coils per set

            final public static int nmse = 100;       // number of mse measurements
            final public static int ntse = 100;       // number of thompson scattering measurements

            // These input variables are mentioned in readin.f, but not in the namelist reading routine anymore.
            // Maybe they have been abandoned in the meantime?
            // ???      phifac:   factor scaling toroidal flux to match apres or limiter
            // ???      pknots:   array of pressure knot values in SQRT(s) space
            // ???      sknots:   array of iota knot values in SQRT(s) space
            // ???     mseprof:   offset array from NAMELIST MSEPROFIL                        \
            // ???                so that the total offset on the i-th MSE data point is      | different namelist:
            // ???                taken to be                                                 | mseprofile --> ignore here
            // ???                = mseangle_offset+mseangle_offsetM*mseprof(i)               /
            // ???        nobd:   number of connected flux loop measurements
            // ???      nobser:   number of individual flux loop positions
            // ???      nbsets:   number of B-coil sets defined in mgrid file
            // ???  nbcoils(n):   number of bfield coils in each set defined in mgrid file
            // ???    nbcoilsn:   total number of bfield coils defined in mgrid file
            // ??? rbcoil(m,n):   R position of the m-th coil in the n-th set from mgrid file
            // ??? zbcoil(m,n):   Z position of the m-th coil in the n-th set from mgrid file
            // ??? abcoil(m,n):   orientation (surface normal wrt R axis; in radians)
            // ???                of the m-th coil in the n-th set from mgrid file.
            // ???       nflxs:   number of flux loop measurements used in matching
            // ???    nbfld(n):   number of selected EXTERNAL bfield measurements in set n from nml file
            // ???      nbfldn:   total number of EXTERNAL bfield measurements used in matching
            // ??? sigma_current:  standard deviation (A) in toroidal current
            // ??? sigma_delphid:  standard deviation (Wb) for diamagnetic match

            // NOTE: FOR STANDARD DEVIATIONS (sigma''s) < 0, INTERPRET
            //       AS PERCENT OF RESPECTIVE MEASUREMENT
            //       THE INPUT DATA FILE WILL ACCEPT BOTH POSITIVE AND NEGATIVE
            //       SIGMAS, WHICH IT INTERPRETS DIFFERENTLY. FOR SIGMA > 0, IT
            //       TAKES SIGMA TO BE THE STANDARD DEVIATION FOR THAT MEASUREMENT
            //       AS DESCRIBED ABOVE. FOR SIGMA < 0, SIGMA IS INTERPRETED AS
            //       THE FRACTION OF THE MEASURED DATA NEEDED TO COMPUTE THE ABSOLUTE
            //       SIGMA, I.E., (-SIGMA * DATA) = ACTUAL SIGMA USED IN CODE.

            // Commends are based on trunk/VMEC2000/Input_Output/readin.f of the ORNL stellinstall repository.
            // Namelist variables and names are from trunk/LIBSTELL/Sources/Modules/vmec_input.f (same repo).
            @namelist_variable  String     mgrid_file;          // full path for vacuum Green's function data
            @namelist_variable  double     time_slice;          // ???
            @namelist_variable     int     nfp;                 // number of toroidal field periods ( =1 for Tokamak)
            @namelist_variable     int     ncurr;               // flux conserving (=0) or prescribed toroidal current (=1)
            @namelist_variable     int     nsin;                // ???
            @namelist_variable     int     niter;               // ??? probably a single limit for number of iterations
            @namelist_variable     int     nstep;               // number of timesteps between printouts on screen
            @namelist_variable     int     nvacskip;            // number of iteration steps between accurate calculation of vacuum response; use fast interpolation scheme in between
            @namelist_variable  double     delt;                // time derivative in force iterations
            @namelist_variable  double     ftol;                // ??? probably a single convergence criterion
            @namelist_variable  double     gamma;               // value of compressibility index (gamma=0 => pressure prescribed)
            @namelist_variable(dim0min=0)
            double[]   am;                  // mass or pressure (gamma=0) expansion coefficients (series in s) in MKS units [NWT/M**2]; Interpretation changes with pmass_type
            @namelist_variable(dim0min=0)
            double[]   ai;                  // expansion coefficients for iota (power series in s) used when ncurr=0; Interpretation changes with piota_type
            @namelist_variable(dim0min=0)
            double[]   ac;                  // expansion coefficients for the normalized (pcurr(s=1) = 1) radial derivative of the flux-averaged toroidal current density (power series in s) used when ncurr=1; Interpretation changes with pcurr_type
            @namelist_variable  double[]   aphi;                // expansion coefficients for radial coordinate; can be used to refine grid in center or at LCFS for example
            @namelist_variable  String     pcurr_type;          // Specifies parameterization type of pcurr function: 'power_series' - I'(s)=Sum[ ac(j) s ** j] - Default; 'gauss_trunc'  - I'(s)=ac(0) (exp(-(s/ac(1)) ** 2) - exp(-(1/ac(1)) ** 2)); others - see function pcurr
            @namelist_variable  String     pmass_type;          // Specifies parameterization type of pmass function: 'power_series' - p(s)=Sum[ am(j) s ** j] - Default; 'gauss_trunc'  - p(s)=am(0) (exp(-(s/am(1)) ** 2) - exp(-(1/am(1)) ** 2)); others - see function pmass
            @namelist_variable  String     piota_type;          // Specifies parameterization type of piota function: 'power_series' - p(s)=Sum[ am(j) s ** j] - Default; others - see function piota
            @namelist_variable  double[]   am_aux_s;            // Auxiliary array for mass profile. Used for splines, s values
            @namelist_variable  double[]   am_aux_f;            // Auxiliary array for mass profile. Used for splines, function values
            @namelist_variable  double[]   ai_aux_s;            // Auxiliary array for iota profile. Used for splines, s values
            @namelist_variable  double[]   ai_aux_f;            // Auxiliary array for iota profile. Used for splines, function values
            @namelist_variable  double[]   ac_aux_s;            // Auxiliary array for current profile. Used for splines, s values
            @namelist_variable  double[]   ac_aux_f;            // Auxiliary array for current profile. Used for splines, function values
            @namelist_variable(dim0min=0)
            double[]   ah;                  // WAC (anisotropic pres)
            @namelist_variable(dim0min=0)
            double[]   at;                  // WAC (anisotropic pres)
            @namelist_variable  double     bcrit;               // WAC (anisotropic pres)
            @namelist_variable(dim0min=-ntord, dim1min=0)
            double[][] rbc;                 // boundary coefficients of COS(m*theta-n*zeta) for R [m]
            @namelist_variable(dim0min=-ntord, dim1min=0)
            double[][] zbs;                 // boundary coefficients of SIN(m*theta-n*zeta) for Z [m]
            @namelist_variable(dim0min=-ntord, dim1min=0)
            double[][] rbs;                 // boundary coefficients of SIN(m*theta-n*zeta) for R [m]
            @namelist_variable(dim0min=-ntord, dim1min=0)
            double[][] zbc;                 // boundary coefficients of COS(m*theta-n*zeta) for Z [m]
            @namelist_variable  double     spres_ped;           // pressure value of pedestal (finite value of pressure at and beyond the LCFS)
            @namelist_variable  double     pres_scale;          // factor used to scale pressure profile (default value = 1) useful so user can fix profile and change beta without having to change all AM coefficients separately
            @namelist_variable(dim0min=0)
            double[]   raxis_cc;            // stellarator-symmetric coefficients of magnetic axis position r
            @namelist_variable(dim0min=0) 
            double[]   zaxis_cs;            // stellarator-symmetric coefficients of magnetic axis position z
            @namelist_variable(dim0min=0) 
            double[]   raxis_cs;            // stellarator-asymmetric coefficients of magnetic axis position r
            @namelist_variable(dim0min=0)
            double[]   zaxis_cc;            // stellarator-asymmetric coefficients of magnetic axis position z
            @namelist_variable     int     mpol;                // upper limit for poloidal mode numbers 
            @namelist_variable     int     ntor;                // upper limit for toroidal mode number range
            @namelist_variable     int     ntheta;              // number of theta grid points (>=2*mpol+6)
            @namelist_variable     int     nzeta;               // number of zeta grid points (=1 IF ntor=0)
            @namelist_variable     int     mfilter_fbdy;        // ???
            @namelist_variable     int     nfilter_fbdy;        // ???
            @namelist_variable     int[]   niter_array;         // array of number of iterations (used to terminate run) at each multigrid iteration
            @namelist_variable     int[]   ns_array;            // array of radial mesh sizes to be used in multigrid sequence
            @namelist_variable  double[]   ftol_array;          // array of value of residual(s) at which each multigrid iteration ends
            @namelist_variable  double     tcon0;               // weight factor for constraint force (=1 by DEFAULT)
            @namelist_variable  String     precon_type;         // specifies type of 2D preconditioner to use ('default', diagonal in m,n, tri-diagonal in s; 'conjugate-gradient', block tri-di, evolve using cg method; 'gmres', block tri-di, generalized minimal residual method; 'tfqmr', block tri-di, transpose-free quasi minimum residual
            @namelist_variable  double     prec2d_threshold;    // value of preconditioned force residuals at which block (2d) tri-di solver is turned on, if requested via type_prec2d
            @namelist_variable  double     curtor;              // value of toroidal current [A]. Used if ncurr = 1 to specify current profile, or IF in data reconstruction mode.
            @namelist_variable  double     sigma_current;       // ???
            @namelist_variable  double[]   extcur;              // array of currents in each external current group. Used to multiply Green's function for fields and loops read in from MGRID file. Should use real current units (A).
            @namelist_variable     int     omp_num_threads;     // ???
            @namelist_variable  double     phiedge;             // toroidal flux enclosed by plasma at edge (in Wb)
            @namelist_variable  double[]   psa;                 // ???
            @namelist_variable  double[]   pfa;                 // ???
            @namelist_variable  double[]   isa;                 // ???
            @namelist_variable  double[]   ifa;                 // ???
            @namelist_variable     int     imatch_phiedge;      // = 1 (default), match value of PHIEDGE in input file; = 0, USE pressure profile width to determine PHIEDGE;  = 2, USE LIMPOS data (in mgrid file) to find PHIEDGE; = 3, USE Ip to find PHIEDGE (fixed-boundary only)
            @namelist_variable     int     opt_raxis;           // ???
            @namelist_variable  double     tensi;               // spline tension for iota
            @namelist_variable  double     tensp;               // spline tension for pressure profile
            @namelist_variable  double     mseangle_offset;     // uniform EXPerimental offset of MSE data (calibration offset) ... PLUS ...
            @namelist_variable  double     mseangle_offsetm;    //  multiplier on mseprof offset array (calibration offset)
            @namelist_variable     int     imse;                //  number of Motional Stark effect data points: >0, USE mse data to find iota; <=0, fixed iota profile ai
            @namelist_variable     int     isnodes;             // number of iota spline points (computed internally unless specified explicitly)
            @namelist_variable  double[]   rstark;              // ???
            @namelist_variable  double[]   datastark;           // pitch angle data from stark measurement
            @namelist_variable  double[]   sigma_stark;         // standard deviation (degrees) in MSE data
            @namelist_variable     int     itse;                // number of pressure profile data points: = 0, no thompson scattering data to READ
            @namelist_variable     int     ipnodes;             // number of pressure spline points (computed internally unless specified explicitly)
            @namelist_variable  double     presfac;             // number by which Thomson scattering data is scaled to get actual pressure
            @namelist_variable  double     pres_offset;         // uniform arbitrary  radial offset of pressure data
            @namelist_variable  double[]   rthom;               // ???
            @namelist_variable  double[]   datathom;            // pressure data from Thompson, CHEERS (Pa)
            @namelist_variable  double[]   sigma_thom;          // standard deviation (Pa) for pressure profile data
            @namelist_variable  double     phidiam;             // diamagnetic toroidal flux (Wb)
            @namelist_variable  double     sigma_delphid;       // ???
            @namelist_variable  double     tensi2;              // vbl spline tension for iota
            @namelist_variable  double     fpolyi;              // vbl spline tension form factor (note: IF tensi!=tensi2 THEN tension(i-th point) = tensi+(tensi2-tensi)*(i/n-1))**fpolyi
            @namelist_variable     int     nflxs;               // number of flux loop measurements used in matching
            @namelist_variable     int[]   indxflx;             // array giving INDEX of flux measurement in iconnect array
            @namelist_variable  double[]   dsiobt;              // measured flux loop signals corresponding to the combination of signals in iconnect array
            @namelist_variable  double[]   sigma_flux;          // standard deviaton (Wb) for EXTERNAL poloidal flux data
            @namelist_variable     int[]   nbfld;               // number of selected EXTERNAL bfield measurements in set n from nml file
            @namelist_variable     int[][] indxbfld;            // array giving INDEX of bfield measurement used in matching
            @namelist_variable  double     bloat;               // scaling factor for phiedge: bloat plasma (bloat>1)
            @namelist_variable(dim0min=0)
            double[]   raxis;               // Backwards compatibility: Obsolete --> see raxis_cc, raxis_cs
            @namelist_variable(dim0min=0)
            double[]   zaxis;               // Backwards compatibility: Obsolete --> see raxis_cc, raxis_cs
            @namelist_variable  double[][] bbc;                 // measured magnetic field at rbcoil(m,n),zbcoil(m,n) at the orientation br*COS(abcoil) + bz*SIN(abcoil)
            @namelist_variable  double[][] sigma_b;             // standard deviation (T) for EXTERNAL magnetic field data
            @namelist_variable boolean     lpofr;               // LOGICAL variable. =.true. IF pressure data are prescribed in REAL space. =.false. IF data in flux space.
            @namelist_variable boolean     lforbal;             // =T, use non-variational forces to ensure <EQUIF> = 0; =F, use variational form of forces, <EQUIF> ~ 0
            @namelist_variable boolean     lfreeb;              // =T, run in free boundary mode if mgrid_file exists
            @namelist_variable boolean     lmove_axis;          // ???
            @namelist_variable boolean     lrecon;              // ??? reconstruction mode ???
            @namelist_variable boolean     lmac;                // ???
            @namelist_variable boolean     lasym;               // =T, run in asymmetric mode; =F, run in stellarator symmetry mode
            @namelist_variable boolean     ledge_dump;          // ???
            @namelist_variable boolean     lspectrum_dump;      // Obsolete
            @namelist_variable boolean     loptim;              // Obsolete
            @namelist_variable boolean     lrfp;                // ??? reverse-field pinch ???
            @namelist_variable boolean     loldout;             // to obtain old fort.8 file
            @namelist_variable boolean     lwouttxt;            // J.Geiger: for txt- and diagno-output
            @namelist_variable boolean     ldiagno;             // J.Geiger: for txt- and diagno-output
            @namelist_variable boolean     lfull3d1out;         // J.Geiger: to force full 3D1-output
            @namelist_variable     int     max_main_iterations; // maximum number of iterations of the main loop, i.e. of all values in ns_array
            @namelist_variable boolean     lgiveup;             // inserted M.Drevlak
            @namelist_variable  double     fgiveup;             // inserted M.Drevlak, giveup-factor for ftolv
            @namelist_variable boolean     lbsubs;              // J Hanson See jxbforce coding

            // Init variables (especially arrays!) to default values/sizes in the constructor.
            // Updated based on trunk/LIBSTELL/Sources/Modules/{vmec_input.f, vparams.f} and trunk/VMEC2000/Sources/Input_Output/readin.f
            // from the stellinstall repository of ORNL at https://github.com/ORNL-Fusion/stellinstall
            public VmecInputNamelist() {
                ns_array    = new int[100];
                niter_array = new int[100];

                rbc = new double[2*ntord+1][mpold+1]; // actual indices are in [-ntord:ntord][0:mpold]
                zbs = new double[2*ntord+1][mpold+1];
                rbs = new double[2*ntord+1][mpold+1];
                zbc = new double[2*ntord+1][mpold+1];

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

                raxis_cc = new double[ntord+1]; // actual indices are in [0:ntord]
                zaxis_cs = new double[ntord+1];
                raxis_cs = new double[ntord+1];
                zaxis_cc = new double[ntord+1];
                raxis = new double[ntord+1];
                zaxis = new double[ntord+1];

                ftol_array = new double[100];

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
                mfilter_fbdy = -1; nfilter_fbdy = -1;
                tcon0 = 1;
                precon_type = "NONE"; prec2d_threshold = 1.e-30;
                curtor = 0; 
                phiedge = 1;
                mgrid_file = "NONE";
                lfreeb = true;
                lmove_axis = true;
                lmac = false;
                lforbal = false;        // SPH: changed 05-14-14
                lasym = false;
                lrfp = false;
                loldout = false;        // obtain old fort.8 file
                ldiagno = false;        // J Geiger 2010-05-04: for txt- and diagno-output
                lgiveup = false;        // inserted M.Drevlak
                fgiveup = 3.E+01;       // inserted M.Drevlak
                lbsubs = false;         // J Hanson. See jxbforce coding
                lfull3d1out = true;     // J Geiger & SPH (5-21-15): to force full 3D1-output
                max_main_iterations = 1;  // to keep a presumably expected standard behavior.
                //              #if defined(NETCDF)
                lwouttxt = false;       // to keep functionality as expected with netcdf
                // J.Geiger: for txt- and diagno-output
                //              #else
                //                    lwouttxt = .true.        ! and without netcdf
                //              #endif
                // J Geiger 2010-05-04 end
                pcurr_type = "power_series";
                piota_type = "power_series";
                pmass_type = "power_series";

                //     ANISTROPY PARAMETERS
                bcrit = 1;
                at[0] = 1;

                //
                // BACKWARDS COMPATIBILITY
                //
                // Work around a bug in gfortran. When performing an optimized build, the WHERE
                // statement would produce incorrect results. Work around this bug by expanding
                // the full WHERE statment. This should have no adverse effects on any other
                // compiler since these statements are equivalent to the older code statement.
                //
                //    WHERE (raxis .ne. 0.0_dp) raxis_cc = raxis
                //    WHERE (zaxis .ne. 0.0_dp) zaxis_cs = zaxis
                //
                // The effect of this bug optimized to code to effectively ignore the WHERE
                // statement and assign all value values of the r/zaxis to the r/zaxis_cc/s
                // arrays. Explicitly adding the r/zaxis .ne. 0.0_dp section prevents this. This
                // bug is known to exist in gfortran 4.9. It may manifest in other versions.
                //    WHERE (raxis .ne. 0.0_dp)
                //       raxis_cc = raxis
                //    ELSEWHERE
                //       raxis_cc = raxis_cc
                //    ENDWHERE
                //    WHERE (zaxis .ne. 0.0_dp)
                //       zaxis_cs = zaxis
                //    ELSEWHERE
                //       zaxis_cs = zaxis_cs
                //    ENDWHERE
                // 
                // ==> If something is given in raxis or zaxis, put it in the corresponding positions of raxis_cc and zaxis_cs.
            }
        }

        String[] testInputs = new String[] {
            "/exampleVmecInput/input.v3fit_jsch_20171207_006_2200ms_pm_10ms_v0_Bprecheck"
            , "/exampleVmecInput/input.minerva-vmec-9aea06f0c9b27276d9a55c2817eda280"
            , "/exampleVmecInput/input.demo_vmec6_90"
            , "/exampleVmecInput/input.demo3d_vmec6_90"
            , "/exampleVmecInput/input.AsMessyAsPossible"
            , "/exampleVmecInput/input.ffhr_d1_nflux_v1.1"
            , "/exampleVmecInput/input.ITER_nflux_verification_0001.0001"
            , "/exampleVmecInput/input.dboe_id_1000_1000_1000_1000_+0000_+0000_v_00_pres_00_it_6" // => forgot value of am
            , "/exampleVmecInput/input.BETA_5_ICUR_5K" // EXTCUR has spaces in indices
        };

        for (String testVmecInput: testInputs) {

            System.out.println("testing \"" + testVmecInput + "\"");
            
            String inputFile = "";
            try {
                inputFile = new String(Files.readAllBytes(Paths.get(getClass().getResource(testVmecInput).toURI())));
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            //System.out.println("read input file:");
            //System.out.println(inputFile);

             

            //        System.out.println("before parsing:");
            //        dumpObject(vmecInput);
            //        System.out.println("\n\n");

            // parse the namelist into a Java object
            VmecInputNamelist vmecInput = new VmecInputNamelist();
            FortranNamelist parser = new FortranNamelist(inputFile, "indata", vmecInput);
            vmecInput = (VmecInputNamelist)parser.getParsed();

            //        System.out.println("\n\n");
            //        System.out.println("after parsing:");
            //        dumpObject(vmecInput);

            assertNotNull(vmecInput);
        }
    }


    public static void dumpObject(Object obj) {
        System.out.println("class dump of \""+obj.toString()+"\" :");
        for (Field field : obj.getClass().getDeclaredFields()) {
            // You might want to set modifier to public first.
            field.setAccessible(true); 

            try {
                Object value = field.get(obj);
                if (value != null) {
                    if (value.getClass().equals(String.class)) {
                        System.out.println(" "+field.getName() + " = \"" + value+"\"");                        
                    } else {
                        System.out.println(" "+field.getName() + " = " + value);
                    }
                } else {
                    System.out.println(" "+field.getName() + " = null");
                }                
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

}
