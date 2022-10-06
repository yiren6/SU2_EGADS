// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#ifndef __NASTRAN_CARDS_H__
#define __NASTRAN_CARDS_H__

#include "capsTypes.h"
#include "feaTypes.h" // Bring in FEA structures
#include "miscUtils.h"

#ifdef __cplusplus
extern "C" {
#endif


/*
 * AELIST
 * 
 * sid     Set identification number (Integer > 0)
 * e       List of aerodynamic boxes generated by CAERO entries (Integer > 0)
 * 
 * numE    Number of `e` values
 */
int nastranCard_aelist(FILE *fp, int *sid, int numE, int *e,
                       feaFileTypeEnum formatType);

/*
 * AERO
 *
 * acsid   Aerodynamic coordinate system identification (Integer >= 0)
 * velocityVelocity for aerodynamic force data recovert and to calculate the
 *         BOV parameter (Real)
 * refc    Reference length for reduced frequency (Real)
 * rhoref  Reference density (Real)
 * symxz   Symmetry key for the aero coordinate x-z plane 
 *         (Integer, one of: (+1, 0, -1), Default = 0)
 * symxy   Symmetry key for the aero coordinate x-y plane 
 *         (Integer, one of: (+1, 0, -1), Default = 0)
 */
int nastranCard_aero(FILE *fp, const int *acsid, const double *velocity, const double *refc,
                     const double *rhoref, const int *symxz, const int *symxy,
                     feaFileTypeEnum formatType);

/*
 * AEROS
 * 
 * acsid   Aerodynamic coordinate system identification (Integer > 0)
 * rcsid   Reference coordinate identification for rigid body motions
 *         (Integer > 0)
 * refc    Reference chord length (Real > 0.0)
 * refb    Reference span (Real > 0.0)
 * refs    Reference wing area (Real > 0.0)
 * symxz   Symmetry key for the aero coordinate x-z plane 
 *         (Integer, one of: (+1, 0, -1), Default = 0)
 * symxy   Symmetry key for the aero coordinate x-y plane 
 *         (Integer, one of: (+1, 0, -1), Default = 0)
 */
int nastranCard_aeros(FILE *fp, const int *acsid, const int *rcsid, const double *refc,
                      const double *refb, const double *refs, const int *symxz,
                      const int *symxy, feaFileTypeEnum formatType);

/*
 * AESURF
 * 
 * id      Controller identification number (Integer > 0) [Ignored]
 * label   Controller name (Character)
 * cid     Identification number of rectangular coordinate system with
 *         a y-axis that defines the hinge line of the control surface
 *         component (Integer > 0)
 * alid    Identification of an AELIST bulk data entry that identifies
 *         all aerodynamic elements that make up the control surface
 *         component (Integer > 0)
 * eff     Control surface effectiveness (Real != 0; Default = 1.0)
 * ldw     Linear downwash flag (Either "LDW" or "NOLDW")
 * crefc   Reference chord length for the control surface 
 *         (Real >= 0.0, Default = 1.0)
 * crefs   Reference surface area for the control surface
 *         (Real >= 0.0, Default = 1.0)
 * pllim   Lower deflection limit for the control surface
 *         in radian (Real, Default = +- pi/2)
 * pulim   Upper deflection limit for the control surface
 *         in radian (Real, Default = +- pi/2)
 * hmllim  Lower hinge moment limit for the control surface
 *         in force-length units (Real, Default = Blank (no limit))
 * hmulim  Upper hinge moment limit for the control surface
 *         in force-length units (Real, Default = Blank (no limit))
 * tqllim  Set identification number of TABLED entry that provides the
 *         lower deflection limits for the control surface as a function
 *         of the dynamic pressure (Integer > 0, Default = Blank (no limit))
 * tqulim  Set identification number of TABLED entry that provides the
 *         upper deflection limits for the control surface as a function
 *         of the dynamic pressure (Integer > 0, Default = Blank (no limit))
 */
int nastranCard_aesurf(FILE *fp, int *id, char *label, int *cid, int *alid,
                       double *eff, char *ldw, int *crefc, int *crefs, 
                       double *pllim, double *pulim, double *hmllim, double *hmulim,
                       int *tqllim, int *tqulim, feaFileTypeEnum formatType);

/*
 * CAERO1
 *
 * eid     Element identification number (Integer > 0)
 * pid     Property identification number of a PAERO1 entry (Integer > 0)
 * cp      Coordinate system for location points 1 and 4 (Integer >= 0,
 *         or blank)
 * nspan   Number of spanwise boxes (Integer >= 0, or blank)
 * nchord  Number of chordwise boxes (Integer >= 0, or blank)
 * lspan   ID of an AEFACT data entry containing a list of division points
 *         for spanwise boxes. Only used if `nspan` is blank. (Integer > 0)
 * lchord  ID of an AEFACT data entry containing a list of division points
 *         for chordwise boxes. Only used if `nchord` is blank. (Integer > 0)
 * igid    Inference group identification (Integer > 0)
 * xyz1    
 * xyz4    Locations of point 1 and 4, in coordinate system `cp` (Real)
 * x12
 * x43     Edge chord lengths in aerodynamic coordinate system (Real >= 0.0,
 *         but not both zero)
 */
int nastranCard_caero1(FILE *fp, const int *eid, const int *pid, const int *cp,
                       const int *nspan, const int *nchord, const int *lspan, const int *lchord,
                       const int *igid, const double xyz1[3], const double xyz4[3],
                       const double *x12, const double *x43,
                       feaFileTypeEnum formatType);

/*
 * CBAR
 */
int nastranCard_cbar(FILE *fp, const int *eid, const int *pid, const int g[2],
                     const double x[3], const int *g0, const int *pa, const int *pb,
                     const double wa[3], const double wb[3],
                     feaFileTypeEnum formatType);

/*
 * CDAMP2
 * 
 * eid     Unique element identification number (Integer > 0)
 * b       Value of the scalar damper (Real)
 * g1
 * g2      Geometric grid point identification number (Integer >= 0,
 *         or blank)
 * c1
 * c2      Component number(0 <= Integer <= 6, or blank)
 */
int nastranCard_cdamp2(FILE *fp, const int *eid, const double *b, const int *g1, const int *g2,
                       const int *c1, const int *c2, feaFileTypeEnum formatType);

/*
 * CELAS2
 *
 * eid     Unique element identification number (Integer > 0)
 * k       Stiffness of the scalar spring (Real)
 * g1
 * g2      Geometric grid point or scalar identification number 
 *         (Integer >= 0, or blank)
 * c1
 * c2      Component number(0 <= Integer <= 6, or blank)
 * ge      Damping coefficient (Real, or blank)
 * s       Stress coefficient (Real)
 */
int nastranCard_celas2(FILE *fp, const int *eid, const double *k, const int *g1, const int *g2,
                       const int *c1, const int *c2, const double *ge, const double *s,
                       feaFileTypeEnum formatType);

/*
 * CMASS2
 *
 * eid     Unique element identification number (Integer > 0)
 * m       Value of scalar mass (Real)
 * g1
 * g2      Geometric grid point or scalar identification number 
 *         (Integer >= 0, or blank)
 * c1
 * c2      Component number(0 <= Integer <= 6, or blank)
 */
int nastranCard_cmass2(FILE *fp, const int *eid, const double *m, const int *g1, const int *g2,
                       const int *c1, const int *c2, feaFileTypeEnum formatType);

/*
 * CONM2
 */
int nastranCard_conm2(FILE *fp, const int *eid, const int *g, const int *cid, const double *m,
                      const double x[3], const double i[6],
                      feaFileTypeEnum formatType);

/*
 * CORD2R
 */
int nastranCard_cord2c(FILE *fp, const int *cid, const int *rid,
                       double a[3], double b[3], double c[3],
                       feaFileTypeEnum formatType);

/*
 * CORD2R
 */
int nastranCard_cord2r(FILE *fp, const int *cid, const int *rid,
                       double a[3], double b[3], double c[3],
                       feaFileTypeEnum formatType);

/*
 * CORD2R
 */
int nastranCard_cord2s(FILE *fp, const int *cid, const int *rid,
                       double a[3], double b[3], double c[3],
                       feaFileTypeEnum formatType);

/*
 * CQUAD4
 */
int nastranCard_cquad4(FILE *fp, int *eid, int *pid, int g[4], 
                       double *theta, int *mcid, double *zoffs,
                       double t[4], feaFileTypeEnum formatType);

/*
 * CQUAD8
 */
int nastranCard_cquad8(FILE *fp, int *eid, int *pid, int g[8], 
                       double *theta, int *mcid, double *zoffs,
                       double t[4], feaFileTypeEnum formatType);

/*
 * CSHEAR
 */
int nastranCard_cshear(FILE *fp, int *eid, int *pid, int g[4],
                       feaFileTypeEnum formatType);

/*
 * CTRIA3
 */
int nastranCard_ctria3(FILE *fp, int *eid, int *pid, int g[3], 
                       double *theta, int *mcid, double *zoffs,
                       double t[3], feaFileTypeEnum formatType);

/*
 * CTRIA6
 */
int nastranCard_ctria6(FILE *fp, int *eid, int *pid, int g[6], 
                       double *theta, int *mcid, double *zoffs,
                       double t[3], feaFileTypeEnum formatType);

/*
 * DCONADD
 */
int nastranCard_dconadd(FILE *fp, const int *dcid, int numDC, const int *dc,
                        feaFileTypeEnum formatType);

/*
 * DCONSTR
 *
 * dcid    Design constraint set identification number (Integer > 0)
 * rid     DRESPi entry identification number (Integer > 0)
 * lallow  Lower bound on the response quantity (Real, Default = -1.E20)
 * uallow  Upper bound on the response quantity (Real, Default =  1.E20)
 */
int nastranCard_dconstr(FILE *fp, const int *dcid, const int *rid,
                        const double *lallow, const double *uallow,
                        feaFileTypeEnum formatType);

/*
 * DDVAL
 */
int nastranCard_ddval(FILE *fp, const int *id, int numDV, const double *dval,
                      feaFileTypeEnum formatType);

/*
 * DEQATN
 * 
 * eqid    Unique equation identification number (Integer > 0)
 * equationEquation(s). Array of equation strings no longer than 64
 *         chars, first one no longer than 56 chars)
 * 
 * numEquation  Number of `equation` values
 */
int nastranCard_deqatn(FILE *fp, const int *eqid, int numEquation,
                       char **equation);

/*
 * DESVAR
 *
 * id      Unique design variable identification number (Integer > 0)
 * label   User-supplied name for printing purposes (Character)
 * xinit   Initial value (Real, `xlb` <= `xinit` <= `xub`)
 * xlb     Lower bound (Real, default = -1.E20)
 * xub     Upper bound (Real, default =  1.E20)
 * delxv   Fractional change allowed for the design variable during
 *         approximate optimization (Real >= 0, or blank)
 * ddval   DDVAL identification number (Integer > 0)
 */
int nastranCard_desvar(FILE *fp, const int *id, const char *label, const double *xinit,
                       const double *xlb, const double *xub, const double *delxv, const int *ddval,
                       feaFileTypeEnum formatType);

/*
 * DLINK
 *
 * id      Unique identification number (Integer > 0)
 * ddvid   Dependent design variable identification (Integer > 0)
 * c0      Constant term (Real, Default = 0.0)
 * cmult   Constant multiplier (Real, default = 1.0)
 * idv     Independent design variable identification numbers (Integer > 0)
 * c       Coefficients corresponding to `idv` (Real)
 *
 * numDV   Number of design variables to link
 */
int nastranCard_dlink(FILE *fp, const int *id, const int *ddvid, const double *c0,
                      const double *cmult, const int numDV, const int *idv, const double *c,
                      feaFileTypeEnum formatType);

/*
 * DMI
 * 
 * name    Name of the matrix (Character)
 * form    Form of matrix (Integer)
 * tin     Type of matrix being input (Integer)
 * tout    Type of matrix being output (Integer)
 * m       Number of rows (Integer > 0)
 * n       Number of columns (Integer > 0)
 * a       Real values (Real)
 * b       Imaginary values (Real)
 */
int nastranCard_dmi(FILE *fp, char *name, int *form, 
                    int *tin, int* tout, int m, int n, 
                    double *a, double *b,
                    feaFileTypeEnum formatType);

/*
 * DOPTPRM
 * 
 * param   Names of the design optimization parameters
 * val     Values of the paremeters
 * 
 * numParam  Number of params
 * paramType  Types of the values (Either 1: Integer or 2: Double)
 */
int nastranCard_doptprm(FILE *fp, int numParam, char **param,
                        int *paramType, void **val,
                        feaFileTypeEnum formatType);

/*
 * DRESP1
 *
 * id      Unique entry identifier (Integer > 0)
 * label   User-defined label (Character)
 * rtype   Response type (Character)
 * ptype   Element flag or property entry name (Character)
 * region  Region identifier for constraint screening (Integer > 0, or blank)
 * atta
 * attb
 * atts     Response attributes (Integer > 0 or Real or blank)
 * 
 * attaType  The type of `atta` entry (1: Integer or 2: Double)
 * attbType  The type of `attb` entry (1: Integer or 2: Double)
 * attsType  The type of `atts` entries (1: Integer or 2: Double)
 * numAtts  The number of `atts` values
 */
int nastranCard_dresp1(FILE *fp, const int *id, const char *label, const char *rtype,
                       const char *ptype, const int *region, int attaType, const void *atta,
                       int attbType, const void *attb,
                       int attsType, int numAtts, const void *atts,
                       feaFileTypeEnum formatType);

/*
 * DRESP2
 * 
 * id      Unique identification number (Integer > 0)
 * label   User-defined label (Character)
 * eqid    DEQATN entry identification number (Integer > 0)
 * region  Region identifier for constant screening (Integer > 0)
 * dvid    DESVAR entry identification numbers
 * labl    Labels for a constant in the DTABLE entry
 * nr      DRESP1 entry identification numbers
 * g       Grid point identification numbers
 * c       Degree of freedom numbers of grid points `g`
 * 
 * numDV   Number of `dvid` values
 * numLabl Number of `labl` values
 * numNR   Number of `nr` values
 * numG    Number of `g` and `c` values
 */
int nastranCard_dresp2(FILE *fp, int *id, char *label, int *eqid,
                       int *region, int numDV, int *dvid, int numLabl,
                       char **labl, int numNR, int *nr, int numG, 
                       int *g, int *c, int numNRR, int *nrr,
                       feaFileTypeEnum formatType);

/*
 * DTABLE
 * 
 * labl    Labels for the constants
 * valu    Values of the constants
 * 
 * numVal  Number of constants
 */
int nastranCard_dtable(FILE *fp, int numVal, char **labl, double *valu,
                       feaFileTypeEnum formatType);

/*
 * DVCREL1
 */
int nastranCard_dvcrel1(FILE *fp, const int *id, const char *type, const int *eid,
                        const char *cpname, const double *cpmin, const double *cpmax,
                        const double *c0, int numDV, const int *dvid, const double *coeff,
                        feaFileTypeEnum formatType);

/*
 * DVMREL1
 */
int nastranCard_dvmrel1(FILE *fp, const int *id, const char *type, const int *mid,
                        const char *mpname, const double *mpmin, const double *mpmax,
                        const double *c0, int numDV, const int *dvid, const double *coeff,
                        feaFileTypeEnum formatType);

/*
 * DVPREL1
 * 
 * id      Unique identification number (Integer > 0)
 * type    Name of property entry (Character)
 * pid     Property entry identification number (Integer > 0)
 * fid     Field position of the property entry, or word position in the
 *         element property table of the analysis model (Integer != 0)
 * pname   Property entry name, used if `fid` is NULL (Character)
 * pmin    Minimum value allowed for this property (Real, or blank)
 * pmax    Maximum value allowed for this property (Real, or blank)
 * c0      Constant term of relation (Real, Default = 1.0E20)
 * dvid    DESVAR entry identification numbers (Integer > 0)
 * coef    Coefficients of linear relation (Real)
 *
 * numDV   Number of design variables
 */
int nastranCard_dvprel1(FILE *fp, const int *id, const char *type, const int *pid,
                        const int *fid, const char *pname, const double *pmin,
                        const double *pmax, const double *c0,
                        int numDV, const int *dvid, const double *coef,
                        feaFileTypeEnum formatType);

/*
 * EIGR
 * 
 * sid     Set identification number (Integer > 0)
 * method  Method of eigenvalue extraction (Character)
 * f1      
 * f2      Frequency range of interest (Real, or blank)
 * ne      Estimate of number of roots in range (Integer > 0, or blank)
 * nd      Desired number of roots (Integer >= 0, or blank)
 * norm    Method of normalizing eigenvectors (Character)
 * g       Grid or scalar point identification number (Integer > 0,
 *         or blank)
 * c       Component number (1 <= Integer <= 6, or blank)
 */
int nastranCard_eigr(FILE *fp, const int *sid, const char *method, const double *f1,
                     const double *f2, const int *ne, const int *nd, const char *norm,
                     const int *g, const int *c, feaFileTypeEnum formatType);

/*
 * EIGRL
 */
int nastranCard_eigrl(FILE *fp, const int *sid, const double *v1, const double *v2,
                      const int *nd, int *msglvl, int *maxset,
                      double *shfscl, char *norm,
                      feaFileTypeEnum formatType);

/*
 * FLFACT
 */
int nastranCard_flfact(FILE *fp, const int *sid, int numF, const double *f,
                       feaFileTypeEnum formatType);

/*
 * FLUTTER
 */
int nastranCard_flutter(FILE *fp, const int *sid, const char *method, const int *dens,
                        const int *mach, const int *rfreq, const char *imeth,
                        const int *nvalue, const double *eps,
                        feaFileTypeEnum formatType);

/*
 * FORCE
 */
int nastranCard_force(FILE *fp, const int *sid, const int *g, const int *cid,
                      const double *f, const double n[3],
                      feaFileTypeEnum formatType);

/*
 * GRAV
 */
int nastranCard_grav(FILE *fp, const int *sid, const int *cid,
                     const double *g, const double n[3],
                     feaFileTypeEnum formatType);

/*
 * LOAD
 */
int nastranCard_load(FILE *fp, int *sid, double *s, int numL, 
                     double *ls, int *l, feaFileTypeEnum formatType);

/*
 * MAT1
 */
int nastranCard_mat1(FILE *fp, const int *mid, const double *e, const double* g,
                     const double *nu, const double *rho, const double *a, const double *tref,
                     const double *ge, const double *st, const double *sc, const double *ss,
                     const int *mcsid, feaFileTypeEnum formatType);

/*
 * MAT8
 */
int nastranCard_mat8(FILE *fp, const int *mid, const double *e1, const double *e2,
                     const double *nu12, const double *g12, const double *g1z,
                     const double *g2z, const double *rho, const double *a1, const double *a2,
                     const double *tref, const double *xt, const double *xc,
                     const double *yt, const double *yc, const double *s, const double *ge,
                     const double *f12, const double *strn,
                     feaFileTypeEnum formatType);

/*
 * MKAERO1
 */
int nastranCard_mkaero1(FILE *fp, int numM, double *m, int numK,
                        double *k, feaFileTypeEnum formatType);

/*
 * MOMENT
 */
int nastranCard_moment(FILE *fp, const int *sid, const int *g, const int *cid, const double *m,
                       const double n[3], feaFileTypeEnum formatType);

/*
 * PAERO1
 */
int nastranCard_paero1(FILE *fp, const int *pid, int numB, const int *b,
                       feaFileTypeEnum formatType);

/*
 * PBAR
 */
int nastranCard_pbar(FILE *fp, const int *pid, const int *mid, const double *a,
                     const double *i1, const double *i2, const double *i12, const double *j,
                     const double *nsm, const double c[2], const double d[2],
                     const double e[2], const double f[2], const double *k1, double *k2,
                     feaFileTypeEnum formatType);

/*
 * PBARL
 */
int nastranCard_pbarl(FILE *fp, const int *pid, const int *mid, const char *type,
                      const double *f0, const int numDim, const double *dim, const double *nsm,
                      feaFileTypeEnum formatType);

/*
 * PCOMP
 */
int nastranCard_pcomp(FILE *fp, const int *pid, const double *z0, const double* nsm,
                      const double *sb, const char *ft, const double *tref, const double *ge,
                      const char *lam, int numLayers, const int *mid, const double *t,
                      const double *theta, const char **sout,
                      feaFileTypeEnum formatType);

/*
 * PLOAD2
 */
int nastranCard_pload2(FILE *fp, const int *sid, const double *p, int numE,
                       const int *eid, feaFileTypeEnum formatType);

/*
 * PLOAD4
 */
int nastranCard_pload4(FILE *fp, const int *sid, const int *eid, const double p[4],
                       const int *g1, const int *g3, const int *cid, const double n[3],
                       feaFileTypeEnum formatType);

/*
 * PROD
 */
int nastranCard_prod(FILE *fp, const int *pid, const int *mid, const double *a, const double *j,
                     const double *c, const double *nsm,
                     feaFileTypeEnum formatType);

/*
 * PSHEAR
 */
int nastranCard_pshear(FILE *fp, const int *pid, const int *mid, const double *t,
                       const double *nsm, const double *f1, const double *f2,
                       feaFileTypeEnum formatType);

/*
 * PSHELL
 */
int nastranCard_pshell(FILE *fp, const int *pid, const int *mid1, const double *t,
                       const int *mid2, const double *i12t3, const int *mid3,
                       const double *tst, const double *nsm, const double *z1,
                       const double *z2, const int *mid4,
                       feaFileTypeEnum formatType);

/*
 * PSOLID
 */
int nastranCard_psolid(FILE *fp, const int *pid, const int *mid, const int *cordm,
                       const char *in, const char *stress, const char *isop, const char *fctn,
                       feaFileTypeEnum formatType);

/*
 * RBE2
 */
int nastranCard_rbe2(FILE *fp, const int *eid, const int *gn, const int *cm,
                     int numGM, const int *gm, feaFileTypeEnum formatType);

/*
 * RBE3
 */
int nastranCard_rbe3(FILE *fp, const int *eid, const int *refgrid, const int *refc,
                     int numG, const double *wt, const int *c, const int *g,
                     int numGM, const int *gm, const int *cm,
                     feaFileTypeEnum formatType);

/*
 * RFORCE
 */
int nastranCard_rforce(FILE *fp, const int *sid, const int *g, const int *cid, const double *a,
                       const double r[3], const int *method, const double *racc,
                       feaFileTypeEnum formatType);

/*
 * SET1
 */
int nastranCard_set1(FILE *fp, const int *sid, int numG, const int *g,
                     feaFileTypeEnum formatType);

/*
 * SPC
 */
int nastranCard_spc(FILE *fp, const int *sid, int numSPC, const int *g, const int *c,
                    const double *d, feaFileTypeEnum formatType);

/*
 * SPC1
 */
int nastranCard_spc1(FILE *fp, const int *sid, const int *c, int numSPC, const int *g,
                     feaFileTypeEnum formatType);

/*
 * SPCADD
 */
int nastranCard_spcadd(FILE *fp, int *sid, int numSPC, int *s,
                       feaFileTypeEnum formatType);

/*
 * SPLINE1
 */
int nastranCard_spline1(FILE *fp, const int *eid, const int *caero, const int *box1,
                        const int *box2, const int *setg, const double *dz,
                        feaFileTypeEnum formatType);

/*
 * SUPORT
 */
int nastranCard_suport(FILE *fp, int numID, const int *id, const int *c,
                       feaFileTypeEnum formatType);

/*
 * SUPORT1
 */
int nastranCard_suport1(FILE *fp, const int *sid, int numID, const int *id, const int *c,
                       feaFileTypeEnum formatType);

/*
 * TEMP
 */
int nastranCard_temp(FILE *fp, const int *sid, int numG, const int *g, const double *t,
                     feaFileTypeEnum formatType);

/*
 * TEMPD
 */
int nastranCard_tempd(FILE *fp, int numSID, const int *sid, const double *t,
                      feaFileTypeEnum formatType);

/*
 * TRIM
 */
int nastranCard_trim(FILE *fp, const int *id, const double *mach, const double *q,
                     int numVar, char **label, const double *ux,
                     feaFileTypeEnum formatType);


#ifdef __cplusplus
}
#endif

#endif // __NASTRAN_CARDS_H__
