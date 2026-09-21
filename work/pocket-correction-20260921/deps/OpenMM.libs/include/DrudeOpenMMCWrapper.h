
#ifndef DRUDE_OPENMM_CWRAPPER_H_
#define DRUDE_OPENMM_CWRAPPER_H_

#ifndef OPENMM_EXPORT_DRUDE
#define OPENMM_EXPORT_DRUDE
#endif
/* Global Constants */


/* Type Declarations */

typedef struct OpenMM_DrudeNoseHooverIntegrator_struct OpenMM_DrudeNoseHooverIntegrator;
typedef struct OpenMM_DrudeForce_struct OpenMM_DrudeForce;
typedef struct OpenMM_DrudeIntegrator_struct OpenMM_DrudeIntegrator;
typedef struct OpenMM_DrudeLangevinIntegrator_struct OpenMM_DrudeLangevinIntegrator;
typedef struct OpenMM_DrudeSCFIntegrator_struct OpenMM_DrudeSCFIntegrator;

typedef struct OpenMM_2D_IntArray_struct OpenMM_2D_IntArray;
typedef struct OpenMM_3D_DoubleArray_struct OpenMM_3D_DoubleArray;

#if defined(__cplusplus)
extern "C" {
#endif

/* OpenMM_3D_DoubleArray */
OPENMM_EXPORT_DRUDE OpenMM_3D_DoubleArray* OpenMM_3D_DoubleArray_create(int size1, int size2, int size3);
OPENMM_EXPORT_DRUDE void OpenMM_3D_DoubleArray_set(OpenMM_3D_DoubleArray* array, int index1, int index2, OpenMM_DoubleArray* values);
OPENMM_EXPORT_DRUDE void OpenMM_3D_DoubleArray_destroy(OpenMM_3D_DoubleArray* array);

/* DrudeNoseHooverIntegrator */
extern OPENMM_EXPORT_DRUDE OpenMM_DrudeNoseHooverIntegrator* OpenMM_DrudeNoseHooverIntegrator_create(double temperature, double collisionFrequency, double drudeTemperature, double drudeCollisionFrequency, double stepSize, int chainLength, int numMTS, int numYoshidaSuzuki);
extern OPENMM_EXPORT_DRUDE void OpenMM_DrudeNoseHooverIntegrator_destroy(OpenMM_DrudeNoseHooverIntegrator* target);
extern OPENMM_EXPORT_DRUDE double OpenMM_DrudeNoseHooverIntegrator_getMaxDrudeDistance(const OpenMM_DrudeNoseHooverIntegrator* target);
extern OPENMM_EXPORT_DRUDE void OpenMM_DrudeNoseHooverIntegrator_setMaxDrudeDistance(OpenMM_DrudeNoseHooverIntegrator* target, double distance);
extern OPENMM_EXPORT_DRUDE double OpenMM_DrudeNoseHooverIntegrator_computeDrudeKineticEnergy(OpenMM_DrudeNoseHooverIntegrator* target);
extern OPENMM_EXPORT_DRUDE double OpenMM_DrudeNoseHooverIntegrator_computeTotalKineticEnergy(OpenMM_DrudeNoseHooverIntegrator* target);
extern OPENMM_EXPORT_DRUDE double OpenMM_DrudeNoseHooverIntegrator_computeSystemTemperature(OpenMM_DrudeNoseHooverIntegrator* target);
extern OPENMM_EXPORT_DRUDE double OpenMM_DrudeNoseHooverIntegrator_computeDrudeTemperature(OpenMM_DrudeNoseHooverIntegrator* target);

/* DrudeForce */
extern OPENMM_EXPORT_DRUDE OpenMM_DrudeForce* OpenMM_DrudeForce_create();
extern OPENMM_EXPORT_DRUDE void OpenMM_DrudeForce_destroy(OpenMM_DrudeForce* target);
extern OPENMM_EXPORT_DRUDE int OpenMM_DrudeForce_getNumParticles(const OpenMM_DrudeForce* target);
extern OPENMM_EXPORT_DRUDE int OpenMM_DrudeForce_getNumScreenedPairs(const OpenMM_DrudeForce* target);
extern OPENMM_EXPORT_DRUDE int OpenMM_DrudeForce_addParticle(OpenMM_DrudeForce* target, int particle, int particle1, int particle2, int particle3, int particle4, double charge, double polarizability, double aniso12, double aniso34);
extern OPENMM_EXPORT_DRUDE void OpenMM_DrudeForce_getParticleParameters(const OpenMM_DrudeForce* target, int index, int* particle, int* particle1, int* particle2, int* particle3, int* particle4, double* charge, double* polarizability, double* aniso12, double* aniso34);
extern OPENMM_EXPORT_DRUDE void OpenMM_DrudeForce_setParticleParameters(OpenMM_DrudeForce* target, int index, int particle, int particle1, int particle2, int particle3, int particle4, double charge, double polarizability, double aniso12, double aniso34);
extern OPENMM_EXPORT_DRUDE int OpenMM_DrudeForce_addScreenedPair(OpenMM_DrudeForce* target, int particle1, int particle2, double thole);
extern OPENMM_EXPORT_DRUDE void OpenMM_DrudeForce_getScreenedPairParameters(const OpenMM_DrudeForce* target, int index, int* particle1, int* particle2, double* thole);
extern OPENMM_EXPORT_DRUDE void OpenMM_DrudeForce_setScreenedPairParameters(OpenMM_DrudeForce* target, int index, int particle1, int particle2, double thole);
extern OPENMM_EXPORT_DRUDE void OpenMM_DrudeForce_updateParametersInContext(OpenMM_DrudeForce* target, OpenMM_Context* context);
extern OPENMM_EXPORT_DRUDE void OpenMM_DrudeForce_setUsesPeriodicBoundaryConditions(OpenMM_DrudeForce* target, OpenMM_Boolean periodic);
extern OPENMM_EXPORT_DRUDE OpenMM_Boolean OpenMM_DrudeForce_usesPeriodicBoundaryConditions(const OpenMM_DrudeForce* target);

/* DrudeIntegrator */
extern OPENMM_EXPORT_DRUDE OpenMM_DrudeIntegrator* OpenMM_DrudeIntegrator_create(double stepSize);
extern OPENMM_EXPORT_DRUDE void OpenMM_DrudeIntegrator_destroy(OpenMM_DrudeIntegrator* target);
extern OPENMM_EXPORT_DRUDE void OpenMM_DrudeIntegrator_step(OpenMM_DrudeIntegrator* target, int steps);
extern OPENMM_EXPORT_DRUDE double OpenMM_DrudeIntegrator_getDrudeTemperature(const OpenMM_DrudeIntegrator* target);
extern OPENMM_EXPORT_DRUDE void OpenMM_DrudeIntegrator_setDrudeTemperature(OpenMM_DrudeIntegrator* target, double temp);
extern OPENMM_EXPORT_DRUDE double OpenMM_DrudeIntegrator_getMaxDrudeDistance(const OpenMM_DrudeIntegrator* target);
extern OPENMM_EXPORT_DRUDE void OpenMM_DrudeIntegrator_setMaxDrudeDistance(OpenMM_DrudeIntegrator* target, double distance);
extern OPENMM_EXPORT_DRUDE void OpenMM_DrudeIntegrator_setRandomNumberSeed(OpenMM_DrudeIntegrator* target, int seed);
extern OPENMM_EXPORT_DRUDE int OpenMM_DrudeIntegrator_getRandomNumberSeed(const OpenMM_DrudeIntegrator* target);

/* DrudeLangevinIntegrator */
extern OPENMM_EXPORT_DRUDE OpenMM_DrudeLangevinIntegrator* OpenMM_DrudeLangevinIntegrator_create(double temperature, double frictionCoeff, double drudeTemperature, double drudeFrictionCoeff, double stepSize);
extern OPENMM_EXPORT_DRUDE void OpenMM_DrudeLangevinIntegrator_destroy(OpenMM_DrudeLangevinIntegrator* target);
extern OPENMM_EXPORT_DRUDE double OpenMM_DrudeLangevinIntegrator_getTemperature(const OpenMM_DrudeLangevinIntegrator* target);
extern OPENMM_EXPORT_DRUDE void OpenMM_DrudeLangevinIntegrator_setTemperature(OpenMM_DrudeLangevinIntegrator* target, double temp);
extern OPENMM_EXPORT_DRUDE double OpenMM_DrudeLangevinIntegrator_getFriction(const OpenMM_DrudeLangevinIntegrator* target);
extern OPENMM_EXPORT_DRUDE void OpenMM_DrudeLangevinIntegrator_setFriction(OpenMM_DrudeLangevinIntegrator* target, double coeff);
extern OPENMM_EXPORT_DRUDE double OpenMM_DrudeLangevinIntegrator_getDrudeFriction(const OpenMM_DrudeLangevinIntegrator* target);
extern OPENMM_EXPORT_DRUDE void OpenMM_DrudeLangevinIntegrator_setDrudeFriction(OpenMM_DrudeLangevinIntegrator* target, double coeff);
extern OPENMM_EXPORT_DRUDE void OpenMM_DrudeLangevinIntegrator_step(OpenMM_DrudeLangevinIntegrator* target, int steps);
extern OPENMM_EXPORT_DRUDE double OpenMM_DrudeLangevinIntegrator_computeSystemTemperature(OpenMM_DrudeLangevinIntegrator* target);
extern OPENMM_EXPORT_DRUDE double OpenMM_DrudeLangevinIntegrator_computeDrudeTemperature(OpenMM_DrudeLangevinIntegrator* target);

/* DrudeSCFIntegrator */
extern OPENMM_EXPORT_DRUDE OpenMM_DrudeSCFIntegrator* OpenMM_DrudeSCFIntegrator_create(double stepSize);
extern OPENMM_EXPORT_DRUDE void OpenMM_DrudeSCFIntegrator_destroy(OpenMM_DrudeSCFIntegrator* target);
extern OPENMM_EXPORT_DRUDE double OpenMM_DrudeSCFIntegrator_getMinimizationErrorTolerance(const OpenMM_DrudeSCFIntegrator* target);
extern OPENMM_EXPORT_DRUDE void OpenMM_DrudeSCFIntegrator_setMinimizationErrorTolerance(OpenMM_DrudeSCFIntegrator* target, double tol);
extern OPENMM_EXPORT_DRUDE void OpenMM_DrudeSCFIntegrator_step(OpenMM_DrudeSCFIntegrator* target, int steps);


#if defined(__cplusplus)
}
#endif

#endif /*DRUDE_OPENMM_CWRAPPER_H_*/
