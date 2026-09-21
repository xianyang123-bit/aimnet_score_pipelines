
MODULE OpenMM_Drude_Types
    implicit none

    ! Global Constants


    ! Type Declarations

    type OpenMM_DrudeNoseHooverIntegrator
        integer*8 :: handle = 0
    end type

    type OpenMM_DrudeForce
        integer*8 :: handle = 0
    end type

    type OpenMM_DrudeIntegrator
        integer*8 :: handle = 0
    end type

    type OpenMM_DrudeLangevinIntegrator
        integer*8 :: handle = 0
    end type

    type OpenMM_DrudeSCFIntegrator
        integer*8 :: handle = 0
    end type

    ! Enumerations

    integer*4, parameter :: OpenMM_Drude_False = 0
    integer*4, parameter :: OpenMM_Drude_True = 1

END MODULE OpenMM_Drude_Types

MODULE OpenMM_Drude
    use OpenMM
    use OpenMM_Drude_Types
    implicit none
    
    interface


    ! OpenMM::DrudeNoseHooverIntegrator
    subroutine OpenMM_DrudeNoseHooverIntegrator_create(result, temperature, &
collisionFrequency, &
drudeTemperature, &
drudeCollisionFrequency, &
stepSize, &
chainLength, &
numMTS, &
numYoshidaSuzuki)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeNoseHooverIntegrator) :: result
        real*8 :: temperature
        real*8 :: collisionFrequency
        real*8 :: drudeTemperature
        real*8 :: drudeCollisionFrequency
        real*8 :: stepSize
        integer*4 :: chainLength
        integer*4 :: numMTS
        integer*4 :: numYoshidaSuzuki
    end subroutine
    subroutine OpenMM_DrudeNoseHooverIntegrator_destroy(destroy)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeNoseHooverIntegrator) :: destroy
    end subroutine
    function OpenMM_DrudeNoseHooverIntegrator_getMaxDrudeDistance(target)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeNoseHooverIntegrator) :: target
        real*8 :: OpenMM_DrudeNoseHooverIntegrator_getMaxDrudeDistance
    end function
    subroutine OpenMM_DrudeNoseHooverIntegrator_setMaxDrudeDistance(target, distance)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeNoseHooverIntegrator) :: target
        real*8 :: distance
    end subroutine
    function OpenMM_DrudeNoseHooverIntegrator_computeDrudeKineticEnergy(target)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeNoseHooverIntegrator) :: target
        real*8 :: OpenMM_DrudeNoseHooverIntegrator_computeDrudeKineticEnergy
    end function
    function OpenMM_DrudeNoseHooverIntegrator_computeTotalKineticEnergy(target)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeNoseHooverIntegrator) :: target
        real*8 :: OpenMM_DrudeNoseHooverIntegrator_computeTotalKineticEnergy
    end function
    function OpenMM_DrudeNoseHooverIntegrator_computeSystemTemperature(target)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeNoseHooverIntegrator) :: target
        real*8 :: OpenMM_DrudeNoseHooverIntegrator_computeSystemTemperature
    end function
    function OpenMM_DrudeNoseHooverIntegrator_computeDrudeTemperature(target)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeNoseHooverIntegrator) :: target
        real*8 :: OpenMM_DrudeNoseHooverIntegrator_computeDrudeTemperature
    end function

    ! OpenMM::DrudeForce
    subroutine OpenMM_DrudeForce_create(result)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeForce) :: result
    end subroutine
    subroutine OpenMM_DrudeForce_destroy(destroy)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeForce) :: destroy
    end subroutine
    function OpenMM_DrudeForce_getNumParticles(target)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeForce) :: target
        integer*4 :: OpenMM_DrudeForce_getNumParticles
    end function
    function OpenMM_DrudeForce_getNumScreenedPairs(target)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeForce) :: target
        integer*4 :: OpenMM_DrudeForce_getNumScreenedPairs
    end function
    function OpenMM_DrudeForce_addParticle(target, particle, &
particle1, &
particle2, &
particle3, &
particle4, &
charge, &
polarizability, &
aniso12, &
aniso34)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeForce) :: target
        integer*4 :: particle
        integer*4 :: particle1
        integer*4 :: particle2
        integer*4 :: particle3
        integer*4 :: particle4
        real*8 :: charge
        real*8 :: polarizability
        real*8 :: aniso12
        real*8 :: aniso34
        integer*4 :: OpenMM_DrudeForce_addParticle
    end function
    subroutine OpenMM_DrudeForce_getParticleParameters(target, index, &
particle, &
particle1, &
particle2, &
particle3, &
particle4, &
charge, &
polarizability, &
aniso12, &
aniso34)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeForce) :: target
        integer*4 :: index
        integer*4 :: particle
        integer*4 :: particle1
        integer*4 :: particle2
        integer*4 :: particle3
        integer*4 :: particle4
        real*8 :: charge
        real*8 :: polarizability
        real*8 :: aniso12
        real*8 :: aniso34
    end subroutine
    subroutine OpenMM_DrudeForce_setParticleParameters(target, index, &
particle, &
particle1, &
particle2, &
particle3, &
particle4, &
charge, &
polarizability, &
aniso12, &
aniso34)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeForce) :: target
        integer*4 :: index
        integer*4 :: particle
        integer*4 :: particle1
        integer*4 :: particle2
        integer*4 :: particle3
        integer*4 :: particle4
        real*8 :: charge
        real*8 :: polarizability
        real*8 :: aniso12
        real*8 :: aniso34
    end subroutine
    function OpenMM_DrudeForce_addScreenedPair(target, particle1, &
particle2, &
thole)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeForce) :: target
        integer*4 :: particle1
        integer*4 :: particle2
        real*8 :: thole
        integer*4 :: OpenMM_DrudeForce_addScreenedPair
    end function
    subroutine OpenMM_DrudeForce_getScreenedPairParameters(target, index, &
particle1, &
particle2, &
thole)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeForce) :: target
        integer*4 :: index
        integer*4 :: particle1
        integer*4 :: particle2
        real*8 :: thole
    end subroutine
    subroutine OpenMM_DrudeForce_setScreenedPairParameters(target, index, &
particle1, &
particle2, &
thole)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeForce) :: target
        integer*4 :: index
        integer*4 :: particle1
        integer*4 :: particle2
        real*8 :: thole
    end subroutine
    subroutine OpenMM_DrudeForce_updateParametersInContext(target, context)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeForce) :: target
        type(OpenMM_Context) :: context
    end subroutine
    subroutine OpenMM_DrudeForce_setUsesPeriodicBoundaryConditions(target, periodic)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeForce) :: target
        integer*4 :: periodic
    end subroutine
    function OpenMM_DrudeForce_usesPeriodicBoundaryConditions(target)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeForce) :: target
        integer*4 :: OpenMM_DrudeForce_usesPeriodicBoundaryConditions
    end function

    ! OpenMM::DrudeIntegrator
    subroutine OpenMM_DrudeIntegrator_create(result, stepSize)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeIntegrator) :: result
        real*8 :: stepSize
    end subroutine
    subroutine OpenMM_DrudeIntegrator_destroy(destroy)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeIntegrator) :: destroy
    end subroutine
    subroutine OpenMM_DrudeIntegrator_step(target, steps)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeIntegrator) :: target
        integer*4 :: steps
    end subroutine
    function OpenMM_DrudeIntegrator_getDrudeTemperature(target)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeIntegrator) :: target
        real*8 :: OpenMM_DrudeIntegrator_getDrudeTemperature
    end function
    subroutine OpenMM_DrudeIntegrator_setDrudeTemperature(target, temp)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeIntegrator) :: target
        real*8 :: temp
    end subroutine
    function OpenMM_DrudeIntegrator_getMaxDrudeDistance(target)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeIntegrator) :: target
        real*8 :: OpenMM_DrudeIntegrator_getMaxDrudeDistance
    end function
    subroutine OpenMM_DrudeIntegrator_setMaxDrudeDistance(target, distance)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeIntegrator) :: target
        real*8 :: distance
    end subroutine
    subroutine OpenMM_DrudeIntegrator_setRandomNumberSeed(target, seed)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeIntegrator) :: target
        integer*4 :: seed
    end subroutine
    function OpenMM_DrudeIntegrator_getRandomNumberSeed(target)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeIntegrator) :: target
        integer*4 :: OpenMM_DrudeIntegrator_getRandomNumberSeed
    end function

    ! OpenMM::DrudeLangevinIntegrator
    subroutine OpenMM_DrudeLangevinIntegrator_create(result, temperature, &
frictionCoeff, &
drudeTemperature, &
drudeFrictionCoeff, &
stepSize)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeLangevinIntegrator) :: result
        real*8 :: temperature
        real*8 :: frictionCoeff
        real*8 :: drudeTemperature
        real*8 :: drudeFrictionCoeff
        real*8 :: stepSize
    end subroutine
    subroutine OpenMM_DrudeLangevinIntegrator_destroy(destroy)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeLangevinIntegrator) :: destroy
    end subroutine
    function OpenMM_DrudeLangevinIntegrator_getTemperature(target)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeLangevinIntegrator) :: target
        real*8 :: OpenMM_DrudeLangevinIntegrator_getTemperature
    end function
    subroutine OpenMM_DrudeLangevinIntegrator_setTemperature(target, temp)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeLangevinIntegrator) :: target
        real*8 :: temp
    end subroutine
    function OpenMM_DrudeLangevinIntegrator_getFriction(target)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeLangevinIntegrator) :: target
        real*8 :: OpenMM_DrudeLangevinIntegrator_getFriction
    end function
    subroutine OpenMM_DrudeLangevinIntegrator_setFriction(target, coeff)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeLangevinIntegrator) :: target
        real*8 :: coeff
    end subroutine
    function OpenMM_DrudeLangevinIntegrator_getDrudeFriction(target)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeLangevinIntegrator) :: target
        real*8 :: OpenMM_DrudeLangevinIntegrator_getDrudeFriction
    end function
    subroutine OpenMM_DrudeLangevinIntegrator_setDrudeFriction(target, coeff)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeLangevinIntegrator) :: target
        real*8 :: coeff
    end subroutine
    subroutine OpenMM_DrudeLangevinIntegrator_step(target, steps)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeLangevinIntegrator) :: target
        integer*4 :: steps
    end subroutine
    function OpenMM_DrudeLangevinIntegrator_computeSystemTemperature(target)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeLangevinIntegrator) :: target
        real*8 :: OpenMM_DrudeLangevinIntegrator_computeSystemTemperature
    end function
    function OpenMM_DrudeLangevinIntegrator_computeDrudeTemperature(target)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeLangevinIntegrator) :: target
        real*8 :: OpenMM_DrudeLangevinIntegrator_computeDrudeTemperature
    end function

    ! OpenMM::DrudeSCFIntegrator
    subroutine OpenMM_DrudeSCFIntegrator_create(result, stepSize)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeSCFIntegrator) :: result
        real*8 :: stepSize
    end subroutine
    subroutine OpenMM_DrudeSCFIntegrator_destroy(destroy)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeSCFIntegrator) :: destroy
    end subroutine
    function OpenMM_DrudeSCFIntegrator_getMinimizationErrorTolerance(target)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeSCFIntegrator) :: target
        real*8 :: OpenMM_DrudeSCFIntegrator_getMinimizationErrorTolerance
    end function
    subroutine OpenMM_DrudeSCFIntegrator_setMinimizationErrorTolerance(target, tol)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeSCFIntegrator) :: target
        real*8 :: tol
    end subroutine
    subroutine OpenMM_DrudeSCFIntegrator_step(target, steps)
        use OpenMM
        use OpenMM_Drude_Types
        implicit none
        type(OpenMM_DrudeSCFIntegrator) :: target
        integer*4 :: steps
    end subroutine


    end interface
END MODULE OpenMM_Drude
