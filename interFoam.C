#include "fvCFD.H"          //在OpenFOAM的求解器中，涉及到的构建时间、组建矩阵、有限体积离散、组建网格、量纲设置等有限体积库
#include "dynamicFvMesh.H"
#include "CMULES.H"         //求解相体积分数的模型
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
int main(int argc, char *argv[])
{
    #include "postProcess.H"
 
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"   //在这里建立mixture模型
            /*Info<< "Reading transportProperties\n" << endl;
            immiscibleIncompressibleTwoPhaseMixture mixture(U, phi); //调用incompressibleTwoPhaseMixture(U, phi),interfaceProperties(alpha1(), U, *this)构造函数
            volScalarField& alpha1(mixture.alpha1());
            volScalarField& alpha2(mixture.alpha2());
            const dimensionedScalar& rho1 = mixture.rho1();
            const dimensionedScalar& rho2 = mixture.rho2();*/


    #include "createAlphaFluxes.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"
 
    turbulence->validate();
 
    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }
 
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;
 
    while (runTime.run())
    {
        #include "readDyMControls.H"
 
        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }
 
        runTime++;
 
        Info<< "Time = " << runTime.timeName() << nl << endl;
 
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.update();
 
                if (mesh.changing())
                {
                    // Do not apply previous time-step mesh compression flux
                    // if the mesh topology changed
                    if (mesh.topoChanging())
                    {
                        talphaPhi1Corr0.clear();
                    }
 
                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;
 
                    MRF.update();
 
                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();
 
                        #include "correctPhi.H"
 
                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
 
                        mixture.correct();
                    }
 
                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }
 
            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"
 
            mixture.correct();
 
            #include "UEqn.H"
 
            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }
 
            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }
 
        runTime.write();
 
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
 
    Info<< "End\n" << endl;
 
    return 0;
}
 
 
// ************************************************************************* //
