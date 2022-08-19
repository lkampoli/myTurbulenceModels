/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "kOmegaSST_BOR_frozen_R.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{
// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //
//#include "evalGeno.H"
template<class BasicTurbulenceModel>
tmp<volScalarField>
kOmegaSST_BOR_frozen_R<BasicTurbulenceModel>::kOmegaSST_BOR_frozen_R::F1
(
    const volScalarField& CDkOmega
) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar(dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
kOmegaSST_BOR_frozen_R<BasicTurbulenceModel>::kOmegaSST_BOR_frozen_R::F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
kOmegaSST_BOR_frozen_R<BasicTurbulenceModel>::kOmegaSST_BOR_frozen_R::F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150*(this->mu()/this->rho_)/(omega_*sqr(y_)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
kOmegaSST_BOR_frozen_R<BasicTurbulenceModel>::kOmegaSST_BOR_frozen_R::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23.ref() *= F3();
    }

    return f23;
}


template<class BasicTurbulenceModel>
void kOmegaSST_BOR_frozen_R<BasicTurbulenceModel>::correctNut
(
    const volScalarField& S2,
    const volScalarField& F2
)
{
    this->nut_ = a1_*k_/max(a1_*omega_, b1_*F2*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kOmegaSST_BOR_frozen_R<BasicTurbulenceModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))), F23());
}


template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> kOmegaSST_BOR_frozen_R<BasicTurbulenceModel>::Pk
(
    const volScalarField::Internal& G
) const
{
    return min(G, (c1_*betaStar_)*this->k_()*this->omega_());
}


template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> kOmegaSST_BOR_frozen_R<BasicTurbulenceModel>::epsilonByk
(
    const volScalarField::Internal& F1,
    const volScalarField::Internal& F2
) const
{
    return betaStar_*omega_();
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kOmegaSST_BOR_frozen_R<BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kOmegaSST_BOR_frozen_R<BasicTurbulenceModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kOmegaSST_BOR_frozen_R<BasicTurbulenceModel>::Qsas
(
    const volScalarField::Internal& S2,
    const volScalarField::Internal& gamma,
    const volScalarField::Internal& beta
) const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaSST_BOR_frozen_R<BasicTurbulenceModel>::kOmegaSST_BOR_frozen_R
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "b1",
            this->coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0
        )
    ),
    F3_
    (
        Switch::lookupOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),
	
	y_(wallDist::New(this->mesh_).y()),
	
    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

////////////////////////////////////////////////////////////////////////
//  added by Yuan Fang --start
////////////////////////////////////////////////////////////////////////    

     nonlinearStress_
    (
        IOobject
        (
            IOobject::groupName("nonlinearStress", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedSymmTensor
        (
            "nonlinearStress",
            sqr(dimVelocity),
            Zero
        )
    ),
     Rij_
    (
        IOobject
        (
            IOobject::groupName("Rij", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedSymmTensor
        (
            "Rij",
            sqr(dimVelocity),
            Zero
        )
    ),
     Rall_
    (
        IOobject
        (
            IOobject::groupName("Rall", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor
        (
            "Rall",
            sqr(dimVelocity),
            Zero
        )
    ),
    xswtch_
    (
        IOobject
        (
            IOobject::groupName("xswtch", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )

////////////////////////////////////////////////////////////////////////
//  added by Yuan Fang --end
////////////////////////////////////////////////////////////////////////    

{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    //ccl: Read the model equation from transportProperties
    //Info<< "Reading tensor expression from GEP \n" << endl;
    //Info<< "Tensor expression of nonlinearStress being run: " <<endl;
    //wordList zeta1
    //(
    //  this->coeffDict_.lookup("zeta1")
    //);
    //wordList zeta2
    //(
    //  this->coeffDict_.lookup("zeta2")
    //);
    //wordList zeta3
    //(
    //  this->coeffDict_.lookup("zeta3")
    //);
    //wordList zeta4
    //(
    //  this->coeffDict_.lookup("zeta4")
    //);
    //wordList zeta5
    //(
    //  this->coeffDict_.lookup("zeta5")
    //);
    //wordList zeta6
    //(
    //  this->coeffDict_.lookup("zeta6")
    //);
	
    //  label i = 0, j = 0, k = 0, e = 0, f=0, g = 0;
    //  Info << "(";
    //  printGeno(zeta1, i);  
    //  Info << ") * T1 + (";
    //  printGeno(zeta2, j); 
    //  Info << ") * T2 + (";
    //  printGeno(zeta3, k);
    //  Info<< ") * T3 \n"<< endl;

    //  Info<< "Tensor expression of Rij being run: " <<endl;
    //  Info << "(";
    //  printGeno(zeta4, e);  
    //  Info << ") * T1 + (";
    //  printGeno(zeta5, f); 
    //  Info << ") * T2 + (";
    //  printGeno(zeta6, g);
    //  Info<< ") * T3 \n"<< endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

////////////////////////////////////////////////////////////////////////
//  added by Yuan Fang --start
////////////////////////////////////////////////////////////////////////

template<class BasicTurbulenceModel>
Foam::tmp<Foam::volSymmTensorField>
kOmegaSST_BOR_frozen_R<BasicTurbulenceModel>::R() const
{
    tmp<volSymmTensorField> tR
    (
        eddyViscosity<RASModel<BasicTurbulenceModel>>::R()
    );
    tR.ref() += nonlinearStress_;
    return tR;
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volSymmTensorField>
kOmegaSST_BOR_frozen_R<BasicTurbulenceModel>::devRhoReff() const
{
    tmp<volSymmTensorField> tdevRhoReff
    (
        eddyViscosity<RASModel<BasicTurbulenceModel>>::devRhoReff()
    );
    tdevRhoReff.ref() += this->rho_*nonlinearStress_;
    return tdevRhoReff;
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::fvVectorMatrix>
kOmegaSST_BOR_frozen_R<BasicTurbulenceModel>::divDevRhoReff
(
    volVectorField& U
) const
{
    return
    (
        fvc::div(this->rho_*nonlinearStress_) +
        eddyViscosity<RASModel<BasicTurbulenceModel>>::divDevRhoReff(U)
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::fvVectorMatrix>
kOmegaSST_BOR_frozen_R<BasicTurbulenceModel>::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    return
    (
        fvc::div(rho*nonlinearStress_) +
        eddyViscosity<RASModel<BasicTurbulenceModel>>::divDevRhoReff(rho, U)
    );
}

////////////////////////////////////////////////////////////////////////
//  added by Yuan Fang --end
////////////////////////////////////////////////////////////////////////

template<class BasicTurbulenceModel>
bool kOmegaSST_BOR_frozen_R<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void kOmegaSST_BOR_frozen_R<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();
	
    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))()()
    );
    dimensionedScalar Nutsmall(
        "0", dimensionSet(0,2,-1,0,0,0,0), 1e-20
    );
    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));
    //volScalarField S = sqrt(2*magSqr(symm(tgradU())));
    volScalarField::Internal GbyNu((dev(twoSymm(tgradU()()))) && tgradU()());
    volScalarField::Internal GbyNuaijx(((dev(twoSymm(tgradU()()))) - (nonlinearStress_/nut())) && tgradU()());
    volScalarField::Internal G(this->GName(), nut()*GbyNuaijx);
    
////////////////////////////////////////////////////////////////////////
//  added by Yuan Fang --start
////////////////////////////////////////////////////////////////////////

    volScalarField S = sqrt(S2);
	
    //volScalarField tau = 1./max( S/a1_ + omegaMin_,omega_ + omegaMin_);
    volScalarField tau =  1./max( S/0.31 + this->omegaMin_,(omega_ + this->omegaMin_));
    volScalarField tau2 = sqr(tau);
	
    volSymmTensorField sij(dev(symm(tgradU()))); 
    volTensorField omegaij((skew(tgradU())));
	
    dimensionedScalar osmall ("0",dimensionSet(0,0,-1,0,0,0,0),0.01);
    dimensionedScalar ksmall ("0",dimensionSet(0,2,-2,0,0,0,0),0.001);

    volSymmTensorField Sijt = sij * tau;
    volTensorField     Oijt = omegaij * tau;

    volScalarField I1 = tr(Sijt & Sijt);
    volScalarField I2 = tr(Oijt & Oijt);
    volScalarField I3 = tr(Sijt & (Sijt & Sijt));
    volScalarField I4 = tr(Oijt & (Oijt & Sijt));
    volScalarField I5 = tr(Oijt & (Oijt & (Sijt & Sijt)));

    volSymmTensorField T1 = Sijt;
    volSymmTensorField T2 = symm((Sijt & Oijt) - (Oijt & Sijt));
    volSymmTensorField T3 = symm( Sijt & Sijt) - scalar(1.0/3.0)*I*I1;
    volSymmTensorField T4 = symm( Oijt & Oijt) - scalar(1.0/3.0)*I*I2;
    volSymmTensorField T5 = symm((Oijt & (Sijt & Sijt)) - ((Sijt & Sijt) & Oijt));
    volSymmTensorField T6 = symm(((Oijt & Oijt) & Sijt) + ((Sijt & Sijt) & Oijt) - (2.0/3.0) * I * tr(Sijt & (Oijt & Oijt)));
    volSymmTensorField T7 = symm((Oijt & (Sijt & (Oijt & Oijt))) - (Oijt & (Oijt & (Sijt & Oijt))));
    volSymmTensorField T8 = symm((Sijt & (Oijt & (Sijt & Sijt))) - (Sijt & (Sijt & (Oijt & Sijt))));
    volSymmTensorField T9 = symm((Oijt & (Oijt & (Sijt & Sijt))) + (Sijt & (Sijt & (Oijt & Oijt))) - (2.0/3.0) * I * tr(Sijt & (Sijt & (Oijt & Oijt))));
    volSymmTensorField T10= symm((Oijt & (Sijt & (Sijt & (Oijt & Oijt)))) -(Oijt & (Oijt & (Sijt & (Sijt & Oijt)))));

    //volScalarField I1(tau2*(sij && sij));
    //volScalarField I2(tau2*(omegaij && omegaij));
	
    //volSymmTensorField T1(tau*sij);
    //volSymmTensorField T2(tau2*(symm((sij & omegaij) - (omegaij & sij))));
    //volSymmTensorField T3(tau2*(symm(sij & sij)) - (scalar(1.0/3.0))*I1*I);   	

    // imported from kOmegaSSTGEP_Rijhat_added.C
    word model (
      this->coeffDict_.lookup ("model")
    );

    word aijx_term (
      this->coeffDict_.lookup("aijx_term")
      //this->coeffDict_.lookupOrDefault<word> (
      //      "aijx_term",
      //      "1.0"
      //  )
    );

    word Rterm_corr (
      this->coeffDict_.lookup("Rterm_corr")
      //this->coeffDict_.lookupOrDefault<word> (
      //      "Rterm_corr",
      //      "1.0"
      //  )
    );

   //ccl: Read the model equation from transportProperties
   // wordList zeta1
   // (
   //   this->coeffDict_.lookup("zeta1")
   // );
   // wordList zeta2
   // (
   //   this->coeffDict_.lookup("zeta2")
   // );
   // wordList zeta3
   // (
   //   this->coeffDict_.lookup("zeta3")
   // );
   // wordList zeta4
   // (
   //   this->coeffDict_.lookup("zeta4")
   // );
   // wordList zeta5
   // (
   //   this->coeffDict_.lookup("zeta5")
   // );
   // wordList zeta6
   // (
   //   this->coeffDict_.lookup("zeta6")
   // );
 
   // Here should be modified:
   // nonlinearStress_ and Rij_ should not be interpreted by evalGeno anymore
   // but read from some other file or directly harcoded from pre-trained 
   // offline model ...
   // label i = 0, j = 0, k = 0, e = 0, f = 0, g = 0;
   // nonlinearStress_=2*k_*((evalGeno(zeta1, i, I1, I2))*T1 + (evalGeno(zeta2, j, I1, I2))*T2 + (evalGeno(zeta3, k, I1, I2))*T3);
   // Rij_=2*k_*((evalGeno(zeta4, e, I1, I2))*T1 + (evalGeno(zeta5, f, I1, I2))*T2 + (evalGeno(zeta6, g, I1, I2))*T3);

    // library of GEP models
    if (model == "linear") {
	nonlinearStress_ = T1 - T1;
    }
    #include "nonLinearModel.H"
    else {
	FatalError << "unknown model " << model << nl << exit(FatalError);
    }

    nonlinearStress_ *= 2*k_;
    nonlinearStress_.correctBoundaryConditions();

    Rij_ *= 2*k_;
    Rij_.correctBoundaryConditions();

    Rall_= ((2.0/3.0)*I)*k_ - this->nut_*dev(twoSymm(tgradU())) + nonlinearStress_;
    Rall_.correctBoundaryConditions();
    
    volScalarField Rterm(Rij_ && symm(tgradU()));	
    tgradU.clear();

    dimensionedScalar small_val
    (
       "small_val",
       dimensionSet(0, 2, -1, 0, 0, 0, 0),
       1e-25
    );

////////////////////////////////////////////////////////////////////////
//  added by Yuan Fang --end
////////////////////////////////////////////////////////////////////////
	
    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    // Access to the density
    //volScalarField& rho = const_cast<volScalarField&>
    //(
    //    this->mesh_.objectRegistry::template
    //    lookupObject<volScalarField>("rho")
    //);
    //Info << " rho = " << rho << endl;

    volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F23(this->F23());

    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
          + fvm::div(alphaRhoPhi, omega_)
          - fvm::laplacian(alpha*rho*DomegaEff(F1), omega_)
         ==
            alpha()*rho()*gamma
           *min
            (
                GbyNuaijx,
                (c1_/a1_)*betaStar_*omega_()
               *max(a1_*omega_(), b1_*F23()*sqrt(S2()))
            )
          + alpha()*rho()*gamma*Rterm/(nut()+small_val)
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, omega_)
          - fvm::Sp(alpha()*rho()*beta*omega_(), omega_)
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/omega_(),
                omega_
            )
          + Qsas(S2(), gamma, beta)
          + omegaSource()
          + fvOptions(alpha, rho, omega_)
        );

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvOptions.correct(omega_);
        bound(omega_, this->omegaMin_);
    }

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(F1), k_)
     ==
        alpha()*rho()*Pk(G)
//      + mag(alpha()*rho()*Rterm)
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*epsilonByk(F1, F23), k_)
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut(S2, F23);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace RASModels
// ************************************************************************* //
