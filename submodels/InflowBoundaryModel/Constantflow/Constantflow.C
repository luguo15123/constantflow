/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "Constantflow.H"
#include "constants.H"
#include "triPointRef.H"
#include "tetIndices.H"
#include "regExp.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::Constantflow<CloudType>::Constantflow
(
    const dictionary& dict,
    CloudType& cloud
)
:
    InflowBoundaryModel<CloudType>(dict, cloud, typeName),
    moleculeTypeIds_(),
    OutgassRate_(),
    outgassFluxAccumulators_(),
    constantFluxAccumulators_(),
    ConstantFlowNo_(),
    OutFlowNo_(),
    Umacro_
    (
      IOobject
      (
          "UMean",
          cloud.mesh().time().timeName(),
          cloud.mesh(),
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
      ),
      cloud.mesh(),
      dimensionedVector
      (
           "UMean",
           dimensionSet(0, 1, -1, 0, 0),
           vector::zero
      )
    ),
    MrhoN_
    (
      IOobject
      (
          "rhoNMean",
          cloud.mesh().time().timeName(),
          cloud.mesh(),
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
      ),
      cloud.mesh(),
      dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), VSMALL)
    ),
    CfNoDens_(),
    PressureFluxAccumulators_(),
    PressureNoDens_()
{
    labelList patchesog;
    labelList patchescf;
    labelList patchesfs;
    labelList patchesps;
    labelList patchespo;
    labelList patchesrf;

    // Temporary storage for the patch types
    DynamicList<label> dpatchesog;
    DynamicList<label> dpatchescf;
    DynamicList<label> dpatchesfs;
    DynamicList<label> dpatchesps;
    DynamicList<label> dpatchespo;
    DynamicList<label> dpatchesrf;

    regExp rog("^outgass");
    regExp rcf("^constant");
    regExp rfs("^freestream");
    regExp rps("^pressure");
    regExp rpo("^pump");
    regExp rrf("^reflect");

    // decide which boundaryMesh is pressure inlet patch
    // which boundary patch is constant flow patch
    // which is the freestream patch
    // which patch is the outgassing patch
    // which patch is the pumping out patch
    forAll(cloud.mesh().boundaryMesh(), p)
    {
        const polyPatch& patch = cloud.mesh().boundaryMesh()[p];

        if (isType<polyPatch>(patch))
        {
            if(rog.search(patch.name()))
                dpatchesog.append(p);

            else if(rcf.search(patch.name()))
                dpatchescf.append(p);

            else if(rfs.search(patch.name()))
		dpatchesfs.append(p);

            else if(rps.search(patch.name()))
		dpatchesps.append(p);
        }
        else if(isType<wallPolyPatch>(patch))
        {
            if(rpo.search(patch.name()))
                dpatchespo.append(p);

            else if(rrf.search(patch.name()))
		dpatchesrf.append(p);
        }
    }

    // transfer the temporary stored patch to the public variables
    patchesog.transfer(dpatchesog);
    patchescf.transfer(dpatchescf);
    patchesfs.transfer(dpatchesfs);
    patchesps.transfer(dpatchesps);
    patchespo.transfer(dpatchespo);
    patchesrf.transfer(dpatchesrf);

    List<word> molecules(cloud.typeIdList());

    moleculeTypeIds_.setSize(molecules.size());

    if(!patchesog.empty())
    {
        // Initialise the outgassFluxAccumulators_
        outgassFluxAccumulators_.setSize(patchesog.size());

        forAll(patchesog, p)
        {
            const polyPatch& patch = cloud.mesh().boundaryMesh()[patchesog[p]];

            outgassFluxAccumulators_[p] = List<Field<scalar> >
            (
                molecules.size(),
                Field<scalar>(patch.size(), 0.0)
            );
        }

        OutgassRate_ = List<List<scalar> >
        (
            patchesog.size(),
            List<scalar>(molecules.size(), 0.0)
        );

        forAll(patchesog, p)
        {
            const polyPatch& patch = cloud.mesh().boundaryMesh()[patchesog[p]];

            const dictionary& patchogdict
            (
                this->dict().subDict("OutgassingRate").subDict(patch.name())
            );

            forAll(molecules, j)
            {
                OutgassRate_[p][j] = readScalar
                (
                    patchogdict.lookup(molecules[j])
                );

                // the equivilent outgass particle number
                OutgassRate_[p][j] /= cloud.nParticle();

                moleculeTypeIds_[j] = findIndex(cloud.typeIdList(), molecules[j]);

                if (moleculeTypeIds_[j] == -1)
                {
                    FatalErrorIn
                    (
                        "Foam::Constantflow<CloudType>::Constantflow"
                        "("
                            "const dictionary&, "
                            "CloudType&"
                        ")"
                    ) << "typeId " << molecules[j] << "not defined in cloud." << nl
                      << abort(FatalError);
                }
            }
        }
    }

    if(!patchescf.empty())
    {

      // molecules mass for constant flow patch: start
      const dictionary& moleculeProperties
      (
          dict.subDict("moleculeProperties")
      );
      //average mass for the constant flow patch:end

      // Initialise the constantFluxAccumulators_
      constantFluxAccumulators_.setSize(patchescf.size());

      MacroVel_.setSize(patchescf.size());

      CfrhoN_.setSize(patchescf.size());

      cfmass_ = List<scalar>(patchescf.size(), 0.0);

      //CfrhoM_.setSize(patchescf.size());

      //CfMom_.setSize(patchescf.size());

      CfNoDens_.setSize(patchescf.size());

      // OutFlowNo_ indicate the number of particles that flow through the constant flow patches
      // Initialising the OutFlowNo_ to zero
      OutFlowNo_ = List<List<scalar> >
      (
          patchescf.size(),
          List<scalar>(molecules.size(), 0.0)
      );

      ConstantFlowNo_ = List<List<scalar> >
      (
          patchescf.size(),
          List<scalar>(molecules.size(), 0.0)
      );

      forAll(patchescf, p)
      {
        const polyPatch& patch = cloud.mesh().boundaryMesh()[patchescf[p]];

        constantFluxAccumulators_[p] = List<Field<scalar> >
        (
            molecules.size(),
            Field<scalar>(patch.size(), 0.0)
        );

        MacroVel_[p] = Field<vector>
        (
          patch.size(),
          vector::zero
        );

        CfrhoN_[p] = Field<scalar>
        (
          patch.size(),
          0.0
        );
		
/*        CfrhoM_[p] = Field<scalar>
        (
          patch.size(),
          0.0
        );
		
	CfMom_[p] = Field<vector>
        (
          patch.size(),
          vector::zero
        );
*/
        CfNoDens_[p] = List<Field<scalar> >
        (
            molecules.size(),
            Field<scalar>(patch.size(), 0.0)
        );

        const dictionary& patchcfrdict
        (
            this->dict().subDict("ConstantflowCoeffs").subDict(patch.name())
        );

        forAll(molecules, j)
        {
            //Info<< molecules[j] << nl <<endl;

            ConstantFlowNo_[p][j] = readScalar
            (
                patchcfrdict.lookup(molecules[j])
            );

            // the equivilent constant flow particle number
            ConstantFlowNo_[p][j] /= cloud.nParticle();

            moleculeTypeIds_[j] = findIndex(cloud.typeIdList(), molecules[j]);

            if (moleculeTypeIds_[j] == -1)
            {
                FatalErrorIn
                (
                    "Foam::Constantflow<CloudType>::Constantflow"
                    "("
                        "const dictionary&, "
                        "CloudType&"
                    ")"
                ) << "typeId " << molecules[j] << "not defined in cloud." << nl
                  << abort(FatalError);
            }
        }

        scalar tocfno = sum(ConstantFlowNo_[p]);

        forAll(molecules, j)
        {
           const dictionary& massdict
           (
              moleculeProperties.subDict(molecules[j])
           );

           scalar mass = readScalar
           (
                massdict.lookup("mass")
           );

           cfmass_[p] += ConstantFlowNo_[p][j]/tocfno*mass;
        }

      }
    }

    if(!patchesps.empty())
    {
      // Initialise the PressureFluxAccumulators_
      PressureFluxAccumulators_.setSize(patchesps.size());

      PressureVol_.setSize(patchesps.size());

      Psmacrovel_.setSize(patchesps.size());

      PressureNoDens_ = List<List<scalar> >
      (
          patchesps.size(),
          List<scalar>(molecules.size(), 0.0)
      );

      forAll(patchesps, p)
      {
        const polyPatch& patch = cloud.mesh().boundaryMesh()[patchesps[p]];

        PressureFluxAccumulators_[p] = List<Field<scalar> >
        (
            molecules.size(),
            Field<scalar>(patch.size(), 0.0)
        );


        Psmacrovel_[p] = Field<vector>(patch.size(), vector::zero);

        PressureVol_[p] = Patchvolume(patch, cloud.mesh());

        const dictionary& patchpnddict
        (
            this->dict().subDict("PressureNoDensity").subDict(patch.name())
        );

        const volScalarField::GeometricBoundaryField& boundaryT
        (
            cloud.boundaryT().boundaryField()
        );

        scalar pstemperature = gAverage(boundaryT[patchesps[p]]);

        forAll(molecules, j)
        {
            PressureNoDens_[p][j] = readScalar
            (
                patchpnddict.lookup(molecules[j])
            );

            Info<<" Pressure for patch " << p << " and molecule "<< molecules[j]
                <<" is " << PressureNoDens_[p][j] << endl;

            // Calculate the corresponding number density 
            if(pstemperature > SMALL)
            {
                PressureNoDens_[p][j] = PressureNoDens_[p][j]
                                    /Foam::constant::physicoChemical::k.value()
                                    /pstemperature;
            }


            // the equivilent constant flow particle number
            PressureNoDens_[p][j] /= cloud.nParticle();

            moleculeTypeIds_[j] = findIndex(cloud.typeIdList(), molecules[j]);

            if (moleculeTypeIds_[j] == -1)
            {
                FatalErrorIn
                (
                    "Foam::Constantflow<CloudType>::Constantflow"
                    "("
                        "const dictionary&, "
                        "CloudType&"
                    ")"
                ) << "typeId " << molecules[j] << "not defined in cloud." << nl
                  << abort(FatalError);
            }
        }
      }

      Pout << " The parameter for the pressure patches " << PressureNoDens_
           << endl;
    }

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::Constantflow<CloudType>::~Constantflow()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::Constantflow<CloudType>::autoMap(const mapPolyMesh& mapper)
{
    CloudType& cloud(this->owner());

    const polyMesh& mesh(cloud.mesh());

    if(cloud.patchescf_.size())
    {
      forAll(cloud.patchescf_, p)
      {
        label patchi = cloud.patchescf_[p];

        const polyPatch& patch = mesh.boundaryMesh()[patchi];
        List<Field<scalar> >& pFA = constantFluxAccumulators_[p];

        forAll(pFA, facei)
        {
            pFA[facei].setSize(patch.size(), 0);
        }
      }
    }

    if(cloud.patchesog_.size())
    {
      forAll(cloud.patchesog_, p)
      {
        label patchi = cloud.patchesog_[p];

        const polyPatch& patch = mesh.boundaryMesh()[patchi];
        List<Field<scalar> >& pFA = outgassFluxAccumulators_[p];

        forAll(pFA, facei)
        {
            pFA[facei].setSize(patch.size(), 0);
        }
      }
    }

    if(cloud.patchesps_.size())
    {
      forAll(cloud.patchesps_, p)
      {
        label patchi = cloud.patchesps_[p];

        const polyPatch& patch = mesh.boundaryMesh()[patchi];
        List<Field<scalar> >& pFA = PressureFluxAccumulators_[p];

        forAll(pFA, facei)
        {
            pFA[facei].setSize(patch.size(), 0);
        }
      }
    }

}


template<class CloudType>
void Foam::Constantflow<CloudType>::inflow()
{
    CloudType& cloud(this->owner());

    const polyMesh& mesh(cloud.mesh());

    const scalar deltaT = mesh.time().deltaTValue();

    Random& rndGen(cloud.rndGen());

    const labelList& patchesog = cloud.patchesog_;

    const labelList& patchescf = cloud.patchescf_;

    const labelList& patchesps = cloud.patchesps_;

    scalar sqrtPi = sqrt(pi);

    const volScalarField::GeometricBoundaryField& boundaryT
    (
        cloud.boundaryT().boundaryField()
    );

    const volVectorField::GeometricBoundaryField& boundaryU
    (
        cloud.boundaryU().boundaryField()
    );

    //Info << 1 << nl << endl;

    if(!patchesog.empty())
    {
    // Insert particles at the outgassing patches
      forAll(patchesog, p)
      {
        label ogparticlesInserted = 0;

        label patchi = patchesog[p];

        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        // Add mass to the accumulators.  negative face area dotted with the
        // velocity to point flux into the domain.

        // Take a reference to the particleFluxAccumulator for this patch
        List<Field<scalar> >& pFA = outgassFluxAccumulators_[p];

        forAll(pFA, i)
        {
            if (min(boundaryT[patchi]) < SMALL)
            {
                FatalErrorIn ("Foam::Constantflow<CloudType>::inflow()")
                    << "Zero boundary temperature detected, check boundaryT "
                    << "condition." << nl
                    << nl << abort(FatalError);
            }
            
            // Read from the dictionnary the numberdensity for the inlet
            scalar nIn = OutgassRate_[p][i];

            pFA[i] += mag(patch.faceAreas())*nIn*deltaT;

            // output the number of particle inserted for the patch
            //Info << "number inserted for patch  <" 
                //<<  patch.name() << ">  is  :" << pFA[i] << endl;

        }

        forAll(patch, pFI)
        {
            // Loop over all faces as the outer loop to avoid calculating
            // geometrical properties multiple times.

            const face& f = patch[pFI];

            label globalFaceIndex = pFI + patch.start();

            label cellI = mesh.faceOwner()[globalFaceIndex];

            const vector& fC = patch.faceCentres()[pFI];

            scalar fA = mag(patch.faceAreas()[pFI]);

            List<tetIndices> faceTets = polyMeshTetDecomposition::faceTetIndices
            (
                mesh,
                globalFaceIndex,
                cellI
            );

            // Cumulative triangle area fractions
            List<scalar> cTriAFracs(faceTets.size(), 0.0);

            scalar previousCummulativeSum = 0.0;

            forAll(faceTets, triI)
            {
                const tetIndices& faceTetIs = faceTets[triI];

                cTriAFracs[triI] =
                    faceTetIs.faceTri(mesh).mag()/fA
                  + previousCummulativeSum;

                previousCummulativeSum = cTriAFracs[triI];
            }

            // Force the last area fraction value to 1.0 to avoid any
            // rounding/non-flat face errors giving a value < 1.0
            cTriAFracs.last() = 1.0;

            // Normal unit vector *negative* so normal is pointing into the
            // domain
            vector n = patch.faceAreas()[pFI];
            n /= -mag(n);

            // Wall tangential unit vector. Use the direction between the
            // face centre and the first vertex in the list
            vector t1 = fC - (mesh.points()[f[0]]);
            t1 /= mag(t1);

            // Other tangential unit vector.  Rescaling in case face is not
            // flat and n and t1 aren't perfectly orthogonal
            vector t2 = n^t1;
            t2 /= mag(t2);

            scalar faceTemperature = boundaryT[patchi][pFI];

            const vector& faceVelocity = boundaryU[patchi][pFI];

            forAll(pFA, i)
            {
                scalar& faceAccumulator = pFA[i][pFI];

                // Number of whole particles to insert, if negative then zero
                label nI = max(label(faceAccumulator), 0);

                // Add another particle with a probability proportional to the
                // remainder of taking the integer part of faceAccumulator
                if ((faceAccumulator - nI) > rndGen.scalar01())
                {
                    nI++;
                }

                faceAccumulator -= nI;

                label typeId = moleculeTypeIds_[i];

                scalar mass = cloud.constProps(typeId).mass();

                for (label i = 0; i < nI; i++)
                {
                    // Choose a triangle to insert on, based on their relative
                    // area

                    scalar triSelection = rndGen.scalar01();

                    // Selected triangle
                    label selectedTriI = -1;

                    forAll(cTriAFracs, triI)
                    {
                        selectedTriI = triI;

                        if (cTriAFracs[triI] >= triSelection)
                        {
                            break;
                        }
                    }

                    // Randomly distribute the points on the triangle.

                    const tetIndices& faceTetIs = faceTets[selectedTriI];

                    point p = faceTetIs.faceTri(mesh).randomPoint(rndGen);

                    // Velocity generation

                    scalar mostProbableSpeed
                    (
                        cloud.maxwellianMostProbableSpeed
                        (
                            faceTemperature,
                            mass
                        )
                    );

                    scalar sCosTheta = (faceVelocity & n)/mostProbableSpeed;

                    // Coefficients required for Bird eqn 12.5
                    scalar uNormProbCoeffA =
                        sCosTheta + sqrt(sqr(sCosTheta) + 2.0);

                    scalar uNormProbCoeffB =
                        0.5*
                        (
                            1.0
                          + sCosTheta*(sCosTheta - sqrt(sqr(sCosTheta) + 2.0))
                        );

                    // Equivalent to the QA value in Bird's DSMC3.FOR
                    scalar randomScaling = 3.0;

                    if (sCosTheta < -3)
                    {
                        randomScaling = mag(sCosTheta) + 1;
                    }

                    scalar P = -1;

                    // Normalised candidates for the normal direction velocity
                    // component
                    scalar uNormal;
                    scalar uNormalThermal;

                    // Select a velocity using Bird eqn 12.5
                    do
                    {
                        uNormalThermal =
                            randomScaling*(2.0*rndGen.scalar01() - 1);

                        uNormal = uNormalThermal + sCosTheta;

                        if (uNormal < 0.0)
                        {
                            P = -1;
                        }
                        else
                        {
                            P = 2.0*uNormal/uNormProbCoeffA
                               *exp(uNormProbCoeffB - sqr(uNormalThermal));
                        }

                    } while (P < rndGen.scalar01());

                    vector U =
                        sqrt(physicoChemical::k.value()*faceTemperature/mass)
                       *(
                            rndGen.GaussNormal()*t1
                          + rndGen.GaussNormal()*t2
                        )
                      + (t1 & faceVelocity)*t1
                      + (t2 & faceVelocity)*t2
                      + mostProbableSpeed*uNormal*n;

                    scalar Ei = cloud.equipartitionInternalEnergy
                    (
                        faceTemperature,
                        cloud.constProps(typeId).internalDegreesOfFreedom()
                    );

                    cloud.addNewParcel
                    (
                        p,
                        U,
                        Ei,
                        cellI,
                        globalFaceIndex,
                        faceTetIs.tetPt(),
                        typeId
                    );

                    ogparticlesInserted++;
                }
            }
        }

        reduce(ogparticlesInserted, sumOp<label>());

        Info<< "    Particles inserted for outgass patch "<< p << " is "
            << ogparticlesInserted << endl;
      }
    }

    // Insert particles at the constant flow rate inlet
    if(!patchescf.empty())
    {
      //const List<scalar>& cfarea = cloud.constantflowarea_;

      //Info<< "cf area      " << cfarea << nl << endl;
      OutFlowNo_ = cloud.outparticleno_;
	  	  	  
      forAll(OutFlowNo_, cfi)
      {
         forAll(OutFlowNo_[cfi], tyj)
         {
             reduce(OutFlowNo_[cfi][tyj], sumOp<scalar>());
         }
      }

      forAll(patchescf, p)
      {
        label cfrparticlesInserted = 0;

        label patchi = patchescf[p];

        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        if(mesh.time().timeIndex() == 1)
        {
            scalar tmpcfnoden = gAverage(Calnodensity(patch, cloud, mesh));

            CfrhoN_[p] = Field<scalar>(patch.size(), tmpcfnoden);
        }
        else if(mesh.time().timeIndex()-mesh.time().startTimeIndex() > 1)
        {
            CfrhoN_[p] = 0.999*CfrhoN_[p]
                       + 0.001*Calnodensity(patch, cloud, mesh);
        }
        else
        {
            CfrhoN_[p] = Calnodensity(patch, cloud, mesh);
        }

        Field<vector> macrovelocity;

        if(mesh.time().timeIndex() > 1)
        {
            macrovelocity = Calmacrovelf(patch, cloud, mesh);
        }
        else
        {
            vector macrotmp = Calmacrovel(patch, cloud, mesh);

            Field<vector> ncf = -patch.faceAreas()/mag(patch.faceAreas());

            macrovelocity = mag(macrotmp & ncf)*ncf;
        }

        if((mesh.time().timeIndex()-mesh.time().startTimeIndex()) >1)
        {
            MacroVel_[p] = 0.999*MacroVel_[p]+0.001*macrovelocity;
        }
        else
        {
            MacroVel_[p] = macrovelocity;
        }

        //Info << "    the modified magmacrovel of the constant patch   "
        //     << MacroVel_[p] << endl;
	   
        const scalar& mass = cfmass_[p];

        if (min(boundaryT[patchi]) < SMALL)
        {
            FatalErrorIn ("Foam::Constantflow<CloudType>::inflow()")
                << "Zero boundary temperature detected, check boundaryT "
                << "condition." << nl
                << nl << abort(FatalError);
        }

        scalarField mostProbableSpeed
        (
            cloud.maxwellianMostProbableSpeed
            (
                boundaryT[patchi],
                mass
            )
        );

        scalarField sCosTheta
        (
           (MacroVel_[p] & -patch.faceAreas()/mag(patch.faceAreas()))
          /mostProbableSpeed
        );

        scalarField tmpdis
        (
            mag(patch.faceAreas())*CfrhoN_[p]
           *mostProbableSpeed
           *(
                exp(-sqr(sCosTheta))
              + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
            )
        );

        scalar dissum = gSum(tmpdis);

        if(dissum != 0)
        {
            tmpdis /= dissum;
        }

        List<Field<scalar> >& pFA = constantFluxAccumulators_[p];

        forAll(moleculeTypeIds_, i)
        {
            scalar nIn = ConstantFlowNo_[p][i]*deltaT + OutFlowNo_[p][i];

            Field<scalar> nInFld = nIn*tmpdis;

            if(mesh.time().timeIndex()-mesh.time().startTimeIndex() > 2)
            {
                CfNoDens_[p][i] = 0.999*CfNoDens_[p][i] + 0.001*nInFld;
            }
            else
            {
                CfNoDens_[p][i] = nInFld;
            }

            pFA[i] += CfNoDens_[p][i];
        }

        // output the number of particle inserted for the patch
        // Info << "number inserted for patch  <" 
        //     <<  patch.name() << ">  is  :" << pFA[i] << endl;

        forAll(patch, pFI)
        {
            // Loop over all faces as the outer loop to avoid calculating
            // geometrical properties multiple times.

            const face& f = patch[pFI];

            label globalFaceIndex = pFI + patch.start();

            label cellI = mesh.faceOwner()[globalFaceIndex];

            const vector& fC = patch.faceCentres()[pFI];

            scalar fA = mag(patch.faceAreas()[pFI]);

            List<tetIndices> faceTets = polyMeshTetDecomposition::faceTetIndices
            (
                mesh,
                globalFaceIndex,
                cellI
            );

            // Cumulative triangle area fractions
            List<scalar> cTriAFracs(faceTets.size(), 0.0);

            scalar previousCummulativeSum = 0.0;

            //Pout << " test test 1" << nl << endl;

            forAll(faceTets, triI)
            {
                const tetIndices& faceTetIs = faceTets[triI];

                cTriAFracs[triI] =
                    faceTetIs.faceTri(mesh).mag()/fA
                  + previousCummulativeSum;

                previousCummulativeSum = cTriAFracs[triI];
            }

            // Force the last area fraction value to 1.0 to avoid any
            // rounding/non-flat face errors giving a value < 1.0
            cTriAFracs.last() = 1.0;

            // Normal unit vector *negative* so normal is pointing into the
            // domain
            vector n = patch.faceAreas()[pFI];
            n /= -mag(n);

            // Wall tangential unit vector. Use the direction between the
            // face centre and the first vertex in the list
            vector t1 = fC - (mesh.points()[f[0]]);
            t1 /= mag(t1);

            // Other tangential unit vector.  Rescaling in case face is not
            // flat and n and t1 aren't perfectly orthogonal
            vector t2 = n^t1;
            t2 /= mag(t2);

            scalar faceTemperature = boundaryT[patchi][pFI];

            const vector& faceVelocity = MacroVel_[p][pFI];

            //Pout << " test test 2 " << nl << endl;

            forAll(pFA, i)
            {
                scalar& faceAccumulator = pFA[i][pFI];

                // Number of whole particles to insert, if negative then zero
                label nI = max(label(faceAccumulator), 0);

                // Add another particle with a probability proportional to the
                // remainder of taking the integer part of faceAccumulator
                if ((faceAccumulator - nI) > rndGen.scalar01())
                {
                    nI++;
                }

                faceAccumulator -= nI;

                label typeId = moleculeTypeIds_[i];

                scalar mass = cloud.constProps(typeId).mass();

                for (label i = 0; i < nI; i++)
                {
                    // Choose a triangle to insert on, based on their relative
                    // area

                    scalar triSelection = rndGen.scalar01();

                    // Selected triangle
                    label selectedTriI = -1;

                    forAll(cTriAFracs, triI)
                    {
                        selectedTriI = triI;

                        if (cTriAFracs[triI] >= triSelection)
                        {
                            break;
                        }
                    }

                    // Randomly distribute the points on the triangle.

                    const tetIndices& faceTetIs = faceTets[selectedTriI];

                    point p = faceTetIs.faceTri(mesh).randomPoint(rndGen);

                    // Velocity generation

                    scalar mostProbableSpeed
                    (
                        cloud.maxwellianMostProbableSpeed
                        (
                            faceTemperature,
                            mass
                        )
                    );

                    scalar sCosTheta = (faceVelocity & n)/mostProbableSpeed;

                    // Coefficients required for Bird eqn 12.5
                    scalar uNormProbCoeffA =
                        sCosTheta + sqrt(sqr(sCosTheta) + 2.0);

                    scalar uNormProbCoeffB =
                        0.5*
                        (
                            1.0
                          + sCosTheta*(sCosTheta - sqrt(sqr(sCosTheta) + 2.0))
                        );

                    // Equivalent to the QA value in Bird's DSMC3.FOR
                    scalar randomScaling = 3.0;

                    if (sCosTheta < -3)
                    {
                        randomScaling = mag(sCosTheta) + 1;
                    }

                    scalar P = -1;

                    // Normalised candidates for the normal direction velocity
                    // component
                    scalar uNormal;
                    scalar uNormalThermal;

                    // Select a velocity using Bird eqn 12.5
                    do
                    {
                        uNormalThermal =
                            randomScaling*(2.0*rndGen.scalar01() - 1);

                        uNormal = uNormalThermal + sCosTheta;

                        if (uNormal < 0.0)
                        {
                            P = -1;
                        }
                        else
                        {
                            P = 2.0*uNormal/uNormProbCoeffA
                               *exp(uNormProbCoeffB - sqr(uNormalThermal));
                        }

                    } while (P < rndGen.scalar01());

                    vector U =
                        sqrt(physicoChemical::k.value()*faceTemperature/mass)
                       *(
                            rndGen.GaussNormal()*t1
                          + rndGen.GaussNormal()*t2
                        )
                      + (t1 & faceVelocity)*t1
                      + (t2 & faceVelocity)*t2
                      + mostProbableSpeed*uNormal*n;

                    scalar Ei = cloud.equipartitionInternalEnergy
                    (
                        faceTemperature,
                        cloud.constProps(typeId).internalDegreesOfFreedom()
                    );

                    cloud.addNewParcel
                    (
                        p,
                        U,
                        Ei,
                        cellI,
                        globalFaceIndex,
                        faceTetIs.tetPt(),
                        typeId
                    );

                    cfrparticlesInserted++;
                }
            }
        }

        reduce(cfrparticlesInserted, sumOp<label>());

        Info<< "    Particles inserted for constant flow patch " << p << " is "
            << cfrparticlesInserted << endl;
      }
    }

    if (!patchesps.empty())
    {
      forAll(patchesps, p)
      {
        label psparticlesInserted = 0;

        label patchi = patchesps[p];

        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        // Add mass to the accumulators.  negative face area dotted with the
        // velocity to point flux into the domain.

        // Take a reference to the particleFluxAccumulator for this patch
        List<Field<scalar> >& pFA = PressureFluxAccumulators_[p];

        Field<vector> psmacrovel = Calmacrovelf(patch, cloud, mesh);

        if(mesh.time().timeIndex() != (mesh.time().startTimeIndex()+1))
        {
            Psmacrovel_[p] = 0.999*Psmacrovel_[p] + 0.001*psmacrovel;
        }
        else
        {
            Psmacrovel_[p] = psmacrovel;
        }

        // Info << Psmacrovel_[p] << nl << endl;

        forAll(pFA, i)
        {
            label typeId = moleculeTypeIds_[i];

            scalar mass = cloud.constProps(typeId).mass();

            // Info << mass << nl << endl;

            if (min(boundaryT[patchi]) < SMALL)
            {
                FatalErrorIn ("Foam::Constantflow<CloudType>::inflow()")
                    << "Zero boundary temperature detected, check boundaryT "
                    << "condition." << nl
                    << nl << abort(FatalError);
            }

            //Info << boundaryT[patchi] << nl << endl;

            scalarField mostProbableSpeed
            (
                cloud.maxwellianMostProbableSpeed
                (
                    boundaryT[patchi],
                    mass
                )
            );

            // Dotting boundary velocity with the face unit normal
            // (which points out of the domain, so it must be
            // negated), dividing by the most probable speed to form
            // molecularSpeedRatio * cosTheta

            scalarField sCosTheta
            (
                (Psmacrovel_[p] & -patch.faceAreas()/mag(patch.faceAreas()))
              / mostProbableSpeed
            );

            // Read from the dictionnary the numberdensity for the pressure inlet
            scalar nIn = PressureNoDens_[p][i];

            //Info << nIn << nl << endl;

            // the type of faceAreas() is subfield
            pFA[i] +=
                    mag(patch.faceAreas())*nIn*deltaT
                   *mostProbableSpeed
                   *(
                        exp(-sqr(sCosTheta))
                      + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                    )
                   /(2.0*sqrtPi);
        }

        forAll(patch, pFI)
        {
            // Loop over all faces as the outer loop to avoid calculating
            // geometrical properties multiple times.

            const face& f = patch[pFI];

            label globalFaceIndex = pFI + patch.start();

            label cellI = mesh.faceOwner()[globalFaceIndex];

            const vector& fC = patch.faceCentres()[pFI];

            scalar fA = mag(patch.faceAreas()[pFI]);

            List<tetIndices> faceTets = polyMeshTetDecomposition::faceTetIndices
            (
                mesh,
                globalFaceIndex,
                cellI
            );

            // Cumulative triangle area fractions
            List<scalar> cTriAFracs(faceTets.size(), 0.0);

            scalar previousCummulativeSum = 0.0;

            forAll(faceTets, triI)
            {
                const tetIndices& faceTetIs = faceTets[triI];

                cTriAFracs[triI] =
                    faceTetIs.faceTri(mesh).mag()/fA
                  + previousCummulativeSum;

                previousCummulativeSum = cTriAFracs[triI];
            }

            // Force the last area fraction value to 1.0 to avoid any
            // rounding/non-flat face errors giving a value < 1.0
            cTriAFracs.last() = 1.0;

            // Normal unit vector *negative* so normal is pointing into the
            // domain
            vector n = patch.faceAreas()[pFI];
            n /= -mag(n);

            // Wall tangential unit vector. Use the direction between the
            // face centre and the first vertex in the list
            vector t1 = fC - (mesh.points()[f[0]]);
            t1 /= mag(t1);

            // Other tangential unit vector.  Rescaling in case face is not
            // flat and n and t1 aren't perfectly orthogonal
            vector t2 = n^t1;
            t2 /= mag(t2);

            scalar faceTemperature = boundaryT[patchi][pFI];

            const vector& faceVelocity = Psmacrovel_[p][pFI];

            forAll(pFA, i)
            {
                scalar& faceAccumulator = pFA[i][pFI];

                // Number of whole particles to insert, if negative then zero
                label nI = max(label(faceAccumulator), 0);

                // Add another particle with a probability proportional to the
                // remainder of taking the integer part of faceAccumulator
                if ((faceAccumulator - nI) > rndGen.scalar01())
                {
                    nI++;
                }

                faceAccumulator -= nI;

                label typeId = moleculeTypeIds_[i];

                scalar mass = cloud.constProps(typeId).mass();

                for (label i = 0; i < nI; i++)
                {
                    // Choose a triangle to insert on, based on their relative
                    // area

                    scalar triSelection = rndGen.scalar01();

                    // Selected triangle
                    label selectedTriI = -1;

                    forAll(cTriAFracs, triI)
                    {
                        selectedTriI = triI;

                        if (cTriAFracs[triI] >= triSelection)
                        {
                            break;
                        }
                    }

                    // Randomly distribute the points on the triangle.

                    const tetIndices& faceTetIs = faceTets[selectedTriI];

                    point p = faceTetIs.faceTri(mesh).randomPoint(rndGen);

                    // Velocity generation

                    scalar mostProbableSpeed
                    (
                        cloud.maxwellianMostProbableSpeed
                        (
                            faceTemperature,
                            mass
                        )
                    );

                    scalar sCosTheta = (faceVelocity & n)/mostProbableSpeed;

                    // Coefficients required for Bird eqn 12.5
                    scalar uNormProbCoeffA =
                        sCosTheta + sqrt(sqr(sCosTheta) + 2.0);

                    scalar uNormProbCoeffB =
                        0.5*
                        (
                            1.0
                          + sCosTheta*(sCosTheta - sqrt(sqr(sCosTheta) + 2.0))
                        );

                    // Equivalent to the QA value in Bird's DSMC3.FOR
                    scalar randomScaling = 3.0;

                    if (sCosTheta < -3)
                    {
                        randomScaling = mag(sCosTheta) + 1;
                    }

                    scalar P = -1;

                    // Normalised candidates for the normal direction velocity
                    // component
                    scalar uNormal;
                    scalar uNormalThermal;

                    // Select a velocity using Bird eqn 12.5
                    do
                    {
                        uNormalThermal =
                            randomScaling*(2.0*rndGen.scalar01() - 1);

                        uNormal = uNormalThermal + sCosTheta;

                        if (uNormal < 0.0)
                        {
                            P = -1;
                        }
                        else
                        {
                            P = 2.0*uNormal/uNormProbCoeffA
                               *exp(uNormProbCoeffB - sqr(uNormalThermal));
                        }

                    } while (P < rndGen.scalar01());

                    vector U =
                        sqrt(physicoChemical::k.value()*faceTemperature/mass)
                       *(
                            rndGen.GaussNormal()*t1
                          + rndGen.GaussNormal()*t2
                        )
                      + (t1 & faceVelocity)*t1
                      + (t2 & faceVelocity)*t2
                      + mostProbableSpeed*uNormal*n;

                    scalar Ei = cloud.equipartitionInternalEnergy
                    (
                        faceTemperature,
                        cloud.constProps(typeId).internalDegreesOfFreedom()
                    );

                    cloud.addNewParcel
                    (
                        p,
                        U,
                        Ei,
                        cellI,
                        globalFaceIndex,
                        faceTetIs.tetPt(),
                        typeId
                    );

                    psparticlesInserted++;
                }
            }
        }

        reduce(psparticlesInserted, sumOp<label>());

        Info << "    Particles inserted for pressure inlet patch " << p
             << " = "<< psparticlesInserted << endl;
      }
    }

}


// Calculate the number density near the patch
template<class CloudType>
Foam::Field<Foam::scalar> Foam::Constantflow<CloudType>::Calnodensity
(
    const polyPatch& nodensitypatch,
    const CloudType& cloud,
    const polyMesh& mesh
)
{
    Field<scalar> tmpdensity = Field<scalar>
    (
        nodensitypatch.size(), 
        0.0
    );

    if
    (
      (mesh.time().timeIndex()-mesh.time().startTimeIndex() > 1) 
      || (mesh.time().timeIndex() ==1)
    )
    {
        const scalarField& rhoN = cloud.rhoN().internalField();

        // Calculate the macrospeed in neighbouring cell
        forAll(nodensitypatch, pFI)
        {
            label globalFaceIndex = pFI + nodensitypatch.start();

            label cellI = mesh.faceOwner()[globalFaceIndex];

            tmpdensity[pFI] = rhoN[cellI];
        }
    }
    else
    {
        const scalarField& MrhoNIn = MrhoN_.internalField();

        forAll(nodensitypatch, pFI)
        {
            label globalFaceIndex = pFI + nodensitypatch.start();

            label cellI = mesh.faceOwner()[globalFaceIndex];

            tmpdensity[pFI] = MrhoNIn[cellI];
        }
    }

    return tmpdensity;
}


// Calculate the number density near the patch
template<class CloudType>
Foam::List<Foam::Field<Foam::scalar> > Foam::Constantflow<CloudType>::Calcfnoden
(
    const polyPatch& nodensitypatch,
    const CloudType& cloud,
    const polyMesh& mesh,
    const label p
)
{
    Field<scalar> tmpmoment = Field<scalar>
    (
        nodensitypatch.size(), 
        0
    );

    const scalarField& rhoN = cloud.rhoN().internalField();

    //const scalarField& rhoM = cloud.rhoM().internalField();

    //const vectorField& momentum = cloud.momentum().internalField();

    // Calculate the macrospeed in neighbouring cell
    forAll(nodensitypatch, pFI)
    {
        label globalFaceIndex = pFI + nodensitypatch.start();

        label cellI = mesh.faceOwner()[globalFaceIndex];

	if(mesh.time().timeIndex() - mesh.time().startTimeIndex() > 1)
	{
	    CfrhoN_[p][pFI] = 0.999*CfrhoN_[p][pFI]+0.001*rhoN[cellI];
	}
        
        tmpmoment[pFI] = mag(MacroVel_[p][pFI] & nodensitypatch.faceAreas()[pFI])
                       *CfrhoN_[p][pFI];
    }

    scalar totalmomentum = gSum(tmpmoment);

    if(mag(totalmomentum) > VSMALL)
    {
        tmpmoment =mag(tmpmoment)/totalmomentum;
    }
    else
    {
       tmpmoment = mag(nodensitypatch.faceAreas())
                  /gSum(mag(nodensitypatch.faceAreas()));
    }

	
    List<Field<scalar> > tmpnoden = List<Field<scalar> >
    (
        moleculeTypeIds_.size(),
        tmpmoment
    );

    return tmpnoden;
}

// Calculate the macro velocity near the patch
template<class CloudType>
Foam::vector Foam::Constantflow<CloudType>::Calmacrovel
(
    const polyPatch& macrovelpatch,
    const CloudType& cloud,
    const polyMesh& mesh
)
{
    // The average number density of particles in the neighbouring cell
    vector macrovelocity = vector::zero;

    double tmpmasstot = 0.0;

    // Calculate the macrospeed in neighbouring cell
    forAll(macrovelpatch, pFI)
    {
        label globalFaceIndex = pFI + macrovelpatch.start();

        label cellI = mesh.faceOwner()[globalFaceIndex];

        const DynamicList<typename CloudType::parcelType* >& particlelist
        (
                cloud.cellOccupancy()[cellI]
        );

        forAll(particlelist, pt)
        {
            typename CloudType::parcelType& ptmp = *particlelist[pt];

            macrovelocity += cloud.constProps(ptmp.typeId()).mass()
                            *ptmp.U();

            tmpmasstot += cloud.constProps(ptmp.typeId()).mass();
        }
    }

    reduce(macrovelocity, sumOp<vector>());

    reduce(tmpmasstot, sumOp<scalar>());
   
    if(tmpmasstot > VSMALL)
    {
        macrovelocity = macrovelocity/tmpmasstot;
    }

    return macrovelocity;
}


// Calculate the macro velocity near the patch face
template<class CloudType>
Foam::Field<Foam::vector> Foam::Constantflow<CloudType>::Calmacrovelf
(
    const polyPatch& macrovelpatch,
    const CloudType& cloud,
    const polyMesh& mesh
)
{
    // The average number density of particles in the neighbouring cell
    Field<vector> macrovelocity = Field<vector>
                                  (
                                     macrovelpatch.size(), vector::zero
                                  );

    //Info << "1" << endl;
    if(mesh.time().timeIndex()-mesh.time().startTimeIndex() > 1)
    {
        const scalarField& rhoM = cloud.rhoM().internalField();

        const vectorField& momentum = cloud.momentum().internalField();

        // Calculate the macrospeed in neighbouring cell
        forAll(macrovelpatch, pFI)
        {
            label globalFaceIndex = pFI + macrovelpatch.start();

            label cellI = mesh.faceOwner()[globalFaceIndex];

            const scalar rhonf = rhoM[cellI];
            const vector momentf = momentum[cellI];
        
            if(rhonf > VSMALL)
            {
                macrovelocity[pFI] = momentf/rhonf;
            }
        }
    }
    else
    {
        const vectorField& UmacroIn = Umacro_.internalField();

        forAll(macrovelpatch, pFI)
        {
          label globalFaceIndex = pFI + macrovelpatch.start();

          label cellI = mesh.faceOwner()[globalFaceIndex];
        
          macrovelocity[pFI] = UmacroIn[cellI];
        }
    }

    return macrovelocity;
}

// Calculate the pressure near the patch
template<class CloudType>
Foam::scalar Foam::Constantflow<CloudType>::Calpressure
(
    const polyPatch& pressurepatch,
    const CloudType& cloud,
    const polyMesh& mesh,
    const scalar& patchvol
)
{
    scalar massvelo2 = 0.0;

    vector massvelo  = vector::zero;

    scalar masstotal = 0.0;

    forAll(pressurepatch, pFI)
    {
        label globalFaceIndex = pFI + pressurepatch.start();

        label cellI = mesh.faceOwner()[globalFaceIndex];

        const DynamicList<typename CloudType::parcelType* >& particlelist
        (
            cloud.cellOccupancy()[cellI]
        );

        forAll(particlelist, j)
        {
            typename CloudType::parcelType& ptmp = *particlelist[j];

            massvelo2 += cloud.constProps(ptmp.typeId()).mass()
                           *(ptmp.U() & ptmp.U());

            massvelo  += cloud.constProps(ptmp.typeId()).mass()*ptmp.U();

            masstotal += cloud.constProps(ptmp.typeId()).mass();
        }
    }

    //Pout<<"massvelo2  "<< massvelo2[i] << nl << "massvelo   "<< massvelo[i]
    //    << nl <<"masstotal  "<< masstotal[i] << nl << endl;

    reduce(massvelo2, sumOp<scalar>());

    reduce(massvelo, sumOp<vector>());

    reduce(masstotal, sumOp<scalar>());

    //Info<< "massvelo2             "<< massvelo2[i] << nl
    //    << "massvelo              "<< massvelo[i]  << nl
    //    << "masstotal             "<< masstotal[i] << nl << endl;

    vector avevelo;

    if (masstotal !=0)
    {
        avevelo = massvelo/masstotal;
    }
    else
    {
        avevelo = vector::zero;
    }

    scalar aveandsum = 2*avevelo & massvelo;

    scalar avevelo2 = avevelo & avevelo;
            
    scalar pretmp = cloud.nParticle()
                   *(massvelo2 - aveandsum + masstotal*avevelo2)/(3*patchvol);

    return pretmp;

}


// Calculate the macro velocity near the patch
template<class CloudType>
Foam::scalar Foam::Constantflow<CloudType>::Patchvolume
(
    const polyPatch& patch,
    const polyMesh& mesh
)
{
    scalar patchvol = 0.0;

    forAll(patch, pFI)
    {
        label globalFaceIndex = pFI + patch.start();

        label cellI = mesh.faceOwner()[globalFaceIndex];

        patchvol += mesh.cellVolumes()[cellI];
    }

    //Pout<< "Constantflow volume " << i << "      "
    //    << constantflowvol_[i] << nl;

    reduce(patchvol, sumOp<scalar>());

    return patchvol;

}
// ************************************************************************* //
