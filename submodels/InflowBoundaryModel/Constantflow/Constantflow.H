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

Class
    Foam::Constantflow

Description
    Inserting new particles across the faces of a all patched of type
    "patch" for a free stream.  Uniform values number density, temperature
    and velocity sourced face-by-face from the boundaryT and boundaryU fields
    of the cloud.

\*---------------------------------------------------------------------------*/

#ifndef Constantflow_H
#define Constantflow_H

#include "InflowBoundaryModel.H"
#include "polyMesh.H"
#include "dsmcParcel.H"
#include "DSMCCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class Constantflow Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class Constantflow
:
    public InflowBoundaryModel<CloudType>
{
    // Private data

        //- The molecule types to be introduced
        List<label> moleculeTypeIds_;

        //- The outgassing rate at the outgassing patch
        List<List<scalar> > OutgassRate_;

        List<List<Field<scalar> > > outgassFluxAccumulators_;

        //- A List of Lists of Fields specifying carry-over of mass flux from
        // one timestep to the next
        // + Outer List - one inner List for each patch
        // + Inner List - one Field for every species to be introduced
        // + Each field entry corresponding to a face to be injected across
        //   with a particular species
        List<List<Field<scalar> > > constantFluxAccumulators_;

        //- The particle No. corresponding to the constant flow rate
        List<List<scalar> > ConstantFlowNo_;

        //- The particle No.that flow out of the inflow boundary at previous
        //  time step
        List<List<scalar> > OutFlowNo_;

        List<Field<vector> > MacroVel_;

        List<Field<scalar> > CfrhoN_;

        List<scalar> cfmass_;

        //List<Field<scalar> > CfrhoM_;

	//List<Field<vector> > CfMom_;

        volVectorField Umacro_;

        volScalarField MrhoN_;

        //- The macro speed of the inflow particles
        List<List<Field<scalar> > > CfNoDens_;

        //- The pressure inlet accumulator
        List<List<Field<scalar> > > PressureFluxAccumulators_;

        //- The number density corresponding to the pressure
        List<List<scalar> > PressureNoDens_;

        //- The patch volume near the pressure patch
        List<scalar> PressureVol_;

        List<Field<vector> > Psmacrovel_;

public:

    //- Runtime type information
    TypeName("Constantflow");


    // Constructors

        //- Construct from dictionary
        Constantflow
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~Constantflow();


    // Member Functions

        // Mapping

            //- Remap the particles to the correct cells following mesh change
            virtual void autoMap(const mapPolyMesh&);

        //- Introduce particles
        virtual void inflow();
        
        //- Calculate the number density near the inflow patch
        Field<scalar> Calnodensity
               (
                   const polyPatch& nodensitypatch,
                   const CloudType& cloud,
                   const polyMesh& mesh
               );

        //- Calculate the number density near the cf inflow patch
        List<Field<scalar> > Calcfnoden
               (
                   const polyPatch& nodensitypatch,
                   const CloudType& cloud,
                   const polyMesh& mesh,
                   const label p
               );

        //- Calculate the macro speed near the inflow patch
        vector Calmacrovel
               (
                   const polyPatch& macrovelpatch,
                   const CloudType& cloud,
                   const polyMesh& mesh
               );

        Field<vector> Calmacrovelf
                      (
                          const polyPatch& macrovelpatch,
                          const CloudType& cloud,
                          const polyMesh& mesh
                      );

        //- Calculate the pressure near the inflow patch
        scalar Calpressure
               (
                   const polyPatch& pressurepatch,
                   const CloudType& cloud,
                   const polyMesh& mesh,
                   const scalar& patchvol
               );

        //- Calculate the cell volume near the patch
        scalar Patchvolume
               (
                   const polyPatch& patch,
                   const polyMesh& mesh
               );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "Constantflow.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
