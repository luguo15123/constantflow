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

#include "DSMCParcel.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
bool Foam::DSMCParcel<ParcelType>::move(TrackData& td, const scalar trackTime)
{
    typename TrackData::cloudType::parcelType& p =
        static_cast<typename TrackData::cloudType::parcelType&>(*this);

    td.switchProcessor = false;
    td.keepParticle = true;

    const polyMesh& mesh = td.cloud().pMesh();
    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();

    scalar tEnd = (1.0 - p.stepFraction())*trackTime;
    const scalar dtMax = tEnd;

    // For reduced-D cases, the velocity used to track needs to be
    // constrained, but the actual U_ of the parcel must not be
    // altered or used, as it is altered by patch interactions an
    // needs to retain its 3D value for collision purposes.
    vector Utracking = U_;
    //vector Ptracking = p.position_;

    while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
    {
        // Apply correction to position for reduced-D cases
        meshTools::constrainToMeshCentre(mesh, p.position());

        Utracking = U_;

        // Apply correction to velocity to constrain tracking for
        // reduced-D cases
        meshTools::constrainDirection(mesh, mesh.solutionD(), Utracking);

        // Set the Lagrangian time-step
        scalar dt = min(dtMax, tEnd);

        dt *= p.trackToFace(p.position() + dt*Utracking, td);

        tEnd -= dt;

        p.stepFraction() = 1.0 - tEnd/trackTime;

        if (p.onBoundary() && td.keepParticle)
        {
            if (isA<processorPolyPatch>(pbMesh[p.patch(p.face())]))
            {
                td.switchProcessor = true;
            }
        }
    }

    /*  //used for counting the angular distribution at
    // the capillary outlet without outlet reservoirs
    if (!td.keepParticle)   
    {
        // the angular statistic data for the particle trave through the outlet
        // and at each time step, the default number of the particle is the zero
    	FixedList<label,91>& AngCountt = td.cloud().AngCount;

        // find the the label index of the patch the particle is on
	label patchInI = pbMesh.whichPatch(p.face());
	// find the patch that the particle is on
	const polyPatch& patchIn = pbMesh[patchInI];
                
	// if the patch is the outlet patch, then analyze the
	// direction of the velocity with an interval of 1 degree
	if (patchIn.name() == "outlet")
	{
	    // the direction of the velocity is calculated on the base
	    // of the normal direction of the patch is parallel with z
	    // - first the polar angle of the velocity is calculated
	    scalar AngDegs = 180/constant::mathematical::pi
                          *atan(sqrt(sqr(p.U()[0]) + sqr(p.U()[1]))/p.U()[2]);
	    label AngDeg = AngDegs;
	    if ((AngDegs - AngDeg) > 0)
	    {
	        AngDeg++;
	    }
		    
            // increase the number of particle by 1 in the direction
            AngCountt[AngDeg] += 1;
	}
    }
    */

    /*  // used for the one with outlet reservoirs
    if ((Ptracking[2] < 0.016) && (p.position_[2] > 0.016)) 
    {
        // the angular statistic data for the particle trave through the outlet
        // and at each time step, the default number of the particle is the zero
    	FixedList<label,91>& AngCountt = td.cloud().AngCount;

        // find the the label index of the patch the particle is on
	//label patchInI = pbMesh.whichPatch(p.face());
	// find the patch that the particle is on
	//const polyPatch& patchIn = pbMesh[patchInI];
                
	// if the patch is the outlet patch, then analyze the
	// direction of the velocity with an interval of 1 degree
	//if (patchIn.name() == "outlet")
	//{
	    // the direction of the velocity is calculated on the base
	    // of the normal direction of the patch is parallel with z
	    // - first the polar angle of the velocity is calculated
	scalar AngDegs = 180/constant::mathematical::pi
                      *atan(sqrt(sqr(p.U()[0]) + sqr(p.U()[1]))/p.U()[2]);
	label AngDeg = AngDegs;
	if ((AngDegs - AngDeg) > 0)
	{
	        AngDeg++;
	}
		    
            // increase the number of particle by 1 in the direction
        AngCountt[AngDeg] += 1;
	//}
    }
    */
    return td.keepParticle;
}


template<class ParcelType>
template<class TrackData>
bool Foam::DSMCParcel<ParcelType>::hitPatch
(
    const polyPatch&,
    TrackData& td,
    const label,
    const scalar,
    const tetIndices&
)
{
    return false;
}


template<class ParcelType>
template<class TrackData>
void Foam::DSMCParcel<ParcelType>::hitProcessorPatch
(
    const processorPolyPatch&,
    TrackData& td
)
{
    td.switchProcessor = true;
}


template<class ParcelType>
template<class TrackData>
void Foam::DSMCParcel<ParcelType>::hitWallPatch
(
    const wallPolyPatch& wpp,
    TrackData& td,
    const tetIndices& tetIs
)
{
    Random& rndGen(td.cloud().rndGen());

    const scalar deltaT = td.cloud().pMesh().time().deltaTValue();

    //const labelList& patchesrf = cloud.patchesrf_;

    const labelList& patchespo = td.cloud().patchespo_;

    label wppIndex = wpp.index();

    label findwppindex = findIndex(patchespo, wppIndex);

    List<scalar>& postpumpno = td.cloud().postpumpno_;

    List<scalar>& prepumpno = td.cloud().prepumpno_;

    List<scalar>& pumpoutno = td.cloud().pumpoutno_;

    List<scalar>& probablepumppo = td.cloud().probablepumppo_;

    // if the wall is a pumping wall - a wall which has a name with pre "pump"
    if ((findwppindex) != -1)
    {
        // counting the particles that collide with the pumping wall
        postpumpno[findwppindex] += 1;

        scalar propump = probablepumppo[findwppindex];

        if
        (
            td.cloud().mesh().time().timeIndex()
             != (td.cloud().mesh().time().startTimeIndex()+1)
        )
        {
            if (rndGen.scalar01() < propump)
            {
                td.keepParticle = false;
                pumpoutno[findwppindex] += 1;
            }
        }
        else
        {
            if(pumpoutno[findwppindex] < prepumpno[findwppindex])
            {
                td.keepParticle = false;
                pumpoutno[findwppindex] += 1;                
            }
        }
    }

    if(td.keepParticle)
    {

        label wppLocalFace = wpp.whichFace(this->face());

        const scalar fA = mag(wpp.faceAreas()[wppLocalFace]);

        const constantProperties& constProps(td.cloud().constProps(typeId_));

        scalar m = constProps.mass();

        vector nw = wpp.faceAreas()[wppLocalFace];
        nw /= mag(nw);

        scalar U_dot_nw = U_ & nw;

        vector Ut = U_ - U_dot_nw*nw;

        scalar invMagUnfA = 1/max(mag(U_dot_nw)*fA, VSMALL);

        td.cloud().rhoNBF()[wppIndex][wppLocalFace] += invMagUnfA;

        td.cloud().rhoMBF()[wppIndex][wppLocalFace] += m*invMagUnfA;

        td.cloud().linearKEBF()[wppIndex][wppLocalFace] +=
            0.5*m*(U_ & U_)*invMagUnfA;

        td.cloud().internalEBF()[wppIndex][wppLocalFace] += Ei_*invMagUnfA;

        td.cloud().iDofBF()[wppIndex][wppLocalFace] +=
            constProps.internalDegreesOfFreedom()*invMagUnfA;

        td.cloud().momentumBF()[wppIndex][wppLocalFace] += m*Ut*invMagUnfA;

        // pre-interaction energy
        scalar preIE = 0.5*m*(U_ & U_) + Ei_;

        // pre-interaction momentum
        vector preIMom = m*U_;

        td.cloud().wallInteraction().correct
        (
            static_cast<DSMCParcel<ParcelType> &>(*this),
            wpp
        );

        U_dot_nw = U_ & nw;

        Ut = U_ - U_dot_nw*nw;

        invMagUnfA = 1/max(mag(U_dot_nw)*fA, VSMALL);

        td.cloud().rhoNBF()[wppIndex][wppLocalFace] += invMagUnfA;

        td.cloud().rhoMBF()[wppIndex][wppLocalFace] += m*invMagUnfA;

        td.cloud().linearKEBF()[wppIndex][wppLocalFace] +=
            0.5*m*(U_ & U_)*invMagUnfA;

        td.cloud().internalEBF()[wppIndex][wppLocalFace] += Ei_*invMagUnfA;

        td.cloud().iDofBF()[wppIndex][wppLocalFace] +=
            constProps.internalDegreesOfFreedom()*invMagUnfA;

        td.cloud().momentumBF()[wppIndex][wppLocalFace] += m*Ut*invMagUnfA;

        // post-interaction energy
        scalar postIE = 0.5*m*(U_ & U_) + Ei_;

        // post-interaction momentum
        vector postIMom = m*U_;

        scalar deltaQ = td.cloud().nParticle()*(preIE - postIE)/(deltaT*fA);

        vector deltaFD = td.cloud().nParticle()*(preIMom - postIMom)/(deltaT*fA);

        td.cloud().qBF()[wppIndex][wppLocalFace] += deltaQ;

        td.cloud().fDBF()[wppIndex][wppLocalFace] += deltaFD;
    }

}


template<class ParcelType>
template<class TrackData>
void Foam::DSMCParcel<ParcelType>::hitPatch(const polyPatch& pp, TrackData& td)
{
//  how to call the member in submodel
//  td.cloud().wallInteraction()
    typename TrackData::cloudType& cloud(td.cloud());
    
    label ppindex = pp.index();

    label findogindex = findIndex(cloud.patchesog_,ppindex);

    if(findogindex != -1)
    {
        td.cloud().wallInteraction().correct
        (
            static_cast<DSMCParcel<ParcelType> &>(*this),
            pp
        );

        td.keepParticle = true;
    }

    else
    {
        label findcfindex = findIndex(cloud.patchescf_,ppindex);
        if(findcfindex !=-1)
        {
           cloud.outparticleno_[findcfindex][typeId_] += 1;
        }

        td.keepParticle = false;
    }

}


template<class ParcelType>
void Foam::DSMCParcel<ParcelType>::transformProperties(const tensor& T)
{
    ParcelType::transformProperties(T);
    U_ = transform(T, U_);
}


template<class ParcelType>
void Foam::DSMCParcel<ParcelType>::transformProperties
(
    const vector& separation
)
{
    ParcelType::transformProperties(separation);
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "DSMCParcelIO.C"

// ************************************************************************* //
