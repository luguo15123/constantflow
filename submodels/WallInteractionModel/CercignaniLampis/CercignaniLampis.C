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

#include "CercignaniLampis.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CercignaniLampis<CloudType>::CercignaniLampis
(
    const dictionary& dict,
    CloudType& cloud
)
:
    WallInteractionModel<CloudType>(dict, cloud, typeName),
    TangentialCoeffs_(readScalar(this->coeffDict().lookup("TangentialCoeffs"))),
    NormalCoeffs_(readScalar(this->coeffDict().lookup("NormalCoeffs")))
{
    //used for debugging
    Info << TangentialCoeffs_<< nl << NormalCoeffs_ << nl << endl; 
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CercignaniLampis<CloudType>::~CercignaniLampis()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::CercignaniLampis<CloudType>::correct
(
    typename CloudType::parcelType& p,
    const wallPolyPatch& wpp
)
{
    vector& U = p.U();

    scalar& Ei = p.Ei();

    label typeId = p.typeId();

    label wppIndex = wpp.index();

    label wppLocalFace = wpp.whichFace(p.face());

    vector nw = p.normal();

    nw /= mag(nw);

    // Normal velocity magnitude
    // is it necessary to substract the velocity component of the wall boundary
    scalar U_dot_nw = U & nw;

    CloudType& cloud(this->owner());

    Random& rndGen(cloud.rndGen());

    // tangential velocity magnitude
    vector Ut = U - U_dot_nw*nw;

    scalar U_dot_tw1 = mag(Ut);
    
    vector tw1 = vector(0.0,0.0,0.0);

    if (U_dot_tw1 > 0.0)
    {
        tw1 = Ut/U_dot_tw1;
    }
    else
    {
        if (nw.x() > 0.0)
        {
            tw1 = vector(-nw.y(),nw.x(),0.0);
        }
        else
        {
            tw1 = vector(0.0,nw.z(),-nw.y());
        }
        tw1 /= mag(tw1);
    }
/*    while (mag(Ut) < SMALL)
    {
        // If the incident velocity is parallel to the face normal, no
        // tangential direction can be chosen.  Add a perturbation to the
        // incoming velocity and recalculate.
        U = vector
            (
                U.x()*(0.8 + 0.2*rndGen.scalar01()),
                U.y()*(0.8 + 0.2*rndGen.scalar01()),
                U.z()*(0.8 + 0.2*rndGen.scalar01())
            );
        U_dot_nw1 = U & nw;
        Ut = U - U_dot_nw1*nw;
    }
*/

    // Other tangential unit vector
    vector tw2 = nw^tw1;

    scalar T = cloud.boundaryT().boundaryField()[wppIndex][wppLocalFace];

    scalar mass = cloud.constProps(typeId).mass();

    direction iDof = cloud.constProps(typeId).internalDegreesOfFreedom();

    //normalization component
    scalar vnormalise = sqrt(2*physicoChemical::k.value()*T/mass);
    //normalized normal component
    scalar tun = mag(U_dot_nw/vnormalise);
    //normalized tangential component
    scalar tut = U_dot_tw1/vnormalise;

    //the mean velocity after reflection
    scalar Unim = sqrt(1-NormalCoeffs_)*tun;
    scalar Utim = (1-TangentialCoeffs_)*tut;
    
    scalar rtn = sqrt((-NormalCoeffs_)*log(rndGen.scalar01()));
    scalar thetan = 2*Foam::constant::mathematical::pi*rndGen.scalar01();
    scalar Unr = sqrt(sqr(Unim)+2*rtn*Unim*cos(thetan)+sqr(rtn));
    
    scalar rtt = sqrt((-TangentialCoeffs_)*(2-TangentialCoeffs_)
              *log(rndGen.scalar01()));
    scalar thetat = 2*Foam::constant::mathematical::pi*rndGen.scalar01();
    scalar Utt1 = Utim + rtt*cos(thetat);
    scalar Utt2 = rtt*sin(thetat);

    U = (-Unr*nw + Utt1*tw1 + Utt2*tw2)*vnormalise;

//    Info << U << nl << endl;
    
//        U =
//            sqrt(physicoChemical::k.value()*T/mass)
//           *(
//                rndGen.GaussNormal()*tw1
//              + rndGen.GaussNormal()*tw2
//              - sqrt(-2.0*log(max(1 - rndGen.scalar01(), VSMALL)))*nw
//            );

    U += cloud.boundaryU().boundaryField()[wppIndex][wppLocalFace];

    Ei = cloud.equipartitionInternalEnergy(T, iDof);

}


// velocity correct function used for 
template<class CloudType>
void Foam::CercignaniLampis<CloudType>::correct
(
    typename CloudType::parcelType& p,
    const polyPatch& wpp
)
{
    vector& U = p.U();

    scalar& Ei = p.Ei();

    label typeId = p.typeId();

    label wppIndex = wpp.index();

    label wppLocalFace = wpp.whichFace(p.face());

    vector nw = p.normal();

    nw /= mag(nw);

    // Normal velocity magnitude
    // is it necessary to substract the velocity component of the wall boundary
    scalar U_dot_nw = U & nw;

    CloudType& cloud(this->owner());

    Random& rndGen(cloud.rndGen());

    // tangential velocity magnitude
    vector Ut = U - U_dot_nw*nw;

    scalar U_dot_tw1 = mag(Ut);
    
    vector tw1 = vector(0.0,0.0,0.0);

    if (U_dot_tw1 > 0.0)
    {
        tw1 = Ut/U_dot_tw1;
    }
    else
    {
        if (nw.x() > 0.0)
        {
            tw1 = vector(-nw.y(),nw.x(),0.0);
        }
        else
        {
            tw1 = vector(0.0,nw.z(),-nw.y());
        }
        tw1 /= mag(tw1);
    }
/*    while (mag(Ut) < SMALL)
    {
        // If the incident velocity is parallel to the face normal, no
        // tangential direction can be chosen.  Add a perturbation to the
        // incoming velocity and recalculate.
        U = vector
            (
                U.x()*(0.8 + 0.2*rndGen.scalar01()),
                U.y()*(0.8 + 0.2*rndGen.scalar01()),
                U.z()*(0.8 + 0.2*rndGen.scalar01())
            );
        U_dot_nw1 = U & nw;
        Ut = U - U_dot_nw1*nw;
    }
*/

    // Other tangential unit vector
    vector tw2 = nw^tw1;

    scalar T = cloud.boundaryT().boundaryField()[wppIndex][wppLocalFace];

    scalar mass = cloud.constProps(typeId).mass();

    direction iDof = cloud.constProps(typeId).internalDegreesOfFreedom();

    //normalization component
    scalar vnormalise = sqrt(2*physicoChemical::k.value()*T/mass);
    //normalized normal component
    scalar tun = mag(U_dot_nw/vnormalise);
    //normalized tangential component
    scalar tut = U_dot_tw1/vnormalise;

    //the mean velocity after reflection
    scalar Unim = sqrt(1-NormalCoeffs_)*tun;
    scalar Utim = (1-TangentialCoeffs_)*tut;
    
    scalar rtn = sqrt((-NormalCoeffs_)*log(rndGen.scalar01()));
    scalar thetan = 2*Foam::constant::mathematical::pi*rndGen.scalar01();
    scalar Unr = sqrt(sqr(Unim)+2*rtn*Unim*cos(thetan)+sqr(rtn));
    
    scalar rtt = sqrt((-TangentialCoeffs_)*(2-TangentialCoeffs_)
              *log(rndGen.scalar01()));
    scalar thetat = 2*Foam::constant::mathematical::pi*rndGen.scalar01();
    scalar Utt1 = Utim + rtt*cos(thetat);
    scalar Utt2 = rtt*sin(thetat);

    U = (-Unr*nw + Utt1*tw1 + Utt2*tw2)*vnormalise;

//    Info << U << nl << endl;
    
//        U =
//            sqrt(physicoChemical::k.value()*T/mass)
//           *(
//                rndGen.GaussNormal()*tw1
//              + rndGen.GaussNormal()*tw2
//              - sqrt(-2.0*log(max(1 - rndGen.scalar01(), VSMALL)))*nw
//            );

    U += cloud.boundaryU().boundaryField()[wppIndex][wppLocalFace];

    Ei = cloud.equipartitionInternalEnergy(T, iDof);

}

// ************************************************************************* //
