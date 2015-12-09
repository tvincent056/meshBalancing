/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Original Copyright (C) 2014 Tyler Voskuilen
     \\/     M anipulation  | Updates Copyright (C) 2015 Timothy Vincent
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "dynamicRefineBalancedFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceInterpolate.H"
#include "volFields.H"
#include "polyTopoChange.H"
#include "surfaceFields.H"
#include "syncTools.H"
#include "pointFields.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicRefineBalancedFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, dynamicRefineBalancedFvMesh, IOobject);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// the PackedBoolList::count method would probably be faster
// since we are only checking for 'true' anyhow
Foam::label Foam::dynamicRefineBalancedFvMesh::count
(
    const PackedBoolList& l,
    const unsigned int val
)
{
    label n = 0;
    forAll(l, i)
    {
        if (l.get(i) == val)
        {
            n++;
        }

        // debug also serves to get-around Clang compiler trying to optimsie
        // out this forAll loop under O3 optimisation
        if (debug)
        {
            Info<< "n=" << n << endl;
        }
    }

    return n;
}


void Foam::dynamicRefineBalancedFvMesh::calculateProtectedCells
(
    PackedBoolList& unrefineableCell
) const
{
    if (protectedCell_.empty())
    {
        unrefineableCell.clear();
        return;
    }

    const labelList& cellLevel = meshCutter_.cellLevel();

    unrefineableCell = protectedCell_;

    // Get neighbouring cell level
    labelList neiLevel(nFaces()-nInternalFaces());

    for (label faceI = nInternalFaces(); faceI < nFaces(); faceI++)
    {
        neiLevel[faceI-nInternalFaces()] = cellLevel[faceOwner()[faceI]];
    }
    syncTools::swapBoundaryFaceList(*this, neiLevel);


    while (true)
    {
        // Pick up faces on border of protected cells
        boolList seedFace(nFaces(), false);

        forAll(faceNeighbour(), faceI)
        {
            label own = faceOwner()[faceI];
            bool ownProtected = unrefineableCell.get(own);
            label nei = faceNeighbour()[faceI];
            bool neiProtected = unrefineableCell.get(nei);

            if (ownProtected && (cellLevel[nei] > cellLevel[own]))
            {
                seedFace[faceI] = true;
            }
            else if (neiProtected && (cellLevel[own] > cellLevel[nei]))
            {
                seedFace[faceI] = true;
            }
        }
        for (label faceI = nInternalFaces(); faceI < nFaces(); faceI++)
        {
            label own = faceOwner()[faceI];
            bool ownProtected = unrefineableCell.get(own);
            if
            (
                ownProtected
             && (neiLevel[faceI-nInternalFaces()] > cellLevel[own])
            )
            {
                seedFace[faceI] = true;
            }
        }

        syncTools::syncFaceList(*this, seedFace, orEqOp<bool>());


        // Extend unrefineableCell
        bool hasExtended = false;

        for (label faceI = 0; faceI < nInternalFaces(); faceI++)
        {
            if (seedFace[faceI])
            {
                label own = faceOwner()[faceI];
                if (unrefineableCell.get(own) == 0)
                {
                    unrefineableCell.set(own, 1);
                    hasExtended = true;
                }

                label nei = faceNeighbour()[faceI];
                if (unrefineableCell.get(nei) == 0)
                {
                    unrefineableCell.set(nei, 1);
                    hasExtended = true;
                }
            }
        }
        for (label faceI = nInternalFaces(); faceI < nFaces(); faceI++)
        {
            if (seedFace[faceI])
            {
                label own = faceOwner()[faceI];
                if (unrefineableCell.get(own) == 0)
                {
                    unrefineableCell.set(own, 1);
                    hasExtended = true;
                }
            }
        }

        if (!returnReduce(hasExtended, orOp<bool>()))
        {
            break;
        }
    }
}


void Foam::dynamicRefineBalancedFvMesh::readDict()
{
    dictionary refineDict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).subDict("dynamicRefineFvMeshCoeffs")
    );

    List<Pair<word> > fluxVelocities = List<Pair<word> >
    (
        refineDict.lookup("correctFluxes")
    );
    // Rework into hashtable.
    correctFluxes_.resize(fluxVelocities.size());
    forAll(fluxVelocities, i)
    {
        correctFluxes_.insert(fluxVelocities[i][0], fluxVelocities[i][1]);
    }

    dumpLevel_ = Switch(refineDict.lookup("dumpLevel"));
}


// Refines cells, maps fields and recalculates (an approximate) flux
Foam::autoPtr<Foam::mapPolyMesh>
Foam::dynamicRefineBalancedFvMesh::refine
(
    const labelList& cellsToRefine
)
{
    // Mesh changing engine.
    polyTopoChange meshMod(*this);

    // Play refinement commands into mesh changer.
    meshCutter_.setRefinement(cellsToRefine, meshMod);

    // Create mesh (with inflation), return map from old to new mesh.
    //autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, true);
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false);

    Info<< "Refined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << globalData().nTotalCells() << " cells." << endl;

    if (debug)
    {
        // Check map.
        for (label faceI = 0; faceI < nInternalFaces(); faceI++)
        {
            label oldFaceI = map().faceMap()[faceI];

            if (oldFaceI >= nInternalFaces())
            {
                FatalErrorIn("dynamicRefineBalancedFvMesh::refine(const labelList&)")
                    << "New internal face:" << faceI
                    << " fc:" << faceCentres()[faceI]
                    << " originates from boundary oldFace:" << oldFaceI
                    << abort(FatalError);
            }
        }
    }

//    // Remove the stored tet base points
//    tetBasePtIsPtr_.clear();
//    // Remove the cell tree
//    cellTreePtr_.clear();

    // Update fields
    updateMesh(map);


    // Move mesh
    /*
    pointField newPoints;
    if (map().hasMotionPoints())
    {
        newPoints = map().preMotionPoints();
    }
    else
    {
        newPoints = points();
    }
    movePoints(newPoints);
    */

    // Correct the flux for modified/added faces. All the faces which only
    // have been renumbered will already have been handled by the mapping.
    {
        const labelList& faceMap = map().faceMap();
        const labelList& reverseFaceMap = map().reverseFaceMap();

        // Storage for any master faces. These will be the original faces
        // on the coarse cell that get split into four (or rather the
        // master face gets modified and three faces get added from the master)
        labelHashSet masterFaces(4*cellsToRefine.size());

        forAll(faceMap, faceI)
        {
            label oldFaceI = faceMap[faceI];

            if (oldFaceI >= 0)
            {
                label masterFaceI = reverseFaceMap[oldFaceI];

                if (masterFaceI < 0)
                {
                    FatalErrorIn
                    (
                        "dynamicRefineBalancedFvMesh::refine(const labelList&)"
                    )   << "Problem: should not have removed faces"
                        << " when refining."
                        << nl << "face:" << faceI << abort(FatalError);
                }
                else if (masterFaceI != faceI)
                {
                    masterFaces.insert(masterFaceI);
                }
            }
        }
        if (debug)
        {
            Pout<< "Found " << masterFaces.size() << " split faces " << endl;
        }

        HashTable<surfaceScalarField*> fluxes
        (
            lookupClass<surfaceScalarField>()
        );
        forAllIter(HashTable<surfaceScalarField*>, fluxes, iter)
        {
            if (!correctFluxes_.found(iter.key()))
            {
                WarningIn("dynamicRefineBalancedFvMesh::refine(const labelList&)")
                    << "Cannot find surfaceScalarField " << iter.key()
                    << " in user-provided flux mapping table "
                    << correctFluxes_ << endl
                    << "    The flux mapping table is used to recreate the"
                    << " flux on newly created faces." << endl
                    << "    Either add the entry if it is a flux or use ("
                    << iter.key() << " none) to suppress this warning."
                    << endl;
                continue;
            }

            const word& UName = correctFluxes_[iter.key()];

            if (UName == "none")
            {
                continue;
            }

            if (UName == "NaN")
            {
                Pout<< "Setting surfaceScalarField " << iter.key()
                    << " to NaN" << endl;

                surfaceScalarField& phi = *iter();

                sigFpe::fillSignallingNan(phi.internalField());

                continue;
            }

            if (debug)
            {
                Pout<< "Mapping flux " << iter.key()
                    << " using interpolated flux " << UName
                    << endl;
            }

            surfaceScalarField& phi = *iter();
            const surfaceScalarField phiU
            (
                fvc::interpolate
                (
                    lookupObject<volVectorField>(UName)
                )
              & Sf()
            );

            // Recalculate new internal faces.
            for (label faceI = 0; faceI < nInternalFaces(); faceI++)
            {
                label oldFaceI = faceMap[faceI];

                if (oldFaceI == -1)
                {
                    // Inflated/appended
                    phi[faceI] = phiU[faceI];
                }
                else if (reverseFaceMap[oldFaceI] != faceI)
                {
                    // face-from-masterface
                    phi[faceI] = phiU[faceI];
                }
            }

            // Recalculate new boundary faces.
            surfaceScalarField::GeometricBoundaryField& bphi =
                phi.boundaryField();
            forAll(bphi, patchI)
            {
                fvsPatchScalarField& patchPhi = bphi[patchI];
                const fvsPatchScalarField& patchPhiU =
                    phiU.boundaryField()[patchI];

                label faceI = patchPhi.patch().start();

                forAll(patchPhi, i)
                {
                    label oldFaceI = faceMap[faceI];

                    if (oldFaceI == -1)
                    {
                        // Inflated/appended
                        patchPhi[i] = patchPhiU[i];
                    }
                    else if (reverseFaceMap[oldFaceI] != faceI)
                    {
                        // face-from-masterface
                        patchPhi[i] = patchPhiU[i];
                    }

                    faceI++;
                }
            }

            // Update master faces
            forAllConstIter(labelHashSet, masterFaces, iter)
            {
                label faceI = iter.key();

                if (isInternalFace(faceI))
                {
                    phi[faceI] = phiU[faceI];
                }
                else
                {
                    label patchI = boundaryMesh().whichPatch(faceI);
                    label i = faceI - boundaryMesh()[patchI].start();

                    const fvsPatchScalarField& patchPhiU =
                        phiU.boundaryField()[patchI];

                    fvsPatchScalarField& patchPhi = bphi[patchI];

                    patchPhi[i] = patchPhiU[i];
                }
            }
        }
    }



    // Update numbering of cells/vertices.
    meshCutter_.updateMesh(map);

    // Update numbering of protectedCell_
    if (protectedCell_.size())
    {
        PackedBoolList newProtectedCell(nCells());

        forAll(newProtectedCell, cellI)
        {
            label oldCellI = map().cellMap()[cellI];
            newProtectedCell.set(cellI, protectedCell_.get(oldCellI));
        }
        protectedCell_.transfer(newProtectedCell);
    }

    // Debug: Check refinement levels (across faces only)
    meshCutter_.checkRefinementLevels(-1, labelList(0));

    return map;
}


// Combines previously split cells, maps fields and recalculates
// (an approximate) flux
Foam::autoPtr<Foam::mapPolyMesh>
Foam::dynamicRefineBalancedFvMesh::unrefine
(
    const labelList& splitPoints
)
{
    polyTopoChange meshMod(*this);

    // Play refinement commands into mesh changer.
    meshCutter_.setUnrefinement(splitPoints, meshMod);


    // Save information on faces that will be combined
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Find the faceMidPoints on cells to be combined.
    // for each face resulting of split of face into four store the
    // midpoint
    Map<label> faceToSplitPoint(3*splitPoints.size());

    {
        forAll(splitPoints, i)
        {
            label pointI = splitPoints[i];

            const labelList& pEdges = pointEdges()[pointI];

            forAll(pEdges, j)
            {
                label otherPointI = edges()[pEdges[j]].otherVertex(pointI);

                const labelList& pFaces = pointFaces()[otherPointI];

                forAll(pFaces, pFaceI)
                {
                    faceToSplitPoint.insert(pFaces[pFaceI], otherPointI);
                }
            }
        }
    }


    // Change mesh and generate map.
    //autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, true);
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false);

    Info<< "Unrefined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << globalData().nTotalCells() << " cells."
        << endl;

    // Update fields
    updateMesh(map);


    // Move mesh
    /*
    pointField newPoints;
    if (map().hasMotionPoints())
    {
        newPoints = map().preMotionPoints();
    }
    else
    {
        newPoints = points();
    }
    movePoints(newPoints);
    */

    // Correct the flux for modified faces.
    {
        const labelList& reversePointMap = map().reversePointMap();
        const labelList& reverseFaceMap = map().reverseFaceMap();

        HashTable<surfaceScalarField*> fluxes
        (
            lookupClass<surfaceScalarField>()
        );
        forAllIter(HashTable<surfaceScalarField*>, fluxes, iter)
        {
            if (!correctFluxes_.found(iter.key()))
            {
                WarningIn("dynamicRefineBalancedFvMesh::refine(const labelList&)")
                    << "Cannot find surfaceScalarField " << iter.key()
                    << " in user-provided flux mapping table "
                    << correctFluxes_ << endl
                    << "    The flux mapping table is used to recreate the"
                    << " flux on newly created faces." << endl
                    << "    Either add the entry if it is a flux or use ("
                    << iter.key() << " none) to suppress this warning."
                    << endl;
                continue;
            }

            const word& UName = correctFluxes_[iter.key()];

            if (UName == "none")
            {
                continue;
            }

            if (debug)
            {
                Info<< "Mapping flux " << iter.key()
                    << " using interpolated flux " << UName
                    << endl;
            }

            surfaceScalarField& phi = *iter();
            surfaceScalarField::GeometricBoundaryField& bphi =
                phi.boundaryField();

            const surfaceScalarField phiU
            (
                fvc::interpolate
                (
                    lookupObject<volVectorField>(UName)
                )
              & Sf()
            );


            forAllConstIter(Map<label>, faceToSplitPoint, iter)
            {
                label oldFaceI = iter.key();
                label oldPointI = iter();

                if (reversePointMap[oldPointI] < 0)
                {
                    // midpoint was removed. See if face still exists.
                    label faceI = reverseFaceMap[oldFaceI];

                    if (faceI >= 0)
                    {
                        if (isInternalFace(faceI))
                        {
                            phi[faceI] = phiU[faceI];
                        }
                        else
                        {
                            label patchI = boundaryMesh().whichPatch(faceI);
                            label i = faceI - boundaryMesh()[patchI].start();

                            const fvsPatchScalarField& patchPhiU =
                                phiU.boundaryField()[patchI];
                            fvsPatchScalarField& patchPhi = bphi[patchI];
                            patchPhi[i] = patchPhiU[i];
                        }
                    }
                }
            }
        }
    }


    // Update numbering of cells/vertices.
    meshCutter_.updateMesh(map);

    // Update numbering of protectedCell_
    if (protectedCell_.size())
    {
        PackedBoolList newProtectedCell(nCells());

        forAll(newProtectedCell, cellI)
        {
            label oldCellI = map().cellMap()[cellI];
            if (oldCellI >= 0)
            {
                newProtectedCell.set(cellI, protectedCell_.get(oldCellI));
            }
        }
        protectedCell_.transfer(newProtectedCell);
    }

    // Debug: Check refinement levels (across faces only)
    meshCutter_.checkRefinementLevels(-1, labelList(0));

    return map;
}


// Get max of connected point
Foam::scalarField
Foam::dynamicRefineBalancedFvMesh::maxPointField(const scalarField& pFld) const
{
    scalarField vFld(nCells(), -GREAT);

    forAll(pointCells(), pointI)
    {
        const labelList& pCells = pointCells()[pointI];

        forAll(pCells, i)
        {
            vFld[pCells[i]] = max(vFld[pCells[i]], pFld[pointI]);
        }
    }
    return vFld;
}


// Get max of connected cell
Foam::scalarField
Foam::dynamicRefineBalancedFvMesh::maxCellField(const volScalarField& vFld) const
{
    scalarField pFld(nPoints(), -GREAT);

    forAll(pointCells(), pointI)
    {
        const labelList& pCells = pointCells()[pointI];

        forAll(pCells, i)
        {
            pFld[pointI] = max(pFld[pointI], vFld[pCells[i]]);
        }
    }
    return pFld;
}


// Simple (non-parallel) interpolation by averaging.
Foam::scalarField
Foam::dynamicRefineBalancedFvMesh::cellToPoint(const scalarField& vFld) const
{
    scalarField pFld(nPoints());

    forAll(pointCells(), pointI)
    {
        const labelList& pCells = pointCells()[pointI];

        scalar sum = 0.0;
        forAll(pCells, i)
        {
            sum += vFld[pCells[i]];
        }
        pFld[pointI] = sum/pCells.size();
    }
    return pFld;
}


// Calculate error. Is < 0 or distance to minLevel, maxLevel
Foam::scalarField Foam::dynamicRefineBalancedFvMesh::error
(
    const scalarField& fld,
    const scalar minLevel,
    const scalar maxLevel
) const
{
    scalarField c(fld.size(), -1);

    forAll(fld, i)
    {
        scalar err = min(fld[i]-minLevel, maxLevel-fld[i]);

        if (err >= 0)
        {
            c[i] = err;
        }
    }
    return c;
}


void Foam::dynamicRefineBalancedFvMesh::selectRefineCandidates
(
    const scalar lowerRefineLevel,
    const scalar upperRefineLevel,
    const scalarField& vFld,
    PackedBoolList& candidateCell
) const
{
    // Get error per cell. Is -1 (not to be refined) to >0 (to be refined,
    // higher more desirable to be refined).
    scalarField cellError
    (
        maxPointField
        (
            error
            (
                cellToPoint(vFld),
                lowerRefineLevel,
                upperRefineLevel
            )
        )
    );

    // Mark cells that are candidates for refinement.
    forAll(cellError, cellI)
    {
        if (cellError[cellI] > 0)
        {
            candidateCell.set(cellI, 1);
        }
    }
}


Foam::labelList Foam::dynamicRefineBalancedFvMesh::selectRefineCells
(
    const label maxCells,
    const label maxRefinement,
    const PackedBoolList& candidateCell
) const
{
    // Every refined cell causes 7 extra cells
    label nTotToRefine = (maxCells - globalData().nTotalCells()) / 7;

    const labelList& cellLevel = meshCutter_.cellLevel();

    // Mark cells that cannot be refined since they would trigger refinement
    // of protected cells (since 2:1 cascade)
    PackedBoolList unrefineableCell;
    calculateProtectedCells(unrefineableCell);

    // Count current selection
    label nLocalCandidates = count(candidateCell, 1);
    label nCandidates = returnReduce(nLocalCandidates, sumOp<label>());

    // Collect all cells
    DynamicList<label> candidates(nLocalCandidates);

    if (nCandidates < nTotToRefine)
    {
        forAll(candidateCell, cellI)
        {
            if
            (
                cellLevel[cellI] < maxRefinement
             && candidateCell.get(cellI)
             && (
                    unrefineableCell.empty()
                 || !unrefineableCell.get(cellI)
                )
            )
            {
                candidates.append(cellI);
            }
        }
    }
    else
    {
        // Sort by error? For now just truncate.
        for (label level = 0; level < maxRefinement; level++)
        {
            forAll(candidateCell, cellI)
            {
                if
                (
                    cellLevel[cellI] == level
                 && candidateCell.get(cellI)
                 && (
                        unrefineableCell.empty()
                     || !unrefineableCell.get(cellI)
                    )
                )
                {
                    candidates.append(cellI);
                }
            }

            if (returnReduce(candidates.size(), sumOp<label>()) > nTotToRefine)
            {
                break;
            }
        }
    }

    // Guarantee 2:1 refinement after refinement
    labelList consistentSet
    (
        meshCutter_.consistentRefinement
        (
            candidates.shrink(),
            true               // Add to set to guarantee 2:1
        )
    );

    Info<< "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
        << " cells for refinement out of " << globalData().nTotalCells()
        << "." << endl;

    return consistentSet;
}


Foam::labelList Foam::dynamicRefineBalancedFvMesh::selectUnrefinePoints
(
    const scalar unrefineLevel,
    const PackedBoolList& markedCell,
    const scalarField& pFld
) const
{
    // All points that can be unrefined
    const labelList splitPoints(meshCutter_.getSplitPoints());

    DynamicList<label> newSplitPoints(splitPoints.size());

    forAll(splitPoints, i)
    {
        label pointI = splitPoints[i];

        if (pFld[pointI] < unrefineLevel)
        {
            // Check that all cells are not marked
            const labelList& pCells = pointCells()[pointI];

            bool hasMarked = false;

            forAll(pCells, pCellI)
            {
                if (markedCell.get(pCells[pCellI]))
                {
                    hasMarked = true;
                    break;
                }
            }

            if (!hasMarked)
            {
                newSplitPoints.append(pointI);
            }
        }
    }


    newSplitPoints.shrink();

    // Guarantee 2:1 refinement after unrefinement
    labelList consistentSet
    (
        meshCutter_.consistentUnrefinement
        (
            newSplitPoints,
            false
        )
    );
    Info<< "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
        << " split points out of a possible "
        << returnReduce(splitPoints.size(), sumOp<label>())
        << "." << endl;

    return consistentSet;
}


void Foam::dynamicRefineBalancedFvMesh::extendMarkedCells
(
    PackedBoolList& markedCell
) const
{
    // Mark faces using any marked cell
    boolList markedFace(nFaces(), false);

    forAll(markedCell, cellI)
    {
        if (markedCell.get(cellI))
        {
            const cell& cFaces = cells()[cellI];

            forAll(cFaces, i)
            {
                markedFace[cFaces[i]] = true;
            }
        }
    }

    syncTools::syncFaceList(*this, markedFace, orEqOp<bool>());

    // Update cells using any markedFace
    for (label faceI = 0; faceI < nInternalFaces(); faceI++)
    {
        if (markedFace[faceI])
        {
            markedCell.set(faceOwner()[faceI], 1);
            markedCell.set(faceNeighbour()[faceI], 1);
        }
    }
    for (label faceI = nInternalFaces(); faceI < nFaces(); faceI++)
    {
        if (markedFace[faceI])
        {
            markedCell.set(faceOwner()[faceI], 1);
        }
    }
}


void Foam::dynamicRefineBalancedFvMesh::checkEightAnchorPoints
(
    PackedBoolList& protectedCell,
    label& nProtected
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const labelList& pointLevel = meshCutter_.pointLevel();

    labelList nAnchorPoints(nCells(), 0);

    forAll(pointLevel, pointI)
    {
        const labelList& pCells = pointCells(pointI);

        forAll(pCells, pCellI)
        {
            label cellI = pCells[pCellI];

            if (pointLevel[pointI] <= cellLevel[cellI])
            {
                // Check if cell has already 8 anchor points -> protect cell
                if (nAnchorPoints[cellI] == 8)
                {
                    if (protectedCell.set(cellI, true))
                    {
                        nProtected++;
                    }
                }

                if (!protectedCell[cellI])
                {
                    nAnchorPoints[cellI]++;
                }
            }
        }
    }


    forAll(protectedCell, cellI)
    {
        if (!protectedCell[cellI] && nAnchorPoints[cellI] != 8)
        {
            protectedCell.set(cellI, true);
            nProtected++;
        }
    }
}


Foam::label Foam::dynamicRefineBalancedFvMesh::topParentID(label p)
{
    label nextP = meshCutter().history().splitCells()[p].parent_.index();
    if( nextP < 0 )
    {
        return p;
    }
    else
    {
        return topParentID(nextP);
    }
}

Foam::Pair<scalar> Foam::dynamicRefineBalancedFvMesh::readRefinementPoints()
{
    dictionary refineDict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).subDict("dynamicRefineFvMeshCoeffs")
    );
    
    return Pair<scalar>
    (
        readScalar(refineDict.lookup("unrefineLevel")), 
        readScalar(refineDict.lookup("lowerRefineLevel"))
    );
}

void Foam::dynamicRefineBalancedFvMesh::updateRefinementField()
{
    Info<< "Calculating internal refinement field" << endl;
    
    volScalarField& intRefFld = *internalRefinementFieldPtr_;
    
    // Set the internal refinement field to zero to start with
    intRefFld = dimensionedScalar("zero",dimless,0.0);
    
    // Get the cell level field from dynamicRefineFvMesh
    const labelList& cellLevel = meshCutter().cellLevel();
    
    // Read the points at which refinement and unrefinement occur from the
    // dynamicMeshDict entries
    Pair<scalar> refinePoints = readRefinementPoints();
    
    // First gradients
    List<word> gradFieldNames = gradFields_.toc();
    
    Field<scalar> cubeRtV = Foam::pow(this->V(),1.0/3.0);
    Field<scalar> refFld(nCells(),0.0);
    
    forAll(gradFieldNames, i)
    {
        word fldName = gradFieldNames[i];
        scalar wgt = gradFields_[fldName].first();
        label maxLevel = static_cast<label>(gradFields_[fldName].second()+0.5);
        
        const volScalarField& fld = this->lookupObject<volScalarField>(fldName);
        
        refFld = wgt * mag(fvc::grad(fld)) * cubeRtV;
        
        // Limit the value of refFld based on its max level
        forAll(refFld, cellI)
        {
            if( cellLevel[cellI] >= maxLevel )
            {
                refFld[cellI] = min
                (
                    refFld[cellI],
                    0.5*(refinePoints.first() + refinePoints.second())
                );
            }
        }
 
        intRefFld.internalField() = max
        (
            intRefFld.internalField(),
            refFld
        );
    }
    
    // Then curls
    List<word> curlFieldNames = curlFields_.toc();
    
    
    forAll(curlFieldNames, i)
    {
        word fldName = curlFieldNames[i];
        scalar wgt = curlFields_[fldName].first();
        label maxLevel = static_cast<label>(curlFields_[fldName].second()+0.5);
        
        const volVectorField& fld = this->lookupObject<volVectorField>(fldName);
        
        refFld = wgt * mag(fvc::curl(fld)) * cubeRtV;
        
        // Limit the value of refFld based on its max level
        forAll(refFld, cellI)
        {
            if( cellLevel[cellI] >= maxLevel )
            {
                refFld[cellI] = min
                (
                    refFld[cellI],
                    0.5*(refinePoints.first() + refinePoints.second())
                );
            }
        }
        
        intRefFld.internalField() = max
        (
            intRefFld.internalField(),
            refFld
        );
    }
    
    // The set refinement physical regions (force the mesh to stay refined
    // near key features)
    forAll(refinedRegions_, regionI)
    {
        const entry& region = refinedRegions_[regionI];
        
        autoPtr<topoSetSource> source = 
            topoSetSource::New(region.keyword(), *this, region.dict());
    
        refFld *= 0.0;
        
        cellSet selectedCellSet
        (
            *this,
            "cellSet",
            nCells()/10+1 //estimate
        );
        
        source->applyToSet
        (
            topoSetSource::NEW,
            selectedCellSet
        );
        
        const labelList cells = selectedCellSet.toc();
        
        label minLevel = readLabel(region.dict().lookup("minLevel"));
        
        forAll(cells, i)
        {
            const label& cellI = cells[i];
            
            if( cellLevel[cellI] < minLevel ) //force refinement
            {
                refFld[cellI] = refinePoints.second() + 1.0;
            }
            else if( cellLevel[cellI] == minLevel ) // keep from coarsening
            {
                refFld[cellI] = 0.5*(refinePoints.first() + refinePoints.second());
            }
            // else do nothing
        }
        
        intRefFld.internalField() = max
        (
            intRefFld.internalField(),
            refFld
        );
    }
    
    intRefFld.correctBoundaryConditions();
    
    Info<<"Min,max refinement field = " << Foam::min(intRefFld).value() << ", "
        << Foam::max(intRefFld).value() << endl;
        
}

Foam::labelList Foam::dynamicRefineBalancedFvMesh::selectProcBoundaryUnrefinePoints
(
    const scalar unrefineLevel,
    const PackedBoolList& markedCell,
    const scalarField& pFld
) const
{
    // All points that can be unrefined
    const labelList splitPoints(meshCutter_.getProcBoundarySplitPoints());

    DynamicList<label> newSplitPoints(splitPoints.size());

    Map<bool> hasMarked;

    forAll(splitPoints, i)
    {
        label pointI = splitPoints[i];

        if (pFld[pointI] < unrefineLevel)
        {
            // Check that all cells are not marked
            const labelList& pCells = pointCells()[pointI];

            if(!hasMarked.found(pointI))
            {
            	hasMarked.insert(pointI, false);
            }

            forAll(pCells, pCellI)
            {
                if (markedCell.get(pCells[pCellI]))
                {
                    hasMarked[pointI] = true;
                    break;
                }
            }
        }
    }

    Pstream::mapCombineGather(hasMarked, orEqOp<bool>());
    Pstream::mapCombineScatter(hasMarked);

    forAll(splitPoints, i)
    {
        label pointI = splitPoints[i];
		if (!hasMarked[pointI])
		{
			newSplitPoints.append(pointI);
		}
    }

    newSplitPoints.shrink();

    // Guarantee 2:1 refinement after unrefinement
    labelList consistentSet
    (
        meshCutter_.consistentUnrefinement
        (
            newSplitPoints,
            false
        )
    );

    return consistentSet;
}

void Foam::dynamicRefineBalancedFvMesh::redistributeUnrefine(labelList& splitPoints)
{
    Map<label> combineCells;
    forAll(splitPoints, i)
    {
        label pointI = splitPoints[i];
        const labelList& pCells = pointCells()[pointI];
        forAll(pCells, j)
        {
        	combineCells[pCells[j]] = meshCutter_.history().parent(pCells[j]).proc();
        }
    }

    Info << "Here1" << endl;
	labelList newDecomp(nCells(),Pstream::myProcNo());
	Map<label>::const_iterator combineIter = combineCells.begin();
	forAll(newDecomp, i)
	{
		if(combineIter == combineCells.end())
			break;

		if(combineIter.key() == i)
		{
			newDecomp[i] = *combineIter;
			combineIter++;
		}
	}
    Info << "Here2" << endl;
    scalar tolDim = globalMeshData::matchTol_ * bounds().mag();

    fvMeshDistribute distributor(*this, tolDim);

    autoPtr<mapDistributePolyMesh> map =
          distributor.distribute( newDecomp );
    Info << "Here3" << endl;

    meshCutter_.distribute(map);

    Info << "Here4" << endl;

    //Correct values on all cyclic patches
    correctBoundaries<scalar>();
    correctBoundaries<vector>();
    correctBoundaries<sphericalTensor>();
    correctBoundaries<symmTensor>();
    correctBoundaries<tensor>();
}

void Foam::dynamicRefineBalancedFvMesh::readRefinementDict()
{
    dictionary dynamicMeshDict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        )
    );
    
    if( dynamicMeshDict.isDict("refinementControls") )
    {
        dictionary refineControlDict = 
            dynamicMeshDict.subDict("refinementControls");
        
        enableRefinementControl_ = 
            Switch(refineControlDict.lookup("enableRefinementControl"));
        
        if( enableRefinementControl_ )
        {
            // Overwrite field name entry in dynamicRefineBalancedFvMeshCoeffs?
            // For now you just have to be smart and enter 
            // 'internalRefinementField' for the field name manually
            
            // Read HashTable of gradient-refinement scalars
            if( refineControlDict.found("gradients") )
            {
                gradFields_ = HashTable< Pair<scalar> >
                (
                    refineControlDict.lookup("gradients")
                );
            }
            
            // Read HashTable of curl-refinement vectors
            if( refineControlDict.found("curls") )
            {
                curlFields_ = HashTable< Pair<scalar> >
                (
                    refineControlDict.lookup("curls")
                );
            }
            
            // Read refinement regions
            if( refineControlDict.found("regions") )
            {
                refinedRegions_ = PtrList<entry>
                (
                    refineControlDict.lookup("regions")
                );
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
            
Foam::dynamicRefineBalancedFvMesh::dynamicRefineBalancedFvMesh
(
    const IOobject& io
)
:
	dynamicFvMesh(io),
	meshCutter_(*this),
	dumpLevel_(false),
	nRefinementIterations_(0),
	protectedCell_(nCells(),0),
    internalRefinementFieldPtr_(NULL),
    gradFields_(),
    curlFields_(),
    refinedRegions_(),
    enableRefinementControl_(false)
{
    setupDynamicRefineFvMesh(io);
    readRefinementDict();
    
    if( enableRefinementControl_ )
    {
        internalRefinementFieldPtr_ = new volScalarField
        (
            IOobject
            (
                "internalRefinementField",
                this->time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            *this,
            dimensionedScalar("zero", dimless, 0.0)
        );
    }
}

void Foam::dynamicRefineBalancedFvMesh::setupDynamicRefineFvMesh(const IOobject& io)
{
	readDict();

    const labelList& cellLevel = meshCutter_.cellLevel();
    const labelList& pointLevel = meshCutter_.pointLevel();

    // Set cells that should not be refined.
    // This is currently any cell which does not have 8 anchor points or
    // uses any face which does not have 4 anchor points.
    // Note: do not use cellPoint addressing

    // Count number of points <= cellLevel
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList nAnchors(nCells(), 0);

    label nProtected = 0;

    forAll(pointCells(), pointI)
    {
        const labelList& pCells = pointCells()[pointI];

        forAll(pCells, i)
        {
            label cellI = pCells[i];

            if (!protectedCell_.get(cellI))
            {
                if (pointLevel[pointI] <= cellLevel[cellI])
                {
                    nAnchors[cellI]++;

                    if (nAnchors[cellI] > 8)
                    {
                        protectedCell_.set(cellI, 1);
                        nProtected++;
                    }
                }
            }
        }
    }


    // Count number of points <= faceLevel
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Bit tricky since proc face might be one more refined than the owner since
    // the coupled one is refined.

    {
        labelList neiLevel(nFaces());

        for (label faceI = 0; faceI < nInternalFaces(); faceI++)
        {
            neiLevel[faceI] = cellLevel[faceNeighbour()[faceI]];
        }
        for (label faceI = nInternalFaces(); faceI < nFaces(); faceI++)
        {
            neiLevel[faceI] = cellLevel[faceOwner()[faceI]];
        }
        syncTools::swapFaceList(*this, neiLevel);


        boolList protectedFace(nFaces(), false);

        forAll(faceOwner(), faceI)
        {
            label faceLevel = max
            (
                cellLevel[faceOwner()[faceI]],
                neiLevel[faceI]
            );

            const face& f = faces()[faceI];

            label nAnchors = 0;

            forAll(f, fp)
            {
                if (pointLevel[f[fp]] <= faceLevel)
                {
                    nAnchors++;

                    if (nAnchors > 4)
                    {
                        protectedFace[faceI] = true;
                        break;
                    }
                }
            }
        }

        syncTools::syncFaceList(*this, protectedFace, orEqOp<bool>());

        for (label faceI = 0; faceI < nInternalFaces(); faceI++)
        {
            if (protectedFace[faceI])
            {
                protectedCell_.set(faceOwner()[faceI], 1);
                nProtected++;
                protectedCell_.set(faceNeighbour()[faceI], 1);
                nProtected++;
            }
        }
        for (label faceI = nInternalFaces(); faceI < nFaces(); faceI++)
        {
            if (protectedFace[faceI])
            {
                protectedCell_.set(faceOwner()[faceI], 1);
                nProtected++;
            }
        }

        // Also protect any cells that are less than hex
        forAll(cells(), cellI)
        {
            const cell& cFaces = cells()[cellI];

            if (cFaces.size() < 6)
            {
                if (protectedCell_.set(cellI, 1))
                {
                    nProtected++;
                }
            }
            else
            {
                forAll(cFaces, cFaceI)
                {
                    if (faces()[cFaces[cFaceI]].size() < 4)
                    {
                        if (protectedCell_.set(cellI, 1))
                        {
                            nProtected++;
                        }
                        break;
                    }
                }
            }
        }

        // Check cells for 8 corner points
        checkEightAnchorPoints(protectedCell_, nProtected);
    }

    if (returnReduce(nProtected, sumOp<label>()) == 0)
    {
        protectedCell_.clear();
    }
    else
    {

        cellSet protectedCells(*this, "protectedCells", nProtected);
        forAll(protectedCell_, cellI)
        {
            if (protectedCell_[cellI])
            {
                protectedCells.insert(cellI);
            }
        }

        Info<< "Detected " << returnReduce(nProtected, sumOp<label>())
            << " cells that are protected from refinement."
            << " Writing these to cellSet "
            << protectedCells.name()
            << "." << endl;

        protectedCells.write();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicRefineBalancedFvMesh::~dynamicRefineBalancedFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicRefineBalancedFvMesh::updateA()
{
    // Re-read dictionary. Choosen since usually -small so trivial amount
    // of time compared to actual refinement. Also very useful to be able
    // to modify on-the-fly.
    dictionary refineDict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).subDict("dynamicRefineFvMeshCoeffs")
    );

    label refineInterval = readLabel(refineDict.lookup("refineInterval"));

    bool hasChanged = false;

    if (refineInterval == 0)
    {
        topoChanging(hasChanged);

        return false;
    }
    else if (refineInterval < 0)
    {
        FatalErrorIn("dynamicRefineFvMesh::update()")
            << "Illegal refineInterval " << refineInterval << nl
            << "The refineInterval setting in the dynamicMeshDict should"
            << " be >= 1." << nl
            << exit(FatalError);
    }




    // Note: cannot refine at time 0 since no V0 present since mesh not
    //       moved yet.

    if (time().timeIndex() > 0 && time().timeIndex() % refineInterval == 0)
    {
        label maxCells = readLabel(refineDict.lookup("maxCells"));

        if (maxCells <= 0)
        {
            FatalErrorIn("dynamicRefineFvMesh::update()")
                << "Illegal maximum number of cells " << maxCells << nl
                << "The maxCells setting in the dynamicMeshDict should"
                << " be > 0." << nl
                << exit(FatalError);
        }

        label maxRefinement = readLabel(refineDict.lookup("maxRefinement"));

        if (maxRefinement <= 0)
        {
            FatalErrorIn("dynamicRefineFvMesh::update()")
                << "Illegal maximum refinement level " << maxRefinement << nl
                << "The maxCells setting in the dynamicMeshDict should"
                << " be > 0." << nl
                << exit(FatalError);
        }

        const word fieldName(refineDict.lookup("field"));

        const volScalarField& vFld = lookupObject<volScalarField>(fieldName);

        const scalar lowerRefineLevel =
            readScalar(refineDict.lookup("lowerRefineLevel"));
        const scalar upperRefineLevel =
            readScalar(refineDict.lookup("upperRefineLevel"));
        const scalar unrefineLevel = refineDict.lookupOrDefault<scalar>
        (
            "unrefineLevel",
            GREAT
        );
        const label nBufferLayers =
            readLabel(refineDict.lookup("nBufferLayers"));

        // Cells marked for refinement or otherwise protected from unrefinement.
        PackedBoolList refineCell(nCells());

        // Determine candidates for refinement (looking at field only)
        selectRefineCandidates
        (
            lowerRefineLevel,
            upperRefineLevel,
            vFld,
            refineCell
        );

        if (globalData().nTotalCells() < maxCells)
        {
            // Select subset of candidates. Take into account max allowable
            // cells, refinement level, protected cells.
            labelList cellsToRefine
            (
                selectRefineCells
                (
                    maxCells,
                    maxRefinement,
                    refineCell
                )
            );

            label nCellsToRefine = returnReduce
            (
                cellsToRefine.size(), sumOp<label>()
            );

            if (nCellsToRefine > 0)
            {
                // Refine/update mesh and map fields
                autoPtr<mapPolyMesh> map = refine(cellsToRefine);

                // Update refineCell. Note that some of the marked ones have
                // not been refined due to constraints.
                {
                    const labelList& cellMap = map().cellMap();
                    const labelList& reverseCellMap = map().reverseCellMap();

                    PackedBoolList newRefineCell(cellMap.size());

                    forAll(cellMap, cellI)
                    {
                        label oldCellI = cellMap[cellI];

                        if (oldCellI < 0)
                        {
                            newRefineCell.set(cellI, 1);
                        }
                        else if (reverseCellMap[oldCellI] != cellI)
                        {
                            newRefineCell.set(cellI, 1);
                        }
                        else
                        {
                            newRefineCell.set(cellI, refineCell.get(oldCellI));
                        }
                    }
                    refineCell.transfer(newRefineCell);
                }

                // Extend with a buffer layer to prevent neighbouring points
                // being unrefined.
                for (label i = 0; i < nBufferLayers; i++)
                {
                    extendMarkedCells(refineCell);
                }

                hasChanged = true;
            }
        }


        {
            // Select unrefineable points that are not marked in refineCell
            labelList procBoundaryUnrefinePoints
            (
                selectProcBoundaryUnrefinePoints
                (
                    unrefineLevel,
                    refineCell,
                    maxCellField(vFld)
                )
            );

            if(!procBoundaryUnrefinePoints.empty())
            	redistributeUnrefine(procBoundaryUnrefinePoints);

            labelList pointsToUnrefine
            (
                selectUnrefinePoints
                (
                    unrefineLevel,
                    refineCell,
                    maxCellField(vFld)
                )
            );

            label nSplitPoints = returnReduce
            (
                pointsToUnrefine.size(),
                sumOp<label>()
            );

            if (nSplitPoints > 0)
            {
                // Refine/update mesh
                unrefine(pointsToUnrefine);

                hasChanged = true;
            }
        }


        if ((nRefinementIterations_ % 10) == 0)
        {
            // Compact refinement history occassionally (how often?).
            // Unrefinement causes holes in the refinementHistory.
            const_cast<refinementHistoryBalanced&>(meshCutter().history()).compact();
        }
        nRefinementIterations_++;
    }

    topoChanging(hasChanged);
    if (hasChanged)
    {
        // Reset moving flag (if any). If not using inflation we'll not move,
        // if are using inflation any follow on movePoints will set it.
        moving(false);
    }

    return hasChanged;
}

bool Foam::dynamicRefineBalancedFvMesh::update()
{
    //Part 0 - Update internally calculated refinement field
    readRefinementDict();
    
    if( enableRefinementControl_ )
    {
        updateRefinementField();
    }
    
    //Part 1 - Call normal update from dynamicRefineFvMesh
    Info << "Running mesh update" << endl;
    bool hasChanged = updateA();
    
    // Part 2 - Load Balancing
//    Info << "Running mesh load balance" << endl;
    dictionary refineDict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).subDict("dynamicRefineFvMeshCoeffs")
    );
    
    Switch enableBalancing = refineDict.lookup("enableBalancing");
    
    if ( Pstream::parRun() && hasChanged )
    {
        const scalar allowableImbalance = 
            readScalar(refineDict.lookup("allowableImbalance"));
            
        //First determine current level of imbalance - do this for all
        // parallel runs with a changing mesh, even if balancing is disabled
        label nGlobalCells = globalData().nTotalCells();
        scalar idealNCells = scalar(nGlobalCells)/scalar(Pstream::nProcs());
        scalar localImbalance = mag(scalar(nCells()) - idealNCells);
        Foam::reduce(localImbalance, maxOp<scalar>());
        scalar maxImbalance = localImbalance/idealNCells;
        
        Info<<"Maximum imbalance = " << 100*maxImbalance << " %" << endl;
        
        //If imbalanced, construct weighted coarse graph (level 0) with node
        // weights equal to their number of subcells. This partitioning works
        // as long as the number of level 0 cells is several times greater than
        // the number of processors.
        if( maxImbalance > allowableImbalance && enableBalancing)
        {
            Info<< "Re-balancing dynamically refined mesh" << endl;
                        
//            const labelIOList& cellLevel = meshCutter().cellLevel();
//            Map<label> coarseIDmap(100);
//            labelList uniqueIndex(nCells(),0);
//
//            label nCoarse = 0;
//
//            forAll(cells(), cellI)
//            {
//                if( cellLevel[cellI] > 0 )
//                {
//                    uniqueIndex[cellI] = nCells() + topParentID
//                    (
//                        meshCutter().history().parent(cellI).index()
//                    );
//                }
//                else
//                {
//                    uniqueIndex[cellI] = cellI;
//                }
//
//                if( coarseIDmap.insert(uniqueIndex[cellI], nCoarse) )
//                {
//                    ++nCoarse;
//                }
//            }
            
//            // Convert to local sequential indexing and calculate coarse
//            // points and weights
//            labelList localIndex(nCells(),0);
//            pointField coarsePoints(nCoarse,vector::zero);
//            scalarField coarseWeights(nCoarse,0.0);
//
//            forAll(uniqueIndex, cellI)
//            {
//                localIndex[cellI] = coarseIDmap[uniqueIndex[cellI]];
//
//                // If 2D refinement (quadtree) is ever implemented, this '3'
//                // should be set in general as the number of refinement
//                // dimensions.
//                label w = (1 << (3*cellLevel[cellI]));
//
//                coarseWeights[localIndex[cellI]] += 1.0;
//                coarsePoints[localIndex[cellI]] += C()[cellI]/w;
//            }

//            Info << "max weight = " << gMax(coarseWeights) << endl;
            
            //Set up decomposer - a separate dictionary is used here so
            // you can use a simple partitioning for decomposePar and
            // ptscotch for the rebalancing (or any chosen algorithms)
            autoPtr<decompositionMethod> decomposer
            (
                decompositionMethod::New
                (
                    IOdictionary
                    (
                        IOobject
                        (
                            "balanceParDict",
                            time().system(),
                            *this,
                            IOobject::MUST_READ_IF_MODIFIED,
                            IOobject::NO_WRITE
                        )
                    )
                )
            );

            //Ideally we should use no agglomeration here
//            labelList finalDecomp = decomposer().decompose
//            (
//                *this,
//                localIndex, //agglomeration vector
//                coarsePoints,
//                coarseWeights
//            );
            labelList finalDecomp = decomposer().decompose
            (
                *this,
                C()
            );

//            List<int> cellCount(Pstream::nProcs(),0);
//            forAll(finalDecomp, cellI)
//            {
//                cellCount[finalDecomp[cellI]]++;
//            }
//            int size = finalDecomp.size();
//            for(size_t i = 0; i < Pstream::nProcs(); i++)
//                sumReduce(cellCount[i],size);
//            Info << "cellCount = " << cellCount << endl;

            scalar tolDim = globalMeshData::matchTol_ * bounds().mag();
            
            
            fvMeshDistribute distributor(*this, tolDim);
            
            autoPtr<mapDistributePolyMesh> map =
                  distributor.distribute(finalDecomp);
                  
//            Info << "Debug point at " << __FILE__ << ':' << __LINE__ << endl;

            meshCutter_.distribute(map);
            
//            Info << "Debug point at " << __FILE__ << ':' << __LINE__ << endl;

            //Correct values on all cyclic patches
            correctBoundaries<scalar>();
            correctBoundaries<vector>();
            correctBoundaries<sphericalTensor>();
            correctBoundaries<symmTensor>();
            correctBoundaries<tensor>();
        }
    }

    return hasChanged;
}

bool Foam::dynamicRefineBalancedFvMesh::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    // Force refinement data to go to the current time directory.
    const_cast<hexRef8Balanced&>(meshCutter_).setInstance(time().timeName());

    bool writeOk =
    (
        dynamicFvMesh::writeObjects(fmt, ver, cmp)
     && meshCutter_.write()
    );

    if (dumpLevel_)
    {
        volScalarField scalarCellLevel
        (
            IOobject
            (
                "cellLevel",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            *this,
            dimensionedScalar("level", dimless, 0)
        );

        const labelList& cellLevel = meshCutter_.cellLevel();

        forAll(cellLevel, cellI)
        {
            scalarCellLevel[cellI] = cellLevel[cellI];
        }

        writeOk = writeOk && scalarCellLevel.write();
    }

    return writeOk;
}


// ************************************************************************* //
