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


Class
    Foam::dynamicRefineBalancedFvMesh

Description
    A fvMesh with built-in refinement and run-time load balancing.

SourceFiles
    dynamicRefineBalancedFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicRefineBalancedFvMesh_H
#define dynamicRefineBalancedFvMesh_H

#include <fstream>

#include "dynamicFvMesh.H"
#include "dynamicRefineFvMesh.H"
#include "hexRef8Balanced.H"
#include "PackedBoolList.H"
#include "Switch.H"
#include "decompositionMethod.H"
#include "fvMeshDistribute.H"
#include "mapDistributePolyMesh.H"
#include "HashTable.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "PtrDictionary.H"
#include "dictionaryEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dynamicRefineBalancedFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicRefineBalancedFvMesh
:
    public dynamicFvMesh
{
protected:

        //- Mesh cutting engine
        hexRef8Balanced meshCutter_;

        //- Dump cellLevel for postprocessing
        Switch dumpLevel_;

        //- Fluxes to map
        HashTable<word> correctFluxes_;

        //- Number of refinement/unrefinement steps done so far.
        label nRefinementIterations_;

        //- Protected cells (usually since not hexes)
        PackedBoolList protectedCell_;


    // Private Member Functions

        //- Count set/unset elements in packedlist.
        static label count(const PackedBoolList&, const unsigned int);

        //- Calculate cells that cannot be refined since would trigger
        //  refinement of protectedCell_ (since 2:1 refinement cascade)
        void calculateProtectedCells(PackedBoolList& unrefineableCell) const;

        //- Read the projection parameters from dictionary
        void readDict();


        //- Refine cells. Update mesh and fields.
        autoPtr<mapPolyMesh> refine(const labelList&);

        //- Unrefine cells. Gets passed in centre points of cells to combine.
        autoPtr<mapPolyMesh> unrefine(const labelList&);

        //- Redistribute cells to their parent processors before unrefining
        void redistributeUnrefine(labelList& splitPoints);


        // Selection of cells to un/refine

            //- Calculates approximate value for refinement level so
            //  we don't go above maxCell
            scalar getRefineLevel
            (
                const label maxCells,
                const label maxRefinement,
                const scalar refineLevel,
                const scalarField&
            ) const;

            //- Get per cell max of connected point
            scalarField maxPointField(const scalarField&) const;

            //- Get point max of connected cell
            scalarField maxCellField(const volScalarField&) const;

            scalarField cellToPoint(const scalarField& vFld) const;

            scalarField error
            (
                const scalarField& fld,
                const scalar minLevel,
                const scalar maxLevel
            ) const;

            //- Select candidate cells for refinement
            virtual void selectRefineCandidates
            (
                const scalar lowerRefineLevel,
                const scalar upperRefineLevel,
                const scalarField& vFld,
                PackedBoolList& candidateCell
            ) const;

            //- Subset candidate cells for refinement
            virtual labelList selectRefineCells
            (
                const label maxCells,
                const label maxRefinement,
                const PackedBoolList& candidateCell
            ) const;

            //- Select points that can be unrefined.
            virtual labelList selectUnrefinePoints
            (
                const scalar unrefineLevel,
                const PackedBoolList& markedCell,
                const scalarField& pFld
            ) const;

            //- Select points that can be unrefined.
            virtual labelList selectProcBoundaryUnrefinePoints
            (
                const scalar unrefineLevel,
                const PackedBoolList& markedCell,
                const scalarField& pFld
            ) const;

            //- Extend markedCell with cell-face-cell.
            void extendMarkedCells(PackedBoolList& markedCell) const;

            //- Check all cells have 8 anchor points
            void checkEightAnchorPoints
            (
                PackedBoolList& protectedCell,
                label& nProtected
            ) const;

private:
        //-
        volScalarField* internalRefinementFieldPtr_;
        
        //-
        HashTable< Pair<scalar> > gradFields_;
        
        //-
        HashTable< Pair<scalar> > curlFields_;
        
        //-
        PtrList<entry> refinedRegions_; 
        
        //-
        Switch enableRefinementControl_;
        
        //-
        void updateRefinementField();
        
        //-
        void readRefinementDict();
        
        //-
        Pair<scalar> readRefinementPoints();
        
        //- Disallow default bitwise copy construct
        dynamicRefineBalancedFvMesh(const dynamicRefineBalancedFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicRefineBalancedFvMesh&);

        label topParentID(label p);
public:

    //- Runtime type information
    TypeName("dynamicRefineBalancedFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit dynamicRefineBalancedFvMesh(const IOobject& io);

        void setupDynamicRefineFvMesh(const IOobject& io);


    //- Destructor
    virtual ~dynamicRefineBalancedFvMesh();


    // Member Functions

        //- Direct access to the refinement engine
        const hexRef8Balanced& meshCutter() const
        {
            return meshCutter_;
        }

        //- Cells which should not be refined/unrefined
        const PackedBoolList& protectedCell() const
        {
            return protectedCell_;
        }

        //- Cells which should not be refined/unrefined
        PackedBoolList& protectedCell()
        {
            return protectedCell_;
        }

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
        virtual bool updateA();
        

    // Writing

        //- Write using given format, version and compression
        virtual bool writeObject
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp
        ) const;

        //- Template to update all volField boundaries
        template<class Type> void correctBoundaries();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dynamicRefineBalancedFvMeshTemplates.C"
#endif

#endif

// ************************************************************************* //
