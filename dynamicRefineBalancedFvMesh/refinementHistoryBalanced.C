/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "DynamicList.H"
#include "refinementHistoryBalanced.H"
#include "ListOps.H"
#include "mapPolyMesh.H"
#include "mapDistributePolyMesh.H"
#include "polyMesh.H"
#include "IPstream.H"
#include "OPstream.H"
#include "PstreamReduceOps.H"
//#include <fstream> //temporary for debugging
//#include "Time.H" //temporary for debugging


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(refinementHistoryBalanced, 0);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::refinementHistoryBalanced::writeEntry
(
    const List<splitCell8>& splitCells,
    const splitCell8& split
)
{
    // Write me:
    if (split.addedCellsPtr_.valid())
    {
        Pout<< "parentProc:" << split.parent_.proc()
            << "parentIndex:" << split.parent_.index()
            << " subCells:" << split.addedCellsPtr_()
            << endl;
    }
    else
    {
        Pout<< "parentProc:" << split.parent_.proc()
            << "parentIndex:" << split.parent_.index()
            << " no subcells"
            << endl;
    }

    if (split.parent_.index() >= 0)
    {
        if (split.parent_.proc() == Pstream::myProcNo())
        {
            Pout<< "parent data:" << endl;
            // Write my parent
            string oldPrefix = Pout.prefix();
            Pout.prefix() = "  " + oldPrefix;
            writeEntry(splitCells, splitCells[split.parent_.index()]);
            Pout.prefix() = oldPrefix;
        }
        else
        {
            //TODO: use gather/scatter to complete the tree
        }
    }
    
}


void Foam::refinementHistoryBalanced::writeDebug
(
    const labelList& visibleCells,
    const List<splitCell8>& splitCells
)
{
    string oldPrefix = Pout.prefix();
    Pout.prefix() = "";

    forAll(visibleCells, cellI)
    {
        label index = visibleCells[cellI];

        if (index >= 0)
        {
            Pout<< "Cell from refinement:" << cellI << " index:" << index
                << endl;

            string oldPrefix = Pout.prefix();
            Pout.prefix() = "  " + oldPrefix;
            writeEntry(splitCells, splitCells[index]);
            Pout.prefix() = oldPrefix;
        }
        else
        {
            Pout<< "Unrefined cell:" << cellI << " index:" << index << endl;
        }
    }
    Pout.prefix() = oldPrefix;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementHistoryBalanced::procIndexPair::procIndexPair()
:
	pair_(Pstream::myProcNo(),-1)
{}

Foam::refinementHistoryBalanced::procIndexPair::procIndexPair(label proc, label index)
:
	pair_(proc,index)
{}

Foam::refinementHistoryBalanced::procIndexPair::procIndexPair(Istream& is)
:
	pair_(is)
{}

//- Construct null
Foam::refinementHistoryBalanced::splitCell8::splitCell8()
:
    parent_(Pstream::myProcNo(),-1),
    addedCellsPtr_(NULL)
{}


//- Construct as child element of parent
Foam::refinementHistoryBalanced::splitCell8::splitCell8(const label parent)
:
    parent_(Pstream::myProcNo(),parent),
    addedCellsPtr_(NULL)
{}

//- Construct as child element of parent
Foam::refinementHistoryBalanced::splitCell8::splitCell8(const label parentProc, const label parent)
:
    parent_(parentProc,parent),
    addedCellsPtr_(NULL)
{}


//- Construct from Istream
Foam::refinementHistoryBalanced::splitCell8::splitCell8(Istream& is)
{
    is >> *this;
}


//- Construct as (deep) copy.
Foam::refinementHistoryBalanced::splitCell8::splitCell8(const splitCell8& sc)
:
    parent_(sc.parent_),
    addedCellsPtr_
    (
        sc.addedCellsPtr_.valid()
      ? new FixedList<procIndexPair, 8>(sc.addedCellsPtr_())
      : NULL
    )
{}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, refinementHistoryBalanced::splitCell8& sc)
{
    List<refinementHistoryBalanced::procIndexPair> addedCells;

    is >> sc.parent_ >> addedCells;

    if (addedCells.size())
    {
        sc.addedCellsPtr_.reset(new FixedList<refinementHistoryBalanced::procIndexPair, 8>(addedCells));
    }
    else
    {
        sc.addedCellsPtr_.reset(NULL);
    }

    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const refinementHistoryBalanced::splitCell8& sc
)
{
    // Output as labelList so we can have 0 sized lists. Alternative is to
    // output as fixedlist with e.g. -1 elements and check for this upon
    // reading. However would cause much more data to be transferred.

    if (sc.addedCellsPtr_.valid())
    {
        return os
            << sc.parent_
            << token::SPACE
            << List<refinementHistoryBalanced::procIndexPair>(sc.addedCellsPtr_());
    }
    else
    {
        return os << sc.parent_ << token::SPACE
        		<< List<refinementHistoryBalanced::procIndexPair>(0);
    }
}

Foam::Istream& Foam::operator>>(Istream& is, refinementHistoryBalanced::procIndexPair& pair)
{
    return is >> pair.pair_;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const refinementHistoryBalanced::procIndexPair& pair)
{
    return os << pair.pair_;
}

void Foam::refinementHistoryBalanced::checkIndices() const
{
    // Check indices.
    forAll(visibleCells_, i)
    {
        if (visibleCells_[i] < 0 && visibleCells_[i] >= splitCells_.size())
        {
            FatalErrorIn("refinementHistoryBalanced::checkIndices() const")
                << "Illegal entry " << visibleCells_[i]
                << " in visibleCells at location" << i << nl
                << "It points outside the range of splitCells : 0.."
                << splitCells_.size()-1
                << abort(FatalError);
        }
    }
}


Foam::label Foam::refinementHistoryBalanced::allocateSplitCell
(
    const label parent,
    const label i
)
{
    label index = -1;

    if (freeSplitCells_.size())
    {
        index = freeSplitCells_.remove();

        splitCells_[index] = splitCell8(parent);
    }
    else
    {
        index = splitCells_.size();

        splitCells_.append(splitCell8(parent));
    }


    // Update the parent field
    if (parent >= 0)
    {
        splitCell8& parentSplit = splitCells_[parent];

        if (parentSplit.addedCellsPtr_.empty())
        {
            // Allocate storage on parent for the 8 subcells.
            parentSplit.addedCellsPtr_.reset(new FixedList<procIndexPair, 8>(procIndexPair(-1,-1)));
        }


        // Store me on my parent
        FixedList<procIndexPair, 8>& parentSplits = parentSplit.addedCellsPtr_();

        parentSplits[i] = procIndexPair(Pstream::myProcNo(),index);
    }

    return index;
}


void Foam::refinementHistoryBalanced::freeSplitCell(const label index)
{
    splitCell8& split = splitCells_[index];

    // Make sure parent does not point to me anymore.
    if (split.parent_.index() >= 0)
    {
    	if(split.parent_.proc() != Pstream::myProcNo())
    	{
            FatalErrorIn("refinementHistoryBalanced::freeSplitCell")
                << "Problem: parent is not on same process as children"
				<< abort(FatalError);
    	}

        autoPtr<FixedList<procIndexPair, 8> >& subCellsPtr =
            splitCells_[split.parent_.index()].addedCellsPtr_;

        if (subCellsPtr.valid())
        {
            FixedList<procIndexPair, 8>& subCells = subCellsPtr();

            label myPos = findIndex(subCells, procIndexPair(Pstream::myProcNo(),index));

            if (myPos == -1)
            {
                FatalErrorIn("refinementHistoryBalanced::freeSplitCell")
                    << "Problem: cannot find myself in"
                    << " parents' children" << abort(FatalError);
            }
            else
            {
                subCells[myPos] = procIndexPair(-1,-1);
            }
        }
    }

    // Mark splitCell as free
    split.parent_.index()= -2;

    // Add to cache of free splitCells
    freeSplitCells_.append(index);
}

// Mark entry in splitCells. Recursively mark its parent and subs.
void Foam::refinementHistoryBalanced::markSplit
(
    const label index,
    labelList& oldToNew,
    DynamicList<splitCell8>& newSplitCells
) const
{
    if (oldToNew[index] == -1)
    {
        // Not yet compacted.

        const splitCell8& split = splitCells_[index];

        oldToNew[index] = newSplitCells.size();
        newSplitCells.append(split);

        if (split.parent_.index() >= 0)
        {
            markSplit(split.parent_.index(), oldToNew, newSplitCells);
        }
        if (split.addedCellsPtr_.valid())
        {
            const FixedList<procIndexPair, 8>& splits = split.addedCellsPtr_();

            forAll(splits, i)
            {
                if (splits[i].index() >= 0)
                {
                    markSplit(splits[i].index(), oldToNew, newSplitCells);
                }
            }
        }
    }
}

void Foam::refinementHistoryBalanced::parallelMarkSplits
(
	boolListList& marked,
    labelList& oldToNew,
    DynamicList<splitCell8>& newSplitCells
) const
{
    //Loop until marking is complete
    bool complete = false;
    while(!complete)
    {
    	complete = true;
    	forAll (marked[Pstream::myProcNo()], index)
		{
    		if (marked[Pstream::myProcNo()][index] && oldToNew[index] == -1)
    		{
    	        // Not yet compacted.
    			complete = false;

    	        const splitCell8& split = splitCells_[index];

    	        oldToNew[index] = newSplitCells.size();
    	        newSplitCells.append(split);

    	        if (split.parent_.index() >= 0)
    	        {
    	            //markSplit(split.parent_, oldToNew, newSplitCells);
    	        	marked[split.parent_.proc()][split.parent_.index()] = true;
    	        }
    	        if (split.addedCellsPtr_.valid())
    	        {
    	            const FixedList<procIndexPair, 8>& splits = split.addedCellsPtr_();

    	            forAll(splits, i)
    	            {
    	                if (splits[i].index() >= 0)
    	                {
    	                    //markSplit(splits[i], oldToNew, newSplitCells);
    	                	marked[splits[i].proc()][splits[i].index()] = true;
    	                }
    	            }
    	        }
    		}
    	}

    	forAll(marked, procI)
    	{
    	    Pstream::listCombineGather(marked[procI],orEqOp<bool>());
    	    Pstream::listCombineScatter(marked[procI]);
    	}

    	reduce(complete,andOp<bool>());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementHistoryBalanced::refinementHistoryBalanced(const IOobject& io)
:
    regIOobject(io)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningIn
        (
            "refinementHistoryBalanced::refinementHistoryBalanced(const IOobject&)"
        )   << "Specified IOobject::MUST_READ_IF_MODIFIED but class"
            << " does not support automatic rereading."
            << endl;
    }

    if
    (
        io.readOpt() == IOobject::MUST_READ
     || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }

    if (debug)
    {
        Pout<< "refinementHistoryBalanced::refinementHistoryBalanced :"
            << " constructed history from IOobject :"
            << " splitCells:" << splitCells_.size()
            << " visibleCells:" << visibleCells_.size()
            << endl;
    }
}


//- Read or construct
Foam::refinementHistoryBalanced::refinementHistoryBalanced
(
    const IOobject& io,
    const List<splitCell8>& splitCells,
    const labelList& visibleCells
)
:
    regIOobject(io),
    splitCells_(splitCells),
    freeSplitCells_(0),
    visibleCells_(visibleCells)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningIn
        (
            "refinementHistoryBalanced::refinementHistoryBalanced"
            "(const IOobject&, const List<splitCell8>&, const labelList&)"
        )   << "Specified IOobject::MUST_READ_IF_MODIFIED but class"
            << " does not support automatic rereading."
            << endl;
    }

    if
    (
        io.readOpt() == IOobject::MUST_READ
     || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }

    // Check indices.
    checkIndices();

    if (debug)
    {
        Pout<< "refinementHistoryBalanced::refinementHistoryBalanced :"
            << " constructed history from IOobject or components :"
            << " splitCells:" << splitCells_.size()
            << " visibleCells:" << visibleCells_.size()
            << endl;
    }
}


// Construct from initial number of cells (all visible)
Foam::refinementHistoryBalanced::refinementHistoryBalanced
(
    const IOobject& io,
    const label nCells
)
:
    regIOobject(io),
    freeSplitCells_(0)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningIn
        (
            "refinementHistoryBalanced::refinementHistoryBalanced"
            "(const IOobject&, const label)"
        )   << "Specified IOobject::MUST_READ_IF_MODIFIED but class"
            << " does not support automatic rereading."
            << endl;
    }

    if
    (
        io.readOpt() == IOobject::MUST_READ
     || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }
    else
    {
        visibleCells_.setSize(nCells);
        splitCells_.setCapacity(nCells);

        for (label cellI = 0; cellI < nCells; cellI++)
        {
            visibleCells_[cellI] = cellI;
            splitCells_.append(splitCell8());
        }
    }

    // Check indices.
    checkIndices();

    if (debug)
    {
        Pout<< "refinementHistoryBalanced::refinementHistoryBalanced :"
            << " constructed history from IOobject or initial size :"
            << " splitCells:" << splitCells_.size()
            << " visibleCells:" << visibleCells_.size()
            << endl;
    }
}


// Construct as copy
Foam::refinementHistoryBalanced::refinementHistoryBalanced
(
    const IOobject& io,
    const refinementHistoryBalanced& rh
)
:
    regIOobject(io),
    splitCells_(rh.splitCells()),
    freeSplitCells_(rh.freeSplitCells()),
    visibleCells_(rh.visibleCells())
{
    if (debug)
    {
        Pout<< "refinementHistoryBalanced::refinementHistoryBalanced : constructed initial"
            << " history." << endl;
    }
}


// Construct from Istream
Foam::refinementHistoryBalanced::refinementHistoryBalanced(const IOobject& io, Istream& is)
:
    regIOobject(io),
    splitCells_(is),
    freeSplitCells_(0),
    visibleCells_(is)
{
    // Check indices.
    checkIndices();

    if (debug)
    {
        Pout<< "refinementHistoryBalanced::refinementHistoryBalanced :"
            << " constructed history from Istream"
            << " splitCells:" << splitCells_.size()
            << " visibleCells:" << visibleCells_.size()
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::refinementHistoryBalanced::resize(const label size)
{
    label oldSize = visibleCells_.size();

    if (debug)
    {
        Pout<< "refinementHistoryBalanced::resize from " << oldSize << " to " << size
            << " cells" << endl;
    }

    visibleCells_.setSize(size);

    // Set additional elements to -1.
    for (label i = oldSize; i < visibleCells_.size(); i++)
    {
        visibleCells_[i] = -1;
    }
}


void Foam::refinementHistoryBalanced::updateMesh(const mapPolyMesh& map)
{
    if (active())
    {
        const labelList& reverseCellMap = map.reverseCellMap();

        // Note that only the live cells need to be renumbered.

        labelList newVisibleCells(map.cellMap().size(), -1);

        forAll(visibleCells_, cellI)
        {
            if (visibleCells_[cellI] != -1)
            {
                label index = visibleCells_[cellI];

                // Check not already set
                if (splitCells_[index].addedCellsPtr_.valid())
                {
                    FatalErrorIn
                    (
                        "refinementHistoryBalanced::updateMesh(const mapPolyMesh&)"
                    )   << "Problem" << abort(FatalError);
                }

                label newCellI = reverseCellMap[cellI];

                if (newCellI >= 0)
                {
                    newVisibleCells[newCellI] = index;
                }
            }
        }

        if (debug)
        {
            Pout<< "refinementHistoryBalanced::updateMesh : from "
                << visibleCells_.size()
                << " to " << newVisibleCells.size()
                << " cells" << endl;
        }

        visibleCells_.transfer(newVisibleCells);
    }
}


// Update numbering for subsetting
void Foam::refinementHistoryBalanced::subset
(
    const labelList& pointMap,
    const labelList& faceMap,
    const labelList& cellMap
)
{
    if (active())
    {
        labelList newVisibleCells(cellMap.size(), -1);

        forAll(newVisibleCells, cellI)
        {
            label oldCellI = cellMap[cellI];

            label index = visibleCells_[oldCellI];

            // Check that cell is live (so its parent has no refinement)
            if (index >= 0 && splitCells_[index].addedCellsPtr_.valid())
            {
                FatalErrorIn
                (
                    "refinementHistoryBalanced::subset"
                    "(const labelList&, const labelList&, const labelList&)"
                )   << "Problem" << abort(FatalError);
            }

            newVisibleCells[cellI] = index;
        }

        if (debug)
        {
            Pout<< "refinementHistoryBalanced::updateMesh : from "
                << visibleCells_.size()
                << " to " << newVisibleCells.size()
                << " cells" << endl;
        }

        visibleCells_.transfer(newVisibleCells);
    }
}


void Foam::refinementHistoryBalanced::countProc
(
    const label index,
    const label newProcNo,
    labelList& splitCellProc,
    labelList& splitCellNum
) const
{
    if (splitCellProc[index] != newProcNo)
    {
        // Different destination processor from other cells using this
        // parent. Reset count.
        splitCellProc[index] = newProcNo;
        splitCellNum[index] = 1;
    }
    else
    {
        splitCellNum[index]++;

        // Increment parent if whole splitCell moves to same processor
        if (splitCellNum[index] == 8)
        {
            if (debug)
            {
                Pout<< "Moving " << splitCellNum[index]
                    << " cells originating from cell " << index
                    << " from processor " << Pstream::myProcNo()
                    << " to processor " << splitCellProc[index]
                    << endl;
            }

            label parent = splitCells_[index].parent_.index();

            if (parent >= 0)
            {
                countProc(parent, newProcNo, splitCellProc, splitCellNum);
            }
        }
    }
}

void Foam::refinementHistoryBalanced::parallelCountProc
(
	const label index,
	const label newProcNo,
    labelListList& splitCellNum
)
{
	if(index >= 0)
	{
		label parentIndex = splitCells_[index].parent_.index();
		if (parentIndex >= 0)
		{
			label parentProc = splitCells_[index].parent_.proc();

			if(parentProc != newProcNo)
				splitCellNum[parentProc][parentIndex]++;
		}
	}
}

void Foam::refinementHistoryBalanced::parallelAddProc
(
	const label index,
	const labelListList& splitCellNum,
	List<labelHashSet>& parents
)
{
	label parentIndex = splitCells_[index].parent_.index();
	label parentProc = splitCells_[index].parent_.proc();
	if (parentIndex >= 0)
	{
		if(splitCellNum[parentProc][parentIndex] == 8)
		{
			parents[parentProc].insert(parentIndex);
		}
	}
}

void Foam::refinementHistoryBalanced::distribute(const mapDistributePolyMesh& map)
{
    if (!active())
    {
        FatalErrorIn
        (
            "refinementHistoryBalanced::distribute(const mapDistributePolyMesh&)"
        )   << "Calling distribute on inactive history" << abort(FatalError);
    }


    if (!Pstream::parRun())
    {
        return;
    }

    // Remove unreferenced history.
    compact();

    //Pout<< nl << "--BEFORE:" << endl;
    //writeDebug();
    //Pout<< "---------" << nl << endl;

    // Determine clusters. This is per every entry in splitCells_ (that is
    // a parent of some refinement) a label giving the processor it goes to
    // if all its children are going to the same processor.

    // Per visible cell the processor it goes to.
    labelList destination(visibleCells_.size());

    const labelListList& subCellMap = map.cellMap().subMap();

    forAll(subCellMap, procI)
    {
        const labelList& newToOld = subCellMap[procI];

        forAll(newToOld, i)
        {
            label oldCellI = newToOld[i];

            destination[oldCellI] = procI;
        }
    }

//Pout<< "refinementHistoryBalanced::distribute :"
//    << " destination:" << destination << endl;

    // Per splitCell entry the processor it moves to
    labelListList splitCellProc(Pstream::nProcs());
	splitCellProc[Pstream::myProcNo()] = labelList(splitCells_.size(), Pstream::myProcNo());
//    Pstream::gatherList(splitCellProc);
//    Pstream::scatterList(splitCellProc);

    // Per splitCell entry the number of live cells that are off the processor
    // (we need to relocate the splitCell if none off its live cells are on this processor)
    labelListList splitCellNum(Pstream::nProcs());
    splitCellNum[Pstream::myProcNo()] = labelList(splitCells_.size(), 0);

    //Debug
//    labelList splitCellsSize(Pstream::nProcs());
//    splitCellsSize[Pstream::myProcNo()] = splitCells_.size();
//    Pstream::gatherList(splitCellsSize);
//    Info << "splitCells.size() = " << splitCellsSize << endl;
//    labelListList allVisibleCells(Pstream::nProcs());
//    allVisibleCells[Pstream::myProcNo()] = visibleCells_;
//    Pstream::gatherList(allVisibleCells);
//    Info << "visibleCells" << allVisibleCells << endl;

    Pstream::gatherList(splitCellNum);
    Pstream::scatterList(splitCellNum);

//    Info << "splitCellNum = " << splitCellNum << endl;

    // Per processor, local splitCell indices that are moving off processor
    List<labelHashSet> parents(Pstream::nProcs());
	forAll(visibleCells_, cellI)
	{
		label index = visibleCells_[cellI];
		parallelCountProc(index, destination[cellI], splitCellNum);
	}

	forAll(splitCellNum,procI)
	{
		Pstream::listCombineGather(splitCellNum[procI],plusEqOp<label>());
		Pstream::listCombineScatter(splitCellNum[procI]);
	}

//	Info << "splitCellNum = " << splitCellNum << endl;

	forAll(visibleCells_, cellI)
	{
		label index = visibleCells_[cellI];
		if(index >= 0)
		{
			splitCellProc[Pstream::myProcNo()][index] = destination[cellI];
			parallelAddProc(index, splitCellNum, parents);
		}
	}

	Pstream::gatherList(splitCellProc);
	Pstream::scatterList(splitCellProc);

	for (label procI = 0; procI < Pstream::nProcs(); procI++)
	{
        // Send to neighbors
        OPstream toNbr(Pstream::blocking, procI);
        toNbr << parents[procI];
        parents[procI].clear();
	}

	for (label procI = 0; procI < Pstream::nProcs(); procI++)
	{
		// Receive from neighbors
        IPstream fromNbr(Pstream::blocking, procI);
        labelHashSet newParents(fromNbr);
        parents[Pstream::myProcNo()] |= newParents;
	}

	forAllConstIter(labelHashSet,parents[Pstream::myProcNo()],iter)
	{
		procIndexPair& firstChild = splitCells_[iter.key()].addedCellsPtr_()[0];
		splitCellProc[Pstream::myProcNo()][iter.key()] = splitCellProc[firstChild.proc()][firstChild.index()];
	}

	Pstream::gatherList(splitCellProc);
	Pstream::scatterList(splitCellProc);

	bool complete = parents[Pstream::myProcNo()].size() == 0;
	reduce(complete, andOp<bool>());

	if(debug)
	{
		labelList parentsSize(Pstream::nProcs());
		parentsSize[Pstream::myProcNo()] = parents.size();
		Pstream::gatherList(parentsSize);
		Info << "parents.size() = " << parentsSize << endl;
	}

	while(!complete)
	{
		forAll(splitCellNum,procI)
		{
			if(procI != Pstream::myProcNo())
			{
				forAll(splitCellNum[procI],i)
				{
					splitCellNum[procI][i] = 0;
				}
			}
		}

		complete = true;

		forAllConstIter(labelHashSet,parents[Pstream::myProcNo()],iter)
		{
			complete = false;
			label parent = iter.key();
			parallelCountProc(parent, splitCellProc[Pstream::myProcNo()][parent], splitCellNum);
		}

		forAll(splitCellNum,procI)
		{
			Pstream::listCombineGather(splitCellNum[procI],plusEqOp<label>());
			Pstream::listCombineScatter(splitCellNum[procI]);
		}

		{
			List<labelHashSet> nextParents(Pstream::nProcs());
			forAllConstIter(labelHashSet,parents[Pstream::myProcNo()],iter)
			{
				complete = false;
				label parent = iter.key();
				parallelAddProc(parent, splitCellNum, nextParents);
			}
			parents = nextParents;
		}

		for (label procI = 0; procI < Pstream::nProcs(); procI++)
		{
	        // Send to neighbors
	        OPstream toNbr(Pstream::blocking, procI);
	        toNbr << parents[procI];
	        parents[procI].clear();
		}

		for (label procI = 0; procI < Pstream::nProcs(); procI++)
		{
			// Receive from neighbors
	        IPstream fromNbr(Pstream::blocking, procI);
	        labelHashSet newParents(fromNbr);
	        parents[Pstream::myProcNo()] |= newParents;
		}

		forAllConstIter(labelHashSet,parents[Pstream::myProcNo()],iter)
		{
			procIndexPair& firstChild = splitCells_[iter.key()].addedCellsPtr_()[0];
			splitCellProc[Pstream::myProcNo()][iter.key()] = splitCellProc[firstChild.proc()][firstChild.index()];
		}
		Pstream::gatherList(splitCellProc);
		Pstream::scatterList(splitCellProc);

		reduce(complete, andOp<bool>());
    }

//	//Debugging
//	std::ostringstream timeInfoFilename;
//	timeInfoFilename << time().timeName() << "_" << Pstream::myProcNo();
//	std::ofstream timeInfo(timeInfoFilename.str().c_str());
//
//	forAll(splitCellProc[Pstream::myProcNo()],i)
//	{
//		bool fail = true;
//		if(splitCells_[i].addedCellsPtr_.valid())
//		{
//			const FixedList<procIndexPair, 8>& splits = splitCells_[i].addedCellsPtr_();
//			forAll(splits, j)
//			{
//				if (splits[j].index() >= 0)
//				{
//					if(splitCellProc[splits[j].proc()][splits[j].index()] == splitCellProc[Pstream::myProcNo()][i])
//						fail = false;
//				}
//			}
//		}
//		else
//			fail = false;
//
//		if(fail)
//		{
//			timeInfo << "(" << Pstream::myProcNo() << "," << i << ") -> " << splitCellProc[Pstream::myProcNo()][i] << std::endl;
//			if(splitCells_[i].addedCellsPtr_.valid())
//			{
//				const FixedList<procIndexPair, 8>& splits = splitCells_[i].addedCellsPtr_();
//				forAll(splits, j)
//				{
//					if (splits[j].index() >= 0)
//					{
//						timeInfo << "    (" << splits[j].proc() << "," << splits[j].index() << ") -> " << splitCellProc[splits[j].proc()][splits[j].index()] << std::endl;
//					}
//				}
//			}
//			else
//			{
//				timeInfo << "    No children" << std::endl;
//			}
//		}
//	}
//	timeInfo.close();

    //Pout<< "refinementHistoryBalanced::distribute :"
    //    << " splitCellProc:" << splitCellProc << endl;
    //
    //Pout<< "refinementHistoryBalanced::distribute :"
    //    << " splitCellNum:" << splitCellNum << endl;


    // Create subsetted refinement tree consisting of all parents that
    // move in their whole to other processor.
	List< DynamicList<splitCell8> > newSplitCells(Pstream::nProcs());
	labelListList oldToNew(Pstream::nProcs());
	oldToNew[Pstream::myProcNo()] = labelList(splitCells_.size(), -1);

    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        //Pout<< "-- Subetting for processor " << procI << endl;

        // Loop over all entries.
        forAll(splitCells_, index)
        {
//            Pout<< "oldCell:" << index
//                << " proc:" << splitCellProc[index]
//                << " nCells:" << splitCellNum[index]
//                << endl;

            if (splitCellProc[Pstream::myProcNo()][index] == procI)
            {
                // Entry moves to procI
                oldToNew[Pstream::myProcNo()][index] = newSplitCells[procI].size();
                newSplitCells[procI].append(splitCells_[index]);

                //Pout<< "Added oldCell " << index
                //    << " info " << newSplitCells.last()
                //    << " at position " << newSplitCells.size()-1
                //    << endl;
            }
        }

        newSplitCells[procI].shrink();
    }

    forAll(newSplitCells,procI)
    {

    	labelList size(Pstream::nProcs());
    	size[Pstream::myProcNo()] = newSplitCells[procI].size();
    	Pstream::gatherList(size);
    	Pstream::scatterList(size);

		label indexOffset = 0;
		for(label procIJ = 0; procIJ < Pstream::myProcNo(); procIJ++)
		{
			indexOffset += size[procIJ];
		}

		forAll(oldToNew[Pstream::myProcNo()], index)
		{
			if(splitCellProc[Pstream::myProcNo()][index] == procI)
				oldToNew[Pstream::myProcNo()][index] += indexOffset;
		}
    }

	Pstream::gatherList(oldToNew);
	Pstream::scatterList(oldToNew);

	forAll(newSplitCells,procI)
	{
    	forAll(newSplitCells[procI],index)
    	{
    		splitCell8& split = newSplitCells[procI][index];
    		label oldParentIndex = split.parent_.index();
    		label oldParentProc = split.parent_.proc();
    		if(oldParentIndex >= 0)
    		{
    			split.parent_.index() = oldToNew[oldParentProc][oldParentIndex];
    			split.parent_.proc() = splitCellProc[oldParentProc][oldParentIndex];
    		}

            if (split.addedCellsPtr_.valid())
            {
            	//Update sub cell entries
                FixedList<procIndexPair, 8>& splits = split.addedCellsPtr_();

                forAll(splits, i)
                {
                    if (splits[i].index() >= 0)
                    {
                    	label oldSplitProc = splits[i].proc();
                    	label oldSplitIndex = splits[i].index();
                    	splits[i].proc() = splitCellProc[oldSplitProc][oldSplitIndex];
                        splits[i].index() = oldToNew[oldSplitProc][oldSplitIndex];
                    }
                }
            }
    	}
    }

    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {

        const labelList& subMap = subCellMap[procI];

        // New visible cells.
        labelList newVisibleCells(subMap.size(), -1);

        forAll(subMap, newCellI)
        {
            label oldCellI = subMap[newCellI];

            label oldIndex = visibleCells_[oldCellI];

            if (oldIndex >= 0) //if it has a splitCell
            {
                newVisibleCells[newCellI] = oldToNew[Pstream::myProcNo()][oldIndex];
            }
        }

        //Pout<< nl << "--Subset for domain:" << procI << endl;
        //writeDebug(newVisibleCells, newSplitCells);
        //Pout<< "---------" << nl << endl;


        // Send to neighbours
        OPstream toNbr(Pstream::blocking, procI);
        toNbr << newSplitCells[procI] << newVisibleCells;
    }


    // Receive from neighbours and merge
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Remove all entries. Leave storage intact.
    splitCells_.clear();

    visibleCells_.setSize(map.mesh().nCells());
    visibleCells_ = -1;

    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        IPstream fromNbr(Pstream::blocking, procI);
        List<splitCell8> newSplitCells(fromNbr);
        labelList newVisibleCells(fromNbr);

        //Pout<< nl << "--Received from domain:" << procI << endl;
        //writeDebug(newVisibleCells, newSplitCells);
        //Pout<< "---------" << nl << endl;


        // newSplitCells contain indices only into newSplitCells so
        // renumbering can be done here.
//        label offset = splitCells_.size();

        //Pout<< "**Renumbering data from proc " << procI << " with offset "
        //    << offset << endl;

//        //Need to do this before we communicate because parents may be off processor
//        forAll(newSplitCells, index)
//        {
//            splitCell8& split = newSplitCells[index];
//
//            if (split.parent_.index()>= 0)
//            {
//                split.parent_.index()+= offset;
//            }
//            if (split.addedCellsPtr_.valid())
//            {
//                FixedList<procIndexPair, 8>& splits = split.addedCellsPtr_();
//
//                forAll(splits, i)
//                {
//                    if (splits[i].index() >= 0)
//                    {
//                        splits[i].index() += offset;
//                    }
//                }
//            }
//
//            splitCells_.append(split);
//        }

        splitCells_.append(newSplitCells);

//        Info << "Debug point at " << __FILE__ << ':' << __LINE__ << endl;

        // Combine visibleCell.
        const labelList& constructMap = map.cellMap().constructMap()[procI];

        forAll(newVisibleCells, i)
        {
            if (newVisibleCells[i] >= 0)
            {
                visibleCells_[constructMap[i]] = newVisibleCells[i];
            }
        }
    }
    splitCells_.shrink();

//    Info << "Debug point at " << __FILE__ << ':' << __LINE__ << endl;

    //Pout<< nl << "--AFTER:" << endl;
    //writeDebug();
    //Pout<< "---------" << nl << endl;
}

void Foam::refinementHistoryBalanced::compact()
{
    if (debug)
    {
        Pout<< "refinementHistoryBalanced::compact() Entering with:"
            << " freeSplitCells_:" << freeSplitCells_.size()
            << " splitCells_:" << splitCells_.size()
            << " visibleCells_:" << visibleCells_.size()
            << endl;

        // Check all free splitCells are marked as such
        forAll(freeSplitCells_, i)
        {
            label index = freeSplitCells_[i];

            if (splitCells_[index].parent_.index() != -2)
            {
                FatalErrorIn("refinementHistoryBalanced::compact()")
                    << "Problem index:" << index
                    << abort(FatalError);
            }
        }

        // Check none of the visible cells are marked as free
        forAll(visibleCells_, cellI)
        {
            if
            (
                visibleCells_[cellI] >= 0
             && splitCells_[visibleCells_[cellI]].parent_.index() == -2
            )
            {
                FatalErrorIn("refinementHistoryBalanced::compact()")
                    << "Problem : visible cell:" << cellI
                    << " is marked as being free." << abort(FatalError);
            }
        }
    }

//    Info << "Debug point at " << __FILE__ << ':' << __LINE__ << endl;

    DynamicList<splitCell8> newSplitCells(splitCells_.size());

    // From uncompacted to compacted splitCells.
    labelList oldToNew(splitCells_.size(), -1);

    // Mark all used splitCell entries. These are either indexed by visibleCells
    // or indexed from other splitCell entries.

    // Mark from visibleCells

    boolListList marked(Pstream::nProcs());
    marked[Pstream::myProcNo()] = boolList(splitCells_.size(),false);

    forAll(visibleCells_, cellI)
    {
        label index = visibleCells_[cellI];

        if (index >= 0)
        {
            // Make sure we only mark visible indices if they either have a
            // parent or subsplits.
            if
            (
                splitCells_[index].parent_.index() != -1
             || splitCells_[index].addedCellsPtr_.valid()
            )
            {
                //markSplit(index, oldToNew, newSplitCells);
            	marked[Pstream::myProcNo()][index] = true;
            }
        }
    }

    Pstream::gatherList(marked);
    Pstream::scatterList(marked);

//    Info << "Debug point at " << __FILE__ << ':' << __LINE__ << endl;

    parallelMarkSplits(marked, oldToNew, newSplitCells);

    // Mark from splitCells
    forAll(splitCells_, index)
    {
        if (splitCells_[index].parent_.index() == -2)
        {
            // freed cell.
        }
        else if
        (
            splitCells_[index].parent_.index() == -1
         && splitCells_[index].addedCellsPtr_.empty()
        )
        {
            // recombined cell. No need to keep since no parent and no subsplits
            // Note that gets marked if reachable from other index!
        }
        else
        {
            // Is used element.
            //markSplit(index, oldToNew, newSplitCells);
        	marked[Pstream::myProcNo()][index] = true;
        }
    }

    Pstream::gatherList(marked);
    Pstream::scatterList(marked);

    parallelMarkSplits(marked, oldToNew, newSplitCells);

    // Now oldToNew is fully complete and compacted elements are in
    // newSplitCells.
    // Renumber contents of newSplitCells and visibleCells.
    labelListList allOldToNew(Pstream::nProcs());
    allOldToNew[Pstream::myProcNo()] = oldToNew;
    Pstream::gatherList(allOldToNew);
    Pstream::scatterList(allOldToNew);
    forAll(newSplitCells, index)
    {
        splitCell8& split = newSplitCells[index];

        if (split.parent_.index() >= 0)
        {
            split.parent_.index() = allOldToNew[split.parent_.proc()][split.parent_.index()];
        }
        if (split.addedCellsPtr_.valid())
        {
            FixedList<procIndexPair, 8>& splits = split.addedCellsPtr_();

            forAll(splits, i)
            {
                if (splits[i].index() >= 0)
                {
                    splits[i].index() = allOldToNew[splits[i].proc()][splits[i].index()];
                }
            }
        }
    }


    if (debug)
    {
        Pout<< "refinementHistoryBalanced::compact : compacted splitCells from "
            << splitCells_.size() << " to " << newSplitCells.size() << endl;
    }

    splitCells_.transfer(newSplitCells);
    freeSplitCells_.clearStorage();


    if (debug)
    {
        Pout<< "refinementHistoryBalanced::compact() NOW:"
            << " freeSplitCells_:" << freeSplitCells_.size()
            << " splitCells_:" << splitCells_.size()
            << " newSplitCells:" << newSplitCells.size()
            << " visibleCells_:" << visibleCells_.size()
            << endl;
    }


    // Adapt indices in visibleCells_
    forAll(visibleCells_, cellI)
    {
        label index = visibleCells_[cellI];

        if (index >= 0)
        {
            // Note that oldToNew can be -1 so it resets newVisibleCells.
            visibleCells_[cellI] = oldToNew[index];
        }
        else
        {
            // Keep -1 value.
        }
    }
}


void Foam::refinementHistoryBalanced::writeDebug() const
{
    writeDebug(visibleCells_, splitCells_);
}


void Foam::refinementHistoryBalanced::storeSplit
(
    const label cellI,
    const labelList& addedCells
)
{
    label parentIndex = -1;

    if (visibleCells_[cellI] != -1)
    {
        // Was already live. The current live cell becomes the
        // parent of the cells split off from it.

        parentIndex = visibleCells_[cellI];

        // It is no longer live (note that actually cellI gets alive
        // again below since is addedCells[0])
        visibleCells_[cellI] = -1;
    }
    else
    {
        // Create 0th level. -1 parent to denote this.
        parentIndex = allocateSplitCell(-1, -1);
    }

    // Create live entries for added cells that point to the
    // cell they were created from (parentIndex)
    forAll(addedCells, i)
    {
        label addedCellI = addedCells[i];

        // Create entries for the split off cells. All of them
        // are visible.
        visibleCells_[addedCellI] = allocateSplitCell(parentIndex, i);
    }
}


void Foam::refinementHistoryBalanced::combineCells
(
    const label masterCellI,
    const labelList& combinedCells
)
{
    // Save the parent structure
    label parentIndex = splitCells_[visibleCells_[masterCellI]].parent_.index();
	label parentProc = splitCells_[visibleCells_[masterCellI]].parent_.proc();

	if(parentProc != Pstream::myProcNo())
	{
		write();
        FatalErrorIn("refinementHistoryBalanced::combineCells")
            << "Problem: parent is not on my process" << endl
			<< "parentIndex: " << parentIndex << endl
			<< "parentProc: " << parentProc << endl
			<< "masterCellI: " << masterCellI << endl
            << abort(FatalError);
	}

    // Remove the information for the combined cells
    forAll(combinedCells, i)
    {
        label cellI = combinedCells[i];

        label index = visibleCells_[cellI];
        splitCell8& split = splitCells_[index];

        // Make sure parent does not point to me anymore.
        if (split.parent_.index() >= 0)
        {

            autoPtr<FixedList<procIndexPair, 8> >& subCellsPtr =
                splitCells_[split.parent_.index()].addedCellsPtr_;

            if (subCellsPtr.valid())
            {
                FixedList<procIndexPair, 8>& subCells = subCellsPtr();

                label myPos = findIndex(subCells, procIndexPair(Pstream::myProcNo(),index));

                if (myPos == -1)
                {
                    FatalErrorIn("refinementHistoryBalanced::freeSplitCell")
                        << "Problem: cannot find myself in"
                        << " parents' children" << abort(FatalError);
                }
                else
                {
                    subCells[myPos] = procIndexPair(-1,-1);
                }
            }
        }

        // Mark splitCell as free
        split.parent_.index()= -2;

        // Add to cache of free splitCells
        freeSplitCells_.append(index);
        visibleCells_[cellI] = -1;
    }

    if (parentProc == Pstream::myProcNo())
    {
		splitCell8& parentSplit = splitCells_[parentIndex];
		parentSplit.addedCellsPtr_.reset(NULL);
		visibleCells_[masterCellI] = parentIndex;
    }
}


bool Foam::refinementHistoryBalanced::readData(Istream& is)
{
    is >> *this;
    return !is.bad();
}


bool Foam::refinementHistoryBalanced::writeData(Ostream& os) const
{
    os << *this;

    return os.good();
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, refinementHistoryBalanced& rh)
{
    rh.freeSplitCells_.clearStorage();

    is >> rh.splitCells_ >> rh.visibleCells_;

    // Check indices.
    rh.checkIndices();

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const refinementHistoryBalanced& rh)
{
    const_cast<refinementHistoryBalanced&>(rh).compact();

    return os   << "// splitCells" << nl
                << rh.splitCells_ << nl
                << "// visibleCells" << nl
                << rh.visibleCells_;
}


// ************************************************************************* //
