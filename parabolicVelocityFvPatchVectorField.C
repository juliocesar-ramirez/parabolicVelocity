/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "addToRunTimeSelectionTable.H"
#include "boundBox.H"
#include "dimensionedType.H"
#include "fvPatchFieldMapper.H"
#include "parabolicVelocityFvPatchVectorField.H"
#include "surfaceFields.H"
#include "vectorField.H"
#include "volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::parabolicVelocityFvPatchVectorField::t() const
{
    return db().time().timeOutputValue();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parabolicVelocityFvPatchVectorField::parabolicVelocityFvPatchVectorField(
    const fvPatch &p, const DimensionedField<vector, volMesh> &iF)
    : fixedValueFvPatchVectorField(p, iF), Vmax_(0.0), y_(Zero), n_(Zero) ,
      wordData_("wordDefault"), labelData_(-1), boolData_(false) {}

Foam::parabolicVelocityFvPatchVectorField::parabolicVelocityFvPatchVectorField(
    const fvPatch &p, const DimensionedField<vector, volMesh> &iF,
    const dictionary &dict)
    : fixedValueFvPatchVectorField(p, iF), Vmax_(dict.lookup<scalar>("Vmax")),
      y_(dict.lookup<vector>("y")), n_(dict.lookup<vector>("n")),
      wordData_(dict.lookupOrDefault<word>("wordName", "wordDefault")),
      labelData_(-1), boolData_(false) {

  fixedValueFvPatchVectorField::evaluate();

  /*
  // Initialise with the value entry if evaluation is not possible
  fvPatchVectorField::operator=
  (
      vectorField("value", dict, p.size())
  );
  */
}

Foam::parabolicVelocityFvPatchVectorField::parabolicVelocityFvPatchVectorField(
    const parabolicVelocityFvPatchVectorField &ptf, const fvPatch &p,
    const DimensionedField<vector, volMesh> &iF,
    const fvPatchFieldMapper &mapper)
    : fixedValueFvPatchVectorField(ptf, p, iF, mapper), Vmax_(ptf.Vmax_),
      y_(ptf.y_), n_(ptf.n_), wordData_(ptf.wordData_), labelData_(-1),
      boolData_(ptf.boolData_) {}

Foam::parabolicVelocityFvPatchVectorField::parabolicVelocityFvPatchVectorField(
    const parabolicVelocityFvPatchVectorField &ptf,
    const DimensionedField<vector, volMesh> &iF)
    : fixedValueFvPatchVectorField(ptf, iF), Vmax_(ptf.Vmax_), y_(ptf.y_),
      n_(ptf.n_), wordData_(ptf.wordData_), labelData_(-1),
      boolData_(ptf.boolData_) {}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::parabolicVelocityFvPatchVectorField::autoMap
          (
    const fvPatchFieldMapper& m
)
          {
  fixedValueFvPatchVectorField::autoMap(m);
}

          void Foam::parabolicVelocityFvPatchVectorField::rmap(
    const fvPatchVectorField &ptf, const labelList &addr) {
  fixedValueFvPatchVectorField::rmap(ptf, addr);

}

          void Foam::parabolicVelocityFvPatchVectorField::updateCoeffs() {
  if (updated()) {
    return;
  }

  boundBox bb(patch().patch().localPoints(), true);

  vector vertice = 0.5 * (bb.max() + bb.min());

  vector radio = 0.5 * (bb.max() - bb.min());

  const vectorField &cellxy = patch().Cf();

  scalarField part1 = ((cellxy - vertice) & y_) / (radio & y_);

  fixedValueFvPatchVectorField::operator==(n_ *Vmax_ *(1 - sqr(part1)));

  fixedValueFvPatchVectorField::updateCoeffs();
}

          void Foam::parabolicVelocityFvPatchVectorField::write
          (
    Ostream& os
) const
          {
  fvPatchVectorField::write(os);
  writeEntry(os, "Vmax", Vmax_);
  writeEntry(os, "y", y_);
  writeEntry(os, "n", n_);
  writeEntry(os, "wordData", wordData_);
  writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

          namespace Foam
          {
  makePatchTypeField(fvPatchVectorField, parabolicVelocityFvPatchVectorField);
}

// ************************************************************************* //
