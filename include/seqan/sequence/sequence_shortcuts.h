// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Shortcuts for certain often-used string types.
// ==========================================================================

#ifndef SEQAN_HEADER_SEQUENCE_SHORTCUTS_H
#define SEQAN_HEADER_SEQUENCE_SHORTCUTS_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

/*!
 * @typedef CharString
 * @headerfile <seqan/sequence.h>
 * @brief An AllocString with char alphabet.
 *
 * You can efficiently cast a CharString to <tt>char *</tt> using <tt>toCString</tt>.
 *
 * @signature typedef String<char> CharString;
 */

/**
.Shortcut.CharString:
..cat:Strings
..summary:A string of $char$.
..signature:CharString
..description:
This is a useful replacement of $std::string$.
If you want a C-style $char *$ string, use @Function.toCString@.
..shortcutfor:Spec.Alloc String
...signature:String<char, Alloc<> >
*/

typedef String<char, Alloc<void> > CharString;

//____________________________________________________________________________

/*!
 * @typedef CharIterator
 * @headerfile <seqan/sequence.h>
 * @brief An iterator intoa @link CharString @endlink.
 *
 * @signature typedef Iterator<CharString, Rooted>::Type CharIterator;
 */

/**
.Shortcut.CharIterator:
..cat:Iterators
..summary:Iterator for @Shortcut.CharString@.
..signature:CharIterator
..shortcutfor:Concept.RootedIteratorConcept
...signature:Iterator<CharString, Rooted>::Type
..see:Shortcut.CharString
*/

typedef Iterator<CharString, Rooted>::Type CharIterator;

//////////////////////////////////////////////////////////////////////////////

/*!
 * @typedef UnicodeString
 * @headerfile <seqan/sequence.h>
 * @brief An AllocString with wchar_t alphabet.
 *
 * @signature typedef String<wchar_t> CharString;
 */

/**
.Shortcut.UnicodeString:
..cat:Strings
..summary:A string of $wchar_t$.
..signature:UnicodeString
..shortcutfor:Spec.Alloc String
...signature:String<wchar_t, Alloc<> >
*/

typedef String<wchar_t, Alloc<void> > UnicodeString;

//____________________________________________________________________________

/**
.Shortcut.UnicodeIterator:
..cat:Iterators
..summary:Iterator for @Shortcut.UnicodeString@.
..signature:UnicodeIterator
..shortcutfor:Concept.RootedIteratorConcept
...signature:Iterator<UnicodeString, Rooted>::Type
..see:Shortcut.UnicodeString
*/

typedef Iterator<UnicodeString, Rooted>::Type UnicodeIterator;

//////////////////////////////////////////////////////////////////////////////

/*!
 * @typedef DnaString
 * @headerfile <seqan/sequence.h>
 * @brief An AllocString with Dna alphabet.
 *
 * @signature typedef String<Dna> CharString;
 */

/**
.Shortcut.DnaString:
..cat:Strings
..summary:A string of @Spec.Dna@.
..signature:DnaString
..shortcutfor:Spec.Alloc String
...signature:String<Dna, Alloc<> >
..see:Spec.Dna
*/

typedef String<Dna, Alloc<void> > DnaString;

//____________________________________________________________________________

/*!
 * @typedef DnaIterator
 * @headerfile <seqan/sequence.h>
 * @brief A rooted iterator over a Dna.
 *
 * @signature typedef Iterator<DnaString, Rooted>::Type DnaIterator;
 */

/**
.Shortcut.DnaIterator:
..cat:Iterators
..summary:Iterator for @Shortcut.DnaString@.
..signature:DnaIterator
..shortcutfor:Concept.RootedIteratorConcept
...signature:Iterator<DnaString, Rooted>::Type
..see:Spec.Dna
..see:Shortcut.DnaString
*/

typedef Iterator<DnaString, Rooted>::Type DnaIterator;

//////////////////////////////////////////////////////////////////////////////

/*!
 * @typedef Dna5String
 * @headerfile <seqan/sequence.h>
 * @brief An AllocString with Dna5 alphabet.
 *
 * @signature typedef String<Dna5> CharString;
 */

/**
.Shortcut.Dna5String:
..cat:Strings
..summary:A string of @Spec.Dna5@.
..signature:Dna5String
..shortcutfor:Spec.Alloc String
...signature:String<Dna5, Alloc<> >
..see:Spec.Dna5
..see:Shortcut.DnaString
*/

typedef String<Dna5, Alloc<void> > Dna5String;

//____________________________________________________________________________

/*!
 * @typedef Dna5Iterator
 * @headerfile <seqan/sequence.h>
 * @brief A rooted iterator over a Dna5.
 *
 * @signature typedef Iterator<Dna5String, Rooted>::Type Dna5Iterator;
 */

/**
.Shortcut.Dna5Iterator:
..cat:Iterators
..summary:Iterator for @Shortcut.Dna5String@.
..signature:Dna5Iterator
..shortcutfor:Concept.RootedIteratorConcept
...signature:Iterator<Dna5String, Rooted>::Type
..see:Spec.Dna5
..see:Shortcut.Dna5String
..see:Shortcut.DnaIterator
*/

typedef Iterator<Dna5String, Rooted>::Type Dna5Iterator;

//////////////////////////////////////////////////////////////////////////////

/*!
 * @typedef RnaString
 * @headerfile <seqan/sequence.h>
 * @brief An AllocString with Rna alphabet.
 *
 * @signature typedef String<Rna> CharString;
 */

/**
.Shortcut.RnaString:
..cat:Strings
..summary:A string of @Spec.Rna@.
..signature:RnaString
..shortcutfor:Spec.Alloc String
...signature:String<Rna, Alloc<> >
..see:Spec.Rna
*/

typedef String<Rna, Alloc<void> > RnaString;

//____________________________________________________________________________

/*!
 * @typedef RnaIterator
 * @headerfile <seqan/sequence.h>
 * @brief A rooted iterator over a Rna.
 *
 * @signature typedef Iterator<RnaString, Rooted>::Type RnaIterator;
 */

/**
.Shortcut.RnaIterator:
..cat:Iterators
..summary:Iterator for @Shortcut.RnaString@.
..signature:RnaIterator
..shortcutfor:Concept.RootedIteratorConcept
...signature:Iterator<RnaString, Rooted>::Type
..see:Spec.Rna
..see:Shortcut.RnaString
*/

typedef Iterator<RnaString, Rooted>::Type RnaIterator;

//////////////////////////////////////////////////////////////////////////////

/*!
 * @typedef Rna5String
 * @headerfile <seqan/sequence.h>
 * @brief An AllocString with Rna5 alphabet.
 *
 * @signature typedef String<Rna5> CharString;
 */

/**
.Shortcut.Rna5String:
..cat:Strings
..summary:A string of @Spec.Rna5@.
..signature:Rna5String
..shortcutfor:Spec.Alloc String
...signature:String<Rna5, Alloc<> >
..see:Spec.Rna5
..see:Shortcut.RnaString
*/

typedef String<Rna5, Alloc<void> > Rna5String;

//____________________________________________________________________________

/*!
 * @typedef Rna5Iterator
 * @headerfile <seqan/sequence.h>
 * @brief A rooted iterator over a Rna5.
 *
 * @signature typedef Iterator<Rna5String, Rooted>::Type Rna5Iterator;
 */

/**
.Shortcut.Rna5Iterator:
..cat:Iterators
..summary:Iterator for @Shortcut.Rna5String@.
..signature:Rna5Iterator
..shortcutfor:Concept.RootedIteratorConcept
...signature:Iterator<Rna5String, Rooted>::Type
..see:Spec.Rna5
..see:Shortcut.Rna5String
..see:Shortcut.RnaIterator
*/

typedef Iterator<Rna5String, Rooted>::Type Rna5Iterator;

//////////////////////////////////////////////////////////////////////////////

/*!
 * @typedef IupacString
 * @headerfile <seqan/sequence.h>
 * @brief An AllocString with Iupac alphabet.
 *
 * @signature typedef String<Iupac> CharString;
 */

/**
.Shortcut.IupacString:
..cat:Strings
..summary:A string of @Spec.Iupac@.
..signature:IupacString
..shortcutfor:Spec.Alloc String
...signature:String<Iupac, Alloc<> >
..see:Spec.Iupac
*/

typedef String<Iupac, Alloc<void> > IupacString;

//____________________________________________________________________________

/*!
 * @typedef IupacIterator
 * @headerfile <seqan/sequence.h>
 * @brief A rooted iterator over a Iupac.
 *
 * @signature typedef Iterator<IupacString, Rooted>::Type IupacIterator;
 */

/**
.Shortcut.IupacIterator:
..cat:Iterators
..summary:Iterator for @Shortcut.IupacString@.
..signature:IupacIterator
..shortcutfor:Concept.RootedIteratorConcept
...signature:Iterator<IupacString, Rooted>::Type
..see:Spec.Iupac
..see:Shortcut.IupacString
*/

typedef Iterator<IupacString, Rooted>::Type IupacIterator;

//////////////////////////////////////////////////////////////////////////////

/*!
 * @typedef Peptide
 * @headerfile <seqan/sequence.h>
 * @brief An AllocString with AminoAcid alphabet.
 *
 * @signature typedef String<AminoAcid> Peptide;
 */

/**
.Shortcut.Peptide:
..cat:Strings
..summary:A string of @Spec.AminoAcid@.
..signature:Peptide
..shortcutfor:Spec.Alloc String
...signature:String<AminoAcid, Alloc<> >
..see:Spec.AminoAcid
*/

typedef String<AminoAcid, Alloc<void> > Peptide;

//____________________________________________________________________________

/*!
 * @typedef PeptideIterator
 * @headerfile <seqan/sequence.h>
 * @brief A rooted iterator over a Peptide.
 *
 * @signature typedef Iterator<Peptide, Rooted>::Type PeptideIterator;
 */

/**
.Shortcut.PeptideIterator:
..cat:Iterators
..summary:Iterator for @Shortcut.Peptide@.
..signature:PeptideIterator
..shortcutfor:Concept.RootedIteratorConcept
...signature:Iterator<Peptide, Rooted>::Type
..see:Spec.AminoAcid
..see:Shortcut.Peptide
*/

typedef Iterator<Peptide, Rooted>::Type PeptideIterator;

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
