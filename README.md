#CompressGV

[Abstract](#abstract) | [Browser & HTTP Interface](#browser) | [Requirements](#req) | [Installation](#install)

##<a id="abstract"></a>Abstract

**Background:** Improvements in speed and cost of genome sequencing are resulting in increasing numbers of novel non-synonymous single nucleotide polymorphisms (nsSNPs) in genes known to be associated with disease. The large number of nsSNPs makes laboratory-based classification infeasible and familial co-segregation with disease is not always possible. *In-silico* methods for classification or triage are thus utilised. A popular tool based on multiple-species sequence alignments (MSAs), Align-GVGD, has been shown to underestimate deleterious effects, particularly as sequence numbers increase.

**Results:** We utilised the DEFLATE compression algorithm to account for expected variation across a number of species. With the adjusted Grantham measure we derived a means of quantitatively clustering known neutral and deleterious nsSNPs from the same gene; this was then used to assign novel variants to the most appropriate cluster as a means of binary classification. Scaling of clusters allows for inter-gene comparison of variants through a single pathogenecity score.

**Conclusions:** The approach improves upon the classification accuracy of Align-GVGD while correcting for sensitivity to large MSAs.

##<a id="browser"></a>Browser & HTTP Interface

A server running CompressGV is made available for anonymous usage; input data is destroyed immediately after classification and output is not stored. A [browser-based interface](http://compressgv.arranschlosberg.com) is available for standard users. Low-level [usage directly over HTTP](http://compressgv.arranschlosberg.com/#help) is encouraged for personal use only.

##<a id="req"></a>Requirements

Versions in brackets refer to those used in development; earlier versions may be compatible but have not been tested.

* [GNU Scientific Library](http://www.gnu.org/software/gsl/) + CBLAS (1.16)
* [zlib](http://www.zlib.net/) (1.2.8)
* [PSO](https://github.com/aschlosberg/pso) implementation ([kkentzo/pso 5a992b312e](https://github.com/kkentzo/pso/commit/5a992b312e21c421b363ed95cf5b0f7dede9890a))

##<a id="install"></a>Installation

A sample Linux makefile is included.

---------------------------------------------------------------------------------------

Copyright 2013 Arran Schlosberg.

This file is part of https://github.com/aschlosberg/CompressGV (CompressGV)

    CompressGV is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CompressGV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CompressGV. If not, see <http://www.gnu.org/licenses/>.

---------------------------------------------------------------------------------------
