**Background:** Improvements in speed and cost of genome sequencing are resulting in increasing numbers of novel non-synonymous single nucleotide polymorphisms (nsSNPs) in genes known to be associated with disease. The large number of nsSNPs makes laboratory-based classification infeasible and familial co-segregation with disease is not always possible. *In-silico* methods for classification or triage are thus utilised. A popular tool based on multiple-species sequence alignments (MSAs), Align-GVGD, has been shown to underestimate deleterious effects, particularly as sequence numbers increase.

**Results:** We utilised the DEFLATE compression algorithm to account for expected variation across a number of species. With the adjusted Grantham measure we derived a means of quantitatively clustering known neutral and deleterious nsSNPs from the same gene; this was then used to assign novel variants to the most appropriate cluster as a means of binary classification. Scaling of clusters allows for inter-gene comparison of variants through a single pathogenecity score.

**Conclusions:** The approach improves upon the classification accuracy of Align-GVGD while correcting for sensitivity to large MSAs.

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
