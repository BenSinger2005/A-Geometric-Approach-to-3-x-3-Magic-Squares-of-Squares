This repository contains the Magma code used to verify the computational results in my paper A Geometric Approach to Magic Squares of Squares (with Asher Auel, arxiv ID to be uploaded). The paper studies the classical problem of 3×3 magic squares of squares using methods from algebraic geometry, including the geometry of quartic del Pezzo surfaces and K3 surfaces. This repository contains all computations required to reproduce the main computational results. The paper is heavily reliant on Magma, particularly for the Picard rank, classifying magic del Pezzo surfaces, and computing Weil polynomials of magic K3 surfaces.
Below is a guide to which files are used to prove which results in the paper.

**THEOREM 3: Picard Rank Bound:**
This is the most computationally intensive result in the paper. It involves three files: Comps.m, ComputeLinePreimages.m, PicardRankCodeSquareOfSquares.m, and IntersectionMatrix.m. PicardRankCodeSquareOfSquares contains code to compute the Picard ranks of both the square of squares variety V and the K3 surface X from section 4.2. 
To run the verification for Theorem 2, paste in PicardRankCodeSquareOfSquares.m up until retrieving the singular points of V. Then, paste in Comps.m. Afterwards, paste the rest of the code up through the first instance of Rank(pairingmat). The computation should take ~2 hours. The raw intersection matrix is also available in the file IntersectionMatrix.m. As far as ComputeLinePreimages.m goes, this is to compute the 532 divisors appearing as preimages of lines on the magic quartic del Pezzo and cubic surfaces derived in Section 4.1. The equations are already in Comps.m, but the curious reader may check that the code gives the divisors we claim to have.

**Hardware note.** Computing the Picard rank of V requires substantially more memory than is typically available on a desktop/laptop installation of Magma. We recommend running this computation on a high-memory computing cluster.

**PROPOSITION 4.2: Classification of Magic Quartic del Pezzo Surfaces:**

The classification of magic quartic del Pezzo surfaces is entirely contained in delPezzoClassification. It is done by computing the Segre symbol of each model, as well as the PGL_2 orbit of the multiset of coefficients in P^1; this completely determines the geometric isomorphism class. The instructions for how to run it are given in the file.

**PROPOSITION 4.8: Geometric Picard Rank of the Magic Octic K3 Surface X:**

This proposition requires 3 files: MagicK3Comps.m, PicardRankCodeSquareOfSquares.m, and MagicK3WeilPolynomials.txt. Verifying the lower bound on the Picard rank of X is similar to the method for Theorem 2. After the end of the code to compute the bound on the Picard rank for V, the code for X begins. Paste it up until the singular points, and then paste MagicK3Comps.m in. Then, paste the rest of the code. The computation takes a few minutes. The Weil polynomials of X are all given in MagicK3WeilPolynomials.txt, along with extra data about the branch sextic and Frobenius traces.

**PROPOSITION 4.10: Classification of Magic Octic K3 Surfaces and Associated Degree 2 Surfaces:**

The classification of magic octic K3 surfaces and their associated double covers is entirely contained in AssociatedDoubleCoverAndLines. We give code to compute the associated double cover of an octic K3 surface, its splitting into lines, and the number of nodes/triple points. More directions are given in the file.

**PROPOSITION 4.11: Verifying Degree of Dominant Rational Map from Bremner's Magic K3:**

In the proof of Proposition 4.11, we show that Bremner's smooth magic octic K3 surface Y is isogenous to a 15 dodecic (degree 12) magic K3 surfaces, enumerated by sorting 6 coordinates into 3 unordered pairs of 2. This involves showing that a projection p from Y is a dominant rational map generically of degree 4, which is reliant on showing that a particular matrix N has nonvanishing 2 x 2 minors for each choice of pairs. This is done in PairingCheck.m. This is a very quick verification; we include it for completeness' sake.

**APPENDIX A: Verifying that V has No Lines:**

We show that V has no lines in Appendix A, building on a Schubert cell stratification approach of Elsenhans--Jahnel. The accompanying file is NoLines.m. This line-counting approach is also used in ComputeLinePreimages.m.

Enjoy! Questions, comments, or bug reports are welcome through GitHub Issues.
