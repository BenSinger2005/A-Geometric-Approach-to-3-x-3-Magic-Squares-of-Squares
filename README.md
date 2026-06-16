This repository contains code supporting the results in my paper with Asher Auel "A Geometric Approach to Magic Squares of Squares" (Arxiv number to be added in here).
The paper is heavily reliant on computations in magma, particularly for the Picard rank, classifying magic del Pezzo surfaces, and computing Weil polynomials of magic K3 surfaces.
Below is a guide to which files are used to prove which results in the paper.

*THEOREM 2: Picard Rank Bound*
4.5
This is the most computationally intensive theorem in the paper. It involves three files: Comps.m, ComputeLinePreimages.m, and PicardRankCodeSquareOfSquares.m. PicardRankCodeSquareOfSquares contains code to compute the Picard ranks of both the square of squares variety V and the K3 surface X from section 4.2. 
To run the verification for Theorem 2, paste in the first few lines of PicardRankCodeSquareOfSquares.m, up until retrieving the singular points of V. Then, paste in Comps.m. Afterwards, paste the rest of the code up through the first instance of Rank(pairingmat). The computation should take ~6 hours. We note that a device running Magma natively will run out of available memory to complete the computation, so we recommend running it on a supercomputing client. As far as ComputeLinePreimages.m goes, this is to compute the 88 divisors appearing as preimages of lines on the magic del Pezzos surfaces derived in Section 4.1. The equations are already in Comps.m, but the curious verifier may check that the code gives the divisors we claim to have.

*PROPOSITION 3.9: Empty Fano Variety of Lines*

To be added

*PROPOSITION 4.2: Classification of Magic Quartic del Pezzo Surfaces*

The classification of magic quartic del Pezzo surfaces is entirely contained in delPezzoClassification. The instructions for how to run it are given in the file.

*PROPOSITION 4.5: Geometric Picard Rank of the Magic Octic K3 Surface X*

This proposition requires 3 files: MagicK3Comps.m, PicardRankCodeSquareOfSquares.m, and MagicK3WeilPolynomials.txt. Verifying the lower bound on the Picard rank of X is similar to the method for Theorem 2. After the end of the code to compute the bound on the Picard rank for V, the code for X begins. Paste it up until the singular points, and then paste MagicK3Comps.m in. Then, paste the rest of the code. The computation takes a few minutes. The Weil polynomials of X are all given in MagicK3WeilPolynomials.txt, along with extra data about the branch sextic and Frobenius traces.

*PROPOSITION 4.7: Classification of Magic Octic K3 Surfaces*

The classification of magic octic K3 surfaces is entirely contained in AssociatedDoubleCoverAndLines. We give code to compute the associated double cover of an octic K3 surface, its splitting into lines, and the number of nodes/triple points. More directions are given in the file.

Enjoy!
