//This is one file to compute the Picard rank of the square of squares variety.

//We need to build our variety and its curves over Q(i, sqrt(2), sqrt(3)). We also need to have Comps.m on hand to produce this.

R<x> := PolynomialRing(Rationals());
L<i, sqrt2, sqrt3> := NumberField([x^2 + 1, x^2 - 2, x^2 - 3]);
P<a,b,c,d,m,e,f,g,h> := ProjectiveSpace(L, 8);
S := Scheme(P,[a^2+b^2+c^2-3*m^2,d^2+e^2-2*m^2,f^2+g^2+h^2-3*m^2,a^2+d^2+f^2-3*m^2,b^2-2*m^2+g^2,c^2+e^2+h^2-3*m^2,a^2-2*m^2+h^2,f^2-2*m^2+c^2]);
pts := Points(SingularSubscheme(S));
r := map<P -> P | [f, d, a, g, m, b, h, e, c]>;
s := map<P -> P | [c, b, a, e, m, d, h, g, f]>;

//Paste in Comps.m here. See Sections 3 and 4 for all of the information. Comps.m is in a state where the curves have not had their Groebner bases taken or had the base field flattened. We do this so that the reader can see all of the explicit equations and observe how the automorphism group/Galois group are acting on the divisors. In order to make sure that the Picard rank code terminates, we have to amend this, as the intersection computations will fail to terminate otherwise.

Labs := OptimizedRepresentation(AbsoluteField(L));
bool, phi := IsIsomorphic(L, Labs);
assert bool;                       // phi : L -> Labs, carries i, sqrt2, sqrt3 over

PK := ProjectiveSpace(Labs, 8);
AssignNames(~PK, ["a","b","c","d","e","f","g","h","m"]);
Rnew := CoordinateRing(PK);
Rold := Universe(Comps[1]);        // Polynomial ring the Comps live in, over L
psi  := hom< Rold -> Rnew | phi, [Rnew.k : k in [1..Rank(Rold)]] >;

//We now rebuild S over the flattened base field, and map the singular points in, and bring Comps over as well. We turn Comps into CurvesOnS, taking Groebner bases of each ideal to speed up later computations.

Snew := Scheme(PK, [psi(f) : f in DefiningEquations(S)]); 

ptsNew := [ PK ! [phi(c) : c in Eltseq(p)] : p in pts ];
CompsNew := [* [ psi(f) : f in comp ] : comp in Comps *];
CurvesOnS := [];
for i in [1..#CompsNew] do
    I := Ideal(Scheme(Snew, CompsNew[i]));
    Groebner(I);                                  // reduce; fast over Labs
    Append(~CurvesOnS, Scheme(Snew, GroebnerBasis(I)));
end for;

//Singular points on each curve, as well as precomputing genera and degrees for later computations.

Cpts  := [ {pt : pt in ptsNew | pt in C} : C in CurvesOnS ];
genus := [ ArithmeticGenus(C) : C in CurvesOnS ];
degs  := [ Degree(C)          : C in CurvesOnS ];
n     := #CurvesOnS;

//We now write the intersection function. See Section 5 for an explanation of how this is computed. We note that the formula subtracts integer multiples of the product of multiplicities at each point, instead of halving them (as is computed in the paper). This is because we are taking the degrees of intersections on the reduced subschemes downstairs, which already adjusts for these.

function intersection(C, j)
    if j le n then
        D  := CurvesOnS[j];
        gC := ArithmeticGenus(C);     // capture before any Difference mutates C
        dC := Degree(C);
        Corig := C;
        m  := 0;
        CC := C meet D;
        if Dimension(CC) eq 1 then
            while IsSubscheme(D, C) do
                m +:= 1;
                C := Difference(C, D);
            end while;
            CC := C meet D;
            assert Dimension(CC) lt 1;
        end if;
        // full node correction: product of both multiplicities, no 1/2
        return (2*gC - 2 - 3*dC)*m + Degree(CC)
               - &+[Integers() | Multiplicity(Corig, pt)*Multiplicity(D, pt)
                                 : pt in Cpts[j] | pt in Corig];
    else
        pt := ptsNew[j - n];
        return pt in C select Multiplicity(C, pt) else 0;
    end if;
end function;

//We now build the upper triangle of the curve pairings, outputting to a separate file to optimize memory/CPU usage.

outname := "MatCC_final.txt";

// It's important to ensure that the file starts empty, or else the Picard rank returned will not be accurate.

SetOutputFile(outname : Overwrite := true); UnsetOutputFile();

for j in [1..n] do
    printf "row %o\n", j;
    Cj := CurvesOnS[j];
    for k in [j..n] do
        if k eq j then
            v := 2*genus[j] - 2 - 3*degs[j];      // diagonal: clean formula
        else
            v := intersection(Cj, k);
        end if;
        fprintf outname, "%o %o %o\n", j, k, v;
    end for;
end for;

//We now assemble the full pairing matrix, starting with the curve pairings, then adding the pairings of the curves with the exceptional divisors and the -2 block of exceptional divisors in separate blocks

bdim := n + #ptsNew;
pairingmat := ZeroMatrix(Integers(), bdim, bdim);

for line in Split(Read(outname), "\n") do
    t := Split(line, " ");
    if #t eq 3 then
        j := StringToInteger(t[1]);
        k := StringToInteger(t[2]);
        v := StringToInteger(t[3]);
        pairingmat[j,k] := v;
        pairingmat[k,j] := v;
    end if;
end for;

// curve <-> exceptional (node) coupling: C~.E = mult (= 1)
MatCP := Matrix(Integers(),
    [ [ ptsNew[k] in ptC select 1 else 0 : k in [1..#ptsNew] ]
      where ptC := Cpts[j] : j in [1..n] ]);
InsertBlock(~pairingmat, MatCP, 1, n+1);
InsertBlock(~pairingmat, Transpose(MatCP), n+1, 1);
// exceptional curves: pairwise disjoint, self-intersection -2
InsertBlock(~pairingmat, DiagonalMatrix(Integers(), [-2 : j in [1..#ptsNew]]), n+1, n+1);

//Check the rank. It should be 518!

Rank(pairingmat);


//This is the code to compute the Picard rank of the Magic K3. Our cycles are defined over Q(i, sqrt(2), sqrt(3)), so we need to define our variety over it. The code is essentially as above, but we need to change the self-intersection formula to line up with the formula for K3 surfaces C^2 = 2g - 2

R<x> := PolynomialRing(Rationals());
L<i, sqrt2, sqrt3> := NumberField([x^2 + 1, x^2 -2, x^2 -3]);
P<a,b,c,m,f,g> := ProjectiveSpace(L, 5);
S := Scheme(P, [a^2 + b^2 + c^2 - 3*m^2, b^2 + g^2 - 2*m^2, c^2 + f^2 - 2*m^2]);
pts := Points(SingularSubscheme(S));

//Copy in MagicK3Comps.m

CurvesOnS := [];
for i in [1..#Comps] do
Append(~CurvesOnS, Scheme(S, Comps[i]));
end for;

Cpts := [{pt : pt in pts | pt in C} : C in CurvesOnS];
//Now that we have imported the curves list, we can create the intersection function

function intersection(C, j)
  // j is index in "CurvesOnS cat pts"
  if j le #CurvesOnS then
    m := 0;
    CC := C meet CurvesOnS[j];
    if Dimension(CC) eq 1 then
      // CurvesOnS[j] is contained in C; find the multiplicity m
      while IsSubscheme(CurvesOnS[j], C) do
        m +:= 1;
        C := Difference(C, CurvesOnS[j]);
      end while;
      CC := C meet CurvesOnS[j];
      assert Dimension(CC) lt 1;
    end if;
//The self-intersection of CurvesOnS[j] is 2g -2 --> (2g - 2)*m
    //The count of intersection points (with multiplicity) --> Degree(CC)
    //Subtract the number of singularities of S-bar among these
    //(blowing them up reduces the intersection number there) --> -&+[...]
    return (2*ArithmeticGenus(C) - 2)*m + Degree(CC) - &+[Integers() | Multiplicity(C, pt) : pt in Cpts[j] | pt in C];
else
    pt := pts[j-#CurvesOnS];
    return pt in C select Multiplicity(C, pt) else 0;
  end if;
end function;

bdim := #CurvesOnS + #pts;
pairingmat := ZeroMatrix(Integers(), bdim, bdim); // initialize the matrix with zeros
//Pairing between curves
MatCC := Matrix(Integers(),
                [[k eq j select 2*ArithmeticGenus(C) - 2 else intersection(C, k) : k in [1..#CurvesOnS]] where C := CurvesOnS[j]
                   : j in [1..#CurvesOnS]]);
//Pairing between curves and singularities
MatCP := Matrix(Integers(),
                [[pts[k] in ptC select 1 else 0 : k in [1..#pts]] where ptC := Cpts[j]
                   : j in [1..#CurvesOnS]]);
//Put the parts together
InsertBlock(~pairingmat, MatCC, 1, 1);
InsertBlock(~pairingmat, MatCP, 1, #CurvesOnS+1);
InsertBlock(~pairingmat, Transpose(MatCP), #CurvesOnS+1, 1);
//The exceptional curves are pairwise disjoint and have self-intersection -2
InsertBlock(~pairingmat, DiagonalMatrix(Integers(), [-2 : j in [1..#pts]]), #CurvesOnS+1, #CurvesOnS+1);

//The rank is 19!
