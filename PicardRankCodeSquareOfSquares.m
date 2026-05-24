//This is the code to compute the Picard rank of the square of squares variety.

//We need to build our variety and its curves over Q(i, sqrt(2), sqrt(3)). We also need to have Comps.m on hand to produce this

R<x> := PolynomialRing(Rationals());
L<i, sqrt2, sqrt3> := NumberField([x^2 + 1, x^2 - 2, x^2 - 3]);
P<a,b,c,d,m,e,f,g,h> := ProjectiveSpace(L, 8);
S := Scheme(P,[a^2+b^2+c^2-3*m^2,d^2+e^2-2*m^2,f^2+g^2+h^2-3*m^2,a^2+d^2+f^2-3*m^2,b^2-2*m^2+g^2,c^2+e^2+h^2-3*m^2,a^2-2*m^2+h^2,f^2-2*m^2+c^2]);

//Get our 256 singular points

pts := Points(SingularSubscheme(S));

//Paste in Comps.m. See Section 4 for all of the information

//Create an empty curves list

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
//The self-intersection of CurvesOnS[j] is 2g -2 - 3deg(C) --> (2g - 2 - 3deg(C))*m
    //The count of intersection points (with multiplicity) --> Degree(CC)
    //Subtract the number of singularities of S-bar among these. Blowing them up reduces the intersection number there) --> -&+[...]
    return (2*ArithmeticGenus(C) - 2 - 3*Degree(C))*m + Degree(CC) - &+[Integers() | Multiplicity(C, pt) : pt in Cpts[j] | pt in C];
else
    pt := pts[j-#CurvesOnS];
    return pt in C select Multiplicity(C, pt) else 0;
  end if;
end function;

bdim := #CurvesOnS + #pts;
pairingmat := ZeroMatrix(Integers(), bdim, bdim); // initialize the matrix with zeros
// pairing between curves
MatCC := Matrix(Integers(),
                [[k eq j select 2*ArithmeticGenus(C) - 2 - 3*Degree(C) else intersection(C, k) : k in [1..#CurvesOnS]] where C := CurvesOnS[j]
                   : j in [1..#CurvesOnS]]);
//Pairing between curves and singularities
MatCP := Matrix(Integers(),
                [[pts[k] in ptC select 1 else 0 : k in [1..#pts]] where ptC := Cpts[j]
                   : j in [1..#CurvesOnS]]);
//Put the parts together
InsertBlock(~pairingmat, MatCC, 1, 1);
InsertBlock(~pairingmat, MatCP, 1, #CurvesOnS+1);
InsertBlock(~pairingmat, Transpose(MatCP), #CurvesOnS+1, 1);
// The exceptional curves are pairwise disjoint and have self-intersection -2
InsertBlock(~pairingmat, DiagonalMatrix(Integers(), [-2 : j in [1..#pts]]), #CurvesOnS+1, #CurvesOnS+1);

//Rank 410 with all components from our Comps file

//The Discriminant: 1192918290276309693519536394422003105074107558973261489180881139895427332790365
3010028225067702372679859379035882355291160065123032782938736332796386587303424
8678330728800058165686009498914845137721557511777653208654461270968000440894114
6565734523931701549451270224820016439751619068663198676599520828331054086951115
1616000000000000000000000000000000000000000. This factors as 2^762*3^101*5^39*7^5*13^20*17^4*19*47^4*61^4*5105753


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
