//Below is the code to verify that the Fano scheme of lines on the 3 \times 3 magic square of squares variety V is empty. Follow Appendix A, we use a method of Elsenhans--Jahnel, enumerating through Schubert cells of the Grassmannian G(1, 8) to show that V has no lines.

Qt := Rationals();

/* We take the coefficients of the six defining quadrics of V in the coordinate order
        1  2  3  4  5  6  7  8  9
        A  B  C  D  M  E  F  G  H
   All six are diagonal, which is what makes the expansion below
   collapse to three bilinear expressions per quadric. */

coeffs := [
  [ 3, 0, 0, 0, 0, 0, -2, -2,  1],   //  3A^2 - 2F^2 - 2G^2 +  H^2
  [ 0, 3, 0, 0, 0, 0, -2,  1, -2],   //  3B^2 - 2F^2 +  G^2 - 2H^2
  [ 0, 0, 3, 0, 0, 0,  1, -2, -2],   //  3C^2 +  F^2 - 2G^2 - 2H^2
  [ 0, 0, 0, 3, 0, 0,  2, -1, -4],   //  3D^2 + 2F^2 -  G^2 - 4H^2
  [ 0, 0, 0, 0, 0, 3, -4, -1,  2],   //  3E^2 - 4F^2 -  G^2 + 2H^2
  [ 0, 0, 0, 0, 3, 0, -1, -1, -1]    //  3M^2 -  F^2 -  G^2 -  H^2
];

//We now make one ambient chart ring with 14 variables, serving every cell. Cells of smaller dimension leave the trailing variables unused, which does not affect emptiness. 

A := PolynomialRing(Qt, 14);
AssignNames(~A, ["a" cat IntegerToString(i) : i in [1..14]]);

//We now generate the affine chart on the Schubert cell with pivots (p,q).

Chart := function(p, q)
    u := [ A | 0 : i in [1..9] ];
    v := [ A | 0 : i in [1..9] ];
    u[p] := 1;
    v[q] := 1;
    n := 0;
    for i in [p+1..9] do
        if i ne q then
            n +:= 1; u[i] := A.n;
        end if;
    end for;
    for i in [q+1..9] do
        n +:= 1; v[i] := A.n;
    end for;
    assert n eq 17 - p - q;
    return u, v, n;
end function;

//The three conditions imposed by one diagonal quadric:

Conds := function(c, u, v)
    return [ &+[ c[i]*u[i]^2       : i in [1..9] ],   // coeff of s^2
             &+[ c[i]*u[i]*v[i]    : i in [1..9] ],   // coeff of s*t
             &+[ c[i]*v[i]^2       : i in [1..9] ] ]; // coeff of t^2
end function;

//We check that a single quadric in P^8 does contain lines as a sanity check.

u, v, n := Chart(1, 2);
Ictl := ideal< A | Conds(coeffs[1], u, v) >;
assert not (A!1 in Ictl);
printf "control: single quadric on the open cell is nonempty, as expected\n\n";

//We now do the main computation over all 36 Schubert cells.

cells := [ [p,q] : p in [1..9], q in [1..9] | p lt q ];
assert #cells eq 36;

dims := [];
for cl in cells do
    p := cl[1]; q := cl[2];
    u, v, n := Chart(p, q);
    eqs := &cat[ Conds(c, u, v) : c in coeffs ];
    assert #eqs eq 18;

    I := ideal< A | eqs >;
    G := GroebnerBasis(I);

    assert A!1 in I;
    Append(~dims, n);
    printf "cell (%o,%o)  dim %2o  #eqs %o  Groebner basis %o  -> empty\n",
              p, q, n, #eqs, G;
end for;

assert Sort(dims) eq Sort([ 17-cl[1]-cl[2] : cl in cells ]);
assert Maximum(dims) eq 14 and Minimum(dims) eq 0;

printf "\nF(V) meets no Schubert cell: F(V) is empty, so V contains no lines.\n";