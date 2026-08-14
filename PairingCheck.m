//This is the code to verify the claims about nonvanishing 2 x 2 minors of the matrix N = (q_{jk}) in Proposition 4.11. Given the smooth octic magic K3 surface Y corresponding to orbit 2, and each of the fifteen pairings of its six coordinates, we exhibit a 2 x 2 minor of N that does not vanish on Y, and we check that the locus {rank(N) \leq 1} is a proper closed subset of Y.

Qt := Rationals();
PP<x0,x1,x2,x3,x4,x5> := ProjectiveSpace(Qt, 5);
R := CoordinateRing(PP);
X := [R.i : i in [1..6]];

//This is the coefficient matrix of the three diagonal quadrics defining Y, in the coordinates [A:B:C:F:G:H] = [x0:x1:x2:x3:x4:x5].
 
Lambda := Matrix(Qt, 3, 6, [
    3,  0,  0, -2, -2,  1,
    0,  3,  0, -2,  1, -2,
    0,  0,  3,  1, -2, -2 ]);

Qs := [ &+[ Lambda[j][i]*X[i]^2 : i in [1..6] ] : j in [1..3] ];
Y  := Scheme(PP, Qs);
IY := Ideal(Y);

assert Dimension(Y) eq 2;
assert IsIrreducible(Y);
assert IsNonsingular(Y);
printf "Y: smooth irreducible surface, degree %o\n\n", Degree(Y);

//Here are the partitions of {1..6} into three unordered pairs.

pairings := [];
for a in [2..6] do
    p1   := [1, a];
    rest := [ i : i in [2..6] | i ne a ];
    b    := rest[1];
    for c in [2..4] do
        p2 := [b, rest[c]];
        p3 := [ i : i in rest | i notin p2 ];
        Append(~pairings, [p1, p2, p3]);
    end for;
end for;
assert #pairings eq 15;

//Here is the matrix N = (q_{jk}) attached to a given pairing.

Nmatrix := function(pr)
    ent := [];
    for j in [1..3] do
        for k in [1..3] do
            Append(~ent, &+[ Lambda[j][i]*X[i]^2 : i in pr[k] ]);
        end for;
    end for;
    return Matrix(R, 3, 3, ent);
end function;

idx2   := [[1,2],[1,3],[2,3]];
Minors := function(N)
    return [ N[rr[1]][cc[1]]*N[rr[2]][cc[2]]
             - N[rr[1]][cc[2]]*N[rr[2]][cc[1]]
             : rr in idx2, cc in idx2 ];
end function;

//Now for the main verification.

for n in [1..#pairings] do
    pr := pairings[n];
    N  := Nmatrix(pr);
    ms := Minors(N);

//We verify that every row of N sums to Q_j, so rank(N) \leq 2 on Y.
    for j in [1..3] do
        assert &+[ N[j][k] : k in [1..3] ] eq Qs[j];
    end for;

//Here is a 2 x 2 minor not vanishing on Y.
    witness := 0;
    for m in ms do
        if m ne 0 and not IsInRadical(m, IY) then
            witness := m; break;
        end if;
    end for;
    assert witness ne 0;

//Now we verify that the rank <= 1 locus is a proper closed subset of Y.
    Zbad := Scheme(PP, Qs cat ms);
    assert Dimension(Zbad) le 1;

    printf "pairing %2o  %o\n", n, pr;
    printf "   witness minor : %o\n", witness;
    printf "   dim{rank<=1}  : %o\n", Dimension(Zbad);
end for;

printf "\nAll fifteen pairings: rank N = 2 generically on Y.\n";