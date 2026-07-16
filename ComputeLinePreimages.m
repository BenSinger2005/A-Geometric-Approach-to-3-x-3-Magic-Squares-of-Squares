//Here is the code to compute the line preimages on V for the 23 del Pezzo coordinate projections. The idea is to do the following:
//(1) Build the 8x9 matrix of the magic relations in the squares;
//(2) For each configurations, derive the 2 diagonal quadrics of the del Pezzo image as the left-nullspace of the complementary columns, read off on the 5 chosen columns;
//(3) find every line on  Q1 = Q2 = 0  in P^4 by the graph-chart method over all 10 pivot pairs (solve the 6-equation system, lines are 2-planes by their echelon form);
//(4) pull each line back to V:  Veqns cat [ 3 linear forms ].
// We work over K = Q(i, sqrt2, sqrt3, sqrt5).

Qx<x> := PolynomialRing(Rationals());
K<i,sqrt2,sqrt3,sqrt5> := NumberField([x^2+1, x^2-2, x^2-3, x^2-5]);

// Ambient coordinate ring of P^8 and the defining equations of V.
P<a,b,c,d,m,e,f,g,h> := PolynomialRing(K, 9);
vars := [a,b,c,d,m,e,f,g,h];
names := ["a","b","c","d","m","e","f","g","h"];
idx := AssociativeArray();
for j in [1..9] do idx[names[j]] := j; end for;

Veqns := [
    a^2 + b^2 + c^2 - 3*m^2,
    d^2 + e^2 - 2*m^2,
    f^2 + g^2 + h^2 - 3*m^2,
    a^2 + d^2 + f^2 - 3*m^2,
    b^2 + g^2 - 2*m^2,
    c^2 + e^2 + h^2 - 3*m^2,
    a^2 + h^2 - 2*m^2,
    c^2 + f^2 - 2*m^2
];

//Here is the 8x9 matrix in the quadratic relations for V.
M := ZeroMatrix(Rationals(), 8, 9);
rels := [
    [<"a",1>,<"b",1>,<"c",1>,<"m",-3>],   // row 1
    [<"d",1>,<"e",1>,<"m",-2>],           // row 2
    [<"f",1>,<"g",1>,<"h",1>,<"m",-3>],   // row 3
    [<"a",1>,<"d",1>,<"f",1>,<"m",-3>],   // col 1
    [<"b",1>,<"g",1>,<"m",-2>],           // col 2
    [<"c",1>,<"e",1>,<"h",1>,<"m",-3>],   // col 3
    [<"a",1>,<"h",1>,<"m",-2>],           // main diagonal
    [<"c",1>,<"f",1>,<"m",-2>]            // anti-diagonal
];
for r in [1..8] do
    for t in rels[r] do M[r, idx[t[1]]] := t[2]; end for;
end for;

//Here are the two diagonal quadrics of a configuration.
//We choose the 5 column indices (in 1..9) of the configuration's variables.
DPQuadrics := function(chosen)
    comp    := [j : j in [1..9] | j notin chosen];
    Mcomp   := Submatrix(M, [1..8], comp);     // 8 x 4
    Mchosen := Submatrix(M, [1..8], chosen);   // 8 x 5
    // left-nullspace rows c with c * Mcomp = 0, then evaluate on chosen columns:
    Kc   := KernelMatrix(Mcomp);               // rows span {c : c*Mcomp = 0}
    Qmat := Kc * Mchosen;                       // (#rels) x 5 diagonal-quadric coeffs
    E    := EchelonForm(Qmat);
    rows := [ [E[r,t] : t in [1..5]] : r in [1..Nrows(E)] | not IsZero(E[r]) ];
    error if #rows ne 2, "expected exactly 2 quadrics";
    return rows[1], rows[2];
end function;

//We now find lines on Q1 = Q2 = 0 in P^4.
FindLines := function(q1, q2)
    aa := [K!q1[t] : t in [1..5]];
    bb := [K!q2[t] : t in [1..5]];
    R  := PolynomialRing(K, 6);
    AssignNames(~R, ["al1","al2","al3","be1","be2","be3"]);
    al := [R.1, R.2, R.3];  be := [R.4, R.5, R.6];
    found := {};            
    lines := [];            
    for p in [1..5] do
      for q in [p+1..5] do
        J := [t : t in [1..5] | t ne p and t ne q];     
        // line is the graph  x_{J[k]} = al[k]*x_p + be[k]*x_q ; impose both quadrics:
        eqs := [
            aa[p] + &+[ aa[J[k]]*al[k]^2      : k in [1..3] ],
                    &+[ aa[J[k]]*al[k]*be[k]  : k in [1..3] ],
            aa[q] + &+[ aa[J[k]]*be[k]^2      : k in [1..3] ],
            bb[p] + &+[ bb[J[k]]*al[k]^2      : k in [1..3] ],
                    &+[ bb[J[k]]*al[k]*be[k]  : k in [1..3] ],
            bb[q] + &+[ bb[J[k]]*be[k]^2      : k in [1..3] ]
        ];
        I := ideal<R | eqs>;
        if Dimension(I) gt 0 then continue; end if;
        for pt in Variety(I) do
            alpha := [pt[1], pt[2], pt[3]];
            beta  := [pt[4], pt[5], pt[6]];
            // canonical key: echelon form of the 2x5 plane matrix
            Pl := ZeroMatrix(K, 2, 5);
            Pl[1,p] := 1;  Pl[2,q] := 1;
            for k in [1..3] do Pl[1,J[k]] := alpha[k]; Pl[2,J[k]] := beta[k]; end for;
            key := Eltseq(EchelonForm(Pl));
            if key notin found then
                Include(~found, key);
                Append(~lines, <p, q, J, alpha, beta>);
            end if;
        end for;
      end for;
    end for;
    return lines;
end function;

//Finally, we pull the line back to 3 linear forms.
ClearDen := function(form)
    den := 1;
    for cc in Coefficients(form) do den := Lcm(den, Denominator(cc)); end for;
    return den * form;
end function;

LineForms := function(cfg, line)
    p := line[1]; q := line[2]; J := line[3]; alpha := line[4]; beta := line[5];
    vp := vars[idx[cfg[p]]];  vq := vars[idx[cfg[q]]];
    return [ ClearDen( vars[idx[cfg[J[k]]]] - alpha[k]*vp - beta[k]*vq ) : k in [1..3] ];
end function;

//We now compute smoothness via finding the number of distinct pencil roots
Kind := function(q1, q2)
    rs := {};
    for t in [1..5] do
        if q2[t] eq 0 then Include(~rs, "inf");
        else Include(~rs, Sprint(q1[t]/q2[t])); end if;
    end for;
    nd := #rs;
    if nd eq 5 then return "smooth dP4 (16 lines)"; end if;
    if nd eq 4 then return "1 node (8 lines)"; end if;
    return "2 nodes (4 lines)";
end function;

//Here are the 23 configurations.
Configs := [
    ["a","b","c","d","m"], ["a","b","c","d","e"], ["a","b","c","d","f"],
    ["a","b","c","d","g"], ["a","b","c","d","h"], ["a","b","c","m","f"],
    ["a","b","c","m","g"], ["a","b","c","f","g"], ["a","b","c","f","h"],
    ["a","b","d","m","e"], ["a","b","d","m","h"], ["a","b","d","e","g"],
    ["a","b","d","e","h"], ["a","b","m","e","f"], ["a","b","m","e","g"],
    ["a","b","m","e","h"], ["a","b","m","f","g"], ["a","b","m","f","h"],
    ["a","b","m","g","h"], ["a","b","m","e","h"], ["a","b","e","f","h"],
    ["a","c","m","f","h"], ["b","d","m","e","g"]
];

//We now print the Veqns cat [...] blocks.
for n in [1..#Configs] do
    cfg    := Configs[n];
    chosen := [ idx[s] : s in cfg ];
    q1, q2 := DPQuadrics(chosen);
    lines  := FindLines(q1, q2);

    cfgstr := cfg[1];
    for t in [2..5] do cfgstr := cfgstr cat ", " cat cfg[t]; end for;
    printf "// ===== Config%o := [%o] -> %o, %o lines =====\n", n, cfgstr, Kind(q1,q2), #lines;
    printf "Config%o_lines := [\n", n;
    for t in [1..#lines] do
        fs  := LineForms(cfg, lines[t]);
        sep := (t lt #lines) select "," else "";
        printf "    Veqns cat [ %o,  %o,  %o ]%o   // Line %o\n", fs[1], fs[2], fs[3], sep, t;
    end for;
    printf "];\n\n";
end for;

//To get the rest of the divisors corresponding to the orbits of the line preimages in each configuration under the action of the automorphism group, add maps corresponding to the rotation and reflection in the automorphism group as outlined in Section 1. Then, loop them over an appended list of all preimages computed as above, adding the orbit of each line preimage. Most of them are equal to components we already have; the new ones take up indices 417..800 in Comps.m.

//Here is the code to compute the preimages of lines on cubic surfaces. The idea is to project each magic dP4 from the point p = [1,1,1,1,1] which is present on all of them, find the lines on all of the cubic surfaces, and then pull back to V.

//We first give the function to compute the projection from a point given the two defining quadrics of a magic dP4.

GetCubic := function(q1, q2)
    A := PolynomialRing(K, 9);          // A.1..A.5 = y (chosen coords), A.6..A.9 = u
    y := [A.t     : t in [1..5]];
    u := [A.(5+t) : t in [1..4]];
    Q1 := &+[ q1[t]*y[t]^2 : t in [1..5] ];
    Q2 := &+[ q2[t]*y[t]^2 : t in [1..5] ];
    proj := [ y[1]-y[2], y[2]-y[3], y[3]-y[4], y[4]-y[5] ];
    I  := ideal<A | Q1, Q2, [ u[t]-proj[t] : t in [1..4] ]>;
    Iu := EliminationIdeal(I, 5);        // eliminate the 5 y's -> ideal in u only
    cub := [b : b in Basis(Iu) | TotalDegree(b) eq 3];
    error if #cub eq 0, "no cubic generator found (p may be a node here)";
    // move the cubic into a clean 4-variable ring
    Ru := PolynomialRing(K, 4);
    toRu := hom< A -> Ru | [Ru | 0,0,0,0,0, Ru.1,Ru.2,Ru.3,Ru.4] >;
    return toRu(cub[1]), Ru;
end function;

// Here is the function to compute the lines on a given cubic surface, following the Elsenhans-Jahnel method similarly to above.

LinesOnCubic := function(F, Ru)
    R := PolynomialRing(K, 4);
    AssignNames(~R, ["al1","al2","be1","be2"]);
    al := [R.1, R.2];  be := [R.3, R.4];
    Rst := PolynomialRing(R, 2);  s := Rst.1;  t := Rst.2;
    found := {};  lines := [];
    for p in [1..4] do
      for q in [p+1..4] do
        J  := [j : j in [1..4] | j ne p and j ne q];
        pt := [Rst| 0,0,0,0];
        pt[p] := s;  pt[q] := t;
        pt[J[1]] := al[1]*s + be[1]*t;
        pt[J[2]] := al[2]*s + be[2]*t;
        Fval := Evaluate(F, pt);                      
        eqs  := [ MonomialCoefficient(Fval, m)
                    : m in [s^3, s^2*t, s*t^2, t^3] ];
        I := ideal<R | eqs>;
        if Dimension(I) gt 0 then continue; end if;
        for sol in Variety(I) do
            a1,a2,b1,b2 := Explode([sol[1],sol[2],sol[3],sol[4]]);
            Pl := ZeroMatrix(K, 2, 4);
            Pl[1,p]:=1; Pl[2,q]:=1;
            Pl[1,J[1]]:=a1; Pl[1,J[2]]:=a2; Pl[2,J[1]]:=b1; Pl[2,J[2]]:=b2;
            key := Eltseq(EchelonForm(Pl));
            if key notin found then
                Include(~found, key);
                Append(~lines, <p, q, J, [a1,a2], [b1,b2]>);
            end if;
        end for;
      end for;
    end for;
    return lines;
end function;

//  We now give the function to compute the preimages of each line. The u-coordinates in the original variables are: uu[t] = C_t - C_{t+1}, C_t = cfg's t-th var.

CubicLineForms := function(cfg, line)
    C  := [ vars[idx[cfg[t]]] : t in [1..5] ];
    uu := [ C[1]-C[2], C[2]-C[3], C[3]-C[4], C[4]-C[5] ];  // pullback of u0..u3
    p := line[1]; q := line[2]; J := line[3];
    a := line[4]; b := line[5];
    L1 := uu[J[1]] - a[1]*uu[p] - b[1]*uu[q];
    L2 := uu[J[2]] - a[2]*uu[p] - b[2]*uu[q];
    return [ ClearDen(L1), ClearDen(L2) ];
end function;

// We now run this over all the configurations, omitting the pullback blocks.

for n in [1..#Configs] do
    cfg    := Configs[n];
    chosen := [ idx[st] : st in cfg ];
    q1, q2 := DPQuadrics(chosen);
    F, Ru  := GetCubic(q1, q2);
    lines  := LinesOnCubic(F, Ru);

    cfgstr := cfg[1]; for st in [2..5] do cfgstr := cfgstr cat ", " cat cfg[st]; end for;
    printf "// ===== CubicConfig%o := [%o]  ->  %o lines on the cubic =====\n",
           n, cfgstr, #lines;
    printf "CubicConfig%o_pullbacks := [\n", n;
    for tt in [1..#lines] do
        L := CubicLineForms(cfg, lines[tt]);
        sep := (tt lt #lines) select "," else "";
        printf "    Veqns cat [ %o,  %o ]%o   // line %o\n", L[1], L[2], sep, tt;
    end for;
    printf "];\n\n";
end for;

PP := Proj(P);   // P^8 with coordinate ring P<a,b,c,d,m,e,f,g,h>

// 1-dimensional irreducible components of a divisor block, with canonical keys.
CurveComponents := function(eqns)
    out := [];   // list of <key, scheme>
    for C in IrreducibleComponents(Scheme(PP, eqns)) do
        if Dimension(C) eq 1 then
            Append(~out, < GroebnerBasis(Ideal(C)), C >);
        end if;
    end for;
    return out;
end function;

//Building the list of cubics. 

CubicLists := [* *];    
for n in [1..#Configs] do
    cfg    := Configs[n];
    chosen := [ idx[st] : st in cfg ];
    q1, q2 := DPQuadrics(chosen);
    F, Ru  := GetCubic(q1, q2);
    lines  := LinesOnCubic(F, Ru);

    blocks := [ Veqns cat CubicLineForms(cfg, lines[tt]) : tt in [1..#lines] ];
    Append(~CubicLists, blocks);

end for;

//Output it! To get the rest of the divisors corresponding to the orbits of the line preimages in each configuration under the action of the automorphism group, add maps corresponding to the rotation and reflection in the automorphism group as outlined in Section 1. Then, loop them over CubicLists, adding the orbit of each line preimage. Most of them are equal to components we already have; the new ones take up indices 801..948 in Comps.m.

CubicLists;
