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
