MW = copy(M);
for j in 1:Int(log2(n))
    p = Int(n/2^(j-1));
    sel = 1:p;
    # average/difference along X
    MW[sel, sel, sel] = cat(1, (MW[1:2:p, sel, sel] + MW[2:2:p, sel, sel])./sqrt(2), (MW[1:2:p, sel, sel] - MW[2:2:p, sel, sel])./sqrt(2) );
    # average/difference along Y
    MW[sel, sel, sel] = cat(2, (MW[sel, 1:2:p, sel] + MW[sel, 2:2:p, sel])./sqrt(2), (MW[sel, 1:2:p, sel] - MW[sel, 2:2:p, sel])./sqrt(2) );
    # average/difference along Z
    MW[sel, sel, sel] = cat(3, (MW[sel, sel, 1:2:p] + MW[sel, sel, 2:2:p])./sqrt(2), (MW[sel, sel, 1:2:p] - MW[sel, sel, 2:2:p])./sqrt(2) );
end
