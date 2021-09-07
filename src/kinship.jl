"""
    function kinship(ped, i, j)
---
This function is handy if just to calculate relationship of a few (pairs of) ID.
It can also speed up by adding `Thread.@threads` before your pair loop.
"""
function kinship(ped, i, j)
    (i == 0 || j == 0) && return 0
    ipa, ima = ped[i, :]          # used both for below and the last
    i == j && (return 1 + .5kinship(ped, ipa, ima))
    if i < j
        jpa, jma = ped[j, :]
        return .5(kinship(ped, i, jpa) + kinship(ped, i, jma))
    end
    return .5(kinship(ped, j, ipa) + kinship(ped, j, ima))
end
