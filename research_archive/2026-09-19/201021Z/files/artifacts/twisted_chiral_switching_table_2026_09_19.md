# Switching versus permutation in the finite twisted family

All matchings d are optimized in each column. Every fixed-child family comparison is exact; these are not global parent optima.

| Child | Untwisted | Switching only | Full family | Permutations indispensable for full minimum |
|---|---|---|---|---|
| order3_class0 | 5 | 5 | 5 | no |
| order4_class0 | 10 | 10 | 10 | no |
| order5_class0 | 13 | 13 | 13 | no |
| order6_class0 | 18 | 18 | 18 | no |
| order7_class0 | 25 | 25 | 21 | yes |
| order7_class1 | 25 | 25 | 21 | yes |
| order7_class2 | 25 | 25 | 25 | no |
| order8_class0 | 32 | 32 | 30 | yes |
| order8_class1 | 32 | 32 | 30 | yes |
| order9_class0 | 39 | 37 | 33 | yes |
| order9_class1 | 41 | 37 | 37 | no |
| order9_class2 | 39 | 37 | 33 | yes |
| order9_class3 | 39 | 37 | 37 | no |
| order9_class4 | 41 | 37 | 37 | no |
| order9_class5 | 39 | 39 | 33 | yes |
| order9_class6 | 39 | 37 | 37 | no |
| order9_class7 | 39 | 33 | 33 | no |
| order9_class8 | 33 | 33 | 33 | no |
| order10_class0 | 44 | 44 | 44 | no |
| order10_class1 | 44 | 44 | 44 | no |
| order10_nonoptimal_conference | 40 | 40 | 40 | no |

Switching-only means p=id and every2^(n-1) switching vector modulo global sign; every2^n matching d is tested.

At J={i:s_i=-1}, rotate (x_i,y_i)=(yprime_i,-xprime_i). The resulting child reverses only edges internal to J, the bridge becomes an untwisted copy of that child, and the matching becomes d*s. This exact signed-permutation matrix identity passed120 independently generated random examples.

All1024 partial reversals of the conference order10 child have Q histogram {15:82,17:762,19:180}; none is an optimal Q13 child.

Source/results: computations/twisted_chiral_switching_audit_2026_09_19.py and computations/results/twisted_chiral_switching_audit_2026_09_19.json.
