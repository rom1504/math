# Exact finite audit of the actual energy-shell retuning split

Status: **complete actual-child enumeration at balanced `4+4` and held-out
`4+5`; finite numerical evidence only**.  The experiment verifies the exact
KL identities in Theorem 37.62 and measures which term carries the averaged
latent-posterior retuning.  It uses contracted-temperature pressure
minimizers selected by exhaustive signing enumeration, not conference or
Paley surrogates.

Reproducible files:

- [`actual_child_energy_shell_retuning.py`](actual_child_energy_shell_retuning.py)
- [`../../computations/results/actual_child_energy_shell_retuning.json`](../../computations/results/actual_child_energy_shell_retuning.json)
- [`../../computations/results/actual_child_energy_shell_retuning_n9.json`](../../computations/results/actual_child_energy_shell_retuning_n9.json)
- [`../../computations/logs/actual_child_energy_shell_retuning.log`](../../computations/logs/actual_child_energy_shell_retuning.log)
- [`../../computations/logs/actual_child_energy_shell_retuning_n9.log`](../../computations/logs/actual_child_energy_shell_retuning_n9.log)

## 1. Evaluated law and identities

For every selected child pair and orientation, the program enumerates all
signed rank-one words `Q`, their exact integer combined-energy labels

```math
E_\epsilon(Q)=H_A(x)+\epsilon H_D(y),
```

and all `2^(mn)` bridges.  It constructs the exact binary-channel likelihood
`p(B)`, the inverse bridge law `q_a(B) proportional p(B)^a` at
`a=-lambda`, and the averaged ordinary forward posterior

```math
\bar\mu(Q)=\sum_Bq_a(B)\mu(Q\mid B).
```

It then evaluates

```math
\begin{aligned}
D_{\rm total}
 &=D(\bar\mu\Vert\mu),\\
D_{\rm shell}
 &=D(\bar p_E\Vert p_E),\\
D_{\rm geom}
 &=\sum_e\bar p_eD(\bar\mu(\cdot\mid e)
                   \Vert U_{\mathcal S_e}).
\end{aligned}
```

The exact chain rule is `D_total=D_shell+D_geom`.  A separate pass evaluates
the posterior mutual informations `I(B;Q)` and `I(B;E)` and checks

```math
E_qD(\mu(\cdot\mid B)\Vert\mu)
=I(B;Q)+D(\bar\mu\Vert\mu).
```

Across every record, the largest absolute residual in the shell identity is
`1.4e-15` and the largest posterior-budget residual is `6.7e-15`.

## 2. Results

The table gives one orientation; the other agrees to displayed precision.
All tested priors have nine nonempty combined-energy shells.

| split | `beta` | `lambda` | `D_total` | `D_shell` | `D_geom` | shell fraction |
|---:|---:|---:|---:|---:|---:|---:|
| `4+4` | 1 | 1 | 0.001594 | 0.001520 | 0.000074 | 0.9534 |
| `4+4` | 1 | 5.3821 | 0.008181 | 0.007693 | 0.000488 | 0.9403 |
| `4+4` | 2 | 1 | 0.170387 | 0.166376 | 0.004012 | 0.9765 |
| `4+4` | 2 | 5.3821 | 0.380210 | 0.362285 | 0.017925 | 0.9529 |
| `4+4` | 4 | 1 | 1.701857 | 1.676889 | 0.024968 | 0.9853 |
| `4+4` | 4 | 5.3821 | 2.621118 | 2.520576 | 0.100541 | 0.9616 |
| `4+5` | 1 | 1 | 0.001071 | 0.000991 | 0.000079 | 0.9258 |
| `4+5` | 1 | 5.3821 | 0.006020 | 0.005488 | 0.000531 | 0.9117 |
| `4+5` | 2 | 1 | 0.119096 | 0.114376 | 0.004720 | 0.9604 |
| `4+5` | 2 | 5.3821 | 0.289268 | 0.267856 | 0.021411 | 0.9260 |
| `4+5` | 4 | 1 | 1.261404 | 1.240518 | 0.020886 | 0.9834 |
| `4+5` | 4 | 5.3821 | 1.928438 | 1.862542 | 0.065895 | 0.9658 |

Thus the scalar shell marginal carries between `91.17%` and `98.53%` of
the finite retuning in this audit.  This is a useful finite target for the
rare-shell branch.  It is not an asymptotic theorem and does not connect the
shell KL to the canonical row lifetime or the optimal reverse-product
dependence.

## 3. Exact scope

The signing cubes, child energies, rank-one words, shell labels, and bridge
cubes are enumerated exactly.  Selection between transcendental child
pressures uses 80-decimal `mpmath`; Gibbs probabilities and information
quantities use double precision.  The `4+5` computation is held out from the
balanced derivation and contains `2^20` bridges.

The observed dominance of `D_shell` does **not** prove that a shell state is
dynamically closed.  Computing the retuned shell marginal still uses the
complete bridge likelihood, and coherent row-factor retuning can occur even
when the averaged latent posterior is unchanged.  The experiment therefore
supports investigating a one-dimensional rare-event mechanism but supplies
no Level-6 recurrence.
