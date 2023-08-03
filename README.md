# Introduction

In this repository, an adaptive version of the restarted GMRES (GMRES(m)) is introduced for the resolution of the ﬁnite difference approximation of the Helmholtz equation. It has been observed that the choice of the restart parameter m strongly affects the convergence of standard GMRES(m). To overcome this problem, the GMRES(m) is formulated as a control problem in order to adaptively combine two strategies: a) the appropriate variation of the restarted parameter m, if a stagnation in the convergence is detected; and b) the augmentation of the search subspace using vectors obtained at previous cycles. The proposal is compared with similar iterative methods of the literature based on standard GMRES(m) with ﬁxed parameters. Numerical results for selected matrices suggest that the switching adaptive proposal method could overcome the stagnation observed in standard methods, and even improve the performance in terms of computational time and memory requirements.

## Instructions

- Run `mi_script_principal.m` for a simple comparative test (compiutational time, number of iterations) between different iterative methods. Line 4: `load 'cavity07.mat';` and 6: `Name_Matrix = 'cavity07';` are suitable for changes (see available `.mat` data).
- Tweak `sistema-lineal.m` to change the discretization parameters.

## Reference

ESPÍNOLA, G., Cabral, J. C., & SCHAERER, C. (2020). Adaptive GMRES (m) for the Electromagnetic Scattering Problem. TEMA (São Carlos), 21, 191-208.

## Online

    https://tema.sbmac.org.br/tema/article/view/1316
