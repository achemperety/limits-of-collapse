# Referee report (Phase I) and response map

*September 30, 2026. Numbers refer to the revised paper (siamart version).*

# Summary and recommendation

The submission assembles classical limits of minimax evaluation and adds a genuine new layer: an exact resolvent formula for the value, its collapse to matching polynomials for undirected geography, the directed resolvent collapse on SCC-symmetric digraphs, and in the addendum a substantial body of results on one-way arcs and set potentials. Most proofs I checked are correct. The weaknesses are of three kinds. First, several statements are broader than their proofs (the band statement for the rescued triangle, the treewidth dichotomy, the “free values” obstacle of §11.7). Second, the organizing conjecture is ambiguous in its field: the definition in §6.5 fixes $\mathbb{Q}(z)$, the addendum silently extends the question to $\overline{\mathbb{Q}(z)}$, and the addendum’s open problem 2 suggests that the complex form is true. It is false. Third, notation, numbering and computational certification fall short of publication standard. **Recommendation: major revision.** The revised paper accompanying this report resolves every item below; its status table and Appendix C give the details.

# Part I. Hidden logical vulnerabilities

- **I.1** **\[critical\]** *The field in Conjecture 6.10, and the truth of its complex form.*

  *Location.* Definition in §6.5 (weights in $\mathbb{Q}(z)$); addendum Theorem 10.8 (“weights anywhere in the algebraic closure”); addendum open problem 2 (“A counterexample would have to realize a complex frozen environment by genuine vertices…”); Proposition 10.10 and Theorem 11.5. The conjecture is stated for $\mathbb{Q}(z)$, but the addendum treats $\overline{\mathbb{Q}(z)}$ as the natural setting and reads the absence of small counterexamples as evidence. Over $\overline{\mathbb{Q}(z)}$ the statement is false: the directed $m$-cycle with $m+1$ isolated vertices has a potential whose cycle weights are $0$ and whose environment weights are the roots of an explicit polynomial $Y_m$; the smallest instance has seven vertices.

  *Resolution.* Theorem 6.10 now names its fields ($\mathbb{Q}(z)$, $\mathbb{R}(z)$, $\mathcal{P}_{\mathbb{R}}$) and is proved (Theorem 9.5); Lemma 10.1 shows that “some field containing $\mathbb{Q}(z)$” and $\overline{\mathbb{Q}(z)}$ are the same question; Theorem 11.6 and Corollary 11.7 disprove the complex form; Proposition 11.9 computes the arithmetic of the weights; the remaining complex question is Open Problem 11.19.

- **I.2** **\[major\]** *Nonvanishing is missing from the definition of a potential.*

  *Location.* Definition (local weighted matching potential), §6.5. The definition requires $R(W,v)=\mu_a(U[W-v])/\mu_a(U[W])$ but not $\mu_a(U[W])\neq0$. The addendum (“the defining identity says $M(W)\neq0$”) and every Gröbner computation (saturation by all $\mu_a(U[W])$) assume it, so the treatise and the addendum work with different definitions.

  *Resolution.* The definition in §6.5 requires $\mu_a(U[W])\neq0$ at every reachable position; Lemma 10.1 encodes it by the Rabinowitsch variable.

- **I.3** **\[major\]** *The band statement of Theorem 11.4 is proved for one graph $U$ only.*

  *Location.* Theorem 11.4 (“for real $z$ outside $\{0,\pm1,\pm\sqrt2\}$, a real frozen potential exists if and only if $|z|<2$”) and the sentence “The band edge is sharp” after Theorem 11.18. The proof solves the frozen identities for the six edges between the triangle and two frozen vertices. Nothing is shown for other environments, although the text uses the statement as if it did.

  *Resolution.* Proposition 9.14 proves, by Brändén’s criterion applied to the path data alone, that at every real $z$ with $|z|>2$ no frozen potential of the directed triangle with real weights exists, for any environment and any $U$; existence inside the band is Theorem 11.3.

- **I.4** **\[major\]** *The “free values” obstacle of §11.7 is not an obstacle.*

  *Location.* §11.7, “What is missing for general environments”: “For $m\ge5$ there are $2^m-m(m-1)-2$ free values, and we have no argument that controls them.” This is the step that kept Conjecture 10.9$'$ open. Real stability controls every coefficient of the generating polynomial to first order at $z=\infty$ from its values on cyclic intervals.

  *Resolution.* Lemma 9.8; with Lemma 9.9 and Lemma 9.10 it gives Theorem 9.5, hence Theorem 6.10 and Theorem 9.12 (formerly Conjecture 10.9$'$) for real weights. Remark 9.11 records the correction.

- **I.5** **\[major\]** *Overlapping path kinds in Lemma 11.9.*

  *Location.* Lemma 11.9 (symmetric branches act as exits), proof of the converse. The kinds overlap: a path $\beta\cdot h$ with $\alpha$ empty is of kinds 1 and 2, and the single vertex $h$ of kinds 0 and 1. The concluding sentence “a vertex set met by paths of two different kinds … the kinds are 2 and 3” is false for the definitions given.

  *Resolution.* Lemma 13.3 defines four disjoint kinds by the vertex set and the start, proves the product formula for each, and treats the only mixed case (kinds 2 and 3) explicitly.

- **I.6** **\[major\]** *The $\Sigma^0_2$-hardness proof uses an ineffective enumeration.*

  *Location.* Theorem 8.3, hardness. The diagonal language is built from “an effective enumeration of polynomial-time machines”, which does not exist; one needs clocked machines, which cover every polynomial-time language.

  *Resolution.* Theorem 8.3 uses an enumeration of clocked polynomial-time machines.

- **I.7** **\[major\]** *Algebraic constants in the real-number argument.*

  *Location.* Theorem 4.5. The result of Allender, Bürgisser, Kjeldgaard-Pedersen and Miltersen places the Boolean part of *constant-free* polynomial-time real computation in $\mathsf{P}^{\mathrm{PosSLP}}$. The theorem is stated for machines with algebraic constants, which needs the elimination of algebraic constants.

  *Resolution.* Theorem 4.5 cites Chapuis–Koiran for that step.

- **I.8** **\[major\]** *Overstated hardness input in Theorem 6.18(1).*

  *Location.* Theorem 6.18(1): “Distinguishing Grundy values 1 and 2 of UVG positions is PSPACE-hard (Burke–Ferland–Teng 2021).” This specific form could not be matched to the cited statement. The proof needs only that Grundy values are hard to compute.

  *Resolution.* Theorem 6.18: the Grundy value is recovered as the unique $m\le n$ with $\mathsf D(G)+\mathsf D(*m)\sim0$, and $\mathsf{PSPACE}$-hardness of computing Grundy values is cited.

- **I.9** **\[major\]** *Informal treewidth dichotomy.*

  *Location.* Addendum, end of “Potentials versus treewidth dynamic programming”: “Any algebraic object computed by a width-bounded dynamic program for DVG faces the same choice: it carries valuations, or it carries exponentially large exact data.” Nothing in the text covers objects indexed by other states, or represented other than by explicit rational functions.

  *Resolution.* §14: Theorem 14.2 (degree $\ge2^{m}-1$ for every exact multiplicative potential, however indexed), Lemma 14.4 (the recursion closes on canonical subgames), Proposition 14.5 ($2^{m}$ distinct resolvents at one vertex at treewidth $3$), Theorem 14.6 (linear-size circuits when components are small), Open Problem 15.9.

- **I.10** **\[minor\]** *Missing hypothesis check in Lemma 11.14(2).* The lemma applies Theorem 11.1 to $\Gamma(K'+u_iH_i;A')$ without checking that the arcs of $A'$ join non-adjacent vertices of the enlarged graph. They do, because the arcs avoid $u_i$ and every added edge contains it.

  *Resolution.* Lemma 12.8 states and proves this, for $p$ tails.

- **I.11** **\[minor\]** *The exclusion of graphs inside the symmetric part in Theorem 10.8.* It is attributed to “Theorems 10.1 and 10.2 … over any field”. Theorem 10.1 needs normal weights and does not apply; Theorem 10.2 alone suffices, applied to the tail of a one-way arc.

  *Resolution.* Proof of Theorem 10.15.

- **I.12** **\[minor\]** *The constant-self-energy remark after Theorem 11.18.* It asserts the extension to $\eta\equiv e$ by “shifting $z$” and then claims Conjecture 10.9$'$ on odd cycles for this class, which also needs Theorem 11.10.

  *Resolution.* Folded into Theorem 11.18 with the substitution $y=z-e$; the conjecture itself is now Theorem 9.12.

- **I.13** **\[minor\]** *Sign-coherence of the potential of Theorem 6.6 is asserted without proof* (Theorem 10.5). It holds because the weights $z-\eta$ have degree $1$, so the potential is normal.

  *Resolution.* Corollary 10.12.

- **I.14** **\[minor\]** *Non-rigorous test points.* Proposition 6.9(3) rests on Gröbner bases “at two test points”; a solution with a pole or a vanishing denominator at a test point is lost by specialization.

  *Resolution.* Theorem 10.15 works over $\mathbb{Q}(z)$ with saturation; the two-point filter is used in the revision only as a pre-filter followed by the exact test (Proposition 11.10).

- **I.15** **\[minor\]** *Wrong labels and names.* The results table of the addendum cites “(Conj 10.15)” for Conjecture 10.9$'$; 10.7$''$ is a Proposition in the table and a Theorem in §11.2; its proof is described as “strongly Rayleigh at three real points”, but it uses Brändén’s Rayleigh-difference criterion for real-stable polynomials, while “strongly Rayleigh” means real stable with nonnegative coefficients.

  *Resolution.* Corollary 10.14; Appendix C.

- **I.16** **\[minor\]** *Bookkeeping.* The verification table of §11.4 lists “Predictions of §11.3 on 5 to 10 vertices: 20 digraphs”; `predictions.py` confirms 17, and the other 3 are in `families_check.py`.

  *Resolution.* Appendix A lists each claim with its script.

- **I.17** **\[minor\]** *Vadhan.* “We did not have Vadhan’s full text, so for maximum matchings we use only the planar bipartite case” is inappropriate in a paper and leaves the bounded-degree case unclear.

  *Resolution.* §6.6 states Vadhan’s theorem (matchings and extremal variants, planar bipartite, bounded degree) and the Miklós–Krész withdrawal.

- **I.18** **\[minor\]** *Quantum attributions* (§3.1 table, §6.4). The $O(\sqrt N)$ bound for balanced formulas is Ambainis–Childs–Reichardt–Špalek–Zhang; $N^{1/2+o(1)}$ for arbitrary formulas is Childs–Cleve–Jordan–Yonge-Mallo and ACRŠZ; Reichardt proves tightness of the general adversary bound.

  *Resolution.* Theorem 3.1 and the quantum paragraph of §6.

- **I.19** **\[minor\]** *Zhuk–Martin and Zhuk.* “Some EGP languages are in P” and “a $\Pi_2^{\mathsf{P}}$-complete language on six elements” are stated without precise pointers. The robust statement is that the three-element classification exhibits EGP languages that are not $\mathsf{PSPACE}$-complete.

  *Resolution.* §7.2.

- **I.20** **\[minor\]** *Alternation and ranking formalities.* Alternation was imposed at all positions, terminals included, where the operators are irrelevant; and in Theorem 2.1(3) the rank $\mathrm{rk}_{+}$, defined only on $\Phi^{-1}(+1)$, is applied to every option of a Min position without noting that these options lie in $\Phi^{-1}(+1)$ by condition (2).

  *Resolution.* §1.1 and Theorem 2.1.

- **I.21** **\[minor\]** *Space bound.* Membership of connected components of succinct state graphs was proved via Savitch ($\mathsf{DSPACE}(N^2)$); Reingold gives $\mathsf{DSPACE}(N)$.

  *Resolution.* Proposition 3.6.

- **I.22** **\[minor\]** *Corollary 6.8 is stated with one maximum-matching computation per vertex*; one Gallai–Edmonds decomposition per component suffices.

  *Resolution.* Corollary 6.8.

# Part II. Missed research opportunities

- **II.1** *A first-order Brändén argument.* The addendum’s real-weight results (Theorems 10.3, 10.5, 10.7$'$, 10.7$''$) all extract information from Rayleigh differences at a few points, one cycle length at a time. Applied to the whole generating polynomial to first order at $z=\infty$, the same criterion settles every digraph.

  *Resolution.* Theorem 9.5, with Lemma 9.7, Lemma 9.8, Lemma 9.9 and Lemma 9.10.

- **II.2** *Generality of the real obstruction.* Nothing in the argument uses matchings beyond the Heilmann–Lieb theorem. The obstruction therefore covers edge-weighted matching potentials, frozen environments and Hermitian determinantal potentials (Borcea–Brändén), turning the “reservoir engineering” analogy of the addendum into a theorem about lossless systems.

  *Resolution.* Definition 9.3, Example 9.4, Corollary 9.6, Remark 9.15(2).

- **II.3** *The sharp band for all environments.* See I.3; Proposition 9.14.

- **II.4** *Pruned cores (§11.3): a continuant criterion for all blow-ups.* The addendum classifies one blown-up vertex and balanced blow-ups separately and attributes the set potentials of $C_4(1,2,1,2)$, $C_4(1,3,1,3)$, $C_4(2,3,2,3)$ and $C_6(1,2,1,2,1,2)$ to a “hidden symmetry”. The game tree from a start in a blow-up is spherically symmetric, so every product $\Gamma(P)$ is a ratio of backward continuants that depends only on the start class and the length.

  *Resolution.* Theorem 13.9 reduces set potentials on every blow-up to finitely many continuant identities; the twin identity Lemma 13.10 is the hidden symmetry; Corollary 13.11 proves the “if” half of Conjecture 11.27(1) for all two-periodic blow-ups. The “only if” half is checked on $331$ size vectors (Conjecture 13.16).

- **II.5** *The two-tail phase (§11.6): $p$ tails and matching-type terminals.* Lemma 11.14 is the case $p=2$ of a general reduction: a play crosses at most $p$ times, and a $p$-tail game is undirected geography with $p$ terminals whose values are $(p-1)$-tail games.

  *Resolution.* Lemma 12.8; Open Problem 15.5; the barrier examples are now regression tests (Proposition 12.6 and Proposition 12.9).

- **II.6** *Treewidth: precise barriers and a construction.* See I.9. Besides the barriers, component-state circuits give Corollary 14.7: DVG is decidable in time $2^{k}\operatorname{poly}(n)$ with $k$ the size of the largest non-symmetric strong component.

- **II.7** *The complex question has a finite answer.* The addendum’s frozen rescues (Theorems 11.4–11.5) and its open problem 2 lead directly to a genuine counterexample once the cycle weights are set to $0$: every identity then prescribes one elementary symmetric function of the environment weights (Vieta), and the starts at the environment vertices collapse to a single condition.

  *Resolution.* Theorem 11.6, Proposition 11.8, Proposition 11.9 and Proposition 11.10, Remark 11.11.

# Part III. Formalization deficits

- **III.1** *Notation collisions.* Symbols carry several meanings, often within one section:

  <div class="center">

  | Symbol    | Meanings in the submission                                                                         | Revision                                                  |
  |:----------|:---------------------------------------------------------------------------------------------------|:----------------------------------------------------------|
  | $V$       | value; vertex set                                                                                  | $\mathcal V$ for the value                                |
  | $U$       | payoff; the graph of a potential                                                                   | $\mathcal U$ for the payoff                               |
  | $M$       | state graph; matchings; potential values $M(S)$                                                    | $\mathcal G$, $M$, $\mathcal M$                           |
  | $m$       | cycle length; frozen potential values $m(W)$; $m^*(G)$                                             | $\mathcal M_F$, $\mathfrak m(G)$                          |
  | $\Gamma$  | constraint language; path products; the game $\Gamma(H;A)$; graphs $\Gamma_n$                      | $\Gamma(P)$ only for products; $\mathfrak G(H;\mathsf A)$ |
  | $K$       | current graph; continuant; field                                                                   | $K$; $\mathcal K$; $\mathbb K$                            |
  | $D$       | digraph; divisor map; a machine; $D_i$ in a proof                                                  | $\mathsf D$ for divisors                                  |
  | $\rho$    | fresh root resolvent; Green’s function                                                             | $R^{\downarrow}$; $\rho(S,s)$                             |
  | $\lambda$ | first-order coefficient $\lambda(P)$; $\lambda=m(\emptyset)$                                       | $\ell(P)$; $\lambda$                                      |
  | $e$       | self-energy constant; $z^{-1}$-coefficients; elementary symmetric functions; truncated exponential | $e$, $e_v$, $e_j(a_F)$, $\exp_m$                          |
  | $\kappa$  | first-order combination; canonical subgame (new)                                                   | $\kappa(S)$; $\operatorname{cs}(W,v)$                     |

  </div>

  *Resolution.* Notation table Appendix B. Local uses of $\tau,\sigma,T,C,X$ inside single proofs are kept and declared there.

- **III.2** *Numbering.* Primes and letters (10.3$'$, 10.7$'$, 10.7$''$, 10.9$'$, 10.11$'$, 11.1$'$, Lemmas 10.A–C), a gap (no 11.7), and theorem-level items that are corollaries of later results.

  *Resolution.* Consecutive numbering per section in both classes; Appendix C maps every old number.

- **III.3** *Undefined or informal notions.* “Real weights” (is $\mathcal{P}_{\mathbb{R}}$ allowed?), “Puiseux series” (which field, which convergence?), “admissible self-energies”, “pruned core”, “occurs” in Conjecture 11.27(2), and “degree” of an algebraic function (depends on the embedding).

  *Resolution.* §9.1, §10.1 (degree for each embedding), Definition 11.1, Lemma 13.3 (pruned core), Conjecture 13.16(2).

- **III.4** *Bibliography.* The addendum cites Brändén’s criterion, Heilmann–Lieb, Metelmann–Clerk and Pólya–Schur informally or not at all; the treatise cites Chapuis–Koiran and Reingold nowhere although it needs them.

  *Resolution.* A complete bibliography with all keys used.

- **III.5** *Certification of computations.* Evaluation at $z=37/11$ is justified as “not an algebraic integer, so never an eigenvalue”; that argument certifies refutations of identities only, not the absence of solutions of a polynomial system. Several rows of the verification tables report counts without the script that produces them.

  *Resolution.* Appendix A separates proofs from computer-assisted results and names a script for every row; the archive adds tests that compare every new formula with brute force.

- **III.6** *Missing response document.* The treatise’s acknowledgment refers to a “Response to review” that is not included. This report, with the concordance, serves as the response for this round.

- **III.7** *Status labeling.* The addendum’s results table lists Theorem 10.9 (shape of a counterexample) as established “by proof” for real weights, although for real weights it is vacuous once Theorem 6.10 holds.

  *Resolution.* Proposition 10.16 restates it over arbitrary fields, where it has content.

# Part IV. Outcome of the mandated research problems

<div class="center">

| Problem                                                                           | Outcome                                                                                                                                                                                                                                                                                                                                                                                      |
|:----------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Conjecture 6.10 over $\mathbb{Q}(z)$, $\mathbb{R}(z)$, $\mathcal{P}_{\mathbb{R}}$ | **Proved**: Theorem 6.10 and Theorem 9.5.                                                                                                                                                                                                                                                                                                                                                    |
| Conjecture 10.9$'$ (frozen environments), real weights                            | **Proved** for every strongly connected non-symmetric digraph and every real self-energy vanishing at infinity: Theorem 9.12.                                                                                                                                                                                                                                                                |
| Conjecture 6.10 over $\overline{\mathbb{Q}(z)}$ (complex weights)                 | **False**: Theorem 11.6 and Corollary 11.7; the weights have monodromy $S_{m+1}$ for $3\le m\le6$ and are unramified at $\infty$ (Proposition 11.9), so no dihedral-monodromy argument can exclude genuine realizations (Remark 11.11). Open over $\overline{\mathbb{Q}(z)}$: strongly connected digraphs (Open Problem 11.19).                                                              |
| Treewidth bridge (“path-dependent bag potential”)                                 | The resolvent recursion closes on canonical subgames (Lemma 14.4); explicit rational potentials have degree $\ge2^{m}-1$ however indexed (Theorem 14.2); exact resolvents are not bag-local even at treewidth $3$ (Proposition 14.5); circuits over component states have linear size for small components (Theorem 14.6). Polynomial-size circuits at bounded treewidth: Open Problem 15.9. |
| Star-tail rule                                                                    | Theorem 12.2 (implemented as `star_tail.star_rule`, tested against exhaustive search).                                                                                                                                                                                                                                                                                                       |
| Cacti taxonomy                                                                    | Theorem 13.13 and Corollary 13.14; blow-ups: Theorem 13.9 and Corollary 13.11.                                                                                                                                                                                                                                                                                                               |

</div>
