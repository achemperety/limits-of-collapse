"""dvgcollapse: exact tools for "Alternation, Invariants, and the Limits of Collapse" (revised version).

Modules
  dvg                    digraphs, DVG positions, outcomes, exact root resolvents over Q(z)
  matching               matching numbers, Gallai-Edmonds classes, weighted matching polynomials
  potentials             Theorem 6.6 potentials; set and valuation potentials (Section 13)
  real_stability         Rayleigh differences, the first-order data of Theorem 9.5, certificates
  groebner_saturation    saturated Groebner tests over Q(z) for (frozen) matching potentials
  complex_counterexample the counterexample D_m to Theorem 6.10 over the algebraic closure (Section 11)
  star_tail              the common-tail evaluator (Theorem 12.2, Corollaries 12.3-12.4)
  blowups                blow-ups of directed cycles: continuant criterion and twin identity (Section 13)
  treewidth              ladders, intrinsic degree, canonical subgames, the memory comb, component-state
                         circuits (Section 14)
  solver                 a unified exact DVG solver dispatching to the methods above (also a CLI)
"""
__version__ = '1.1.0'
