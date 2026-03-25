# Questions for the Reviewer

## Critical checks

1. **Reproduce all numbers.** Run CIR_pipeline.py and compare every number in
   VERIFICATION_CHECKLIST.md against your output. Any discrepancy = bug.

2. **Check the chemical index formulas.** Are DBE, AI_mod, and NOSC computed
   correctly? Compare against Koch & Dittmar (2006) and LaRowe & Van Cappellen (2011).

3. **Check compound class boundaries.** Do they match Kim et al. (2003)?
   Note: the priority order matters — AI_mod thresholds are checked before H/C and O/C.

4. **Verify the worked example.** Manually compute CIR for C15H20O7 using the
   formulas in CIR_METHODS_DETAIL.md. You should get log_CIR ≈ 1.172.

## Scientific questions

5. **Is the CIR scoring function well-motivated?**
   The five components are heuristic. Could the weights (0.8, 0.7, 0.5, 0.12, 0.05)
   be derived more rigorously? Would a sensitivity analysis change the conclusions?

6. **Is the TIMS calibration valid?**
   The linear rescaling to log10(6)–log10(40) assumes the score distribution maps
   uniformly to the TIMS range. Is there a better calibration approach?

7. **Could the CIR-richness anti-correlation be artifactual?**
   If samples with more formulae tend to include more low-intensity (noisy) peaks,
   and those peaks happen to be from classes with lower CIR, the anti-correlation
   could be spurious. Check: does the correlation survive after filtering to only
   formulae detected in all three replicates of a site?

8. **Is the CRAM-latitude correlation robust?**
   r = -0.74 is very strong for environmental data. Check:
   - Are there high-leverage outlier sites driving the correlation?
   - Does it survive leave-one-out cross-validation?
   - Is it confounded with another variable (e.g., DOC, which also varies with latitude)?

9. **Are the constraint rules biased toward certain classes?**
   CRAM has 1 required + 4 forbidden = 5 constraints.
   Lipid-like has 0-1 required + 4 forbidden = 4-5 constraints.
   Carbohydrate-like has 1-2 required + 5 forbidden = 6-7 constraints.
   Does the number of constraints drive CIR differences, or does the underlying
   structural chemistry dominate? Test: run CIR with uniform constraints across
   all classes and see if the class ordering persists.

10. **What is the actual scientific value?**
    Be honest: CIR is a heuristic proxy. The CRAM-latitude correlation is the
    strongest empirical result. But does CIR tell us anything we couldn't learn
    from simpler metrics (e.g., mean DBE, mean AI_mod, or mean molecular weight
    of the CRAM fraction)? Compute those and compare correlations with latitude.
