# Anti-invariant weave verification record

2026-09-07. All twenty fixed-density bounds and all ten density-chord intervals PASS. The completed directed run checked 100466 boxes, accepted 50243 leaves, and took 157.13 seconds summed over its twenty points. The largest certified chord cap is below `0.498795641`, hence strictly below `499/1000`.

The canonical full rational results are `computations/results/principle_synthesis_2026_09_07_anti_weave_interval.json`. Every scalar rectangle was outward-rounded at 40 decimal digits. The independent source audit and independently implemented entropy-tangent chord check are in `principle_construct_2026_09_07_skew_lift_and_paired_weave_audit.md` and its linked computations.

## Preserved exploratory and operational history

- The first floating grid used `p=.96,t=4.85`; it did NOT certify a subhalf theorem. Its largest grid cap was approximately `.50147717452` at `theta=0`, while the midpoint was approximately `.49357512388`. This unsuccessful parameter choice is reproducible with the default diagnostic script.
- The next floating grid used `p=.96,t=3.5`, and its largest grid cap was approximately `.49864889282`. These scans selected rational upper targets; neither scan was used to accept interval rectangles.
- A preliminary directed run at indices `1,10,20` passed (4465, 3041, 26097 boxes). The full twenty-point run then repeated these points and checked all the remaining points.
- A polling wrapper attempted to store an undefined completed-session ID and reported a serialization error AFTER it had stored the returned output. Inspection recovered every complete `POINT_RESULT` line, indices1 through20. The tool truncated part of the redundant aggregate `FINAL_RESULT`; all point records remained intact. The cheap chord computation was rerun and parsed completely. The canonical JSON is reconstructed from those verified complete records, not from the truncated duplicate line.

No solver or interval budget stop occurred in the successful parameter certificate. The family result does not improve the original all-signing upper endpoint and does not establish a seed-sensitive recurrence.

