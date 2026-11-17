# Spatial repair implementation record

Date: 17 November 2026. Status: COMPLETE, verified for the stated model and selected numerical cases.

Authorized task: implement the two-complex spatial maintenance model, audit originality first, compare uniform and localized starts with a calibrated well-mixed comparator, retain rigorous derivations, visuals and verification receipts.

## Gates and receipts

1. PASS: originality audit inspected primary equations before authorizing model implementation as replication/targeted extension, without a novelty claim. See spatial-repair-originality-audit.md and spatial-repair-search-receipt.json (28 queries, 20 source records).
2. PASS: periodic relative Brownian position, reversible damage and contact repair, finite fuel/waste, product equilibrium and exact event ledger. Full equilibrium and first-collapse survival are explicitly separate.
3. PASS: implemented model plus nine independent new tests; full suite 51 tests passed. A separately assembled labelled periodic-grid generator, a dense exponential and a spherical Doi reference exercise independent reductions and limits.
4. PASS for declared numerical acceptance criteria: grids 32/48/64/96, 1601 observation times, every forward fuel-edge mean calibrated before survival comparison. Finest adjacent curve change at most 3.30115e-5; no rigorous continuum discretization enclosure is claimed. Exact-arithmetic polynomial truncation bounds and calibration solve residuals are stored separately.
5. PASS: 200000 independent continuous-space Brownian paths, explicit event ledgers, seeds and sampled curves. Familywise 99% DKW sampling envelope is separate from mesh error. See results/spatial-repair-verification.json.
6. PASS: model report with proofs, three PNG/PDF figure pairs, data and source hashes. Desktop and mobile browser checks render 104 equations with no math/image/link errors or horizontal page overflow. See results/spatial-repair-ui-validation.json.

## Deliverables

- src/purposeful/spatial_repair.py: model, exact finite-grid symmetry reduction, calibration, killed semigroup.
- scripts/run_spatial_repair.py: deterministic experiment and two figure pairs.
- tests/test_spatial_repair.py and scripts/verify_spatial_repair.py: independent verification and third figure pair.
- docs/SPATIAL-REPAIR.html and docs/spatial-repair.md: scientific report and source.
- research/spatial-repair-*.md/json: audit and verification rationale.
- results/spatial-repair*.json: numerical, statistical, UI and integrity receipts.

Parallel ownership was separated between audit, model, independent verification, and integration/documentation. Independent final review found no unresolved critical, high or medium scientific issue. Project-private MEMORY.md was retained locally and is outside the delivery manifest. Earlier workspace changes were preserved.

## Scientific limits

This is maintenance of two preassembled point complexes, with fast-mixed fuel and assumed dimensionless parameters. It is not spatial self-assembly or experimentally calibrated biology. Localized survival effects persist against both calibration choices, but no new general concept or theorem is established. A chemostatted/inventory ablation would be required to attribute the gap specifically to fuel depletion. Statistical bounds compare the true Brownian survival law to fixed finest-grid-calibrated count coefficients; exact continuum calibration error remains unenclosed.
