# GFR Benchmark Suite

This directory contains benchmarks that provide empirical evidence for the GFR paper.
All benchmarks are designed to be reproducible with fixed random seeds.

## Building

```bash
cd benchmark/build
cmake ..
cmake --build . --config Release
```

## Paper-Critical Benchmarks

These benchmarks produce data cited in the paper:

### multi_mesh_benchmark
**Paper section:** Results, Executive Summary, Comparison with FlipOut

Compares GFR vs FlipOut across multiple meshes, producing paper-ready markdown tables.
This is the primary benchmark for speedup and path quality claims.

```bash
./bin/Release/multi-mesh-benchmark --paths 1000
```

**Outputs:**
- Per-mesh timing comparison (GFR vs FlipOut)
- Path quality comparison (% difference in total path length)
- Iteration statistics (avg iterations, iters/face ratio)
- Cache hit rates

---

### initial_distance_test
**Paper section:** Sleeve Property Verification

Verifies that GFR's initial funnel distance is never longer than the discrete edge path.
This validates correctness of sleeve construction and the funnel algorithm.

```bash
./bin/Release/initial-distance-test --paths 5000
```

**Expected result:** Zero violations across all tested paths.

---

### fan_direction_test
**Paper section:** Fan Direction Selection

Validates that the angle-sum criterion for fan direction selection matches
minimum-length sleeves. Also measures the performance impact of suboptimal
fan direction (before the try-both fix was implemented).

```bash
./bin/Release/fan-direction-test --paths 1000
```

**Key finding:** Optimal fan direction improved GFR speedup from 1.57x to 1.90x.

---

### epsilon_test
**Paper section:** Note on acceptance threshold (Algorithm section)

Tests the effect of the acceptance threshold on convergence and path quality.
Compares strict improvement (>0) vs accepting worsening flips.

```bash
./bin/Release/epsilon-test --paths 200
```

**Key findings:**
- strict (>0): ~25 iterations, 0.39% improvement
- relaxed (<0): ~40 iterations, 0.41% improvement
- Accepting worsening flips adds 40-60% iterations for marginal quality gain

---

### epsilon_vs_flipout_test
**Paper section:** Note on acceptance threshold (Algorithm section)

Tests how acceptance threshold affects win/loss rate against FlipOut.
Counts per-path wins/ties/losses with different thresholds.

```bash
./bin/Release/epsilon-vs-flipout-test --paths 200
```

**Key findings:**
- stanford-bunny: strict 37% win rate → relaxed 43.5% (+13 paths flipped to wins)
- armadillo: strict 63.5% win rate → relaxed 61% (slightly worse)
- Relaxed threshold helps on some meshes, hurts on others
- Strict threshold provides consistent cross-mesh behavior

---

### cold_cache_benchmark
**Paper section:** Single-Path vs Batch Performance

Measures performance with cold cache (no pre-warming) to understand
cache buildup behavior during batch runs.

```bash
./bin/Release/cold-cache-benchmark --paths 100
```

**Key finding:** Both algorithms start cold; GFR's cache builds naturally during run.

---

## Diagnostic Tools

These tools are for debugging specific paths or behaviors:

### debug_path
Detailed comparison of GFR vs FlipOut on a single vertex pair.

```bash
./bin/Release/debug-path mesh.obj 123 456
```

---

### diagnose_facestrip
Compares face strips from VeryDiscreteGeodesic vs Dijkstra for debugging
sleeve construction issues.

```bash
./bin/Release/diagnose-facestrip mesh.obj 123 456
```

---

### geodesic_comparison
Single-mesh comparison of GFR vs FlipOut with detailed output.

```bash
./bin/Release/geodesic-comparison mesh.obj 100
```

---

### gfr_strategy_benchmark
Compares different GFR strategies (CornerGreedy, WedgeGreedy, CoherentMiniWedge).

```bash
./bin/Release/gfr-strategy-benchmark mesh.obj 100
```

---

## Reproducing Paper Results

To reproduce the main benchmark results from the paper:

```bash
# Main benchmark (15,000 paths across 15 meshes)
./bin/Release/multi-mesh-benchmark --paths 1000 > results.md

# Sleeve property verification (30,000 paths)
./bin/Release/initial-distance-test --paths 5000

# Fan direction validation
./bin/Release/fan-direction-test --paths 1000

# Epsilon threshold analysis
./bin/Release/epsilon-test --paths 200
./bin/Release/epsilon-vs-flipout-test --paths 200
./bin/Release/epsilon-vs-flipout-test --mesh armadillo.obj --paths 200
```

## Test Meshes

Default mesh directory: `C:/Dev/Colonel/Data/Meshes`

Meshes used in paper benchmarks:
- woody.obj (694 vertices)
- pig.obj (1,843 vertices)
- bunny-small.obj (2,503 vertices)
- spot.obj (2,930 vertices)
- alligator.obj (3,208 vertices)
- homer.obj (6,002 vertices)
- fandisk.obj (6,475 vertices)
- cheburashka.obj (6,669 vertices)
- rocker-arm.obj (10,044 vertices)
- stanford-bunny.obj (34,834 vertices)
- horse.obj (48,485 vertices)
- happy-buddha.obj (49,251 vertices)
- armadillo.obj (49,990 vertices)
- max-planck.obj (50,077 vertices)
- bimba.obj (112,455 vertices)
