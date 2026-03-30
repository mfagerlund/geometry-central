# GFR vs FlipOut Benchmark Results
Generated: 2026-03-30 14:48
Seed: 42

## Executive Summary

| Mesh | Vertices | Paths | GFR Time | FlipOut Time | Speedup | vs FlipOut |
|------|----------|-------|----------|--------------|---------|------------|
| woody | 694 | 1,000 | 201ms | 874ms | **4.4x** | **0.00%** |
| pig | 1,843 | 1,000 | 890ms | 2.3s | **2.6x** | **-0.16%** |
| bunny-small | 2,503 | 1,000 | 1.5s | 3.2s | **2.2x** | **-0.05%** |
| spot | 2,930 | 1,000 | 1.4s | 3.6s | **2.5x** | **-0.20%** |
| alligator | 3,208 | 1,000 | 2.7s | 4.2s | **1.6x** | **0.00%** |
| homer | 6,002 | 1,000 | 3.6s | 7.9s | **2.2x** | **-0.14%** |
| fandisk | 6,475 | 1,000 | 3.9s | 8.0s | **2.1x** | **-0.26%** |
| cheburashka | 6,669 | 1,000 | 3.7s | 9.3s | **2.5x** | **-0.08%** |
| rocker-arm | 10,044 | 1,000 | 5.7s | 13.2s | **2.3x** | **-0.03%** |
| stanford-bunny | 34,834 | 1,000 | 24.2s | 52.3s | **2.2x** | **-0.18%** |
| horse | 48,485 | 1,000 | 28.1s | 75.5s | **2.7x** | **-0.03%** |
| happy-buddha | 49,251 | 1,000 | 44.9s | 76.3s | **1.7x** | **-0.38%** |
| armadillo | 49,990 | 1,000 | 32.4s | 89.9s | **2.8x** | **-0.11%** |
| max-planck | 50,077 | 1,000 | 31.6s | 78.8s | **2.5x** | **-0.02%** |
| bimba | 112,455 | 1,000 | 101.4s | 169.1s | **1.7x** | **-0.08%** |

**Average across 15 meshes (15,000 paths):** GFR is 2.38x faster and produces 0.12% shorter paths.

## Detailed Statistics

### woody (694 vertices, 1,000 paths)
- GFR wins: 0 (0.0%)
- FlipOut wins: 2 (0.2%)
- Ties: 998 (99.8%)
- GFR avg iterations: 3.8
- FlipOut avg iterations: 16.8

### pig (1,843 vertices, 1,000 paths)
- GFR wins: 132 (13.2%)
- FlipOut wins: 96 (9.6%)
- Ties: 772 (77.2%)
- GFR avg iterations: 4.9
- FlipOut avg iterations: 22.0

### bunny-small (2,503 vertices, 1,000 paths)
- GFR wins: 248 (24.8%)
- FlipOut wins: 270 (27.0%)
- Ties: 482 (48.2%)
- GFR avg iterations: 7.3
- FlipOut avg iterations: 22.9

### spot (2,930 vertices, 1,000 paths)
- GFR wins: 154 (15.4%)
- FlipOut wins: 124 (12.4%)
- Ties: 722 (72.2%)
- GFR avg iterations: 7.6
- FlipOut avg iterations: 33.8

### alligator (3,208 vertices, 1,000 paths)
- GFR wins: 0 (0.0%)
- FlipOut wins: 50 (5.0%)
- Ties: 950 (95.0%)
- GFR avg iterations: 25.9
- FlipOut avg iterations: 116.8

### homer (6,002 vertices, 1,000 paths)
- GFR wins: 218 (21.8%)
- FlipOut wins: 147 (14.7%)
- Ties: 635 (63.5%)
- GFR avg iterations: 11.0
- FlipOut avg iterations: 42.4

### fandisk (6,475 vertices, 1,000 paths)
- GFR wins: 107 (10.7%)
- FlipOut wins: 35 (3.5%)
- Ties: 858 (85.8%)
- GFR avg iterations: 12.7
- FlipOut avg iterations: 64.8

### cheburashka (6,669 vertices, 1,000 paths)
- GFR wins: 240 (24.0%)
- FlipOut wins: 214 (21.4%)
- Ties: 546 (54.6%)
- GFR avg iterations: 10.7
- FlipOut avg iterations: 46.1

### rocker-arm (10,044 vertices, 1,000 paths)
- GFR wins: 226 (22.6%)
- FlipOut wins: 250 (25.0%)
- Ties: 524 (52.4%)
- GFR avg iterations: 17.8
- FlipOut avg iterations: 57.6

### stanford-bunny (34,834 vertices, 1,000 paths)
- GFR wins: 434 (43.4%)
- FlipOut wins: 515 (51.5%)
- Ties: 51 (5.1%)
- GFR avg iterations: 25.6
- FlipOut avg iterations: 224.9

### horse (48,485 vertices, 1,000 paths)
- GFR wins: 213 (21.3%)
- FlipOut wins: 774 (77.4%)
- Ties: 13 (1.3%)
- GFR avg iterations: 17.2
- FlipOut avg iterations: 318.5

### happy-buddha (49,251 vertices, 1,000 paths)
- GFR wins: 620 (62.0%)
- FlipOut wins: 367 (36.7%)
- Ties: 13 (1.3%)
- GFR avg iterations: 20.6
- FlipOut avg iterations: 122.5

### armadillo (49,990 vertices, 1,000 paths)
- GFR wins: 665 (66.5%)
- FlipOut wins: 259 (25.9%)
- Ties: 76 (7.6%)
- GFR avg iterations: 32.6
- FlipOut avg iterations: 136.2

### max-planck (50,077 vertices, 1,000 paths)
- GFR wins: 155 (15.5%)
- FlipOut wins: 568 (56.8%)
- Ties: 277 (27.7%)
- GFR avg iterations: 61.4
- FlipOut avg iterations: 248.5

### bimba (112,455 vertices, 1,000 paths)
- GFR wins: 383 (38.3%)
- FlipOut wins: 595 (59.5%)
- Ties: 22 (2.2%)
- GFR avg iterations: 47.4
- FlipOut avg iterations: 280.4
