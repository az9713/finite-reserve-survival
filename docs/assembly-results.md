## Certified numerical conclusions

For the two-pair, four-fuel benchmark with formation deadline one and required post-assembly lifetime four,

\[
0.07486803450\le V^*\le0.07496800805.
\]

The decimal endpoints above are rounded outward. The full rational endpoints are retained in the certificate. The absolute optimization gap is below \(10^{-4}\), or 0.01 percentage points. The design

\[
\alpha=529/1024,\qquad z=63/256,\qquad
\theta=\log(256/63)
\]

attains the certified lower bound to its probability-enclosure precision. Every boundary point of the compact square has

\[
V\le0.07296461445<0.07486803450.
\]

Thus **every global maximizer is strictly interior**. It uses positive assembly and maintenance capacity and a finite, positive binding strength in this declared family. The result holds for every real parameter pair in the square, including the limiting infinite-binding boundary. It certifies existence of an interior optimum without asserting a unique maximizer or an exact symbolic value for it.

The cover has 46,655 tree nodes, 23,328 leaf rectangles, and 23,969 evaluated rational points. The replay script recomputes the point intervals, verifies every partition and interpolation bound, and compares all points with an independent matrix-exponential evaluation. Its output also supplies an exact positive-vector enclosure for the reachable spectral exponent at the certified design and the conditional resource allocation.

The first run rebuilt parameter-independent curvature coefficients for every rectangle and took about 468 seconds. Caching those exact coefficients reduced the repeat optimization to about 141 seconds and reproduced the same rational bounds and construction. The optimization is over the specified two-parameter static family. Arbitrary reaction mechanisms and time-varying controllers require different comparisons.

[Full rational certificate](../results/assembly.json) · [Complete certificate replay receipt](../results/assembly-certificate-check.json) · [Independent stochastic evidence](../results/assembly-evidence.json)
