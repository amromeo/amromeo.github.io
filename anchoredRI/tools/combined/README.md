# Task G package regression

The regression fits all 749 prespecified audit datasets through the public
`fit_anchored_ri()` entry point. It uses the observed routine vector and
guidance panel from the Stage 3 export, with the recorded rounding unit,
Stage 2 fit seed and Stage 1 PIT seed. It checks the Arm B decision, reason
and model choice, the Stage 1 decision and reason when the fallback runs,
the reporting method and the released limits. Endpoints must match within
an absolute tolerance of 1e-10. A nonreturn must keep both endpoints missing.
[Task G decision](https://github.com/amromeo/amromeo.github.io/blob/bd0af4d5eb7b0c7c6657b9389a26bf32377088d3/anchoredRI/tools/combined/DECISION.md)

The Actions workflow builds the package once in the container from the
audited Task G run. It uses Stage 3 dependencies from the frozen library bundle. The package check includes the external refineR comparator. The audit runs with the external refineR package absent to exercise the vendored fallback. Each of the
12 jobs installs that same compiled library and pins the BLAS settings.
`prepare-audit.py` verifies the data against the published checksums.
`audit-package.R` compares the public function's results with the recorded
routes and endpoints. The final job checks the 749 unique IDs and writes
the route counts and largest absolute endpoint difference.
[Audited Task G run](https://github.com/amromeo/amromeo.github.io/actions/runs/45614043039)
