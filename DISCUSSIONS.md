# GitHub Discussions Guide

GitHub Discussions is the forum for questions, proposed tests, reproductions, failed reproductions, and interpretation of this repository's atomic-spectra analysis.

## Recommended categories

- **Announcements**: Maintainer updates, releases, frozen results, and changes to canonical status.
- **General**: Broad conversation that does not yet require a structured technical report.
- **Q&A**: Installation, implementation, statistical, data-provenance, and interpretation questions. Mark a response as the accepted answer when resolved.
- **Ideas**: Proposed methods, controls, datasets, or reporting improvements with explicit success and failure criteria.
- **Show and tell**: Exact reproductions, failed reproductions, independent implementations, sensitivity analyses, and exploratory extensions.
- **Polls**: Use sparingly for community priorities, never as evidence for a scientific conclusion.

## Evidence language

Use the narrowest accurate status:

- **Reproduced**: The declared computation was rerun and agreed within a stated tolerance.
- **Independent implementation**: Different code reproduced the declared calculation.
- **Robustness supported**: A stated sensitivity or negative-control battery was passed.
- **Failed reproduction**: The declared result did not reproduce under the reported environment and procedure.
- **Incomplete**: Required data, artifacts, runtime, or provenance were unavailable.
- **Exploratory**: The analysis was not preregistered or is outside the canonical verdict path.
- **Empirically supported**: Reserved for evidence that goes beyond software reproduction and survives appropriate experimental or observational controls.

Computational agreement does not by itself establish a physical interpretation. Public NIST Atomic Spectra Database exports are the data source; this project is not endorsed, certified, or validated by NIST.

## Minimum standard for technical claims

Include the repository commit, environment, exact command, input provenance, changed parameters, random seed where relevant, machine-readable outputs, and a clear comparison against the canonical result. Report null and negative results as prominently as favorable results.

When a discussion identifies a concrete implementation task, convert or restate it as a GitHub issue. Keep open scientific interpretation, questions, and research proposals in Discussions.
