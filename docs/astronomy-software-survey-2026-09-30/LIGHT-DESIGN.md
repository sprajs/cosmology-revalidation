# Cosmology design relocated

The expanded design now lives in the Cosmology project:

- [Main plan](/home/szymon/Projects/cosmology/Plan/README.md)
- [Sequential build order](/home/szymon/Projects/cosmology/Plan/BUILD-ORDER.md)
- [Agent development and analysis instructions](/home/szymon/Projects/cosmology/AGENTS.md)
- [Historical light design](/home/szymon/Projects/cosmology/Plan/LIGHT-DESIGN.md)

The current design uses Rust for agent/data orchestration and C++20 for owned science and numerical work, with measured CPU SIMD/assembly, optional CUDA and later cluster backends. The earlier light design is preserved as historical reference; active requirements are in the new plan.

The software survey remains in this directory and is also snapshotted under the new plan's reference folder. This relocation did not move scientific datasets or rerun research.
