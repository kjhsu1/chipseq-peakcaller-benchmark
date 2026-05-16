# Goal

Complete the revised ChIP-seq benchmark implementation wave: replace the naive
real-study-to-PMF plan with a stronger simulator-based plan, implement the
chosen simulator workflow, prepare the EPIC2 sweep, improve repo workflow
tracking, and clean the active pipeline structure.

## Context

- Project: ChIP-seq peak-caller benchmark focused on how control/input coverage affects precision, recall, F1, and failure modes across different simulated backgrounds.
- Stack: Python, Snakemake, pandas/numpy/matplotlib, MACS2, EPIC2, Slurm, shell scripts, public ChIP-seq metadata/data, external ChIP-seq simulators such as ChIPs.
- Current state: `chipseq_pipeline_v2` is the active controlled truth-aware pipeline; `chipseq_pipeline_v2_realstudy` is the active realstudy/realism-oriented sibling pipeline; the selected clean worm broad-mark study is downloaded locally; controlled-sweep bug fixes were implemented, but old sweep plots/results need rerun/audit before being trusted.
- Working dir: `/Users/kentahsu/Code/KorfLab/Background_Forked/.codex_worktrees/feature-realstudy-benchmark-prototype`
- Constraints: do not directly treat one real pileup as the true generative PMF; keep ontology windowing; prefer a published external ChIP-seq simulator for realistic read generation; EPIC2 only needs local non-execution validation because it runs on cluster, not local MacBook; avoid helper-script bloat and unnecessary config parameters.
- Audience: Kenta and professor/lab stakeholders reviewing whether the benchmark is scientifically credible and execution-ready.
- Research front 1, controlled benchmark repair: the existing controlled simulator is valuable because it is truth-aware and isolates bias effects, but recent bug hunting found issues in the old sweep interpretation. The current work must repair aggregation consistency, treatment/control randomization, broad-truth evaluation, peakcaller path handling, and auditability before old PR/F1/control-coverage conclusions are trusted.
- Research front 2, revised realistic benchmark plan: the original real-study-conditioned plan was conceptually useful but scientifically weak if it directly converted an observed real pileup into a "true" read-generation PMF. The revised plan should keep real-study grounding and ontology evaluation, but use a stronger published ChIP-seq simulator where possible to model realistic experimental processes.
- Research front 3, simulator search and selection: before building the new realistic pipeline, define what a simulator must support for this project, search the literature/web for ChIP-seq simulators, compare candidates, and justify the best choice. Compatibility with control-depth sweeps and ontology-based failure analysis matters more than simulator popularity alone.
- Research front 4, simulator workflow execution: after choosing the simulator, build or adapt a Snakemake workflow that can produce simulated treatment/control reads, preserve truth or truth-proxy information, feed existing alignment/peakcalling where practical, and produce outputs that the ontology evaluation can consume.
- Research front 5, EPIC2 comparison sweep: after the controlled-sweep bug fixes are in place, prepare a matched EPIC2 version of the existing controlled sweep. The purpose is to test whether the control-depth conclusions are MACS2-specific or persist under a broad-domain-oriented peak caller. EPIC2 execution is expected to be cluster-only, so local work should validate configs, run tables, paths, and scripts.
- Research front 6, ontology-based evaluation: keep the region ontology idea as the interpretive core. The final benchmark should not only ask whether more control helps, but where it helps, where it hurts, and which background/foreground/confounder classes drive failures.
- Research front 7, real-study assets and manifests: the new clean worm broad-mark study has been selected and downloaded locally, and realstudy manifests exist. These assets should support the revised realistic benchmark, but should not be treated as direct biological truth without a defensible simulator/truth-proxy strategy.
- Research front 8, repo organization and file-system hygiene: the repo was recently flattened and legacy materials were archived. The next cleanup front is to keep active workflows obvious, prune low-utility scripts from active paths, add a plain-language scripts index, and avoid growing many helper scripts unless they clearly reduce complexity.
- Research front 9, goal/task operating system: create a durable `/goals` workflow inside the repo using `GOAL.md`, `PLAN.md`, `PROGRESS.md`, and `LOG.md`, and update `AGENTS.md` so future Codex runs know how to report progress, log blockers, preserve context, and avoid repeating earlier confusion.

## Success Criteria

All must be true.

Research front 1, controlled benchmark repair:

1. Controlled-sweep bug fixes are treated as prerequisite context for all new sweep interpretation.
2. The EPIC2 follow-up sweep is designed after the current bug-fix refactor, not against the old scientifically suspect evaluation path.
3. Local validation confirms the repaired controlled pipeline still passes its test suite.
4. Old six TF-clean sweep outputs are treated as baseline comparison artifacts, not trusted final evidence, until rerun or audited under the repaired logic.

Research front 2, revised realistic benchmark plan:

5. The 11-step plan is revised from the naive "real study pileup -> direct PMF conversion" approach into a scientifically stronger plan that uses real data for parameterization/reference context but avoids claiming the observed pileup is biological truth.
6. The revised 11-step plan still preserves the ontology windowing idea, including region classes/failure-mode evaluation as a central endpoint rather than a side analysis.
7. The revised plan explicitly explains the first-principles split between the controlled truth-aware benchmark and the more realistic simulator-based benchmark.
8. The revised plan explicitly states where an external ChIP-seq simulator replaces custom direct read sampling from real pileups.

Research front 3, simulator search and selection:

9. Before choosing a simulator, a clear compatibility checklist is written for this project, including support or workaround paths for control-depth sweeps, treatment/control generation, known/estimated truth regions, ontology-window evaluation, Snakemake execution, Slurm execution, reproducibility, and realistic ChIP-seq experimental noise.
10. A web/literature search is performed for candidate ChIP-seq simulators, including ChIPs and other credible alternatives.
11. A simulator comparison table is produced that lists each candidate, evidence/source, strengths, weaknesses, installation/runtime concerns, compatibility with the ontology framework, compatibility with control-depth sweeps, and final recommendation.
12. The chosen simulator is justified in writing, with a clear explanation of why it is better for this benchmark than direct real-pileup resampling.

Research front 4, simulator workflow execution:

13. A Snakemake workflow is built or adapted to run the chosen simulator path end to end at workflow level.
14. The simulator workflow has explicit inputs and outputs for real-study/reference inputs, simulated treatment reads, simulated control reads, coverage-depth parameters, seeds, metadata/manifests, and downstream evaluation handoff.
15. The simulator workflow connects to the existing alignment and peak-calling structure where practical instead of inventing a new orchestration model unnecessarily.
16. The simulator workflow is compatible with ontology-based evaluation, either by preserving simulator truth outputs or by producing a documented truth-proxy/region-label handoff.
17. Slurm assets are created for the chosen simulator workflow, including at least one `.sbatch` file and one shell submitter/launcher script.
18. The Slurm/shell assets document expected working directory, environment setup, config path, output path, resources, and how to launch the workflow on the cluster.

Research front 5, EPIC2 comparison sweep:

19. The EPIC2 sweep reuses the same parameter family as the current controlled sweep, changing the peak caller to EPIC2 while keeping the comparison as apples-to-apples as possible.
20. The EPIC2 sweep design targets the intended approximately `1700` run scale or explains the exact resolved run count if the Cartesian product differs.
21. The EPIC2 workflow/config changes are made only where necessary; no unnecessary new sweep parameters are introduced.
22. Slurm and shell scripts are created for the EPIC2 sweep, with cluster execution clearly documented.
23. Local validation confirms EPIC2 run-table generation, config parsing, path expansion, and non-EPIC2 workflow code paths, while explicitly documenting that actual EPIC2 execution is cluster-only.
24. The EPIC2 deliverable includes the workflow/config changes needed, Slurm scripts, shell scripts, expected run table, and commands needed to compute the final performance table and six PR/F1/control-coverage plots after cluster execution.
25. The final EPIC2 target deliverables are documented as: computed performance table, screenshots or generated files for six PR/F1 vs control-coverage plots, and a short note comparing the EPIC2 setup against the MACS2 setup.

Research front 6, ontology-based evaluation:

26. The revised benchmark keeps ontology-based evaluation as the interpretive endpoint, not just a visualization layer.
27. The simulator/evaluation handoff preserves enough truth or truth-proxy information to classify windows and evaluate performance by ontology class.
28. The final evaluation plan can answer where more control helps, where it hurts, and which background/foreground/confounder classes drive failures.

Research front 7, real-study assets and manifests:

29. Existing realstudy manifests and downloaded clean worm broad-mark inputs are incorporated into the revised plan as reference/parameterization assets.
30. The realstudy assets are not described as direct truth unless a defensible simulator/truth-proxy mechanism supports that claim.
31. Any needed manifest or path updates are made without breaking `chipseq_pipeline_v2_realstudy`.

Research front 8, repo organization and file-system hygiene:

32. The active `scripts/` folder is reviewed and pruned.
33. Low-utility scripts that have not been used recently are removed from the active script path or archived, not left mixed into the active workflow.
34. The remaining active script set is minimal enough that future work can understand which scripts matter.
35. A simple README/index is added to the active `scripts/` folder explaining each remaining script in plain language.
36. Cleanup does not break active folders, especially `chipseq_pipeline_v2`, `chipseq_pipeline_v2_realstudy`, and their Snakemake workflows.
37. Any path changed by cleanup is patched in docs, configs, Slurm scripts, shell scripts, or workflow rules.

Research front 9, goal/task operating system:

38. A dedicated folder is created for goal/tracking markdown files, containing templates or active files for `PLAN.md`, `PROGRESS.md`, `LOG.md`, and `GOAL.md`.
39. Root `AGENTS.md` remains outside that goal/tracking folder.
40. The tracking markdown files explain their roles clearly enough for future Codex runs: `PLAN.md` for task strategy, `PROGRESS.md` for checklist/status, `LOG.md` for attempts/findings/blockers/assumptions, and `GOAL.md` for the current full objective.
41. `AGENTS.md` is updated to clarify what progress updates should look like during `/goals` tasks.
42. The progress-update rule says updates should summarize the current `PLAN.md` step, `PROGRESS.md` checklist state, and important findings/blockers from `LOG.md`.
43. `AGENTS.md` is updated with a repo hygiene rule that prevents helper-script bloat in Snakemake workflows unless there is a clear need.
44. `AGENTS.md` is updated with a rule to avoid adding new config parameters for sweeps unless necessary.
45. `AGENTS.md` is updated with a rule that tests should live in appropriate test folders rather than being scattered through active scripts.

Final validation across fronts:

46. The final deliverable runs local checks without errors.
47. Proof is shown through test output, dry-run output where applicable, generated comparison tables, generated run tables, and generated or documented plot commands.
48. Goal score/checklist reaches all items checked, except cluster-only execution items may be marked "ready but not locally executable" with a clear reason and command to run on cluster.

## Feedback Loop

- Primary verification command: `cd chipseq_pipeline_v2 && pytest tests`
- Fast feedback command: `cd chipseq_pipeline_v2 && python -m py_compile scripts/*.py`
- Workflow shape command: `cd chipseq_pipeline_v2 && snakemake -s Snakefile.py --configfile <config> --dry-run`
- Realstudy verification command: `cd chipseq_pipeline_v2_realstudy && pytest tests`
- Full validation command: `cd chipseq_pipeline_v2 && pytest tests && cd ../chipseq_pipeline_v2_realstudy && pytest tests`
- Scoring method: checklist of all success criteria above, plus pass/fail for tests, dry-runs, run-table generation, and artifact existence.
- Baseline to compare against: current feature branch after controlled-sweep bug fixes; old six TF-clean outputs under `chipseq_pipeline_v2/analysis_outputs/tfclean_balanced_288_config_report_20260505/`; current realstudy manifests and downloaded worm study files.

## Tracking Files

Create or update these files while working:

- `PLAN.md`: high-level plan, task breakdown, and current strategy.
- `PROGRESS.md`: progress on `PLAN.md`, short checklist with minor notes.
- `LOG.md`: every important attempt, finding, assumption, edge case, and context for future runs.

## Operating Rules

1. PLAN FIRST. Output a numbered task list before writing any code. Update `PLAN.md` when needed.
2. WORK AUTONOMOUSLY. Don't ask clarifying questions unless genuinely blocked.
3. SELF-VERIFY. After every step: run tests, inspect output, confirm it worked.
4. DEBUG YOURSELF. If it fails, diagnose and fix. Don't hand it back.
5. USE EVERY TOOL. MCPs, terminal, web, code execution, pull real data.
6. NO PLACEHOLDERS. No TODOs, no stubs, real components and real states.
7. PROGRESS LOG. Track completed, in-flight, decisions, and blockers in `LOG.md`.
8. STAY ON GOAL. Discoveries off-spec? Note them and keep moving.
9. IF BLOCKED. Log the wall and continue everything parallelizable.
10. CHECK SUCCESS BEFORE STOPPING. Re-read criteria and confirm each is met.
11. KEEP THE LOOP TIGHT. Prefer the fastest useful test before running slow full validation.
12. MAKE PROGRESS MEASURABLE. Convert vague improvements into numbers, checklists, or pass/fail criteria.
13. PERSIST CONTEXT. Keep tracking files updated so work can continue across long runs.
14. AVOID FLAILING. If an approach fails twice, log why, change strategy, and continue.

## Quality Bar

- Code: clean, typed where useful, follows project conventions, and avoids unnecessary helper-script bloat.
- Design: workflow structure is simple, explicit, and easy for future Codex runs or lab members to resume.
- Output: survives a senior code review and a scientific-method review.
- Docs: every new pattern, env var, simulator choice, workflow assumption, and blocked cluster-only step is logged.
- Verification: fast checks pass first, then full local validation passes before stopping; cluster-only EPIC2 execution is documented with dry-run proof where local execution is impossible.

## Final Deliverable

- Confirmation each criterion is satisfied.
- Final score/checklist status.
- Every file created or modified.
- How to run, test, and deploy.
- Proof: test output, dry-run output, generated run tables, simulator comparison table, and plot artifacts or plot commands.
- Decisions made and anything to know.
- Known limitations and follow-ups, especially cluster-only EPIC2 execution and simulator limitations.
- Summary of failed attempts and why the final approach worked.

Begin by outputting your plan. Then execute end-to-end without checking in until done or genuinely blocked.
