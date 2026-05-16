# Contribution Guidelines

## Coding Standards
- Only use the libraries that the codebase already uses. You may use other
  libraries only when explicitly instructed to do so in the prompt.
- Do not use the logging module or any logging calls unless explicitly directed.

## Restricted Areas
- Do not edit anything outside of this repository root unless explicitly instructed.

## Environment Setup
- Run `source snakemake_stuff/setup.sh` before executing tests or scripts.
- On the local Mac, if `setup.sh` points to the cluster Conda path, use
  `eval "$(/opt/anaconda3/bin/conda shell.zsh hook)"` followed by
  `conda activate background_project`.
- For local Snakemake commands inside the sandbox, set
  `HOME=/private/tmp/chipseq_snakemake_home` so Snakemake writes its source
  cache somewhere writable.

## Style Guidelines
- At start of each program, ensure sections are marked with triple-quoted comments and leave one blank line above and below.
  - Example: `"""Imports"""`
  - Example: `"""Functions"""`
- Use one line docstrings for functions and explain inputs and outputs when complex.
- Use `argparse` to handle user arguments.
- Never add dunder name/main checks (`if __name__ == '__main__'`) unless explicitly requested.

## Documentation
- Update the relevant `README.md` whenever a PR requires changes to its content.

## Interaction Instructions
- Make a reasonable assumption and keep moving unless a decision would create
  hidden risk or destructive changes.

## Goal Tracking
- Use `goals/current/GOAL.md` as the full objective, `goals/current/PLAN.md`
  as the working strategy, `goals/current/PROGRESS.md` as the checklist, and
  `goals/current/LOG.md` as the durable record of attempts, findings,
  assumptions, blockers, and validation proof.
- Keep `AGENTS.md` at the repo root. Do not move it into the goal tracking
  folder.
- During `/goal` work, progress updates should state the current `PLAN.md`
  step, the relevant `PROGRESS.md` checklist state, and any important
  finding/blocker from `LOG.md`.
- Update `LOG.md` whenever an approach fails, an assumption matters, or a
  cluster-only/local-only limitation is discovered.

## Repo Hygiene
- Keep active workflows obvious. Archive low-utility historical scripts instead
  of leaving them mixed into active workflow paths.
- Avoid helper-script bloat in Snakemake workflows. Add a new helper script
  only when it removes real repetition or makes a workflow easier to verify.
- Avoid adding new sweep config parameters unless they are needed for the
  scientific question or required by a tool.
- Put tests in the relevant `tests/` folder. Do not scatter ad hoc test scripts
  through active script directories.
