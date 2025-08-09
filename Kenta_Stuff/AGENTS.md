# Contribution Guidelines

## Coding Standards
- Only use the libraries that the codebase already uses. You may use other
  libraries only when explicitly instructed to do so in the prompt.
- Do not use the logging module or any logging calls unless explicitly directed.

## Restricted Areas
- Do not edit anything outside of `/Kenta_Stuff`.

## Environment Setup
- Run `source Kenta_Stuff/snakemake_stuff/setup.sh` before executing tests or scripts.

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
- Always ask clarifying questions if unsure.
