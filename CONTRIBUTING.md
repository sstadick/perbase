# Contributing to Perbase

## Common Commands

`perbase` follows standard Rust project structure and commands:

```
# Run tests
cargo test

# Build binary in debug mode
cargo build

# Build in release mode
cargo build --release

# Format code
cargo fmt

# Lint code 
cargo clippy
```

## Submitting PRs

For all PRs, please include an update to the changelog (add a new section if needed).
Any significant changes should include tests, if the change is not already covered by tests.
Any significant changes to usage or output should be checked against the README docs to make sure no updates are needed there as well.

## AI Policy

This project adopts the policy outlined in [ghotstty](https://github.com/ghostty-org/ghostty/blob/7c7d5421f9380f7ca65a0b413d3117d5a88e2dd1/AI_POLICY.md).

In short:

- AI usage fine, but must be disclosed.
- AI usage for discussions must have humans in the loop.
- Any AI generated code must be reviewed in full by the submitter.
