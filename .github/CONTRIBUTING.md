
# Contributing to snpick

Thank you for your interest in contributing to **snpick**! We welcome contributions that help improve the tool, fix bugs, or enhance the documentation. This document outlines the process for contributing and the guidelines to follow.

## Table of Contents
- [Project Description](#project-description)
- [How to Contribute](#how-to-contribute)
  - [Reporting Issues](#reporting-issues)
  - [Feature Requests](#feature-requests)
  - [Contributing Code](#contributing-code)
- [Coding Guidelines](#coding-guidelines)
- [Setting Up the Development Environment](#setting-up-the-development-environment)
- [Submitting a Pull Request](#submitting-a-pull-request)
- [Code of Conduct](#code-of-conduct)

## Project Description
`snpick` is a tool designed to identify **Multi-Nucleotide Variants (MNVs)** within the same codon in genomic sequences. MNVs occur when multiple Single Nucleotide Variants (SNVs) are present within the same codon, leading to the translation of a different amino acid. This tool addresses limitations in current annotation programs like **ANNOVAR** or **SnpEff**, which are primarily designed to work with individual SNVs and might overlook the actual amino acid changes resulting from MNVs.

### Current Limitations
**IMPORTANT**: This script currently works only with **SNVs** against a reference genome. Insertions and deletions that modify the reading frame are not supported yet.

## How to Contribute
We appreciate all contributions, whether it’s fixing bugs, proposing new features, improving the documentation, or suggesting a new direction for the tool.

### Reporting Issues
If you encounter a bug, have a question, or want to request a feature, please [open an issue](https://github.com/PathoGenOmics-Lab/snpick/issues) on our GitHub repository. When reporting an issue, please include:
- A detailed description of the problem.
- Steps to reproduce the issue.
- Any relevant logs or screenshots.
- Version information of `snpick` and your operating system.

### Feature Requests
We welcome suggestions for new features and improvements! Please open an issue labeled **Feature Request** and provide as much detail as possible regarding your suggestion and its potential use cases.

### Contributing Code
If you’d like to contribute code, follow these steps:
1. Fork the repository.
2. Create a new branch for your feature or bugfix (`git checkout -b feature/new-feature`).
3. Make your changes.
4. Test your code.
5. Submit a pull request following the guidelines in the **Submitting a Pull Request** section.

## Coding Guidelines
- **Code Style**: Follow Rust’s official [Rustfmt style guide](https://github.com/rust-lang/rustfmt) for formatting.
- **Testing**: Ensure that your code changes include relevant tests. Use `cargo test` to run tests locally before submitting a pull request.
- **Documentation**: Document all public methods, structs, and modules using doc comments (`///`).

## Setting Up the Development Environment
1. Clone the repository:
   ```bash
   git clone https://github.com/PathoGenOmics-Lab/snpick.git
   cd snpick
   ```

2. Install the necessary dependencies using `cargo`:
   ```bash
   cargo build
   ```

3. Run the tests to ensure everything is set up correctly:
   ```bash
   cargo test
   ```

4. You’re ready to start contributing!

## Submitting a Pull Request
1. Ensure that your code follows the coding guidelines and passes all tests.
2. Write a clear commit message detailing what your change does.
3. Submit a pull request (PR) and fill in the PR template. Include a summary of the changes, why they are necessary, and any relevant issue numbers.
4. A project maintainer will review your PR and provide feedback if needed.

## Code of Conduct
This project follows the [Contributor Covenant Code of Conduct](https://www.contributor-covenant.org/version/2/0/code_of_conduct/). By participating, you agree to abide by its terms. Please be respectful and professional in all interactions.

We look forward to your contributions!
