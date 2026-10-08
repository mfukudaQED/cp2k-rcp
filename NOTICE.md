# Attribution and release notes

## Upstream CP2K

This distribution provides modifications against **CP2K 2026.2**, derived from
[the CP2K source repository](https://github.com/cp2k/cp2k). CP2K source headers
carry the following license designation:

- `SPDX-License-Identifier: GPL-2.0-or-later`
- Copyright statements and author notices remain in the patch and in the
  corresponding patched source files.

The [LICENSE](LICENSE) file is copied verbatim from the CP2K source checkout
used to prepare this research extension. No additional proprietary permission
is claimed for upstream CP2K code.

## RCP extension

The RCP patch, manuals, inputs, and supporting utility scripts are provided
for reproducibility and scientific use without any warranty of correctness
for unsupported methods or untested materials.

This is an **independent research extension**: the existence of this repository
does not imply that the RCP extension is part of the official CP2K release.

The isolated experimental Pauli-spinor unit-test kernel in the complete patch
does **not** provide integrated noncollinear or relativistic RCP functionality.

## Research-data separation

Only curated, small input decks and supporting text/code are included.
Calculation outputs and numerical research datasets are not part of this
distribution. Large or provenance-dependent surface/chain simulations should
be evaluated and released separately as needed.
