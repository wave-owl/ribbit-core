# Source SBOM

Source repository inventory only; completeness is UNKNOWN. Syft inventories tracked source/package manifests. Gitlinks and declared vendored/model evidence are inspected at immutable pins. CMake declarations are conditional requirements, not resolved installations. This is NOT proof of operating-system packages, static-link closure, build-time resolution, or the contents of release binaries.

Generation: **success**

## Source identity

wave-owl/ribbit-core `1490b479f05d6f72b740da82f947353192261740`

- sbom:role: source-root
- sbom:git-commit: 1490b479f05d6f72b740da82f947353192261740
- sbom:git-tree: 43d7ebadbec603b68509be901d42a2789856aac3
- sbom:source-repository: https://github.com/wave-owl/ribbit-core
- sbom:release-tags: ["v0.1.0"]
- sbom:license-evidence: [{"path": "LICENSE", "sha256": "7ce48f319a0642e30f2f6eca7c9987e5ecb2114bdd49298fa531fc81c9d8b3f8"}]
- sbom:license-resolution: unresolved; manual review required

## Generation errors

None.

## Components and unresolved licenses

| Component | Version | Role | License |
| --- | --- | --- | --- |
| wave-owl/ribbit-core | 1490b479f05d6f72b740da82f947353192261740 | source-root | unresolved; manual review required |
| wave-owl/rattlegram | 4df0f3af3c6da63f0d70512392e80c7a31b72099 | reference-only | unresolved; manual review required |
| src | 4df0f3af3c6da63f0d70512392e80c7a31b72099 | declared-vendored-source | unresolved; manual review required |

License evidence paths/hashes, source manifests, model hashes and conditional
CMake declarations are recorded in SBOM.cdx.json. Declared upstream pins do not
prove vendored equivalence. No dependency installation or build was performed.
