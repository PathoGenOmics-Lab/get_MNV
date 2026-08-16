<!--
Fill in what applies and delete what does not. A short pull request needs a
short description; nobody is asking for an essay to fix a typo.
-->

## What this changes

<!-- One or two sentences. What is different after this is merged? -->

## Why

<!-- The problem it solves. Link an issue with "Closes #123" if there is one. -->

## How it was verified

<!--
The important part, and the one a reviewer cannot reconstruct.

Not "it should work", but what you ran and what it printed. For example: the
test you added and the fact that it fails without the fix, the VCF you ran
through it and the annotation you got before and after, or the comparison run
whose output changed.

If it is a change that CI already covers, say which check covers it.
-->

## Checklist

- [ ] `cargo test` passes and `cargo clippy` is clean.
- [ ] New behaviour has a test, or I have said above why it does not.
- [ ] Documentation under `docs/` is updated if the change is user-visible.
- [ ] The CLI and the desktop app still agree on defaults, if I touched either.
